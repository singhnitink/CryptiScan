#!/usr/bin/env python3
"""
cryptic_pocket_analysis.py  --  CryptiScan Stage 9
==================================================
Quantify the CRYPTIC POCKET ONLY across generative (aSAM) ensembles, and report
how much its volume changes in the mutant relative to the wild type.

The pocket of interest is pinned to the predicted cryptic residues in two ways:

  1. fpocket is run with -P, which guides alpha-sphere clustering to collect the
     spheres in contact with those residues. The pocket is therefore DEFINED by
     the cryptic site rather than found generically and matched afterwards.
     (Disable with --no-guided if your fpocket build predates -P.)

  2. Among the pockets fpocket reports, the cryptic pocket is identified by
     residue overlap with the cryptic site AND proximity of its centroid to the
     cryptic-residue centroid. Both criteria must pass, so a large neighbouring
     pocket that merely brushes one cryptic residue is rejected.

A frame in which no pocket satisfies both criteria is a genuinely CLOSED frame
and enters the statistics with volume 0 -- not as missing data. That is what
makes the open fraction meaningful.

fpocket treats every structure independently, so NO TRAJECTORY ALIGNMENT IS
NEEDED. aSAM frames are independent samples, which matches this exactly. The
frame index is a label, not a time: results are reported as distributions, and
because the samples are independent the standard error is sigma/sqrt(N) with no
block averaging or autocorrelation correction.

Usage
-----
Mutant vs wild type (what you want):

    python cryptic_pocket_analysis.py \
        --mutant-top mut/ensemble_output.top.pdb \
        --mutant-dcd mut/ensemble_output.traj.dcd \
        --wt-top      wt/ensemble_output.top.pdb \
        --wt-dcd      wt/ensemble_output.traj.dcd \
        --residues 200 205 233 235 255 \
        --outdir pocket_analysis --nproc 8

Single ensemble (no comparison):

    python cryptic_pocket_analysis.py \
        --mutant-top ... --mutant-dcd ... --residues ... --outdir ...

RESIDUE NUMBERING: aSAM renumbers residues from 0, so the numbers in
ensemble_output.top.pdb are NOT the original PDB numbers. This script prints the
residue NAME at every number you pass -- verify them before trusting anything.
"""

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import mdtraj as md


# ===========================================================================
# fpocket output parsing
# ===========================================================================
def parse_info(info_path):
    """fpocket <name>_info.txt -> {pocket_number: {descriptor: value}}."""
    pockets, cur = {}, None
    with open(info_path) as fh:
        for line in fh:
            line = line.strip()
            m = re.match(r"^Pocket\s+(\d+)\s*:", line)
            if m:
                cur = int(m.group(1))
                pockets[cur] = {}
                continue
            if cur is not None and ":" in line:
                key, val = line.split(":", 1)
                try:
                    pockets[cur][key.strip()] = float(val.strip())
                except ValueError:
                    pass
    return pockets


def read_target_residues(residues_json):
    """Which residues to analyse, in ORIGINAL PDB numbering.

    Accepts either file the pipeline writes, and prefers the manifest:

      <PDB>_manifest.json      the residues the run actually EXAMINED, each
                               carrying target_mut when it was mutated. This is
                               the correct source: since the pipeline walks down
                               the crypticity ranking past residues whose wild
                               type is already charged/polar, the residues that
                               were mutated are NOT necessarily the top-N by
                               score.

      <PDB>_<chain>_cryptic.json   the raw prediction, top-N by score only.
                               Fine when nothing was skipped, wrong otherwise,
                               so it warns.
    """
    data = json.loads(Path(residues_json).read_text())

    if "cryptic_residues" in data:
        entries = data["cryptic_residues"]
        mutated = [r for r in entries if r.get("target_mut")]
        if mutated:
            skipped = len(entries) - len(mutated)
            note = (f"manifest: {len(mutated)} mutated residue(s)"
                    + (f", {skipped} examined-but-skipped excluded" if skipped else ""))
            return [(r["resid"], r.get("resname_3letter", "?")) for r in mutated], note
        return ([(r["resid"], r.get("resname_3letter", "?")) for r in entries],
                "manifest: no residue carries target_mut, using all examined")

    if "top_residues" in data:
        return ([(r["resid"], r.get("resname_3letter", "?"))
                 for r in data["top_residues"]],
                "cryptic prediction: top-N by score. WARNING -- if the run "
                "skipped any residue these are NOT the residues that were "
                "mutated; pass <PDB>_manifest.json instead")

    sys.exit(f"ERROR: {residues_json} has neither 'cryptic_residues' "
             f"(manifest) nor 'top_residues' (cryptic prediction).")


def map_pdb_to_asam_resids(residues_json, source_pdb, chain, force=False):
    """Translate the pipeline's cryptic residues into aSAM numbering.

    The pipeline reports residues by their ORIGINAL PDB number; aSAM renumbers
    residues from 0 in sequence order. The i-th residue of the chain that was
    handed to aSAM therefore becomes resid i. This walks the input chain in
    order and returns the index of each cryptic residue, so the pipeline never
    has to hard-code an offset.
    """
    # Walk the chain in file order. Residues are keyed on (number, insertion
    # code): 100, 100A and 100B are three residues, and collapsing them would
    # shift every later index by one and silently measure the wrong pocket.
    order, seen = [], set()
    with open(source_pdb) as fh:
        for line in fh:
            if not line.startswith("ATOM") or len(line) < 27:
                continue
            ch = line[21]
            if ch != chain and ch != " ":
                continue
            try:
                rid = int(line[22:26])
            except ValueError:
                continue
            icode = line[26].strip()
            if (rid, icode) not in seen:
                seen.add((rid, icode))
                order.append((rid, icode, line[17:20].strip()))

    if not order:
        sys.exit(f"ERROR: no residues read from {source_pdb} chain '{chain}'.")

    # First occurrence of each residue NUMBER wins; the cryptic predictor
    # reports bare numbers, so an insertion-coded duplicate is ambiguous.
    index, ambiguous = {}, set()
    for i, (rid, icode, _) in enumerate(order):
        if rid not in index:
            index[rid] = i
        else:
            ambiguous.add(rid)

    wanted, source_note = read_target_residues(residues_json)
    if not wanted:
        sys.exit(f"ERROR: {residues_json} lists no residues.")
    print(f"Residue source: {source_note}")

    mapped, report, mismatched = [], [], []
    for rid, name in wanted:
        if rid not in index:
            sys.exit(f"ERROR: cryptic residue {name}{rid} is not present in "
                     f"{source_pdb} chain {chain}.")
        i = index[rid]
        found = order[i][2]
        flag = ""
        if found != name:
            mismatched.append((name, rid, found, i))
            flag = f"   <-- structure has {found}"
        if rid in ambiguous:
            flag += "   <-- insertion codes share this number; used the first"
        mapped.append(i)
        report.append(f"{name}{rid} -> aSAM resid {i}{flag}")

    print(f"Cryptic residue mapping  ({Path(source_pdb).name} -> aSAM numbering):")
    for line in report:
        print(f"    {line}")

    if mismatched:
        print("\nERROR: residue identities do not match the cryptic prediction.")
        print("  The usual cause is passing the MUTANT structure as --source-pdb.")
        print("  Mapping must use the WILD-TYPE chain (the Stage 2 output,")
        print("  e.g. 1jwp.pdb), whose residue names match the Stage 3 JSON.")
        print("  MODELLER preserves numbering and residue order, so the")
        print("  wild-type chain maps correctly for both ensembles.")
        print("  Override with --force-mapping only if you are certain.")
        if not force:
            sys.exit(1)
        print("  --force-mapping given: continuing anyway.\n")

    return sorted(mapped)


def lining_residues_and_centroid(atm_pdb):
    """Residue numbers lining a pocket, and the centroid of its lining atoms."""
    resids, xyz = set(), []
    with open(atm_pdb) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 54:
                try:
                    resids.add(int(line[22:26]))
                    xyz.append((float(line[30:38]),
                                float(line[38:46]),
                                float(line[46:54])))
                except ValueError:
                    pass
    centroid = np.array(xyz).mean(axis=0) if xyz else None
    return resids, centroid


# ===========================================================================
# per-frame worker
# ===========================================================================
def analyse_frame(task):
    """Run fpocket on one frame and return descriptors of the cryptic pocket."""
    (idx, pdb_path, cryptic, site_centroid, min_overlap,
     max_dist, fpocket_bin, guide_str) = task

    pdb_path = Path(pdb_path)
    out_dir = pdb_path.with_name(pdb_path.stem + "_out")

    closed = {"frame": idx, "status": "closed", "volume": 0.0, "drug_score": 0.0,
              "n_alpha": 0, "total_sasa": 0.0, "apolar_sasa": 0.0,
              "hydrophobicity": np.nan, "overlap": 0, "dist": np.nan,
              # specificity control: every other pocket fpocket found in this
              # frame. Free -- fpocket reports them all whether we ask or not.
              "n_pockets": 0, "other_volume": 0.0, "lining": []}

    # Absolute path, no cwd override. This is the invocation empirically shown
    # to work on this fpocket build (1458 pockets found across 2000 frames);
    # a cwd-relative variant was tried and is NOT used, because the known-good
    # form should not be replaced on the strength of a theory about output-
    # directory resolution.
    cmd = [fpocket_bin, "-f", str(pdb_path)]
    if guide_str:
        cmd += ["-P", guide_str]

    # Capture WHY a call failed. Swallowing fpocket's stderr turns every
    # distinct failure -- bad flag, missing file, unwritable directory,
    # segfault -- into the same opaque "fpocket_failed", which is unactionable
    # and makes remote debugging a guessing game.
    try:
        subprocess.run(cmd, capture_output=True, timeout=600, check=True)
    except subprocess.CalledProcessError as exc:
        err = (exc.stderr or b"").decode("utf-8", "replace").strip()
        why = f"exit {exc.returncode}: {err.splitlines()[-1][:160]}" if err \
              else f"exit {exc.returncode} (no stderr)"
        shutil.rmtree(out_dir, ignore_errors=True)
        return {**closed, "status": "fpocket_failed", "volume": np.nan,
                "error": why}
    except Exception as exc:
        shutil.rmtree(out_dir, ignore_errors=True)
        return {**closed, "status": "fpocket_failed", "volume": np.nan,
                "error": f"{type(exc).__name__}: {str(exc)[:160]}"}

    info = out_dir / f"{pdb_path.stem}_info.txt"
    if not info.is_file():
        shutil.rmtree(out_dir, ignore_errors=True)
        return {**closed, "status": "no_output", "volume": np.nan}

    pockets = parse_info(info)
    all_vol = {p: d.get("Volume", 0.0) for p, d in pockets.items()}
    n_pockets = len(pockets)

    # ---- identify the cryptic pocket: overlap AND proximity ---------------
    best, best_ov, best_dist, best_lining = None, 0, np.nan, set()
    for pnum in pockets:
        # info.txt numbers pockets from 1; the files on disk from 0.
        atm = out_dir / "pockets" / f"pocket{pnum - 1}_atm.pdb"
        if not atm.is_file():
            continue
        resids, centroid = lining_residues_and_centroid(atm)
        ov = len(resids & cryptic)
        if ov < min_overlap or centroid is None:
            continue
        d = float(np.linalg.norm(centroid - site_centroid))
        if d > max_dist:
            continue                      # near the right residues but elsewhere
        # more shared lining residues wins; ties broken by proximity
        if ov > best_ov or (ov == best_ov and best is not None and d < best_dist):
            best, best_ov, best_dist, best_lining = pnum, ov, d, resids

    if best is None:
        shutil.rmtree(out_dir, ignore_errors=True)
        # Every pocket in this frame is a non-cryptic one.
        return {**closed, "n_pockets": n_pockets,
                "other_volume": float(sum(all_vol.values()))}

    d = pockets[best]
    result = {
        "frame": idx,
        "status": "open",
        "volume": d.get("Volume", np.nan),
        "drug_score": d.get("Druggability Score", np.nan),
        "n_alpha": d.get("Number of Alpha Spheres", 0),
        "total_sasa": d.get("Total SASA", np.nan),
        "apolar_sasa": d.get("Apolar SASA", np.nan),
        "hydrophobicity": d.get("Hydrophobicity score", np.nan),
        "overlap": best_ov,
        "dist": best_dist,
        # SPECIFICITY CONTROL: summed volume of every pocket that is NOT the
        # cryptic one. If the mutation opens the cryptic site specifically,
        # this stays flat between wild type and mutant while the cryptic
        # volume rises. Costs nothing -- fpocket already found these.
        "n_pockets": n_pockets,
        "other_volume": float(sum(v for p, v in all_vol.items() if p != best)),
        # which residues line the pocket, for the per-residue persistence map
        "lining": sorted(best_lining),
    }
    shutil.rmtree(out_dir, ignore_errors=True)
    return result


# ===========================================================================
# global integrity controls
# ===========================================================================
def structural_metrics(traj, cv_a=None, cv_b=None):
    """Whole-protein descriptors that answer the obvious objection.

    Burying five polar groups is a deliberate destabilisation, so a reviewer
    will ask whether the 'pocket' is simply local unfolding. These are the
    cheapest possible answers: if the mutant's radius of gyration, backbone
    RMSD and secondary-structure content match the wild type's, the fold is
    intact and the volume change is a real cavity, not disorder.

    Optionally also computes a collective variable: the distance between the
    centroids of two residue groups (the inter-helix distance used as a CV in
    the underlying mutagenesis study), letting the generative ensemble be
    compared directly against published molecular-dynamics results.
    """
    bb = traj.top.select("backbone")
    ca = traj.top.select("name CA")
    out = {}

    # Radius of gyration -- global compactness. A mutant that unfolded would
    # show a shifted, broadened distribution.
    out["rg"] = md.compute_rg(traj) * 10.0                       # nm -> A

    # Backbone RMSD from the input conformation (md.rmsd superposes
    # internally and does not modify the trajectory).
    out["rmsd"] = md.rmsd(traj, traj, frame=0, atom_indices=bb) * 10.0

    # Per-residue RMSF. Computed on a copy so the superposition cannot
    # perturb the coordinates fpocket is about to be given.
    tr = md.Trajectory(traj.xyz.copy(), traj.top)
    tr.superpose(tr, frame=0, atom_indices=bb)
    mean_xyz = tr.xyz[:, ca, :].mean(axis=0)
    out["rmsf"] = np.sqrt(
        ((tr.xyz[:, ca, :] - mean_xyz) ** 2).sum(axis=2).mean(axis=0)) * 10.0
    out["rmsf_resseq"] = np.array([traj.top.atom(i).residue.resSeq for i in ca])

    # Secondary structure content. mdtraj's DSSP needs no external binary.
    try:
        dssp = md.compute_dssp(traj, simplified=True)             # H / E / C
        out["frac_helix"] = (dssp == "H").mean(axis=1)
        out["frac_sheet"] = (dssp == "E").mean(axis=1)
    except Exception as e:
        print(f"    (DSSP unavailable: {type(e).__name__}; skipping "
              f"secondary structure)")
        out["frac_helix"] = np.full(traj.n_frames, np.nan)
        out["frac_sheet"] = np.full(traj.n_frames, np.nan)

    # Optional collective variable: centroid separation of two residue groups.
    if cv_a and cv_b:
        ia = traj.top.select("resSeq " + " ".join(str(r) for r in cv_a))
        ib = traj.top.select("resSeq " + " ".join(str(r) for r in cv_b))
        if ia.size and ib.size:
            ca_ = traj.xyz[:, ia, :].mean(axis=1)
            cb_ = traj.xyz[:, ib, :].mean(axis=1)
            out["cv"] = np.linalg.norm(ca_ - cb_, axis=1) * 10.0
        else:
            print("    (CV groups matched no atoms; skipping)")
            out["cv"] = np.full(traj.n_frames, np.nan)
    else:
        out["cv"] = np.full(traj.n_frames, np.nan)

    return out


# ===========================================================================
# one ensemble
# ===========================================================================
def analyse_ensemble(label, top_pdb, dcd, residues, args, chain="A"):
    traj = md.load(dcd, top=top_pdb, stride=args.stride)
    print(f"\n[{label}] {traj.n_frames} frames / {traj.n_atoms} atoms "
          f"(stride {args.stride})")

    cryptic = set(residues)
    present = {r.resSeq: r.name for r in traj.top.residues if r.resSeq in cryptic}
    print(f"[{label}] cryptic sites in this topology:")
    for r in sorted(cryptic):
        print(f"           resid {r:<6} {present.get(r, '*** NOT FOUND ***')}")
    missing = cryptic - set(present)
    if missing:
        sys.exit(f"ERROR [{label}]: residues {sorted(missing)} absent from "
                 f"{top_pdb}. aSAM renumbers from 0 -- map via sequence position.")

    # Per-frame centroid of the cryptic residues (Angstrom), used to keep the
    # identified pocket physically at the cryptic site.
    idx_cryptic = traj.top.select(
        "resSeq " + " ".join(str(r) for r in sorted(cryptic)))
    centroids = traj.xyz[:, idx_cryptic, :].mean(axis=1) * 10.0

    # SASA of the cryptic residues: cheap, and a pocket-detector-free cross-check
    sasa_all = md.shrake_rupley(traj, mode="residue")
    cols = [i for i, r in enumerate(traj.top.residues) if r.resSeq in cryptic]
    sasa = sasa_all[:, cols].sum(axis=1) * 100.0            # nm^2 -> A^2

    # Global integrity controls (seconds, and they pre-empt the main objection)
    print(f"[{label}] computing structural controls (Rg, RMSD, RMSF, DSSP) ...")
    struct = structural_metrics(traj, args.cv_group_a, args.cv_group_b)

    # fpocket -P string: resnum:icode:chain joined by '.'
    guide = None
    if not args.no_guided:
        guide = ".".join(f"{r}::{chain}" for r in sorted(cryptic))
    if not hasattr(args, "guided_used"):
        args.guided_used = bool(guide)

    tmp = Path(tempfile.mkdtemp(prefix=f"fp_{label}_"))
    try:
        def make_task(i):
            return (i, str(tmp / f"f{i:06d}.pdb"), cryptic, centroids[i],
                    args.min_overlap, args.max_dist, args.fpocket, guide)

        def write_frame(i):
            path = tmp / f"f{i:06d}.pdb"
            try:
                traj[i].save_pdb(str(path))
            except OSError as exc:
                sys.exit(
                    f"\nERROR [{label}]: could not write frame {i} to {tmp}:\n"
                    f"  {exc}\n"
                    f"Scratch space is the usual cause. Set TMPDIR to a "
                    f"filesystem with several GB free and re-run; check with "
                    f"'df -h {tmp.parent}'.")
            return path

        # Frames are written and analysed in BATCHES, not all at once.
        # Writing all N frames up front makes peak scratch usage scale with
        # the ensemble size (2000 frames is of order a gigabyte once fpocket's
        # per-frame output directories are counted). When that filesystem
        # fills, fpocket exits without writing its _info.txt and every
        # remaining frame is recorded as a failure -- which looks identical to
        # a pocket that never opens. Batching caps peak usage at roughly
        # batch_size frames regardless of ensemble length.
        batch = max(args.nproc * 4, 50)
        tasks = [make_task(i) for i in range(traj.n_frames)]

        # ------------------------------------------------------------------
        # Probe -P on ONE frame before committing the whole ensemble to it.
        #
        # fpocket 4.0 advertises -P in its help text but does not implement it:
        # it exits 0, prints "No pocket to refine!" to stderr, and writes no
        # output directory at all. Run unchecked, that turns every frame into
        # status=no_output -- which reads exactly like a protein whose pocket
        # never opens, except that a genuinely closed frame is status=closed
        # with volume 0. A silent tool failure must never be presentable as a
        # null biological result, so probe once and fall back. The working -P
        # arrived in fpocket 4.2.
        # ------------------------------------------------------------------
        if guide and tasks:
            write_frame(0)
            probe = analyse_frame(tasks[0])
            Path(tasks[0][1]).unlink(missing_ok=True)
            if probe["status"] in ("fpocket_failed", "no_output"):
                print(f"[{label}] WARNING: fpocket produced no output with -P on "
                      f"the first frame; this build does not support guided "
                      f"detection (4.0 advertises -P without implementing it, "
                      f"4.2+ works).")
                print(f"[{label}]          Falling back to UNGUIDED detection. "
                      f"The cryptic pocket is still identified by the overlap "
                      f"and distance criteria, so the results remain valid -- "
                      f"-P only narrowed the search.")
                guide = None
                tasks = [t[:-1] + (None,) for t in tasks]
                # The report must state what was actually run, not what was
                # requested, or a fallback run is indistinguishable from a
                # guided one in the archived output.
                args.guided_used = False

        print(f"[{label}] running fpocket "
              f"({'guided -P' if guide else 'unguided'}, {args.nproc} workers, "
              f"batches of {batch}) ...")
        results = []
        with Pool(args.nproc) as pool:
            for start in range(0, len(tasks), batch):
                sl = tasks[start:start + batch]
                for t in sl:
                    write_frame(t[0])
                results.extend(pool.map(analyse_frame, sl, chunksize=4))
                for t in sl:
                    Path(t[1]).unlink(missing_ok=True)
                done = min(start + batch, len(tasks))
                n_ok = sum(1 for r in results
                           if r["status"] in ("open", "closed"))
                print(f"[{label}]   {done}/{len(tasks)} frames   "
                      f"{n_ok} usable, {len(results) - n_ok} failed",
                      flush=True)
                # Fail fast: if the first full batch produced nothing usable,
                # the remaining frames will not either, and an hour spent
                # confirming that helps no one.
                if start == 0 and n_ok == 0:
                    print(f"[{label}]   ABORTING this ensemble: the first "
                          f"{len(sl)} frames all failed. Check scratch space "
                          f"and fpocket stderr before re-running.")
                    break
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    results.sort(key=lambda r: r["frame"])

    # Report the distinct failure reasons, once each, with counts.
    reasons = {}
    for r in results:
        if r.get("error"):
            reasons[r["error"]] = reasons.get(r["error"], 0) + 1
    if reasons:
        print(f"[{label}] fpocket failure reasons:")
        for why, n in sorted(reasons.items(), key=lambda kv: -kv[1])[:5]:
            print(f"[{label}]   {n:>6} x  {why}")
        print(f"[{label}]   command was: {args.fpocket} -f <frame.pdb>"
              + (" -P <residues>" if guide else ""))

    for i, r in enumerate(results):
        r["cryptic_sasa"] = float(sasa[i])
        r["rg"] = float(struct["rg"][i])
        r["rmsd"] = float(struct["rmsd"][i])
        r["frac_helix"] = float(struct["frac_helix"][i])
        r["frac_sheet"] = float(struct["frac_sheet"][i])
        r["cv"] = float(struct["cv"][i])
    return results, struct


def write_csv(results, path):
    cols = ["frame", "status", "volume", "drug_score", "n_alpha", "total_sasa",
            "apolar_sasa", "hydrophobicity", "overlap", "dist", "cryptic_sasa",
            "n_pockets", "other_volume", "rg", "rmsd", "frac_helix",
            "frac_sheet", "cv"]
    with open(path, "w") as fh:
        fh.write("frame,status,volume_A3,drug_score,n_alpha_spheres,"
                 "pocket_total_sasa,pocket_apolar_sasa,hydrophobicity,"
                 "lining_overlap,centroid_dist_A,cryptic_sasa_A2,"
                 "n_pockets_total,noncryptic_volume_A3,rg_A,rmsd_A,"
                 "frac_helix,frac_sheet,cv_A\n")
        for r in results:
            fh.write(",".join(
                f"{r[c]:.3f}" if isinstance(r[c], float) else str(r[c])
                for c in cols) + "\n")


def write_rmsf(struct, path):
    """Per-residue RMSF, for the supplementary flexibility figure."""
    with open(path, "w") as fh:
        fh.write("resSeq,rmsf_A\n")
        for r, v in zip(struct["rmsf_resseq"], struct["rmsf"]):
            fh.write(f"{int(r)},{v:.3f}\n")


def write_persistence(results, cryptic, path):
    """Fraction of frames in which each residue lines the cryptic pocket.

    Descriptive, NOT a validation of the cryptic prediction: the pocket was
    selected by overlap with those residues, so they are guaranteed to score
    highly. Its value is in showing which OTHER residues join the pocket when
    it opens.
    """
    counts, n = {}, len(results)
    for r in results:
        for rid in r.get("lining", []):
            counts[rid] = counts.get(rid, 0) + 1
    with open(path, "w") as fh:
        fh.write("resSeq,frames_lining,fraction,is_cryptic_site\n")
        for rid in sorted(counts, key=lambda k: -counts[k]):
            fh.write(f"{rid},{counts[rid]},{counts[rid] / n:.4f},"
                     f"{int(rid in cryptic)}\n")
    return counts


# ===========================================================================
# statistics
# ===========================================================================
def summarise(results):
    """Volume statistics over the whole ensemble. Closed frames count as 0."""
    status = np.array([r["status"] for r in results])
    vol = np.array([r["volume"] for r in results], dtype=float)
    ok = ~np.isin(status, ["fpocket_failed", "no_output"])
    vol, status = vol[ok], status[ok]
    openmask = status == "open"
    openv = vol[openmask]
    n = vol.size

    s = {
        "n_frames": int(n),
        "n_failed": int((~ok).sum()),
        "n_open": int(openmask.sum()),
        "open_fraction": float(openmask.mean()) if n else np.nan,
        # Mean over ALL frames, closed counted as zero volume. This is the
        # headline quantity: it folds together how often the pocket opens and
        # how large it is when open.
        "mean_volume_all": float(vol.mean()) if n else np.nan,
        "sem_volume_all": float(vol.std(ddof=1) / np.sqrt(n)) if n > 1 else np.nan,
        "mean_volume_open": float(openv.mean()) if openv.size else np.nan,
        "sem_volume_open": (float(openv.std(ddof=1) / np.sqrt(openv.size))
                            if openv.size > 1 else np.nan),
        "median_volume_open": float(np.median(openv)) if openv.size else np.nan,
        "p95_volume_open": float(np.percentile(openv, 95)) if openv.size else np.nan,
        "max_volume": float(openv.max()) if openv.size else np.nan,
        "mean_drug_open": float(np.nanmean([r["drug_score"] for r in results
                                            if r["status"] == "open"]))
                          if openmask.any() else np.nan,
        "frac_druggable": (float(np.nanmean(
            [r["drug_score"] > 0.5 for r in results if r["status"] == "open"]))
            if openmask.any() else np.nan),
        "mean_cryptic_sasa": float(np.mean([r["cryptic_sasa"] for r in results])),
        # --- specificity control -------------------------------------------
        "mean_n_pockets": float(np.mean([r["n_pockets"] for r in results])),
        "mean_noncryptic_volume": float(np.mean([r["other_volume"]
                                                 for r in results])),
        "sem_noncryptic_volume": float(np.std([r["other_volume"]
                                               for r in results], ddof=1)
                                       / np.sqrt(len(results)))
                                 if len(results) > 1 else np.nan,
        # --- global integrity controls --------------------------------------
        "mean_rg": float(np.mean([r["rg"] for r in results])),
        "std_rg": float(np.std([r["rg"] for r in results], ddof=1))
                  if len(results) > 1 else np.nan,
        "mean_rmsd": float(np.mean([r["rmsd"] for r in results])),
        "mean_frac_helix": float(np.nanmean([r["frac_helix"] for r in results])),
        "mean_frac_sheet": float(np.nanmean([r["frac_sheet"] for r in results])),
        "mean_cv": float(np.nanmean([r["cv"] for r in results])),
        "_vol_all": vol,
        "_vol_open": openv,
        "_other": np.array([r["other_volume"] for r in results], dtype=float),
        "_rg": np.array([r["rg"] for r in results], dtype=float),
        "_cv": np.array([r["cv"] for r in results], dtype=float),
    }
    return s


def histogram(values, width=46, bins=12):
    if values.size == 0:
        return ["    (no open frames)"]
    counts, edges = np.histogram(values, bins=bins)
    top = counts.max() if counts.max() else 1
    return [f"    {lo:7.0f}-{hi:7.0f} | {'#' * int(width * c / top)} {c}"
            for c, lo, hi in zip(counts, edges[:-1], edges[1:])]


def compare(wt, mut, outdir):
    """Mutant vs wild type. Frames are independent, so plain two-sample tests
    are valid without correlation-time corrections."""
    dV = mut["mean_volume_all"] - wt["mean_volume_all"]
    fold = (mut["mean_volume_all"] / wt["mean_volume_all"]
            if wt["mean_volume_all"] > 0 else np.inf)
    d_open = mut["open_fraction"] - wt["open_fraction"]

    out = {"delta_mean_volume_A3": dV, "fold_change_volume": fold,
           "delta_open_fraction": d_open,
           # Specificity: the cryptic pocket should move while everything else
           # stays put. A comparable shift here means the mutation loosened the
           # whole surface rather than opening this site.
           "delta_noncryptic_volume_A3": (mut["mean_noncryptic_volume"]
                                          - wt["mean_noncryptic_volume"]),
           "specificity_ratio": (abs(dV) / abs(mut["mean_noncryptic_volume"]
                                               - wt["mean_noncryptic_volume"])
                                 if abs(mut["mean_noncryptic_volume"]
                                        - wt["mean_noncryptic_volume"]) > 1e-6
                                 else np.inf),
           # Integrity: these should be ~0. A large shift means the fold moved.
           "delta_rg_A": mut["mean_rg"] - wt["mean_rg"],
           "delta_frac_helix": mut["mean_frac_helix"] - wt["mean_frac_helix"],
           "delta_frac_sheet": mut["mean_frac_sheet"] - wt["mean_frac_sheet"],
           "delta_cv_A": mut["mean_cv"] - wt["mean_cv"]}

    # A total detection failure leaves nothing to compare. Guard here rather
    # than letting scipy raise: the traceback that results is about array
    # sizes and says nothing about the real problem, which is that fpocket
    # returned no usable frame in one or both ensembles.
    if mut["_vol_all"].size == 0 or wt["_vol_all"].size == 0:
        empty = [n for n, s_ in (("mutant", mut), ("wild type", wt))
                 if s_["_vol_all"].size == 0]
        print()
        print("  " + "!" * 66)
        print(f"  NO USABLE FRAMES for: {', '.join(empty)}")
        print("  Every frame was fpocket_failed or no_output, so fpocket never")
        print("  returned a parseable result. This is a TOOL FAILURE, not a")
        print("  closed pocket: a genuinely closed frame has status 'closed'")
        print("  and volume 0. Statistical comparison is skipped; the exposure,")
        print("  flexibility and fold-integrity columns in the CSVs remain valid.")
        print("  Check a frame by hand with --stride 500 and read fpocket stderr.")
        print("  " + "!" * 66)
        out["comparison_skipped"] = f"no usable frames: {', '.join(empty)}"
        return out

    try:
        from scipy import stats
        u = stats.mannwhitneyu(mut["_vol_all"], wt["_vol_all"],
                               alternative="greater")
        ks = stats.ks_2samp(mut["_vol_all"], wt["_vol_all"])
        # Cliff's delta from the Mann-Whitney U statistic: effect size that does
        # not assume normality, which these skewed volume distributions are not.
        n1, n2 = mut["_vol_all"].size, wt["_vol_all"].size
        cliffs = 2.0 * u.statistic / (n1 * n2) - 1.0
        out.update({"mannwhitney_U": float(u.statistic),
                    "mannwhitney_p_greater": float(u.pvalue),
                    "ks_statistic": float(ks.statistic),
                    "ks_p": float(ks.pvalue),
                    "cliffs_delta": float(cliffs)})
    except ImportError:
        print("  (scipy unavailable -- statistical tests skipped)")
    return out


def make_plot(wt, mut, path):
    # Nothing to draw when an ensemble produced no usable frame: open_fraction
    # is NaN and matplotlib rejects NaN axis limits. The run has already told
    # the user what went wrong by this point, so a traceback here only buries
    # that message under an unrelated one about axis limits.
    if (wt["_vol_all"].size == 0 or mut["_vol_all"].size == 0
            or not np.isfinite(wt["open_fraction"])
            or not np.isfinite(mut["open_fraction"])):
        print(f"  (no usable frames -- {path.name} not written)")
        return False

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print(f"  (matplotlib unavailable -- no plot written; add it to the "
              f"container to enable {path.name})")
        return False

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))

    hi = max(mut["_vol_all"].max(initial=0), wt["_vol_all"].max(initial=0), 1.0)
    bins = np.linspace(0, hi * 1.05, 36)
    ax[0].hist(wt["_vol_all"], bins=bins, alpha=0.65, label="wild type",
               color="#2E7D8C", edgecolor="none")
    ax[0].hist(mut["_vol_all"], bins=bins, alpha=0.65, label="mutant",
               color="#D97D28", edgecolor="none")
    ax[0].set_xlabel("cryptic pocket volume ($\\AA^3$)")
    ax[0].set_ylabel("frames")
    ax[0].set_title("Volume distribution (closed frames at 0)")
    ax[0].legend(frameon=False)

    labels = ["wild type", "mutant"]
    openf = [100 * wt["open_fraction"], 100 * mut["open_fraction"]]
    meanv = [wt["mean_volume_all"], mut["mean_volume_all"]]
    x = np.arange(2)
    b = ax[1].bar(x, openf, color=["#2E7D8C", "#D97D28"], width=0.55)
    ax[1].set_xticks(x); ax[1].set_xticklabels(labels)
    ax[1].set_ylabel("frames with pocket open (%)")
    ax[1].set_title("Open fraction")
    for rect, o, m in zip(b, openf, meanv):
        ax[1].text(rect.get_x() + rect.get_width() / 2, o,
                   f"{o:.1f}%\n$\\langle V\\rangle$={m:.0f} $\\AA^3$",
                   ha="center", va="bottom", fontsize=9)
    ax[1].set_ylim(0, max(openf) * 1.35 + 1)

    for a in ax:
        a.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)
    return True


# ===========================================================================
def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mutant-top", required=True)
    ap.add_argument("--mutant-dcd", required=True)
    ap.add_argument("--wt-top", help="wild-type topology (enables comparison)")
    ap.add_argument("--wt-dcd", help="wild-type trajectory")
    ap.add_argument("--residues", type=int, nargs="*", default=None,
                    metavar="RESID",
                    help="cryptic residue numbers AS THEY APPEAR IN the aSAM "
                         "topology (aSAM renumbers from 0). Omit this and use "
                         "--residues-json/--source-pdb to map automatically.")
    ap.add_argument("--residues-json",
                    help="<PDB>_manifest.json (PREFERRED -- gives the residues "
                         "actually mutated) or <PDB>_<chain>_cryptic.json "
                         "(top-N by score, which differs from the mutated set "
                         "whenever the run skipped a residue). Numbering is "
                         "mapped to aSAM's automatically.")
    ap.add_argument("--source-pdb",
                    help="the WILD-TYPE single-chain PDB from Stage 2 (e.g. "
                         "1jwp.pdb) -- its residue names match the Stage 3 JSON. "
                         "Required with --residues-json.")
    ap.add_argument("--force-mapping", action="store_true",
                    help="continue even if mapped residue identities disagree "
                         "with the cryptic prediction (not recommended)")
    ap.add_argument("--cv-group-a", type=int, nargs="*", default=None,
                    metavar="RESID",
                    help="optional collective variable: residues (aSAM "
                         "numbering) forming the first group. The distance "
                         "between the two group centroids is reported, so the "
                         "generative ensemble can be compared against a "
                         "published MD collective variable.")
    ap.add_argument("--cv-group-b", type=int, nargs="*", default=None,
                    metavar="RESID", help="second CV residue group")
    ap.add_argument("--outdir", default="pocket_analysis")
    ap.add_argument("--chain", default="A", help="chain id in the aSAM PDB")
    ap.add_argument("--stride", type=int, default=1,
                    help="use every Nth frame (try 50 for a quick test)")
    ap.add_argument("--min-overlap", type=int, default=2,
                    help="min cryptic residues lining the pocket (default 2)")
    ap.add_argument("--max-dist", type=float, default=12.0,
                    help="max Angstrom from pocket centroid to cryptic-site "
                         "centroid (default 12)")
    ap.add_argument("--no-guided", action="store_true",
                    help="do not pass -P to fpocket (for builds without it)")
    ap.add_argument("--nproc", type=int, default=4)
    # fpocket lives in its own conda env inside the container and is NOT on
    # PATH, so fall back to the FPOCKET_BIN the image exports. Without this,
    # running the stage by hand fails with "fpocket not on PATH" even though
    # the binary is present.
    ap.add_argument("--fpocket",
                    default=os.environ.get("FPOCKET_BIN", "fpocket"),
                    help="fpocket binary (default: $FPOCKET_BIN, else PATH)")
    args = ap.parse_args()

    if not shutil.which(args.fpocket):
        sys.exit(f"ERROR: '{args.fpocket}' not on PATH.")
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    # Residues: given directly in aSAM numbering, or mapped from the pipeline's
    # Stage 3 output, which uses the original PDB numbering.
    if args.residues_json:
        if not args.source_pdb:
            sys.exit("ERROR: --residues-json requires --source-pdb "
                     "(the structure handed to aSAM).")
        args.residues = map_pdb_to_asam_resids(args.residues_json,
                                               args.source_pdb, args.chain,
                                               force=args.force_mapping)
    if not args.residues:
        sys.exit("ERROR: supply --residues, or --residues-json with --source-pdb.")

    cryptic_set = set(args.residues)

    mut_res, mut_struct = analyse_ensemble("mutant", args.mutant_top,
                                           args.mutant_dcd, args.residues,
                                           args, args.chain)
    write_csv(mut_res, outdir / "mutant_pocket.csv")
    write_rmsf(mut_struct, outdir / "mutant_rmsf.csv")
    write_persistence(mut_res, cryptic_set, outdir / "mutant_persistence.csv")
    mut = summarise(mut_res)

    wt = None
    if args.wt_top and args.wt_dcd:
        wt_res, wt_struct = analyse_ensemble("wild type", args.wt_top,
                                             args.wt_dcd, args.residues,
                                             args, args.chain)
        write_csv(wt_res, outdir / "wildtype_pocket.csv")
        write_rmsf(wt_struct, outdir / "wildtype_rmsf.csv")
        write_persistence(wt_res, cryptic_set,
                          outdir / "wildtype_persistence.csv")
        wt = summarise(wt_res)

    # ---------------------------------------------------------------- report
    lines = []
    A = lines.append
    A("=" * 70)
    A(" CRYPTIC POCKET VOLUME  --  CryptiScan Stage 9")
    A("=" * 70)
    A(f" cryptic residues (aSAM numbering) : "
      f"{', '.join(str(r) for r in sorted(args.residues))}")
    _guided = getattr(args, "guided_used", not args.no_guided)
    A(f" pocket definition                 : "
      f"{'fpocket -P (guided)' if _guided else 'unguided + matching'}"
      + ("" if _guided or args.no_guided
         else "   [-P unsupported by this fpocket build; fell back]"))
    A(f" identification                    : >= {args.min_overlap} lining "
      f"residues and centroid within {args.max_dist:.0f} A")
    A("")

    def block(tag, s):
        A(f" {tag}")
        A(f"   frames analysed            : {s['n_frames']}"
          + (f"   ({s['n_failed']} fpocket failures excluded)"
             if s['n_failed'] else ""))
        A(f"   pocket open                : {s['n_open']}  "
          f"({100 * s['open_fraction']:.1f} %)")
        A(f"   mean volume, ALL frames    : {s['mean_volume_all']:.1f} "
          f"+/- {s['sem_volume_all']:.1f} A^3    <-- headline")
        A(f"   mean volume, open frames   : {s['mean_volume_open']:.1f} "
          f"+/- {s['sem_volume_open']:.1f} A^3")
        A(f"   median / 95th pct (open)   : {s['median_volume_open']:.1f} / "
          f"{s['p95_volume_open']:.1f} A^3")
        A(f"   largest pocket seen        : {s['max_volume']:.1f} A^3")
        A(f"   mean druggability (open)   : {s['mean_drug_open']:.3f}  "
          f"({100 * s['frac_druggable']:.0f} % of open frames > 0.5)")
        A(f"   cryptic-residue SASA       : {s['mean_cryptic_sasa']:.1f} A^2")
        A("   volume distribution (open frames, A^3)")
        lines.extend(histogram(s["_vol_open"]))
        A("   controls")
        A(f"     pockets per frame        : {s['mean_n_pockets']:.1f}")
        A(f"     non-cryptic pocket vol   : {s['mean_noncryptic_volume']:.1f} "
          f"+/- {s['sem_noncryptic_volume']:.1f} A^3")
        A(f"     radius of gyration       : {s['mean_rg']:.2f} "
          f"+/- {s['std_rg']:.2f} A")
        A(f"     backbone RMSD from input : {s['mean_rmsd']:.2f} A")
        A(f"     helix / sheet content    : {100 * s['mean_frac_helix']:.1f} % / "
          f"{100 * s['mean_frac_sheet']:.1f} %")
        if not np.isnan(s["mean_cv"]):
            A(f"     collective variable      : {s['mean_cv']:.2f} A")
        A("")

    if wt:
        block("WILD TYPE", wt)
    block("MUTANT", mut)

    comp = None
    if wt:
        comp = compare(wt, mut, outdir)
        A("-" * 70)
        A(" MUTANT vs WILD TYPE")
        A("-" * 70)
        A(f"   delta mean volume          : {comp['delta_mean_volume_A3']:+.1f} A^3")
        A(f"   fold change                : {comp['fold_change_volume']:.2f} x")
        A(f"   delta open fraction        : "
          f"{100 * comp['delta_open_fraction']:+.1f} percentage points")
        if "mannwhitney_p_greater" in comp:
            A(f"   Mann-Whitney U (mut>wt)    : p = "
              f"{comp['mannwhitney_p_greater']:.3e}")
            A(f"   Kolmogorov-Smirnov         : D = {comp['ks_statistic']:.3f}, "
              f"p = {comp['ks_p']:.3e}")
            A(f"   Cliff's delta (effect size): {comp['cliffs_delta']:+.3f}")
            A("")
            A("   Frames are independent samples, so these two-sample tests")
            A("   apply directly -- no block averaging or correlation-time")
            A("   correction is needed, unlike a molecular-dynamics trajectory.")
        A("")
        A("   SPECIFICITY  (did only the cryptic pocket change?)")
        A(f"     delta non-cryptic volume : "
          f"{comp['delta_noncryptic_volume_A3']:+.1f} A^3")
        A(f"     specificity ratio        : "
          f"{comp['specificity_ratio']:.1f} x   "
          f"(cryptic change / background change; >5 is clean)")
        A("")
        A("   FOLD INTEGRITY  (is the mutant still folded?)")
        A(f"     delta radius of gyration : {comp['delta_rg_A']:+.2f} A")
        A(f"     delta helix content      : "
          f"{100 * comp['delta_frac_helix']:+.1f} pp")
        A(f"     delta sheet content      : "
          f"{100 * comp['delta_frac_sheet']:+.1f} pp")
        A("     These should all be near zero. A large shift means the")
        A("     substitutions perturbed the fold, and the volume change")
        A("     may be local unfolding rather than a genuine cavity.")
        if not np.isnan(comp["delta_cv_A"]):
            A("")
            A(f"   COLLECTIVE VARIABLE")
            A(f"     delta CV                 : {comp['delta_cv_A']:+.2f} A")
        A("")
        A("   NOTE: open fractions from a generative model are Boltzmann-")
        A("   weighted only insofar as the model reproduces its training")
        A("   distribution. Report these as a relative mutant-vs-wild-type")
        A("   comparison, not as thermodynamic populations.")
    else:
        A("-" * 70)
        A(" No wild-type ensemble supplied, so there is no reference. A volume")
        A(" on its own carries no meaning -- rerun aSAM on the pre-mutation")
        A(" structure and pass --wt-top/--wt-dcd.")
    A("=" * 70)

    report = "\n".join(lines)
    print("\n" + report)
    (outdir / "pocket_report.txt").write_text(report + "\n")

    payload = {"cryptic_residues": sorted(args.residues),
               "guided": getattr(args, "guided_used", not args.no_guided),
               "min_overlap": args.min_overlap, "max_dist_A": args.max_dist,
               "mutant": {k: v for k, v in mut.items() if not k.startswith("_")}}
    if wt:
        payload["wildtype"] = {k: v for k, v in wt.items() if not k.startswith("_")}
        payload["comparison"] = comp
    (outdir / "pocket_summary.json").write_text(json.dumps(payload, indent=2))

    if wt:
        make_plot(wt, mut, outdir / "pocket_volume_comparison.png")

    print(f"\nWritten to {outdir}/:")
    print("   mutant_pocket.csv          per-frame volume, controls, structure")
    print("   mutant_rmsf.csv            per-residue flexibility")
    print("   mutant_persistence.csv     residues lining the pocket, by frequency")
    if wt:
        print("   wildtype_pocket.csv        per-frame, reference ensemble")
        print("   wildtype_rmsf.csv          per-residue flexibility")
        print("   wildtype_persistence.csv   pocket lining, reference")
        print("   pocket_volume_comparison.png")
    print("   pocket_report.txt          the text above")
    print("   pocket_summary.json        machine-readable, for batch screening")


if __name__ == "__main__":
    main()
