#!/usr/bin/env python3
"""
diagnose_fpocket.py -- find out why Stage 9 returned no_output on every frame
============================================================================
Run this ON THE MACHINE THAT HAS THE ENSEMBLES, inside the container:

    apptainer exec pipeline_local.sif /opt/conda/bin/python3 \
        diagnose_fpocket.py \
            --top ensemble_output.top.pdb \
            --dcd ensemble_output.traj.dcd \
            --chain A \
            --residues 200 205 233 235 255

It writes ONE frame, runs fpocket on it four ways, and prints the exit status,
stderr and the files produced each time. Takes about ten seconds and tells you
exactly which of the four candidate causes it is.

Pass --residues in aSAM numbering -- the numbers Stage 9 printed under
"Cryptic residue mapping", not the original PDB numbers.
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def run(cmd, cwd=None):
    p = subprocess.run(cmd, capture_output=True, text=True, cwd=cwd, timeout=600)
    return p.returncode, p.stdout, p.stderr


def show(tag, rc, out, err, workdir, stem):
    print(f"\n--- {tag}")
    print(f"    exit code : {rc}")
    if out.strip():
        print("    stdout    : " + out.strip().splitlines()[0][:150])
    if err.strip():
        for line in err.strip().splitlines()[:6]:
            print(f"    stderr    : {line[:150]}")
    outdir = Path(workdir) / f"{stem}_out"
    if outdir.is_dir():
        files = sorted(p.name for p in outdir.iterdir())
        print(f"    {stem}_out/ : {len(files)} entries -> {files[:6]}")
        info = outdir / f"{stem}_info.txt"
        # Record the result BEFORE deleting the directory. Reading is_file()
        # after rmtree always returns False and makes every probe look failed.
        ok = info.is_file()
        print(f"    info.txt  : {'PRESENT' if ok else '*** MISSING ***'}")
        if ok:
            n = sum(1 for l in info.read_text().splitlines()
                    if l.strip().startswith("Pocket"))
            print(f"    pockets   : {n}")
        shutil.rmtree(outdir, ignore_errors=True)
        return ok
    else:
        print(f"    {stem}_out/ : *** NOT CREATED ***")
        # fpocket may have written somewhere unexpected
        stray = [p.name for p in Path(workdir).iterdir() if p.name != f"{stem}.pdb"]
        if stray:
            print(f"    other     : {stray[:8]}")
        return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--top", required=True)
    ap.add_argument("--dcd", required=True)
    ap.add_argument("--chain", default="A")
    ap.add_argument("--residues", type=int, nargs="+", required=True,
                    help="cryptic residues in aSAM numbering")
    ap.add_argument("--frame", type=int, default=0)
    ap.add_argument("--fpocket", default=os.environ.get("FPOCKET_BIN", "fpocket"))
    a = ap.parse_args()

    fp = shutil.which(a.fpocket) or a.fpocket
    print("=" * 72)
    print(" fpocket diagnostic")
    print("=" * 72)
    print(f" binary : {fp}")
    if not Path(fp).exists() and not shutil.which(a.fpocket):
        sys.exit("ERROR: fpocket not found. Set FPOCKET_BIN or --fpocket.")

    # ---- 1. version and whether -P is advertised --------------------------
    rc, out, err = run([fp, "--help"])
    blob = (out + err)
    if not blob.strip():
        rc, out, err = run([fp])
        blob = out + err
    ver = [l for l in blob.splitlines() if "ersion" in l or "fpocket" in l.lower()][:2]
    print(" version banner:")
    for l in ver:
        print("   " + l.strip()[:120])
    has_P = "-P" in blob
    print(f" '-P' advertised in help text : {'YES' if has_P else 'NO'}")
    if not has_P:
        print("   >>> If NO, this build predates the -P option and the guided call")
        print("   >>> fails. Re-run Stage 9 with --no-guided, or upgrade fpocket.")

    # ---- 2. write one frame ----------------------------------------------
    import mdtraj as md
    traj = md.load_frame(a.dcd, a.frame, top=a.top)
    tmp = Path(tempfile.mkdtemp(prefix="fpdiag_"))
    stem = "frame"
    pdb = tmp / f"{stem}.pdb"
    traj.save_pdb(str(pdb))
    print(f"\n wrote {pdb}  ({traj.n_atoms} atoms, {traj.n_residues} residues)")

    # ---- 3. what the written PDB actually contains ------------------------
    chains, resids, names = set(), [], {}
    for line in pdb.read_text().splitlines():
        if line.startswith("ATOM") and len(line) > 26:
            chains.add(line[21])
            try:
                r = int(line[22:26])
            except ValueError:
                continue
            if r not in names:
                names[r] = line[17:20].strip()
                resids.append(r)
    print(f" chain IDs in file  : {sorted(repr(c) for c in chains)}")
    print(f" residue numbering  : {min(resids)} .. {max(resids)}  ({len(resids)} residues)")
    print(" requested cryptic residues:")
    missing = []
    for r in a.residues:
        ok = r in names
        if not ok:
            missing.append(r)
        print(f"    {r:>5}  {'-> ' + names[r] if ok else '*** NOT IN FILE ***'}")
    bad_chain = a.chain not in chains
    if bad_chain:
        print(f"\n >>> CHAIN MISMATCH: you passed --chain {a.chain!r} but the frame")
        print(f" >>> contains {sorted(chains)!r}. The -P string names a chain that")
        print(" >>> does not exist, so fpocket cannot resolve the site.")
    if missing:
        print(f"\n >>> RESIDUES ABSENT: {missing} are not in the frame. The mapping")
        print(" >>> from PDB to aSAM numbering is wrong; re-check --source-pdb.")

    # ---- 4. four invocations ---------------------------------------------
    guide = ".".join(f"{r}::{a.chain}" for r in sorted(a.residues))
    print(f"\n -P string: {guide}")

    results = {}
    results["plain"] = show("A. plain:  fpocket -f frame.pdb",
                            *run([fp, "-f", str(pdb)], cwd=tmp), tmp, stem)
    results["guided"] = show(f"B. guided: fpocket -f frame.pdb -P {guide[:40]}...",
                             *run([fp, "-f", str(pdb), "-P", guide], cwd=tmp),
                             tmp, stem)
    # chain-less form: some builds want resnum::  with no chain
    guide_nochain = ".".join(f"{r}::" for r in sorted(a.residues))
    results["guided_nochain"] = show("C. guided, no chain id",
                                     *run([fp, "-f", str(pdb), "-P", guide_nochain],
                                          cwd=tmp), tmp, stem)
    # relative path, in case an absolute path confuses the output-dir logic
    results["relative"] = show("D. plain, relative path (cwd = frame dir)",
                               *run([fp, "-f", f"{stem}.pdb"], cwd=tmp), tmp, stem)

    # ---- 5. verdict -------------------------------------------------------
    print("\n" + "=" * 72)
    print(" VERDICT")
    print("=" * 72)
    if results["plain"] and not results["guided"]:
        print(" fpocket works WITHOUT -P and fails WITH it.")
        print(" -> The guided call is the problem.")
        if results["guided_nochain"]:
            print(" -> The chain-less -P form works: your frames carry a chain ID")
            print("    that does not match --chain. Fix the chain, or use form C.")
        else:
            print(" -> Neither -P form works. fpocket 4.0 advertises -P but its")
            print("    implementation is incomplete ('No pocket to refine!'); the")
            print("    working one arrived in 4.2.")
            print()
            print("    FIX: re-run Stage 9 with --no-guided. Detection is then")
            print("    unguided and the cryptic pocket is identified by the")
            print("    overlap + distance test alone, which is what --min-overlap")
            print("    and --max-dist are for. Nothing in the analysis depends on")
            print("    -P; it only narrowed the search. Optionally also upgrade:")
            print("      conda install -y -n analysis -c conda-forge 'fpocket>=4.2'")
    elif not results["plain"] and results["relative"]:
        print(" fpocket works only when called with a RELATIVE path.")
        print(" -> Absolute-path handling in this build is broken; have the")
        print("    worker chdir into the frame directory before calling it.")
    elif not any(results.values()):
        print(" fpocket produced no parseable output in ANY form.")
        print(" -> The build itself is broken, or the frame PDB is unreadable by it.")
        print("    Reinstall: conda install -y -c conda-forge fpocket")
        print("    Then re-run case A by hand and read the stderr above.")
    else:
        print(" fpocket works in the guided form on this frame.")
        print(" -> The per-frame failure is elsewhere: check that the worker's")
        print("    temp directory is writable and not being cleaned concurrently,")
        print("    and re-run Stage 9 with --nproc 1 to rule out a race.")
    shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    main()
