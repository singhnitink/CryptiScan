#!/usr/bin/env python3
"""
CryptoBank -> ESM-Scan -> MODELLER Pipeline Orchestrator
=========================================================
Stages:
    1.  clean_structure.py  : Fetch the entry from RCSB and split it with
                              MDAnalysis. Keeps BOTH the full protein
                              (<pdb>_protein.pdb, all chains) and one file per
                              ligand under <pdb>_ligands/. Waters are dropped.
    2.  clean_structure.py  : Select the user-requested chain -> <pdb>.pdb.
                              Every later stage operates on this single chain.
    3.  scrape_cryptobank.py / predict_cryptic_local.py
                            : Find the top cryptic pocket residues.
    4.  esm_scanner.py      : Run ESM-Scan (ESM-1v zero-shot variant scoring) on
                              the top cryptic residues to find the best amino acid
                              mutation for each position (highest fitness score).
    5.  modeller_mutate.py  : Build the relaxed 3D mutant structure using MODELLER.
    5b. clean_structure.py  : Validate the mutant for SAM2 -> sam2_input.pdb.
    6.  manifest.json       : Manifest describing the run, ready for the
                              SAM2 / aSAM ensemble stage.

Usage:
------
    python3 run_pipeline.py --pdb 1JWP --chain A --top 5 \
        --esm_model /opt/models/esm1v_t33_650M_UR90S_1.pt \
        --esm_python /opt/conda/bin/python3 \
        --mode combined \
        --check_numbering
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

VALID_AA3 = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
}

# Auto-detect MODELLER install directory if not explicitly passed
if "MODINSTALL10v8" not in os.environ:
    import glob
    _mod_dirs = glob.glob("/opt/conda/lib/modeller-*")
    if _mod_dirs:
        os.environ["MODINSTALL10v8"] = _mod_dirs[0]


def log(msg):
    print(f"[PIPELINE] {msg}", file=sys.stderr)


# ---------------------------------------------------------------------------
# STAGE 3: cryptic residue prediction (CryptoBank web or local ProtT5)
# ---------------------------------------------------------------------------
def get_cryptic_residues(pdb, chain, top, score_type, workdir, scraper, python_exe):
    """Queries CryptoBank API to identify top cryptic residues."""
    out_json = workdir / f"{pdb}_{chain}_cryptic.json"
    cmd = [python_exe, str(scraper), "--pdb", pdb, "--chain", chain,
           "--top", str(top), "--score_type", score_type, "--output", str(out_json)]
    log(f"Stage 3: cryptic residues -> {out_json.name}")
    log(f"         scraper python: {python_exe}")
    subprocess.run(cmd, check=True)
    data = json.loads(out_json.read_text())
    residues = data.get("top_residues", [])
    log(f"         got {len(residues)} cryptic residues: "
        + ", ".join(f"{r['resname_3letter']}{r['resid']}" for r in residues))
    return residues


# ---------------------------------------------------------------------------
# STAGES 1 & 2: fetch, split protein/ligands, then select the chain
# ---------------------------------------------------------------------------
def fetch_and_split(pdb, workdir, do_clean, prep_script, prep_python):
    """STAGE 1. Fetches the entry from RCSB and splits it with MDAnalysis into
    workdir/<pdb>_protein.pdb  (all chains, standard residues) and
    workdir/<pdb>_ligands/<RES>_<CH><N>.pdb (one file per ligand, any resname).
    Waters are dropped. Returns the all-chain protein file."""
    base = pdb.lower()
    protein = workdir / f"{base}_protein.pdb"
    if protein.exists():
        log(f"Stage 1: using existing {protein.name} (not overwriting)")
        return protein

    raw = workdir / f"{base}_raw.pdb"
    src = Path(pdb)
    if src.is_file():
        shutil.copy(str(src), str(raw))
        log(f"Stage 1: using local structure {src.name} -> {raw.name}")
    elif raw.exists():
        log(f"Stage 1: reusing already-downloaded {raw.name}")
    else:
        url = f"https://files.rcsb.org/download/{pdb.upper()}.pdb"
        log(f"Stage 1: fetching {url}")
        try:
            urllib.request.urlretrieve(url, raw)
        except Exception as e:
            log(f"  ERROR downloading PDB: {e}")
            log(f"  Place the file manually at: {raw}")
            sys.exit(1)

    if not do_clean:
        shutil.copy(str(raw), str(protein))
        log(f"Stage 1: skipping cleanup -> {protein.name} (--no_clean)")
        return protein

    log(f"Stage 1: splitting protein/ligands with MDAnalysis")
    subprocess.run([prep_python, str(prep_script),
                    "--mode", "split",
                    "--pdb", str(raw),
                    "--outdir", str(workdir),
                    "--basename", base], check=True)
    if not protein.exists():
        log(f"  ERROR: structure preparation did not produce {protein.name}")
        sys.exit(1)
    return protein


def select_chain(protein_all, chain, pdb, workdir, do_clean, prep_script, prep_python):
    """STAGE 2. Reduces the structure to the single chain the user selected and
    validates it (standard residues only, every residue has a CA). This one
    chain is what CryptoBank, ESM-Scan, MODELLER and SAM2 all operate on."""
    dest = workdir / f"{pdb.lower()}.pdb"
    if dest.exists():
        log(f"Stage 2: using existing {dest.name} (not overwriting)")
        return dest

    if not do_clean:
        shutil.copy(str(protein_all), str(dest))
        log(f"Stage 2: skipping chain selection -> {dest.name} (--no_clean)")
        return dest

    log(f"Stage 2: selecting chain {chain} -> {dest.name}")
    subprocess.run([prep_python, str(prep_script),
                    "--mode", "extract",
                    "--pdb", str(protein_all),
                    "--chain", chain,
                    "--output", str(dest)], check=True)
    return dest


def validate_for_sam2(mutant_pdb, chain, workdir, prep_script, prep_python):
    """Re-checks the MODELLER output against SAM2's requirements before the
    ensemble stage. The structure is already single-chain by this point, so
    this is a validation pass rather than a filter."""
    dest = workdir / "sam2_input.pdb"
    log(f"Stage 5b: validating mutant for SAM2 -> {dest.name}")
    subprocess.run([prep_python, str(prep_script),
                    "--mode", "extract",
                    "--pdb", str(mutant_pdb),
                    "--chain", chain,
                    "--output", str(dest)], check=True)
    return dest


def check_numbering(pdb_file, chain, residues):
    """Verifies that residue IDs from CryptoBank match residues present in the PDB."""
    present = {}
    for line in Path(pdb_file).read_text().splitlines():
        if line.startswith("ATOM") and len(line) > 25 and line[21] == chain:
            try:
                present[int(line[22:26])] = line[17:20].strip()
            except ValueError:
                continue
    mism = []
    for r in residues:
        got = present.get(r["resid"])
        if got is None:
            mism.append(f"  resid {r['resid']} ({r['resname_3letter']}) NOT FOUND in chain {chain}")
        elif got != r["resname_3letter"]:
            mism.append(f"  resid {r['resid']}: CryptoBank={r['resname_3letter']}, PDB={got}")
    if mism:
        log("WARNING: residue-numbering mismatches:")
        for m in mism:
            log(m)
        log("  -> numbering may be offset; resolve before trusting mutations.")
        return False
    log(f"Stage 3: numbering check OK ({len(residues)} residues matched in chain {chain})")
    return True


# ---------------------------------------------------------------------------
# STAGE 4: ESM-Scan mutation scoring (select best mutation per residue)
# ---------------------------------------------------------------------------
def run_esm_scan(pdb_file, chain, residues, workdir, esm_script, esm_model, esm_python,
                 mutation_set="charged_polar"):
    """
    Uses ESM-Scan to score substitutions at each cryptic position and select the
    highest-scoring one from --mutation_set (charged/polar by default, so the
    mutation actually destabilises the closed state).

    A residue whose wild type is already charged/polar is marked skipped: it is
    scored for the report but not mutated.
    """
    res_list_str = ",".join(str(r["resid"]) for r in residues)
    out_json = workdir / f"{pdb_file.stem}_esm_results.json"
    out_csv = workdir / f"{pdb_file.stem}_esm_summary.csv"

    cmd = [
        esm_python, str(esm_script),
        "--pdb", str(pdb_file),
        "--chain", chain,
        "--residues", res_list_str,
        "--model_path", str(esm_model),
        "--mutation_set", str(mutation_set),
        "--output_json", str(out_json),
        "--output_csv", str(out_csv)
    ]
    log("Stage 4: Running ESM-Scan (predicting best mutation for each residue)...")
    log(f"         Model: {esm_model}")
    log(f"         Mutation set: {mutation_set}")
    log(f"         ESM Python: {esm_python}")
    subprocess.run(cmd, check=True)

    esm_data = json.loads(out_json.read_text())
    esm_map = {item["resid"]: item for item in esm_data}

    # Attach the best mutation to each residue in the list
    for r in residues:
        resid = r["resid"]
        info = esm_map.get(resid)

        if info is None:
            # Never invent a substitution here. This used to default to ALA,
            # which is nonpolar -- the exact opposite of what the restricted
            # set is for -- and it happened silently.
            log(f"WARNING: no ESM score for resid {resid}; it will NOT be mutated.")
            r["target_mut"] = None
            r["esm_score"] = None
            r["skip_mutation"] = True
            r["skip_reason"] = "no ESM score returned for this residue"
            continue

        r["esm_candidates"] = info.get("candidates_ranked", [])
        r["esm_all_candidates"] = info.get("all_candidates_ranked", [])
        r["selection_set"] = info.get("selection_set")

        if info.get("skip_mutation"):
            reason = ("wild type is already in the selection set"
                      if info.get("wt_in_selection_set")
                      else "no candidate available in the selection set")
            log(f"         Residue {r['resname_3letter']}{resid} -> SKIPPED ({reason})")
            r["target_mut"] = None
            r["esm_score"] = None
            r["skip_mutation"] = True
            r["skip_reason"] = reason
            continue

        r["target_mut"] = info["best_mut_aa3"]
        r["esm_score"] = info["best_score"]
        r["skip_mutation"] = False
        log(f"         Residue {r['resname_3letter']}{resid} -> Best ESM mutation: "
            f"{r['target_mut']} (score: {r['esm_score']:+.4f})")

    return residues


# ---------------------------------------------------------------------------
# STAGE 4: MODELLER mutations
# ---------------------------------------------------------------------------
def run_modeller(modeller_script, base, respos, restyp, chain, workdir, modeller_python):
    """Executes modeller_mutate.py to introduce a single mutation and refine sidechain."""
    cmd = [modeller_python, str(modeller_script), base, str(respos), restyp, chain]
    logf = workdir / f"modeller_{restyp}{respos}.log"
    log(f"  mutate {chain}/{respos} -> {restyp}   (log: {logf.name})")
    with open(logf, "w") as lf:
        subprocess.run(cmd, cwd=str(workdir), check=True, stdout=lf, stderr=subprocess.STDOUT)
    out = workdir / f"{base}{restyp}{respos}.pdb"
    if not out.exists():
        raise FileNotFoundError(f"MODELLER did not produce {out} (see {logf})")
    return out


def mutate_independent(modeller_script, pdb_base, residues, chain, workdir, modeller_python):
    """Builds N separate single-point mutant structures."""
    mutants = []
    for r in residues:
        respos, wt, target_mut = r["resid"], r["resname_3letter"], r["target_mut"]
        out = run_modeller(modeller_script, pdb_base, respos, target_mut, chain,
                           workdir, modeller_python)
        canonical = workdir / f"{pdb_base}_{wt}{respos}{target_mut}.pdb"
        shutil.move(str(out), str(canonical))
        mutants.append({
            "residue": f"{wt}{respos}",
            "resid": respos,
            "wt": wt,
            "mut": target_mut,
            "chain": chain,
            "pdb_file": canonical.name,
            "crypticity_score": r.get("score"),
            "esm_score": r.get("esm_score")
        })
    return mutants


def mutate_combined(modeller_script, pdb_base, residues, chain, workdir, modeller_python, n, suffix="esm"):
    """Applies ALL mutations sequentially into ONE combined structure."""
    current_base = pdb_base
    applied = []
    for i, r in enumerate(residues, 1):
        respos, wt, target_mut = r["resid"], r["resname_3letter"], r["target_mut"]
        out = run_modeller(modeller_script, current_base, respos, target_mut, chain,
                           workdir, modeller_python)
        rolling = workdir / f"{pdb_base}_step{i}.pdb"
        shutil.move(str(out), str(rolling))
        current_base = rolling.stem
        applied.append({
            "residue": f"{wt}{respos}",
            "resid": respos,
            "wt": wt,
            "mut": target_mut,
            "esm_score": r.get("esm_score")
        })

    final = workdir / f"{pdb_base}_top{n}_{suffix}_mutant.pdb"
    shutil.move(str(workdir / f"{current_base}.pdb"), str(final))

    # Clean intermediate step files
    for i in range(1, len(applied)):
        p = workdir / f"{pdb_base}_step{i}.pdb"
        if p.exists():
            p.unlink()

    log(f"         combined mutant with {len(applied)} mutations -> {final.name}")
    return [{
        "residue": "combined",
        "chain": chain,
        "mutations": applied,
        "pdb_file": final.name
    }]


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(
        description="CryptoBank -> ESM-Scan -> MODELLER pipeline orchestrator",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pdb", required=True, help="Target PDB code (e.g. 1JWP)")
    p.add_argument("--chain", default="A", help="Target chain ID (default: A)")
    p.add_argument("--top", type=int, default=5, help="Number of top cryptic residues to mutate")
    p.add_argument("--score_type", choices=["raw", "normalized"], default="raw")
    
    # Mutation strategy: ESM-Scan (default) or manual fixed amino acid
    p.add_argument("--strategy", choices=["esm", "manual"], default="esm",
                   help="Mutation strategy: 'esm' (scores and picks best mutation via ESM-1v) "
                        "or 'manual' (mutates all residues to --mutate_to)")
    p.add_argument("--mutate_to", default="GLU",
                   help="Target residue type (3-letter) if strategy is 'manual'.")
    p.add_argument("--mutation_set", default="charged_polar",
                   help="Which amino acids ESM-Scan may choose from: charged, polar, "
                        "charged_polar (default), all, or a custom comma list. "
                        "Charged/polar substitutions destabilise the closed state so "
                        "the cryptic pocket can open; 'all' restores the unrestricted "
                        "19-way scan, which tends to pick conservative nonpolar swaps.")
    
    p.add_argument("--mode", choices=["independent", "combined"], default="combined",
                   help="combined = ONE file with all mutations (default); "
                        "independent = N separate mutants")
    p.add_argument("--no_clean", action="store_true",
                   help="Feed the raw downloaded PDB to MODELLER (skip auto-cleaning)")
    p.add_argument("--workdir", default="pipeline_out")
    
    # Scripts & Models
    p.add_argument("--scraper", default="scrape_cryptobank.py")
    p.add_argument("--esm_script", default="esm_scanner.py")
    p.add_argument("--modeller_script", default="modeller_mutate.py")
    p.add_argument("--prep_script", default="clean_structure.py",
                   help="MDAnalysis structure-preparation script")
    p.add_argument("--esm_model",
                   default=os.environ.get("ESM_MODEL", "/opt/models/esm1v_t33_650M_UR90S_1.pt"),
                   help="Path to pre-trained ESM-1v model checkpoint")
    
    # Python Environments
    p.add_argument("--scraper_python", default=sys.executable,
                   help="Python interpreter for CryptoBank scraper")
    p.add_argument("--esm_python",
                   default="/opt/conda/bin/python3",
                   help="Python interpreter with torch & esm")
    p.add_argument("--modeller_python", default="/opt/conda/bin/python3",
                   help="Python interpreter with MODELLER")
    p.add_argument("--prep_python", default="/opt/conda/envs/prep/bin/python3",
                   help="Python interpreter with MDAnalysis (separate conda env)")
    
    p.add_argument("--check_numbering", action="store_true")
    p.add_argument("--skip_modeller", action="store_true")
    args = p.parse_args()

    # Resolve paths
    scraper_path = Path(args.scraper).resolve()
    esm_script_path = Path(args.esm_script).resolve()
    modeller_path = Path(args.modeller_script).resolve()
    prep_path = Path(args.prep_script).resolve()
    workdir = Path(args.workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    # Stage 1: fetch from RCSB and split protein / ligands
    protein_all = fetch_and_split(args.pdb, workdir,
                                  do_clean=not args.no_clean,
                                  prep_script=prep_path,
                                  prep_python=args.prep_python)

    # Stage 2: chain selection. Everything downstream -- CryptoBank, ESM-Scan,
    # MODELLER and SAM2 -- works on this single chain.
    pdb_file = select_chain(protein_all, args.chain, args.pdb, workdir,
                            do_clean=not args.no_clean,
                            prep_script=prep_path,
                            prep_python=args.prep_python)
    pdb_base = pdb_file.stem

    # Stage 3: CryptoBank / local ProtT5 on the selected chain
    residues = get_cryptic_residues(args.pdb, args.chain, args.top, args.score_type,
                                    workdir, scraper_path, args.scraper_python)
    if not residues:
        log("No cryptic residues returned; stopping.")
        sys.exit(1)

    if args.check_numbering:
        check_numbering(pdb_file, args.chain, residues)

    # Stage 4: Decide mutations (ESM-Scan or Manual)
    if args.strategy == "esm":
        if not os.path.exists(args.esm_model):
            log(f"WARNING: ESM model not found at {args.esm_model}. Falling back to manual '{args.mutate_to}'")
            for r in residues:
                r["target_mut"] = args.mutate_to.upper()
                r["esm_score"] = None
            suffix = args.mutate_to.lower()
        else:
            residues = run_esm_scan(pdb_file, args.chain, residues, workdir,
                                    esm_script_path, args.esm_model, args.esm_python,
                                    mutation_set=args.mutation_set)
            suffix = "esm"
    else:
        mutate_to = args.mutate_to.upper()
        if mutate_to not in VALID_AA3:
            log(f"ERROR: --mutate_to '{mutate_to}' is not a valid amino acid.")
            sys.exit(1)
        for r in residues:
            r["target_mut"] = mutate_to
            r["esm_score"] = None
        suffix = mutate_to.lower()

    # Residues whose wild type is already charged/polar are scored but not
    # mutated -- swapping one polar residue for another will not open a pocket.
    to_mutate = [r for r in residues if r.get("target_mut")]
    n_skipped = len(residues) - len(to_mutate)
    if n_skipped:
        log(f"Stage 4: {len(to_mutate)} of {len(residues)} residues will be mutated; "
            f"{n_skipped} skipped:")
        for r in residues:
            if not r.get("target_mut"):
                log(f"         {r['resname_3letter']}{r['resid']} -- "
                    f"{r.get('skip_reason', 'no mutation selected')}")
    if not to_mutate:
        log("ERROR: every cryptic residue was skipped, so there is nothing to mutate.")
        log("       All of them are already charged/polar. Consider a different "
            "--mutation_set, or raise --top to reach more residues.")
        sys.exit(1)

    # Stage 5: MODELLER mutagenesis
    mutants = []
    if not args.skip_modeller:
        log(f"Stage 5: MODELLER ({args.mode} mode)")
        if args.mode == "combined":
            mutants = mutate_combined(modeller_path, pdb_base, to_mutate,
                                      args.chain, workdir, args.modeller_python,
                                      len(to_mutate), suffix=suffix)
        else:
            mutants = mutate_independent(modeller_path, pdb_base, to_mutate,
                                         args.chain, workdir, args.modeller_python)
    else:
        log("Stage 5: skipped (--skip_modeller)")

    # Stage 5b: validate the mutant before the ensemble stage.
    # In "combined" mode there is one entry; in "independent" mode the last
    # mutant is used, matching the launcher's previous "first *_mutant.pdb" pick.
    sam2_input = None
    if mutants:
        final_mutant = workdir / mutants[-1]["pdb_file"]
        if final_mutant.exists():
            sam2_input = validate_for_sam2(final_mutant, args.chain, workdir,
                                            prep_path, args.prep_python)
        else:
            log(f"Stage 5b: mutant {final_mutant.name} not found, skipping SAM2 prep")

    # Write final manifest
    manifest = {
        "pdb": args.pdb,
        "chain": args.chain,
        "mode": args.mode,
        "strategy": args.strategy,
        "mutation_set": args.mutation_set if args.strategy == "esm" else None,
        "score_type": args.score_type,
        "cleaned": not args.no_clean,
        "protein_all_chains": protein_all.name,
        "source_structure": pdb_file.name,
        "n_residues_scored": len(residues),
        "n_residues_mutated": len(to_mutate),
        "n_residues_skipped": n_skipped,
        # Every scored residue is kept, skipped ones included, so the report can
        # show their charged/polar scores.
        "cryptic_residues": residues,
        "mutants": mutants,
        "sam2_input": sam2_input.name if sam2_input else None,
        "next_stage": "ensemble generation (aSAM / BioEmu)",
    }
    manifest_path = workdir / f"{args.pdb}_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2))
    log(f"Wrote manifest -> {manifest_path}")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
