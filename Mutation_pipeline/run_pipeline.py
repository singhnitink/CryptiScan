#!/usr/bin/env python3
"""
CryptoBank -> MODELLER Pipeline Orchestrator
=============================================
    scrape_cryptobank.py  (top cryptic residues)
            |
            v
    [fetch PDB from RCSB + auto-clean: keep chain, drop waters/hetero/altloc]
            |
            v
    modeller_mutate.py    (build mutant structure)
            |
            v
    manifest.json         (ready for the aSAM / BioEmu ensemble stage)

You do NOT activate any environment. Run this with no env active and point it
at each environment's python:

    python3 run_pipeline.py --pdb 1JWP --chain A --top 5 --mutate_to GLU \
        --mode combined \
        --scraper_python  /full/path/to/cryptobank_env/bin/python3 \
        --modeller_python python3 \
        --check_numbering

--mode combined   -> ONE output file with all top-N mutations applied  (what you want)
--mode independent-> N separate single-point mutants
--no_clean        -> feed the raw downloaded PDB to MODELLER (skip auto-cleaning)
"""

import argparse
import json
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

VALID_AA3 = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
}


def log(msg):
    print(f"[PIPELINE] {msg}", file=sys.stderr)


# ---------------------------------------------------------------------------
# STAGE 1: scrape CryptoBank
# ---------------------------------------------------------------------------
def get_cryptic_residues(pdb, chain, top, score_type, workdir, scraper, python_exe):
    out_json = workdir / f"{pdb}_{chain}_cryptic.json"
    cmd = [python_exe, str(scraper), "--pdb", pdb, "--chain", chain,
           "--top", str(top), "--score_type", score_type, "--output", str(out_json)]
    log(f"Stage 1: scraping CryptoBank -> {out_json.name}")
    log(f"         scraper python: {python_exe}")
    subprocess.run(cmd, check=True)
    data = json.loads(out_json.read_text())
    residues = data.get("top_residues", [])
    log(f"         got {len(residues)} cryptic residues: "
        + ", ".join(f"{r['resname_3letter']}{r['resid']}" for r in residues))
    return residues


# ---------------------------------------------------------------------------
# STAGE 2: fetch + (optionally) clean the PDB
# ---------------------------------------------------------------------------
def clean_pdb_text(raw_text, chain):
    """Keep only the target chain's ATOM records; drop HETATM/waters; keep the
    first alt-loc (blank or 'A') and blank the alt-loc indicator. Returns text.
    NOTE: modified residues stored as HETATM (e.g. MSE selenomethionine) are
    dropped — if your structure has those, use --no_clean or convert them first."""
    out = []
    for line in raw_text.splitlines():
        if line.startswith("ATOM"):
            if len(line) < 22 or line[21] != chain:
                continue
            altloc = line[16]
            if altloc not in (" ", "A"):
                continue
            out.append(line[:16] + " " + line[17:])  # blank the alt-loc
    out.append("TER")
    out.append("END")
    return "\n".join(out) + "\n"


def fetch_and_prepare_pdb(pdb, chain, workdir, do_clean):
    """Produce workdir/<pdb>.pdb ready for MODELLER.
    Priority: (1) a file you already placed there is used as-is;
              (2) otherwise download raw from RCSB and clean (unless --no_clean)."""
    dest = workdir / f"{pdb.lower()}.pdb"
    if dest.exists():
        log(f"Stage 2: using existing {dest.name} (not overwriting)")
        return dest

    raw = workdir / f"{pdb.lower()}_raw.pdb"
    url = f"https://files.rcsb.org/download/{pdb.upper()}.pdb"
    log(f"Stage 2: fetching {url}")
    try:
        urllib.request.urlretrieve(url, raw)
    except Exception as e:
        log(f"  ERROR downloading PDB: {e}")
        log(f"  Place the file manually at: {dest}")
        sys.exit(1)

    if do_clean:
        cleaned = clean_pdb_text(raw.read_text(), chain)
        n_atoms = sum(1 for l in cleaned.splitlines() if l.startswith("ATOM"))
        if n_atoms == 0:
            log(f"  ERROR: cleaning left 0 atoms for chain {chain}. "
                f"Is '{chain}' the right chain? Use --no_clean to inspect the raw file.")
            sys.exit(1)
        dest.write_text(cleaned)
        log(f"Stage 2: cleaned -> {dest.name} (chain {chain}, {n_atoms} atoms, "
            f"waters/hetero/alt-locs removed)")
    else:
        shutil.move(str(raw), str(dest))
        log(f"Stage 2: using raw structure -> {dest.name} (--no_clean)")
    return dest


def check_numbering(pdb_file, chain, residues):
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
    log(f"Stage 2: numbering check OK ({len(residues)} residues matched in chain {chain})")
    return True


# ---------------------------------------------------------------------------
# STAGE 3: MODELLER mutations
# ---------------------------------------------------------------------------
def run_modeller(modeller_script, base, respos, restyp, chain, workdir, modeller_python):
    """modeller_mutate.py args: modelname respos restyp chain
    reads <base>.pdb, writes <base><restyp><respos>.pdb in workdir."""
    cmd = [modeller_python, str(modeller_script), base, str(respos), restyp, chain]
    logf = workdir / f"modeller_{restyp}{respos}.log"
    log(f"  mutate {chain}/{respos} -> {restyp}   (log: {logf.name})")
    with open(logf, "w") as lf:
        subprocess.run(cmd, cwd=str(workdir), check=True, stdout=lf, stderr=subprocess.STDOUT)
    out = workdir / f"{base}{restyp}{respos}.pdb"
    if not out.exists():
        raise FileNotFoundError(f"MODELLER did not produce {out} (see {logf})")
    return out


def mutate_independent(modeller_script, pdb_base, residues, mutate_to, chain,
                       workdir, modeller_python):
    mutants = []
    for r in residues:
        respos, wt = r["resid"], r["resname_3letter"]
        out = run_modeller(modeller_script, pdb_base, respos, mutate_to, chain,
                           workdir, modeller_python)
        canonical = workdir / f"{pdb_base}_{wt}{respos}{mutate_to}.pdb"
        shutil.move(str(out), str(canonical))
        mutants.append({"residue": f"{wt}{respos}", "resid": respos, "wt": wt,
                        "mut": mutate_to, "chain": chain,
                        "pdb_file": canonical.name, "crypticity_score": r.get("score")})
    return mutants


def mutate_combined(modeller_script, pdb_base, residues, mutate_to, chain,
                    workdir, modeller_python, n):
    """Apply ALL mutations sequentially into ONE structure."""
    current_base = pdb_base
    applied = []
    for i, r in enumerate(residues, 1):
        respos, wt = r["resid"], r["resname_3letter"]
        out = run_modeller(modeller_script, current_base, respos, mutate_to, chain,
                           workdir, modeller_python)
        rolling = workdir / f"{pdb_base}_step{i}.pdb"
        shutil.move(str(out), str(rolling))
        current_base = rolling.stem
        applied.append({"residue": f"{wt}{respos}", "resid": respos,
                        "wt": wt, "mut": mutate_to})
    final = workdir / f"{pdb_base}_top{n}_{mutate_to}_mutant.pdb"
    shutil.move(str(workdir / f"{current_base}.pdb"), str(final))
    # tidy intermediate step files
    for i in range(1, len(applied)):
        p = workdir / f"{pdb_base}_step{i}.pdb"
        if p.exists():
            p.unlink()
    log(f"         combined mutant with {len(applied)} mutations -> {final.name}")
    return [{"residue": "combined", "chain": chain, "mutations": applied,
             "mut": mutate_to, "pdb_file": final.name}]


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(
        description="CryptoBank -> MODELLER pipeline orchestrator",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pdb", required=True)
    p.add_argument("--chain", default="A")
    p.add_argument("--top", type=int, default=5)
    p.add_argument("--score_type", choices=["raw", "normalized"], default="raw")
    p.add_argument("--mutate_to", default="ALA",
                   help="Target residue type (3-letter). Set to your method's rule.")
    p.add_argument("--mode", choices=["independent", "combined"], default="combined",
                   help="combined = ONE file with all mutations (default); "
                        "independent = N separate mutants")
    p.add_argument("--no_clean", action="store_true",
                   help="Feed the raw downloaded PDB to MODELLER (skip auto-cleaning)")
    p.add_argument("--workdir", default="pipeline_out")
    p.add_argument("--scraper", default="scrape_cryptobank.py")
    p.add_argument("--modeller_script", default="modeller_mutate.py")
    p.add_argument("--scraper_python", default=sys.executable,
                   help="Full path to the python that has gradio_client (cryptobank_env)")
    p.add_argument("--modeller_python", default="python3",
                   help="Python that has MODELLER (usually your base env)")
    p.add_argument("--check_numbering", action="store_true")
    p.add_argument("--skip_modeller", action="store_true")
    args = p.parse_args()

    mutate_to = args.mutate_to.upper()
    if mutate_to not in VALID_AA3:
        log(f"ERROR: --mutate_to '{mutate_to}' is not a valid amino acid.")
        sys.exit(1)

    # MODELLER runs with cwd=workdir, so script paths must be absolute
    scraper_path = Path(args.scraper).resolve()
    modeller_path = Path(args.modeller_script).resolve()

    workdir = Path(args.workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    residues = get_cryptic_residues(args.pdb, args.chain, args.top, args.score_type,
                                    workdir, scraper_path, args.scraper_python)
    if not residues:
        log("No cryptic residues returned; stopping.")
        sys.exit(1)

    pdb_file = fetch_and_prepare_pdb(args.pdb, args.chain, workdir, do_clean=not args.no_clean)
    pdb_base = pdb_file.stem

    if args.check_numbering:
        check_numbering(pdb_file, args.chain, residues)

    mutants = []
    if not args.skip_modeller:
        log(f"Stage 3: MODELLER ({args.mode} mode, all -> {mutate_to})")
        if args.mode == "combined":
            mutants = mutate_combined(modeller_path, pdb_base, residues,
                                      mutate_to, args.chain, workdir,
                                      args.modeller_python, args.top)
        else:
            mutants = mutate_independent(modeller_path, pdb_base, residues,
                                         mutate_to, args.chain, workdir, args.modeller_python)
    else:
        log("Stage 3: skipped (--skip_modeller)")

    manifest = {
        "pdb": args.pdb, "chain": args.chain, "mode": args.mode,
        "mutate_to": mutate_to, "score_type": args.score_type,
        "cleaned": not args.no_clean,
        "source_structure": pdb_file.name,
        "cryptic_residues": residues, "mutants": mutants,
        "next_stage": "ensemble generation (aSAMt / BioEmu)",
    }
    manifest_path = workdir / f"{args.pdb}_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2))
    log(f"Wrote manifest -> {manifest_path}")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
