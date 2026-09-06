#!/usr/bin/env python3
"""
Builds the human-readable run report for a CryptiScan job.

Reads the structured outputs the pipeline already writes -- <PDB>_manifest.json,
<PDB>_<chain>_cryptic.json, <pdb>_esm_results.json and <pdb>_prep.json -- rather
than scraping stderr, so nothing depends on log formatting.

Also writes the wild-type and mutant sequences as FASTA.

Deliberately dependency-free (standard library only) so it runs in either conda
env, or on the host. mdtraj is used only if present, to count ensemble frames.

Usage:
    python3 make_report.py --workdir <run dir> [--jobname NAME] [--container FILE]
"""

import argparse
import datetime as _dt
import json
import re
import sys
from pathlib import Path

THREE_TO_ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}

W = 80
def rule(c="-"):
    return c * W

def head(title):
    return f"{rule()}\n {title}\n{rule()}"


def read_json(path):
    try:
        return json.loads(Path(path).read_text())
    except Exception:
        return None


def chain_sequence(pdb_path, chain=None):
    """(sequence, [resid...]) from CA atoms, in file order. Pure text parsing."""
    seq, ids = [], []
    try:
        lines = Path(pdb_path).read_text().splitlines()
    except Exception:
        return "", []
    for line in lines:
        if not line.startswith("ATOM") or len(line) < 27:
            continue
        if line[12:16].strip() != "CA":
            continue
        if chain and line[21] != chain:
            continue
        resname = line[17:20].strip()
        seq.append(THREE_TO_ONE.get(resname, "X"))
        try:
            ids.append(int(line[22:26]))
        except ValueError:
            ids.append(None)
    return "".join(seq), ids


def wrap_fasta(seq, width=60):
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))


def ensemble_frames(workdir):
    dcd = workdir / "ensemble_output.traj.dcd"
    top = workdir / "ensemble_output.top.pdb"
    if not dcd.exists():
        return None
    try:
        import mdtraj
        return mdtraj.load(str(dcd), top=str(top)).n_frames
    except Exception:
        # DCD header: frame count is a 4-byte int at offset 8
        try:
            import struct
            with open(dcd, "rb") as fp:
                fp.seek(8)
                return struct.unpack("<i", fp.read(4))[0]
        except Exception:
            return None


def build(workdir, jobname, container):
    workdir = Path(workdir)
    out = []

    manifest_files = sorted(workdir.glob("*_manifest.json"))
    manifest = read_json(manifest_files[0]) if manifest_files else {}
    if not manifest:
        print(f"ERROR: no *_manifest.json in {workdir}", file=sys.stderr)
        sys.exit(1)

    pdb_id = manifest.get("pdb", "?")
    chain = manifest.get("chain", "?")
    base = pdb_id.lower()

    cryptic = read_json(workdir / f"{pdb_id}_{chain}_cryptic.json") or {}
    esm = read_json(workdir / f"{base}_esm_results.json") or []
    prep = read_json(workdir / f"{base}_prep.json") or {}

    # timestamp from the <jobname>_<YYYYmmdd>_<HHMMSS> directory name
    stamp = "unknown"
    m = re.search(r"_(\d{8})_(\d{6})$", workdir.name)
    if m:
        try:
            stamp = _dt.datetime.strptime(m.group(1) + m.group(2),
                                          "%Y%m%d%H%M%S").strftime("%Y-%m-%d %H:%M:%S")
        except ValueError:
            pass
    if not jobname:
        jobname = re.sub(r"_\d{8}_\d{6}$", "", workdir.name)

    engine = cryptic.get("source", "unknown")
    mode_label = "local (ProtT5, offline)" if "Local" in engine else \
                 "web (CryptoBank HF Space)" if "hf.space" in engine else "unknown"

    # ---------------------------------------------------------------- header
    out.append(rule("="))
    out.append(" CryptiScan - Run Report")
    out.append(rule("="))
    out.append(f" Job name       : {jobname}")
    out.append(f" Run directory  : {workdir.name}")
    out.append(f" Run started    : {stamp}")
    out.append(f" Report written : {_dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    out.append(f" Pipeline       : {mode_label}")
    if container:
        out.append(f" Container      : {container}")
    out.append("")

    # ------------------------------------------------------- target structure
    out.append(head("1. TARGET STRUCTURE"))
    out.append(f" PDB ID                : {pdb_id}")
    out.append(f" Selected chain        : {chain}")
    out.append(f" Source                : https://files.rcsb.org/download/{pdb_id.upper()}.pdb")
    if prep:
        out.append(f" Chains in entry       : {', '.join(prep.get('chains', [])) or 'n/a'}")
        out.append(f" Protein residues      : {prep.get('n_protein_residues', '?')} (all chains)")
        out.append(f" Waters dropped        : {prep.get('n_waters_dropped', '?')}")
    wt_seq, wt_ids = chain_sequence(workdir / manifest.get("source_structure", f"{base}.pdb"), chain)
    out.append(f" Selected chain length : {len(wt_seq)} residues")
    out.append(f" Cleaned               : {manifest.get('cleaned')}")
    out.append("")
    ligs = prep.get("ligands", []) if prep else []
    out.append(f" Ligands extracted     : {len(ligs)}")
    for lg in ligs:
        out.append(f"     {lg['resname']:<4} chain {lg.get('chain') or '?':<2} "
                   f"resid {lg['resid']:<5} {lg['n_atoms']:>3} atoms  -> {lg['file']}")
    out.append("")

    # ---------------------------------------------------- cryptic prediction
    out.append(head("2. CRYPTIC POCKET PREDICTION"))
    out.append(f" Engine       : {engine}")
    out.append(f" Score type   : {manifest.get('score_type', '?')}")
    tops = cryptic.get("top_residues", [])
    out.append(f" Residues     : top {len(tops)} of {cryptic.get('n_total_residues', '?')}")
    out.append("")
    out.append("   Rank  Residue      Crypticity")
    out.append("   ----  -----------  ----------")
    for i, r in enumerate(tops, 1):
        res = f"{r['resname_3letter']}{r['resid']}"
        sc = r.get("score")
        out.append(f"   {i:>4}  {res:<11}  {sc:.4f}" if isinstance(sc, (int, float))
                   else f"   {i:>4}  {res:<11}  {sc}")
    out.append("")

    # ------------------------------------------------------------- ESM-Scan
    out.append(head("3. ESM-SCAN MUTATION SELECTION"))
    out.append(" Score = logP(mutant) - logP(wild-type) from ESM-1v; higher is better.")
    out.append(" ESM chooses WHICH amino acid to place at each cryptic position;")
    out.append(" the positions themselves come from the cryptic prediction above.")
    out.append("")
    out.append("   Residue      Chosen    Score     Next best alternatives")
    out.append("   -----------  --------  --------  ------------------------------------")
    for e in esm:
        res = f"{e.get('wt_aa3','?')}{e.get('resid','?')}"
        alts = "; ".join(f"{c['mut_aa3']}({c['esm_score']:+.2f})"
                         for c in e.get("candidates_ranked", [])[1:4])
        out.append(f"   {res:<11}  {e.get('best_mut_aa3','?'):<8}  "
                   f"{e.get('best_score', 0):>+8.4f}  {alts}")
    out.append("")
    out.append(f" Full ranking of all 19 substitutions per site: {base}_esm_summary.csv")
    out.append("")

    # ----------------------------------------------------------- mutagenesis
    out.append(head("4. IN SILICO MUTAGENESIS (MODELLER)"))
    mode = manifest.get("mode", "?")
    out.append(f" Mode : {mode}")
    if mode == "combined":
        out.append("        all mutations applied sequentially into ONE structure;")
        out.append("        each mutation is relaxed in the presence of the previous ones.")
    else:
        out.append("        each mutation applied separately to the wild-type structure.")
    out.append("")
    mutants = manifest.get("mutants", [])
    applied = []
    for mu in mutants:
        if mu.get("residue") == "combined":
            applied = mu.get("mutations", [])
            out.append(" Applied in this order:")
            for i, mm in enumerate(applied, 1):
                out.append(f"   {i}. {mm['wt']}{mm['resid']} -> {mm['mut']}"
                           f"   (ESM {mm.get('esm_score', 0):+.4f})")
            out.append("")
            out.append(f" Output structure : {mu.get('pdb_file')}")
        else:
            out.append(f"   {mu.get('residue')} -> {mu.get('mut')}   ->  {mu.get('pdb_file')}")
            applied.append({"wt": mu.get("wt"), "resid": mu.get("resid"), "mut": mu.get("mut")})
    out.append("")

    # -------------------------------------------------------------- sequences
    out.append(head("5. SEQUENCES"))
    mut_pdb = None
    for mu in mutants:
        if mu.get("pdb_file"):
            mut_pdb = workdir / mu["pdb_file"]
    mut_seq, mut_ids = chain_sequence(mut_pdb, chain) if mut_pdb else ("", [])

    out.append(f" Wild-type (chain {chain}, {len(wt_seq)} aa)")
    out.append("")
    out.append(f" >{pdb_id}_{chain}|wild-type|{len(wt_seq)}aa")
    out.append(wrap_fasta(wt_seq))
    out.append("")
    if mut_seq:
        out.append(f" Mutant ({len(applied)} mutations, {len(mut_seq)} aa)")
        out.append("")
        out.append(f" >{pdb_id}_{chain}|mutant|{len(applied)}mut|{len(mut_seq)}aa")
        out.append(wrap_fasta(mut_seq))
        out.append("")
        if len(wt_seq) == len(mut_seq):
            diffs = [(i, a, b) for i, (a, b) in enumerate(zip(wt_seq, mut_seq)) if a != b]
            out.append(f" Differences vs wild-type : {len(diffs)}")
            for i, a, b in diffs:
                rid = wt_ids[i] if i < len(wt_ids) else "?"
                out.append(f"   seq position {i+1:<5} (PDB resid {rid})  {a} -> {b}")
            if len(diffs) != len(applied):
                out.append(f"   NOTE: {len(applied)} mutations were requested but "
                           f"{len(diffs)} sequence differences are present.")
        else:
            out.append(f" NOTE: length mismatch, wild-type {len(wt_seq)} vs mutant {len(mut_seq)};"
                       " sequence diff skipped.")
    else:
        out.append(" Mutant sequence unavailable (no mutant structure found).")
    out.append("")

    # --------------------------------------------------------------- ensemble
    out.append(head("6. ENSEMBLE GENERATION (aSAM / SAM2)"))
    n = ensemble_frames(workdir)
    if n is not None:
        out.append(f" Frames generated : {n}")
        out.append(f" Input structure  : {manifest.get('sam2_input')}")
        out.append(" Files            : ensemble_output.top.pdb, ensemble_output.traj.dcd")
        out.append("")
        out.append(" NOTE: SAM2 renumbers residues from 0 and drops the C-terminal OXT")
        out.append("       atom, so load the trajectory with ensemble_output.top.pdb as")
        out.append("       the topology, not the input PDB.")
    else:
        out.append(" No ensemble found (SAM2 stage did not run or produced no output).")
    out.append("")

    # ---------------------------------------------------------------- files
    out.append(head("7. FILES IN THIS ARCHIVE"))
    for f in sorted(workdir.rglob("*")):
        if f.is_file():
            out.append(f"   {str(f.relative_to(workdir)):<44} {f.stat().st_size:>12,} B")
    out.append("   (plus this report and the two .fasta files, written afterwards)")
    out.append("")
    out.append(rule("="))
    out.append(" End of report")
    out.append(rule("="))

    return "\n".join(out) + "\n", (pdb_id, chain, wt_seq, mut_seq, len(applied))


def main():
    p = argparse.ArgumentParser(description="Build the CryptiScan run report")
    p.add_argument("--workdir", required=True)
    p.add_argument("--jobname", default=None)
    p.add_argument("--container", default=None)
    p.add_argument("--output", default=None, help="Report path (default <workdir>/report.txt)")
    a = p.parse_args()

    workdir = Path(a.workdir)
    text, (pdb_id, chain, wt_seq, mut_seq, n_mut) = build(workdir, a.jobname, a.container)

    out = Path(a.output) if a.output else workdir / "report.txt"
    out.write_text(text)
    print(f"[REPORT] wrote {out}", file=sys.stderr)

    if wt_seq:
        f = workdir / f"{pdb_id}_{chain}_wildtype.fasta"
        f.write_text(f">{pdb_id}_{chain}|wild-type|{len(wt_seq)}aa\n{wrap_fasta(wt_seq)}\n")
        print(f"[REPORT] wrote {f}", file=sys.stderr)
    if mut_seq:
        f = workdir / f"{pdb_id}_{chain}_mutant.fasta"
        f.write_text(f">{pdb_id}_{chain}|mutant|{n_mut}mut|{len(mut_seq)}aa\n{wrap_fasta(mut_seq)}\n")
        print(f"[REPORT] wrote {f}", file=sys.stderr)


if __name__ == "__main__":
    main()
