#!/usr/bin/env python3
"""
CryptoBank Terminal Client
===========================
Get per-residue cryptic pocket predictions via the CryptoBank Gradio API.

Usage:
    python3 scrape_cryptobank.py --pdb 1JWP
    python3 scrape_cryptobank.py --pdb 1JWP --chain A --top 10
    python3 scrape_cryptobank.py --pdb 1JWP --score_type normalized
    python3 scrape_cryptobank.py --pdb 1JWP --output results.json

Setup:  pip install gradio_client
"""

import argparse
import json
import sys
import time
import re

GRADIO_URL = "https://thorbenf-cryptobank.hf.space"

THREE_TO_ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}


def predict_cryptobank(pdb_id, chain="A", score_type="Raw Scores"):
    from gradio_client import Client

    print(f"[1/3] Connecting to CryptoBank ...", file=sys.stderr)
    client = Client(GRADIO_URL)

    print(f"[2/3] Fetching PDB {pdb_id} ...", file=sys.stderr)
    client.predict("PDB ID", pdb_id, None, api_name="/fetch_interface")

    print(f"[3/3] Running prediction (chain {chain}, {score_type}) ...", file=sys.stderr)
    result = client.predict(
        "PDB ID", pdb_id, None, chain, score_type,
        api_name="/process_interface"
    )

    # result is a tuple:
    #   [0] = PyMOL commands text
    #   [1] = HTML viewer
    #   [2] = list of file paths [residues.txt, scores.pdb]
    # Read the residues txt file from result[2]
    prediction_text = ""
    if isinstance(result, tuple) and len(result) > 2:
        file_list = result[2]
        if isinstance(file_list, list):
            for fpath in file_list:
                if 'binding_site_residues' in str(fpath):
                    print(f"  Reading: {fpath}", file=sys.stderr)
                    with open(fpath, 'r') as f:
                        prediction_text = f.read()
                    break

    return prediction_text


def parse_prediction_text(text, chain):
    residues = []
    for line in text.split('\n'):
        line = line.strip()
        m = re.match(r'^([A-Z]{3})\s+(\d+)\s+([A-Z])\s+([\d.]+)$', line)
        if m:
            resname_3, resnum, resname_1, score = m.groups()
            if resname_3 in THREE_TO_ONE:
                residues.append({
                    "resid": int(resnum),
                    "resname_3letter": resname_3,
                    "resname_1letter": resname_1,
                    "chain": chain,
                    "score": round(float(score), 4),
                })
    return residues


def print_results(pdb_id, chain, residues, top_n, score_type, fmt="text"):
    residues_sorted = sorted(residues, key=lambda x: x["score"], reverse=True)
    top = residues_sorted[:top_n]

    output = {
        "query": {"pdb_id": pdb_id, "chain": chain, "score_type": score_type},
        "source": "CryptoBank PLM (thorbenf-cryptobank.hf.space)",
        "top_residues": top,
        "all_residues": residues_sorted,
        "n_total_residues": len(residues),
    }

    if fmt == "json":
        return json.dumps(output, indent=2)

    lines = []
    lines.append("")
    lines.append("=" * 65)
    lines.append("  CryptoBank PLM — Top Cryptic Residues")
    lines.append(f"  PDB: {pdb_id}  Chain: {chain}  Scores: {score_type}")
    lines.append("=" * 65)

    if not top:
        lines.append("")
        lines.append("  No residue-level predictions obtained.")
        lines.append("")
    else:
        lines.append("")
        lines.append(f"  Top {len(top)} residues by crypticity score")
        lines.append(f"  (out of {len(residues)} total residues)")
        lines.append("")
        lines.append(f"  {'Rank':<6}{'Residue':<14}{'ResNum':<8}{'Chain':<7}{'Score':<10}")
        lines.append(f"  {'─'*6}{'─'*14}{'─'*8}{'─'*7}{'─'*10}")

        for i, r in enumerate(top, 1):
            label = f"{r['resname_3letter']} ({r['resname_1letter']})"
            lines.append(
                f"  {i:<6}{label:<14}{r['resid']:<8}{r['chain']:<7}{r['score']:<10.4f}"
            )

        lines.append("")
        res_list = ", ".join(f"{r['resname_3letter']}{r['resid']}" for r in top)
        lines.append(f"  Top residues: {res_list}")
        lines.append("")

    lines.append("=" * 65)
    lines.append("")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="CryptoBank Terminal Client — Per-residue cryptic pocket prediction",
        epilog="Setup:  pip install gradio_client"
    )
    parser.add_argument("--pdb", type=str, required=True, help="PDB ID")
    parser.add_argument("--chain", type=str, default="A", help="Chain (default: A)")
    parser.add_argument("--top", type=int, default=5, help="Top N residues (default: 5)")
    parser.add_argument("--score_type", choices=["raw", "normalized"], default="raw")
    parser.add_argument("--output", type=str, help="Save JSON results to file")
    parser.add_argument("--format", choices=["text", "json"], default="text")

    args = parser.parse_args()
    score_type_val = "Raw Scores" if args.score_type == "raw" else "Normalized Scores"

    t0 = time.time()
    prediction_text = predict_cryptobank(args.pdb, args.chain, score_type_val)
    print(f"[INFO] Done in {time.time()-t0:.0f}s", file=sys.stderr)

    residues = parse_prediction_text(prediction_text, args.chain)
    print(f"[INFO] Parsed {len(residues)} residues", file=sys.stderr)

    fmt = "json" if args.output or args.format == "json" else "text"
    output = print_results(args.pdb, args.chain, residues, args.top, args.score_type, fmt)

    if args.output:
        with open(args.output, 'w') as f:
            f.write(print_results(args.pdb, args.chain, residues, args.top, args.score_type, "json"))
        print(f"[INFO] Saved to: {args.output}", file=sys.stderr)

    print(output)


if __name__ == "__main__":
    main()
