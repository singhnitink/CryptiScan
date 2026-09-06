#!/usr/bin/env python3
"""
ESM-Scan Mutation Selector for CryptiScan
=========================================
This script uses Facebook Research's ESM-1v protein language model
to evaluate and score all possible amino acid mutations at target
residue positions (e.g. cryptic residues identified by CryptoBank).

How ESM Scoring Works:
----------------------
1. The wild-type (WT) sequence is passed through the ESM-1v model.
2. For each position, the model computes log-probabilities for all 20 amino acids.
3. For any mutation (WT -> Mutant):
       score = log_prob(Mutant) - log_prob(WT)
   - A score near 0 or positive means the mutation is favored / tolerated.
   - A negative score means the mutation is disfavored / deleterious.
4. For each target residue, all 19 alternative amino acids are scored.
   The mutation with the HIGHEST (best) score is selected.

Usage (Standalone):
-------------------
    python3 esm_scanner.py \
        --pdb 1JWP.pdb \
        --chain A \
        --residues 26,45,102 \
        --model_path /opt/models/esm1v_t33_650M_UR90S_1.pt \
        --output esm_results.json
"""

import argparse
import csv
import json
import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple

# Standard 20 amino acids
THREE_TO_ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}
ONE_TO_THREE = {v: k for k, v in THREE_TO_ONE.items()}
ALL_AMINO_ACIDS = sorted(list(ONE_TO_THREE.keys()))  # 20 standard 1-letter codes


def log(msg: str):
    """Prints formatted progress messages to standard error."""
    print(f"[ESM-SCAN] {msg}", file=sys.stderr)


def extract_chain_sequence_from_pdb(pdb_path: str, chain_id: str = "A") -> Tuple[str, Dict[int, int]]:
    """
    Extracts the sequence of a single chain from a PDB file.

    Parameters:
        pdb_path: Path to the .pdb file.
        chain_id: Target chain identifier (default: "A").

    Returns:
        sequence (str): 1-letter amino acid sequence of the chain.
        resid_to_seqidx (dict): Mapping from PDB residue number (resid)
                                to 0-based index in the sequence string.
    """
    seen_resids = set()
    ordered_residues = []

    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith("ATOM") and len(line) >= 26:
                c = line[21]
                if c != chain_id:
                    continue
                try:
                    resid = int(line[22:26])
                except ValueError:
                    continue

                # Record each residue only once (first atom encountered)
                if resid not in seen_resids:
                    seen_resids.add(resid)
                    resname_3 = line[17:20].strip()
                    resname_1 = THREE_TO_ONE.get(resname_3, "X")
                    ordered_residues.append((resid, resname_1))

    sequence = "".join(r[1] for r in ordered_residues)
    resid_to_seqidx = {r[0]: idx for idx, r in enumerate(ordered_residues)}
    return sequence, resid_to_seqidx


def score_residue_mutations(
    sequence: str,
    target_resids: List[int],
    resid_to_seqidx: Dict[int, int],
    model_path: str,
    device: str = None,
) -> List[Dict]:
    """
    Runs ESM-1v on the sequence and calculates mutation scores for each target residue.

    Parameters:
        sequence: Wild-type protein sequence (1-letter code).
        target_resids: List of PDB residue numbers to evaluate.
        resid_to_seqidx: Mapping from PDB resid to 0-based index in sequence.
        model_path: Path to pre-trained ESM-1v weights (.pt file).
        device: 'cuda' or 'cpu' (auto-detects GPU if available).

    Returns:
        results: List of dictionaries containing the best mutation and
                 ranked candidates for each target residue.
    """
    try:
        import torch
        import esm
    except ImportError as e:
        log(f"ERROR: Missing PyTorch or fair-esm library ({e}).")
        log("This script is meant to run inside the CryptiScan container.")
        sys.exit(1)

    # 1. Device selection (GPU if available, else CPU)
    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"
    log(f"Using device: {device}")

    # 2. Load the ESM-1v model and alphabet
    log(f"Loading ESM-1v model from: {model_path}")
    if not os.path.exists(model_path):
        raise FileNotFoundError(f"ESM model file not found at: {model_path}")

    model, alphabet = esm.pretrained.load_model_and_alphabet(model_path)
    model.eval()
    model = model.to(device)
    log("ESM-1v model loaded successfully.")

    # 3. Tokenize sequence and run a single forward pass
    batch_converter = alphabet.get_batch_converter()
    data = [("target_protein", sequence)]
    _, _, batch_tokens = batch_converter(data)
    batch_tokens = batch_tokens.to(device)

    log("Computing sequence log-probabilities with ESM-1v...")
    with torch.no_grad():
        token_logits = model(batch_tokens)["logits"]
        token_probs = torch.log_softmax(token_logits, dim=-1).cpu()

    results = []

    # 4. Score each target residue
    for resid in target_resids:
        if resid not in resid_to_seqidx:
            log(f"Warning: Residue {resid} not found in chain sequence. Skipping.")
            continue

        seq_idx = resid_to_seqidx[resid]
        wt_aa1 = sequence[seq_idx]
        wt_aa3 = ONE_TO_THREE.get(wt_aa1, "UNK")

        if wt_aa1 not in ALL_AMINO_ACIDS:
            log(f"Warning: Non-standard residue {wt_aa1} at resid {resid}. Skipping.")
            continue

        wt_tok = alphabet.get_idx(wt_aa1)
        # Note: In ESM tokens, token 0 is BOS (beginning-of-sequence),
        # so sequence index 0 corresponds to token position 1.
        wt_prob = token_probs[0, 1 + seq_idx, wt_tok].item()

        # Score all 19 alternative amino acids
        candidates = []
        for mt_aa1 in ALL_AMINO_ACIDS:
            if mt_aa1 == wt_aa1:
                continue  # Skip wild-type

            mt_tok = alphabet.get_idx(mt_aa1)
            mt_prob = token_probs[0, 1 + seq_idx, mt_tok].item()

            # Zero-shot score: difference in log-probability (mutant - WT)
            esm_score = mt_prob - wt_prob

            candidates.append({
                "mut_aa1": mt_aa1,
                "mut_aa3": ONE_TO_THREE[mt_aa1],
                "esm_score": round(float(esm_score), 4),
            })

        # Sort candidates descending: highest score = most favored mutation
        candidates.sort(key=lambda x: x["esm_score"], reverse=True)
        best = candidates[0]

        log(f"Residue {wt_aa3}{resid} -> Best ESM mutation: {best['mut_aa3']} (score: {best['esm_score']:+.4f})")

        results.append({
            "resid": resid,
            "wt_aa1": wt_aa1,
            "wt_aa3": wt_aa3,
            "best_mut_aa1": best["mut_aa1"],
            "best_mut_aa3": best["mut_aa3"],
            "best_score": best["esm_score"],
            "candidates_ranked": candidates,
        })

    return results


def save_reports(results: List[Dict], output_json: str, output_csv: str):
    """Saves formatted JSON and CSV summaries of ESM predictions."""
    # 1. Save detailed JSON
    with open(output_json, "w") as f:
        json.dump(results, f, indent=2)
    log(f"Saved JSON report -> {output_json}")

    # 2. Save human-readable CSV summary
    with open(output_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["resid", "wt", "best_mutant", "best_score", "top3_options"])
        for r in results:
            top3 = "; ".join(f"{c['mut_aa3']}({c['esm_score']:+.2f})" for c in r["candidates_ranked"][:3])
            writer.writerow([r["resid"], r["wt_aa3"], r["best_mut_aa3"], r["best_score"], top3])
    log(f"Saved CSV report  -> {output_csv}")


def main():
    parser = argparse.ArgumentParser(
        description="ESM-Scan: Zero-shot variant scoring for CryptiScan cryptic residues.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--pdb", required=True, help="Path to PDB file (e.g. 1JWP.pdb)")
    parser.add_argument("--chain", default="A", help="Target chain ID")
    parser.add_argument("--residues", required=True, help="Comma-separated PDB residue numbers (e.g. 26,45,102)")
    parser.add_argument("--model_path",
                        default=os.environ.get("ESM_MODEL", "/opt/models/esm1v_t33_650M_UR90S_1.pt"),
                        help="Path to pre-trained ESM-1v model checkpoint")
    parser.add_argument("--device", default=None, choices=["cuda", "cpu"], help="Compute device (default: auto)")
    parser.add_argument("--output_json", default="esm_scan_results.json", help="Output JSON path")
    parser.add_argument("--output_csv", default="esm_scan_summary.csv", help="Output CSV path")

    args = parser.parse_args()

    # Parse target residue numbers
    target_resids = [int(r.strip()) for r in args.residues.split(",") if r.strip()]

    # Extract sequence and residue mapping from PDB
    sequence, resid_to_seqidx = extract_chain_sequence_from_pdb(args.pdb, args.chain)
    log(f"Extracted chain {args.chain} sequence: {len(sequence)} amino acids.")

    # Run ESM scoring
    results = score_residue_mutations(
        sequence=sequence,
        target_resids=target_resids,
        resid_to_seqidx=resid_to_seqidx,
        model_path=args.model_path,
        device=args.device,
    )

    # Save reports
    save_reports(results, args.output_json, args.output_csv)
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
