#!/usr/bin/env python3
"""
CryptoBank Local Inference Predictor
====================================
Predicts cryptic pocket residues entirely offline using the local ProtT5
model weights (prot_t5_xl_uniref50_full_v2). Zero web/API dependency!

Usage:
    python3 predict_cryptic_local.py \
        --pdb 1JWP \
        --chain A \
        --top 5 \
        --model_dir /opt/models/prot_t5_xl_uniref50_full_v2 \
        --output 1jwp_A_cryptic.json
"""

import argparse
import json
import os
import sys
from pathlib import Path
import numpy as np
from scipy.special import expit

import torch
from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1

# Import the model loader from the local directory
try:
    from cryptobank_model_loader import load_model
except ImportError:
    from model_loader import load_model

THREE_TO_ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}


def log(msg: str):
    print(f"[CryptoBank-Local] {msg}", file=sys.stderr)


def normalize_scores(scores):
    min_score = np.min(scores)
    max_score = np.max(scores)
    if max_score > min_score:
        return (scores - min_score) / (max_score - min_score)
    return scores


def resolve_pdb_file(pdb_id_or_file: str, workdir: Path = None) -> str:
    """Finds existing PDB file locally or downloads it from RCSB."""
    p = Path(pdb_id_or_file)
    if p.is_file():
        return str(p.resolve())

    # Check workdir if given
    if workdir:
        candidate = workdir / f"{pdb_id_or_file.lower()}.pdb"
        if candidate.is_file():
            return str(candidate.resolve())
        candidate_raw = workdir / f"{pdb_id_or_file.lower()}_raw.pdb"
        if candidate_raw.is_file():
            return str(candidate_raw.resolve())

    # Download from RCSB PDB. Saved with a "_raw" suffix -- NOT the filename
    # run_pipeline.py's Stage 2 writes for its cleaned output -- so Stage 2
    # still performs its own chain-filter / HETATM-strip pass instead of
    # mistaking this raw download for an already-cleaned structure.
    out_dir = workdir if workdir else Path(".")
    dest = out_dir / f"{pdb_id_or_file.lower()}_raw.pdb"
    url = f"https://files.rcsb.org/download/{pdb_id_or_file.upper()}.pdb"
    log(f"Downloading PDB from RCSB: {url} -> {dest}")
    import urllib.request
    urllib.request.urlretrieve(url, dest)
    return str(dest.resolve())


def extract_chain_data(pdb_path: str, target_chain: str):
    """Extracts amino acid sequence and residue IDs for the target chain."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)

    # Find the target chain
    model = structure[0]
    if target_chain not in model:
        avail = [c.id for c in model.get_chains()]
        raise ValueError(f"Chain '{target_chain}' not found in {pdb_path}. Available: {avail}")

    chain = model[target_chain]
    protein_residues = [res for res in chain if is_aa(res)]
    if not protein_residues:
        raise ValueError(f"No amino acid residues found in chain '{target_chain}' of {pdb_path}")

    sequence = "".join(seq1(res.resname) for res in protein_residues)
    sequence_id = [res.id[1] for res in protein_residues]
    resnames_3 = [res.resname for res in protein_residues]

    return sequence, sequence_id, resnames_3


def resolve_model_dir(model_dir: str) -> str:
    """Finds the local prot_t5 weights directory with intelligent fallback."""
    candidates = [
        Path(model_dir) if model_dir else None,
        Path(os.environ.get("PROT_T5_MODEL", "")) if os.environ.get("PROT_T5_MODEL") else None,
        Path(__file__).parent / "prot_t5_xl_uniref50_full_v2",
        Path(__file__).parent / "required_files_weights" / "prot_t5_xl_uniref50_full_v2",
        Path("/opt/models/prot_t5_xl_uniref50_full_v2"),
    ]
    for c in candidates:
        if c and c.is_dir() and (c / "cpt.pth").is_file():
            return str(c.resolve())
    # Return original if none matched to let downstream error report it
    return model_dir


def predict_local(pdb_input: str, chain_id: str, top_n: int, score_type: str, model_dir: str,
                  device: str = None, workdir: Path = None):
    pdb_path = resolve_pdb_file(pdb_input, workdir)
    log(f"Using structure: {pdb_path}")

    sequence, sequence_id, resnames_3 = extract_chain_data(pdb_path, chain_id)
    log(f"Chain {chain_id}: length = {len(sequence)} residues")

    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"
    log(f"Using device: {device}")

    resolved_dir = resolve_model_dir(model_dir)
    log(f"Loading local ProtT5 model from: {resolved_dir}")
    model, tokenizer = load_model(resolved_dir, max_length=1500)
    model.to(device)
    model.eval()
    log("ProtT5 model loaded successfully.")

    log("Running local cryptic residue token classification...")
    spaced_seq = " ".join(sequence)
    input_ids = tokenizer(spaced_seq, return_tensors="pt").input_ids.to(device)

    with torch.no_grad():
        outputs = model(input_ids).logits
        outputs_cpu = outputs.detach().cpu().numpy().squeeze()

        del outputs, input_ids
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

    # If length has BOS/EOS or batch dimension, handle shape
    if outputs_cpu.ndim == 2:
        # Expected shape (seq_len, 2)
        # Tokenizer might add an end token (EOS)
        if len(outputs_cpu) == len(sequence) + 1:
            outputs_cpu = outputs_cpu[:len(sequence)]
        elif len(outputs_cpu) != len(sequence):
            outputs_cpu = outputs_cpu[:len(sequence)]

    raw_scores = expit(outputs_cpu[:, 1] - outputs_cpu[:, 0])
    normalized_scores = normalize_scores(raw_scores)

    scores = normalized_scores if score_type.lower() == "normalized" else raw_scores

    residues = []
    for i in range(len(sequence)):
        r3 = resnames_3[i]
        r1 = sequence[i]
        resid = sequence_id[i]
        sc = round(float(scores[i]), 4)
        residues.append({
            "resid": int(resid),
            "resname_3letter": r3,
            "resname_1letter": r1,
            "chain": chain_id,
            "score": sc,
        })

    residues_sorted = sorted(residues, key=lambda x: x["score"], reverse=True)
    top_residues = residues_sorted[:top_n]

    pdb_id = Path(pdb_input).stem.upper().split("_")[0]

    output_data = {
        "query": {
            "pdb_id": pdb_id,
            "chain": chain_id,
            "score_type": score_type,
        },
        "source": "CryptoBank PLM Local (prot_t5_xl_uniref50_full_v2)",
        "top_residues": top_residues,
        "all_residues": residues_sorted,
        "n_total_residues": len(residues),
    }

    log(f"Top {len(top_residues)} cryptic residues identified:")
    for r in top_residues:
        log(f"  Residue {r['resname_3letter']}{r['resid']} (Chain {r['chain']}) -> score: {r['score']:.4f}")

    return output_data


def main():
    parser = argparse.ArgumentParser(description="CryptoBank Local Offline Predictor")
    parser.add_argument("--pdb", required=True, help="PDB ID or path to .pdb file")
    parser.add_argument("--chain", default="A", help="Target chain ID (default: A)")
    parser.add_argument("--top", type=int, default=5, help="Number of top cryptic residues (default: 5)")
    parser.add_argument("--score_type", default="normalized", choices=["normalized", "raw", "Raw Scores"],
                        help="Score type (default: normalized)")
    parser.add_argument("--model_dir", default=os.environ.get("PROT_T5_MODEL", "/opt/models/prot_t5_xl_uniref50_full_v2"),
                        help="Path to prot_t5_xl_uniref50_full_v2 directory")
    parser.add_argument("--output", default=None, help="Output JSON path")
    parser.add_argument("--device", default=None, help="Compute device ('cuda' or 'cpu')")

    args = parser.parse_args()

    # Keep any RCSB download beside the results rather than in the current directory
    workdir = Path(args.output).parent if args.output else None
    if workdir:
        workdir.mkdir(parents=True, exist_ok=True)

    results = predict_local(
        pdb_input=args.pdb,
        chain_id=args.chain,
        top_n=args.top,
        score_type=args.score_type,
        model_dir=args.model_dir,
        device=args.device,
        workdir=workdir,
    )

    if args.output:
        out_path = Path(args.output)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, "w") as fp:
            json.dump(results, fp, indent=2)
        log(f"Results saved to: {out_path}")
    else:
        print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
