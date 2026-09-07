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
4. All 19 alternatives are scored at every target residue, but the mutation is
   chosen from a restricted set (--mutation_set, default charged_polar): the one
   with the HIGHEST score within that set wins.

Why the set is restricted:
--------------------------
ESM-1v measures evolutionary fitness, so unrestricted it picks whatever best
preserves the fold -- for a buried hydrophobic residue that means another
hydrophobic one (ILE -> VAL). That defeats the purpose here. These mutations
exist to destabilise the closed state so the cryptic pocket opens during MD,
which needs a charged or polar residue buried in a hydrophobic environment.
Scores will therefore be markedly negative; that is expected, not a failure.

A residue whose wild type is already charged/polar is scored but NOT mutated
(GLU -> ASP would not open anything). Such residues are reported with
skip_mutation: true and are excluded from the MODELLER stage.

Usage (Standalone):
-------------------
    python3 esm_scanner.py \
        --pdb 1JWP.pdb \
        --chain A \
        --residues 26,45,102 \
        --mutation_set charged_polar \
        --model_path /opt/models/esm1v_t33_650M_UR90S_1.pt \
        --output_json esm_results.json
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

# ---------------------------------------------------------------------------
# Substitution sets
# ---------------------------------------------------------------------------
# ESM-1v is an evolutionary fitness model: unrestricted, it favours whatever
# best preserves the fold, which for a buried hydrophobic residue means another
# hydrophobic one (ILE -> VAL). That is useless here -- the point of mutating a
# cryptic site is to DESTABILISE the closed state so the pocket opens in MD.
# Burying a charge or a polar group costs desolvation energy and adds steric
# strain, which is what actually drives opening. Restricting the candidates
# lets ESM still rank within that set, so we take the least evolutionarily
# disruptive of the destabilising options rather than an arbitrary one.
CHARGED = "DERKH"   # ASP GLU ARG LYS HIS (His is titratable; polar regardless)
POLAR = "NQSTY"     # ASN GLN SER THR TYR
# CYS is deliberately absent: a free cysteine invites spurious disulfides and MD
# artefacts. TRP is aromatic/nonpolar. Both remain reachable via a custom list.
MUTATION_SETS = {
    "charged": CHARGED,
    "polar": POLAR,
    "charged_polar": CHARGED + POLAR,
    "all": "".join(ALL_AMINO_ACIDS),
}


def resolve_mutation_set(spec: str) -> List[str]:
    """Turn a preset name or a custom comma list into 1-letter codes.

    Accepts either form, so '--mutation_set ASP,GLU' and '--mutation_set D,E'
    are equivalent.
    """
    spec = (spec or "").strip()
    if spec.lower() in MUTATION_SETS:
        return sorted(set(MUTATION_SETS[spec.lower()]))

    out = set()
    for tok in spec.replace(" ", "").split(","):
        if not tok:
            continue
        t = tok.upper()
        if len(t) == 3 and t in THREE_TO_ONE:
            out.add(THREE_TO_ONE[t])
        elif len(t) == 1 and t in ONE_TO_THREE:
            out.add(t)
        else:
            log(f"ERROR: '{tok}' is not a standard amino acid (1- or 3-letter).")
            log(f"       Presets: {', '.join(sorted(MUTATION_SETS))}")
            sys.exit(1)
    if not out:
        log(f"ERROR: --mutation_set '{spec}' resolved to no amino acids.")
        sys.exit(1)
    return sorted(out)


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
    selection_set: List[str] = None,
) -> List[Dict]:
    """
    Runs ESM-1v on the sequence and calculates mutation scores for each target residue.

    Parameters:
        sequence: Wild-type protein sequence (1-letter code).
        target_resids: List of PDB residue numbers to evaluate.
        resid_to_seqidx: Mapping from PDB resid to 0-based index in sequence.
        model_path: Path to pre-trained ESM-1v weights (.pt file).
        device: 'cuda' or 'cpu' (auto-detects GPU if available).
        selection_set: 1-letter codes the mutation may be chosen from. A residue
                       whose wild type is already in this set is scored but
                       flagged skip_mutation, since swapping like for like
                       (e.g. GLU->ASP) will not open a pocket.

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

        # Score every alternative once; the selection set only decides which of
        # them may be chosen, so the full 19 stay available as a record.
        def score_of(mt_aa1):
            mt_tok = alphabet.get_idx(mt_aa1)
            mt_prob = token_probs[0, 1 + seq_idx, mt_tok].item()
            # Zero-shot score: difference in log-probability (mutant - WT)
            return round(float(mt_prob - wt_prob), 4)

        all_candidates, candidates = [], []
        for mt_aa1 in ALL_AMINO_ACIDS:
            if mt_aa1 == wt_aa1:
                continue  # Skip wild-type
            entry = {
                "mut_aa1": mt_aa1,
                "mut_aa3": ONE_TO_THREE[mt_aa1],
                "esm_score": score_of(mt_aa1),
            }
            all_candidates.append(entry)
            if mt_aa1 in selection_set:
                candidates.append(dict(entry))

        # Sort descending: highest score = most favored mutation
        all_candidates.sort(key=lambda x: x["esm_score"], reverse=True)
        candidates.sort(key=lambda x: x["esm_score"], reverse=True)

        # A wild type that is already charged/polar cannot be usefully swapped
        # for another one -- GLU->ASP will not open anything -- so score it for
        # the report but do not hand it to MODELLER.
        wt_in_set = wt_aa1 in selection_set
        if wt_in_set:
            log(f"Residue {wt_aa3}{resid} -> SKIPPED "
                f"(wild-type is already in the selection set)")
            best = None
        elif not candidates:
            # Cannot happen with the built-in presets, but a custom set of one
            # residue equal to the wild type would empty the list.
            log(f"Residue {wt_aa3}{resid} -> SKIPPED (no candidates in the selection set)")
            best = None
        else:
            best = candidates[0]
            log(f"Residue {wt_aa3}{resid} -> Best ESM mutation: "
                f"{best['mut_aa3']} (score: {best['esm_score']:+.4f})")

        results.append({
            "resid": resid,
            "wt_aa1": wt_aa1,
            "wt_aa3": wt_aa3,
            "best_mut_aa1": best["mut_aa1"] if best else None,
            "best_mut_aa3": best["mut_aa3"] if best else None,
            "best_score": best["esm_score"] if best else None,
            "skip_mutation": best is None,
            "wt_in_selection_set": wt_in_set,
            "selection_set": list(selection_set),
            "candidates_ranked": candidates,
            "all_candidates_ranked": all_candidates,
        })

    return results


def save_reports(results: List[Dict], output_json: str, output_csv: str):
    """Saves formatted JSON and CSV summaries of ESM predictions."""
    # 1. Save detailed JSON
    with open(output_json, "w") as f:
        json.dump(results, f, indent=2)
    log(f"Saved JSON report -> {output_json}")

    # 2. Save human-readable CSV summary: one column per candidate in the
    #    selection set, so every charged/polar score is visible at a glance
    #    -- including for residues that were skipped.
    sel = []
    for r in results:
        for a in r.get("selection_set", []):
            if a not in sel:
                sel.append(a)
    sel_cols = [ONE_TO_THREE[a] for a in sorted(sel)]

    with open(output_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["resid", "wt", "status", "best_mutant", "best_score"] + sel_cols)
        for r in results:
            scores = {c["mut_aa3"]: c["esm_score"] for c in r["all_candidates_ranked"]}
            status = "skipped" if r.get("skip_mutation") else "mutated"
            writer.writerow(
                [r["resid"], r["wt_aa3"], status,
                 r["best_mut_aa3"] or "", r["best_score"] if r["best_score"] is not None else ""]
                # blank where the candidate IS the wild type, so the gap is visible
                + [scores.get(c, "") for c in sel_cols]
            )
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
    parser.add_argument("--mutation_set", default="charged_polar",
                        help="Amino acids the mutation may be chosen from: a preset "
                             "(" + ", ".join(sorted(MUTATION_SETS)) + ") or a custom "
                             "comma list in 1- or 3-letter form (e.g. 'D,E' or 'ASP,GLU'). "
                             "Restricting to charged/polar residues is what destabilises "
                             "the closed state so the pocket can open; 'all' restores the "
                             "unrestricted 19-way scan.")
    parser.add_argument("--output_json", default="esm_scan_results.json", help="Output JSON path")
    parser.add_argument("--output_csv", default="esm_scan_summary.csv", help="Output CSV path")

    args = parser.parse_args()

    selection_set = resolve_mutation_set(args.mutation_set)
    log(f"Mutation set '{args.mutation_set}': "
        f"{', '.join(ONE_TO_THREE[a] for a in selection_set)}")

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
        selection_set=selection_set,
    )

    n_skip = sum(1 for r in results if r.get("skip_mutation"))
    if n_skip:
        log(f"{n_skip} of {len(results)} residues skipped "
            f"(wild type already in the selection set); they are scored in the "
            f"reports but will not be mutated.")

    # Save reports
    save_reports(results, args.output_json, args.output_csv)
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
