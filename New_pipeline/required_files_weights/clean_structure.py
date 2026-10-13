#!/usr/bin/env python3
"""
Structure preparation for CryptiScan (MDAnalysis).

Two modes:

  --mode split     Split a raw RCSB entry into
                     <base>_protein.pdb             protein, ALL chains, standard residues
                     <base>_ligands/<RES>_<CH><N>.pdb   one file per ligand (any resname)
                   Waters are dropped. Alternate locations are collapsed to one.

  --mode extract   Pull a single chain out of a structure and validate that it is
                   safe to feed to SAM2 (standard residues only, every residue has CA).

Ligand handling is generic: anything that is not one of the 20 standard amino
acids and not water is treated as a ligand, whatever its residue name.
"""

import argparse
import json
import sys
import warnings
from pathlib import Path

import MDAnalysis as mda

# The 20 standard amino acids. Deliberately NOT MDAnalysis' broader "protein"
# selection: SAM2's atom14 table keys on exactly these names, so protonation
# variants (HSD/HSE/HID/HIP...) would raise KeyError there.
STANDARD_AA = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
]
WATER = ["HOH", "WAT", "TIP3", "TIP4", "SOL", "DOD"]

PROTEIN_SEL = "resname " + " ".join(STANDARD_AA)
WATER_SEL = "resname " + " ".join(WATER)


def log(msg):
    print(f"[PREP] {msg}", file=sys.stderr)


def load(pdb_path):
    """Load a structure, keeping only the first model of a multi-model entry."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        u = mda.Universe(str(pdb_path))
    if len(u.trajectory) > 1:
        u.trajectory[0]
        log(f"multi-model entry: using model 1 of {len(u.trajectory)}")
    return u


def drop_altlocs(atoms):
    """Keep only the blank or first-listed alternate location for each atom."""
    if not hasattr(atoms, "altLocs"):
        return atoms
    keep = (atoms.altLocs == "") | (atoms.altLocs == "A")
    n_drop = int((~keep).sum())
    if n_drop:
        log(f"dropped {n_drop} atoms in alternate locations (kept blank/'A')")
    return atoms[keep]


def residue_tag(res):
    """Filesystem-safe identifier: RESNAME_<chain><resid><icode>."""
    chain = (getattr(res.atoms, "chainIDs", [""])[0] or "X").strip() or "X"
    icode = ""
    try:
        icode = (res.atoms.icodes[0] or "").strip()
    except (AttributeError, IndexError):
        pass
    return f"{res.resname.strip()}_{chain}{res.resid}{icode}"


def write(atoms, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        atoms.write(str(path))


def do_split(raw_pdb, outdir, basename):
    u = load(raw_pdb)
    outdir = Path(outdir)

    protein = drop_altlocs(u.select_atoms(PROTEIN_SEL))
    if protein.n_atoms == 0:
        log("ERROR: no standard amino-acid atoms found in the structure.")
        sys.exit(1)

    protein_path = outdir / f"{basename}_protein.pdb"
    write(protein, protein_path)
    chains = sorted({c for c in protein.chainIDs if c.strip()})
    log(f"protein  -> {protein_path.name} "
        f"({protein.n_atoms} atoms, {protein.n_residues} residues, chains: {','.join(chains) or 'n/a'})")

    # Anything that is neither a standard amino acid nor water is a ligand,
    # regardless of its residue name.
    other = u.select_atoms(f"not ({PROTEIN_SEL}) and not ({WATER_SEL})")
    n_water = u.select_atoms(WATER_SEL).n_residues

    ligands = []
    if other.n_atoms:
        lig_dir = outdir / f"{basename}_ligands"
        for res in other.residues:
            tag = residue_tag(res)
            lig_path = lig_dir / f"{tag}.pdb"
            write(res.atoms, lig_path)
            ligands.append({
                "resname": res.resname.strip(),
                "resid": int(res.resid),
                "chain": (getattr(res.atoms, "chainIDs", [""])[0] or "").strip(),
                "n_atoms": int(res.atoms.n_atoms),
                "file": str(lig_path.relative_to(outdir)),
            })
            log(f"ligand   -> {lig_path.relative_to(outdir)} ({res.atoms.n_atoms} atoms)")

    if n_water:
        log(f"waters   -> dropped ({n_water} residues)")
    if not ligands:
        log("ligands  -> none found")

    manifest = {
        "source": str(raw_pdb),
        "protein_file": protein_path.name,
        "chains": chains,
        "n_protein_residues": int(protein.n_residues),
        "n_waters_dropped": int(n_water),
        "ligands": ligands,
    }
    man_path = outdir / f"{basename}_prep.json"
    man_path.write_text(json.dumps(manifest, indent=2))
    log(f"manifest -> {man_path.name}")
    return manifest


def do_extract(in_pdb, chain, out_pdb):
    """Single chain for SAM2, with the checks SAM2 itself would fail on."""
    u = load(in_pdb)
    sel = u.select_atoms(f"({PROTEIN_SEL}) and chainID {chain}")
    if sel.n_atoms == 0:
        avail = sorted({c for c in u.atoms.chainIDs if c.strip()})
        log(f"ERROR: chain '{chain}' has no standard protein atoms. Available: {avail}")
        sys.exit(1)
    sel = drop_altlocs(sel)

    # SAM2 raises KeyError on any non-standard residue and ValueError on any
    # residue lacking CA. Fail here instead, with a message that says why.
    bad_names = sorted({r.resname.strip() for r in sel.residues
                        if r.resname.strip() not in STANDARD_AA})
    if bad_names:
        log(f"ERROR: non-standard residues would break SAM2: {bad_names}")
        sys.exit(1)

    missing_ca = [f"{r.resname.strip()}{r.resid}" for r in sel.residues
                  if "CA" not in set(r.atoms.names)]
    if missing_ca:
        log(f"ERROR: {len(missing_ca)} residue(s) have no CA atom, SAM2 requires one: "
            f"{missing_ca[:10]}{' ...' if len(missing_ca) > 10 else ''}")
        sys.exit(1)

    write(sel, Path(out_pdb))
    log(f"SAM2 input -> {Path(out_pdb).name} "
        f"(chain {chain}, {sel.n_residues} residues, {sel.n_atoms} atoms; validated)")


def main():
    p = argparse.ArgumentParser(description="CryptiScan structure preparation (MDAnalysis)")
    p.add_argument("--mode", choices=["split", "extract"], required=True)
    p.add_argument("--pdb", required=True, help="Input structure")
    p.add_argument("--outdir", default=".", help="Output directory (split mode)")
    p.add_argument("--basename", help="Output basename (split mode)")
    p.add_argument("--chain", help="Chain to extract (extract mode)")
    p.add_argument("--output", help="Output file (extract mode)")
    a = p.parse_args()

    if a.mode == "split":
        if not a.basename:
            p.error("--basename is required for --mode split")
        do_split(a.pdb, a.outdir, a.basename)
    else:
        if not (a.chain and a.output):
            p.error("--chain and --output are required for --mode extract")
        do_extract(a.pdb, a.chain, a.output)


if __name__ == "__main__":
    main()
