#!/usr/bin/env python3
"""
Point mutagenesis for CryptiScan with PDBFixer + OpenMM.

Drop-in replacement for modeller_mutate.py: same arguments, same output file.

    python openmm_mutate.py modelname respos restyp chain

    Example: python openmm_mutate.py 2wgb 121 ARG A > 2wgb.log

Reads <modelname>.pdb, mutates residue <respos> of <chain> to <restyp> and
writes <modelname><restyp><respos>.pdb to the current directory. Chain IDs and
residue numbers are copied from the input unchanged; Stage 10 maps the mutant
ensemble back to the wild-type chain by them.

Procedure, mirroring what the MODELLER script did:
  1. PDBFixer swaps the residue and builds the new side chain from its
     template. Heavy atoms missing anywhere else are rebuilt as well, as
     MODELLER's build() did. Missing residues (gaps) are never modelled.
  2. Hydrogens are added at pH 7 so the Amber force field can be applied.
  3. Only the side chain of the mutated residue (and of any residue that had
     atoms rebuilt) moves, plus hydrogens. Every other heavy atom, the whole
     backbone included, is held fixed by giving it zero mass.
  4. Minimise, anneal 300 -> 600 -> 300 K, minimise again, in Amber14 with
     GBn2 implicit solvent. The anneal lets the side chain leave the template
     rotamer; MODELLER did the same by randomising it and running MD.
  5. Hydrogens are stripped again, so the output is heavy-atom only, like the
     MODELLER output that clean_structure.py and SAM2 were built around.

PDBFixer and OpenMM are MIT licensed, so unlike MODELLER they can be shipped
inside the image and need no licence key.
"""

import argparse
import math
import os
import sys

import numpy as np
import openmm
from openmm import app, unit
from pdbfixer import PDBFixer

STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}
# Held fixed even inside a mobile residue, so the backbone SAM2 starts from is
# exactly the wild-type backbone.
BACKBONE = {"N", "CA", "C", "O", "OXT"}

ANNEAL_K = (300, 400, 500, 600, 500, 400, 300)
STEPS_PER_TEMP = 500       # x 1 fs. No constraints: they cannot involve a
TIMESTEP_FS = 1.0          # massless (frozen) atom, so the step stays short.

# Fixed seed, like MODELLER's rand_seed: rerunning gives the same model.
SEED = 49837


def log(msg):
    print(f"[MUTATE] {msg}", file=sys.stderr)


def res_key(res):
    """Identifies a residue across topology rebuilds. PDBFixer replaces the
    Topology on every edit, so Residue objects themselves go stale."""
    return (res.chain.id, res.id, res.insertionCode)


def find_residue(topology, chain, respos):
    hits = [r for r in topology.residues()
            if r.chain.id.strip() == chain.strip()
            and r.id == str(respos)
            and not r.insertionCode.strip()]
    if len(hits) != 1:
        log(f"ERROR: expected one residue {respos} in chain '{chain}', found {len(hits)}")
        sys.exit(1)
    return hits[0]


def energy(sim):
    e = sim.context.getState(getEnergy=True).getPotentialEnergy()
    return e.value_in_unit(unit.kilojoule_per_mole)


def main():
    p = argparse.ArgumentParser(
        description="CryptiScan point mutagenesis (PDBFixer + OpenMM)")
    p.add_argument("modelname", help="input is <modelname>.pdb")
    p.add_argument("respos", help="residue number to mutate")
    p.add_argument("restyp", help="target residue, 3-letter code")
    p.add_argument("chain", help="chain ID ('' if the file has none)")
    p.add_argument("--platform", default="CPU",
                   help="OpenMM platform (default CPU: deterministic, needs no GPU)")
    p.add_argument("--seed", type=int, default=SEED)
    a = p.parse_args()

    restyp = a.restyp.upper()
    if restyp not in STANDARD_AA:
        log(f"ERROR: '{a.restyp}' is not a standard amino acid")
        sys.exit(1)

    src = a.modelname + ".pdb" if os.path.exists(a.modelname + ".pdb") else a.modelname
    out = f"{a.modelname}{restyp}{a.respos}.pdb"
    log(f"input  {src}")

    fixer = PDBFixer(filename=src)
    target = find_residue(fixer.topology, a.chain, a.respos)
    wt, target_key = target.name, res_key(target)

    # --- 1. mutate and rebuild missing heavy atoms ------------------------
    if wt == restyp:
        log(f"{target.chain.id}/{wt}{a.respos} is already {restyp}; relaxing only")
    else:
        fixer.applyMutations([f"{wt}-{a.respos}-{restyp}"], target.chain.id)
        log(f"mutate {target.chain.id}/{wt}{a.respos} -> {restyp}")

    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findNonstandardResidues()
    if fixer.nonstandardResidues:
        log("replacing non-standard residues: " + ", ".join(
            f"{r.name}{r.id}->{new}" for r, new in fixer.nonstandardResidues))
    fixer.replaceNonstandardResidues()
    # Only reached with --no_clean input; the force field has no templates
    # for ligands, and the SAM2 input is protein-only anyway.
    fixer.removeHeterogens(False)
    fixer.findMissingAtoms()

    # The mutated residue is listed here too: its new side chain is "missing".
    mobile_res = {res_key(r) for r in fixer.missingAtoms} | {target_key}
    others = len(mobile_res) - 1
    if others:
        log(f"rebuilding missing heavy atoms in {others} other residue(s); "
            f"their side chains are relaxed as well")
    fixer.addMissingAtoms(seed=a.seed)

    # --- 2. hydrogens, force field ----------------------------------------
    fixer.addMissingHydrogens(7.0)
    top = fixer.topology
    top.setPeriodicBoxVectors(None)   # CRYST1 from RCSB is meaningless here

    ff = app.ForceField("amber14-all.xml", "implicit/gbn2.xml")
    system = ff.createSystem(top, nonbondedMethod=app.CutoffNonPeriodic,
                             nonbondedCutoff=1.2 * unit.nanometer,
                             constraints=None, removeCMMotion=False)

    # --- 3. freeze everything but the mobile side chains and hydrogens ----
    atoms = list(top.atoms())
    mobile = {at.index for at in atoms
              if at.element == app.element.hydrogen
              or (res_key(at.residue) in mobile_res and at.name not in BACKBONE)}
    for i in range(system.getNumParticles()):
        if i not in mobile:
            system.setParticleMass(i, 0.0)
    n_heavy = sum(1 for i in mobile if atoms[i].element != app.element.hydrogen)
    log(f"atoms  {len(atoms)} total, {len(mobile)} mobile ({n_heavy} heavy)")

    # --- 4. minimise, anneal, minimise ------------------------------------
    integrator = openmm.LangevinMiddleIntegrator(
        ANNEAL_K[0] * unit.kelvin, 1.0 / unit.picosecond, TIMESTEP_FS * unit.femtosecond)
    integrator.setRandomNumberSeed(a.seed)
    platform = openmm.Platform.getPlatformByName(a.platform)
    sim = app.Simulation(top, system, integrator, platform)
    sim.context.setPositions(fixer.positions)
    start = sim.context.getState(getPositions=True).getPositions(asNumpy=True)

    log(f"energy {energy(sim):14.1f} kJ/mol  (as built)")
    sim.minimizeEnergy()
    log(f"energy {energy(sim):14.1f} kJ/mol  (minimised)")
    minimised = sim.context.getState(getPositions=True).getPositions()

    try:
        sim.context.setVelocitiesToTemperature(ANNEAL_K[0] * unit.kelvin, a.seed)
        for t in ANNEAL_K:
            integrator.setTemperature(t * unit.kelvin)
            sim.step(STEPS_PER_TEMP)
        sim.minimizeEnergy()
        if not math.isfinite(energy(sim)):
            raise openmm.OpenMMException("energy is not finite")
        log(f"energy {energy(sim):14.1f} kJ/mol  (annealed {'-'.join(map(str, ANNEAL_K))} K, "
            f"minimised)")
    except openmm.OpenMMException as exc:
        log(f"WARNING: annealing failed ({exc}); keeping the minimised structure")
        sim.context.setPositions(minimised)

    final = sim.context.getState(getPositions=True).getPositions(asNumpy=True)
    frozen = np.array([i for i in range(system.getNumParticles()) if i not in mobile])
    # Strip units before the numpy maths: np.abs drops them from a Quantity.
    drift = np.abs(final.value_in_unit(unit.angstrom)[frozen]
                   - start.value_in_unit(unit.angstrom)[frozen]).max()
    log(f"max displacement of fixed atoms: {drift:.2e} A")

    # --- 5. strip hydrogens and write -------------------------------------
    # openmm.app.Modeller is OpenMM's topology editor, unrelated to MODELLER.
    edit = app.Modeller(top, sim.context.getState(getPositions=True).getPositions())
    edit.delete([at for at in edit.topology.atoms()
                 if at.element == app.element.hydrogen])
    with open(out, "w") as fh:
        app.PDBFile.writeFile(edit.topology, edit.positions, fh, keepIds=True)
    log(f"output {out}")


if __name__ == "__main__":
    main()
