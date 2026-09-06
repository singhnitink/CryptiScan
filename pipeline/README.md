# CryptiScan Pipeline

This directory contains the unified, end-to-end CryptiScan workflow.

```
                      PDB ID  or  local .pdb
                                  │
                                  ▼
  ┌───────────────────────────────────────────────────────────────┐
  │ 1. Fetch from RCSB + Cleanup        [MDAnalysis]              │
  │    <pdb>_raw.pdb       entry as downloaded                    │
  │    <pdb>_protein.pdb   full protein, ALL chains               │
  │    <pdb>_ligands/      one file per ligand, any resname       │
  │    waters dropped, alternate locations collapsed              │
  └───────────────────────────────┬───────────────────────────────┘
                                  ▼
  ┌───────────────────────────────────────────────────────────────┐
  │ 2. Chain Selection                  [MDAnalysis]              │
  │    <pdb>.pdb   the chain set by CHAIN in the launcher.        │
  │    Standard residues + CA verified. Every stage below         │
  │    operates on this single chain.                             │
  └───────────────────────────────┬───────────────────────────────┘
                                  ▼
  ┌───────────────────────────────────────────────────────────────┐
  │ 3. Cryptic Residue Prediction                                 │
  │    Local ProtT5 (offline)  or  CryptoBank Space (web)         │
  │    -> <PDB>_<chain>_cryptic.json                              │
  └───────────────────────────────┬───────────────────────────────┘
                                  ▼
  ┌───────────────────────────────────────────────────────────────┐
  │ 4. ESM-Scan                             [ESM-1v]              │
  │    Scores all 19 substitutions at each cryptic site,          │
  │    keeps the highest-scoring mutation per position            │
  └───────────────────────────────┬───────────────────────────────┘
                                  ▼
  ┌───────────────────────────────────────────────────────────────┐
  │ 5. In Silico Mutagenesis              [MODELLER]              │
  │    Mutates and relaxes sidechains                             │
  │    -> <pdb>_top<N>_esm_mutant.pdb                             │
  └───────────────────────────────┬───────────────────────────────┘
                                  ▼
  ┌───────────────────────────────────────────────────────────────┐
  │ 5b. Validate for SAM2               [MDAnalysis]              │
  │     -> sam2_input.pdb                                         │
  └───────────────────────────────┬───────────────────────────────┘
                                  ▼
  ┌───────────────────────────────────────────────────────────────┐
  │ 6. Ensemble Generation             [aSAM / SAM2]              │
  │    GPU conformational sampling -> ensemble_output.*           │
  └───────────────────────────────┬───────────────────────────────┘
                                  ▼
  ┌───────────────────────────────────────────────────────────────┐
  │ 7. Report + Sequences                                         │
  │    report.txt, <PDB>_<chain>_wildtype.fasta,                  │
  │    <PDB>_<chain>_mutant.fasta                                 │
  └───────────────────────────────┬───────────────────────────────┘
                                  ▼
  ┌───────────────────────────────────────────────────────────────┐
  │ 8. Archive  ->  <jobname>_<timestamp>.zip                     │
  └───────────────────────────────────────────────────────────────┘
```

---

## Two Versions

| | **Local Version** | **Web Version** |
|---|---|---|
| **Launcher script** | `launch_pipeline_local.sh` | `launch_pipeline_web.sh` |
| **Apptainer def** | `pipeline_local.def` | `pipeline_web.def` |
| **Container image** | `pipeline_local.sif` (~20 GB) | `pipeline_web.sif` (~12 GB) |
| **Stage 1 engine** | Local ProtT5 (`predict_cryptic_local.py`) | CryptoBank HF Space (`scrape_cryptobank.py`) |
| **Internet needed** | Only to fetch the structure from RCSB (see note) | Yes — stage 1 calls `thorbenf-cryptobank.hf.space` |
| **Baked-in weights** | ESM-1v, ProtT5-XL + head, SAM2 | ESM-1v, SAM2 |

> **Note on "offline".** The local version performs cryptic-pocket *inference*
> locally, with no HuggingFace Space call. Stage 2 still downloads the target
> structure from `files.rcsb.org`. To run with no network at all, set `PDB` in
> the launcher to the path of a local `.pdb` file.

Both versions are fully self-contained. **The launchers pass no `--bind` flags,
by design** — every script and weight is read from inside the `.sif`. A bind
mount would silently mask a stale image with host files that whoever you share
the image with will not have. Please do not add them back.

---

## How to Run

Everything is driven by the bash launcher. Open it, edit the configuration block
at the top, and run it.

```bash
# 1. Edit PDB, CHAIN, TOP, MODE and MODELLER_KEY at the top of the script
# 2. Run:
bash launch_pipeline_local.sh    # local version
bash launch_pipeline_web.sh      # web version
```

The configuration block looks like this:

```bash
PDB="1JWP"                     # PDB accession code, or a path to a local .pdb
CHAIN="A"                      # Target chain identifier
TOP=5                          # Number of top cryptic residues to mutate
MODE="combined"                # "combined" = one structure; "independent" = N mutants
STRATEGY="esm"                 # "esm" = pick substitutions with ESM-Scan
JOBNAME="${JOBNAME:-1jwp}"     # Names the output folder and zip
MODELLER_KEY="${MODELLER_KEY:-MODELIRANJE}"
```

`JOBNAME` and `MODELLER_KEY` can also be set from the environment without
editing the file:

```bash
JOBNAME=1lzt MODELLER_KEY=your_key_here bash launch_pipeline_local.sh
```

Before stage 1, the launcher verifies that the image contains and can read every
file it needs. If the `.sif` is stale or was built from an older definition, it
stops immediately and tells you to rebuild.

### Output

Everything is packaged into a single timestamped archive:

```
<jobname>_<timestamp>.zip
├── 1jwp_raw.pdb                     entry as downloaded from RCSB
├── 1jwp_protein.pdb                 full protein, ALL chains
├── 1jwp_ligands/PO4_A1.pdb          one file per ligand, whatever its resname
├── 1jwp_prep.json                   what was kept, split out, and dropped
├── 1jwp.pdb                         SELECTED CHAIN - used by every stage below
├── 1JWP_A_cryptic.json              cryptic residue predictions
├── 1jwp_esm_summary.csv             all 19 substitutions scored per site
├── 1jwp_top5_esm_mutant.pdb         MODELLER mutant
├── sam2_input.pdb                   mutant, validated for SAM2
├── ensemble_output.*                SAM2 / aSAM conformational ensemble
├── 1JWP_manifest.json               run manifest
├── report.txt                       human-readable summary of the whole run
├── 1JWP_A_wildtype.fasta            wild-type sequence of the selected chain
└── 1JWP_A_mutant.fasta              mutant sequence
```

**Stage 1 — cleanup.** MDAnalysis splits the raw entry: the full protein (all
chains, the 20 standard residues) goes to `<pdb>_protein.pdb`, every
non-protein non-water residue is written to its own file under
`<pdb>_ligands/` regardless of residue name, and waters are dropped.
Alternate locations are collapsed to one. Ligands are kept for your own
analysis and are not fed downstream.

**Stage 2 — chain selection.** `<pdb>.pdb` is reduced to the chain named by
`CHAIN` in the launcher, and validated: standard residues only, every residue
has a CA. Cryptic-residue prediction, ESM-Scan, MODELLER and SAM2 all run on
this single chain. The full protein stays available in `<pdb>_protein.pdb` if
you want to push the whole structure to CryptoBank instead.

Selecting one chain is also what keeps SAM2 correct: it flattens residues
across chains with no error, so more than one chain yields a quietly wrong
ensemble rather than a crash. Stage 5b re-runs the same validation on the
MODELLER output before the ensemble stage.

---

## Prerequisites

- **Linux** (x86_64)
- **NVIDIA GPU** with a standard driver (version ≥ 525)
- **Apptainer** (or Singularity)
- A free **MODELLER license key** (academic registration at
  [salilab.org](https://salilab.org/modeller/registration.html))

Nothing else. CUDA 12.5, PyTorch 2.4, ESM-1v, ProtT5, MODELLER 10.8 and SAM2 all
live inside the container — no host Python, no conda environment, no `pip
install`, no weights to download.

---

## Files in this Directory

| File | Purpose |
|---|---|
| `launch_pipeline_local.sh` | Launcher for the local version. |
| `launch_pipeline_web.sh` | Launcher for the web version. |
| `pipeline_local.def` | Apptainer definition for the local image. |
| `pipeline_web.def` | Apptainer definition for the web image. |
| `required_files_weights/` | Everything baked into the images (see manifest below). |

### `required_files_weights/` manifest

This directory is the single source of truth for the container contents. Check
these sizes after copying the folder from OneDrive — a truncated download is the
most common cause of a failed build.

| Path | Size |
|---|---|
| `run_pipeline.py` | 16 KB |
| `predict_cryptic_local.py` | 8.6 KB |
| `cryptobank_model_loader.py` | 27 KB |
| `scrape_cryptobank.py` | 6.9 KB |
| `esm_scanner.py` | 9.8 KB |
| `modeller_mutate.py` | 5.7 KB |
| `clean_structure.py` | 7.0 KB |
| `make_report.py` | 11 KB |
| `esm1v_t33_650M_UR90S_1.pt` | 7.3 GB |
| `prot_t5_xl_uniref50_full_v2/` | 6.8 GB |
| `sam2_weights/weights/mdcath_1.0/` | 1.0 GB |
| **total** | **~15 GB** |

Quick check:

```bash
du -sh required_files_weights          # expect ~15 G
find required_files_weights -type l    # expect no output (no symlinks)
```

Script roles: `run_pipeline.py` orchestrates stages 1–4;
`predict_cryptic_local.py` runs local ProtT5 inference;
`cryptobank_model_loader.py` is the ProtT5 classification-head architecture;
`scrape_cryptobank.py` queries the HF Space; `esm_scanner.py` runs ESM-1v;
`modeller_mutate.py` does the point mutagenesis and sidechain minimization;
`clean_structure.py` does the MDAnalysis structure preparation;
`make_report.py` builds the run report and FASTA files.

### Two conda environments inside each image

| Env | Path | Holds |
|---|---|---|
| main | `/opt/conda` | torch 2.4.0+cu124, ESM-1v, MODELLER 10.8, SAM2, mdtraj, (local only: transformers/ProtT5) |
| prep | `/opt/conda/envs/prep` | MDAnalysis 2.7.0 only |

MDAnalysis is deliberately isolated. The main env pins `numpy==1.26.4` because
`torch 2.4.0` and `mdtraj 1.10.3` require it, and a shared install risks a
dependency pulling numpy to 2.x and silently breaking SAM2. `run_pipeline.py`
calls into the prep env through `--prep_python`, the same way it already uses
separate interpreters for the scraper, ESM and MODELLER stages. `apptainer test`
asserts the main env's numpy is still 1.26.4 after the prep env is created.

---

## How ESM-Scan Scoring Works

Rather than mutating cryptic residues to a generic amino acid (like `GLU` or
`ALA`), **ESM-Scan** leverages Facebook Research's pre-trained **ESM-1v** protein
language model:

1. The target chain sequence is passed through ESM-1v.
2. For each position, the model computes log-probabilities for all 20 standard amino acids.
3. For any mutation $\text{WT} \to \text{Mutant}$:
   $$\text{Score} = \log P(\text{Mutant}) - \log P(\text{WT})$$
4. All 19 alternative amino acids are evaluated at each cryptic site. The amino
   acid with the **highest (best) score** is chosen for the mutation.
5. All scores are saved to `<PDB>_esm_summary.csv` for full transparency.

---

## Building the Container Images

Build from inside this directory, so the relative `%files` paths resolve:

```bash
cd pipeline/
sudo apptainer build pipeline_web.sif   pipeline_web.def     # ~12 GB, faster
sudo apptainer build pipeline_local.sif pipeline_local.def   # ~20 GB
chmod 0644 pipeline_local.sif pipeline_web.sif
```

Building needs internet access (base image, conda, pip, the SAM2 clone) and
roughly 60 GB of free scratch space.

Verify an image before using or sharing it:

```bash
apptainer test pipeline_local.sif        # checks packages, and that weights are readable
apptainer run-help pipeline_local.sif    # what the image is and how to drive it
apptainer inspect --deffile pipeline_local.sif | diff - pipeline_local.def   # provenance
```

That last command is worth the habit: it proves the `.sif` was actually built
from the `.def` sitting next to it.

---

## Sharing this Folder

Copy the whole `pipeline/` directory. A recipient can either run the prebuilt
`.sif` files or rebuild them from the `.def` files.

Two things to do after downloading from OneDrive:

```bash
chmod +x launch_pipeline_local.sh launch_pipeline_web.sh   # OneDrive strips the executable bit
du -sh required_files_weights                        # confirm the transfer completed
```

The images are built with the baked weights world-readable, so they work
regardless of which user account runs them.

> **MODELLER license.** The launchers ship with a placeholder key. Sali Lab
> licenses are issued per user — register for your own free academic key at
> [salilab.org](https://salilab.org/modeller/registration.html) and set it via
> `MODELLER_KEY` rather than reusing someone else's.
