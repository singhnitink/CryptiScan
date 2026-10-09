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
  │    Scores every substitution, then picks the best CHARGED     │
  │    or POLAR one. Conservative nonpolar swaps (ILE->VAL)       │
  │    would leave the pocket shut.                               │
  └───────────────────────────────┬───────────────────────────────┘
                                  ▼
  ┌───────────────────────────────────────────────────────────────┐
  │ 5. In Silico Mutagenesis              [MODELLER]              │
  │    Mutates and relaxes sidechains                             │
  │    -> <pdb>_top<N>_esm_mutant.pdb                             │
  └───────────────────────────────┬───────────────────────────────┘
                                  ▼
  ┌───────────────────────────────────────────────────────────────┐
  │ 5b. Validate for aSAM               [MDAnalysis]              │
  │     -> sam2_input.pdb                                         │
  └───────────────────────────────┬───────────────────────────────┘
                                  ▼
  ┌───────────────────────────────────────────────────────────────┐
  │ 6. Ensemble Generation                     [aSAM]             │
  │    GPU conformational sampling -> ensemble_output.*           │
  └───────────────────────────────┬───────────────────────────────┘
                                  ▼
  ┌───────────────────────────────────────────────────────────────┐
  │ 6b. WILD-TYPE Ensemble                     [aSAM]             │
  │     The SAME sampling, run on the UNMUTATED chain.            │
  │     -> wt_ensemble.*    NOT optional -- see "Why the          │
  │     wild-type ensemble is mandatory" below.                   │
  └───────────────────────────────┬───────────────────────────────┘
                                  ▼
  ┌───────────────────────────────────────────────────────────────┐
  │ 7. Ensemble Clustering          [MDAnalysis + scipy]          │
  │     PCA + k-means on CA, one frame per cluster                │
  │     -> ensemble_clusters/                                     │
  └───────────────────────────────┬───────────────────────────────┘
                                  ▼
  ┌───────────────────────────────────────────────────────────────┐
  │ 8. Report + Sequences                                         │
  │    <jobname>_report.txt, <PDB>_<chain>_wildtype.fasta,        │
  │    <PDB>_<chain>_mutant.fasta                                 │
  └───────────────────────────────┬───────────────────────────────┘
                                  ▼
  ┌───────────────────────────────────────────────────────────────┐
  │ 9. Quantification         [mdtraj + fpocket]                  │
  │    Per frame, in BOTH ensembles:                              │
  │      * SASA of the predicted residues   (detector-free)       │
  │      * cryptic-pocket volume + druggability  (fpocket)        │
  │      * specificity + fold-integrity controls                  │
  │    -> pocket_analysis/                                        │
  └───────────────────────────────┬───────────────────────────────┘
                                  ▼
  ┌───────────────────────────────────────────────────────────────┐
  │ 9b. Figures                      [matplotlib]                 │
  │     -> figures/*.png + *.pdf + stats_summary.{csv,txt}        │
  └───────────────────────────────┬───────────────────────────────┘
                                  ▼
  ┌───────────────────────────────────────────────────────────────┐
  │ 10. Archive  ->  <jobname>_<timestamp>.zip                    │
  └───────────────────────────────────────────────────────────────┘
```

---

## What changed in this revision

If you are coming from an earlier copy of this folder, these are the differences
that matter.

**A wild-type ensemble is now generated and is mandatory.** Stage 6b runs the
same sampling on the unmutated chain. This doubles GPU time and is the single
most important change in the pipeline.

**Stage 9 quantifies the predicted site in both ensembles.** Two independent
measures: the solvent-accessible surface area of the predicted residues, which
needs no cavity detector and has no threshold parameters, and the fpocket
cavity volume with its druggability score. Specificity and fold-integrity
controls come free with the same pass.

**Stage 9b writes publication figures and a statistics table** from the
per-frame CSVs.

**`--no-guided` is now the default** (`POCKET_GUIDED=false`). fpocket 4.0
advertises the `-P` guided-detection flag but does not implement it: it exits 0,
writes no output directory, and every frame is recorded as a failure — which is
indistinguishable from a pocket that never opens unless you look at the status
column. Stage 9 probes `-P` on one frame and falls back by itself, but leaving
the flag off skips the wasted probe. Nothing in the analysis depends on `-P`; it
only narrowed the search, and the overlap-plus-distance test does the real work.

**Frames are written and analysed in batches**, so peak scratch usage no longer
scales with ensemble size, and a per-batch progress line makes a failing run
obvious within a minute instead of at the end.

---

## Two Versions

| | **Local Version** | **Web Version** |
|---|---|---|
| **Launcher script** | `launch_pipeline_local.sh` | `launch_pipeline_web.sh` |
| **Apptainer def** | `pipeline_local.def` | `pipeline_web.def` |
| **Container image** | `pipeline_local.sif` (~15 GB) | `pipeline_web.sif` (~10 GB) |
| **Stage 3 engine** | Local ProtT5 (`predict_cryptic_local.py`) | CryptoBank HF Space (`scrape_cryptobank.py`) |
| **Internet needed** | Only to fetch the structure from RCSB (see note) | Yes — stage 3 calls `thorbenf-cryptobank.hf.space` |
| **Baked-in weights** | ESM-1v, ProtT5-XL + head, aSAM | ESM-1v, aSAM |

> **Note on "offline".** The local version performs cryptic-pocket *inference*
> locally, with no HuggingFace Space call. Stage 1 still downloads the target
> structure from `files.rcsb.org`. To run with no network at all, set `PDB` in
> the launcher to the path of a local `.pdb` file.

Both versions are fully self-contained. **The launchers pass no `--bind` flags,
by design** — every script and weight is read from inside the `.sif`. A bind
mount would silently mask a stale image with host files that whoever you share
the image with will not have. Please do not add them back for production runs;
the one legitimate exception is iterating on a single script during development,
and you must rebuild before sharing.

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
TOP=5                          # Number of mutations DELIVERED; a charged/polar wild type is
                               # skipped and the next-ranked cryptic residue is used instead
MODE="combined"                # "combined" = one structure; "independent" = N mutants
STRATEGY="esm"                 # "esm" = pick substitutions with ESM-Scan
MUTATION_SET="charged_polar"   # charged_polar | charged | polar | all | "ASP,GLU,..."

# --- Stage 7: clustering ---
RUN_CLUSTERING=true
N_CLUSTERS=10                  # representative conformations to extract
PCA_DIM=10                     # PCA components used for clustering

# --- Stages 6b and 9: the wild-type reference and quantification ---
RUN_WT_ENSEMBLE=true           # generate the matched wild-type ensemble. Required.
POCKET_ANALYSIS=true           # measure the predicted site in both ensembles
POCKET_MIN_OVERLAP=2           # min predicted residues lining an accepted cavity
POCKET_MAX_DIST=12.0           # max Angstrom, cavity centroid to site centroid
POCKET_GUIDED=false            # fpocket -P. Needs fpocket >= 4.2; leave false otherwise
POCKET_NPROC=16                # parallel fpocket workers
POCKET_SCRATCH=""              # optional TMPDIR for per-frame PDBs

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
├── 1JWP_A_cryptic.json              cryptic residue predictions, ALL positions
├── 1jwp_esm_summary.csv             ESM score for every charged/polar option
├── 1jwp_esm_results.json            full per-position candidate matrix
├── 1jwp_top5_esm_mutant.pdb         MODELLER mutant
├── sam2_input.pdb                   mutant, validated for aSAM
├── ensemble_output.top.pdb          MUTANT ensemble topology
├── ensemble_output.traj.dcd         MUTANT ensemble trajectory
├── wt_ensemble.top.pdb              WILD-TYPE ensemble topology
├── wt_ensemble.traj.dcd             WILD-TYPE ensemble trajectory
├── ensemble_clusters/
│   ├── clustering/cluster_representatives.dcd   one frame per cluster
│   ├── clustering/cluster_representatives.json  ensemble frame and size of each cluster
│   └── cluster_meta.sh              frame and cluster counts
├── pocket_analysis/
│   ├── pocket_report.txt            human-readable summary and the comparison
│   ├── pocket_summary.json          machine-readable, for batch screening
│   ├── mutant_pocket.csv            18 per-frame descriptors (see below)
│   ├── wildtype_pocket.csv          the same for the reference
│   ├── mutant_rmsf.csv              per-residue RMSF
│   ├── wildtype_rmsf.csv
│   ├── mutant_persistence.csv       which residues line the pocket, and how often
│   ├── wildtype_persistence.csv
│   └── pocket_volume_comparison.png
├── figures/
│   ├── fig2_pocket_opening.png/.pdf   exposure distribution + effect sizes
│   ├── figS_integrity.png/.pdf        Rg and DSSP controls
│   ├── figS_rmsf.png/.pdf             per-residue flexibility
│   ├── figS_clusters.png/.pdf         cluster populations
│   └── stats_summary.csv / .txt       every comparison with effect sizes
├── 1JWP_manifest.json               run manifest
├── 1jwp_local_report.txt            human-readable summary of the whole run
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
has a CA. Cryptic-residue prediction, ESM-Scan, MODELLER and aSAM all run on
this single chain. The full protein stays available in `<pdb>_protein.pdb` if
you want to push the whole structure to CryptoBank instead.

Selecting one chain is also what keeps aSAM correct: it flattens residues
across chains with no error, so more than one chain yields a quietly wrong
ensemble rather than a crash. Stage 5b re-runs the same validation on the
MODELLER output before the ensemble stage.

---

## Why the wild-type ensemble is mandatory

An absolute pocket volume, or an absolute solvent exposure, carries no
information. Its magnitude depends on the protein, on the generative model, and
— for cavity volume — on the detection thresholds. "The mutant's pocket is open
in 73% of frames" is uninterpretable on its own: if the wild type is also 70%,
nothing has been shown.

The only interpretable quantity is the difference between two ensembles that
differ in exactly one respect: the substitutions. That is why `RUN_WT_ENSEMBLE`
exists and why the launcher derives both runs from a single configuration block.
A difference in temperature, frame count or seed between them confounds the
comparison. Do not disable it to save GPU time; the result is not a weaker
result, it is no result.

---

## Reading the Stage 9 output

`pocket_report.txt` is laid out as a gated decision, and should be read
top-down rather than as a list of numbers.

**1. The headline.** Change in mean pocket volume over all frames (closed frames
count as zero), change in open fraction, and Cliff's δ. With thousands of
independent frames a *p*-value is driven by sample size, so report δ.

**2. Gate — specificity.** The report prints a specificity ratio. Below ~2, the
mutation loosened the whole surface and the cryptic number is not about the
cryptic site. Check this *before* believing the headline.

**3. Gate — fold integrity.** ΔRg, Δhelix and Δsheet should sit near zero. A
large shift means you measured partial unfolding, not a cavity.

**4. Cross-check — exposure.** `cryptic_sasa_A2` is computed geometrically and
never touches fpocket:

| Pocket volume | Cryptic SASA | Reading |
|---|---|---|
| ↑ | ↑ | Real opening. The strong result. |
| flat | ↑ | Widened but not enclosed — a groove, not a cavity. Thresholds may be too strict. |
| ↑ | flat | Suspect. Check the specificity ratio; a neighbouring cavity may be being assigned. |

**Status values in the CSVs.** `open` and `closed` are results. `no_output` and
`fpocket_failed` are tool failures — if you see those, the numbers are not a
null result and must not be reported as one. `plot_cryptiscan.py` refuses to
draw volume panels in that case and says why.

---

## Prerequisites

- **Linux** (x86_64)
- **NVIDIA GPU** with a standard driver (version ≥ 525)
- **Apptainer** (or Singularity)
- A free **MODELLER license key** (academic registration at
  [salilab.org](https://salilab.org/modeller/registration.html))

Nothing else. CUDA 12.5, PyTorch 2.4, ESM-1v, ProtT5, MODELLER 10.8, aSAM and
fpocket all live inside the container — no host Python, no conda environment, no
`pip install`, no weights to download.

**Disk.** Two 2000-frame ensembles plus their per-frame scratch want a few GB of
free space on whatever `TMPDIR` points at. Set `POCKET_SCRATCH` if `/tmp` is
small.

---

## Files in this Directory

| File | Purpose |
|---|---|
| `launch_pipeline_local.sh` | Launcher for the local version. |
| `launch_pipeline_web.sh` | Launcher for the web version. |
| `pipeline_local.def` | Apptainer definition for the local image. |
| `pipeline_web.def` | Apptainer definition for the web image. |
| `download_weights.sh` | Downloads the ~15 GB of weights. **Run this first.** |
| `required_files_weights/` | Everything baked into the images (see manifest below). |

### `required_files_weights/` manifest

This directory is the single source of truth for the container contents. A
truncated weight file is the most common cause of a broken build, so after copying
the folder from elsewhere (e.g. OneDrive) run `bash download_weights.sh` once: it
checks every weight against its SHA-256 and re-downloads only what is missing or
damaged.

| Path | Size |
|---|---|
| `run_pipeline.py` | 24 KB |
| `predict_cryptic_local.py` | 8.6 KB |
| `cryptobank_model_loader.py` | 27 KB |
| `scrape_cryptobank.py` | 6.9 KB |
| `esm_scanner.py` | 17 KB |
| `modeller_mutate.py` | 5.7 KB |
| `clean_structure.py` | 7.0 KB |
| `make_report.py` | 18 KB |
| `asam_cluster_analysis.py` | 4.2 KB |
| `cryptic_pocket_analysis.py` | 49 KB |
| `plot_cryptiscan.py` | 24 KB |
| `diagnose_fpocket.py` | 9.0 KB |
| `esm1v_t33_650M_UR90S_1.pt` | 7.3 GB |
| `prot_t5_xl_uniref50_full_v2/` | 6.8 GB |
| `sam2_weights/weights/mdcath_1.0/` | 1.0 GB |
| **total** | **~15 GB** |

Populate it with `bash download_weights.sh`. Quick check:

```bash
du -sh required_files_weights          # expect ~15 G
find required_files_weights -type l    # expect no output (no symlinks)
```

Script roles: `run_pipeline.py` orchestrates stages 1–5;
`predict_cryptic_local.py` runs local ProtT5 inference;
`cryptobank_model_loader.py` is the ProtT5 classification-head architecture;
`scrape_cryptobank.py` queries the HF Space; `esm_scanner.py` runs ESM-1v;
`modeller_mutate.py` does the point mutagenesis and sidechain minimization;
`clean_structure.py` does the MDAnalysis structure preparation;
`asam_cluster_analysis.py` clusters the ensemble;
`make_report.py` builds the run report and FASTA files;
`cryptic_pocket_analysis.py` is Stage 9;
`plot_cryptiscan.py` is Stage 9b;
`diagnose_fpocket.py` is a standalone troubleshooter (see below).

### Three conda environments inside each image

| Env | Path | Holds |
|---|---|---|
| main | `/opt/conda` | torch 2.4.0+cu124, ESM-1v, MODELLER 10.8, aSAM, mdtraj, pandas, scipy, matplotlib, (local only: transformers/ProtT5) |
| prep | `/opt/conda/envs/prep` | MDAnalysis 2.7.0 (and its scipy), used for structure prep and ensemble clustering |
| analysis | `/opt/conda/envs/analysis` | fpocket. A standalone binary, invoked by subprocess |

MDAnalysis and fpocket are deliberately isolated. The main env pins
`numpy==1.26.4` because `torch 2.4.0` and `mdtraj 1.10.3` require it, and
resolving either of the others against conda-forge in the main env risks pulling
numpy to 2.x and silently breaking aSAM. `run_pipeline.py` calls into the prep
env through `--prep_python`; Stage 9 reaches fpocket through `$FPOCKET_BIN`,
which the image exports. `apptainer test` asserts the main env's numpy is still
1.26.4 after the other envs are created, and logs the fpocket version.

---

## Running a single stage

Stage 9 can be re-run on ensembles you already have, without repeating the
sampling. This is the entry point to use when re-analysing a deposited ensemble
or testing threshold sensitivity:

```bash
apptainer exec pipeline_local.sif /opt/conda/bin/python3 \
  /opt/pipeline/cryptic_pocket_analysis.py \
    --mutant-top ensemble_output.top.pdb \
    --mutant-dcd ensemble_output.traj.dcd \
    --wt-top     wt_ensemble.top.pdb \
    --wt-dcd     wt_ensemble.traj.dcd \
    --residues-json 1JWP_manifest.json \
    --source-pdb 1jwp.pdb --chain A \
    --outdir pocket_analysis \
    --no-guided --min-overlap 2 --max-dist 12.0 \
    --fpocket /opt/conda/envs/analysis/bin/fpocket \
    --nproc 16
```

Add `--stride 100` for a one-minute smoke test before committing to the full
ensemble. Two arguments are easy to get wrong:

- `--source-pdb` must be the **wild-type** chain from Stage 2, not the mutant.
  The mutant's residue names no longer match the Stage 3 prediction, so the
  identity check would correctly reject it. MODELLER preserves numbering, so the
  wild-type chain maps both ensembles.
- `--residues-json` should be the **manifest**, not the raw cryptic JSON. The
  ESM walk-down means the mutated set is not in general the top *N* by score.

Then the figures:

```bash
apptainer exec pipeline_local.sif /opt/conda/bin/python3 \
  /opt/pipeline/plot_cryptiscan.py \
    --indir pocket_analysis \
    --clusters ensemble_clusters/clustering/cluster_representatives.json \
    --outdir figures --cryptic 200 205 233 235 255
```

`--cryptic` takes the residues in **aSAM numbering** — the numbers Stage 9 prints
under "Cryptic residue mapping", not the original PDB numbers.

---

## Troubleshooting

**Every frame is `no_output` or `fpocket_failed`.** A tool failure, not a closed
pocket. Run the diagnostic — it writes one frame, calls fpocket four ways, prints
exit codes and stderr for each, and gives a verdict:

```bash
apptainer exec pipeline_local.sif /opt/conda/bin/python3 \
  /opt/pipeline/diagnose_fpocket.py \
    --top ensemble_output.top.pdb --dcd ensemble_output.traj.dcd \
    --chain A --residues 200 205 233 235 255
```

The usual answer is `POCKET_GUIDED=true` on an fpocket older than 4.2. The second
is scratch space.

**`ERROR: 'fpocket' not on PATH`.** fpocket lives in its own conda env, which is
not on `PATH`. The launcher passes `--fpocket` explicitly; if you are invoking
the stage by hand, pass it too, or rely on `$FPOCKET_BIN`.

**The run used the wrong residues.** Check the "Cryptic residue mapping" block
Stage 9 prints. It maps positionally, keyed on residue number *and* insertion
code, because PDB numbering frequently has gaps (Ambler numbering in
β-lactamases) or insertion codes that an arithmetic offset would mis-resolve.
The mapping aborts on any residue-identity mismatch.

**Stale image.** The launchers pass no `--bind`, so a `.sif` built from an older
definition fails the payload check at startup. To confirm provenance:

```bash
diff <(apptainer exec pipeline_local.sif cat /opt/pipeline/run_pipeline.py) \
     required_files_weights/run_pipeline.py
```

No output means the image matches its sources.

---

## How ESM-Scan Scoring Works

Rather than mutating cryptic residues to a generic amino acid (like `GLU` or
`ALA`), **ESM-Scan** leverages Facebook Research's pre-trained **ESM-1v** protein
language model:

1. The target chain sequence is passed through ESM-1v.
2. For each position, the model computes log-probabilities for all 20 standard amino acids.
3. For any mutation $\text{WT} \to \text{Mutant}$:
   $$\text{Score} = \log P(\text{Mutant}) - \log P(\text{WT})$$
4. Candidates are restricted to `MUTATION_SET`. The restriction is the point:
   an unconstrained scan favours conservative hydrophobic swaps that preserve
   core packing, which is the opposite of the intended effect. Within a set
   chosen to be destabilising, the score still picks the variant least
   disruptive to the fold.
5. A residue whose wild type is already charged/polar is scored but not
   mutated, and the next residue down the crypticity ranking takes its place.
   `TOP=5` therefore always yields 5 mutations; the report lists every
   residue examined and why any were skipped. No substitution is ever supplied
   by default — a position with no candidate is recorded as skipped, with its
   reason, in the manifest.
6. All scores are saved to `<PDB>_esm_summary.csv` for full transparency, and the
   full per-position matrix over all 20 amino acids to `<PDB>_esm_results.json`.

---

## Downloading the Weights

The weight files are ~15 GB in total, far too large for GitHub, so they are not
in the repository. After cloning, fetch them with:

```bash
cd pipeline/
bash download_weights.sh
```

That places everything in `required_files_weights/`, which is where the `.def`
files expect it. Every file is checked against a pinned SHA-256, so a truncated
or corrupted download is caught here instead of ending up inside an image.
Re-running is safe: verified files are skipped and interrupted downloads resume.
If you only intend to build the **web** image, skip the 6.8 GB ProtT5 download
with `bash download_weights.sh --skip-prot-t5`.

| Weight | Size | Source |
|---|---|---|
| `esm1v_t33_650M_UR90S_1.pt` | 7.3 GB | `dl.fbaipublicfiles.com` (Facebook Research) |
| `prot_t5_xl_uniref50_full_v2/` | 6.8 GB | HuggingFace `ThorbenF/prot_t5_xl_uniref50_full_v2` |
| `sam2_weights/` | 1.0 GB | GitHub release `giacomo-janson/sam2` `data-1.0` |

ProtT5 is pulled file by file rather than with `git clone`, because cloning that
repository also drags in a ~6.8 GB `.git` LFS cache holding a second copy of the
same weights.

---

## Building the Container Images

Build from inside this directory, so the relative `%files` paths resolve:

```bash
cd pipeline/
sudo apptainer build pipeline_web.sif   pipeline_web.def     # ~10 GB, faster
sudo apptainer build pipeline_local.sif pipeline_local.def   # ~15 GB
chmod 0644 pipeline_local.sif pipeline_web.sif
```

Building needs internet access (base image, conda, pip, the aSAM clone) and
about 70 GB of free disk for the local image (about 45 GB for the web image).
At its peak the build holds the downloaded weights, the unpacked image, a
temporary squashfs and the finished `.sif` all at once. Point `APPTAINER_TMPDIR`
at a disk with room if `/tmp` is small (with `sudo`, pass it through using
`sudo -E`).

Verify an image before using or sharing it:

```bash
apptainer test pipeline_local.sif        # checks packages, weights, and logs the fpocket version
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

---

## Citing

If CryptiScan is useful in your work, please cite the resource paper and the
underlying mutagenesis protocol, together with the models the pipeline calls:
CryptoBank, ESM-1v, aSAM and mdCATH, MODELLER, fpocket, MDAnalysis and MDTraj.
Full references are in the manuscript's bibliography.
