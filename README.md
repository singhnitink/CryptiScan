# CryptiScan

CryptiScan predicts cryptic binding pockets from sequence, destabilises the
closed state with a small number of selected substitutions, generates matched
mutant and wild-type conformational ensembles, and measures whether the
predicted site actually opens.

Every stage runs inside one Apptainer image. The host needs Linux, an NVIDIA
driver, and Apptainer. No host Python, no conda environment, no `pip install`.

**You edit two lines: the PDB accession and the chain.** Nothing asks for
residue numbers at any point.

---

## Quick start

```bash
cd pipeline/

bash download_weights.sh                                     # ~15 GB, once
sudo apptainer build pipeline_local.sif pipeline_local.def    # ~2 h

export MODELLER_KEY='your_academic_key'
# edit PDB and CHAIN at the top of launch_pipeline_local.sh
nohup bash launch_pipeline_local.sh > run.log 2>&1 &
tail -f run.log
```

A worked example with its inputs and expected outputs is in
[`sample_output/`](sample_output/).

---

## The twelve stages

```
                      PDB ID  or  local .pdb
                                  |
  A. TARGET ACQUISITION, CURATION AND CRYPTIC SITE DISCOVERY
  ---------------------------------------------------------------------
   1  Target ingestion          RCSB API / local file
   2  Structure curation        MDAnalysis 2.7     waters dropped,
                                                   ligands split out,
                                                   altlocs collapsed
   3  Chain selection           standard residues and backbone verified
   4  Cryptic pocket scoring    CryptoBank ProtT5  per-residue ranking
                                  |
  B. VARIANT SELECTION, MUTAGENESIS AND GENERATIVE SAMPLING
  ---------------------------------------------------------------------
   5  ESM-scan selection        ESM-1v 650M        charged/polar set
   6  In silico mutagenesis     MODELLER 10.8
   7  Mutant ensemble           aSAM (mdCATH)      2000 states, 300 K
   8  Wild-type ensemble        aSAM (mdCATH)      identical settings
                                  |
  C. QUANTIFICATION, CONTROLS AND REPORTING
  ---------------------------------------------------------------------
   9  Ensemble clustering       PCA + k-means      10 representatives
  10  Pocket quantification     fpocket -P -u 2    per-frame volume, SASA
  11  Statistics and figures    Cliff's delta, KS, MWU
  12  Report and archive        <job>_<timestamp>.zip
```

Stages 7 and 8 are the GPU cost. Stage 10 parallelises across cores and is
the only stage that scales with your core count.

---

## Configuration

Open `launch_pipeline_local.sh`. Section 1 is the whole user interface:

```bash
PDB="2WGB"                     # PDB accession, or a path to a local .pdb
CHAIN="A"                      # chain to analyse
MODELLER_KEY="${MODELLER_KEY:-MODELIRANJE}"
```

Section 2 holds the values used in the paper. Change them only if you intend
to, and report what you changed:

```bash
TOP=5                 # cryptic residues mutated. A position whose wild type is
                      # already charged or polar is scored but skipped, and the
                      # scan walks further down the ranking, so TOP mutations
                      # are always delivered.
MUTATION_SET="charged_polar"
FRAMES=2000           # per ensemble
TEMPERATURE=300       # kelvin
POCKET_MIN_EXPLICIT=2 # fpocket -u, see below
POCKET_SPECIFICITY=true
POCKET_NPROC=$(nproc)
POCKET_SCRATCH="${TMPDIR:-/tmp}"
```

`MODELLER_KEY` and `JOBNAME` can be set from the environment instead of edited
into the file:

```bash
JOBNAME=2wgb_rep2 MODELLER_KEY=your_key bash launch_pipeline_local.sh
```

Point `POCKET_SCRATCH` at local disk. Stage 10 writes one PDB per frame, and
on network storage that dominates the runtime of the whole pipeline.

---

## Why residue numbers are never an input

The site comes from the Stage 4 crypticity ranking, and Stage 5 records which
positions were actually mutated and which were skipped. Stage 10 reads that
manifest and maps the positions into aSAM's topology by walking the chain in
file order, keyed on residue number and insertion code together.

This matters because aSAM renumbers residues from zero. Any arithmetic offset
between the deposited numbering and the ensemble numbering breaks on PDB
numbering gaps and on insertion codes, and it breaks silently: you get a volume
for the wrong five residues rather than an error. The mapping aborts on any
residue-identity mismatch instead, and prints the correspondence it used.

---

## Two parameters that will cost you a run if you change them carelessly

**`POCKET_MIN_EXPLICIT` (fpocket `-u`), default 2.** An alpha sphere is kept
only if at least this many of its four contact atoms belong to the listed
cryptic residues. fpocket's own default is 4, meaning all four, which for a
site of five residues retains essentially nothing and reports every frame as
`No Pockets Found while refining`. That reads as a closed pocket, not as a
parameter error. The value changes absolute volumes, so use one setting across
everything you intend to compare, and state it. The sensitivity sweep is
available with `POCKET_USWEEP=1`.

**`RUN_WT_ENSEMBLE`, default true.** An absolute pocket volume carries no
information. Its magnitude depends on the protein, on the generative model and
on the detection parameters. "Open in 94% of frames" means nothing until you
know the wild type is 54%. The only interpretable quantity is the difference
between two ensembles differing in exactly one respect, which is why both runs
come from one configuration block: a difference in temperature, frame count or
seed between them confounds the comparison. Disabling this does not give a
weaker result, it gives no result.

---

## A note on guided detection

Guided detection with fpocket's `-P` flag works, and is what the pipeline uses.
If you are coming from an earlier copy of this folder that set
`POCKET_GUIDED=false` and said the flag was unimplemented, that was wrong. The
flag was being handed a malformed residue specification.

`-P` takes `resnum:icode:chain`. fpocket splits each field with `strtok_r`,
which collapses an empty field, and the switch consuming the pieces falls
through. Passing `127::A` for a residue with no insertion code therefore leaves
the chain letter sitting in the insertion-code slot, and the match fails for
every residue in every frame. The correct form is `127:-:A`, with an explicit
placeholder.

Five further defects had the same signature, each of which independently
returned no pocket for all 4000 frames. They are asserted against at build time
by `check_stage10.py`, which the `%test` block runs: the image fails to build
rather than let one of them return silently. If you modify
`cryptic_pocket_analysis.py`, run that script against your copy before
trusting the output.

---

## Output

```
<jobname>_<timestamp>/
├── 2wgb_raw.pdb                  entry as downloaded
├── 2wgb_protein.pdb              full protein, all chains
├── 2wgb_ligands/                 one file per ligand
├── 2wgb_prep.json                what was kept, split and dropped
├── 2wgb.pdb                      SELECTED CHAIN, used by every later stage
├── 2WGB_A_cryptic.json           crypticity, all positions
├── 2wgb_esm_summary.csv          ESM score for every candidate
├── 2wgb_esm_results.json         full per-position matrix, all 20 residues
├── 2wgb_top5_esm_mutant.pdb      MODELLER mutant
├── sam2_input.pdb                mutant, validated for aSAM
├── ensemble_output.top.pdb       MUTANT ensemble
├── ensemble_output.traj.dcd
├── wt_ensemble.top.pdb           WILD-TYPE ensemble
├── wt_ensemble.traj.dcd
├── ensemble_clusters/
│   └── clustering/cluster_representatives.{dcd,json}
├── pocket_analysis/
│   ├── pocket_report.txt         read this first, it is a gated decision
│   ├── pocket_summary.json       machine-readable, for batch screening
│   ├── mutant_pocket.csv         per-frame descriptors
│   ├── wildtype_pocket.csv
│   ├── mutant_rmsf.csv           per-residue flexibility
│   ├── wildtype_rmsf.csv
│   ├── mutant_persistence.csv    which residues line the pocket, how often
│   └── wildtype_persistence.csv
├── figures/
│   ├── fig2_opening.{png,pdf}        exposure, volume, ECDF, effect sizes
│   ├── figS1_integrity.{png,pdf}     Rg and secondary-structure controls
│   ├── figS2_rmsf.{png,pdf}
│   ├── figS3_persistence.{png,pdf}
│   ├── figS4_clusters.{png,pdf}
│   └── stats_summary.{csv,txt}
├── 2WGB_manifest.json            run manifest, including mutations applied
├── 2wgb_local_report.txt
├── 2WGB_A_wildtype.fasta
└── 2WGB_A_mutant.fasta
```

`pocket_summary.json` records `fpocket_version`,
`fpocket_min_explicit_u` and `rmsd_reference`, so a result can be traced to the
parameters that produced it.

---

## Reading `pocket_report.txt`

It is laid out as a gated decision. Read it top down.

**1. Headline.** Change in mean pocket volume over all frames, closed frames
counted as zero; change in open fraction; Cliff's delta. With thousands of
independent frames a *p*-value is driven by sample size, so delta is the
quantity to report.

**2. Gate, specificity.** Below a ratio of about 2, the mutation loosened the
whole surface and the cryptic number is not about the cryptic site. Check this
before believing the headline. It comes from a separate unguided pass, which is
why `POCKET_SPECIFICITY` doubles the Stage 10 runtime.

**3. Gate, fold integrity.** Delta Rg, delta helix and delta sheet should sit
near zero. A large shift means you measured partial unfolding, not a cavity.

**4. Cross-check, exposure.** `cryptic_sasa_A2` is geometric and never touches
fpocket.

| Pocket volume | Cryptic SASA | Reading |
|---|---|---|
| up | up | Real opening. The strong result. |
| flat | up | Widened but not enclosed. A groove, not a cavity. |
| up | flat | Suspect. Check specificity; a neighbouring cavity may be assigned. |

**Status values.** `open` and `closed` are results. `no_output` and
`fpocket_failed` are tool failures. If you see those, the numbers are not a
null result and must not be reported as one. `plot_cryptiscan.py` refuses to
draw volume panels in that case and says why.

---

## Prerequisites

- Linux x86_64
- NVIDIA GPU, driver >= 525
- Apptainer or Singularity
- A free academic MODELLER key from
  [salilab.org](https://salilab.org/modeller/registration.html). Licenses are
  issued per user; register for your own rather than reusing someone else's.

Disk: about 70 GB free to build the image, a few GB of scratch per run.

CUDA 12.4, PyTorch 2.4, ESM-1v, ProtT5, MODELLER 10.8, aSAM and fpocket are all
inside the image.

---

## Files in `pipeline/`

| File | Purpose |
|---|---|
| `launch_pipeline_local.sh` | The launcher. The only file you edit. |
| `pipeline_local.def` | Apptainer definition for the local image. |
| `download_weights.sh` | Fetches the ~15 GB of weights. Run first. |
| `required_files_weights/` | Everything baked into the image. |
| `sample_output/` | A complete 2WGB run: inputs, key outputs, figures. |

### `required_files_weights/` manifest

| Path | Size | Role |
|---|---|---|
| `run_pipeline.py` | 24 KB | orchestrates Stages 1 to 6 |
| `clean_structure.py` | 7.0 KB | MDAnalysis curation and chain selection |
| `predict_cryptic_local.py` | 8.6 KB | local ProtT5 crypticity inference |
| `cryptobank_model_loader.py` | 27 KB | ProtT5 classification-head architecture |
| `esm_scanner.py` | 17 KB | ESM-1v substitution scoring |
| `modeller_mutate.py` | 5.7 KB | point mutagenesis and sidechain relaxation |
| `asam_cluster_analysis.py` | 4.2 KB | Stage 9 clustering |
| `cryptic_pocket_analysis.py` | 49 KB | Stage 10 quantification |
| `plot_cryptiscan.py` | 32 KB | Stage 11 figures and statistics |
| `check_stage10.py` | 3.7 KB | build-time regression guard |
| `make_report.py` | 22 KB | Stage 12 report and FASTA |
| `esm1v_t33_650M_UR90S_1.pt` | 7.3 GB | ESM-1v weights |
| `prot_t5_xl_uniref50_full_v2/` | 6.8 GB | ProtT5 weights plus head |
| `sam2_weights/weights/mdcath_1.0/` | 1.0 GB | aSAM weights |

A truncated weight file is the commonest cause of a broken build.
`download_weights.sh` checks every file against a pinned SHA-256 and
re-downloads only what is missing or damaged, so re-running it is safe. Use
`--skip-prot-t5` if you only need the web image.

```bash
du -sh required_files_weights          # expect ~15 G
find required_files_weights -type l    # expect no output
```

### Three conda environments inside the image

| Env | Path | Holds |
|---|---|---|
| main | `/opt/conda` | torch 2.4.0+cu124, ESM-1v, MODELLER 10.8, aSAM, mdtraj, transformers, matplotlib |
| prep | `/opt/conda/envs/prep` | MDAnalysis 2.7.0, used for curation and clustering |
| analysis | `/opt/conda/envs/analysis` | fpocket, a standalone binary invoked by subprocess |

MDAnalysis and fpocket are isolated deliberately. The main environment pins
`numpy==1.26.4` because torch 2.4.0 and mdtraj 1.10.3 both require it, and
resolving either of the others in the main environment risks pulling numpy to
2.x, which breaks aSAM quietly rather than loudly. `apptainer test` asserts the
main environment's numpy after the other two are created.

`pyarrow` is pinned to 17.0.0 for the same reason and is the harder case: it
declares no numpy constraint in its wheel metadata, so pip cannot catch the
conflict and nothing warns at install time. From pyarrow 18 onward the compiled
extension asserts numpy >= 2.0 at import, and `import datasets` fails during the
build with an error that names numpy. Do not relax that pin without rebuilding
and running `apptainer test`.

---

## Building

Build from inside `pipeline/`, so the relative `%files` paths resolve:

```bash
cd pipeline/
sudo apptainer build pipeline_local.sif pipeline_local.def
chmod 0644 pipeline_local.sif
```

Building needs internet access and about 70 GB of free disk. At its peak the
build holds the weights, the unpacked image, a temporary squashfs and the
finished `.sif` at once. If `/tmp` is small, point `APPTAINER_TMPDIR` at a
larger disk, passing it through with `sudo -E`.

Verify before using or sharing an image:

```bash
apptainer test pipeline_local.sif
apptainer run-help pipeline_local.sif
apptainer inspect --deffile pipeline_local.sif | diff - pipeline_local.def
```

`apptainer test` checks every package and weight, asserts the numpy and pyarrow
pins, logs the fpocket version, and runs the Stage 10 guards. The last command
proves the `.sif` was built from the `.def` sitting next to it.

The launcher passes no `--bind` flags, by design. Every script and weight is
read from inside the `.sif`. A bind mount would mask a stale image with host
files that whoever you share the image with will not have, and the launcher's
own payload check would pass while the run used different code. The one
legitimate exception is iterating on a single script during development, and
you must rebuild before sharing.

---

## Re-analysing ensembles you already have

Stage 10 and Stage 11 can be re-run on existing ensembles without repeating the
sampling. This is the entry point for re-analysing a deposited ensemble or
testing parameter sensitivity:

```bash
apptainer exec pipeline_local.sif /opt/conda/bin/python3 \
  /opt/pipeline/cryptic_pocket_analysis.py \
    --mutant-top ensemble_output.top.pdb \
    --mutant-dcd ensemble_output.traj.dcd \
    --wt-top     wt_ensemble.top.pdb \
    --wt-dcd     wt_ensemble.traj.dcd \
    --residues-json 2WGB_manifest.json \
    --source-pdb 2wgb.pdb --chain A \
    --outdir pocket_analysis \
    --min-explicit 2 \
    --fpocket "$FPOCKET_BIN" \
    --nproc "$(nproc)"
```

Add `--stride 100` for a one-minute smoke test before committing to the full
ensemble. Two arguments are easy to get wrong:

- `--source-pdb` must be the **wild-type** chain from Stage 3, not the mutant.
  The mutant's residue names no longer match the Stage 4 prediction, so the
  identity check would correctly reject it. MODELLER preserves numbering, so
  the wild-type chain maps both ensembles.
- `--residues-json` must be the **manifest**, not the raw cryptic JSON. The ESM
  walk-down means the mutated set is not in general the top *N* by score.

Then the figures, which read the cryptic residues out of
`pocket_summary.json` and need no residue argument:

```bash
apptainer exec pipeline_local.sif /opt/conda/bin/python3 \
  /opt/pipeline/plot_cryptiscan.py \
    --indir pocket_analysis \
    --clusters ensemble_clusters/clustering/cluster_representatives.json \
    --outdir figures --label 2WGB
```

---

## Troubleshooting

**Every frame is `no_output` or `fpocket_failed`.** A tool failure, not a closed
pocket. First confirm the guards are intact:

```bash
apptainer exec pipeline_local.sif /opt/conda/bin/python3 \
  /opt/pipeline/check_stage10.py
```

If those pass, the usual remaining cause is scratch space. Check that
`POCKET_SCRATCH` points at local disk with room for one PDB per worker.

**Only wild-type frames fail.** A path-sanitisation problem, which looks like a
one-sided biological result. fpocket cannot write into a directory whose name
contains a space, and the ensemble label is "wild type". The fix is in the
shipped code; if you see this, you are running a patched copy.

**`'fpocket' not on PATH`.** fpocket lives in its own environment, which is not
on `PATH`. The launcher passes `--fpocket` explicitly. If you are invoking the
stage by hand, pass it too, or use `$FPOCKET_BIN`, which the image exports.

**The run used the wrong residues.** Read the cryptic residue mapping block
Stage 10 prints. It maps positionally on residue number and insertion code and
aborts on any identity mismatch, so a wrong mapping should not be reachable.
Report it if it is.

**Stale image.** The launcher checks the payload before Stage 1 and stops in
seconds rather than three hours in. To confirm provenance:

```bash
diff <(apptainer exec pipeline_local.sif cat /opt/pipeline/cryptic_pocket_analysis.py) \
     required_files_weights/cryptic_pocket_analysis.py
```

No output means the image matches its sources.

**The build fails at `import datasets` with a numpy error.** The pyarrow pin was
relaxed. See the environments section above.

---

## How ESM-scan selection works

1. The target chain sequence is passed through ESM-1v.
2. At each position the model gives log-probabilities for all 20 amino acids.
3. For a substitution WT -> Mut, the score is
   `log P(Mut) - log P(WT)`.
4. Candidates are restricted to `MUTATION_SET`, by default charged and polar.
   The restriction is the point: an unconstrained scan favours conservative
   hydrophobic swaps that preserve core packing, which is the opposite of the
   intended effect. Within a set chosen to be destabilising, the score still
   picks the variant least disruptive to the fold.
5. A position whose wild type is already charged or polar is scored but not
   mutated, and the next position down the crypticity ranking takes its place.
   `TOP=5` therefore always yields five mutations. No substitution is ever
   supplied by default: a position with no candidate is recorded as skipped,
   with its reason, in the manifest.
6. Every score goes to `<pdb>_esm_summary.csv`, and the full per-position matrix
   over all 20 amino acids to `<pdb>_esm_results.json`.

---

## Runtime

Dominated by Stages 7 and 8, which are GPU-bound, and Stage 10, which
parallelises across cores. Stage 10 with `POCKET_SPECIFICITY=true` makes two
passes over 4000 frames.

Measured timings for the 2WGB benchmark, with the hardware they were measured
on, are in `sample_output/timings.txt`.

---

## Citing

If CryptiScan is useful in your work, please cite the resource paper and the
underlying mutagenesis protocol, together with the models the pipeline calls:
CryptoBank, ESM-1v, aSAM and mdCATH, MODELLER, fpocket, MDAnalysis and MDTraj.
Full references are in the manuscript bibliography.
