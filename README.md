# CryptiScan

Find a cryptic binding pocket from sequence, open it, and measure how far it
opened. Everything runs in one Apptainer image on one GPU.

**You edit two lines: the PDB accession and the chain.** No residue numbers,
no structural alignment, no manual site definition at any point.

---

## The pipeline

```
  ┌────────────────────────────────────────────────────────────┐
  │                     PDB ID   +   CHAIN                     │
  └────────────────────────────────────────────────────────────┘
                                 │   the only two inputs
                                 ▼
  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
   A   DISCOVERY                                                
  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  ┌────────────────────────────────────────────────────────────┐
  │  1  Target ingestion            RCSB API / local .pdb      │
  ├────────────────────────────────────────────────────────────┤
  │  2  Structure curation          MDAnalysis 2.7             │
  ├────────────────────────────────────────────────────────────┤
  │  3  Chain selection             MDAnalysis 2.7             │
  ├────────────────────────────────────────────────────────────┤
  │  4  Crypticity scoring          CryptoBank ProtT5          │
  └────────────────────────────────────────────────────────────┘
                                 │
                                 ▼
  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
   B   PERTURBATION  AND  GENERATIVE SAMPLING                   
  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  ┌────────────────────────────────────────────────────────────┐
  │  5  Variant selection           ESM-1v (650M)              │
  ├────────────────────────────────────────────────────────────┤
  │  6  Mutagenesis                 MODELLER 10.8              │
  ├────────────────────────────────────────────────────────────┤
  │  7  Mutant ensemble             aSAM  2000 @ 300 K         │
  ├────────────────────────────────────────────────────────────┤
  │  8  Wild-type ensemble          aSAM  2000 @ 300 K         │
  └────────────────────────────────────────────────────────────┘
                                 │
                                 ▼
  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
   C   QUANTIFICATION,  CONTROLS  AND  REPORTING                
  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  ┌────────────────────────────────────────────────────────────┐
  │  9  Clustering                  PCA + k-means              │
  ├────────────────────────────────────────────────────────────┤
  │ 10  Pocket detection            fpocket  -P  -u 2          │
  ├────────────────────────────────────────────────────────────┤
  │ 11  Statistics, figures         Cliff's d, KS, MWU         │
  ├────────────────────────────────────────────────────────────┤
  │ 12  Report and archive          <job>.zip                  │
  └────────────────────────────────────────────────────────────┘
                                 │
                                 ▼
  ┌────────────────────────────────────────────────────────────┐
  │             EFFECT SIZE  vs  MATCHED WILD TYPE             │
  └────────────────────────────────────────────────────────────┘
```

Stage 7 and 8 are the GPU cost. Stage 10 parallelises across CPU cores.

---

## Quick start

```bash
cd pipeline/

bash download_weights.sh                                    # ~15 GB, once
sudo apptainer build pipeline_local.sif pipeline_local.def  # ~2 h

export MODELLER_KEY='your_academic_key'
# edit PDB and CHAIN at the top of launch_pipeline_local.sh
nohup bash launch_pipeline_local.sh > run.log 2>&1 &
tail -f run.log
```

A complete worked example is in [`sample_output/`](sample_output/).

---

## Configuration

Section 1 of `launch_pipeline_local.sh` is the whole user interface:

```bash
PDB="2WGB"      # PDB accession, or a path to a local .pdb
CHAIN="A"       # chain to analyse
```

Section 2 holds the values used in the paper. Change them only deliberately,
and report what you changed.

| Setting | Default | What it does |
|---|---|---|
| `TOP` | 5 | mutations delivered. A position whose wild type is already charged or polar is skipped and the next-ranked one is used, so you always get this many. |
| `MUTATION_SET` | `charged_polar` | the candidate set. An unrestricted scan picks conservative hydrophobic swaps and leaves the pocket shut. |
| `FRAMES` | 2000 | conformations per ensemble. 500 reproduces the effect sizes at a quarter of the cost, with weaker bounds on the controls. |
| `TEMPERATURE` | 300 | kelvin |
| `POCKET_MIN_EXPLICIT` | 2 | fpocket `-u`. See below. |
| `POCKET_SPECIFICITY` | true | run the unguided pass that gives the specificity control. Doubles Stage 10. |
| `POCKET_SCRATCH` | `$TMPDIR` | **point this at local disk.** Stage 10 writes one file per frame, and on network storage that dominates the whole runtime. |

`MODELLER_KEY` and `JOBNAME` can be set from the environment instead:

```bash
JOBNAME=2wgb_rep2 MODELLER_KEY=your_key bash launch_pipeline_local.sh
```

---

## Two settings worth understanding

**`POCKET_MIN_EXPLICIT`, fpocket's `-u`.** An alpha sphere is kept only if at
least this many of its four contact atoms belong to the cryptic residues.
fpocket's own default is 4, meaning all four, which for a five-residue site
keeps almost nothing and reports every frame as closed. That looks like a
negative result, not a parameter error. Use one value across everything you
intend to compare, and state it. `POCKET_USWEEP=1` writes the sensitivity
table.

**`RUN_WT_ENSEMBLE`, default true.** An absolute pocket volume carries no
information: "open in 95% of frames" means nothing until you know the wild
type is 54%. Both ensembles come from one configuration block because a
difference in temperature, frame count or seed between them confounds the
comparison. Turning this off does not give a weaker result, it gives no
result.

---

## What you get

```
<jobname>_<timestamp>/
├── 2wgb.pdb                      curated chain, used by every later stage
├── 2wgb_esm_summary.csv          score for every candidate substitution
├── 2wgb_top5_esm_mutant.pdb      the mutant
├── ensemble_output.{top.pdb,traj.dcd}   MUTANT ensemble
├── wt_ensemble.{top.pdb,traj.dcd}       WILD-TYPE ensemble
├── ensemble_clusters/            10 representatives, weighted, for docking
├── pocket_analysis/
│   ├── pocket_report.txt         read this first
│   ├── pocket_summary.json       machine-readable
│   ├── {mutant,wildtype}_pocket.csv        per-frame descriptors
│   ├── {mutant,wildtype}_rmsf.csv          per-residue flexibility
│   └── {mutant,wildtype}_persistence.csv   which residues line the pocket
├── figures/
│   ├── fig2_opening.{png,pdf}        exposure, volume, effect sizes
│   ├── figS1_integrity.{png,pdf}     Rg and secondary structure
│   ├── figS2_rmsf.{png,pdf}
│   ├── figS3_persistence.{png,pdf}
│   ├── figS4_clusters.{png,pdf}
│   └── stats_summary.{csv,txt}
├── 2WGB_manifest.json            run manifest, mutations applied
└── 2wgb_local_report.txt
```

---

## Reading `pocket_report.txt`

It is a gated decision. Read it top down and stop at the first gate that fails.

1. **Headline.** Change in mean pocket volume, change in open fraction, and
   Cliff's delta. With thousands of independent frames a *p*-value is driven by
   sample size, so delta is the quantity to report.
2. **Gate: specificity.** Below a ratio of about 2, the mutation loosened the
   whole surface and the cryptic number is not about the cryptic site.
3. **Gate: fold integrity.** Delta Rg, helix and sheet should sit near zero. A
   large shift means you measured partial unfolding, not a cavity.
4. **Cross-check: exposure.** `cryptic_sasa_A2` is geometric and never touches
   fpocket. Volume and SASA both up is the strong result. SASA up with flat
   volume is a groove, not a cavity.

`open` and `closed` are results. `no_output` and `fpocket_failed` are tool
failures, not null results, and the figure stage refuses to draw volume panels
when they are present.

---

## Requirements

- Linux x86_64, NVIDIA GPU with driver 525 or later, Apptainer or Singularity
- A free academic MODELLER key from [salilab.org](https://salilab.org/modeller/registration.html)
- About 70 GB free to build, a few GB of local scratch per run

CUDA 12.4, PyTorch 2.4, ESM-1v, ProtT5, MODELLER 10.8, aSAM and fpocket are all
inside the image. No host Python, no conda environment, no `pip install`.

Verify an image before using or sharing it:

```bash
apptainer test pipeline_local.sif
apptainer inspect --deffile pipeline_local.sif | diff - pipeline_local.def
```

The second command proves the image was built from the definition beside it.
The launcher passes no bind mounts by design, so a stale image fails the
payload check in seconds rather than three hours in.

---

## Re-analysing an ensemble you already have

Stages 10 and 11 run on existing ensembles without repeating the sampling:

```bash
apptainer exec pipeline_local.sif /opt/conda/bin/python3 \
  /opt/pipeline/cryptic_pocket_analysis.py \
    --mutant-top ensemble_output.top.pdb --mutant-dcd ensemble_output.traj.dcd \
    --wt-top wt_ensemble.top.pdb --wt-dcd wt_ensemble.traj.dcd \
    --residues-json <PDB>_manifest.json --source-pdb <pdb>.pdb --chain A \
    --outdir pocket_analysis --min-explicit 2 \
    --fpocket "$FPOCKET_BIN" --nproc "$(nproc)"
```

Pass the **wild-type** chain as `--source-pdb`, not the mutant, and the
**manifest** as `--residues-json`, not the raw crypticity output. Add
`--stride 100` for a one-minute check first.

---

## Files

| File | Purpose |
|---|---|
| `launch_pipeline_local.sh` | the launcher, the only file you edit |
| `pipeline_local.def` | Apptainer definition |
| `download_weights.sh` | fetches the weights, checksum-verified. Run first. |
| `required_files_weights/` | everything baked into the image |
| `sample_output/` | a complete 2WGB run |

---

## Citing

Please cite the resource paper and the underlying mutagenesis protocol, with
the models the pipeline calls: CryptoBank, ESM-1v, aSAM and mdCATH, MODELLER,
fpocket, MDAnalysis and MDTraj. Full references are in the manuscript.
