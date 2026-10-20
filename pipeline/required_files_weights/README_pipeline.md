# CryptiScan v1.1 — what changed and how to install it

## Install

```bash
cd /data/Routh/DeepCrypt/Nitin/CryptiScan/Refined/pipeline

cp launch_pipeline_local.sh launch_pipeline_local.sh.bak
cp pipeline_local.def       pipeline_local.def.bak

# replaced
cp <downloads>/cryptic_pocket_analysis.py required_files_weights/
cp <downloads>/make_report.py             required_files_weights/
cp <downloads>/launch_pipeline_local.sh   .
cp <downloads>/pipeline_local.def         .

# new
cp <downloads>/plot_cryptiscan.py         required_files_weights/
cp <downloads>/check_stage10.py           required_files_weights/

sudo apptainer build pipeline_local.sif pipeline_local.def
apptainer test pipeline_local.sif
```

`apptainer test` now fails the build if any of the six Stage 10 regressions
below is reintroduced.

## Run

Edit two lines at the top of `launch_pipeline_local.sh`:

```bash
PDB="2WGB"
CHAIN="A"
```

Then `bash launch_pipeline_local.sh`. Nothing else needs touching, and nothing
asks for residue numbers at any point.

## Six bugs fixed in Stage 10

All six made guided detection return nothing while looking like a genuine
negative result. Each is now asserted at build time by `check_stage10.py`.

1. **`-P` spec format.** The code built `127::A`. fpocket splits on `:` with
   `strtok_r`, which collapses the empty field, and the switch consuming the
   pieces has no `break` statements, so the insertion code ended up holding the
   chain letter and the match failed for every residue. Must be `127:-:A`.

2. **`-u` was never passed.** An alpha sphere is kept only if at least
   `min_n_explicit_pocket_atoms` of its four contact atoms belong to the listed
   residues. fpocket's default is 4, i.e. all four, which for a five-residue
   site retains almost nothing. Now 2, exposed as `POCKET_MIN_EXPLICIT`.

3. **Off-by-one in the pocket filename.** `info.txt` and the files on disk both
   number from 1: `Pocket 1 :` is `pockets/pocket1_atm.pdb`. The code looked for
   `pocket0_atm.pdb`, so every pocket was skipped and every frame was reported
   closed — even after fixes 1 and 2.

4. **The temp path contained a space.** The per-ensemble scratch directory was
   named from the label, so the wild type got `fp_wild type_XXXX/`. fpocket
   builds output paths with C string operations and cannot write there: it
   prints "could not be opened" and leaves no `info.txt`. **Every wild-type
   frame failed, every run**, while the mutant succeeded. That reads as a
   one-sided biological result rather than a path bug.

5. **The specificity control was always zero.** It summed "every other pocket"
   over the *guided* run, which by construction contains only the cryptic
   pocket. A separate unguided pass per frame now supplies the real background.

6. **Backbone RMSD used each ensemble's own frame 0.** That makes the two
   numbers incomparable. Both now use one common reference; the old value is
   kept as `rmsd_ownref_A`.

## Other changes

**Stage numbering.** Twelve stages, numbered consistently in the launcher, the
report, the `.def` comments and the new flowchart. Wild-type ensemble, pocket
quantification and figures are now numbered stages rather than bolt-ons.

**`TEMPERATURE` and `FRAMES` are variables.** The old launcher hard-coded
`-T 300` and `-n 1000` in two places each. **The manuscript says 320 K and the
2000-frame runs you analysed; the launcher said 300 K.** Defaults are now 320
and 2000, in one place, used by both ensembles. Confirm which temperature the
2WGB run actually used before the Methods section goes out.

**`POCKET_SPECIFICITY`** toggles the second unguided pass. On by default; off
halves Stage 10 runtime and loses the control.

**`POCKET_SCRATCH`** points per-frame temporaries at local disk. On network
storage this dominates the runtime.

**Stage 11 figures** are new to the container: `plot_cryptiscan.py` writes
`fig2_opening`, `figS1_integrity`, `figS2_rmsf`, `figS3_persistence`,
`figS4_clusters` and `stats_summary`. Filenames match the manuscript's
numbering. It reads the cryptic residues from `pocket_summary.json` itself.

**Stage 10 in the report.** `make_report.py` gains section 7: open fractions,
volumes, the specificity control, effect sizes, and the fpocket build and `-u`
value that produced them.

**Provenance.** `pocket_summary.json` and the report now record the fpocket
version string, the `-u` value and which RMSD reference was used.

## Figure 1

`flowchart.pdf` and `.png` replace the old pipeline image; `make_flowchart.py`
regenerates them if you want to edit wording. Put the PDF in `Figs/`.
