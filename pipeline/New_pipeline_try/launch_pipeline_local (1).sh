#!/bin/bash

# ==============================================================================
# CryptiScan Pipeline Launcher (Local / Offline Version)
# ==============================================================================
# Workflow:
#   1.  Fetch the entry from RCSB and clean it with MDAnalysis. Keeps the full
#       protein (<pdb>_protein.pdb, all chains) and each ligand separately
#       under <pdb>_ligands/. Waters dropped, alt-locs collapsed.
#   2.  Select chain $CHAIN -> <pdb>.pdb. Every stage below uses this chain.
#   3.  Predict cryptic residues locally with ProtT5 (no HuggingFace call)
#   4.  Run ESM-Scan (ESM-1v) to find the best-scoring mutation per residue
#   5.  Run MODELLER to build the relaxed 3D mutant structure
#   5b. Validate the mutant for SAM2 -> sam2_input.pdb
#   6.  Run SAM2 / aSAM ensemble generation on that chain
#   7.  Write <jobname>_report.txt plus wild-type / mutant FASTA sequences
#   8.  Compress all results into a .zip archive
#
# Requirements on the host PC:
#   - Linux (Ubuntu, Debian, RHEL, Arch, ...)
#   - NVIDIA GPU with standard drivers (driver version >= 525)
#   - Apptainer (or Singularity) installed
#   * No host Python, CUDA toolkit, or scientific packages needed.
#   * Cryptic-pocket inference is local - no HuggingFace Space is contacted.
#     Stage 2 still fetches the structure from files.rcsb.org, so set PDB to a
#     local .pdb path below if this machine has no internet access.
#
# Everything else lives inside pipeline_local.sif. This script passes NO --bind flags
# on purpose: a bind mount would silently mask a stale image with host files.
#
# For the lighter web-scraping version, use: bash launch_pipeline_web.sh
# ==============================================================================

# Automatically detect the directory where this script lives
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# ------------------------------------------------------------------------------
# 1. User Configuration (Edit these for your target protein)
# ------------------------------------------------------------------------------
PDB="1JWP"                     # PDB accession code (4 characters) or local PDB path
CHAIN="A"                      # Target protein chain identifier
TOP=5                          # Number of top cryptic residues to mutate
MODE="combined"                # "combined" = all mutations in 1 structure; "independent" = N separate mutants
STRATEGY="esm"                 # "esm" = evaluates substitutions via ESM-Scan protein language model
MUTATION_SET="charged_polar"   # Which residues ESM may choose from:
                               #   charged_polar (default) | charged | polar | all | "ASP,GLU,..."
                               # Charged/polar substitutions destabilise the closed state so the
                               # cryptic pocket can open in MD. "all" reverts to the unrestricted
                               # 19-way scan, which tends to pick conservative nonpolar swaps
                               # (ILE->VAL) that leave the pocket shut.
JOBNAME="${JOBNAME:-1jwp}_local"     # Job name for folder and zip naming (e.g. "1lzt", "1jwp")

# === STAGE 6: ensemble clustering ===
RUN_CLUSTERING=true            # PCA + k-means on the mutant ensemble, one
                               # representative frame per cluster. make_report.py
                               # reads its cluster_meta.sh, so leaving this off
                               # makes the report say "Clusters: not run".
N_CLUSTERS=10                  # representative conformations to extract
PCA_DIM=10                     # principal components used for clustering

# === STAGE 9: cryptic pocket quantification ===
RUN_WT_ENSEMBLE=true           # also generate a WILD-TYPE ensemble. Required for
                               # a meaningful result: a mutant pocket volume with
                               # no reference cannot be interpreted. Costs a second
                               # aSAM run (~2x the GPU time).
POCKET_ANALYSIS=true           # measure the cryptic pocket in both ensembles
POCKET_MIN_OVERLAP=2           # min cryptic residues lining the pocket
POCKET_MAX_DIST=12.0           # max Angstrom, pocket centroid to cryptic centroid
POCKET_GUIDED=false            # fpocket -P (guided detection). REQUIRES fpocket >= 4.2.
                               # 4.0 advertises -P but does not implement it: it exits 0,
                               # writes nothing, and every frame is recorded as a failure
                               # -- which reads exactly like a pocket that never opens.
                               # Stage 9 probes -P on one frame and falls back on its own,
                               # but leaving this false skips the wasted probe entirely.
                               # Nothing in the analysis depends on -P; it only narrowed
                               # the search. The overlap + distance test does the work.
POCKET_NPROC=16                # parallel fpocket workers. 16 is plenty -- the stage is
                               # I/O-bound, and more workers mainly multiply scratch use.
POCKET_SCRATCH=""              # optional: TMPDIR for the per-frame PDBs. Leave empty to
                               # use the system default. Point it at a filesystem with a
                               # few GB free if /tmp is small.

# ------------------------------------------------------------------------------
# 2. MODELLER License Key (Required)
# ------------------------------------------------------------------------------
# MODELLER requires a free academic license key from Sali Lab.
# Register for free at: https://salilab.org/modeller/registration.html
# Paste your key below, or export MODELLER_KEY in your shell to override it.
MODELLER_KEY="${MODELLER_KEY:-MODELIRANJE}"

# ------------------------------------------------------------------------------
# 3. Job Directory & Container Setup
# ------------------------------------------------------------------------------
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
WORKDIR="${JOBNAME}_${TIMESTAMP}"
CONTAINER="${SCRIPT_DIR}/pipeline_local.sif"

# ------------------------------------------------------------------------------
# 4. Preflight Checks
# ------------------------------------------------------------------------------
# Check for Apptainer or Singularity on the host
if command -v apptainer &>/dev/null; then
    CONTAINER_RUNTIME="apptainer"
elif command -v singularity &>/dev/null; then
    CONTAINER_RUNTIME="singularity"
else
    echo "======================================================================"
    echo "ERROR: Apptainer or Singularity is not installed on this system!"
    echo "Please install Apptainer: https://apptainer.org/docs/admin/main/installation.html"
    echo "======================================================================"
    exit 1
fi

# Check MODELLER license key
if [ -z "$MODELLER_KEY" ]; then
    echo "======================================================================"
    echo "ERROR: MODELLER license key is missing!"
    echo "Please set MODELLER_KEY at the top of this script:"
    echo "    MODELLER_KEY=\"your_key_here\""
    echo "Get a free key at: https://salilab.org/modeller/registration.html"
    echo "======================================================================"
    exit 1
fi

# Check that the container image exists
if [ ! -f "$CONTAINER" ]; then
    echo "======================================================================"
    echo "ERROR: Container image not found: $CONTAINER"
    echo "Build it first using:"
    echo "    cd ${SCRIPT_DIR} && sudo apptainer build pipeline_local.sif pipeline_local.def"
    echo "======================================================================"
    exit 1
fi

# Verify the image actually carries (and can read) its baked payload. Because
# this script uses no bind mounts, a stale or mis-permissioned .sif fails here
# rather than halfway through a long run.
echo "Verifying container payload..."
for REQUIRED in \
    /opt/pipeline/run_pipeline.py \
    /opt/pipeline/predict_cryptic_local.py \
    /opt/pipeline/cryptobank_model_loader.py \
    /opt/pipeline/esm_scanner.py \
    /opt/pipeline/modeller_mutate.py \
    /opt/pipeline/clean_structure.py \
    /opt/pipeline/make_report.py \
    /opt/pipeline/asam_cluster_analysis.py \
    /opt/pipeline/cryptic_pocket_analysis.py \
    /opt/pipeline/plot_cryptiscan.py \
    /opt/pipeline/diagnose_fpocket.py \
    /opt/conda/envs/analysis/bin/fpocket \
    /opt/conda/envs/prep/bin/python3 \
    /opt/models/esm1v_t33_650M_UR90S_1.pt \
    /opt/models/prot_t5_xl_uniref50_full_v2/cpt.pth \
    /opt/sam2_weights/weights/mdcath_1.0/nn.gen.pt
do
    if ! $CONTAINER_RUNTIME exec "$CONTAINER" test -r "$REQUIRED" 2>/dev/null; then
        echo "======================================================================"
        echo "ERROR: $CONTAINER is missing or cannot read: $REQUIRED"
        echo "The image is stale or was built from an older definition file."
        echo "Rebuild it with:"
        echo "    cd ${SCRIPT_DIR} && sudo apptainer build pipeline_local.sif pipeline_local.def"
        echo "======================================================================"
        exit 1
    fi
done
echo "Container payload OK."

mkdir -p "$WORKDIR"

echo "=========================================================="
echo " CryptiScan Pipeline (Local / Offline Version)"
echo " PDB: $PDB | Chain: $CHAIN | Top: $TOP | Strategy: $STRATEGY"
echo " Mode: $MODE | Workdir: $WORKDIR"
echo " Container: $CONTAINER"
echo " Runtime: $CONTAINER_RUNTIME"
echo "=========================================================="

# ------------------------------------------------------------------------------
# STAGES 1-4: Local CryptoBank -> Clean PDB -> ESM-Scan -> MODELLER
# ------------------------------------------------------------------------------
echo "Running Stages 1-4 inside container (local ProtT5 inference)..."

$CONTAINER_RUNTIME exec --nv \
    --env KEY_MODELLER="$MODELLER_KEY" \
    --env KEY_MODELLER10v8="$MODELLER_KEY" \
    --env MODINSTALL10v8="/opt/conda/lib/modeller-10.8" \
    --env PROT_T5_MODEL="/opt/models/prot_t5_xl_uniref50_full_v2" \
    "$CONTAINER" /opt/conda/bin/python3 /opt/pipeline/run_pipeline.py \
        --pdb "$PDB" \
        --chain "$CHAIN" \
        --top "$TOP" \
        --strategy "$STRATEGY" \
        --mode "$MODE" \
        --workdir "$WORKDIR" \
        --esm_model /opt/models/esm1v_t33_650M_UR90S_1.pt \
        --scraper /opt/pipeline/predict_cryptic_local.py \
        --esm_script /opt/pipeline/esm_scanner.py \
        --modeller_script /opt/pipeline/modeller_mutate.py \
        --prep_script /opt/pipeline/clean_structure.py \
        --mutation_set "$MUTATION_SET" \
        --esm_python /opt/conda/bin/python3 \
        --prep_python /opt/conda/envs/prep/bin/python3 \
        --modeller_python /opt/conda/bin/python3 \
        --scraper_python /opt/conda/bin/python3 \
        --check_numbering

if [ $? -ne 0 ]; then
    echo "ERROR: Stages 1-4 failed. Please review the output above."
    exit 1
fi

# ------------------------------------------------------------------------------
# STAGE 5: SAM2 Ensemble Generation
# ------------------------------------------------------------------------------
# The structure has been a single chain since Stage 2; run_pipeline.py Stage 5b
# re-validated the MODELLER output against SAM2's requirements (standard
# residues only, every residue has a CA).
if [ ! -f "$WORKDIR/sam2_input.pdb" ]; then
    echo "ERROR: $WORKDIR/sam2_input.pdb was not produced. Check pipeline output above."
    exit 1
fi
cp "$WORKDIR/sam2_input.pdb" "$WORKDIR/protein.pdb"
echo "SAM2 input ready (chain ${CHAIN}): $WORKDIR/protein.pdb"

# Run SAM2 ensemble generation using container
echo "Running SAM2 ensemble generation inside container..."
$CONTAINER_RUNTIME exec --nv \
    --env SAM_WEIGHTS_PATH="/opt/sam2_weights" \
    "$CONTAINER" /opt/conda/bin/python3 -c "import logging; logging.Logger.infox = getattr(logging.Logger, 'info', None); from sam.scripts.generate_ensemble import main; main()" \
        --config mdcath \
        -i "$WORKDIR/protein.pdb" \
        -o "$WORKDIR/ensemble_output" \
        --data_dir "/opt/sam2_weights" \
        -n 1000 \
        -b 4 \
        -T 320 \
        -d cuda

if [ -f "$WORKDIR/ensemble_output.top.pdb" ]; then
    echo "aSAM ensemble generation completed successfully!"
    echo "Ensemble files created: $WORKDIR/ensemble_output.*"
else
    echo "NOTE: aSAM ensemble generation skipped. Stages 1-4 results are preserved."
fi

# ------------------------------------------------------------------------------
# === ADDED ===  STAGE 5b: WILD-TYPE reference ensemble
# ------------------------------------------------------------------------------
# The wild type is the chain as it left Stage 2, before MODELLER. It goes through
# the same validation the mutant did, or aSAM rejects it. Settings MUST match the
# mutant run exactly -- a different temperature or frame count invalidates the
# comparison.
PDB_LOWER=$(echo "$PDB" | tr '[:upper:]' '[:lower:]')
WT_CHAIN_PDB="$WORKDIR/${PDB_LOWER}.pdb"

if [ "$RUN_WT_ENSEMBLE" = true ]; then
    if [ ! -f "$WT_CHAIN_PDB" ]; then
        echo "WARNING: $WT_CHAIN_PDB not found; skipping the wild-type ensemble."
        echo "         Stage 9 will have no reference to compare against."
    else
        echo "Validating wild-type structure for aSAM..."
        $CONTAINER_RUNTIME exec \
            "$CONTAINER" /opt/conda/envs/prep/bin/python3 \
            /opt/pipeline/clean_structure.py \
                --mode extract \
                --pdb "$WT_CHAIN_PDB" \
                --chain "$CHAIN" \
                --output "$WORKDIR/wt_input.pdb"

        echo "Running aSAM ensemble generation on the WILD TYPE..."
        $CONTAINER_RUNTIME exec --nv \
            --env SAM_WEIGHTS_PATH="/opt/sam2_weights" \
            "$CONTAINER" /opt/conda/bin/python3 -c "from sam.scripts.generate_ensemble import main; main()" \
                --config mdcath \
                -i "$WORKDIR/wt_input.pdb" \
                -o "$WORKDIR/wt_ensemble" \
                --data_dir "/opt/sam2_weights" \
                -n 1000 \
                -b 4 \
                -T 320 \
                -d cuda

        if [ -f "$WORKDIR/wt_ensemble.top.pdb" ]; then
            echo "Wild-type reference ensemble complete."
        else
            echo "NOTE: wild-type ensemble failed. Stage 9 will report the mutant"
            echo "      alone, which has no reference to compare against."
        fi
    fi
fi

# ------------------------------------------------------------------------------
# === STAGE 6 ===  Ensemble clustering (PCA + k-means on CA coordinates)
# ------------------------------------------------------------------------------
# Extracts N representative conformations from the mutant ensemble. Must run
# BEFORE Stage 7: make_report.py looks for ensemble_clusters/cluster_meta.sh and
# reports "Clusters: not run" if it is absent.
#
# Runs in the prep env (MDAnalysis + scipy). The third argument is a reference
# structure the script accepts but does not use -- it aligns to the trajectory's
# own first frame -- so the topology is passed to satisfy the positional count.
if [ "$RUN_CLUSTERING" = true ] && [ -f "$WORKDIR/ensemble_output.top.pdb" ]; then
    echo "Stage 6: clustering the mutant ensemble (k=$N_CLUSTERS, $PCA_DIM PCs)..."
    $CONTAINER_RUNTIME exec \
        "$CONTAINER" /opt/conda/envs/prep/bin/python3 \
        /opt/pipeline/asam_cluster_analysis.py \
            "$WORKDIR/ensemble_output.top.pdb" \
            "$WORKDIR/ensemble_output.traj.dcd" \
            "$WORKDIR/ensemble_output.top.pdb" \
            "$WORKDIR/ensemble_clusters" \
            "$N_CLUSTERS" \
            "$PCA_DIM" \
        || echo "WARNING: clustering failed; the run continues without it."

    if [ -f "$WORKDIR/ensemble_clusters/cluster_meta.sh" ]; then
        echo "Stage 6: representatives -> $WORKDIR/ensemble_clusters/clustering/"
    fi
fi

# ------------------------------------------------------------------------------
# STAGE 7: Human-readable report + FASTA sequences
# ------------------------------------------------------------------------------
echo "Building run report..."
$CONTAINER_RUNTIME exec "$CONTAINER" /opt/conda/bin/python3 /opt/pipeline/make_report.py \
    --workdir "$WORKDIR" \
    --jobname "$JOBNAME" \
    --container "$(basename "$CONTAINER")"

# ------------------------------------------------------------------------------
# === STAGE 9 ===  Cryptic pocket volume: mutant vs wild type
# ------------------------------------------------------------------------------
# Runs BEFORE the archive step so the results land inside the .zip.
#
# Residue numbers are mapped from the original PDB numbering to aSAM's 0-based
# numbering automatically, by walking the chain in file order. Never substitute a
# fixed offset: PDB numbering often has gaps (Ambler numbering in beta-lactamases)
# or insertion codes (antibodies), and subtraction then silently selects the wrong
# residues.
#
# --source-pdb is the WILD-TYPE chain, NOT sam2_input.pdb. The mutant's residue
# names no longer match the Stage 3 prediction, so the identity check would
# correctly reject it. MODELLER preserves numbering and residue order, so the
# wild-type chain maps both ensembles.
if [ "$POCKET_ANALYSIS" = true ] && [ -f "$WORKDIR/ensemble_output.top.pdb" ]; then
    echo "=========================================================="
    echo " Stage 9: quantifying the cryptic pocket"
    echo "=========================================================="

    # The MANIFEST, not the raw cryptic prediction. run_pipeline.py walks down
    # the crypticity ranking past residues whose wild type is already
    # charged/polar, so the residues actually mutated are NOT necessarily the
    # top-N by score. The manifest records which ones were mutated; the cryptic
    # JSON does not, and using it would analyse the wrong site whenever the run
    # skipped anything.
    RESIDUES_JSON="$WORKDIR/${PDB}_manifest.json"
    [ -f "$RESIDUES_JSON" ] || RESIDUES_JSON="$WORKDIR/${PDB_LOWER}_manifest.json"
    if [ ! -f "$RESIDUES_JSON" ]; then
        # fall back to the prediction; the script warns that it may not match
        RESIDUES_JSON="$WORKDIR/${PDB}_${CHAIN}_cryptic.json"
        [ -f "$RESIDUES_JSON" ] || RESIDUES_JSON="$WORKDIR/${PDB_LOWER}_${CHAIN}_cryptic.json"
    fi

    if [ ! -f "$RESIDUES_JSON" ]; then
        echo "WARNING: neither manifest nor cryptic JSON found; skipping Stage 9."
    elif [ ! -f "$WT_CHAIN_PDB" ]; then
        echo "WARNING: $WT_CHAIN_PDB not found; skipping Stage 9."
    else
        GUIDED_FLAG=""
        [ "$POCKET_GUIDED" = true ] || GUIDED_FLAG="--no-guided"
        [ -n "$POCKET_SCRATCH" ] && { mkdir -p "$POCKET_SCRATCH"; export TMPDIR="$POCKET_SCRATCH"; }

        WT_ARGS=""
        if [ -f "$WORKDIR/wt_ensemble.top.pdb" ]; then
            WT_ARGS="--wt-top $WORKDIR/wt_ensemble.top.pdb --wt-dcd $WORKDIR/wt_ensemble.traj.dcd"
        fi

        $CONTAINER_RUNTIME exec \
            "$CONTAINER" /opt/conda/bin/python3 /opt/pipeline/cryptic_pocket_analysis.py \
                --mutant-top "$WORKDIR/ensemble_output.top.pdb" \
                --mutant-dcd "$WORKDIR/ensemble_output.traj.dcd" \
                $WT_ARGS \
                --residues-json "$RESIDUES_JSON" \
                --source-pdb "$WT_CHAIN_PDB" \
                --chain "$CHAIN" \
                --outdir "$WORKDIR/pocket_analysis" \
                --min-overlap "$POCKET_MIN_OVERLAP" \
                --max-dist "$POCKET_MAX_DIST" \
                --fpocket /opt/conda/envs/analysis/bin/fpocket \
                $GUIDED_FLAG \
                --nproc "$POCKET_NPROC" \
            || echo "WARNING: Stage 9 failed; all earlier results are preserved."

        # --- Stage 9b: figures -------------------------------------------
        # Publication figures from the per-frame tables. Panels whose input is
        # missing are skipped with a printed reason rather than drawn empty, so
        # a failed detection stage cannot masquerade as a null result.
        if [ -d "$WORKDIR/pocket_analysis" ]; then
            $CONTAINER_RUNTIME exec \
                "$CONTAINER" /opt/conda/bin/python3 /opt/pipeline/plot_cryptiscan.py \
                    --indir "$WORKDIR/pocket_analysis" \
                    --clusters "$WORKDIR/ensemble_clusters/clustering/cluster_representatives.json" \
                    --outdir "$WORKDIR/figures" \
                || echo "NOTE: figure generation failed; the CSVs are unaffected." 
        fi
    fi
fi

# ------------------------------------------------------------------------------
# STAGE 8: Archive results and clean up directory
# ------------------------------------------------------------------------------
ZIP_FILE="${WORKDIR}.zip"
echo "Compressing all run results into ${ZIP_FILE}..."

# Use host zip if installed; otherwise use Python's built-in zipfile from container
if command -v zip &>/dev/null; then
    zip -r "${ZIP_FILE}" "${WORKDIR}"
else
    $CONTAINER_RUNTIME exec "$CONTAINER" /opt/conda/bin/python3 -m zipfile -c "${ZIP_FILE}" "${WORKDIR}"
fi

if [ -f "${ZIP_FILE}" ]; then
    echo "Archive created successfully: ${ZIP_FILE}"
    rm -rf "${WORKDIR}"
    echo "Cleaned up temporary directory: ${WORKDIR}/"
    echo "Result archive is ready: ${ZIP_FILE}"
else
    echo "WARNING: Could not create zip archive. Keeping raw folder: ${WORKDIR}/"
fi

echo "=========================================================="
echo " All pipeline steps completed successfully!"
echo " Output archive: $(pwd)/${ZIP_FILE}"
echo "=========================================================="
