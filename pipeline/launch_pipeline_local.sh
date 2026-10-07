#!/bin/bash

# CryptiScan Pipeline Launcher (Local / Offline Version)
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
#   6b. Cluster the ensemble and keep one representative frame per cluster
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

# Automatically detect the directory where this script lives
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# 1. User Configuration (Edit these for your target protein)
PDB="1JWP"                     # PDB accession code (4 characters) or local PDB path
CHAIN="A"                      # Target protein chain identifier
TOP=5                          # Number of residues to mutate; a charged/polar wild type is
                               # skipped and the next-ranked cryptic residue is used instead
MODE="combined"                # "combined" = all mutations in 1 structure; "independent" = N separate mutants
STRATEGY="esm"                 # "esm" = evaluates substitutions via ESM-Scan protein language model
MUTATION_SET="charged_polar"   # Which residues ESM may choose from:
                               #   charged_polar (default) | charged | polar | all | "ASP,GLU,..."
                               # Charged/polar substitutions destabilise the closed state so the
                               # cryptic pocket can open in MD. "all" reverts to the unrestricted
                               # 19-way scan, which tends to pick conservative nonpolar swaps
                               # (ILE->VAL) that leave the pocket shut.
N_CLUSTERS=10                  # Max clusters for the SAM2 ensemble (one representative frame each)
PCA_DIM=10                     # PCA components used for clustering
JOBNAME="${JOBNAME:-1jwp}_local"     # Job name for folder and zip naming (e.g. "1lzt", "1jwp")

# 2. MODELLER License Key (Required)
# MODELLER requires a free academic license key from Sali Lab.
# Register for free at: https://salilab.org/modeller/registration.html
# Paste your key below, or export MODELLER_KEY in your shell to override it.
MODELLER_KEY="${MODELLER_KEY:-MODELIRANJE}"

# 3. Job Directory & Container Setup
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
WORKDIR="${JOBNAME}_${TIMESTAMP}"
CONTAINER="${SCRIPT_DIR}/pipeline_local.sif"

# 4. Preflight Checks
# Check for Apptainer or Singularity on the host
if command -v apptainer &>/dev/null; then
    CONTAINER_RUNTIME="apptainer"
elif command -v singularity &>/dev/null; then
    CONTAINER_RUNTIME="singularity"
else
    echo "ERROR: Apptainer or Singularity is not installed on this system!"
    echo "Please install Apptainer: https://apptainer.org/docs/admin/main/installation.html"
    exit 1
fi

# Check MODELLER license key
if [ -z "$MODELLER_KEY" ]; then
    echo "ERROR: MODELLER license key is missing!"
    echo "Please set MODELLER_KEY at the top of this script:"
    echo "    MODELLER_KEY=\"your_key_here\""
    echo "Get a free key at: https://salilab.org/modeller/registration.html"
    exit 1
fi

# Check that the container image exists
if [ ! -f "$CONTAINER" ]; then
    echo "ERROR: Container image not found: $CONTAINER"
    echo "Build it first using:"
    echo "    cd ${SCRIPT_DIR} && sudo apptainer build pipeline_local.sif pipeline_local.def"
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
    /opt/conda/envs/prep/bin/python3 \
    /opt/models/esm1v_t33_650M_UR90S_1.pt \
    /opt/models/prot_t5_xl_uniref50_full_v2/cpt.pth \
    /opt/sam2_weights/weights/mdcath_1.0/nn.gen.pt
do
    if ! $CONTAINER_RUNTIME exec "$CONTAINER" test -r "$REQUIRED" 2>/dev/null; then
        echo "ERROR: $CONTAINER is missing or cannot read: $REQUIRED"
        echo "The image is stale or was built from an older definition file."
        echo "Rebuild it with:"
        echo "    cd ${SCRIPT_DIR} && sudo apptainer build pipeline_local.sif pipeline_local.def"
        exit 1
    fi
done
echo "Container payload OK."

mkdir -p "$WORKDIR"

echo "CryptiScan Pipeline (Local / Offline Version)"
echo "PDB: $PDB | Chain: $CHAIN | Top: $TOP | Strategy: $STRATEGY"
echo "Mode: $MODE | Workdir: $WORKDIR"
echo "Container: $CONTAINER"
echo "Runtime: $CONTAINER_RUNTIME"

# STAGES 1-4: Local CryptoBank -> Clean PDB -> ESM-Scan -> MODELLER
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

# STAGE 5: SAM2 Ensemble Generation
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
    echo "SAM2 ensemble generation completed successfully!"
    echo "Ensemble files created: $WORKDIR/ensemble_output.*"
else
    echo "NOTE: SAM2 ensemble generation skipped. Stages 1-4 results are preserved."
fi

# STAGE 6: Cluster the SAM2 ensemble
# PCA + k-means on CA coordinates; the frame nearest each cluster centre is kept
# as that cluster's representative. Runs in the prep env (MDAnalysis + scipy).
if [ -f "$WORKDIR/ensemble_output.traj.dcd" ]; then
    echo "Clustering SAM2 ensemble (up to ${N_CLUSTERS} clusters)..."
    $CONTAINER_RUNTIME exec "$CONTAINER" /opt/conda/envs/prep/bin/python3 /opt/pipeline/asam_cluster_analysis.py \
        "$WORKDIR/ensemble_output.top.pdb" \
        "$WORKDIR/ensemble_output.traj.dcd" \
        "$WORKDIR/sam2_input.pdb" \
        "$WORKDIR/ensemble_clusters" \
        "$N_CLUSTERS" \
        "$PCA_DIM" \
        || echo "WARNING: Ensemble clustering failed. The full ensemble is still in the archive."
else
    echo "NOTE: No SAM2 trajectory found, skipping ensemble clustering."
fi

# STAGE 7: Human-readable report + FASTA sequences
echo "Building run report..."
$CONTAINER_RUNTIME exec "$CONTAINER" /opt/conda/bin/python3 /opt/pipeline/make_report.py \
    --workdir "$WORKDIR" \
    --jobname "$JOBNAME" \
    --container "$(basename "$CONTAINER")"

# STAGE 8: Archive results and clean up directory
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

echo "All pipeline steps completed successfully!"
echo "Output archive: $(pwd)/${ZIP_FILE}"
