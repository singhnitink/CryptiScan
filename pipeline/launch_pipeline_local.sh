#!/bin/bash

# ==============================================================================
# CryptiScan Pipeline Launcher (Local / Offline Version)
# ==============================================================================
# EDIT ONLY SECTION 1. Everything below it has working defaults.
#
# Workflow (twelve stages, all automatic):
#    1. Fetch the entry from RCSB and split protein / ligands (waters dropped)
#    2. Select chain $CHAIN; every later stage uses this chain
#    3. Verify backbone continuity and standard residues
#    4. Score per-residue crypticity locally with ProtT5 (no network call)
#    5. Score substitutions with ESM-1v, restricted to charged/polar residues
#    6. Build the relaxed mutant with PDBFixer + OpenMM, then validate it for aSAM
#    7. Generate the MUTANT ensemble with aSAM
#    8. Generate the matched WILD-TYPE ensemble with aSAM
#    9. Cluster the mutant ensemble (PCA + k-means)
#   10. Quantify the cryptic pocket in BOTH ensembles (fpocket, guided)
#   11. Draw the figures and write the statistics table
#   12. Build the run report and archive everything
#
# NOTHING asks you for residue numbers. The cryptic residues come from the
# crypticity ranking, and Stage 10 maps them into aSAM's renumbered topology by
# walking the chain in file order, so no offset is ever assumed.
#
# Host requirements: Linux, an NVIDIA driver >= 525, and Apptainer. No host
# Python, CUDA toolkit or scientific packages. Everything else is in the .sif.
# Stage 1 still reaches files.rcsb.org, so point PDB at a local .pdb path if
# this machine has no internet access.
# ==============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# ==============================================================================
# 1. USER CONFIGURATION  --  this is the only block you need to edit
# ==============================================================================
PDB="1JWP"                     # PDB accession code, or a path to a local .pdb
CHAIN="A"                      # chain to analyse

# ==============================================================================
# 2. ADVANCED SETTINGS  --  defaults are the ones used in the paper
# ==============================================================================
JOBNAME="${JOBNAME:-${PDB,,}}_local"

# --- design ---
TOP=5                          # cryptic residues to MUTATE. Residues whose wild
                               # type is already charged/polar are scored but
                               # skipped, and the scan walks further down the
                               # ranking, so this many are always delivered.
MODE="combined"                # combined = all mutations in one structure
STRATEGY="esm"                 # esm = ESM-1v chooses the substitution
MUTATION_SET="charged_polar"   # charged_polar | charged | polar | all | "ASP,GLU,..."
                               # Burying a charge or polar group is what
                               # destabilises the closed state. "all" reverts to
                               # the unrestricted scan, which picks conservative
                               # nonpolar swaps and leaves the pocket shut.

# --- ensembles (Stages 7 and 8). BOTH ensembles use these values; a difference
#     in frames, temperature or seed between them invalidates the comparison. ---
FRAMES=2000                    # conformations per ensemble
TEMPERATURE=300                # kelvin, passed to aSAM as -T
BATCH=4

RUN_WT_ENSEMBLE=true           # the wild-type reference. NOT optional for a
                               # meaningful result: a pocket volume with no
                               # reference cannot be interpreted. Costs a
                               # second aSAM run.

# --- clustering (Stage 9) ---
RUN_CLUSTERING=true
N_CLUSTERS=10
PCA_DIM=10

# --- pocket quantification (Stage 10) ---
POCKET_ANALYSIS=true
POCKET_MIN_EXPLICIT=2          # fpocket -u: how many of an alpha sphere's four
                               # contact atoms must belong to the cryptic
                               # residues. fpocket's own default is 4, i.e. all
                               # four, which returns NO pocket at all for a site
                               # of five residues. This value changes the
                               # absolute volume, so keep one setting across
                               # every run you intend to compare and report it.
POCKET_MIN_OVERLAP=2           # min cryptic residues lining the pocket
POCKET_MAX_DIST=12.0           # max A, pocket centroid to cryptic-site centroid
POCKET_SPECIFICITY=true        # also run unguided detection, giving the summed
                               # volume of every OTHER pocket. This is what
                               # separates a site-specific opening from a
                               # general loosening of the surface. Doubles the
                               # Stage 10 runtime.
POCKET_NPROC=$(nproc)
POCKET_SCRATCH="${TMPDIR:-/tmp}"   # LOCAL disk. Per-frame temporaries on
                                   # network storage dominate the runtime.

# --- figures (Stage 11) ---
MAKE_FIGURES=true

# ==============================================================================
# 3. Job directory and container
# ==============================================================================
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
WORKDIR="${JOBNAME}_${TIMESTAMP}"
CONTAINER="${SCRIPT_DIR}/pipeline_local_openmm.sif"
PDB_LOWER=$(echo "$PDB" | tr '[:upper:]' '[:lower:]')

# ==============================================================================
# 4. Preflight
# ==============================================================================
if command -v apptainer &>/dev/null; then
    CONTAINER_RUNTIME="apptainer"
elif command -v singularity &>/dev/null; then
    CONTAINER_RUNTIME="singularity"
else
    echo "ERROR: neither Apptainer nor Singularity is installed."
    echo "       https://apptainer.org/docs/admin/main/installation.html"
    exit 1
fi

if [ ! -f "$CONTAINER" ]; then
    echo "ERROR: container image not found: $CONTAINER"
    echo "       cd ${SCRIPT_DIR} && sudo apptainer build pipeline_local.sif pipeline_local.def"
    exit 1
fi

# A stale .sif fails here rather than three hours in. This script passes no
# --bind flags on purpose: a bind mount would silently mask a stale image with
# host files that a recipient of the .sif will not have.
echo "Verifying container payload..."
for REQUIRED in \
    /opt/pipeline/run_pipeline.py \
    /opt/pipeline/predict_cryptic_local.py \
    /opt/pipeline/cryptobank_model_loader.py \
    /opt/pipeline/esm_scanner.py \
    /opt/pipeline/openmm_mutate.py \
    /opt/pipeline/clean_structure.py \
    /opt/pipeline/make_report.py \
    /opt/pipeline/asam_cluster_analysis.py \
    /opt/pipeline/cryptic_pocket_analysis.py \
    /opt/pipeline/plot_cryptiscan.py \
    /opt/conda/envs/analysis/bin/fpocket \
    /opt/conda/envs/prep/bin/python3 \
    /opt/conda/envs/mutate/bin/python3 \
    /opt/models/esm1v_t33_650M_UR90S_1.pt \
    /opt/models/prot_t5_xl_uniref50_full_v2/cpt.pth \
    /opt/sam2_weights/weights/mdcath_1.0/nn.gen.pt
do
    if ! $CONTAINER_RUNTIME exec "$CONTAINER" test -r "$REQUIRED" 2>/dev/null; then
        echo "ERROR: $CONTAINER cannot read: $REQUIRED"
        echo "       The image is stale or was built from an older definition."
        echo "       cd ${SCRIPT_DIR} && sudo apptainer build pipeline_local.sif pipeline_local.def"
        exit 1
    fi
done
echo "Container payload OK."

mkdir -p "$WORKDIR"

echo "=========================================================="
echo " CryptiScan (local / offline)"
echo " PDB: $PDB   Chain: $CHAIN   Mutations: $TOP   Set: $MUTATION_SET"
echo " Ensembles: 2 x $FRAMES frames at ${TEMPERATURE} K"
echo " Workdir: $WORKDIR"
echo "=========================================================="

# ==============================================================================
# STAGES 1-6: curate -> crypticity -> ESM-1v -> mutagenesis -> validate
# ==============================================================================
echo
echo "### Stages 1-6: curation, crypticity, ESM-1v scan, mutagenesis"
$CONTAINER_RUNTIME exec --nv \
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
        --mutate_script /opt/pipeline/openmm_mutate.py \
        --prep_script /opt/pipeline/clean_structure.py \
        --mutation_set "$MUTATION_SET" \
        --esm_python /opt/conda/bin/python3 \
        --prep_python /opt/conda/envs/prep/bin/python3 \
        --mutate_python /opt/conda/envs/mutate/bin/python3 \
        --scraper_python /opt/conda/bin/python3 \
        --check_numbering

if [ $? -ne 0 ]; then
    echo "ERROR: stages 1-6 failed. See the output above."
    exit 1
fi

# ==============================================================================
# STAGE 7: mutant ensemble
# ==============================================================================
if [ ! -f "$WORKDIR/sam2_input.pdb" ]; then
    echo "ERROR: $WORKDIR/sam2_input.pdb was not produced."
    exit 1
fi
cp "$WORKDIR/sam2_input.pdb" "$WORKDIR/protein.pdb"

echo
echo "### Stage 7: aSAM ensemble, MUTANT ($FRAMES frames at ${TEMPERATURE} K)"
$CONTAINER_RUNTIME exec --nv \
    --env SAM_WEIGHTS_PATH="/opt/sam2_weights" \
    "$CONTAINER" /opt/conda/bin/python3 -c "from sam.scripts.generate_ensemble import main; main()" \
        --config mdcath \
        -i "$WORKDIR/protein.pdb" \
        -o "$WORKDIR/ensemble_output" \
        --data_dir "/opt/sam2_weights" \
        -n "$FRAMES" -b "$BATCH" -T "$TEMPERATURE" -d cuda

if [ -f "$WORKDIR/ensemble_output.top.pdb" ]; then
    echo "Mutant ensemble complete."
else
    echo "NOTE: mutant ensemble failed. Stages 1-6 results are preserved."
fi

# ==============================================================================
# STAGE 8: wild-type reference ensemble
# ==============================================================================
# The wild type is the chain as it left Stage 2, before mutagenesis. Settings MUST
# match Stage 7 exactly or the comparison is void.
WT_CHAIN_PDB="$WORKDIR/${PDB_LOWER}.pdb"

if [ "$RUN_WT_ENSEMBLE" = true ]; then
    echo
    echo "### Stage 8: aSAM ensemble, WILD TYPE (same settings)"
    if [ ! -f "$WT_CHAIN_PDB" ]; then
        echo "WARNING: $WT_CHAIN_PDB not found; skipping the wild-type ensemble."
        echo "         Stage 10 will then have no reference to compare against."
    else
        $CONTAINER_RUNTIME exec \
            "$CONTAINER" /opt/conda/envs/prep/bin/python3 \
            /opt/pipeline/clean_structure.py \
                --mode extract --pdb "$WT_CHAIN_PDB" --chain "$CHAIN" \
                --output "$WORKDIR/wt_input.pdb"

        $CONTAINER_RUNTIME exec --nv \
            --env SAM_WEIGHTS_PATH="/opt/sam2_weights" \
            "$CONTAINER" /opt/conda/bin/python3 -c "from sam.scripts.generate_ensemble import main; main()" \
                --config mdcath \
                -i "$WORKDIR/wt_input.pdb" \
                -o "$WORKDIR/wt_ensemble" \
                --data_dir "/opt/sam2_weights" \
                -n "$FRAMES" -b "$BATCH" -T "$TEMPERATURE" -d cuda

        [ -f "$WORKDIR/wt_ensemble.top.pdb" ] \
            && echo "Wild-type reference ensemble complete." \
            || echo "NOTE: wild-type ensemble failed; Stage 10 has no reference."
    fi
fi

# ==============================================================================
# STAGE 9: clustering
# ==============================================================================
# Must run BEFORE the report: make_report.py reads ensemble_clusters/cluster_meta.sh.
if [ "$RUN_CLUSTERING" = true ] && [ -f "$WORKDIR/ensemble_output.top.pdb" ]; then
    echo
    echo "### Stage 9: clustering the mutant ensemble (k=$N_CLUSTERS, $PCA_DIM PCs)"
    $CONTAINER_RUNTIME exec \
        "$CONTAINER" /opt/conda/envs/prep/bin/python3 \
        /opt/pipeline/asam_cluster_analysis.py \
            "$WORKDIR/ensemble_output.top.pdb" \
            "$WORKDIR/ensemble_output.traj.dcd" \
            "$WORKDIR/ensemble_output.top.pdb" \
            "$WORKDIR/ensemble_clusters" \
            "$N_CLUSTERS" "$PCA_DIM" \
        || echo "WARNING: clustering failed; the run continues without it."
fi

# ==============================================================================
# STAGE 10: cryptic pocket quantification
# ==============================================================================
# Runs BEFORE the report and the archive, so both pick up its output.
#
# Residues come from the MANIFEST, not the raw crypticity JSON. The pipeline
# walks down the ranking past residues whose wild type is already charged or
# polar, so the residues actually mutated are NOT the top-N by score. The
# manifest records which ones were mutated; the cryptic JSON does not.
#
# --source-pdb is the WILD-TYPE chain, never sam2_input.pdb: the mutant's
# residue names no longer match the Stage 4 prediction, so the identity check
# would correctly reject it. Mutagenesis preserves numbering and residue order, so
# the wild-type chain maps both ensembles.
if [ "$POCKET_ANALYSIS" = true ] && [ -f "$WORKDIR/ensemble_output.top.pdb" ]; then
    echo
    echo "### Stage 10: quantifying the cryptic pocket"

    RESIDUES_JSON="$WORKDIR/${PDB}_manifest.json"
    [ -f "$RESIDUES_JSON" ] || RESIDUES_JSON="$WORKDIR/${PDB_LOWER}_manifest.json"
    if [ ! -f "$RESIDUES_JSON" ]; then
        RESIDUES_JSON="$WORKDIR/${PDB}_${CHAIN}_cryptic.json"
        [ -f "$RESIDUES_JSON" ] || RESIDUES_JSON="$WORKDIR/${PDB_LOWER}_${CHAIN}_cryptic.json"
    fi

    SPEC_FLAG=""
    [ "$POCKET_SPECIFICITY" = true ] || SPEC_FLAG="--no-specificity"

    if [ ! -f "$RESIDUES_JSON" ]; then
        echo "WARNING: neither manifest nor cryptic JSON found; skipping Stage 10."
    elif [ ! -f "$WT_CHAIN_PDB" ]; then
        echo "WARNING: $WT_CHAIN_PDB not found; skipping Stage 10."
    else
        WT_ARGS=""
        if [ -f "$WORKDIR/wt_ensemble.top.pdb" ]; then
            WT_ARGS="--wt-top $WORKDIR/wt_ensemble.top.pdb --wt-dcd $WORKDIR/wt_ensemble.traj.dcd"
        fi

        $CONTAINER_RUNTIME exec \
            --env TMPDIR="$POCKET_SCRATCH" \
            "$CONTAINER" /opt/conda/bin/python3 /opt/pipeline/cryptic_pocket_analysis.py \
                --mutant-top "$WORKDIR/ensemble_output.top.pdb" \
                --mutant-dcd "$WORKDIR/ensemble_output.traj.dcd" \
                $WT_ARGS \
                --residues-json "$RESIDUES_JSON" \
                --source-pdb "$WT_CHAIN_PDB" \
                --chain "$CHAIN" \
                --outdir "$WORKDIR/pocket_analysis" \
                --min-explicit "$POCKET_MIN_EXPLICIT" \
                --min-overlap "$POCKET_MIN_OVERLAP" \
                --max-dist "$POCKET_MAX_DIST" \
                $SPEC_FLAG \
                --fpocket /opt/conda/envs/analysis/bin/fpocket \
                --nproc "$POCKET_NPROC" \
            || echo "WARNING: Stage 10 failed; all earlier results are preserved."
    fi
fi

# ==============================================================================
# STAGE 11: figures and statistics
# ==============================================================================
if [ "$MAKE_FIGURES" = true ] && [ -f "$WORKDIR/pocket_analysis/mutant_pocket.csv" ]; then
    echo
    echo "### Stage 11: figures"
    CLUSTERS_JSON="$WORKDIR/ensemble_clusters/clustering/cluster_representatives.json"
    CLUS_ARG=""
    [ -f "$CLUSTERS_JSON" ] && CLUS_ARG="--clusters $CLUSTERS_JSON"

    $CONTAINER_RUNTIME exec \
        "$CONTAINER" /opt/conda/bin/python3 /opt/pipeline/plot_cryptiscan.py \
            --indir "$WORKDIR/pocket_analysis" \
            $CLUS_ARG \
            --outdir "$WORKDIR/figures" \
            --label "$PDB chain $CHAIN" \
        || echo "WARNING: plotting failed; the CSVs are intact."
fi

# ==============================================================================
# STAGE 12: report and archive
# ==============================================================================
echo
echo "### Stage 12: report"
$CONTAINER_RUNTIME exec "$CONTAINER" /opt/conda/bin/python3 /opt/pipeline/make_report.py \
    --workdir "$WORKDIR" \
    --jobname "$JOBNAME" \
    --container "$(basename "$CONTAINER")"

ZIP_FILE="${WORKDIR}.zip"
echo "Compressing into ${ZIP_FILE}..."
if command -v zip &>/dev/null; then
    zip -r "${ZIP_FILE}" "${WORKDIR}" > /dev/null
else
    $CONTAINER_RUNTIME exec "$CONTAINER" /opt/conda/bin/python3 -m zipfile -c "${ZIP_FILE}" "${WORKDIR}"
fi

if [ -f "${ZIP_FILE}" ]; then
    echo "Archive created: ${ZIP_FILE}"
    rm -rf "${WORKDIR}"
else
    echo "WARNING: could not create the archive. Keeping ${WORKDIR}/"
fi

echo "=========================================================="
echo " Done."
echo " Archive: $(pwd)/${ZIP_FILE}"
echo "=========================================================="
