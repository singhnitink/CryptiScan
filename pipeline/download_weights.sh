#!/bin/bash

# ==============================================================================
# CryptiScan - Download the model weights
# ==============================================================================
# The weight files are far too large for GitHub (~15 GB in total), so they are
# not in the repository. Run this script once after cloning to fetch them all
# into required_files_weights/, which is where the .def files expect them.
#
#     bash download_weights.sh
#
# It is safe to run again: anything already downloaded is skipped. If a download
# was interrupted, just run it again -- wget resumes where it left off.
#
# What gets downloaded:
#
#   ESM-1v      7.3 GB   Facebook Research   scores substitutions (Stage 4)
#   ProtT5      6.8 GB   HuggingFace         local cryptic prediction (Stage 3)
#   aSAM/SAM2   1.0 GB   GitHub release      conformational ensembles (Stage 6)
#
# ProtT5 is only needed by the LOCAL image. If you are only building the web
# image, run:  bash download_weights.sh --skip-prot-t5
#
# Requirements: wget and unzip. On Ubuntu/Debian:
#     sudo apt-get install wget unzip
# ==============================================================================

set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEST="${SCRIPT_DIR}/required_files_weights"

SKIP_PROT_T5=0
for arg in "$@"; do
    case "$arg" in
        --skip-prot-t5) SKIP_PROT_T5=1 ;;
        -h|--help)
            sed -n '3,26p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
            exit 0 ;;
        *)
            echo "Unknown option: $arg  (try --help)"
            exit 1 ;;
    esac
done

# ------------------------------------------------------------------------------
# Checks before we start
# ------------------------------------------------------------------------------
for tool in wget unzip; do
    if ! command -v "$tool" &>/dev/null; then
        echo "ERROR: '$tool' is not installed."
        echo "       On Ubuntu/Debian:  sudo apt-get install wget unzip"
        exit 1
    fi
done

mkdir -p "$DEST"

echo "=========================================================="
echo " CryptiScan - downloading model weights"
echo " Destination: $DEST"
echo "=========================================================="
echo
echo "This downloads about 15 GB and can take a while on a slow"
echo "connection. You can stop it (Ctrl+C) and re-run it later;"
echo "it picks up where it left off."
echo

# Free space check: 16 GB for the files, plus room to unzip aSAM.
AVAIL_GB=$(df -BG --output=avail "$DEST" | tail -1 | tr -dc '0-9')
if [ -n "$AVAIL_GB" ] && [ "$AVAIL_GB" -lt 18 ]; then
    echo "WARNING: only ${AVAIL_GB} GB free here; about 18 GB is needed."
    echo "         Free some space, or Ctrl+C now."
    echo
fi

# ------------------------------------------------------------------------------
# Helper: download one file, skipping it if it is already complete
# ------------------------------------------------------------------------------
fetch() {
    local url="$1" out="$2" label="$3"

    # Compare against the size the server reports, so a half-finished file is
    # resumed rather than mistaken for a good one.
    local remote_size
    remote_size=$(wget --spider --server-response "$url" 2>&1 \
        | tr -d '\r' | awk 'BEGIN{IGNORECASE=1} /^ *Content-Length:/{s=$2} END{print s+0}')

    if [ -f "$out" ]; then
        local local_size
        local_size=$(stat -c%s "$out" 2>/dev/null || echo 0)
        if [ "$remote_size" -gt 0 ] && [ "$local_size" -eq "$remote_size" ]; then
            echo "  [skip] $label -- already downloaded"
            return 0
        fi
        echo "  [resume] $label -- incomplete, continuing"
    else
        echo "  [get]  $label"
    fi

    mkdir -p "$(dirname "$out")"
    if ! wget --continue --show-progress --progress=bar:force -O "$out" "$url"; then
        echo
        echo "ERROR: download failed: $label"
        echo "       URL: $url"
        echo "       Re-run this script to try again."
        exit 1
    fi

    # Confirm the finished file is exactly the size the server advertised. A
    # server that ignores resume requests can leave a file that looks complete
    # but is not; deleting it here turns that into an obvious error rather than
    # a corrupt weight file that only fails much later inside the container.
    if [ "$remote_size" -gt 0 ]; then
        local final_size
        final_size=$(stat -c%s "$out" 2>/dev/null || echo 0)
        if [ "$final_size" -ne "$remote_size" ]; then
            rm -f "$out"
            echo
            echo "ERROR: $label finished at ${final_size} bytes but should be ${remote_size}."
            echo "       The partial file has been deleted. Re-run this script."
            exit 1
        fi
    fi
}

# ------------------------------------------------------------------------------
# 1. ESM-1v  (Stage 4: scores substitutions at each cryptic residue)
# ------------------------------------------------------------------------------
echo "----------------------------------------------------------"
echo " 1/3  ESM-1v  (7.3 GB)"
echo "----------------------------------------------------------"
fetch "https://dl.fbaipublicfiles.com/fair-esm/models/esm1v_t33_650M_UR90S_1.pt" \
      "${DEST}/esm1v_t33_650M_UR90S_1.pt" \
      "esm1v_t33_650M_UR90S_1.pt"
echo

# ------------------------------------------------------------------------------
# 2. ProtT5 + CryptoBank classifier head  (Stage 3, LOCAL image only)
#    From https://huggingface.co/ThorbenF/prot_t5_xl_uniref50_full_v2
#    Fetched file by file on purpose: "git clone" of this repo drags in a
#    ~6.8 GB .git/lfs cache holding a second copy of the same weights.
# ------------------------------------------------------------------------------
echo "----------------------------------------------------------"
echo " 2/3  ProtT5 + classifier head  (6.8 GB)"
echo "----------------------------------------------------------"
if [ "$SKIP_PROT_T5" -eq 1 ]; then
    echo "  [skip] --skip-prot-t5 given (web image does not need it)"
else
    HF_BASE="https://huggingface.co/ThorbenF/prot_t5_xl_uniref50_full_v2/resolve/main"
    PT5="${DEST}/prot_t5_xl_uniref50_full_v2"
    # NOTE: the directory name must keep the substring "prot_t5" --
    # cryptobank_model_loader.py picks its code path from the path string.
    for f in model.safetensors cpt.pth spiece.model config.json \
             tokenizer_config.json special_tokens_map.json added_tokens.json; do
        fetch "${HF_BASE}/${f}" "${PT5}/${f}" "prot_t5_xl_uniref50_full_v2/${f}"
    done
fi
echo

# ------------------------------------------------------------------------------
# 3. aSAM / SAM2 weights  (Stage 6: conformational ensembles)
#    Published as a GitHub release zip. Extracting it into sam2_weights/weights/
#    produces weights/mdcath_1.0/, which is the layout SAM2 looks for under
#    --data_dir.
# ------------------------------------------------------------------------------
echo "----------------------------------------------------------"
echo " 3/3  aSAM / SAM2 weights  (1.0 GB)"
echo "----------------------------------------------------------"
SAM_DIR="${DEST}/sam2_weights"
if [ -f "${SAM_DIR}/weights/mdcath_1.0/nn.gen.pt" ]; then
    echo "  [skip] sam2_weights/weights/mdcath_1.0 -- already present"
else
    ZIP="${SAM_DIR}/mdcath_1.0.zip"
    fetch "https://github.com/giacomo-janson/sam2/releases/download/data-1.0/mdcath_1.0.zip" \
          "$ZIP" "mdcath_1.0.zip"
    echo "  [unzip] into sam2_weights/weights/"
    mkdir -p "${SAM_DIR}/weights"
    unzip -q -o "$ZIP" -d "${SAM_DIR}/weights"
    rm -f "$ZIP"
fi
echo

# ------------------------------------------------------------------------------
# Make everything world-readable.
# %files in the .def keeps the permissions it finds, and files that arrive as
# 0600 would be baked into the image unreadable by anyone but you.
# ------------------------------------------------------------------------------
chmod -R u+rwX,go+rX,go-w "$DEST" 2>/dev/null || true

# ------------------------------------------------------------------------------
# Check what we ended up with
# ------------------------------------------------------------------------------
echo "=========================================================="
echo " Checking the downloads"
echo "=========================================================="

MISSING=0
check_file() {
    if [ -s "$1" ]; then
        printf "  OK       %-46s %6s\n" "$2" "$(du -h "$1" | cut -f1)"
    else
        printf "  MISSING  %-46s\n" "$2"
        MISSING=1
    fi
}

check_file "${DEST}/esm1v_t33_650M_UR90S_1.pt" "esm1v_t33_650M_UR90S_1.pt"
if [ "$SKIP_PROT_T5" -eq 0 ]; then
    check_file "${DEST}/prot_t5_xl_uniref50_full_v2/model.safetensors" "prot_t5_xl_uniref50_full_v2/model.safetensors"
    check_file "${DEST}/prot_t5_xl_uniref50_full_v2/cpt.pth"           "prot_t5_xl_uniref50_full_v2/cpt.pth"
    check_file "${DEST}/prot_t5_xl_uniref50_full_v2/spiece.model"      "prot_t5_xl_uniref50_full_v2/spiece.model"
fi
for f in nn.gen.pt nn.enc.pt nn.dec.pt enc_std_scaler.pt; do
    check_file "${SAM_DIR}/weights/mdcath_1.0/${f}" "sam2_weights/weights/mdcath_1.0/${f}"
done

echo
echo "  Total in required_files_weights: $(du -sh "$DEST" 2>/dev/null | cut -f1)"
echo

if [ "$MISSING" -ne 0 ]; then
    echo "=========================================================="
    echo " Some files are missing. Re-run this script to retry."
    echo "=========================================================="
    exit 1
fi

echo "=========================================================="
echo " All weights downloaded."
echo
echo " Next, build the container images:"
echo "     cd ${SCRIPT_DIR}"
echo "     sudo apptainer build pipeline_web.sif   pipeline_web.def"
if [ "$SKIP_PROT_T5" -eq 0 ]; then
echo "     sudo apptainer build pipeline_local.sif pipeline_local.def"
fi
echo "=========================================================="
