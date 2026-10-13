#!/bin/bash
# Download CryptiScan model weights (~15 GB total).
# Run once after cloning:  bash download_weights.sh
# Use --skip-prot-t5 if you only need the web image.
# Requires: wget, unzip

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEST="${SCRIPT_DIR}/required_files_weights"

SKIP_PROT_T5=0
for arg in "$@"; do
    case "$arg" in
        --skip-prot-t5) SKIP_PROT_T5=1 ;;
        -h|--help)
            echo "Usage: bash download_weights.sh [--skip-prot-t5]"
            echo "Downloads model weights into required_files_weights/."
            echo "  --skip-prot-t5   Skip ProtT5 weights (not needed for web image)."
            exit 0 ;;
        *) echo "Unknown option: $arg (try --help)"; exit 1 ;;
    esac
done

for tool in wget unzip; do
    command -v "$tool" &>/dev/null || { echo "ERROR: '$tool' is not installed."; exit 1; }
done

mkdir -p "$DEST"

# Download a file. Skips if it already exists and is non-empty.
fetch() {
    local url="$1" out="$2" label="$3"
    mkdir -p "$(dirname "$out")"

    if [ -f "$out" ] && [ -s "$out" ]; then
        echo "  [skip] $label (already exists)"
        return 0
    fi

    echo "  [downloading] $label"
    local attempt wait_time
    for attempt in 1 2 3; do
        if wget --show-progress --max-redirect=10 -O "$out" "$url"; then
            break
        fi
        wait_time=$((attempt * 20))
        echo "  [retry $attempt/3] waiting ${wait_time}s before retrying $label..."
        sleep "$wait_time"
    done

    if [ ! -s "$out" ]; then
        rm -f "$out"
        echo "  ERROR: failed to download $label after 3 attempts."
        echo "  URL: $url"
        exit 1
    fi
    # Brief pause between downloads to avoid rate-limiting (HuggingFace CDN).
    sleep 3
}

echo "========================================="
echo " CryptiScan — downloading model weights"
echo " Destination: $DEST"
echo "========================================="

# --- 1. ESM-1v (7.3 GB) ---
echo ""
echo "[1/3] ESM-1v (7.3 GB)"
fetch "https://dl.fbaipublicfiles.com/fair-esm/models/esm1v_t33_650M_UR90S_1.pt" \
      "${DEST}/esm1v_t33_650M_UR90S_1.pt" \
      "esm1v_t33_650M_UR90S_1.pt"

# --- 2. ProtT5 + classifier head (6.8 GB) ---
echo ""
echo "[2/3] ProtT5 + classifier head (6.8 GB)"
if [ "$SKIP_PROT_T5" -eq 1 ]; then
    echo "  [skip] --skip-prot-t5 (not needed for web image)"
else
    HF_BASE="https://huggingface.co/ThorbenF/prot_t5_xl_uniref50_full_v2/resolve/main"
    PT5="${DEST}/prot_t5_xl_uniref50_full_v2"
    for f in model.safetensors cpt.pth spiece.model config.json \
             tokenizer_config.json special_tokens_map.json added_tokens.json; do
        fetch "${HF_BASE}/${f}" "${PT5}/${f}" "prot_t5/${f}"
    done
fi

# --- 3. aSAM / SAM2 weights (1.0 GB) ---
echo ""
echo "[3/3] aSAM / SAM2 weights (1.0 GB)"
SAM_DIR="${DEST}/sam2_weights"
if [ -f "${SAM_DIR}/weights/mdcath_1.0/nn.gen.pt" ]; then
    echo "  [skip] sam2_weights (already extracted)"
else
    ZIP="${SAM_DIR}/mdcath_1.0.zip"
    fetch "https://github.com/giacomo-janson/sam2/releases/download/data-1.0/mdcath_1.0.zip" \
          "$ZIP" "mdcath_1.0.zip"
    echo "  [unzipping] sam2_weights..."
    mkdir -p "${SAM_DIR}/weights"
    unzip -q -o "$ZIP" -d "${SAM_DIR}/weights"
    rm -f "$ZIP"
fi

# Make everything readable (needed for container builds).
chmod -R u+rwX,go+rX,go-w "$DEST" 2>/dev/null || true

echo ""
echo "========================================="
echo " Done! All weights downloaded."
echo ""
echo " Next: build the container images:"
echo "   cd ${SCRIPT_DIR}"
echo "   sudo apptainer build pipeline_web.sif pipeline_web.def"
if [ "$SKIP_PROT_T5" -eq 0 ]; then
echo "   sudo apptainer build pipeline_local.sif pipeline_local.def"
fi
echo "========================================="
