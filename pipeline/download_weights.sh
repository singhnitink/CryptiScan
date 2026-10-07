#!/bin/bash
# Download CryptiScan model weights (~15 GB total).
# Run once after cloning:  bash download_weights.sh
# Use --skip-prot-t5 if you only need the web image.
# Requires: wget, unzip, sha256sum
#
# Every file is checked against the SHA-256 below and only counts as downloaded
# once it matches, so a truncated or corrupted file is never baked into an image.
# Re-running is safe: verified files are skipped, interrupted downloads resume.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEST="${SCRIPT_DIR}/required_files_weights"

ESM_URL="https://dl.fbaipublicfiles.com/fair-esm/models/esm1v_t33_650M_UR90S_1.pt"
ESM_SHA256="9519ee60f1cddad3c101afb1f42612499e188534969c3f682e94850870f70433"

# Pinned to a revision so the checksums below stay valid if the repo changes.
PT5_URL="https://huggingface.co/ThorbenF/prot_t5_xl_uniref50_full_v2/resolve/c12bb96dfec9fe25bbb287ec0896680371b762d4"
PT5_FILES=(model.safetensors cpt.pth spiece.model config.json
           tokenizer_config.json special_tokens_map.json added_tokens.json)
declare -A PT5_SHA256=(
    [model.safetensors]="8a9fbcd82de95fc796a9474b0355b2a20952fe5329daa7c5ad2912df3462d668"
    [cpt.pth]="64d44ac3e046adf07214902ae0347bd6eaef0548000824e5674351de05cbced8"
    [spiece.model]="74da7b4afcde53faa570114b530c726135bdfcdb813dec3abfb27f9d44db7324"
    [config.json]="97a252d48d11df51a69fd17a2fd29ec762a7b7f007021bdd1059ee9615614625"
    [tokenizer_config.json]="007a5d1675f05df316c83813d5324904a3841db75a3e0ac684ec76ec1f4ef925"
    [special_tokens_map.json]="7a1985a994c41886db38c719d2a3d2f40606663cc19d7c5d6a85d349320e06d2"
    [added_tokens.json]="3893cf0a0185c9f84d45e39ace998d9a3f1103714361b61b715237123671a99f"
)

SAM_URL="https://github.com/giacomo-janson/sam2/releases/download/data-1.0/mdcath_1.0.zip"
SAM_ZIP_SHA256="c17593a10a832cffa2a3fbeb40032fdd174ccc25c1c12ac609457e5fb7ff1abb"
SAM_FILES=(nn.gen.pt nn.enc.pt nn.dec.pt enc_std_scaler.pt)
declare -A SAM_SHA256=(
    [nn.gen.pt]="0207bd155f5fdc1c95dea8629637f7662c4baea4934830c98429b429ca41d3c4"
    [nn.enc.pt]="23d2d641c972d463663173e1290b9a1e223484ac51886af796d7eff0bd042a3d"
    [nn.dec.pt]="ce9b7f6a7eab9607048ec3efbfde09e5906338158252d067637564599fbbfb29"
    [enc_std_scaler.pt]="3d4604c3a992c41f51f045a34ce0bede68eaed3aeb9dae87cb0922c28d0ceac0"
)

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

for tool in wget unzip sha256sum; do
    command -v "$tool" &>/dev/null || { echo "ERROR: '$tool' is not installed."; exit 1; }
done

mkdir -p "$DEST"

# True if the file exists and matches the expected SHA-256.
sha_ok() {
    [ -f "$1" ] && echo "$2  $1" | sha256sum --check --status
}

# Download a file until it matches its SHA-256. A partial file left by an
# interrupted run is resumed; a file that downloads completely but fails the
# check is deleted and fetched again from scratch.
fetch() {
    local url="$1" out="$2" label="$3" sha="$4"
    mkdir -p "$(dirname "$out")"

    if [ -f "$out" ]; then
        echo "  [checking] $label"
        if sha_ok "$out" "$sha"; then
            echo "  [skip] $label (already downloaded, checksum OK)"
            return 0
        fi
    fi

    local attempt wait_time
    for attempt in 1 2 3; do
        if [ -s "$out" ]; then
            echo "  [resuming] $label"
        else
            echo "  [downloading] $label"
        fi
        if wget --continue --show-progress --max-redirect=10 -O "$out" "$url"; then
            if sha_ok "$out" "$sha"; then
                # Brief pause between downloads to avoid rate-limiting (HuggingFace CDN).
                sleep 3
                return 0
            fi
            echo "  [checksum mismatch] $label, deleting and downloading again"
            rm -f "$out"
        fi
        if [ "$attempt" -lt 3 ]; then
            wait_time=$((attempt * 20))
            echo "  [retry $attempt/3] waiting ${wait_time}s before retrying $label..."
            sleep "$wait_time"
        fi
    done

    echo "  ERROR: could not download a verified copy of $label after 3 attempts."
    echo "  URL: $url"
    if [ -s "$out" ]; then
        echo "  The partial file was kept; re-run this script to resume."
    fi
    exit 1
}

echo "CryptiScan: downloading model weights"
echo "Destination: $DEST"

# 1. ESM-1v (7.3 GB)
echo ""
echo "[1/3] ESM-1v (7.3 GB)"
fetch "$ESM_URL" "${DEST}/esm1v_t33_650M_UR90S_1.pt" "esm1v_t33_650M_UR90S_1.pt" "$ESM_SHA256"

# 2. ProtT5 + classifier head (6.8 GB)
echo ""
echo "[2/3] ProtT5 + classifier head (6.8 GB)"
if [ "$SKIP_PROT_T5" -eq 1 ]; then
    echo "  [skip] --skip-prot-t5 (not needed for web image)"
else
    PT5="${DEST}/prot_t5_xl_uniref50_full_v2"
    for f in "${PT5_FILES[@]}"; do
        fetch "${PT5_URL}/${f}" "${PT5}/${f}" "prot_t5/${f}" "${PT5_SHA256[$f]}"
    done
fi

# 3. aSAM / SAM2 weights (1.0 GB)
echo ""
echo "[3/3] aSAM / SAM2 weights (1.0 GB)"
SAM_DIR="${DEST}/sam2_weights"
SAM_WEIGHTS="${SAM_DIR}/weights/mdcath_1.0"
sam_ok() {
    local f
    for f in "${SAM_FILES[@]}"; do
        sha_ok "${SAM_WEIGHTS}/${f}" "${SAM_SHA256[$f]}" || return 1
    done
}
if [ -d "$SAM_WEIGHTS" ] && { echo "  [checking] sam2_weights"; sam_ok; }; then
    echo "  [skip] sam2_weights (already extracted, checksums OK)"
else
    ZIP="${SAM_DIR}/mdcath_1.0.zip"
    fetch "$SAM_URL" "$ZIP" "mdcath_1.0.zip" "$SAM_ZIP_SHA256"
    echo "  [unzipping] sam2_weights..."
    mkdir -p "${SAM_DIR}/weights"
    unzip -q -o "$ZIP" -d "${SAM_DIR}/weights"
    if ! sam_ok; then
        echo "  ERROR: extracted sam2_weights do not match their checksums."
        exit 1
    fi
    rm -f "$ZIP"
fi

# Make everything readable (needed for container builds).
chmod -R u+rwX,go+rX,go-w "$DEST" 2>/dev/null || true

echo ""
echo "All weights downloaded and verified."
echo ""
echo "Next: build the container images:"
echo "  cd ${SCRIPT_DIR}"
echo "  sudo apptainer build pipeline_web.sif pipeline_web.def"
if [ "$SKIP_PROT_T5" -eq 0 ]; then
echo "  sudo apptainer build pipeline_local.sif pipeline_local.def"
fi
