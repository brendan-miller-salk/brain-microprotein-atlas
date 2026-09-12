#!/usr/bin/env bash
#
# Upload the three figure libraries to the Hugging Face *dataset*
# brmiller/brain-microprotein-atlas.
#
# Runs the four steps end to end:
#   1. Authenticate with Hugging Face (skipped if already logged in)
#   2. Check the three figure directories are present and report their size
#   3. Pick a staging parent OUTSIDE the Box-synced tree
#   4. Run upload_to_hf.py (stage -> upload -> restore -> clean up)
#
# Step 3 is the reason this wrapper exists. upload_large_folder writes resumable
# bookkeeping (.cache/.huggingface/ with a .lock per file) *inside* folder_path.
# With ~28,000 files that is thousands of lock files, and if folder_path sits in
# the Box tree, Box syncs every one of them to the cloud. upload_to_hf.py moves
# the directories to $HF_STAGE_DIR (default ~/Desktop) to dodge this; the check
# here fails loudly if that path is itself inside a synced folder, which the
# Python script only warns about.
#
# This is NOT part of the dashboard deploy. launch_hf_dashboard.sh ships the app
# code to the Space and deliberately excludes these images; the deployed
# dashboard streams them from the dataset this script populates. Re-run this only
# when the figure libraries themselves are regenerated — it is a multi-GB upload.
#
# Usage (from anywhere):
#   Results/upload_figures_to_hf.sh
#   Results/upload_figures_to_hf.sh --yes             # skip the confirmation
#   HF_STAGE_DIR=/tmp/stage Results/upload_figures_to_hf.sh
#
# Non-interactive login: export HF_TOKEN=hf_xxx before running.

set -euo pipefail

DATASET_REPO="brmiller/brain-microprotein-atlas"
DIRS=(mirror_plots expression_profiles smorf_cartoon_figures)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"   # Github/Results

ASSUME_YES=0
for arg in "$@"; do
    case "$arg" in
        -y|--yes) ASSUME_YES=1 ;;
        -h|--help) sed -n '3,/^set -euo/p' "${BASH_SOURCE[0]}" | sed -e 's/^#//' -e 's/^ //' -e '$d'; exit 0 ;;
        *) echo "Unknown option: $arg (try --help)" >&2; exit 1 ;;
    esac
done

# Optionally activate a conda env that has huggingface_hub installed.
if [[ -n "${DASHBOARD_CONDA_ENV:-}" ]]; then
    CONDA_SH="${CONDA_PROFILE:-$HOME/miniconda/etc/profile.d/conda.sh}"
    if [[ -f "$CONDA_SH" ]]; then
        # shellcheck disable=SC1090
        source "$CONDA_SH"
        conda activate "$DASHBOARD_CONDA_ENV"
    else
        echo "Warning: conda profile not found at $CONDA_SH; using current environment." >&2
    fi
fi

echo "==> Dataset:  $DATASET_REPO"
echo "==> Source:   $SCRIPT_DIR"
echo

# ---------------------------------------------------------------------------
# Step 1 — Authenticate with Hugging Face
# ---------------------------------------------------------------------------
# huggingface_hub >= 1.0 ships the `hf` CLI (`hf auth ...`); older versions use
# `huggingface-cli`. Pick whichever is available.
echo "==> [1/4] Hugging Face authentication"
if command -v hf >/dev/null 2>&1; then
    HF_CLI="hf auth"
elif command -v huggingface-cli >/dev/null 2>&1; then
    HF_CLI="huggingface-cli"
else
    echo "    ERROR: neither 'hf' nor 'huggingface-cli' found on PATH." >&2
    echo "    Install with: pip install -U huggingface_hub" >&2
    exit 1
fi

if [[ -n "${HF_TOKEN:-}" ]]; then
    echo "    Using HF_TOKEN from the environment."
    $HF_CLI login --token "$HF_TOKEN" --add-to-git-credential >/dev/null
elif $HF_CLI whoami >/dev/null 2>&1; then
    echo "    Already logged in as: $($HF_CLI whoami)"
else
    echo "    Not logged in — launching interactive login."
    $HF_CLI login
fi
echo

# ---------------------------------------------------------------------------
# Step 2 — Confirm the payload is on disk
# ---------------------------------------------------------------------------
# upload_to_hf.py only warns on a missing directory and uploads whatever is
# left, which would silently publish a partial set. Fail here instead.
echo "==> [2/4] Checking figure directories"
MISSING=()
TOTAL_FILES=0
for d in "${DIRS[@]}"; do
    if [[ -d "$SCRIPT_DIR/$d" ]]; then
        n=$(find "$SCRIPT_DIR/$d" -type f | wc -l | tr -d ' ')
        sz=$(du -sh "$SCRIPT_DIR/$d" 2>/dev/null | cut -f1)
        TOTAL_FILES=$((TOTAL_FILES + n))
        printf '    %-24s %8s  %6s files\n' "$d" "$sz" "$n"
    else
        MISSING+=("$d")
        printf '    %-24s %8s\n' "$d" "ABSENT"
    fi
done

if (( ${#MISSING[@]} )); then
    echo >&2
    echo "    ERROR: missing figure director$( ((${#MISSING[@]}>1)) && echo ies || echo y ): ${MISSING[*]}" >&2
    echo "    Regenerate them, or pull them from the dataset, before uploading." >&2
    echo "    Uploading now would publish an incomplete set." >&2
    exit 1
fi
echo "    Total: $TOTAL_FILES files"
echo

# ---------------------------------------------------------------------------
# Step 3 — Verify the staging parent is outside the Box-synced tree
# ---------------------------------------------------------------------------
# The staging dir must not be under CloudStorage: upload_large_folder writes a
# .lock per file there, and Box would sync all of them. upload_to_hf.py only
# prints a warning in that case; here it is a hard stop.
echo "==> [3/4] Staging location"
STAGE_PARENT="${HF_STAGE_DIR:-$HOME/Desktop}"
echo "    HF_STAGE_DIR = $STAGE_PARENT"

if [[ ! -d "$STAGE_PARENT" ]]; then
    echo "    ERROR: staging parent does not exist: $STAGE_PARENT" >&2
    echo "    Create it, or set HF_STAGE_DIR to an existing local path." >&2
    exit 1
fi

case "$STAGE_PARENT" in
    *CloudStorage*|*"Box Sync"*|*Dropbox*|*"Google Drive"*)
        echo "    ERROR: staging parent is inside a cloud-synced tree." >&2
        echo "    upload_large_folder writes a .lock per file there — with" >&2
        echo "    $TOTAL_FILES files that means syncing thousands of lock files." >&2
        echo "    Set HF_STAGE_DIR to a local path (e.g. \$HOME/Desktop)." >&2
        exit 1
        ;;
esac

# Same filesystem => the stage/restore moves are instant renames, not copies.
DEV_SRC=$(df -P "$SCRIPT_DIR" | awk 'NR==2 {print $1}')
DEV_STAGE=$(df -P "$STAGE_PARENT" | awk 'NR==2 {print $1}')
if [[ "$DEV_SRC" == "$DEV_STAGE" ]]; then
    echo "    Same filesystem as the figures — staging is a rename, no data copy."
else
    echo "    WARNING: different filesystem from the figures ($DEV_SRC vs $DEV_STAGE)." >&2
    echo "    Staging will COPY ~5 GB instead of renaming. Consider another path." >&2
fi
echo

if (( ! ASSUME_YES )); then
    echo "About to upload $TOTAL_FILES files to $DATASET_REPO (multi-GB, slow)."
    read -r -p "Continue? [y/N] " reply
    [[ "$reply" =~ ^[Yy]$ ]] || { echo "Aborted."; exit 0; }
    echo
fi

# ---------------------------------------------------------------------------
# Step 4 — Upload
# ---------------------------------------------------------------------------
# upload_to_hf.py stages the dirs, uploads, then restores them in a finally
# block — the restore runs even if the upload fails partway.
echo "==> [4/4] Uploading figure libraries"
HF_STAGE_DIR="$STAGE_PARENT" python "$SCRIPT_DIR/upload_to_hf.py"
echo

echo "Done."
echo "  Dataset: https://huggingface.co/datasets/${DATASET_REPO}"
