#!/usr/bin/env bash
#
# Publish the UCSC browser tracks to the Hugging Face *dataset*
# brmiller/brain-microprotein-atlas, under browser_tracks/.
#
# Runs the three steps end to end:
#   1. Authenticate with Hugging Face (skipped if already logged in)
#   2. Check the track directory exists and report what is on disk
#   3. Run push_browser_tracks.py (diff -> upload -> verify by hash)
#
# This is the step that actually makes a track rebuild visible in UCSC. The
# tracks the browser loads live in the dataset repo, NOT in the Space and NOT in
# the shared Box folder — bigDataUrl points at the dataset. A rebuild that stops
# at Box leaves the browser serving the previous build while every local check
# reports success.
#
# The Python script does its own preflight (bigDataUrl resolution, stale-legend
# warning, post-upload hash verification) and uploads only files whose published
# copy differs, so re-running it is cheap and safe.
#
# Usage (from anywhere) — all arguments are passed through:
#   Results/push_browser_tracks.sh                    # changed files only
#   Results/push_browser_tracks.sh --dry-run          # report, upload nothing
#   Results/push_browser_tracks.sh -m "message"
#   Results/push_browser_tracks.sh --all              # ignore the diff
#
# Non-interactive login: export HF_TOKEN=hf_xxx before running.

set -euo pipefail

DATASET_REPO="brmiller/brain-microprotein-atlas"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"   # Github/Results
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"                    # Github/
TRACKS_DIR="$REPO_ROOT/Code/data/browser_tracks"

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

echo "==> Dataset: $DATASET_REPO / browser_tracks/"
echo "==> Source:  $TRACKS_DIR"
echo

# ---------------------------------------------------------------------------
# Step 1 — Authenticate with Hugging Face
# ---------------------------------------------------------------------------
# huggingface_hub >= 1.0 ships the `hf` CLI (`hf auth ...`); older versions use
# `huggingface-cli`. Pick whichever is available.
echo "==> [1/3] Hugging Face authentication"
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
# Step 2 — Report what is on disk
# ---------------------------------------------------------------------------
# The Python script skips absent files with a warning rather than failing (not
# every rebuild produces every artifact), so show the inventory up front —
# that warning is easy to miss in the middle of its output.
echo "==> [2/3] Track directory"
if [[ ! -d "$TRACKS_DIR" ]]; then
    echo "    ERROR: track directory not found: $TRACKS_DIR" >&2
    echo "    Note: brain_tracks.txt moved into Code/data/browser_tracks/;" >&2
    echo "    it is no longer at the Code/data/ top level." >&2
    exit 1
fi

n=$(find "$TRACKS_DIR" -maxdepth 1 -type f | wc -l | tr -d ' ')
sz=$(du -sh "$TRACKS_DIR" 2>/dev/null | cut -f1)
echo "    $n files, $sz"

if [[ ! -f "$TRACKS_DIR/brain_tracks.txt" ]]; then
    echo "    ERROR: brain_tracks.txt is missing — that is the hub config UCSC reads." >&2
    exit 1
fi
echo

# ---------------------------------------------------------------------------
# Step 3 — Push
# ---------------------------------------------------------------------------
# push_browser_tracks.py compares local against published by sha256, uploads
# only what differs in a single commit, then re-downloads and re-hashes each
# uploaded file. It exits non-zero if any verification fails.
echo "==> [3/3] Pushing browser tracks"
python "$SCRIPT_DIR/push_browser_tracks.py" "$@"
echo

echo "Done."
echo "  Dataset: https://huggingface.co/datasets/${DATASET_REPO}/tree/main/browser_tracks"
