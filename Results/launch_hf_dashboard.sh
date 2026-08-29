#!/usr/bin/env bash
#
# Deploy the Brain Microprotein Atlas dashboard to its public Hugging Face
# Space and sync the same content to GitHub / Streamlit Cloud.
#
# Runs the five steps end to end:
#   1. Authenticate with Hugging Face (skipped if already logged in)
#   2. Stage the working tree so `git ls-files` sees newly added files
#   3. Push the git-tracked file set to the HF Streamlit Space
#   4. Set the DASHBOARD_PUBLIC=1 Space variable (password-free mirror)
#   5. Re-push the repo to GitHub (Streamlit Cloud auto-redeploys)
#
# Step 2 exists because push_to_hf_space.py enumerates files with `git ls-files`,
# which reads the index. Step 5 (clean_and_push_repo.sh) is what normally stages
# new files, so without this the HF upload silently omits anything added since
# the last deploy — it would only reach the Space one cycle late.
#
# Usage (from anywhere):
#   Results/launch_hf_dashboard.sh
#
# Non-interactive login: export HF_TOKEN=hf_xxx before running.

set -euo pipefail

SPACE_REPO="brmiller/brain-microprotein-atlas-app"

# Resolve directories relative to this script.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"   # Github/Results
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"                    # Github/

# Optionally activate a conda env that has streamlit + huggingface_hub installed.
# Set DASHBOARD_CONDA_ENV to activate one; otherwise the current environment is used.
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

echo "==> Space:     $SPACE_REPO"
echo "==> Repo root: $REPO_ROOT"
echo

# ---------------------------------------------------------------------------
# Step 1 — Authenticate with Hugging Face
# ---------------------------------------------------------------------------
# huggingface_hub >= 1.0 ships the `hf` CLI (`hf auth ...`); older versions use
# `huggingface-cli`. Pick whichever is available.
echo "==> [1/5] Hugging Face authentication"
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
# Step 2 — Stage the working tree so the HF upload sees new files
# ---------------------------------------------------------------------------
# push_to_hf_space.py lists files with `git ls-files`, which reads the index —
# untracked files are invisible to it. Staging here (rather than relying on the
# `git add` inside clean_and_push_repo.sh at step 5) is what makes a file added
# this session actually reach the Space this cycle. Harmless to re-run: step 5
# wipes .git and re-stages from scratch anyway.
echo "==> [2/5] Staging working tree for the HF file enumeration"
git -C "$REPO_ROOT" add -A
echo "    Staged $(git -C "$REPO_ROOT" ls-files | wc -l | tr -d ' ') tracked files."
echo

# ---------------------------------------------------------------------------
# Step 3 — Deploy to the Hugging Face Space
# ---------------------------------------------------------------------------
echo "==> [3/5] Pushing dashboard to the Hugging Face Space"
python "$SCRIPT_DIR/push_to_hf_space.py"
echo

# ---------------------------------------------------------------------------
# Step 4 — Make the Space public (bypass the password gate)
# ---------------------------------------------------------------------------
echo "==> [4/5] Setting DASHBOARD_PUBLIC=1 Space variable"
python - "$SPACE_REPO" <<'PY'
import sys
from huggingface_hub import HfApi

repo_id = sys.argv[1]
HfApi().add_space_variable(repo_id=repo_id, key="DASHBOARD_PUBLIC", value="1")
print(f"    DASHBOARD_PUBLIC=1 set on {repo_id}")
PY
echo

# ---------------------------------------------------------------------------
# Step 5 — Sync to GitHub (Streamlit Cloud redeploys from this)
# ---------------------------------------------------------------------------
echo "==> [5/5] Re-pushing repo to GitHub"
( cd "$REPO_ROOT" && ./clean_and_push_repo.sh )
echo

echo "Done."
echo "  HF Space:        https://huggingface.co/spaces/${SPACE_REPO}"
echo "  Streamlit Cloud: https://brain-microprotein-atlas.streamlit.app/"
