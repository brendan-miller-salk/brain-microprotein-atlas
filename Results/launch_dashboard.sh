#!/usr/bin/env bash

set -euo pipefail

# Launch the dashboard from this script's directory.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Optionally activate a conda env. Set DASHBOARD_CONDA_ENV to use one; otherwise
# the script runs in whatever environment already has the deps installed
# (e.g. `pip install -r requirements.txt` in a venv, per the README).
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

cd "$SCRIPT_DIR"

if [[ ! -f "microproteins_dashboard.py" ]]; then
    echo "microproteins_dashboard.py not found in $SCRIPT_DIR"
    exit 1
fi

echo "Launching dashboard at http://localhost:8505"

streamlit run microproteins_dashboard.py \
    --server.port 8505 \
    --server.headless true \
    --server.runOnSave true