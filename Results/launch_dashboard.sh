#!/usr/bin/env bash

set -euo pipefail

# Launch the dashboard from this script's directory.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

source /Users/brendanmiller/miniconda/etc/profile.d/conda.sh
conda activate github

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