#!/bin/bash
#SBATCH --job-name=prosit_pipeline
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --time=6:00:00
#SBATCH --output=prosit_pipeline_%j.log

# -------------------------------------------------------
# PROSIT TMT Pipeline driver
# -------------------------------------------------------
# Usage:
#   sbatch run_prosit_pipeline.sh                        # full run
#   sbatch run_prosit_pipeline.sh --test                 # 20-peptide test
#   sbatch run_prosit_pipeline.sh --phase2 --no-plots    # SA only, skip plots
#
# Edit WORKSPACE / INPUT_CSV / OUTPUT_DIR for your project.
# -------------------------------------------------------

set -euo pipefail

# --- User configuration ---------------------------------------------------
WORKSPACE="${PROSIT_WORKSPACE:-$PWD}"
INPUT_CSV="${PROSIT_INPUT_CSV:-$WORKSPACE/cleaned_tryptic_peptides_detailed_under_151aa.csv}"
OUTPUT_DIR="${PROSIT_OUTPUT_DIR:-$WORKSPACE/prosit}"
CONDA_ENV="${PROSIT_CONDA_ENV:-prosit}"
# --------------------------------------------------------------------------

cd "$WORKSPACE"

echo "=== PROSIT Pipeline ==="
echo "Job ID:    ${SLURM_JOB_ID:-local}"
echo "Node:      $(hostname)"
echo "CPUs:      ${SLURM_CPUS_PER_TASK:-1}"
echo "Date:      $(date)"
echo "Workspace: $WORKSPACE"
echo "Input:     $INPUT_CSV"
echo "Output:    $OUTPUT_DIR"
echo "========================"

# Activate conda (edit path for local installs)
if [[ -f /apps/conda/miniforge3/24.3.0/etc/profile.d/conda.sh ]]; then
    source /apps/conda/miniforge3/24.3.0/etc/profile.d/conda.sh
fi
conda activate "$CONDA_ENV"

echo "Python: $(which python)"
echo "Args:   $*"
echo "========================"

python "$(dirname "$0")/prosit_pipeline.py" \
    --workspace  "$WORKSPACE" \
    --input-csv  "$INPUT_CSV" \
    --output-dir "$OUTPUT_DIR" \
    --workers    "${SLURM_CPUS_PER_TASK:-4}" \
    "$@"

echo ""
echo "=== Finished at $(date) ==="
