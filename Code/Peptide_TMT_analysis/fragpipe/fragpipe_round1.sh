#!/bin/bash
# -----------------------------------------------------------------------------
# FragPipe TMT (round 1, 10-plex) array job
# -----------------------------------------------------------------------------
# Generalized SLURM driver for the round-1 TMT-ROSMAP FragPipe searches.
# Each array task processes one batch directory `${BATCH_PREFIX}${TASK_ID}`
# under $BASE_DIR. Override any of the variables in the "User configuration"
# block below via the environment, e.g.:
#
#   BASE_DIR=/path/to/tmt_rosmap \
#   FRAGPIPE_BIN=/opt/fragpipe/bin \
#   WORKFLOW_FILE=/path/to/tmt_shortstop_appended_rescored.workflow \
#   sbatch --array=1-50 fragpipe_round1.sh
#
# Edit the #SBATCH lines (account / partition) for your cluster.
# -----------------------------------------------------------------------------
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=5:00:00
#SBATCH --mem-per-cpu=24GB
#SBATCH --job-name=fragpipe_tmt_r1
#SBATCH --mail-type=END
#SBATCH --output=fragpipe_r1_%A_%a.log

set -euo pipefail

# --- User configuration -------------------------------------------------------
BASE_DIR="${BASE_DIR:-$PWD}"
BATCH_PREFIX="${BATCH_PREFIX:-b}"
TASK_ID="${SLURM_ARRAY_TASK_ID:-${TASK_ID:-1}}"

FRAGPIPE_BIN="${FRAGPIPE_BIN:-}"                  # e.g. /opt/fragpipe/bin
FRAGPIPE_TOOLS_DIR="${FRAGPIPE_TOOLS_DIR:-}"      # --config-tools-folder
FRAGPIPE_PYTHON="${FRAGPIPE_PYTHON:-$(command -v python3 || true)}"

WORKFLOW_FILE="${WORKFLOW_FILE:-${BASE_DIR}/fragpipe_processing_files/tmt_shortstop_appended_rescored.workflow}"
RESULTS_SUBDIR="${RESULTS_SUBDIR:-shortstop_proteogenomics_appended_results_cpm05}"

JAVA_MEM="${JAVA_MEM:-150G}"
START_DELAY_PER_TASK="${START_DELAY_PER_TASK:-5}"  # seconds * TASK_ID
# -----------------------------------------------------------------------------

DIR="${BASE_DIR}/${BATCH_PREFIX}${TASK_ID}"
OUTPUT_FILE="${DIR}/${BATCH_PREFIX}${TASK_ID}_tmt.out"

# Stagger array task starts to avoid IO/license storms.
DELAY=$(( TASK_ID * START_DELAY_PER_TASK ))
echo "Delaying start by ${DELAY} seconds"
sleep "${DELAY}"

# Load Java (cluster-specific; comment out if not using environment modules).
if command -v module >/dev/null 2>&1; then
    module load openjdk || true
fi

if [[ -n "${FRAGPIPE_BIN}" ]]; then
    export PATH="${FRAGPIPE_BIN}:${PATH}"
fi

export JAVA_OPTS="-Xmx${JAVA_MEM}"
echo "Capping Java memory to: ${JAVA_MEM}"

# 10-plex TMT annotation (round 1 channel layout).
ANNOTATION_CONTENT=$(cat << EOF
126     ${BATCH_PREFIX}${TASK_ID}.126
127C    ${BATCH_PREFIX}${TASK_ID}.127C
127N    ${BATCH_PREFIX}${TASK_ID}.127N
128C    ${BATCH_PREFIX}${TASK_ID}.128C
128N    ${BATCH_PREFIX}${TASK_ID}.128N
129C    ${BATCH_PREFIX}${TASK_ID}.129C
129N    ${BATCH_PREFIX}${TASK_ID}.129N
130C    ${BATCH_PREFIX}${TASK_ID}.130C
130N    ${BATCH_PREFIX}${TASK_ID}.130N
131N    ${BATCH_PREFIX}${TASK_ID}.131
EOF
)

if [[ ! -d "${DIR}" ]]; then
    echo "Directory ${DIR} does not exist." | tee "${OUTPUT_FILE}"
    exit 1
fi

echo "${ANNOTATION_CONTENT}" > "${DIR}/annotation.txt"

MANIFEST="${DIR}/${BATCH_PREFIX}${TASK_ID}.manifest"
rm -f "${MANIFEST}"
for file in "${DIR}"/*${BATCH_PREFIX}${TASK_ID}*mzML; do
    [[ -e "${file}" ]] || continue
    printf '%s\tDDA\n' "${file}" >> "${MANIFEST}"
done

mkdir -p "${DIR}/${RESULTS_SUBDIR}"
if [[ -n "$(ls -A "${DIR}/${RESULTS_SUBDIR}" 2>/dev/null)" ]]; then
    rm -rf "${DIR}/${RESULTS_SUBDIR}"/*
fi

cp "${WORKFLOW_FILE}" "${DIR}/"
LOCAL_WORKFLOW="${DIR}/$(basename "${WORKFLOW_FILE}")"

FRAGPIPE_ARGS=(
    --headless
    --workflow  "${LOCAL_WORKFLOW}"
    --manifest  "${MANIFEST}"
    --workdir   "${DIR}/${RESULTS_SUBDIR}"
)
[[ -n "${FRAGPIPE_TOOLS_DIR}" ]] && FRAGPIPE_ARGS+=( --config-tools-folder "${FRAGPIPE_TOOLS_DIR}" )
[[ -n "${FRAGPIPE_PYTHON}"    ]] && FRAGPIPE_ARGS+=( --config-python       "${FRAGPIPE_PYTHON}"    )

fragpipe "${FRAGPIPE_ARGS[@]}"
