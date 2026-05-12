#!/bin/bash

# === Set paths ===
ANNOTATOR_SCRIPT="scripts/Annotator.py"
SMORF_GTF="data/shortstop/final_microproteins.gtf"
ENSEMBL_GTF="data/ensembl/ensembl_hg38_filtered.gtf"
OUTPUT_DIR="results/espresso"
INTERSECT_OUT="${OUTPUT_DIR}/intersect.gtf"
NON_INTERSECT_OUT="${OUTPUT_DIR}/nonintersect.gtf"
OUTPUT_FILE="${OUTPUT_DIR}/annotated_smorf_types.txt"

# === Create output directory if it doesn't exist ===
mkdir -p "${OUTPUT_DIR}"

# === Run Annotator ===
python "${ANNOTATOR_SCRIPT}" smorf_types \
  --smorf_gtf "${SMORF_GTF}" \
  --ensembl_gtf "${ENSEMBL_GTF}" \
  --outdir "${OUTPUT_DIR}" \
  --intersect_output "${INTERSECT_OUT}" \
  --non_intersect_output "${NON_INTERSECT_OUT}" \
  --output_file "${OUTPUT_FILE}"