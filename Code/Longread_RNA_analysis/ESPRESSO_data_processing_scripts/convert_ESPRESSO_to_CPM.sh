#!/bin/bash

# Inputs are author-local (controlled-access ESPRESSO output; not shipped in the
# repo). Override these via environment variables, e.g.:
#   ESPRESSO_ABUNDANCE=/path/to/abundance.esp \
#   ESPRESSO_OUTPUT=/path/to/out/prefix bash convert_ESPRESSO_to_CPM.sh
ESPRESSO_ABUNDANCE="${ESPRESSO_ABUNDANCE:-/Users/brendanmiller/Library/CloudStorage/Box-Box/brain_smorfs/espresso/ESPRESSO_syn52047893_hg38NCBI_N2_R0_abundance.esp}"
ESPRESSO_OUTPUT="${ESPRESSO_OUTPUT:-/Users/brendanmiller/Library/CloudStorage/Box-Box/brain_smorfs/espresso/cleaned_files/brain_espresso}"

# Run the Python script with specified arguments
python convert_ESPRESSO_to_CPM_and_filter.py \
    --abundance-esp "$ESPRESSO_ABUNDANCE" \
    --output-path "$ESPRESSO_OUTPUT"