# RP3 (Ribosome Profiling) Analysis

This module summarizes ribosome-profiling evidence (RiboCode / RP3) for the
brain microprotein atlas. The raw RP3 / RiboCode pipeline (read alignment,
P-site calling, ORF prediction) is run on dbGaP-controlled FASTQ outside of
this repository; only the post-RiboCode summarization step lives here.

## Overview
- **Input**: `Code/data/microprotein_master.csv` (master annotation table
  with RP3 / RiboCode evidence columns already merged in).
- **Output**: `Results/RP3/RP3_Results_summary.csv` and `RP3_psORFs.csv` -
  per-microprotein RiboCode evidence (ORF type, RPKM, classification) used
  by the dashboard and the manuscript supplemental tables.
- **Main Script**: `RP3_Results_summary.py`.

The RiboCode reference outputs themselves (BED / GTF / TXT and the RPKM
mapping-group files) are also kept under `Results/RP3/` so the dashboard
can resolve coordinates and per-mapping-group expression without
re-running RP3.

## Usage

```bash
# from this directory
python RP3_Results_summary.py
```

The script is also called automatically by
`bash run_all_analyses.sh --mode=run` from the repo root.

## Dependencies
- Python: `pandas`, `numpy`, `os`.
- Shared filter: `Code/gold_standard_filtering_criteria.py`.

## Outputs
- `../../Results/RP3/RP3_Results_summary.csv`
- `../../Results/RP3/RP3_psORFs.csv`
- (Already present) `ribocode_results.{bed,gtf,txt}`,
  `ribocode_results_collapsed.*`, `mapping_groups_rpkm*.txt`,
  `ribocode_results_ORFs_category.pdf` -- copied/curated RiboCode outputs.
