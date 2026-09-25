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
- **P-site frames**: `Code/data/psite_frame_counts_per_orf.tsv` holds P-site counts per
  codon position for every CDS in the combined GENCODE v43 + microprotein GTF
  (41 adult Ribo-seq libraries; made by `psite_frame_counts_per_orf.py` from
  the `Code/data/browser_tracks` bigwigs). `Psite_frame_mapping.py` attaches
  `Psites_pct_frame0/1/2` and `Psites_frame0_RPKM` to each microprotein:
  unreviewed smORFs by sequence; TrEMBL by the coordinates in
  `GTF_and_BED_files/Unreviewed_Brain_Microproteins_mapping_coordinates_to_sequences.tsv`
  matched to a GENCODE CDS span; Swiss-Prot by gene name. TrEMBL and
  Swiss-Prot also require codon count = protein length.

The RiboCode reference outputs themselves (BED / GTF / TXT and the RPKM
mapping-group files) are also kept under `Results/RP3/` so the dashboard
can resolve coordinates and per-mapping-group expression without
re-running RP3.

## Usage

```bash
# from this directory
# (optional) regenerate the P-site table; needs the adult_psite bigwigs in
# Code/data/browser_tracks/ (Hugging Face dataset) and pyBigWig
python psite_frame_counts_per_orf.py --total-psites 116287730
python Psite_frame_mapping.py      # writes Results/RP3/Psite_frame_by_sequence.csv
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
- `../../Results/RP3/Psite_frame_by_sequence.csv`
- `../../supplementary/Table5_RP3_Results_summary.csv`
- (Already present) `ribocode_results.{bed,gtf,txt}`,
  `ribocode_results_collapsed.*`, `mapping_groups_rpkm*.txt`,
  `ribocode_results_ORFs_category.pdf` -- copied/curated RiboCode outputs.
