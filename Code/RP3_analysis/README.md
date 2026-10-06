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
- **Flanking Ribo-seq**: `smorf_flank_psites.py` counts P-sites (same bigwigs,
  all three frames, same strand) in the CDS and in 50/100/250/1000-nt windows
  up- and downstream of every Salk and TrEMBL smORF. Windows walk along the
  spliced host transcript (smORF GTF for Salk; the smORF Rules BED for TrEMBL)
  and stop at its ends. The stop codon always counts as CDS. Frame 0/1/2 %
  is relative to the smORF start and uses only positions outside annotated
  GENCODE CDS, skipping the 6 nt just upstream of the start codon, where the
  initiation P-site peak spills over. Writes `Code/data/smorf_flank_psites.csv`
  (per smORF), `smorf_flank_structure.csv` (exons, CDS and annotated-CDS spans,
  which the dashboard uses to plot the profile live) and
  `smorf_flank_data_dictionary.csv`. Needs the git-ignored bigwigs and the
  933 MB combined GTF, so it is run by hand and not by `run_all_analyses.sh`.

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
# (optional) flanking P-sites; same bigwigs + GTF_and_BED_files/Ensembl_and_Unreviewed_Brain_Microproteins.gtf
python smorf_flank_psites.py       # writes Code/data/smorf_flank_{psites,structure,data_dictionary}.csv (~20 s)
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
- `../data/smorf_flank_psites.csv`, `../data/smorf_flank_structure.csv`,
  `../data/smorf_flank_data_dictionary.csv`
- `../../supplementary/Table5_RP3_Results_summary.csv`
- (Already present) `ribocode_results.{bed,gtf,txt}`,
  `ribocode_results_collapsed.*`, `mapping_groups_rpkm*.txt`,
  `ribocode_results_ORFs_category.pdf` -- copied/curated RiboCode outputs.
