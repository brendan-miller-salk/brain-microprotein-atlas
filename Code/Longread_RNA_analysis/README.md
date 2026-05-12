# Long-read RNA Analysis

Processes long-read (Oxford Nanopore) brain RNA-seq from the
ROSMAP / Heberle et al. cohort using the ESPRESSO transcript-quantitation
pipeline followed by DESeq2 differential expression. Raw FASTQ and `.esp`
abundance files are not redistributed here (see
`../../DATA_AVAILABILITY.md`); only the summary script is run by
`run_all_analyses.sh`.

## Layout
```
Longread_RNA_analysis/
├── README.md
├── Long-Read_Transcriptomics_Results_summary.py
└── ESPRESSO_data_processing_scripts/
    ├── convert_ESPERSSO_to_CPM.sh           # bash wrapper
    ├── convert_ESPRESSO_to_CPM_and_filter.py # ESPRESSO -> CPM matrix + filter
    └── deseq_brain_espresso.r                # DESeq2 + PCA / heatmaps / volcano
```

## Pipeline

1. **ESPRESSO -> CPM** (`convert_ESPRESSO_to_CPM_and_filter.py`,
   `convert_ESPERSSO_to_CPM.sh`): converts raw ESPRESSO counts to CPM and
   applies the median-CPM filter that produces
   `Code/data/brain_espresso_medianCPM05.csv`.
2. **DESeq2** (`deseq_brain_espresso.r`): performs differential expression
   on the filtered CPM matrix using `Code/data/nanopore_metadata.csv` for
   sample design; emits PCA, heatmap, volcano, and a results table.
3. **Summary** (`Long-Read_Transcriptomics_Results_summary.py`): collates
   the DESeq2 results into
   `../../Results/Transcriptomics/Long-Read_Transcriptomics_Results_summary.csv`,
   keyed against `Code/data/microprotein_master.csv`.

## Usage

```bash
# (Optional - requires raw ESPRESSO + nanopore data)
bash ESPRESSO_data_processing_scripts/convert_ESPERSSO_to_CPM.sh
Rscript ESPRESSO_data_processing_scripts/deseq_brain_espresso.r

# Summary (works against shipped data)
python Long-Read_Transcriptomics_Results_summary.py

# Or as part of the full summary pipeline:
bash ../../run_all_analyses.sh --mode=run
```

## Dependencies
- Python: `pandas`, `numpy`, `argparse`.
- R: `DESeq2`, `dplyr`, `tibble`, `ggplot2`, `biomaRt`.

## Inputs
- ESPRESSO `*_abundance.esp` files (e.g.
  `ESPRESSO_syn52047893_hg38NCBI_N2_R0_abundance.esp`) - external.
- `Code/data/nanopore_metadata.csv` (shipped).
- `Code/data/brain_espresso_medianCPM05.csv` (shipped).
- `Code/data/microprotein_master.csv` (shipped).

## Outputs
- `../../Results/Transcriptomics/Long-Read_Transcriptomics_Results_summary.csv`
- (When raw pipeline is run) CPM matrices + DESeq2 results, PCA / heatmap /
  volcano figures.
