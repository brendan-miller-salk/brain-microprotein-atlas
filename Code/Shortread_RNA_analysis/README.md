# Short-read RNA Analysis

This module processes short-read RNA-seq data (ROSMAP DLPFC and MSBB) for
differential expression and main-ORF / smORF co-expression analyses.

## Overview
- **Input**: RNA-seq count matrices generated with FeatureCounts on BAMs
  downloaded from Synapse (ROSMAP, MSBB).
- **Output**: DESeq2 differential expression results, batch-corrected VST
  expression, per-pair main-ORF / smORF co-expression statistics and figures.
- **Summary script**: `Short-Read_Transcriptomics_Results_summary.py`.

## Directory layout
```
Shortread_RNA_analysis/
├── README.md
├── Short-Read_Transcriptomics_Results_summary.py
└── Shortread_deseq_processing_scripts/
    ├── RNA_differential_expression.R
    ├── Main_smORF_LRT_analysis.R
    └── Main_smORF_coexpression_analysis.R
```

## Analysis pipeline

### 1. Differential expression — `RNA_differential_expression.R`
Generalizable DESeq2 workflow. Edit the `CONFIG` block at the top of the
script to point at either ROSMAP or MSBB inputs (counts matrix, sample
metadata, design covariates, group column / levels, contrast).

Performs:
- Counts + metadata loading and harmonisation.
- Diagnosis classification (ROSMAP: CERAD + Braak + MMSE/cogdx;
  MSBB: CERAD + Braak + CDR).
- DESeq2 fit with `~ covariates + group` design.
- VST and `limma::removeBatchEffect`-corrected VST.
- AD vs. Control results table with ENSG → SYMBOL mapping.

Outputs (written under `OUT_DIR`, prefixed with `tolower(DATASET)`):
- `<prefix>_dds.rds`
- `<prefix>_vsd.rds`
- `<prefix>_vsd_corrected.rds`  ← consumed by the two scripts below
- `<prefix>_deseq2_<A>_vs_<B>_results.csv`

### 2. Main-ORF / smORF LRT — `Main_smORF_LRT_analysis.R`
For every smORF / main-ORF pair in the MASTER annotation table, fits per-pair
linear models on the batch-corrected VST values and tests whether the
relationship between the smORF and its main ORF differs by diagnosis group.

For each pair the script reports:
- Pearson r in each group and Δr.
- Fisher Z test of the two correlations (`psych::paired.r`).
- LRT additive model (`smORF ~ main + group`) F / p / β / SE for the AD term.
- LRT interaction model (`smORF ~ main * group`) F / p / β / SE for the
  `main:group` term.

Output: `<OUT_DIR>/<prefix>_main_smorf_LRT_results.csv`.

### 3. Main-ORF / smORF co-expression — `Main_smORF_coexpression_analysis.R`
Generates a per-pair triptych figure
`[scatter w/ group-coloured fits | main-ORF violin | smORF violin]`
and routes each PDF to `coupled/` or `non_coupled/` based on whether
`|Δr| > DELTA_R_ROUTE_THRESHOLD` (default 0.1).

Outputs (under `OUT_ROOT`):
- `coupled/*.pdf`, `non_coupled/*.pdf` — per-pair triptych figures
- `<prefix>_main_smorf_coexpression_results.csv` — Pearson r per group, Δr, route
- `<prefix>_global_lines.rds` — per-pair fit-line data for rebuilding the
  global "landscape" figure without re-running all pairs

### 4. Summary — `Short-Read_Transcriptomics_Results_summary.py`
Collates the DESeq2 result tables across cohorts into the manuscript
supplemental tables.

## Usage

All three R scripts use repo-relative paths configurable via environment
variables (or by editing the `CONFIG` block at the top of each script):

| Variable                     | Default                                | Used by                       |
|------------------------------|----------------------------------------|-------------------------------|
| `AD_PAPER_DATA_DIR`          | `./data`                               | All three (counts, metadata, MASTER table) |
| `AD_PAPER_OUT_DIR`           | `./outputs/deseq`                      | DESeq output dir              |
| `AD_PAPER_DESEQ_DIR`         | `./outputs/deseq`                      | LRT + co-expression inputs    |
| `AD_PAPER_LRT_OUT_DIR`       | `./outputs/main_smorf_LRT`             | LRT output dir                |
| `AD_PAPER_COEXPR_OUT_DIR`    | `./outputs/main_smorf_coexpression`    | Co-expression output dir      |

```bash
# 1. Differential expression (edit CONFIG to switch ROSMAP <-> MSBB)
Rscript Shortread_deseq_processing_scripts/RNA_differential_expression.R

# 2. Main-ORF / smORF LRT analysis
Rscript Shortread_deseq_processing_scripts/Main_smORF_LRT_analysis.R

# 3. Main-ORF / smORF co-expression triptychs + routing
Rscript Shortread_deseq_processing_scripts/Main_smORF_coexpression_analysis.R

# 4. Final summary tables
python Short-Read_Transcriptomics_Results_summary.py
```

## Dependencies
- **R**: DESeq2, dplyr, tibble, limma, org.Hs.eg.db, AnnotationDbi, psych,
  ggplot2, patchwork, scales
- **Python**: pandas, numpy

## Inputs
- Counts matrix (genes × samples) per cohort
- Sample metadata (cohort-specific clinical + RNA-seq metadata)
- MASTER microprotein annotation table with `RNA_Ensembl_Parent` and a
  smORF coordinate / symbol column

## Outputs
- DESeq2 objects, VST and batch-corrected VST matrices
- AD vs. Control DE results (CSV)
- Per-pair main-ORF / smORF LRT statistics (CSV)
- Per-pair triptych co-expression figures routed by Δr (PDF)
