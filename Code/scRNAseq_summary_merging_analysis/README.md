# Single-cell RNA-seq Analysis

Integrates published per-cell-type differential-expression statistics from
the Mathys et al. (2024) ROSMAP snRNA-seq atlas
(<https://www.nature.com/articles/s41586-024-07606-7>) with our
microprotein annotation table to produce the cell-type enrichment view
used in the dashboard and the manuscript.

## Overview
- **Input**:
  - `Code/data/microprotein_master.csv` (master annotation + evidence).
  - Per-cell-type DE tables published as supplementary data with Mathys
    et al. 2024 (downloaded separately; not redistributed).
- **Output** (written to `../../Results/scRNA_Enrichment/`):
  - `scRNA_Enrichment_summary.csv` - tidy per-microprotein x cell-type
    enrichment statistics.
  - `heatmap_PSM.pdf`, `heatmap_log2FC.pdf`,
    `cell_type_smorf_type_heatmap.pdf` - cell-type heatmaps.
  - `volcano_all_celltypes.pdf` - per-cell-type volcano panel.
  - `UpSet_*.pdf`, `UpSet_input_matrix_microproteins.csv`,
    `UpSet_sharedness_summary.csv` - cross-cell-type overlap analysis.
- **Main script**: `scRNAseq_summary.R`.

## Usage

```bash
Rscript scRNAseq_summary.R
```

Also runs as Phase 4 of `bash run_all_analyses.sh --mode=run`.

## Dependencies
- R: `dplyr`, `tibble`, `ggplot2`, `patchwork`, `ComplexHeatmap`,
  `UpSetR`, `org.Hs.eg.db`, `AnnotationDbi`.

## Notes
This module does **not** re-process raw single-cell counts; the upstream
clustering / DE was performed by the original Mathys et al. authors and we
consume their published statistics.
