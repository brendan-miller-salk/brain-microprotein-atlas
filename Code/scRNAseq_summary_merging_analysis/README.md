# Single-cell RNA-seq Analysis

Integrates published per-cell-type differential-expression statistics from
the Mathys et al. (2024) ROSMAP snRNA-seq atlas
(<https://www.nature.com/articles/s41586-024-07606-7>) with our
microprotein annotation table to produce the cell-type enrichment view
used in the dashboard and the manuscript.

## Overview
- **Input**:
  - `Code/data/microprotein_master.csv` (master annotation + evidence).
  - `Code/data/combined_gpath_results.csv` — per-cell-type DE statistics
    derived from the supplementary data published with Mathys et al. 2024.
    Too large for git; download it from
    [Zenodo](https://doi.org/10.5281/zenodo.20045161)
    (`large_data_tables.tar.gz`) into `Code/data/` before running.
  - `Code/data/ac_list_mapping.csv` — cell-type / accession lookup (shipped).
- **Output** (written to `../../Results/scRNA_Enrichment/`):
  - `scRNA_Enrichment_summary.csv` - tidy per-microprotein x cell-type
    enrichment statistics, **filtered to `p_adj.glb < 0.05`** (this is the
    table the heatmaps, volcanoes and the dashboard's scRNA view are built
    from).
  - `scRNA_Enrichment_all_celltypes.csv.gz` - the same statistics *before*
    that filter: every tested (microprotein, cell type) pair, with a
    `significant` flag stamped pre-rounding. The dashboard's per-microprotein
    cell-type panel reads this so it can draw non-significant cell types
    instead of leaving "tested but not significant" indistinguishable from
    "not tested". ~149k rows, sorted by sequence so it gzips to <3 MB.
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
