# `Code/` — analysis modules

Each subdirectory mirrors one experimental platform and has two halves:

- a **processing** half, which needs the raw instrument or sequencing data. Those
  inputs are under controlled access and are not shipped here, so these scripts
  are documented for reference but are not runnable from a clone.
- a **summary** half, which runs entirely on `data/microprotein_master.zip` and
  writes the CSVs under `Results/` that the dashboard and the manuscript figures
  consume. These *are* runnable from a clone.

`bash run_all_analyses.sh --mode=run` (from the repo root) runs every summary
script in order. `--mode=docs` prints the step list, including the processing
steps, without executing anything.

## Modules

| Directory | What it covers | Summary script |
|-----------|----------------|----------------|
| [`Microprotein_annotation_summary/`](Microprotein_annotation_summary/) | smORF discovery and classification (`Annotator/`); exports the BED / GTF / FASTA bundles in `GTF_and_BED_files/` | `Brain_Microproteins_Discovery_summary.py`, `ShortStop_Microproteins_summary.py`, `Create_BED_GTF_FASTA_files.py`, `Create_Alt_Proteoform_Files.py` |
| [`Codon_context/`](Codon_context/) | Start-codon context: Kozak strength and non-ATG initiation | `initiation_pipeline.py`, `kozak_pipeline.py` |
| [`RP3_analysis/`](RP3_analysis/) | Ribo-seq translation evidence (RiboCode) | `RP3_Results_summary.py` |
| [`Peptide_TMT_analysis/`](Peptide_TMT_analysis/) | FragPipe TMT-MS processing, TAMPOR batch correction, and PROSIT spectral-angle validation | `Proteomics_Results_summary.py` |
| [`Shortread_RNA_analysis/`](Shortread_RNA_analysis/) | ROSMAP / MSBB short-read DESeq2 differential expression (R) | `Short-Read_Transcriptomics_Results_summary.py` |
| [`Longread_RNA_analysis/`](Longread_RNA_analysis/) | Nanopore ESPRESSO isoform differential expression (R) | `Long-Read_Transcriptomics_Results_summary.py` |
| [`scRNAseq_summary_merging_analysis/`](scRNAseq_summary_merging_analysis/) | Mathys et al. 2024 cell-type enrichment (R) | `scRNAseq_summary.R` |
| [`Miscellaneous/`](Miscellaneous/) | Actin fluorescence-microscopy quantitation — standalone, not in the pipeline | — |

Each module's own README documents its expected inputs and outputs.

## Shared across modules

| File | Role |
|------|------|
| [`gold_standard_filtering_criteria.py`](gold_standard_filtering_criteria.py) | **The canonical filter.** Turns `data/microprotein_master.csv` into the deduplicated, evidence-gated microprotein set used everywhere — figures, summary CSVs, and the dashboard. Every summary script adds `Code/` to `sys.path` and imports `load_and_filter_master` from it; the R modules use the equivalent `R_TEMPLATE` embedded in the same file. Change filtering criteria here and nowhere else. |
| [`data/`](data/) | Shared input tables. See [`data/README.md`](data/README.md). |
| `verify_gtf_bed_fasta_files.py` | Checks the generated `GTF_and_BED_files/` bundles for internal consistency. |

The dashboard (`Results/microproteins_dashboard.py`) carries an inlined copy of
the gold-standard filter rather than importing it, so it can be deployed
standalone. If you change the filter, change the dashboard's copy to match.

One divergence between the two is deliberate: the dashboard counts only DDA
evidence toward `has_MS`, while the manuscript figures use the DIA-inclusive
variant. This accounts for a known 14-microprotein difference between the
dashboard's totals and the figures'.
