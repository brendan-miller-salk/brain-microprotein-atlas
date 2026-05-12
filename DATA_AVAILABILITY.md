# Data Availability

This repository contains the analysis code, the processed master annotation
table (`Code/data/microprotein_master.csv`), the per-analysis summary CSVs
under `Results/`, and the interactive dashboard
(`Results/microproteins_dashboard.py`). Raw input files are **not**
redistributed here because most are protected (consent / data-use
agreements) or impractically large.

This document lists where each raw data class can be obtained, what is
included in the repository as a substitute, and what is needed to fully
re-run the upstream processing.

---

## What ships with the repository

| File / Directory                                    | Purpose                                                            |
|-----------------------------------------------------|--------------------------------------------------------------------|
| `Code/data/microprotein_master.csv` (+ `.zip`)      | Master annotation + evidence table that drives every summary script |
| `Code/data/microprotein_master_first_submission.csv`| Archived version from the initial manuscript submission            |
| `Code/data/brain_espresso_medianCPM05.csv`          | Long-read ESPRESSO transcript abundances (CPM ≥ 0.5)               |
| `Code/data/counts.csv`                              | Example bulk RNA-seq count matrix                                  |
| `Code/data/cleaned_tryptic_peptides_under_151aa.csv`| Microprotein tryptic peptide list (input to PROSIT)                |
| `Code/data/cleaned_tryptic_peptides_detailed_under_151aa.csv` | Per-peptide metadata (PROSIT input)                       |
| `Code/data/combined_gpath_results.csv`              | Pathway analysis results                                           |
| `Code/data/ac_list_mapping.csv`                     | Actin-related microprotein ID mapping                              |
| `Code/data/nanopore_metadata.csv`                   | Long-read sample metadata                                          |
| `Code/data/rbp_splice_per_junction.csv`             | RBP splice-junction features                                       |
| `Code/data/tss_motif_prevalence.csv`                | TSS motif frequencies                                              |
| `GTF_and_BED_files/`                                | Generated genomic coordinates + FASTAs for unreviewed microproteins |
| `Results/Annotations/`, `Results/Proteomics/`, `Results/RP3/`, `Results/Transcriptomics/`, `Results/scRNA_Enrichment/`, `Results/ShortStop/` | Summary CSVs consumed by the dashboard |
| `Results/expression_profiles/`                      | Per-pair main-ORF / smORF triptych figures (PDF/PNG)               |
| `Results/mirror_plots/`                             | PROSIT 3-panel mirror plots stratified by `Confidence` tier        |
| `Results/smorf_cartoon_figures/`                    | Per-locus smORF cartoons (PDF/PNG)                                 |
| `supplementary/`                                    | Final supplemental tables (xlsx + csv)                             |

---

## Raw data sources (not redistributed)

### 1. Long-read RNA-seq (ESPRESSO output)
- **Source**: Heberle et al., ROSMAP nanopore cohort.
- **Repository**: AD Knowledge Portal (Synapse).
- **Inputs needed to re-run**: `.esp` ESPRESSO abundance files; raw FASTQ via
  Synapse (data-use agreement required).
- **Processing scripts**: `Code/Longread_RNA_analysis/ESPRESSO_data_processing_scripts/`.

### 2. Short-read RNA-seq (ROSMAP DLPFC, MSBB)
- **Source**: ROSMAP and MSBB consortia.
- **Repository**: AD Knowledge Portal (Synapse) — restricted access.
- **Inputs needed to re-run**: FeatureCounts matrices generated against a
  custom GTF that combines GENCODE v43 with our unreviewed microproteins
  (DDA, DIA, RiboCode, ShortStop sources); ROSMAP / MSBB clinical and assay
  metadata.
- **Processing scripts**: `Code/Shortread_RNA_analysis/Shortread_deseq_processing_scripts/`
  (DESeq2, main-ORF / smORF LRT, co-expression).

### 3. TMT Proteomics (FragPipe output, two rounds)
- **Source**: ROSMAP TMT-MS proteomics.
- **Repository**: AD Knowledge Portal (Synapse) — restricted access.
- **Inputs needed to re-run**:
  - FragPipe `peptide.tsv` / `psm.tsv` for round 1 (`b1`–`b50`) and round 2
    (`round2/b1`–`round2/b14`).
  - Source `*.mzML` spectra (for PROSIT spectral-angle validation).
  - TMT channel metadata.
- **Processing scripts**: `Code/Peptide_TMT_analysis/fragpipe_results_processing_scripts/`
  (PSM aggregation, batch correction, TAMPOR rounds 1/2/combined, ANOVA).
- **PROSIT validation**: `Code/Peptide_TMT_analysis/prosit/prosit_pipeline.py`
  uses the `Prosit_2020_intensity_TMT` model via the
  [Koina](https://koina.wilhelmlab.org) API; requires network access to
  `koina.wilhelmlab.org:443`.

### 4. Ribosome Profiling (RP3 / RiboCode)
- **Source**: Duffy et al. brain Ribo-seq.
- **Repository**: dbGaP — controlled access.
- **Inputs needed to re-run**: RiboCode outputs (`.bed`, `.gtf`, `.txt`,
  RPKM mapping group files). Raw FASTQ requires dbGaP approval.
- **Processing scripts**: available on request (the RP3 module here only
  consumes pre-computed RiboCode tables).

### 5. Single-cell RNA-seq
- **Source**: Mathys et al. (2024) ROSMAP snRNA-seq atlas
  (https://www.nature.com/articles/s41586-024-07606-7) and other published
  datasets cited in the manuscript.
- **Inputs needed to re-run**: per-cell-type DE statistics published as
  supplementary data with the original studies.
- **Processing scripts**: `Code/scRNAseq_summary_merging_analysis/scRNAseq_summary.R`.

### 6. Microscopy (actin imaging)
- **Source**: in-house F-actin fluorescence microscopy (multi-channel
  ND2 / TIFF). Available on request.
- **Processing script**: `Code/Miscellanous/actin_quant_pipeline.py`.

---

## Reproducing summary results without raw data

Everything under `Results/` and `supplementary/` can be regenerated from the
shipped master table:

```bash
bash setup_environment.sh
bash run_all_analyses.sh --mode=run
```

The dashboard
([Results/microproteins_dashboard.py](Results/microproteins_dashboard.py))
reads only those summary CSVs plus the figure directories, so it is fully
runnable from a fresh clone.

---

## Data Access & Collaboration

For access to processed intermediates that are not in the repository (TMT
intensity matrices, batch-corrected VST matrices, FragPipe protein-ID
tables, etc.) or for collaboration:

- **Email**: `brmiller@salk.edu`
- **Institution**: Salk Institute
- **Lab**: PBL-A

We are happy to share processed data under appropriate data-sharing
agreements and AD Knowledge Portal data-access approvals.

