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
| `Code/data/microprotein_master.zip`                 | Compressed master annotation + evidence table that drives every summary script (full uncompressed `microprotein_master.csv` is on Zenodo) |
| `Code/data/brain_espresso_medianCPM05.csv`          | Long-read ESPRESSO transcript abundances (CPM ≥ 0.5)               |
| `Code/data/cleaned_tryptic_peptides_under_151aa.csv`| Microprotein tryptic peptide list (input to PROSIT)                |
| `Code/data/cleaned_tryptic_peptides_detailed_under_151aa.csv` | Per-peptide metadata (PROSIT input)                       |
| `Code/data/ac_list_mapping.csv`                     | Actin-related microprotein ID mapping                              |
| `Code/data/nanopore_metadata.csv`                   | Long-read sample metadata                                          |
| `Code/data/rbp_splice_per_junction.csv`             | RBP splice-junction features                                       |
| `Code/data/alt_proteoform_table.csv`                | Alternative-initiation proteoform records                          |
| `Code/data/codon_context_sequences.csv.gz`          | Start-codon sequence context derived from the ESPRESSO assembly (input to `Code/Codon_context/`) |
| `Code/data/uniprotkb_proteome_UP000005640_2026_07_13.tsv` | UniProt human proteome annotation scores, read by the dashboard |
| `Code/gold_standard_filtering_criteria.py`          | The canonical master-table filter. Every summary script and the dashboard apply it — nothing downstream reproduces without it |
| `GTF_and_BED_files/`                                 | Generated hg38 genomic coordinates + FASTAs for unreviewed microproteins (the ~890 MB combined Ensembl + microprotein GTF is on Zenodo) |
| `Results/Annotations/`, `Results/Codon_context/`, `Results/Proteomics/`, `Results/RP3/`, `Results/Transcriptomics/`, `Results/scRNA_Enrichment/`, `Results/ShortStop/` | Per-analysis summary CSVs (and small PDFs) consumed by the dashboard |

> **Not shipped in the git repo:** the full figure libraries (mirror plots,
> expression-profile triptychs, smORF cartoons), the uncompressed master table,
> and several large intermediate tables. The dashboard streams the figures on
> demand from the
> [Hugging Face dataset](https://huggingface.co/datasets/brmiller/brain-microprotein-atlas),
> and everything large is archived on Zenodo (next section).

---

## Available via Zenodo

Files that exceed GitHub's size limits are archived on Zenodo
([10.5281/zenodo.20045161](https://doi.org/10.5281/zenodo.20045161)) as six
tarballs:

| Archive                            | Contents                                                                                       |
|------------------------------------|------------------------------------------------------------------------------------------------|
| `large_data_tables.tar.gz`         | Full uncompressed `microprotein_master.csv`; `combined_gpath_results.csv` (pathway analysis); `tss_motif_prevalence.csv` (TSS motif frequencies) |
| `large_reference_files.tar.gz`     | Complete `GTF_and_BED_files/` including the combined Ensembl + microprotein GTF and all PROSIT-tier / Ribo-ShortStop / `Alt_Proteoforms/` bundles, plus the raw RiboCode outputs (`Results/RP3/ribocode_results*.{gtf,bed,txt}`) |
| `figures_mirror_plots.tar.gz`      | PROSIT MS2 mirror plots (~11k PDFs/PNGs)                                                        |
| `figures_expression_profiles.tar.gz` | Main-ORF / smORF coupling triptychs (~8.6k figures)                                          |
| `figures_smorf_cartoons.tar.gz`    | Per-locus smORF cartoons (~8.7k figures)                                                        |
| `espresso_sequence_files.tar.gz`   | ESPRESSO long-read transcript assembly: `.gtf` (transcript models), `.nuc` (transcript sequences), `.split_nuc` (per-ORF spliced CDS) and `_ORFs.gtf` (per-ORF exon blocks). Upstream inputs to `Code/Codon_context/build_sequence_context.py`; the analysis pipelines read only its 352 KB derived output, which ships in the repo |

---

## Raw data sources (not redistributed)

### 1. Long-read RNA-seq (ESPRESSO output)
- **Source**: Heberle et al., ROSMAP nanopore cohort.
- **Repository**: AD Knowledge Portal (Synapse).
- **Inputs needed to re-run**: `.esp` ESPRESSO abundance files; raw FASTQ via
  Synapse (data-use agreement required).
- **Processing scripts**: `Code/Longread_RNA_analysis/ESPRESSO_data_processing_scripts/`.
- **The assembled output is redistributed**: the ESPRESSO transcript models and
  their nucleotide sequences are in `espresso_sequence_files.tar.gz` on Zenodo,
  so the assembly itself can be inspected without re-running ESPRESSO or
  obtaining the raw FASTQ.

### 2. Short-read RNA-seq (ROSMAP DLPFC, MSBB)
- **Source**: ROSMAP and MSBB consortia.
- **Repository**: AD Knowledge Portal (Synapse) — restricted access.
- **Inputs needed to re-run**: FeatureCounts matrices generated against a
  custom GTF that combines GENCODE v43 with our unreviewed microproteins
  (DDA, DIA, RiboCode, ShortStop sources); ROSMAP / MSBB clinical and assay
  metadata.
- **Restricted count matrix**: the bulk RNA-seq count matrix
  (`Code/data/counts.csv`) is **not** redistributed — it is derived from
  restricted ROSMAP / MSBB individual-level data and requires AMP-AD / Synapse
  data-access approval. The DESeq scripts expect it at `Code/data/counts.csv`.
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
- **Processing script**: `Code/Miscellaneous/actin_quant_pipeline.py`.

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

