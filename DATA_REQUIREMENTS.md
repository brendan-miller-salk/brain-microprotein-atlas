# Data Requirements and File Formats

This document specifies the input data files expected by the analysis
pipelines in this repository, and the formats they must follow.

For *where* to obtain the raw inputs (Synapse, dbGaP, etc.) see
[DATA_AVAILABILITY.md](DATA_AVAILABILITY.md). For *what* the summary
pipeline outputs, see [run_all_analyses.sh](run_all_analyses.sh) and the
per-module README files.

---

## 1. Files shipped in `Code/data/`

These files are included in the repository and are sufficient to run the
summary pipeline (`bash run_all_analyses.sh --mode=run`) and the dashboard.

| File | Description |
|------|-------------|
| `microprotein_master.csv` (`.zip`) | Master annotation + evidence table. Loaded by every summary script via `gold_standard_filtering_criteria.py`. Contains classification, sequences, MS / Ribo-Seq / RNA-seq evidence columns. |
| `microprotein_master_first_submission.csv` (`.zip`) | Archived snapshot from the initial submission (kept for provenance). |
| `brain_espresso_medianCPM05.csv` | Long-read ESPRESSO transcript abundances filtered at median CPM ≥ 0.5. |
| `counts.csv` | Example bulk RNA-seq count matrix (genes × samples). |
| `cleaned_tryptic_peptides_under_151aa.csv` | List of unique tryptic peptides per microprotein (≤ 151 aa). |
| `cleaned_tryptic_peptides_detailed_under_151aa.csv` | Detailed per-peptide metadata; primary input to the PROSIT pipeline. |
| `combined_gpath_results.csv` | Combined pathway analysis results. |
| `ac_list_mapping.csv` | Mapping between accession identifiers and Ensembl gene IDs (used by the actin/microprotein analysis). |
| `nanopore_metadata.csv` | Long-read sample metadata (used by the ESPRESSO DESeq2 step). |
| `rbp_splice_per_junction.csv` | Per-junction RBP feature table. |
| `tss_motif_prevalence.csv` | TSS motif frequency table. |

The `gold_standard_filtering_criteria.py` helper at the repo root and at
`Code/gold_standard_filtering_criteria.py` defines the canonical filter
applied to `microprotein_master.csv` before any per-module summary.

---

## 2. Inputs required to re-run the upstream processing

These are **not** included in the repository (controlled access). Each
module's README documents its expected inputs in more detail.

### Long-read RNA-seq (`Code/Longread_RNA_analysis/`)
- ESPRESSO `*_abundance.esp` files
  (e.g. `ESPRESSO_syn52047893_hg38NCBI_N2_R0_abundance.esp`).
- `nanopore_metadata.csv` (shipped).

### Short-read RNA-seq (`Code/Shortread_RNA_analysis/`)
- Per-cohort counts matrices (genes × samples) generated with FeatureCounts
  against the combined GENCODE v43 + unreviewed-microprotein GTF.
- ROSMAP metadata: `synapse_download_rnaseq.csv`,
  `ROSMAP_assay_rnaSeq_metadata.csv`, `ROSMAP_clinical.csv`.
- MSBB metadata: cohort-specific clinical + RNA-seq metadata.
- The `RNA_differential_expression.R` `CONFIG` block selects ROSMAP vs MSBB.

### TMT Proteomics (`Code/Peptide_TMT_analysis/`)
- FragPipe output laid out as
  `{batch}/shortstop_proteogenomics_appended_results_cpm05/DDA/{peptide.tsv,psm.tsv}`
  for round 1 batches `b1`–`b50` and round 2 batches `round2/b1`–`round2/b14`.
- Source `*.mzML` spectra in each batch directory (PROSIT spectral-angle phase).
- TMT channel metadata.

### Ribosome Profiling (`Code/RP3_analysis/`)
- RiboCode output files (`.bed`, `.gtf`, `.txt`) and RPKM mapping group files.

### Annotation pipeline (`Code/Microprotein_annotation_summary/`)
- Input smORF GTF files.
- Ensembl reference annotations (GTF).
- ShortStop pipeline outputs.

### Single-cell RNA-seq (`Code/scRNAseq_summary_merging_analysis/`)
- Per-cell-type DE tables from Mathys et al. (2024) supplementary data.

### Microscopy (`Code/Miscellanous/actin_quant_pipeline.py`)
- Multi-channel ND2 / TIFF fluorescence stacks (DAPI + phalloidin + signal).

---

## 3. File format examples

### `microprotein_master.csv` (excerpt of key columns)
```csv
Database,protein_class_length,discovery_origin,Global.PG.Q.Value,gene_symbol,sequence
Salk,Microprotein,proteogenomics,0.001,GENE1,MATKL...
Salk,Microprotein,ribocode_shortstop,,GENE2,MVKLS...
```

### Counts / expression matrix
```csv
gene_id,sample1,sample2,sample3,...
ENSG00000001,12.5,8.2,15.1,...
ENSG00000002,0.0,2.1,1.8,...
```

### PROSIT input (`cleaned_tryptic_peptides_detailed_under_151aa.csv`)
- Required columns: `protein_id`, `peptide_sequence` (Python list-style:
  `['AAAPQAPAAR', 'RAAAPQAPAAR']`).

---

## 4. Summary outputs (consumed by the dashboard)

After running `bash run_all_analyses.sh --mode=run`, the dashboard expects
the following files (written by the pipeline):

```
Results/Annotations/Brain_Microproteins_Discovery_summary.csv
Results/Annotations/ShortStop_Microproteins_summary.csv
Results/Annotations/smORF_type_definitions.csv
Results/Proteomics/Proteomics_Results_summary.csv
Results/RP3/RP3_Results_summary.csv
Results/Transcriptomics/Short-Read_Transcriptomics_Results_summary.csv
Results/Transcriptomics/Long-Read_Transcriptomics_Results_summary.csv
Results/scRNA_Enrichment/scRNA_Enrichment_summary.csv
GTF_and_BED_files/Unreviewed_Brain_Microproteins_CDS_absent_from_UniProt.bed
supplementary/Supplemental_Tables.xlsx
```

Optional figure directories used by the dashboard:
- `Results/expression_profiles/{coupled,non_coupled}/` — main-ORF / smORF
  triptychs (PNG/PDF, named by `GENE_chrX_start-end`).
- `Results/mirror_plots/{Strong,Moderate,Weak,Insufficient}/` — PROSIT
  mirror plots (PNG/PDF, named `PEPTIDE_zCHARGE_SCAN`).
- `Results/smorf_cartoon_figures/` — per-locus smORF cartoons keyed by
  genomic coordinate.

---

## Getting Help

For clarification on file formats or help adapting the pipeline to your data:
**Brendan A. Miller** — `brmiller@salk.edu`.
