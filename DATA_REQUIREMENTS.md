# Data Requirements and File Formats

This document specifies the input data files expected by the analysis
pipelines in this repository, and the formats they must follow.

For *where* to obtain the raw inputs (Synapse, dbGaP, etc.) see
[DATA_AVAILABILITY.md](DATA_AVAILABILITY.md). For *what* the summary
pipeline outputs, see [run_all_analyses.sh](run_all_analyses.sh) and the
per-module README files.

---

## 1. Files shipped in `Code/data/`

These files come with a clone and are sufficient to run the summary pipeline
(`bash run_all_analyses.sh --mode=run`) and the dashboard.

| File | Description |
|------|-------------|
| `microprotein_master.zip` | Master annotation + evidence table, compressed. `run_all_analyses.sh` unzips it to `microprotein_master.csv` on first run. Loaded by every summary script via `gold_standard_filtering_criteria.py`. Contains classification, sequences, MS / Ribo-Seq / RNA-seq evidence columns. |
| `brain_espresso_medianCPM05.csv` | Long-read ESPRESSO transcript abundances filtered at median CPM ≥ 0.5. |
| `cleaned_tryptic_peptides_under_151aa.csv` | List of unique tryptic peptides per microprotein (≤ 151 aa). |
| `cleaned_tryptic_peptides_detailed_under_151aa.csv` | Detailed per-peptide metadata; primary input to the PROSIT pipeline. |
| `ac_list_mapping.csv` | Mapping between accession identifiers and Ensembl gene IDs (used by the actin/microprotein analysis). |
| `alt_proteoform_table.csv` | Alternative-initiation proteoform records, used to build `GTF_and_BED_files/Alt_Proteoforms/`. |
| `nanopore_metadata.csv` | Long-read sample metadata (used by the ESPRESSO DESeq2 step). |
| `rbp_splice_per_junction.csv` | Per-junction RBP feature table. |
| `uniprotkb_proteome_UP000005640_2026_07_13.tsv` | UniProt human proteome export; supplies annotation scores to the dashboard. |
| `codon_context_sequences.csv.gz` | Per-smORF nucleotide sequence, exon blocks and 12 nt of flanking transcript either side. Derived from the ESPRESSO long-read sequence files by `Code/Codon_context/build_sequence_context.py`; the only input the two start-codon pipelines need beyond the master table. |

### Referenced by the code but *not* in a clone

These four live under `Code/data/` in the authors' working tree but are excluded
from git — either because they exceed GitHub's size limits or because they derive
from controlled-access individual-level data. Download them from
[Zenodo](https://doi.org/10.5281/zenodo.20045161) (`large_data_tables.tar.gz`),
except `counts.csv`, which requires AMP-AD / Synapse approval.

| File | Why it is absent | Needed by |
|------|------------------|-----------|
| `microprotein_master.csv` | 41 MB uncompressed; ships as `.zip` instead | everything (unzipped automatically) |
| `counts.csv` | Restricted ROSMAP / MSBB individual-level data | `Code/Shortread_RNA_analysis/` DESeq2 scripts |
| `combined_gpath_results.csv` | Size | `Code/scRNAseq_summary_merging_analysis/scRNAseq_summary.R` |
| `tss_motif_prevalence.csv` | Size | TSS motif figures |

The canonical filter applied to `microprotein_master.csv` before any per-module
summary is defined once, in `Code/gold_standard_filtering_criteria.py`. Each
summary script adds `Code/` to `sys.path` and imports `load_and_filter_master`
from it; the dashboard carries an inlined copy of the same logic. Change the
criteria in that one file.

---

## 2. Inputs required to re-run the upstream processing

These are **not** included in the repository (controlled access). Each
module's README documents its expected inputs in more detail.

### Long-read RNA-seq (`Code/Longread_RNA_analysis/`)
- ESPRESSO `*_abundance.esp` files
  (e.g. `ESPRESSO_syn52047893_hg38NCBI_N2_R0_abundance.esp`).
- `nanopore_metadata.csv` (shipped).

### Start-codon context (`Code/Codon_context/`)
- The three ESPRESSO sequence files sharing one path prefix:
  `<prefix>.split_nuc` (spliced ORF CDS), `<prefix>_ORFs.gtf` (exon blocks) and
  `<prefix>.nuc` (assembled transcripts), ~1.5 GB total.
- Needed only by `build_sequence_context.py`, to regenerate
  `codon_context_sequences.csv.gz`. Both analysis pipelines read that derived
  table and never touch these.
- Available from Zenodo as `espresso_sequence_files.tar.gz` (~158 MB); extract
  anywhere and pass the shared path prefix to `--espresso-prefix`.

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

### Microscopy (`Code/Miscellaneous/actin_quant_pipeline.py`)
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
Results/Codon_context/kozak/kozak_context.csv
Results/Codon_context/initiation/initiation_summary.csv
Results/Proteomics/Proteomics_Results_summary.csv
Results/RP3/RP3_Results_summary.csv
Results/ShortStop/ShortStop_Microproteins_summary.csv
Results/Transcriptomics/Short-Read_Transcriptomics_Results_summary.csv
Results/Transcriptomics/Long-Read_Transcriptomics_Results_summary.csv
Results/scRNA_Enrichment/scRNA_Enrichment_summary.csv
Results/scRNA_Enrichment/scRNA_Enrichment_all_celltypes.csv.gz
GTF_and_BED_files/Unreviewed_Brain_Microproteins_CDS_Absent_from_UniProt.bed
supplementary/Supplemental_Tables.xlsx
```

`Results/Annotations/smORF_type_definitions.csv` is also read by the dashboard,
but it is a hand-curated lookup shipped in the repo — no script generates it.

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
