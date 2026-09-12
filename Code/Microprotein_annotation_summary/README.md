# smORF Annotation & Discovery Summaries

This module provides:

1. The **Annotator pipeline** (`Annotator/`) used to classify candidate
   smORF GTFs against Ensembl / UniProt references.
2. The **summary scripts** that turn the master annotation table into the
   manuscript / dashboard discovery tables.
3. The **genomic-coordinate file generator** that produces the FASTA, BED,
   and GTF distributions in `GTF_and_BED_files/` for the unreviewed brain
   microproteins.

## Files

| File | Purpose |
|------|---------|
| `run_Annotator.sh` | Bash wrapper that runs the full Annotator workflow on a smORF GTF + Ensembl reference. |
| `Annotator/Annotator.py` | Main Python annotation engine (smORF type classification). |
| `Annotator/check_genes.py` | Utility for extracting and validating gene CDS coordinates from GTF. |
| `Annotator/src/` | Supporting package: `annotation/smorf_annotator.py`, `annotation/bedtools_smorf_intersect.py`, `pipeline/pipeline.py`, `pipeline/pipeline_structure.py`. |
| `Brain_Microproteins_Discovery_summary.py` | Filters the master table to Salk-discovered microproteins and writes the discovery summary used by the dashboard. |
| `ShortStop_Microproteins_summary.py` | Builds the ShortStop ML-classification summary. |
| `Create_BED_GTF_FASTA_files.py` | Generates the BED / GTF / FASTA / coordinate-mapping files in `GTF_and_BED_files/` for unreviewed brain microproteins, including the PROSIT confidence tiers and the Ribo-ShortStop bundle. |
| `Create_Alt_Proteoform_Files.py` | Generates the `GTF_and_BED_files/Alt_Proteoforms/` bundle for alternative-initiation proteoforms. |
| `alt_proteoform_records.py` | Record-building helpers imported by `Create_Alt_Proteoform_Files.py`. |
| `compute_nterm_peptide_substitutions.py` | BLASTs each N-terminal (aa 1-2) tryptic peptide against its matched UniProt isoform and records amino-acid substitutions, for the dashboard's N-terminus filter (`Code/data/blast_nterm_peptide_substitutions.csv`). |

## Usage

```bash
# Full annotator pipeline (requires raw smORF GTF + Ensembl reference)
bash run_Annotator.sh

# Summary scripts (run from this directory; reads Code/data/microprotein_master.csv)
python Brain_Microproteins_Discovery_summary.py
python ShortStop_Microproteins_summary.py
python Create_BED_GTF_FASTA_files.py
python Create_Alt_Proteoform_Files.py

# Dashboard N-terminus filter support (needs network access to UniProt only
# for accessions not already cached in Code/data/blast_isoform_reference_sequences.fasta)
python compute_nterm_peptide_substitutions.py
```

All four summary scripts are also invoked by
`bash run_all_analyses.sh --mode=run` from the repo root. The Annotator
pipeline and `compute_nterm_peptide_substitutions.py` are not — the former
needs raw data not shipped here, the latter needs network access on first run.

The ribosome-profiling summary lives in `Code/RP3_analysis/`, not this module.

## Dependencies
- Python: `pandas`, `numpy`, `pybedtools`.
- `bedtools` (binary on `PATH`) for genomic intersections.
- Shared filter: `Code/gold_standard_filtering_criteria.py`.

## Inputs
- Raw smORF GTF files (Annotator pipeline only - not in repo).
- Ensembl reference GTF (Annotator pipeline only - not in repo).
- `Code/data/microprotein_master.csv` (shipped).

## Outputs
- `../../Results/Annotations/Brain_Microproteins_Discovery_summary.csv`
- `../../Results/Annotations/ShortStop_Microproteins_summary.csv`
- `../data/blast_isoform_reference_sequences.fasta`, `../data/blast_nterm_peptide_substitutions.csv`
  (from `compute_nterm_peptide_substitutions.py`)
- `../../GTF_and_BED_files/`:
  - `Unreviewed_Brain_Microproteins.fasta`
  - `Unreviewed_Brain_Microproteins_Absent_from_UniProt.gtf`
  - `Unreviewed_Brain_Microproteins_CDS_Absent_from_UniProt.bed`
  - `Unreviewed_Brain_Microproteins_All_Non_SwissProt.bed` (the same CDS records
    plus the TrEMBL entries as single-interval CDS envelopes; `itemRgb`
    colour-coded by record type)
  - `Unreviewed_Brain_Microproteins_IDs.txt`
  - `Unreviewed_Brain_Microproteins_genomic_coordinates.txt`
  - `Unreviewed_Brain_Microproteins_mapping_coordinates_to_sequences.tsv`
  - `Ensembl_and_Unreviewed_Brain_Microproteins.gtf` (only when the combined
    source GTF is present locally; too large for git)
  - `PROSIT_confidence_tiers/` — 16 files, four per Strong / Moderate / Weak /
    Insufficient tier (FASTA, GTF, CDS BED, ID list)
  - `Ribo_ShortStop/` — 4 files for the 1,592 microproteins with RiboCode +
    ShortStop translation support but no DDA mass-spec detection
  - `Alt_Proteoforms/` — 5 files describing alternative-initiation proteoforms,
    from `Create_Alt_Proteoform_Files.py`

`Results/Annotations/smORF_type_definitions.csv` sits alongside the two summary
CSVs above and is read by the dashboard, but it is a hand-curated lookup table
shipped in the repo — this module does not generate it.

See [`../../GTF_and_BED_files/README.md`](../../GTF_and_BED_files/README.md) for
per-file record counts and column definitions.
