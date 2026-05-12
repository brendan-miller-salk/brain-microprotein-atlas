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
| `Brain_Microproteins_Discovery_summary.py` | Filters the master table to Salk-discovered microproteins and writes the discovery summary used by the dashboard. |
| `ShortStop_Microproteins_summary.py` | Builds the ShortStop ML-classification summary. |
| `RP3_Results_summary.py` | Annotation-side complement to the RP3 module's ribosome-profiling summary. |
| `Create_BED_GTF_FASTA_files.py` | Generates the BED / GTF / FASTA / coordinate-mapping files in `GTF_and_BED_files/` for unreviewed brain microproteins. |

## Usage

```bash
# Full annotator pipeline (requires raw smORF GTF + Ensembl reference)
bash run_Annotator.sh

# Summary scripts (run from this directory; reads Code/data/microprotein_master.csv)
python Brain_Microproteins_Discovery_summary.py
python ShortStop_Microproteins_summary.py
python RP3_Results_summary.py
python Create_BED_GTF_FASTA_files.py
```

All four summary scripts are also invoked by
`bash run_all_analyses.sh --mode=run` from the repo root.

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
- `../../Results/Annotations/smORF_type_definitions.csv`
- `../../GTF_and_BED_files/`:
  - `Unreviewed_Brain_Microproteins.fasta`
  - `Unreviewed_Brain_Microproteins_absent_from_UniProt.gtf`
  - `Unreviewed_Brain_Microproteins_CDS_absent_from_UniProt.bed`
  - `Unreviewed_Brain_Microproteins_IDs.txt`
  - `Unreviewed_Brain_Microproteins_genomic_coordinates.txt`
  - `Unreviewed_Brain_Microproteins_mapping_coordinates_to_sequences.tsv`
  - `Ensembl_and_Unreviewed_Brain_Microproteins.gtf`
