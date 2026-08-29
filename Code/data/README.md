# `Code/data/` — shared inputs

Every analysis module reads from this directory. Paths are resolved relative to
it, so the files must keep these exact names.

## The master table

`microprotein_master.zip` is the single annotation + evidence table that drives
everything downstream. `run_all_analyses.sh` unzips it to
`microprotein_master.csv` (39 MB) on first run; only the compressed form is in
git.

Never read it raw. Load it through the canonical filter, which deduplicates and
applies the evidence gates used throughout the manuscript:

```python
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))   # Code/
from gold_standard_filtering_criteria import load_and_filter_master

mp = load_and_filter_master('Code/data/microprotein_master.csv')
```

`Code/gold_standard_filtering_criteria.py` also embeds an equivalent R version
as the `R_TEMPLATE` string, for the R-based modules.

## Shipped in git

| File | Size | Used by |
|------|-----:|---------|
| `microprotein_master.zip` | 15 M | everything |
| `brain_espresso_medianCPM05.csv` | 3.9 M | `Longread_RNA_analysis/` |
| `nanopore_metadata.csv` | 397 B | `Longread_RNA_analysis/` |
| `cleaned_tryptic_peptides_under_151aa.csv` | 639 K | `Results/generate_supplemental_tables.py` |
| `cleaned_tryptic_peptides_detailed_under_151aa.csv` | 714 K | `Peptide_TMT_analysis/prosit/` |
| `cleaned_tryptic_peptides_detailed_under_151aa_with_SA.csv` | 997 K | `Microprotein_annotation_summary/alt_proteoform_records.py`; PROSIT-annotated output of the above |
| `ac_list_mapping.csv` | 87 K | `scRNAseq_summary_merging_analysis/` |
| `alt_proteoform_table.csv` | 550 K | `Microprotein_annotation_summary/Create_Alt_Proteoform_Files.py` |
| `codon_context_sequences.csv.gz` | 343 K | `Codon_context/` |
| `rbp_splice_per_junction.csv` | 9.0 M | RBP splice-junction figures |
| `uniprotkb_proteome_UP000005640_2026_07_13.tsv` | 14 M | `Results/microproteins_dashboard.py` |
| `blast_isoform_reference_sequences.fasta` | small | `Microprotein_annotation_summary/compute_nterm_peptide_substitutions.py` (cache; fetched from UniProt) |
| `blast_nterm_peptide_substitutions.csv` | small | `Results/microproteins_dashboard.py` (N-terminus filter) |
| `browser_tracks/` | 51 M | UCSC custom-track definition (`brain_tracks.txt`, not read by any script) plus the six bigWig/bigBed binaries it loads |

## UCSC genome-browser tracks

[`browser_tracks/brain_tracks.txt`](browser_tracks/brain_tracks.txt) is a UCSC
custom-track definition file — paste its contents into the browser's "add custom
tracks" box (hg38) to load six DLPFC tracks:

| Track | Type | What it shows |
|-------|------|---------------|
| DLPFC Ribo (Fw) / (Rev) | bigWig | adult-brain ribosome-profiling coverage, per strand |
| DLPFC ESPRESSO (CPM) | bigWig | long-read isoform coverage, mean CPM over 12 ROSMAP samples (sums overlapping isoforms per base) |
| DLPFC ESPRESSO Isoforms | bigBed | the same isoforms individually, coloured on log1p(mean CPM) |
| smORF PROSIT Peptides | bigGenePred | tryptic peptides projected to the genome, coloured by PROSIT confidence tier |
| smORF Rules | bigGenePred | smORFs coloured by decay class (NMD 50-nt rule, or non-stop decay) × main-ORF status |

Nothing in the pipeline reads this file — it is a browser configuration, kept
here so the track definitions are versioned alongside the data. The comments in
the file itself carry the display settings, colour-ramp rationale, and per-tier
feature counts.

Track data streams from the public dataset
[`brmiller/brain-microprotein-atlas`](https://huggingface.co/datasets/brmiller/brain-microprotein-atlas),
so nothing needs downloading first. The same six binaries are also kept in
[`browser_tracks/`](browser_tracks/); note that having them locally does **not**
make the track file work offline, since UCSC fetches `bigDataUrl` byte-ranges
server-side and a local path can never satisfy it. See
[`browser_tracks/README.md`](browser_tracks/README.md) for the hosting
requirements.

`smorf_trx_priority_color_class.bed` is the source BED behind the smORF Rules
track and is shipped here, alongside
`smorf_trx_priority_color_class_assignments.tsv` and
`smorf_nmd_mainorf_summary.tsv` (the per-smORF NMD/NSD × main-ORF calls the
colouring is derived from). The **assignments** TSV — not the summary — is what
the dashboard's ORF Rules column reads, via `load_decay_classes()`; its
`category` column carries the six classes, and `nsd_confident_tx` /
`nsd_gate_reason` record why each smORF did or did not qualify as non-stop
decay. Their generating scripts (`transcript_rules/`, `espresso/`,
`nmd_analysis/`) live outside this repo, in Box.

## Not in a clone

Excluded by `.gitignore` — either too large for GitHub or derived from
controlled-access data. Fetch them before running the modules that need them.

| File | Size | Where to get it |
|------|-----:|-----------------|
| `microprotein_master.csv` | 39 M | unzipped from `microprotein_master.zip` automatically |
| `combined_gpath_results.csv` | 84 M | Zenodo, `large_data_tables.tar.gz` |
| `tss_motif_prevalence.csv` | 75 M | Zenodo, `large_data_tables.tar.gz` |
| `counts.csv` | 21 M | **Restricted.** ROSMAP / MSBB individual-level data — requires AMP-AD / Synapse approval |

`swissprot_razor_peptides_detailed_under_151aa.csv` is an intermediate working
file kept out of the public repo; nothing in the pipeline reads it.

`blast_isoform_reference_sequences.fasta` and
`blast_nterm_peptide_substitutions.csv` are generated by
`Microprotein_annotation_summary/compute_nterm_peptide_substitutions.py`,
which BLASTs each N-terminal (aa 1-2) tryptic peptide against its matched
UniProt isoform's sequence and records any amino-acid substitutions — a
peptide that differs from the canonical protein at that position can't be a
bare Met-excision fragment of it, so the dashboard's N-terminus filter treats
1-2 substitutions as real evidence. The FASTA is a UniProt-fetched cache —
network access is only needed the first time a given accession is seen;
re-running the script only fetches what's missing, so a full cache makes both
it and the dashboard network-free.

See [../../DATA_AVAILABILITY.md](../../DATA_AVAILABILITY.md) and
[../../DATA_REQUIREMENTS.md](../../DATA_REQUIREMENTS.md) for full provenance and
column-level format specifications.
