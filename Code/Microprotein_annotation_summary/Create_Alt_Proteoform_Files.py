"""
Standalone export of alternative proteoforms (non-ATG alternative initiation).

Writes GTF_and_BED_files/Alt_Proteoforms/ -- the alternative proteoforms on
their own, including a BED12 with real splice-block structure that the merged
release files (BED6) cannot represent.

These same records are ALSO merged into the main release file set by
Create_BED_GTF_FASTA_files.py. Both scripts build them through
alt_proteoform_records.build_alt_records, so the two cannot drift.

NOTE: do not `import` Create_BED_GTF_FASTA_files -- it has no
`if __name__ == "__main__"` guard, so importing it re-runs the entire exporter
and rewrites every file in GTF_and_BED_files/.
"""

import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from gold_standard_filtering_criteria import load_and_filter_master

import alt_proteoform_records as apr

CODE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))

# Author-local combined source GTF (not committed, too large). Override with
# the COMBINED_SOURCE_GTF env var -- same contract as the main exporter.
COMBINED_SOURCE_GTF = os.environ.get(
    "COMBINED_SOURCE_GTF",
    "/Users/brendanmiller/Library/CloudStorage/Box-Box/brain_smorfs/annotation/"
    "2025-08-28_annotation_for_ribocode_plaqueomics_proteogenomics_DIA_shortstop/"
    "gencodev43_shortstop_plaqueomics_proteomics_ribocode_DIA_appended_brain_microproteins.gtf",
)

outdir = os.path.join(REPO_ROOT, 'GTF_and_BED_files', 'Alt_Proteoforms')
os.makedirs(outdir, exist_ok=True)

print("=" * 70)
print("Alternative proteoform export (non-ATG alternative initiation)")
print("=" * 70)

if not os.path.exists(COMBINED_SOURCE_GTF):
    print(f"\nWARNING: Source GTF not found at {COMBINED_SOURCE_GTF}")
    print("Coordinates cannot be derived without it. Nothing written.")
    sys.exit(0)

mp = load_and_filter_master(os.path.join(CODE_DIR, 'data', 'microprotein_master.csv'))
salk = mp[mp['Database'] == 'Salk'].copy()
salk['Annotation Status'] = salk['has_MS'].map({True: 'MS', False: 'RiboCode_SAM'})

records, failures, dropped = apr.build_alt_records(
    code_dir=CODE_DIR,
    source_gtf=COMBINED_SOURCE_GTF,
    gold_standard_ids=set(salk['gene_symbol'].dropna().unique()),
    parent_meta=salk.drop_duplicates(subset='gene_symbol').set_index('gene_symbol'),
)

spliced = sum(1 for r in records if len(r['blocks']) > 1)
print(f"\nRecords built: {len(records)}  ({spliced} spliced, "
      f"{len(records) - spliced} single-block)")
print(f"Total CDS lines: {sum(len(r['blocks']) for r in records):,}")
print(f"Dropped (not gold standard): {len(dropped)} rows / "
      f"{dropped['orf_id'].nunique() if len(dropped) else 0} parents")

if failures:
    print(f"\nVALIDATION FAILURES: {len(failures)} (records excluded)")
    for pid, check, exp, obs in failures[:20]:
        print(f"  {pid}  {check}: expected {exp}, observed {obs}")
else:
    print("Validation: all records passed (envelope, coding length, "
          "chrom/strand, block order)")

if not records:
    print("\nNo records to write. Exiting.")
    sys.exit(0)

# --- FASTA ---
fasta_file = os.path.join(outdir, 'Alt_Proteoforms.fasta')
with open(fasta_file, 'w') as f:
    for r in records:
        f.write(f"{apr.fasta_header(r)}\n{r['sequence']}\n")
print(f"\nFASTA written to: {fasta_file} ({len(records)} records)")

# --- ID list ---
ids_file = os.path.join(outdir, 'Alt_Proteoform_IDs.txt')
with open(ids_file, 'w') as f:
    for r in records:
        f.write(f"{r['alt_protein_id']}\n")
print(f"IDs written to: {ids_file}")

# --- mapping TSV ---
mapping_file = os.path.join(outdir, 'Alt_Proteoforms_mapping.tsv')
pd.DataFrame([{
    'alt_protein_id': r['alt_protein_id'],
    'parent_orf_id': r['parent_orf_id'],
    'parent_gene_id': r['parent_gene_id'],
    'gene_name': r['gene_name'],
    'chromosome': r['chrom'],
    'strand': r['strand'],
    'cds_start': r['blocks'][0][0],
    'cds_end': r['blocks'][-1][1],
    'n_cds_blocks': len(r['blocks']),
    'aa_position': r['aa_position'],
    'aa_length': r['aa_length'],
    'codon': r['codon'],
    'initiation_tier': r['initiation_tier'],
    'cognate_status': r['cognate_status'],
    'context_strength': r['context_strength'],
    'initiation_variant': r['initiation_variant'],
    'includes_stop': r['includes_stop'],
    'prosit_grade': r['prosit_grade'],
    'n_support_peptides': r['n_support_peptides'],
    'sequence': r['sequence'],
} for r in records]).to_csv(mapping_file, sep='\t', index=False)
print(f"Mapping written to: {mapping_file}")

# --- GTF (CDS features only) ---
gtf_file = os.path.join(outdir, 'Alt_Proteoforms.gtf')
gtf_line_count = 0
with open(gtf_file, 'w') as f:
    for r in records:
        for line in apr.gtf_lines(r):
            f.write(line)
            gtf_line_count += 1
print(f"GTF written to: {gtf_file} ({gtf_line_count:,} CDS lines)")

# --- BED12 with real block structure ---
bed_file = os.path.join(outdir, 'Alt_Proteoforms_CDS.bed')
with open(bed_file, 'w') as f:
    for r in records:
        f.write(apr.bed12_line(r))
print(f"BED12 written to: {bed_file} ({len(records)} records)")

# The former separate .bed12 is now redundant -- .bed carries the blocks.
_legacy = os.path.join(outdir, 'Alt_Proteoforms_CDS.bed12')
if os.path.exists(_legacy):
    os.remove(_legacy)
    print(f"Removed superseded {os.path.basename(_legacy)} (blocks now in .bed)")

# ============================================================
# Summary
# ============================================================
print("\n" + "=" * 70)
print("Summary")
print("=" * 70)
print(f"Alternative proteoforms written : {len(records)}")
print(f"  parent ORFs                   : {len(set(r['parent_orf_id'] for r in records))}")
print(f"  spliced (>1 CDS block)        : {spliced}")
print(f"  GTF CDS lines                 : {gtf_line_count:,}")
print(f"Validation failures             : {len(failures)}")
print(f"\nOutput directory: {outdir}")

graded = pd.Series([r['prosit_grade'] for r in records])
print("\nPROSIT grade (derived from each proteoform's own downstream peptides):")
print(graded.value_counts(dropna=False).to_string())
print("\nInitiation tier (from alt_proteoform_table):")
print(pd.Series([r['initiation_tier'] for r in records]).value_counts().sort_index().to_string())
