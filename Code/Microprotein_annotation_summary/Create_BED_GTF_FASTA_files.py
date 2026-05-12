import pandas as pd
import os
import sys

# === 1. Load and filter MASTER file using gold standard criteria ===
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from gold_standard_filtering_criteria import load_and_filter_master

mp = load_and_filter_master(os.path.join(os.path.dirname(__file__), '..', 'data', 'microprotein_master.csv'))

# === 2. Keep only Salk (novel) microproteins ===
df = mp[mp['Database'] == 'Salk'].copy()

# === 3. Set annotation label based on which evidence is present ===
df['Annotation Status'] = df['has_MS'].map({True: 'MS', False: 'RiboCode_SAM'})

# Path to the full Ensembl+unreviewed combined source GTF (not committed to Github due to size)
COMBINED_SOURCE_GTF = (
    "/Users/brendanmiller/Library/CloudStorage/Box-Box/brain_smorfs/annotation/"
    "2025-08-28_annotation_for_ribocode_plaqueomics_proteogenomics_DIA_shortstop/"
    "gencodev43_shortstop_plaqueomics_proteomics_ribocode_DIA_appended_brain_microproteins.gtf"
)

outdir = '/Users/brendanmiller/Library/CloudStorage/Box-Box/brendan_alan_shared_folder/AD_paper/Github/GTF_and_BED_files'
os.makedirs(outdir, exist_ok=True)

# ============================================================
# BED6 — unreviewed microprotein CDS coordinates
# NOTE: gene_symbol encodes GTF (1-based) coordinates.
# BED is 0-based half-open, so start = GTF_start - 1.
# ============================================================

bed = df['gene_symbol'].str.extract(
    r'^[^+-]+(?P<strand>[+-])(?P<chr>chr[\w]+):(?P<start>\d+)-(?P<end>\d+)'
)
bed['name'] = df['gene_symbol'].values
bed['score'] = '0'
bed = bed[['chr', 'start', 'end', 'name', 'score', 'strand']]
bed = bed.dropna(subset=['chr', 'start', 'end', 'strand'])
bed['start'] = bed['start'].astype(int) - 1   # convert GTF 1-based → BED 0-based
bed['end']   = bed['end'].astype(int)

bed_file = os.path.join(outdir, 'Unreviewed_Brain_Microproteins_CDS_absent_from_UniProt.bed')
bed.to_csv(bed_file, sep='\t', header=False, index=False)
print(f"BED file written to: {bed_file}")
print(f"Total entries: {len(bed)}")
print("\nCounts of Annotation Status:")
print(df['Annotation Status'].value_counts())
print(f"\nFirst few BED entries:")
print(bed.head())

# ============================================================
# Gene IDs text file
# ============================================================
gene_ids_file = os.path.join(outdir, 'Unreviewed_Brain_Microproteins_IDs.txt')
df['gene_id'].dropna().unique().tofile(gene_ids_file, sep='\n')
with open(gene_ids_file, "r") as f:
    lines = f.readlines()
with open(gene_ids_file, "w") as f:
    for line in lines:
        f.write(line.replace('"', '').replace("'", "").strip() + "\n")
print(f"Gene IDs file written to: {gene_ids_file}")

# ============================================================
# Genomic coordinates text file
# ============================================================
coords = df['genomic_coordinates'].dropna().unique()
coords_file = os.path.join(outdir, 'Unreviewed_Brain_Microproteins_genomic_coordinates.txt')
with open(coords_file, 'w') as f:
    for c in coords:
        f.write(f"{c}\n")
print(f"Genomic coordinates written to: {coords_file} ({len(coords)} entries)")

# ============================================================
# Coordinate → sequence mapping
# ============================================================
mapping = df[['genomic_coordinates', 'sequence']].dropna(subset=['genomic_coordinates', 'sequence']).drop_duplicates()
mapping_file = os.path.join(outdir, 'Unreviewed_Brain_Microproteins_mapping_coordinates_to_sequences.tsv')
mapping.to_csv(mapping_file, sep='\t', index=False)
print(f"Coordinate-to-sequence mapping written to: {mapping_file} ({len(mapping)} entries)")

# ============================================================
# FASTA — protein sequences
# ============================================================
fasta_file = os.path.join(outdir, 'Unreviewed_Brain_Microproteins.fasta')
with open(fasta_file, 'w') as f:
    for _, row in df.iterrows():
        header = f">{row['gene_id']}|{row['gene_name']}|{row['Database']}|{row['Annotation Status']}"
        f.write(f"{header}\n{row['sequence']}\n")
print(f"FASTA file written to: {fasta_file}")

# ============================================================
# GTF — unreviewed CDS filtered to gold standard
# ============================================================
# The transcript_id in GTF2FastaPatched entries matches gene_symbol in master CSV.
# Build a lookup set of all gold-standard gene_symbols.
gold_standard_ids = set(df['gene_symbol'].dropna().unique())

unreviewed_gtf_file = os.path.join(outdir, 'Unreviewed_Brain_Microproteins_absent_from_UniProt.gtf')

if not os.path.exists(COMBINED_SOURCE_GTF):
    print(f"\nWARNING: Source GTF not found at {COMBINED_SOURCE_GTF}")
    print("Skipping GTF generation. Provide the combined source GTF to regenerate.")
else:
    kept = 0
    with open(COMBINED_SOURCE_GTF, 'r') as src, open(unreviewed_gtf_file, 'w') as out:
        for line in src:
            if line.startswith('#'):
                continue
            fields = line.split('\t')
            if len(fields) < 9:
                continue
            source = fields[1]
            if source != 'GTF2FastaPatched':
                continue
            # Extract transcript_id from attributes field
            attrs = fields[8]
            tid_start = attrs.find('transcript_id "')
            if tid_start == -1:
                continue
            tid_start += len('transcript_id "')
            tid_end = attrs.find('"', tid_start)
            transcript_id = attrs[tid_start:tid_end]
            if transcript_id in gold_standard_ids:
                out.write(line)
                kept += 1
    print(f"\nUnreviewed GTF written to: {unreviewed_gtf_file} ({kept:,} lines)")

    # ============================================================
    # Combined GTF — Ensembl (HAVANA + ENSEMBL) + filtered unreviewed CDS
    # ============================================================
    combined_gtf_file = os.path.join(outdir, 'Ensembl_and_Unreviewed_Brain_Microproteins.gtf')
    ensembl_lines = 0
    unreviewed_lines = 0

    with open(COMBINED_SOURCE_GTF, 'r') as src, open(combined_gtf_file, 'w') as out:
        out.write('# Combined GTF: GENCODEv43/Ensembl annotations + gold-standard unreviewed brain microprotein features (transcript, exon, CDS)\n')
        out.write('# Ensembl sources: HAVANA, ENSEMBL\n')
        out.write('# Unreviewed source: GTF2FastaPatched (filtered to gold standard evidence criteria)\n')
        for line in src:
            if line.startswith('#'):
                continue
            fields = line.split('\t')
            if len(fields) < 9:
                continue
            source = fields[1]
            if source in ('HAVANA', 'ENSEMBL'):
                out.write(line)
                ensembl_lines += 1
            elif source == 'GTF2FastaPatched':
                attrs = fields[8]
                tid_start = attrs.find('transcript_id "')
                if tid_start == -1:
                    continue
                tid_start += len('transcript_id "')
                tid_end = attrs.find('"', tid_start)
                transcript_id = attrs[tid_start:tid_end]
                if transcript_id in gold_standard_ids:
                    out.write(line)
                    unreviewed_lines += 1

    print(f"Combined GTF written to: {combined_gtf_file}")
    print(f"  Ensembl (HAVANA + ENSEMBL) lines: {ensembl_lines:,}")
    print(f"  Unreviewed CDS lines:            {unreviewed_lines:,}")
    print(f"  Total lines:                       {ensembl_lines + unreviewed_lines:,}")
