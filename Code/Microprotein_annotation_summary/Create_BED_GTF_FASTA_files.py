import pandas as pd
import os
import re
import sys
import ast

# === 1. Load and filter MASTER file using gold standard criteria ===
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from gold_standard_filtering_criteria import load_and_filter_master

mp = load_and_filter_master(os.path.join(os.path.dirname(__file__), '..', 'data', 'microprotein_master.csv'))

# === 2. Keep only Salk (novel) microproteins ===
df = mp[mp['Database'] == 'Salk'].copy()

# === 3. Set annotation label based on which evidence is present ===
df['Annotation Status'] = df['has_MS'].map({True: 'MS', False: 'RiboCode_SAM'})

# === 3b. Derive best PROSIT confidence grade per microprotein ===
# The 'Confidence' column holds a bracketed list of per-peptide PROSIT
# spectral-angle tiers (e.g. "['Strong', 'Moderate']"). Collapse to the
# single best (highest-confidence) tier observed for each microprotein.
_confidence_rank = {'Strong': 0, 'Moderate': 1, 'Weak': 2, 'Insufficient': 3}


def best_confidence(val):
    if pd.isna(val):
        return None
    try:
        items = ast.literal_eval(val)
    except (ValueError, SyntaxError):
        return None
    if not isinstance(items, list):
        return None
    valid = [x for x in items if x in _confidence_rank]
    if not valid:
        return None
    return min(valid, key=lambda x: _confidence_rank[x])


df['PROSIT Grade'] = df['Confidence'].apply(best_confidence)

# Path to the full Ensembl+unreviewed combined source GTF (author-local; not
# committed to Github due to size). Override with the COMBINED_SOURCE_GTF env var.
COMBINED_SOURCE_GTF = os.environ.get(
    "COMBINED_SOURCE_GTF",
    "/Users/brendanmiller/Library/CloudStorage/Box-Box/brain_smorfs/annotation/"
    "2025-08-28_annotation_for_ribocode_plaqueomics_proteogenomics_DIA_shortstop/"
    "gencodev43_shortstop_plaqueomics_proteomics_ribocode_DIA_appended_brain_microproteins.gtf",
)

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
outdir = os.path.join(REPO_ROOT, 'GTF_and_BED_files')
os.makedirs(outdir, exist_ok=True)

# ============================================================
# Alternative proteoforms (non-ATG downstream initiation)
# ============================================================
# Merged into every file below alongside the canonical stop-to-stop products.
# They stay separable by the GTF source field (AltInitiation vs
# GTF2FastaPatched), the ALT_INIT FASTA field, and the _alt_initiation_N id
# suffix. Requires the combined source GTF, since real (spliced) coordinates
# are derived by trimming the parent's CDS blocks -- see alt_proteoform_records.
import alt_proteoform_records as apr

# One streaming pass over the source GTF serves both the BED12 block structure
# for annotated-start records and the alternative proteoform coordinates.
cds_index = {}
alt_records = []
if os.path.exists(COMBINED_SOURCE_GTF):
    print("\nIndexing CDS blocks from source GTF...")
    cds_index = apr.index_parent_cds(COMBINED_SOURCE_GTF,
                                     set(df['gene_symbol'].dropna().unique()))
    spliced_parents = sum(1 for v in cds_index.values() if len(v['blocks']) > 1)
    print(f"  {len(cds_index):,} records indexed, {spliced_parents:,} spliced "
          f"({100 * spliced_parents / max(len(cds_index), 1):.1f}%)")

    print("Building alternative proteoform records...")
    alt_records, alt_failures, alt_dropped = apr.build_alt_records(
        code_dir=os.path.join(os.path.dirname(__file__), '..'),
        source_gtf=COMBINED_SOURCE_GTF,
        gold_standard_ids=set(df['gene_symbol'].dropna().unique()),
        parent_meta=df.drop_duplicates(subset='gene_symbol').set_index('gene_symbol'),
        cds_index=cds_index,
    )
    print(f"  alternative proteoforms: {len(alt_records)} "
          f"({sum(1 for r in alt_records if len(r['blocks']) > 1)} spliced); "
          f"dropped {len(alt_dropped)} non-gold-standard; "
          f"{len(alt_failures)} validation failures")
    for pid, check, exp, obs in alt_failures[:10]:
        print(f"    FAIL {pid} {check}: expected {exp}, observed {obs}")
    graded = sum(1 for r in alt_records if r['prosit_grade'])
    print(f"  PROSIT grade from own downstream peptides: {graded} graded, "
          f"{len(alt_records) - graded} ungraded (upstream-only peptide evidence)")
else:
    print(f"\nWARNING: Source GTF not found at {COMBINED_SOURCE_GTF}")
    print("Alternative proteoforms need it for coordinates -- they will be OMITTED.")

alt_by_tier = {}
for _record in alt_records:
    alt_by_tier.setdefault(_record['prosit_grade'], []).append(_record)


# gene_symbol encodes the locus, e.g.
# "ENST00000448905.6+chrX:150983407-150985604_F:1_P:0". Used only as a fallback
# when the source GTF is unavailable -- it yields the outer envelope, which for a
# spliced record spans introns.
BED_COORD_RE = re.compile(r'^[^+-]+(?P<strand>[+-])(?P<chr>chr[\w]+):(?P<start>\d+)-(?P<end>\d+)')


def has_derivable_locus(symbol):
    """True if write_bed12 can place this record on the genome."""
    return symbol in cds_index or bool(BED_COORD_RE.match(str(symbol)))


def recolor(line, rgb):
    """Swap a BED12 line's itemRgb field (column 9). No-op when rgb is None."""
    if not rgb:
        return line
    fields = line.rstrip('\n').split('\t')
    fields[8] = rgb
    return '\t'.join(fields) + '\n'


def write_bed12(path, gene_symbols, alt_subset, extra_rows=(), colors=None,
                track_line=None):
    """Write a BED12 file with real CDS block structure.

    29.3% of these records are spliced, so the previous BED6 envelope described
    intervals that were ~93% intron by base count despite being labelled CDS.
    Blocks come from the source GTF; a record missing from it falls back to the
    single-interval envelope parsed out of gene_symbol, and is counted so the
    caller can warn.

    extra_rows appends already-formatted BED12 lines (see locus_rows) for
    records that carry no CDS blocks at all.

    colors is an optional {'cds': 'R,G,B', 'alt': 'R,G,B'} that fills the itemRgb
    column; extra_rows arrive pre-coloured. A track_line is written first when
    given -- UCSC ignores itemRgb without itemRgb="On" on the track.
    """
    colors = colors or {}
    written = fallback = 0
    seen = set()
    with open(path, 'w') as out:
        if track_line:
            out.write(track_line)
        for symbol in gene_symbols:
            if pd.isna(symbol) or symbol in seen:
                continue
            entry = cds_index.get(symbol)
            if entry is not None:
                out.write(recolor(apr.bed12_line_from_blocks(
                    entry['chrom'], entry['blocks'], symbol, entry['strand']),
                    colors.get('cds')))
            else:
                m = BED_COORD_RE.match(str(symbol))
                if not m:
                    continue          # no locus derivable -- FASTA/IDs only
                out.write(recolor(apr.bed12_line_from_blocks(
                    m.group('chr'), [(int(m.group('start')), int(m.group('end')))],
                    symbol, m.group('strand')), colors.get('cds')))
                fallback += 1
            seen.add(symbol)
            written += 1
        # Alternative-start records carry their own trimmed blocks -- never parse
        # their id, which embeds the PARENT's span.
        for record in alt_subset:
            out.write(recolor(apr.bed12_line(record), colors.get('alt')))
            written += 1
        for line in extra_rows:
            out.write(line)
            written += 1
    return written, fallback

# ============================================================
# BED12 — unreviewed microprotein CDS coordinates, with splice blocks
# Blocks come from the source GTF (1-based); BED is 0-based half-open,
# so start = GTF_start - 1.
# ============================================================

bed_file = os.path.join(outdir, 'Unreviewed_Brain_Microproteins_CDS_Absent_from_UniProt.bed')
n_bed, n_fallback = write_bed12(bed_file, df['gene_symbol'], alt_records)
print(f"BED12 file written to: {bed_file}")
print(f"Total entries: {n_bed} "
      f"({n_bed - len(alt_records)} annotated-start + {len(alt_records)} alternative-start)")
if n_fallback:
    print(f"  WARNING: {n_fallback} records fell back to a single-interval envelope "
          f"(no CDS blocks in source GTF)")
print("\nCounts of Annotation Status:")
print(df['Annotation Status'].value_counts())

# ============================================================
# BED12 — every non-Swiss-Prot microprotein, TrEMBL included
# ============================================================
# The BED above holds only records with block structure, which leaves out the
# TrEMBL half of the unreviewed set (Database "TrEMBL" in the master, folded into
# "Salk" by the canonical filter and marked smorf_type == "TrEMBL"). Those rows
# arrive with a plain gene name as gene_symbol and no transcript id, so nothing
# ties them to CDS blocks in the source GTF.
#
# What they do carry is genomic_coordinates: one chrN:start-end interval, no
# strand. That interval is a CDS *envelope* -- start codon to stop codon with the
# introns left in, not a gene body (only 7 of 1,115 equal an Ensembl gene span,
# while 536 equal some Ensembl transcript's CDS envelope exactly). Splicing is
# why it runs a median ~15x the coding length.
#
# The envelope is drawn as a full-width thick block, and the three record types
# are separated by colour instead -- itemRgb plus the track line below, so a
# browser distinguishes them at a glance and the colour survives into any tool
# that reads column 9. Colour is the only marker of "this thick block spans
# introns"; nothing in the geometry says so.

GENE_NAME_RE = re.compile(r'gene_name "([^"]+)"')
LOCUS_RE = re.compile(r'^(?P<chr>chr[\w.]+):(?P<start>\d+)-(?P<end>\d+)$')


def index_gene_strands(source_gtf, gene_names):
    """gene_name -> [(chrom, start, end, strand)] from the Ensembl half of the GTF."""
    index = {}
    with open(source_gtf) as src:
        for line in src:
            if line.startswith('#'):
                continue
            fields = line.split('\t')
            if (len(fields) < 9 or fields[2] != 'gene'
                    or fields[1] not in ('HAVANA', 'ENSEMBL')):
                continue
            m = GENE_NAME_RE.search(fields[8])
            if m and m.group(1) in gene_names:
                index.setdefault(m.group(1), []).append(
                    (fields[0], int(fields[3]), int(fields[4]), fields[6]))
    return index


# itemRgb per record type. Purple matches the BED badge in the repo README.
RGB_CDS = '31,78,121'        # annotated-start, real CDS blocks
RGB_ALT = '130,80,223'       # alternative-start proteoform, real CDS blocks
RGB_TREMBL = '217,119,6'     # TrEMBL, single-interval CDS envelope (introns in)

TRACK_LINE = (
    'track name="Brain microproteins (non-Swiss-Prot)" '
    'description="Unreviewed brain microproteins: blue=annotated-start CDS, '
    'purple=alternative-start CDS, amber=TrEMBL CDS envelope (spans introns)" '
    'itemRgb="On" visibility=pack\n'
)


def locus_bed12_line(chrom, start, end, name, strand, rgb=RGB_TREMBL):
    """BED12 for a CDS envelope: one block, thick across it, coloured as such."""
    chrom_start = start - 1                      # master is 1-based, BED is 0-based
    return '\t'.join([
        chrom, str(chrom_start), str(end), name, '0', strand,
        str(chrom_start), str(end), rgb, '1',
        f'{end - chrom_start},', '0,',
    ]) + '\n'


def locus_rows(rows, gene_strands):
    """BED12 lines for microproteins carrying only a single-interval CDS envelope."""
    lines, unplaced, unstranded = [], 0, 0
    for _, row in rows.iterrows():
        m = LOCUS_RE.match(str(row['genomic_coordinates']))
        if not m:
            unplaced += 1                        # no coordinates in the master at all
            continue
        chrom, start, end = m.group('chr'), int(m.group('start')), int(m.group('end'))
        # gene_name is the only handle on strand. Require the locus to sit on the
        # matching gene so a symbol reused elsewhere cannot mislabel it.
        strand = '.'
        for g_chrom, g_start, g_end, g_strand in gene_strands.get(row['gene_name'], []):
            if g_chrom == chrom and g_start <= end and start <= g_end:
                strand = g_strand
                break
        if strand == '.':
            unstranded += 1
        lines.append(locus_bed12_line(chrom, start, end, row['gene_id'], strand))
    return lines, unplaced, unstranded


no_locus = df[~df['gene_symbol'].fillna('').map(has_derivable_locus)].copy()
gene_strands = {}
if os.path.exists(COMBINED_SOURCE_GTF) and not no_locus.empty:
    print("\nIndexing Ensembl gene strands for the block-less records...")
    gene_strands = index_gene_strands(COMBINED_SOURCE_GTF,
                                      set(no_locus['gene_name'].dropna()))
locus_lines, n_unplaced, n_unstranded = locus_rows(no_locus, gene_strands)

all_bed_file = os.path.join(outdir, 'Unreviewed_Brain_Microproteins_All_Non_SwissProt.bed')
n_all, _ = write_bed12(all_bed_file, df['gene_symbol'], alt_records,
                       extra_rows=locus_lines,
                       colors={'cds': RGB_CDS, 'alt': RGB_ALT},
                       track_line=TRACK_LINE)
print(f"\nAll-non-Swiss-Prot BED12 written to: {all_bed_file}")
print(f"Total entries: {n_all} "
      f"({n_all - len(alt_records) - len(locus_lines)} annotated-start CDS + "
      f"{len(alt_records)} alternative-start CDS + {len(locus_lines)} CDS envelopes "
      f"[{(no_locus['smorf_type'] == 'TrEMBL').sum()} TrEMBL rows in])")
if not gene_strands and not no_locus.empty:
    print("  WARNING: source GTF unavailable — every CDS envelope is strandless ('.')")
elif n_unstranded:
    print(f"  {n_unstranded} CDS envelopes left strandless ('.') — no matching "
          f"Ensembl gene at that locus")
if n_unplaced:
    print(f"  {n_unplaced} rows omitted — no genomic_coordinates in the master")

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
    for record in alt_records:
        f.write(f"{record['alt_protein_id']}\n")
print(f"Gene IDs file written to: {gene_ids_file} "
      f"(+{len(alt_records)} alternative proteoforms)")

# ============================================================
# Genomic coordinates text file
# ============================================================
coords = df['genomic_coordinates'].dropna().unique()
# The 186 native/M-substituted variant pairs share one locus, so dedupe.
alt_coords = [c for c in dict.fromkeys(r['genomic_coordinates'] for r in alt_records)
              if c not in set(coords)]
coords_file = os.path.join(outdir, 'Unreviewed_Brain_Microproteins_genomic_coordinates.txt')
with open(coords_file, 'w') as f:
    for c in list(coords) + alt_coords:
        f.write(f"{c}\n")
print(f"Genomic coordinates written to: {coords_file} "
      f"({len(coords) + len(alt_coords)} entries, +{len(alt_coords)} alternative)")

# ============================================================
# Coordinate → sequence mapping
# ============================================================
mapping = df[['genomic_coordinates', 'sequence']].dropna(subset=['genomic_coordinates', 'sequence']).drop_duplicates()
if alt_records:
    mapping = pd.concat([mapping, pd.DataFrame(
        [{'genomic_coordinates': r['genomic_coordinates'], 'sequence': r['sequence']}
         for r in alt_records]
    )], ignore_index=True).drop_duplicates()
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
    for record in alt_records:
        f.write(f"{apr.fasta_header(record)}\n{record['sequence']}\n")
print(f"FASTA file written to: {fasta_file} "
      f"({len(df) + len(alt_records)} records, +{len(alt_records)} alternative proteoforms)")

# ============================================================
# GTF — unreviewed CDS filtered to gold standard
# ============================================================
# The transcript_id in GTF2FastaPatched entries matches gene_symbol in master CSV.
# Build a lookup set of all gold-standard gene_symbols.
gold_standard_ids = set(df['gene_symbol'].dropna().unique())

unreviewed_gtf_file = os.path.join(outdir, 'Unreviewed_Brain_Microproteins_Absent_from_UniProt.gtf')

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
        alt_gtf_lines = 0
        for record in alt_records:
            for gtf_line in apr.gtf_lines(record):
                out.write(gtf_line)
                alt_gtf_lines += 1
    print(f"\nUnreviewed GTF written to: {unreviewed_gtf_file} "
          f"({kept + alt_gtf_lines:,} lines: {kept:,} canonical + "
          f"{alt_gtf_lines:,} alternative proteoform CDS)")

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
        out.write('# Alternative proteoform source: AltInitiation (non-ATG downstream initiation; CDS features only)\n')
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
        alt_combined_lines = 0
        for record in alt_records:
            for gtf_line in apr.gtf_lines(record):
                out.write(gtf_line)
                alt_combined_lines += 1

    print(f"Combined GTF written to: {combined_gtf_file}")
    print(f"  Ensembl (HAVANA + ENSEMBL) lines: {ensembl_lines:,}")
    print(f"  Unreviewed CDS lines:            {unreviewed_lines:,}")
    print(f"  Alternative proteoform CDS lines: {alt_combined_lines:,}")
    print(f"  Total lines:                       "
          f"{ensembl_lines + unreviewed_lines + alt_combined_lines:,}")

# ============================================================
# Per-PROSIT-confidence-tier files (Strong / Moderate / Weak / Insufficient)
# Splits the unreviewed microproteins into individual FASTA/BED/GTF files
# based on their best PROSIT spectral-angle confidence grade.
# ============================================================
print("\n" + "=" * 60)
print("Per-PROSIT-confidence-tier files")
print("=" * 60)
print("PROSIT Grade counts:")
print(df['PROSIT Grade'].value_counts(dropna=False))

tier_dir = os.path.join(outdir, 'PROSIT_confidence_tiers')
os.makedirs(tier_dir, exist_ok=True)

TIERS = ['Strong', 'Moderate', 'Weak', 'Insufficient']
tier_to_ids = {}

for tier in TIERS:
    tier_df = df[df['PROSIT Grade'] == tier].copy()
    tier_alt = alt_by_tier.get(tier, [])
    if tier_df.empty and not tier_alt:
        print(f"\n[{tier}] no microproteins — skipping")
        continue

    # FASTA
    tier_fasta = os.path.join(tier_dir, f'Unreviewed_Brain_Microproteins_{tier}.fasta')
    with open(tier_fasta, 'w') as f:
        for _, row in tier_df.iterrows():
            header = (f">{row['gene_id']}|{row['gene_name']}|{row['Database']}|"
                      f"{row['Annotation Status']}|PROSIT={tier}")
            f.write(f"{header}\n{row['sequence']}\n")
        for record in tier_alt:
            f.write(f"{apr.fasta_header(record)}|PROSIT={tier}\n{record['sequence']}\n")

    # BED12 with real CDS block structure (same convention as the main BED above)
    tier_bed = os.path.join(tier_dir, f'Unreviewed_Brain_Microproteins_CDS_{tier}.bed')
    n_tbed, _ = write_bed12(tier_bed, tier_df['gene_symbol'], tier_alt)

    # ID list
    tier_ids = os.path.join(tier_dir, f'Unreviewed_Brain_Microproteins_IDs_{tier}.txt')
    with open(tier_ids, 'w') as f:
        for gid in tier_df['gene_id'].dropna().unique():
            f.write(f"{str(gid).replace(chr(34), '').replace(chr(39), '').strip()}\n")
        for record in tier_alt:
            f.write(f"{record['alt_protein_id']}\n")

    tier_to_ids[tier] = set(tier_df['gene_symbol'].dropna().unique())
    print(f"\n[{tier}] {len(tier_df) + len(tier_alt)} microproteins "
          f"({len(tier_df)} canonical + {len(tier_alt)} alternative proteoforms)")
    print(f"  FASTA → {tier_fasta}")
    print(f"  BED12 → {tier_bed} ({n_tbed} rows)")
    print(f"  IDs   → {tier_ids}")

# Per-tier GTF — single pass over the combined source GTF
if tier_to_ids and os.path.exists(COMBINED_SOURCE_GTF):
    id_to_tier = {gid: tier for tier, ids in tier_to_ids.items() for gid in ids}
    # A tier may have alternative proteoforms but no canonical members, so open
    # a handle for every tier that has either.
    tiers_with_output = set(tier_to_ids) | {t for t in alt_by_tier if t in TIERS}
    handles = {
        tier: open(
            os.path.join(tier_dir, f'Unreviewed_Brain_Microproteins_{tier}.gtf'), 'w'
        )
        for tier in tiers_with_output
    }
    counts = {tier: 0 for tier in tiers_with_output}
    alt_counts = {tier: 0 for tier in tiers_with_output}
    with open(COMBINED_SOURCE_GTF, 'r') as src:
        for line in src:
            if line.startswith('#'):
                continue
            fields = line.split('\t')
            if len(fields) < 9 or fields[1] != 'GTF2FastaPatched':
                continue
            attrs = fields[8]
            tid_start = attrs.find('transcript_id "')
            if tid_start == -1:
                continue
            tid_start += len('transcript_id "')
            tid_end = attrs.find('"', tid_start)
            tier = id_to_tier.get(attrs[tid_start:tid_end])
            if tier:
                handles[tier].write(line)
                counts[tier] += 1
    for tier in tiers_with_output:
        for record in alt_by_tier.get(tier, []):
            for gtf_line in apr.gtf_lines(record):
                handles[tier].write(gtf_line)
                alt_counts[tier] += 1
    for tier, h in handles.items():
        h.close()
        print(f"[{tier}] GTF → Unreviewed_Brain_Microproteins_{tier}.gtf "
              f"({counts[tier] + alt_counts[tier]:,} lines: {counts[tier]:,} canonical + "
              f"{alt_counts[tier]:,} alternative)")
elif tier_to_ids:
    print(f"\nWARNING: Source GTF not found at {COMBINED_SOURCE_GTF}")
    print("Skipping per-tier GTF generation. Provide the combined source GTF to regenerate.")

# ============================================================
# Ribo-ShortStop file set (RiboCode + ShortStop SAM, no DDA mass-spec)
# Matches the "Ribo-Seq Only" group in the Figure 7 Panel J detection-
# evidence breakdown: microproteins with RiboSAM translation support that
# were not detected by DDA mass spectrometry (has_RiboSAM & ~DDA_evidence).
# ============================================================
print("\n" + "=" * 60)
print("Ribo-ShortStop file set")
print("=" * 60)

ribo_dir = os.path.join(outdir, 'Ribo_ShortStop')
os.makedirs(ribo_dir, exist_ok=True)

# No alternative proteoforms appear here: every alt parent is DDA-detected
# (DDA_evidence True, has_RiboSAM False for all of them), so the
# has_RiboSAM & ~DDA_evidence gate excludes them by construction.
ribo_df = df[df['has_RiboSAM'] & ~df['DDA_evidence']].copy()
print(f"Ribo-ShortStop microproteins: {len(ribo_df)} "
      f"(alternative proteoforms excluded by definition — all are DDA-detected)")

if not ribo_df.empty:
    # FASTA
    ribo_fasta = os.path.join(ribo_dir, 'Unreviewed_Brain_Microproteins_Ribo_ShortStop.fasta')
    with open(ribo_fasta, 'w') as f:
        for _, row in ribo_df.iterrows():
            header = (f">{row['gene_id']}|{row['gene_name']}|{row['Database']}|"
                      f"{row['Annotation Status']}")
            f.write(f"{header}\n{row['sequence']}\n")

    # BED12 with real CDS block structure (no alternative-start records here)
    ribo_bed = os.path.join(ribo_dir, 'Unreviewed_Brain_Microproteins_CDS_Ribo_ShortStop.bed')
    n_rbed, _ = write_bed12(ribo_bed, ribo_df['gene_symbol'], [])

    # ID list
    ribo_ids_file = os.path.join(ribo_dir, 'Unreviewed_Brain_Microproteins_IDs_Ribo_ShortStop.txt')
    with open(ribo_ids_file, 'w') as f:
        for gid in ribo_df['gene_id'].dropna().unique():
            f.write(f"{str(gid).replace(chr(34), '').replace(chr(39), '').strip()}\n")

    print(f"  FASTA → {ribo_fasta}")
    print(f"  BED12 → {ribo_bed} ({n_rbed} rows)")
    print(f"  IDs   → {ribo_ids_file}")

    # GTF — single pass over the combined source GTF
    ribo_ids = set(ribo_df['gene_symbol'].dropna().unique())
    if os.path.exists(COMBINED_SOURCE_GTF):
        ribo_gtf = os.path.join(ribo_dir, 'Unreviewed_Brain_Microproteins_Ribo_ShortStop.gtf')
        ribo_lines = 0
        with open(COMBINED_SOURCE_GTF, 'r') as src, open(ribo_gtf, 'w') as out:
            for line in src:
                if line.startswith('#'):
                    continue
                fields = line.split('\t')
                if len(fields) < 9 or fields[1] != 'GTF2FastaPatched':
                    continue
                attrs = fields[8]
                tid_start = attrs.find('transcript_id "')
                if tid_start == -1:
                    continue
                tid_start += len('transcript_id "')
                tid_end = attrs.find('"', tid_start)
                if attrs[tid_start:tid_end] in ribo_ids:
                    out.write(line)
                    ribo_lines += 1
        print(f"  GTF   → {ribo_gtf} ({ribo_lines:,} lines)")
    else:
        print(f"\nWARNING: Source GTF not found at {COMBINED_SOURCE_GTF}")
        print("Skipping Ribo-ShortStop GTF generation. Provide the combined source GTF to regenerate.")
