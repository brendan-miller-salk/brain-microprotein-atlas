"""
Verify GTF_and_BED_files/ against the master table and the alternative
proteoform table.

Checks that every released GTF/BED/FASTA file contains exactly the records the
canonical filter says it should, that the canonical and alternative-proteoform
halves are internally consistent, and that the alternative records' coordinates
actually encode their sequences.

Run from anywhere:
    python Code/verify_gtf_bed_fasta_files.py

Exits 0 if every check passes, 1 otherwise. The 1.1 GB combined source GTF is
NOT required -- the deep coordinate checks are self-contained (GTF blocks are
validated against FASTA sequence length and against the BED12 block structure).
Pass --source-gtf to additionally re-derive alternative coordinates from the
parent CDS blocks.

Companion to verify_dashboard_vs_figures.py, which reconciles a different pair
(dashboard filtering vs manuscript Source Data).
"""

import argparse
import os
import re
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gold_standard_filtering_criteria import load_and_filter_master

CODE_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(CODE_DIR, '..'))
FILES_DIR = os.path.join(REPO_ROOT, 'GTF_and_BED_files')
TIER_DIR = os.path.join(FILES_DIR, 'PROSIT_confidence_tiers')
RIBO_DIR = os.path.join(FILES_DIR, 'Ribo_ShortStop')
ALT_DIR = os.path.join(FILES_DIR, 'Alt_Proteoforms')

TIERS = ['Strong', 'Moderate', 'Weak', 'Insufficient']
ALT_MARKER = '_alt_initiation_'

# gene_symbol encodes the locus as e.g. "ENST00000448905.6+chrX:150983407-150985604_F:1_P:0".
# Entries whose gene_symbol is a plain gene name cannot yield coordinates.
COORD_RE = re.compile(r'^[^+-]+(?P<strand>[+-])(?P<chr>chr[\w]+):(?P<start>\d+)-(?P<end>\d+)')

# TrEMBL rows carry only the gene-level interval in genomic_coordinates.
LOCUS_RE = re.compile(r'^(?P<chr>chr[\w.]+):(?P<start>\d+)-(?P<end>\d+)$')

results = []


def check(name, passed, detail=''):
    results.append((name, passed, detail))
    print(f"  [{'PASS' if passed else 'FAIL'}] {name}" + (f" — {detail}" if detail else ''))
    return passed


def info(name, detail):
    print(f"  [INFO] {name} — {detail}")


# ============================================================
# Readers
# ============================================================

def read_fasta(path):
    """Return {record_id: (full_header, sequence)}. Sequences are unwrapped."""
    records, header = {}, None
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                header = line
            elif header is not None:
                records[header[1:].split('|')[0]] = (header, line)
                header = None
    return records


def read_bed6(path):
    rows = []
    with open(path) as f:
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) >= 6:
                rows.append({'chrom': p[0], 'start': int(p[1]), 'end': int(p[2]),
                             'name': p[3], 'strand': p[5]})
    return rows


def read_bed12(path):
    rows = []
    with open(path) as f:
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) < 12:
                continue
            sizes = [int(x) for x in p[10].rstrip(',').split(',')]
            starts = [int(x) for x in p[11].rstrip(',').split(',')]
            rows.append({'chrom': p[0], 'start': int(p[1]), 'end': int(p[2]),
                         'name': p[3], 'strand': p[5], 'block_count': int(p[9]),
                         'sizes': sizes, 'starts': starts})
    return rows


def read_gtf(path):
    """Return {transcript_id: {'source', 'chrom', 'strand', 'blocks', 'attrs'}}."""
    entries = {}
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            p = line.rstrip('\n').split('\t')
            if len(p) < 9:
                continue
            attrs = dict(re.findall(r'(\w+) "([^"]*)";', p[8]))
            tid = attrs.get('transcript_id')
            if tid is None:
                continue
            e = entries.setdefault(tid, {'source': p[1], 'chrom': p[0], 'strand': p[6],
                                         'blocks': [], 'attrs': attrs, 'features': set()})
            e['features'].add(p[2])
            if p[2] == 'CDS':
                e['blocks'].append((int(p[3]), int(p[4]), p[7]))
    return entries


def read_ids(path):
    with open(path) as f:
        return [l.strip() for l in f if l.strip()]


def split_alt(items, key=lambda x: x):
    alt = [i for i in items if ALT_MARKER in key(i)]
    canonical = [i for i in items if ALT_MARKER not in key(i)]
    return canonical, alt


# ============================================================
# Expected state, recomputed from source of truth
# ============================================================

def build_expectations():
    mp = load_and_filter_master(os.path.join(CODE_DIR, 'data', 'microprotein_master.csv'))
    df = mp[mp['Database'] == 'Salk'].copy()

    coords = df['gene_symbol'].str.extract(COORD_RE)
    has_coords = coords[['chr', 'start', 'end', 'strand']].notna().all(axis=1)

    alt = pd.read_csv(os.path.join(CODE_DIR, 'data', 'alt_proteoform_table.csv'))
    alt = alt[alt['alt_protein'] == 'Yes']
    alt = alt[alt['orf_id'].isin(set(df['gene_symbol']))]

    # Rows the CDS BED cannot place — TrEMBL entries — reach the all-non-Swiss-Prot
    # BED as gene-level loci keyed by gene_id.
    locus_only = df.loc[~has_coords]
    locus_m = locus_only['genomic_coordinates'].astype(str).str.extract(LOCUS_RE)
    locus_ok = locus_m.notna().all(axis=1)

    return {
        'df': df,
        'n_rows': len(df),
        'locus_coords': {
            gid: (m_chr, int(m_start), int(m_end))
            for gid, m_chr, m_start, m_end in zip(
                locus_only.loc[locus_ok, 'gene_id'], locus_m.loc[locus_ok, 'chr'],
                locus_m.loc[locus_ok, 'start'], locus_m.loc[locus_ok, 'end'])
        },
        'n_locus_unplaced': int((~locus_ok).sum()),
        'gene_ids': set(df['gene_id'].dropna()),
        'gene_symbols_with_coords': set(df.loc[has_coords, 'gene_symbol'].dropna()),
        'n_no_coords': int((~has_coords).sum()),
        'sequences': set(df['sequence'].dropna()),
        'alt_ids': set(alt['alt_protein_id']),
        'alt_seqs': dict(zip(alt['alt_protein_id'], alt['sequence'])),
        'ribo_symbols': set(df.loc[df['has_RiboSAM'] & ~df['DDA_evidence'], 'gene_symbol'].dropna()),
    }


# ============================================================
# Checks
# ============================================================

def check_main_files(exp):
    print("\nMain release files")

    fasta = read_fasta(os.path.join(FILES_DIR, 'Unreviewed_Brain_Microproteins.fasta'))
    can_f, alt_f = split_alt(list(fasta), key=str)
    check("FASTA canonical records == gold-standard Salk rows",
          len(can_f) == exp['n_rows'], f"{len(can_f)} vs {exp['n_rows']} expected")
    check("FASTA canonical ids == master gene_id set",
          set(can_f) == exp['gene_ids'],
          f"{len(set(can_f) ^ exp['gene_ids'])} symmetric difference")
    check("FASTA alternative records == expected alt proteoforms",
          set(alt_f) == exp['alt_ids'], f"{len(alt_f)} vs {len(exp['alt_ids'])} expected")
    check("FASTA record ids unique", len(fasta) == len(can_f) + len(alt_f))
    check("FASTA alternative sequences match alt_proteoform_table",
          all(fasta[i][1] == exp['alt_seqs'][i] for i in alt_f))
    check("FASTA canonical headers still 4-field (backward compatible)",
          all(fasta[i][0].count('|') == 3 for i in can_f))
    check("FASTA alternative headers tagged ALT_INIT",
          all('|ALT_INIT|' in fasta[i][0] for i in alt_f))

    ids = read_ids(os.path.join(FILES_DIR, 'Unreviewed_Brain_Microproteins_IDs.txt'))
    can_i, alt_i = split_alt(ids)
    check("IDs canonical == master gene_id set", set(can_i) == exp['gene_ids'],
          f"{len(can_i)} entries")
    check("IDs alternative == expected alt proteoforms", set(alt_i) == exp['alt_ids'])
    check("IDs file has no duplicates", len(ids) == len(set(ids)))

    bed = read_bed12(os.path.join(FILES_DIR,
                                  'Unreviewed_Brain_Microproteins_CDS_Absent_from_UniProt.bed'))
    can_b, alt_b = split_alt(bed, key=lambda r: r['name'])
    check("BED annotated-start == gene_symbols carrying coordinates",
          {r['name'] for r in can_b} == exp['gene_symbols_with_coords'],
          f"{len(can_b)} entries")
    check("BED alternative-start == expected alt proteoforms",
          {r['name'] for r in alt_b} == exp['alt_ids'])
    check("BED intervals well-formed (start < end, start >= 0)",
          all(0 <= r['start'] < r['end'] for r in bed))
    info("spliced records in main BED",
         f"{sum(1 for r in bed if r['block_count'] > 1)} of {len(bed)} span >1 CDS block")
    info("microproteins with no BED/GTF locus",
         f"{exp['n_no_coords']} of {exp['n_rows']} gold-standard rows have a plain "
         f"gene name as gene_symbol, so no coordinates can be derived "
         f"(pre-existing; they appear in FASTA/IDs only)")

    return fasta, bed


def check_all_non_swissprot(exp, main_bed):
    """The all-non-Swiss-Prot BED = the CDS BED plus TrEMBL gene-level loci."""
    print("\nAll-non-Swiss-Prot BED (CDS records + TrEMBL loci)")
    path = os.path.join(FILES_DIR, 'Unreviewed_Brain_Microproteins_All_Non_SwissProt.bed')
    if not os.path.exists(path):
        info("skipped", "file absent — rerun Create_BED_GTF_FASTA_files.py to build it")
        return

    bed = read_bed12(path)
    names = [r['name'] for r in bed]
    expected = (exp['gene_symbols_with_coords'] | exp['alt_ids']
                | set(exp['locus_coords']))
    check("all-non-Swiss-Prot BED names == CDS records + TrEMBL loci",
          set(names) == expected,
          f"{len(bed)} rows, {len(set(names) ^ expected)} symmetric difference")
    check("all-non-Swiss-Prot BED has no duplicate names",
          len(names) == len(set(names)))

    # The CDS half must be the main BED verbatim — same blocks, same strands.
    main_by_name = {r['name']: r for r in main_bed}
    drifted = [r['name'] for r in bed if r['name'] in main_by_name
               and (r['chrom'], r['start'], r['end'], r['strand'],
                    r['sizes'], r['starts']) !=
                   (main_by_name[r['name']]['chrom'], main_by_name[r['name']]['start'],
                    main_by_name[r['name']]['end'], main_by_name[r['name']]['strand'],
                    main_by_name[r['name']]['sizes'], main_by_name[r['name']]['starts'])]
    check("CDS records identical to the main BED",
          not drifted, f"{len(drifted)} drifted" if drifted
          else f"{len(main_by_name)} compared")

    # itemRgb is what separates the three record types, so it is load-bearing:
    # colour is the only marker that a TrEMBL block spans introns.
    with open(path) as f:
        header = f.readline()
        cols = {}
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) >= 9:
                cols[p[3]] = (int(p[6]), int(p[7]), p[8])
    check("track line declares itemRgb=\"On\" (UCSC ignores colour without it)",
          header.startswith('track') and 'itemRgb="On"' in header)

    RGB_CDS, RGB_ALT, RGB_TREMBL = '31,78,121', '130,80,223', '217,119,6'
    wrong_colour = [
        n for n, (_, _, rgb) in cols.items()
        if rgb != (RGB_TREMBL if n in exp['locus_coords']
                   else RGB_ALT if ALT_MARKER in n else RGB_CDS)]
    check("itemRgb matches record type (blue CDS / purple alt / amber TrEMBL)",
          not wrong_colour, f"{len(wrong_colour)} miscoloured" if wrong_colour
          else f"{len(cols)} rows")

    # Locus records: one block spanning exactly the master interval, drawn thick.
    bad_interval, bad_thick = [], []
    for r in bed:
        want = exp['locus_coords'].get(r['name'])
        if want is None:
            continue
        if (r['chrom'], r['start'] + 1, r['end']) != want or r['block_count'] != 1:
            bad_interval.append(r['name'])
        if cols[r['name']][:2] != (r['start'], r['end']):
            bad_thick.append(r['name'])
    check("TrEMBL envelopes match master genomic_coordinates (single block)",
          not bad_interval, f"{len(bad_interval)} wrong"
          if bad_interval else f"{len(exp['locus_coords'])} checked")
    check("TrEMBL envelopes drawn thick across the whole interval",
          not bad_thick, f"{len(bad_thick)} wrong" if bad_thick else "all full-width")

    stranded = sum(1 for r in bed if r['name'] in exp['locus_coords'] and r['strand'] != '.')
    info("TrEMBL loci with a strand",
         f"{stranded} of {len(exp['locus_coords'])} matched an Ensembl gene; "
         f"the rest are '.'")
    if exp['n_locus_unplaced']:
        info("TrEMBL rows with no locus at all",
             f"{exp['n_locus_unplaced']} have no genomic_coordinates in the master "
             f"and appear in FASTA/IDs only")


def check_main_gtf(exp, fasta):
    print("\nMain GTF")
    gtf = read_gtf(os.path.join(FILES_DIR,
                                'Unreviewed_Brain_Microproteins_Absent_from_UniProt.gtf'))
    canonical = {k: v for k, v in gtf.items() if v['source'] == 'GTF2FastaPatched'}
    alt = {k: v for k, v in gtf.items() if v['source'] == 'AltInitiation'}

    check("GTF sources are exactly GTF2FastaPatched + AltInitiation",
          {v['source'] for v in gtf.values()} == {'GTF2FastaPatched', 'AltInitiation'})
    check("GTF canonical transcript_ids subset of coordinate-bearing gene_symbols",
          set(canonical) <= exp['gene_symbols_with_coords'],
          f"{len(canonical)} entries")
    check("GTF alternative transcript_ids == expected alt proteoforms",
          set(alt) == exp['alt_ids'], f"{len(alt)} entries")
    check("no id appears under both sources", not (set(canonical) & set(alt)))

    # Alternative records carry CDS only -- this is what keeps RSEM/salmon/
    # StringTie/featureCounts from ever seeing them as transcripts.
    check("alternative records emit CDS features only (no exon/transcript)",
          all(v['features'] == {'CDS'} for v in alt.values()))
    check("alternative records carry protein_id",
          all('protein_id' in v['attrs'] for v in alt.values()))
    check("alternative protein_id == transcript_id",
          all(v['attrs']['protein_id'] == k for k, v in alt.items()))
    check("alternative gene_id differs from transcript_id (parent linkage kept)",
          all(v['attrs'].get('gene_id') != k for k, v in alt.items()))

    # Coordinates must actually encode the sequence, and phase must be right.
    bad_len, bad_phase = [], []
    for tid, entry in alt.items():
        blocks = sorted((s, e) for s, e, _ in entry['blocks'])
        coding = sum(e - s + 1 for s, e in blocks)
        seq_len = len(fasta[tid][1])
        if coding not in (3 * seq_len, 3 * seq_len + 3):
            bad_len.append(tid)
        ordered = sorted(entry['blocks']) if entry['strand'] == '+' \
            else sorted(entry['blocks'], reverse=True)
        cum = 0
        for s, e, phase in ordered:
            if phase != str((3 - cum % 3) % 3):
                bad_phase.append(tid)
                break
            cum += e - s + 1
    check("alternative CDS length == 3 x protein length (+ optional stop)",
          not bad_len, f"{len(bad_len)} bad" if bad_len else f"{len(alt)} checked")
    check("alternative CDS phase computed correctly",
          not bad_phase, f"{len(bad_phase)} bad" if bad_phase else f"{len(alt)} checked")
    check("alternative CDS blocks disjoint and ascending",
          all(all(a[1] < b[0] for a, b in zip(sorted((s, e) for s, e, _ in v['blocks']),
                                              sorted((s, e) for s, e, _ in v['blocks'])[1:]))
              for v in alt.values()))
    info("spliced alternative proteoforms",
         f"{sum(1 for v in alt.values() if len(v['blocks']) > 1)} of {len(alt)} "
         f"span >1 CDS block")
    return alt


def check_tiers(exp, fasta, alt_gtf):
    print("\nPROSIT confidence tiers")
    tier_ids, all_alt_in_tiers = {}, set()
    for tier in TIERS:
        ids = read_ids(os.path.join(TIER_DIR,
                                    f'Unreviewed_Brain_Microproteins_IDs_{tier}.txt'))
        can, alt = split_alt(ids)
        tier_ids[tier] = (set(can), set(alt))
        all_alt_in_tiers |= set(alt)

        gtf = read_gtf(os.path.join(TIER_DIR, f'Unreviewed_Brain_Microproteins_{tier}.gtf'))
        gtf_alt = {k for k, v in gtf.items() if v['source'] == 'AltInitiation'}
        check(f"[{tier}] GTF alternative ids == tier ID list", gtf_alt == set(alt),
              f"{len(gtf_alt)} in GTF vs {len(alt)} in IDs")

        fa = read_fasta(os.path.join(TIER_DIR, f'Unreviewed_Brain_Microproteins_{tier}.fasta'))
        check(f"[{tier}] FASTA ids == tier ID list", set(fa) == set(ids),
              f"{len(fa)} vs {len(ids)}")
        check(f"[{tier}] every record tagged PROSIT={tier}",
              all(f'PROSIT={tier}' in h for h, _ in fa.values()))

    for a, b in [(x, y) for i, x in enumerate(TIERS) for y in TIERS[i + 1:]]:
        overlap = (tier_ids[a][0] | tier_ids[a][1]) & (tier_ids[b][0] | tier_ids[b][1])
        check(f"tiers {a}/{b} disjoint", not overlap, f"{len(overlap)} shared")

    check("tier alternative records subset of all alternative records",
          all_alt_in_tiers <= exp['alt_ids'])
    ungraded = exp['alt_ids'] - all_alt_in_tiers
    info("ungraded alternative proteoforms",
         f"{len(ungraded)} have no PROSIT grade (supporting peptides lie upstream "
         f"of the alternative start, so they are parent-only evidence) and are "
         f"correctly absent from all tier files")

    # Cross-check tier assignment against the standalone mapping table.
    mapping_path = os.path.join(ALT_DIR, 'Alt_Proteoforms_mapping.tsv')
    if os.path.exists(mapping_path):
        m = pd.read_csv(mapping_path, sep='\t')
        grade = dict(zip(m['alt_protein_id'], m['prosit_grade']))
        mismatched = [t for t in TIERS
                      if tier_ids[t][1] != {k for k, v in grade.items() if v == t}]
        check("tier membership matches derived prosit_grade in mapping TSV",
              not mismatched, f"mismatched tiers: {mismatched}" if mismatched else '')


def check_all_bed_files(alt_gtf):
    """Every .bed must be BED12 with internally consistent blocks."""
    print("\nBED12 block integrity (all .bed files)")
    paths = []
    for root, _, names in os.walk(FILES_DIR):
        paths.extend(os.path.join(root, n) for n in sorted(names) if n.endswith('.bed'))

    bad_cols, bad_blocks = [], []
    for path in sorted(paths):
        rel = os.path.relpath(path, FILES_DIR)
        with open(path) as f:
            first = f.readline().rstrip('\n')
            while first.startswith(('track', 'browser')):   # UCSC track headers
                first = f.readline().rstrip('\n')
        if first and len(first.split('\t')) != 12:
            bad_cols.append(f"{rel} ({len(first.split(chr(9)))} cols)")
            continue
        for r in read_bed12(path):
            if (r['starts'][0] != 0
                    or r['block_count'] != len(r['sizes']) != len(r['starts'])
                    or r['start'] + r['starts'][-1] + r['sizes'][-1] != r['end']
                    or any(b <= 0 for b in r['sizes'])
                    or any(r['starts'][i] + r['sizes'][i] > r['starts'][i + 1]
                           for i in range(len(r['sizes']) - 1))):
                bad_blocks.append(f"{rel}:{r['name']}")

    check("every .bed file is 12-column BED12", not bad_cols,
          '; '.join(bad_cols) if bad_cols else f"{len(paths)} files")
    check("all BED12 blocks internally consistent "
          "(ascending, disjoint, ending at chromEnd)", not bad_blocks,
          f"{len(bad_blocks)} bad" if bad_blocks else "all records")

    # The main BED's blocks must agree with the GTF's CDS blocks.
    main = read_bed12(os.path.join(FILES_DIR,
                                   'Unreviewed_Brain_Microproteins_CDS_Absent_from_UniProt.bed'))
    mismatched = []
    for r in main:
        if r['name'] not in alt_gtf:
            continue
        derived = sorted((r['start'] + st + 1, r['start'] + st + sz)
                         for st, sz in zip(r['starts'], r['sizes']))
        if derived != sorted((s, e) for s, e, _ in alt_gtf[r['name']]['blocks']):
            mismatched.append(r['name'])
    check("alternative-start BED blocks match their GTF CDS blocks",
          not mismatched, f"{len(mismatched)} mismatched" if mismatched
          else f"{sum(1 for r in main if r['name'] in alt_gtf)} checked")


def check_ribo(exp):
    print("\nRibo-ShortStop set")
    ids = read_ids(os.path.join(RIBO_DIR,
                                'Unreviewed_Brain_Microproteins_IDs_Ribo_ShortStop.txt'))
    _, alt = split_alt(ids)
    check("Ribo-ShortStop contains no alternative proteoforms "
          "(all alt parents are DDA-detected)", not alt, f"{len(alt)} found")
    gtf = read_gtf(os.path.join(RIBO_DIR,
                                'Unreviewed_Brain_Microproteins_Ribo_ShortStop.gtf'))
    check("Ribo-ShortStop GTF ids subset of expected has_RiboSAM & ~DDA set",
          set(gtf) <= exp['ribo_symbols'], f"{len(gtf)} entries")


def check_standalone(exp, alt_gtf):
    print("\nStandalone Alt_Proteoforms/ set")
    if not os.path.isdir(ALT_DIR):
        info("Alt_Proteoforms/", "directory absent — skipping")
        return
    fa = read_fasta(os.path.join(ALT_DIR, 'Alt_Proteoforms.fasta'))
    check("standalone FASTA == expected alt proteoforms", set(fa) == exp['alt_ids'],
          f"{len(fa)} records")

    gtf = read_gtf(os.path.join(ALT_DIR, 'Alt_Proteoforms.gtf'))
    same_blocks = all(
        sorted((s, e) for s, e, _ in gtf[k]['blocks'])
        == sorted((s, e) for s, e, _ in alt_gtf[k]['blocks'])
        for k in gtf if k in alt_gtf)
    check("standalone GTF blocks identical to merged GTF blocks", same_blocks,
          "shared builder has not drifted")

    bed12 = read_bed12(os.path.join(ALT_DIR, 'Alt_Proteoforms_CDS.bed'))
    check("standalone BED12 record set == expected alt proteoforms",
          {r['name'] for r in bed12} == exp['alt_ids'], f"{len(bed12)} records")
    bad = []
    for r in bed12:
        if r['starts'][0] != 0 or r['block_count'] != len(r['sizes']):
            bad.append(r['name'])
        elif r['start'] + r['starts'][-1] + r['sizes'][-1] != r['end']:
            bad.append(r['name'])
        else:
            derived = sorted((r['start'] + st + 1, r['start'] + st + sz)
                             for st, sz in zip(r['starts'], r['sizes']))
            if derived != sorted((s, e) for s, e, _ in alt_gtf[r['name']]['blocks']):
                bad.append(r['name'])
    check("BED12 blocks internally consistent and match GTF", not bad,
          f"{len(bad)} bad" if bad else f"{len(bed12)} checked")


def check_against_source_gtf(exp, alt_gtf, source_gtf):
    """Optional: re-derive alternative coordinates from parent CDS blocks."""
    print("\nDeep check against combined source GTF")
    sys.path.insert(0, os.path.join(CODE_DIR, 'Microprotein_annotation_summary'))
    import alt_proteoform_records as apr

    alt_table = pd.read_csv(os.path.join(CODE_DIR, 'data', 'alt_proteoform_table.csv'))
    alt_table = alt_table[alt_table['alt_protein_id'].isin(exp['alt_ids'])]
    index = apr.index_parent_cds(source_gtf, set(alt_table['orf_id']))
    bad = []
    for _, row in alt_table.iterrows():
        parent = index.get(row['orf_id'])
        if parent is None:
            bad.append(row['alt_protein_id'])
            continue
        blocks, _ = apr.trim_5prime(parent['blocks'], parent['strand'],
                                    3 * (int(row['aa_position']) - 1))
        if blocks != sorted((s, e) for s, e, _ in alt_gtf[row['alt_protein_id']]['blocks']):
            bad.append(row['alt_protein_id'])
    check("released coordinates reproduce trim of parent CDS blocks", not bad,
          f"{len(bad)} mismatched" if bad else f"{len(alt_table)} checked")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source-gtf', default=os.environ.get('COMBINED_SOURCE_GTF'),
                    help='optional combined source GTF for the deep coordinate check')
    args = ap.parse_args()

    print("=" * 70)
    print("Verifying GTF_and_BED_files/")
    print("=" * 70)

    exp = build_expectations()
    print(f"\nExpected from master + alt table:")
    print(f"  gold-standard Salk rows (== unique gene_id, == unique sequence): {exp['n_rows']}")
    print(f"  of those, carrying derivable coordinates                      : "
          f"{len(exp['gene_symbols_with_coords'])}")
    print(f"  alternative proteoforms                                       : "
          f"{len(exp['alt_ids'])}")

    fasta, main_bed = check_main_files(exp)
    check_all_non_swissprot(exp, main_bed)
    alt_gtf = check_main_gtf(exp, fasta)
    check_all_bed_files(alt_gtf)
    check_tiers(exp, fasta, alt_gtf)
    check_ribo(exp)
    check_standalone(exp, alt_gtf)
    if args.source_gtf and os.path.exists(args.source_gtf):
        check_against_source_gtf(exp, alt_gtf, args.source_gtf)
    else:
        print("\nDeep check against combined source GTF")
        info("skipped", "pass --source-gtf or set COMBINED_SOURCE_GTF to enable")

    failed = [n for n, ok, _ in results if not ok]
    print("\n" + "=" * 70)
    print(f"{len(results) - len(failed)}/{len(results)} checks passed")
    if failed:
        print("\nFAILED:")
        for n in failed:
            print(f"  - {n}")
    print("=" * 70)
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
