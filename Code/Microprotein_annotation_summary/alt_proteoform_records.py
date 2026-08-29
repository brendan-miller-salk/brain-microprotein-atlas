"""
Shared builder for alternative proteoform (non-ATG alternative initiation) records.

Imported by both Create_BED_GTF_FASTA_files.py (which merges these records into
the main release file set) and Create_Alt_Proteoform_Files.py (which writes the
standalone Alt_Proteoforms/ set). Pure functions only -- importing this module
has no side effects.

Two things here are easy to get wrong and are handled explicitly:

Coordinates. `CDS_start`/`CDS_end` in alt_proteoform_table.csv are the genomic
ENVELOPE, not a contiguous CDS -- 27% of alternative proteoforms are spliced
(worst case a 73-aa protein spanning 109,663 bp). Real coordinates come from
trimming the parent's CDS blocks in TRANSCRIPT space; the CSV envelope is used
only as a cross-check.

PROSIT grade. An alternative proteoform has no `Confidence` of its own in the
master -- that column lives on the parent. The grade is therefore derived from
the parent's tryptic peptides that fall entirely DOWNSTREAM of the alternative
initiation site, since only those peptides are evidence for the shorter
proteoform. Peptides lying upstream are parent-only evidence and must not
confer a grade (34 gold-standard proteoforms are correctly left ungraded on
this basis). Note this is a different axis from the table's own `tier` column,
which scores initiation-site confidence, not spectral-angle match quality.
"""

import ast
import os

import pandas as pd

# PROSIT spectral-angle tiers, best first.
CONFIDENCE_RANK = {'Strong': 0, 'Moderate': 1, 'Weak': 2, 'Insufficient': 3}

GTF_SOURCE = 'AltInitiation'   # distinct from GTF2FastaPatched so alternative
                               # records stay separable with a single awk


def _parse_list(value):
    """Parse a bracketed-list cell (e.g. "['Strong', 'Moderate']") to a list."""
    if pd.isna(value):
        return []
    try:
        parsed = ast.literal_eval(value) if isinstance(value, str) else value
    except (ValueError, SyntaxError):
        return []
    return parsed if isinstance(parsed, list) else []


def best_confidence(values):
    """Collapse PROSIT tiers to the single best (highest-confidence) one."""
    valid = [v for v in values if v in CONFIDENCE_RANK]
    return min(valid, key=lambda v: CONFIDENCE_RANK[v]) if valid else None


def trim_5prime(blocks, strand, n_nt):
    """Remove n_nt from the 5' end of a CDS block list, in transcript space.

    Walks blocks left-to-right on '+' and right-to-left on '-', dropping blocks
    fully consumed and partially trimming the block that straddles the offset.
    Returns (blocks sorted ascending by coordinate, bases actually consumed).
    """
    ordered = sorted(blocks) if strand == '+' else sorted(blocks, reverse=True)
    out, remaining = [], n_nt
    for start, end in ordered:
        block_len = end - start + 1
        if remaining >= block_len:
            remaining -= block_len
            continue
        out.append((start + remaining, end) if strand == '+' else (start, end - remaining))
        remaining = 0
    return sorted(out), n_nt - remaining


def cds_phases(blocks, strand):
    """GTF phase per block, keyed by (start, end).

    Phase is the bases to remove from the start of the feature to reach the
    first base of the next codon -- (3 - cumulative % 3) % 3, walking in
    TRANSLATION order. The 5'-most block is always 0 because translation begins
    exactly at the alternative initiation codon.
    """
    ordered = sorted(blocks) if strand == '+' else sorted(blocks, reverse=True)
    phases, cumulative = {}, 0
    for start, end in ordered:
        phases[(start, end)] = (3 - cumulative % 3) % 3
        cumulative += end - start + 1
    return phases


def load_peptide_prosit(code_dir):
    """Map protein_id -> list of {start, end, conf} for every tryptic peptide."""
    path = os.path.join(code_dir, 'data',
                        'cleaned_tryptic_peptides_detailed_under_151aa_with_SA.csv')
    peptides = {}
    table = pd.read_csv(path)
    for _, row in table.iterrows():
        seqs = _parse_list(row['peptide_sequence'])
        if not seqs:
            continue
        starts, ends = _parse_list(row['start']), _parse_list(row['end'])
        confs = _parse_list(row['Confidence'])
        peptides.setdefault(row['protein_id'], []).extend(
            {'start': starts[i], 'end': ends[i],
             'conf': confs[i] if i < len(confs) else None}
            for i in range(len(seqs))
        )
    return peptides


def index_parent_cds(source_gtf, wanted_ids):
    """Stream the combined source GTF once, indexing GTF2FastaPatched CDS blocks."""
    index = {}
    with open(source_gtf, 'r') as src:
        for line in src:
            if line.startswith('#'):
                continue
            fields = line.split('\t')
            # Cheapest filters first -- this is a single pass over ~1.1 GB.
            if len(fields) < 9 or fields[2] != 'CDS' or fields[1] != 'GTF2FastaPatched':
                continue
            attrs = fields[8]
            tid_start = attrs.find('transcript_id "')
            if tid_start == -1:
                continue
            tid_start += len('transcript_id "')
            transcript_id = attrs[tid_start:attrs.find('"', tid_start)]
            if transcript_id not in wanted_ids:
                continue
            entry = index.setdefault(
                transcript_id, {'chrom': fields[0], 'strand': fields[6], 'blocks': []}
            )
            entry['blocks'].append((int(fields[3]), int(fields[4])))
    return index


def build_alt_records(code_dir, source_gtf, gold_standard_ids=None,
                      parent_meta=None, verbose=True, cds_index=None):
    """Build validated alternative proteoform records.

    gold_standard_ids : restrict to alt proteoforms whose parent survives the
                        canonical filter. None keeps all.
    parent_meta       : DataFrame indexed by parent gene_symbol supplying
                        'gene_id' and 'Annotation Status'. Optional.
    cds_index         : prebuilt index from index_parent_cds, to avoid a second
                        streaming pass over the ~1.1 GB source GTF when the
                        caller already needs one. Built here if omitted.

    Returns (records, failures). Records failing validation are excluded and
    reported rather than raising, matching the existing exporters' behaviour.
    """
    alt = pd.read_csv(os.path.join(code_dir, 'data', 'alt_proteoform_table.csv'))
    alt = alt[alt['alt_protein'] == 'Yes'].copy()

    # The 'length' column disagrees with len(sequence) on the 337 M-substituted
    # rows, so it is dropped outright rather than trusted.
    alt = alt.drop(columns=['length'])
    alt['aa_length'] = alt['sequence'].str.len()

    assert alt['alt_protein_id'].is_unique, "alt_protein_id is not unique"
    assert alt['sequence'].notna().all(), "missing sequence in alt table"
    assert (alt['aa_position'] >= 2).all(), "aa_position < 2 implies no trim"

    # Where a parent has two rows they are the same initiation site written
    # twice: once with the native residue and once M-substituted. Deliberate
    # sequence variants, not duplicates -- keep both, separated by protein_id.
    paired = alt.groupby(['orf_id', 'aa_position'])['alt_protein_id'].transform('size') > 1
    alt['initiation_variant'] = 'single'
    alt.loc[paired & (alt['sequence'].str[0] == 'M'), 'initiation_variant'] = 'met_substituted'
    alt.loc[paired & (alt['sequence'].str[0] != 'M'), 'initiation_variant'] = 'native'

    dropped = pd.DataFrame()
    if gold_standard_ids is not None:
        in_gs = alt['orf_id'].isin(gold_standard_ids)
        dropped = alt[~in_gs]
        alt = alt[in_gs].copy()

    peptides = load_peptide_prosit(code_dir)
    if cds_index is None:
        cds_index = index_parent_cds(source_gtf, set(alt['orf_id'].unique()))
    if verbose:
        wanted = set(alt['orf_id'].unique())
        print(f"  alternative-start parents resolved: "
              f"{len(wanted & set(cds_index))}/{len(wanted)}")

    records, failures = [], []
    for _, row in alt.iterrows():
        parent = cds_index.get(row['orf_id'])
        if parent is None:
            failures.append((row['alt_protein_id'], 'parent_missing_from_gtf', '', ''))
            continue

        trim_nt = 3 * (int(row['aa_position']) - 1)
        blocks, consumed = trim_5prime(parent['blocks'], parent['strand'], trim_nt)

        checks = []
        if consumed != trim_nt:
            checks.append(('trim_underconsumed', trim_nt, consumed))
        if not blocks:
            checks.append(('empty_after_trim', '>0 blocks', 0))
        else:
            coding_len = sum(e - s + 1 for s, e in blocks)
            expected = 3 * int(row['aa_length'])
            if coding_len not in (expected, expected + 3):
                checks.append(('coding_length', f'{expected} or {expected + 3}', coding_len))
            if pd.notna(row['CDS_start']):
                csv_envelope = (int(row['CDS_start']), int(row['CDS_end']))
                if (blocks[0][0], blocks[-1][1]) != csv_envelope:
                    checks.append(('envelope_vs_csv', csv_envelope,
                                   (blocks[0][0], blocks[-1][1])))
            if parent['chrom'] != row['chromosome'] or parent['strand'] != row['strand']:
                checks.append(('chrom_strand', f"{row['chromosome']}{row['strand']}",
                               f"{parent['chrom']}{parent['strand']}"))
            prev_end = None
            for s, e in blocks:
                if prev_end is not None and s <= prev_end:
                    checks.append(('blocks_overlap', 'ascending disjoint', str(blocks)))
                    break
                prev_end = e

        if checks:
            failures.extend((row['alt_protein_id'], n, e, o) for n, e, o in checks)
            continue

        # PROSIT grade from peptides lying entirely downstream of the new start.
        start_aa = int(row['aa_position'])
        downstream = [p for p in peptides.get(row['orf_id'], []) if p['start'] >= start_aa]
        grade = best_confidence([p['conf'] for p in downstream])

        coding_len = sum(e - s + 1 for s, e in blocks)
        record = {
            'alt_protein_id': row['alt_protein_id'],
            'parent_orf_id': row['orf_id'],
            'gene_name': row['gene_name'],
            'chrom': parent['chrom'],
            'strand': parent['strand'],
            'blocks': blocks,
            'sequence': row['sequence'],
            'aa_length': int(row['aa_length']),
            'aa_position': start_aa,
            'codon': row['codon'],
            'initiation_tier': row['tier'],
            'cognate_status': row['cognate_status'],
            'context_strength': row['context_strength'],
            'n_peptides': row['n_peptides'],
            'n_support_peptides': len(downstream),
            'smorf_type': row['smorf_type'],
            'initiation_variant': row['initiation_variant'],
            'includes_stop': coding_len == 3 * int(row['aa_length']) + 3,
            'prosit_grade': grade,
            'Database': 'Salk',
            'genomic_coordinates': f"{parent['chrom']}:{blocks[0][0]}-{blocks[-1][1]}",
        }
        if parent_meta is not None and row['orf_id'] in parent_meta.index:
            record['parent_gene_id'] = parent_meta.at[row['orf_id'], 'gene_id']
            record['annotation_status'] = parent_meta.at[row['orf_id'], 'Annotation Status']
        else:
            record['parent_gene_id'] = row['orf_id']
            record['annotation_status'] = 'MS'
        records.append(record)

    return records, failures, dropped


# ============================================================
# Rendering helpers -- shared so both exporters emit identical formats
# ============================================================

def fasta_header(record):
    """>{id}|{gene_name}|{Database}|{Annotation Status}|ALT_INIT|codon=..|tier=..

    Extends the canonical 4-field header in field order, so any consumer doing
    split('|')[:4] keeps working.
    """
    return (f">{record['alt_protein_id']}|{record['gene_name']}|{record['Database']}|"
            f"{record['annotation_status']}|ALT_INIT|codon={record['codon']}|"
            f"init_tier={record['initiation_tier']}")


def gtf_lines(record):
    """Yield GTF CDS lines for one record (CDS features only -- see module docstring)."""
    phases = cds_phases(record['blocks'], record['strand'])
    attrs = ' '.join(f'{k} "{v}";' for k, v in [
        ('gene_id', record['parent_gene_id']),
        ('transcript_id', record['alt_protein_id']),
        ('protein_id', record['alt_protein_id']),
        ('gene_name', record['gene_name']),
        ('parent_orf_id', record['parent_orf_id']),
        ('start_codon_seq', record['codon']),
        ('aa_position', record['aa_position']),
        ('initiation_tier', record['initiation_tier']),
        ('cognate_status', record['cognate_status']),
        ('context_strength', record['context_strength']),
        ('initiation_variant', record['initiation_variant']),
        ('includes_stop', str(record['includes_stop']).lower()),
        ('prosit_grade', record['prosit_grade'] if record['prosit_grade'] else 'NA'),
        ('smorf_type', record['smorf_type']),
    ])
    for start, end in record['blocks']:
        yield (f"{record['chrom']}\t{GTF_SOURCE}\tCDS\t{start}\t{end}\t1000\t"
               f"{record['strand']}\t{phases[(start, end)]}\t{attrs}\n")


def bed6_line(record):
    """BED6 envelope, matching the existing release files' shape."""
    return (f"{record['chrom']}\t{record['blocks'][0][0] - 1}\t{record['blocks'][-1][1]}\t"
            f"{record['alt_protein_id']}\t0\t{record['strand']}\n")


def bed12_line_from_blocks(chrom, blocks, name, strand):
    """BED12 from 1-based inclusive GTF CDS blocks.

    Blocks always ascend by genomic coordinate regardless of strand, per the BED
    spec. thickStart/thickEnd span the whole feature because every block here is
    coding. A single-block feature degrades to the same interval a BED6 would
    carry, so this is a strict superset of the old format.
    """
    blocks = sorted(blocks)
    chrom_start = blocks[0][0] - 1               # GTF 1-based -> BED 0-based
    chrom_end = blocks[-1][1]
    sizes = [e - s + 1 for s, e in blocks]
    starts = [(s - 1) - chrom_start for s, _ in blocks]
    assert starts[0] == 0
    assert chrom_start + starts[-1] + sizes[-1] == chrom_end
    return '\t'.join([
        chrom, str(chrom_start), str(chrom_end), name, '0', strand,
        str(chrom_start), str(chrom_end), '0', str(len(blocks)),
        ','.join(map(str, sizes)) + ',',
        ','.join(map(str, starts)) + ',',
    ]) + '\n'


def bed12_line(record):
    """BED12 for one alternative proteoform record."""
    return bed12_line_from_blocks(record['chrom'], record['blocks'],
                                  record['alt_protein_id'], record['strand'])
