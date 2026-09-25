"""Attach per-ORF P-site frame statistics to the gold-standard microprotein set.

Input is Code/data/psite_frame_counts_per_orf.tsv (one row per CDS in the combined
GENCODE v43 + brain microprotein GTF; see psite_frame_counts_per_orf.py).

Matching:
  - Unreviewed (Salk smORFs): by amino-acid sequence. The GTF smORF records
    carry their translated sequence, so this is exact.
  - TrEMBL: sequence -> genomic coordinates from
    GTF_and_BED_files/Unreviewed_Brain_Microproteins_mapping_coordinates_to_sequences.tsv,
    then to the GENCODE CDS with that span (with or without the stop codon)
    and n_codons == protein_length.
  - Swiss-Prot-MP: by gene_name, restricted to the gene's CDSs with
    n_codons == protein_length.
  A gene can have several CDSs that pass. They are used only if their P-site
  values agree; if not, the CDS whose span matches the entry's coordinates is
  preferred, and anything still ambiguous is left unmatched.

Output: ../../Results/RP3/Psite_frame_by_sequence.csv, keyed by sequence, read
by RP3_Results_summary.py and the dashboard.

usage: python Psite_frame_mapping.py
"""
import os
import sys

import pandas as pd

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, '..'))
from gold_standard_filtering_criteria import load_and_filter_master

psite_file = os.path.join(here, '..', 'data', 'psite_frame_counts_per_orf.tsv')
master_file = os.path.join(here, '..', 'data', 'microprotein_master.csv')
coord_file = os.path.join(here, '..', '..', 'GTF_and_BED_files',
                          'Unreviewed_Brain_Microproteins_mapping_coordinates_to_sequences.tsv')
out_file = os.path.join(here, '..', '..', 'Results', 'RP3', 'Psite_frame_by_sequence.csv')

value_cols = ['Psites_pct_frame0', 'Psites_pct_frame1', 'Psites_pct_frame2', 'Psites_frame0_RPKM']


def add_spans(psites):
    """CDS span as chr:start-end, and the same span extended over the stop codon."""
    plus = psites['strand'] == '+'
    start, end = psites['cds_start'], psites['cds_end']
    psites['span'] = psites['chrom'] + ':' + start.astype(str) + '-' + end.astype(str)
    psites['span_stop'] = (psites['chrom'] + ':' + (start - 3 * ~plus).astype(str)
                           + '-' + (end + 3 * plus).astype(str))
    return psites


def resolve(cands, method):
    """One row per sequence from candidate CDSs (see module docstring)."""
    cands = cands.copy()
    cands['span_match'] = ((cands['genomic_coordinates'] == cands['span'])
                           | (cands['genomic_coordinates'] == cands['span_stop']))

    def n_distinct(d):
        return d.groupby('sequence')[value_cols].nunique(dropna=False).max(axis=1)

    n = n_distinct(cands)
    ambiguous = n[n > 1].index
    by_span = cands[cands['sequence'].isin(ambiguous) & cands['span_match']]
    n_span = n_distinct(by_span)
    by_span = by_span[by_span['sequence'].isin(n_span[n_span <= 1].index)]
    kept = pd.concat([cands[~cands['sequence'].isin(ambiguous)], by_span])

    m = (kept.groupby('sequence')
             .agg({'transcript_id': lambda s: ';'.join(sorted(set(s))),
                   **{c: 'first' for c in value_cols}})
             .reset_index())
    m['Psite_match'] = method
    return m, len(ambiguous) - by_span['sequence'].nunique()


def match_salk(mp, smorfs):
    salk = mp[(mp['Database'] == 'Salk') & (mp['smorf_type'] != 'TrEMBL')]
    m = salk[['sequence']].merge(smorfs[['sequence', 'transcript_id'] + value_cols], on='sequence')
    m['Psite_match'] = 'sequence'
    return m


def match_trembl(mp, ensembl):
    coords = pd.read_csv(coord_file, sep='\t')
    tr = mp.loc[mp['smorf_type'] == 'TrEMBL', ['sequence', 'protein_length']].merge(coords, on='sequence')
    spans = pd.concat([ensembl.assign(key=ensembl['span']),
                       ensembl.assign(key=ensembl['span_stop'])]).drop_duplicates(['source', 'transcript_id', 'key'])
    cands = tr.merge(spans, left_on='genomic_coordinates', right_on='key')
    cands = cands[cands['n_codons'] == cands['protein_length']]
    return resolve(cands, 'coordinates')


def match_swissprot(mp, ensembl):
    sp = mp.loc[mp['Database'] == 'Swiss-Prot-MP',
                ['sequence', 'gene_name', 'protein_length', 'genomic_coordinates']]
    cands = sp.merge(ensembl, on='gene_name')
    cands = cands[cands['n_codons'] == cands['protein_length']]
    return resolve(cands, 'gene_name')


def main():
    mp = load_and_filter_master(master_file)
    psites = add_spans(pd.read_csv(psite_file, sep='\t', low_memory=False))
    smorfs = psites[psites['is_smorf'] & psites['sequence'].notna()]
    ensembl = psites[~psites['is_smorf']].drop(columns=['sequence', 'smorf_type'])

    by_trembl, amb_trembl = match_trembl(mp, ensembl)
    by_sp, amb_sp = match_swissprot(mp, ensembl)
    out = pd.concat([match_salk(mp, smorfs), by_trembl, by_sp], ignore_index=True)
    out = out.rename(columns={'transcript_id': 'Psite_transcript_id'})
    out = out[['sequence'] + value_cols + ['Psite_match', 'Psite_transcript_id']]
    assert not out['sequence'].duplicated().any()

    groups = mp['smorf_type'].where(mp['smorf_type'] == 'TrEMBL', mp['Database'])
    matched = mp['sequence'].isin(out['sequence'])
    print('Matched / total by group:')
    for g in groups.unique():
        print(f"  {g:14s} {matched[groups == g].sum():5d} / {(groups == g).sum()}")
    print(f"  left unmatched as ambiguous: TrEMBL {amb_trembl}, Swiss-Prot-MP {amb_sp}")

    os.makedirs(os.path.dirname(out_file), exist_ok=True)
    out.to_csv(out_file, index=False, float_format='%.4g')
    print('wrote', os.path.normpath(out_file), len(out), 'rows')


if __name__ == '__main__':
    main()
