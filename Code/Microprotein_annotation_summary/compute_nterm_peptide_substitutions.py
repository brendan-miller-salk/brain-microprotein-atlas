"""
Detect amino-acid substitutions in N-terminal tryptic peptides vs. the
matched UniProt isoform, for the dashboard's N-terminus filter.

The dashboard's N-terminus filter (`NTERM_MODE_NON_NTERM` /
`NTERM_MODE_ACETYL_OR_NON_NTERM` in Results/microproteins_dashboard.py)
discards a microprotein whose only tryptic evidence is a peptide at aa 1 or
2, on the theory that it's indistinguishable from Met-excision of the parent
ORF -- unless the microprotein is Nt-acetylated.

But if that exact peptide's own sequence carries 1-2 amino-acid
substitutions relative to the matched canonical/isoform protein at the same
position, it CANNOT be a plain Met-excision fragment of that canonical
protein -- it's a materially different sequence, and its own N-terminal
residue was captured by MS with a real variant. That's peptide-level
evidence, independent of where in the isoform the alignment starts (a
population dominated by legitimately downstream-start "Short-Isoform"
smORFs, where BLAST naturally reports a large offset -- see the earlier,
now-abandoned position-offset approach this script replaces).

This script BLASTs each N-terminal peptide itself (not the whole
microprotein) against its matched isoform's full-length sequence, and
records which residues differ.

Outputs (both in Code/data/, both git-shipped once generated):
  - blast_isoform_reference_sequences.fasta: cached full-length UniProt
    sequences for every matched accession seen so far (shared with/reused
    from the earlier position-offset script; only re-fetches accessions not
    already cached).
  - blast_nterm_peptide_substitutions.csv: one row per (microprotein,
    N-terminal peptide) pair, with the resolved substitutions (if any).

Network access to https://rest.uniprot.org is required only when the FASTA
cache is missing an accession; a full local cache makes subsequent runs
network-free, consistent with the shipped-data-only philosophy documented in
the repo root CLAUDE.md.
"""

import ast
import os
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.request

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from gold_standard_filtering_criteria import load_and_filter_master

CODE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
DATA_DIR = os.path.join(CODE_DIR, 'data')
MASTER_CSV = os.path.join(DATA_DIR, 'microprotein_master.csv')
REFERENCE_FASTA = os.path.join(DATA_DIR, 'blast_isoform_reference_sequences.fasta')
SUBSTITUTIONS_CSV = os.path.join(DATA_DIR, 'blast_nterm_peptide_substitutions.csv')

BLASTP_BIN = '/usr/local/ncbi-blast/bin/blastp'

# A partial/fragmentary hit says nothing about substitutions -- require the
# alignment to cover most of the peptide before trusting its mismatch count.
MIN_COVERAGE = 0.8

UNIPROT_FASTA_URL = "https://rest.uniprot.org/uniprotkb/{accession}.fasta"


def _has_only_nterm_peptides(pos_str):
    """True if every tryptic peptide starts at aa 1 or 2 (the artifact zone).

    Mirrors `_has_non_nterm_peptide` in Results/microproteins_dashboard.py --
    that function asks "is there a peptide at aa>=3?"; this asks the
    complementary question restricted to rows that actually have peptide
    data, i.e. exactly the population this script needs to re-examine.
    """
    if not pos_str or pd.isna(pos_str):
        return False
    try:
        positions = ast.literal_eval(str(pos_str))
    except (ValueError, SyntaxError):
        return False
    if isinstance(positions, (int, float)):
        positions = [positions]
    try:
        positions = [int(p) for p in positions]
    except (TypeError, ValueError):
        return False
    if not positions:
        return False
    return all(p <= 2 for p in positions)


def _as_list(val):
    if val is None or (isinstance(val, float) and pd.isna(val)):
        return []
    try:
        parsed = ast.literal_eval(str(val))
    except (ValueError, SyntaxError):
        return []
    return parsed if isinstance(parsed, list) else [parsed]


def _select_candidates():
    """(microprotein, N-terminal peptide) pairs to examine: every peptide
    starting at aa <=2 belonging to a microprotein whose tryptic evidence is
    entirely confined to aa 1-2, with a BLAST match to re-check against."""
    mp = load_and_filter_master(MASTER_CSV)
    only_nterm = mp['start'].map(_has_only_nterm_peptides)
    has_match = mp['blastp_uniprot_accession_match'].notna()
    cand_proteins = mp.loc[only_nterm & has_match,
                           ['sequence', 'peptide_sequence', 'start',
                            'blastp_uniprot_accession_match']].drop_duplicates(subset='sequence')

    rows = []
    for row in cand_proteins.itertuples(index=False):
        peptides = _as_list(row.peptide_sequence)
        starts = _as_list(row.start)
        for pep, st in zip(peptides, starts):
            try:
                if int(st) <= 2:
                    rows.append({
                        'sequence': row.sequence,
                        'peptide_sequence': pep,
                        'peptide_start': int(st),
                        'uniprot_accession': row.blastp_uniprot_accession_match,
                    })
            except (TypeError, ValueError):
                continue
    return pd.DataFrame(rows)


def _read_cached_fasta(path):
    """Parse a FASTA file into {accession: sequence}, keyed by the first
    pipe-or-whitespace-delimited token after '>' (works for both bare
    accessions and UniProt's 'sp|ACCESSION|NAME' headers)."""
    cache = {}
    if not os.path.exists(path):
        return cache
    header, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if header is not None:
                    cache[header] = ''.join(chunks)
                raw = line[1:].split()[0]
                parts = raw.split('|')
                header = parts[1] if len(parts) >= 2 else parts[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if header is not None:
        cache[header] = ''.join(chunks)
    return cache


def _fetch_uniprot_fasta(accession):
    """Fetch one accession's FASTA record from UniProt's REST API.

    Isoform accessions (e.g. 'P12345-2') are fetched individually via the
    single-entry endpoint rather than batched, because the isoform suffix
    denotes a *specific* sequence variant that a bulk accession-list request
    is not guaranteed to resolve correctly -- this only runs at most once per
    unique accession, ever, since results are cached to REFERENCE_FASTA.
    """
    url = UNIPROT_FASTA_URL.format(accession=accession)
    try:
        with urllib.request.urlopen(url, timeout=30) as resp:
            text = resp.read().decode('utf-8')
    except urllib.error.HTTPError as e:
        print(f"  [WARN] {accession}: HTTP {e.code}")
        return None
    except urllib.error.URLError as e:
        print(f"  [WARN] {accession}: {e}")
        return None
    lines = text.splitlines()
    if not lines or not lines[0].startswith('>'):
        return None
    seq = ''.join(l.strip() for l in lines[1:])
    return seq or None


def _ensure_reference_sequences(accessions):
    """Make sure REFERENCE_FASTA has an entry for every accession, fetching
    only what's missing. Returns {accession: sequence}."""
    cache = _read_cached_fasta(REFERENCE_FASTA)
    missing = sorted(set(accessions) - set(cache))
    if missing:
        print(f"Fetching {len(missing)} new isoform sequence(s) from UniProt "
              f"({len(accessions) - len(missing)} already cached)...")
        with open(REFERENCE_FASTA, 'a') as fh:
            for i, acc in enumerate(missing, 1):
                seq = _fetch_uniprot_fasta(acc)
                if seq:
                    fh.write(f">{acc}\n{seq}\n")
                    cache[acc] = seq
                if i % 25 == 0 or i == len(missing):
                    print(f"  {i}/{len(missing)} fetched")
                time.sleep(0.1)
    else:
        print(f"All {len(accessions)} isoform sequence(s) already cached in "
              f"{os.path.basename(REFERENCE_FASTA)}.")
    return cache


def _align_peptide(peptide_seq, subject_seq):
    """BLAST a single short peptide against one subject sequence.

    Returns the first HSP's (qseq, sseq, qstart, qend, pident, mismatch,
    length), or None if there's no hit. Uses -task blastp-short, NCBI's
    recommended mode for queries under ~30 residues.
    """
    with tempfile.NamedTemporaryFile('w', suffix='.fasta', delete=False) as qf, \
         tempfile.NamedTemporaryFile('w', suffix='.fasta', delete=False) as sf:
        qf.write(f">query\n{peptide_seq}\n")
        sf.write(f">subject\n{subject_seq}\n")
        qpath, spath = qf.name, sf.name
    try:
        result = subprocess.run(
            [BLASTP_BIN, '-task', 'blastp-short',
             '-query', qpath, '-subject', spath,
             '-outfmt', '6 qseq sseq qstart qend pident mismatch length'],
            capture_output=True, text=True, timeout=30,
        )
    finally:
        os.unlink(qpath)
        os.unlink(spath)
    if result.returncode != 0 or not result.stdout.strip():
        return None
    first_line = result.stdout.strip().splitlines()[0]
    qseq, sseq, qstart, qend, pident, mismatch, length = first_line.split('\t')
    return {
        'qseq': qseq, 'sseq': sseq,
        'qstart': int(qstart), 'qend': int(qend),
        'pident': float(pident), 'mismatch': int(mismatch),
        'length': int(length),
    }


def _enumerate_substitutions(hsp, peptide_len):
    """Walk the aligned qseq/sseq strings and list every position where both
    sides are real (non-gap) residues and differ.

    Returns (substitutions_str, mismatch_count, coverage). Positions in the
    output string are 1-based and relative to the peptide itself (not the
    microprotein), since qstart marks where in the *peptide* this HSP began.
    """
    q_pos = hsp['qstart']  # 1-based position within the peptide
    subs = []
    for qc, sc in zip(hsp['qseq'], hsp['sseq']):
        if qc != '-' and sc != '-':
            if qc != sc:
                subs.append(f"{q_pos}:{qc}>{sc}")
            q_pos += 1
        elif qc != '-':
            q_pos += 1
    coverage = hsp['length'] / peptide_len if peptide_len else 0.0
    return ';'.join(subs), len(subs), coverage


def main():
    print("=" * 70)
    print("N-terminal peptide substitution computation")
    print("=" * 70)

    candidates = _select_candidates()
    print(f"N-terminal peptides to examine (from microproteins whose tryptic "
          f"evidence is confined to aa 1-2, with a BLAST match): {len(candidates)}")
    if candidates.empty:
        print("Nothing to do.")
        return

    accessions = candidates['uniprot_accession'].unique()
    ref_seqs = _ensure_reference_sequences(accessions)

    rows = []
    n_comparable, n_one_or_two, n_more = 0, 0, 0
    for i, row in enumerate(candidates.itertuples(index=False), 1):
        subject_seq = ref_seqs.get(row.uniprot_accession)
        record = {
            'sequence': row.sequence,
            'peptide_sequence': row.peptide_sequence,
            'peptide_start': row.peptide_start,
            'uniprot_accession': row.uniprot_accession,
            'pident': None, 'coverage': None, 'mismatch_count': None,
            'substitutions': '', 'comparable': False,
        }
        if subject_seq:
            hsp = _align_peptide(row.peptide_sequence, subject_seq)
            if hsp:
                subs_str, n_subs, coverage = _enumerate_substitutions(
                    hsp, len(row.peptide_sequence))
                comparable = coverage >= MIN_COVERAGE
                record.update({
                    'pident': hsp['pident'], 'coverage': coverage,
                    'mismatch_count': n_subs, 'substitutions': subs_str,
                    'comparable': comparable,
                })
                if comparable:
                    n_comparable += 1
                    if 1 <= n_subs <= 2:
                        n_one_or_two += 1
                    elif n_subs > 2:
                        n_more += 1
        rows.append(record)
        if i % 100 == 0 or i == len(candidates):
            print(f"  aligned {i}/{len(candidates)}")

    out = pd.DataFrame(rows)
    out.to_csv(SUBSTITUTIONS_CSV, index=False)

    print("-" * 70)
    print(f"Peptides examined:                  {len(candidates)}")
    print(f"Comparable (coverage>={MIN_COVERAGE:.0%}):        {n_comparable}")
    print(f"  with 1-2 substitutions (rescues):  {n_one_or_two}")
    print(f"  with >2 substitutions:              {n_more}")
    print(f"Written: {SUBSTITUTIONS_CSV}")


if __name__ == "__main__":
    main()
