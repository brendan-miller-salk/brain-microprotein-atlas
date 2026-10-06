# Peptide Specificity

Do the MS peptides of an unreviewed microprotein identify *that* entry, or could
the same spectra come from another human protein?

For the 4,814 unreviewed entries (Salk smORFs + TrEMBL), every detected peptide is
compared against all human UniProt proteins (UP000005640: Swiss-Prot + TrEMBL +
isoforms). A peptide is flagged when another protein contains the same peptide
(I = L), a same-mass rearrangement, or a deamidation-like variant, including
stretches trypsin would *not* release from that protein (possible proteolytic
fragments). Swiss-Prot entries are intentionally not flagged.

This is an annotation layer only. It does not change the gold-standard filter or
any figure count.

## Overview

- **Input**: `Code/data/microprotein_master.zip` (`peptide_sequence`, `sequence`,
  `gene_name`) and the UniProt human proteome FASTA (not shipped; download below).
- **Output**: `Code/data/unreviewed_microproteins_peptide_flags.csv` (4,814 rows,
  joined to the master by `gene_id`, one-to-one) and
  `Code/data/peptide_flags_data_dictionary.csv` (column definitions).
- **Consumer**: the dashboard's *Peptide specificity* facet, table column and
  entry-page twin table (`Results/microproteins_dashboard.py`).

| script | role |
|---|---|
| `make_peptide_flags.py` | Base flag table: per-entry verdict, match type, matching protein(s), trypsin status, distinguishing fragment ions. ~1.5 h, one core. |
| `add_twin_peptides.py` | Adds the side-by-side columns (`entry_peptides`, `matched_entry_peptide`, `twin_peptide`, `differences`, `twin_peptide_protein`, `twin_peptide_trypsin`). ~2 min, 8 processes. |

## Usage

```bash
curl -L -o up5640_iso.fasta.gz \
  "https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000005640&format=fasta&includeIsoform=true&compressed=true"
python make_peptide_flags.py --fasta up5640_iso.fasta.gz        # --limit N for a test run
python add_twin_peptides.py  --fasta up5640_iso.fasta.gz        # updates the flag table in place
```

Outputs default to `../data/`; `--master`, `--outdir`, `--flags` and `--out` override paths.
Requires Python ≥3.9, pandas, numpy.

**Reproducibility caveat.** The shipped CSV was produced by the original analysis
code, which these two scripts consolidate. `add_twin_peptides.py` reproduced the
shipped comparison columns. `make_peptide_flags.py` has not been rerun end-to-end
against the shipped table. Matches also depend on the UniProt release: the shipped
table used a copy with 169,651 sequences (169,217 distinct). Record the release
number when rerunning.

## Verdicts (`interpretation`)

| verdict | all peptides | only some peptides |
|---|---|---|
| Supported: peptides unique to this entry | 1,917 | — |
| Ambiguous: could be a fragment of the full-length protein from the same gene | 627 | 17 |
| Ambiguous: same-mass peptide in another protein; fragment ions may separate them | 306 | 55 |
| Ambiguous: could be a fragment of a different protein | 191 | 32 |
| Ambiguous: peptide shared with another protein | 41 | 2 |
| Mostly supported: differs only by a deamidation-like mass match | 13 | 2 |
| No MS evidence | 1,611 | — |

Verdicts with only some matching peptides carry the suffix "(only some peptides; others unique)".

**"Ambiguous" means unsupported by peptide specificity alone, not disproved.**

**How to read "could be a fragment".** The atlas search was fully tryptic. When a
peptide sits inside a larger protein at a position trypsin does not cut, the search
engine never generated it for that protein. The smORF's own start or stop supplies
the "tryptic" end, so the peptide was credited to the smORF.

- **Same gene** (mostly N-truncated isoform and TrEMBL fragments): the alternatives are
  a downstream-initiated proteoform or a proteolytic fragment of the full-length
  protein. N-terminal acetylation, Ribo-seq initiation or targeted MS can discriminate
  them.
- **Different protein**: the alternative is a non-specific cleavage product of an
  unrelated, often abundant protein (e.g. tubulin, actin, MIF, ALDOA).

## Rules

- An entry is never matched against itself. Matches identical to the entry's own
  sequence are excluded, so TrEMBL entries do not match their own UniProt record.
- **Same mass**: 1- and 2-residue units with identical elemental composition
  (C,H,N,O,S), plus I = L. N→D and Q→E deamidation are reported separately.
- **Trypsin**: cleavage after K/R, with protein termini counting as cleavage sites.
  No proline rule is applied, which matches the atlas search.
- **Twin shown**: one representative twin per peptide, namely its most tryptic
  occurrence, with ties broken by FASTA order. `matching_protein` lists up to 8
  genes.

## Limitations

1. **Theoretical, not spectral.** Whether the distinguishing b/y ions
   (`n_fragment_ions_that_distinguish`) were actually observed needs per-PSM spectra,
   which are not in the public release.
2. **Reference-bound.** Only human UniProt proteins are considered. Other atlas
   smORFs and non-human contaminants are not.
3. **Fully-tryptic search assumption.** The direct test of the "fragment" categories
   is a semi-specific re-search against canonical proteins + contaminants.
