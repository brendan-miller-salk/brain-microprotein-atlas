# Start-Codon Context Analysis

One pipeline, run in two steps: where microprotein translation begins, and how good the
sequence context around that start codon is.

smORFs here are called with **GTFtoFASTA** (Martinez et al., 2019). An ORF containing an
in-frame ATG takes that ATG as its start. An ORF that does not is reported
**stop-codon-to-stop-codon**, so its annotated N-terminus is where the calling method had
to stop rather than a measured initiation site.

The pipeline exists to **infer where initiation could actually occur in that non-ATG
subset** — and, as the Kozak step shows, to establish that the same stop-to-stop
boundary *manufactures* the weak upstream context those ORFs appear to have.

| step | question |
|---|---|
| `initiation_pipeline.py` | For each non-ATG smORF: where **could** translation start, what protein would that give, and is the observed tryptic peptide compatible with it? |
| `kozak_pipeline.py` | For **every** smORF, ATG and non-ATG alike: how good is the initiation context of the annotated start codon? |
| `build_sequence_context.py` | Processing step: extracts the per-ORF sequence the other two need from the ESPRESSO long-read files. |

The ATG smORFs are what make the non-ATG numbers interpretable: same transcripts, same
ORF caller, same sequence files, but a start codon that is not an artefact of the
stop-to-stop definition.

## Overview

- **Input**: `Code/data/microprotein_master.csv` (cohort, protein sequences, observed
  peptides) and `Code/data/codon_context_sequences.csv.gz` (per-ORF nucleotide
  sequence, exon blocks and flanking transcript window).
- **Output**: `Results/Codon_context/kozak/` and `Results/Codon_context/initiation/` —
  CSVs, a `report.md`, a generated `COLUMNS.md` data dictionary, figures (Kozak only)
  and a run log.
- **Main scripts**: `initiation_pipeline.py`, `kozak_pipeline.py`. Shared configuration,
  path resolution and the codon tier tables are in `codon_context_common.py`.

## Usage

```bash
# from this directory - defaults resolve to the shipped master table
python initiation_pipeline.py
python kozak_pipeline.py
```

Run the initiation step **first**: the Kozak step cross-checks its own context scores
against `initiation_candidates.csv` and merely warns when that file is absent.

Both are called automatically, in that order, by `bash run_all_analyses.sh --mode=run`
from the repo root.

Every path can be overridden, so the same code runs against the full pre-release
annotation table without editing:

```bash
python kozak_pipeline.py --master /path/to/full_annotation.csv \
                         --sequences /path/to/codon_context_sequences.csv.gz \
                         --outdir /path/to/output
```

or, for a whole session, via `$MICROPROTEIN_MASTER` and `$CODON_CONTEXT_SEQUENCES`.
`--help` lists the defaults.

### Rebuilding the sequence-context table

`codon_context_sequences.csv.gz` is derived from three ESPRESSO long-read files
(`.split_nuc`, `_ORFs.gtf`, `.nuc`, ~1.5 GB) that are too large for GitHub. They are
archived on Zenodo as `espresso_sequence_files.tar.gz` (~158 MB, alongside the `.gtf`
transcript models) — see [DATA_AVAILABILITY.md](../../DATA_AVAILABILITY.md). The derived
table ships in this repository instead, so neither pipeline needs them. To rebuild it
(only necessary when the cohort changes), extract the bundle anywhere and pass the shared
path prefix:

```bash
python build_sequence_context.py \
  --espresso-prefix .../ESPRESSO_syn52047893_hg38NCBI_N2_R0_updated.sorted_medianCPM05_GENE_CLEANED
```

Each source file supplies something the others cannot: `.split_nuc` the ORF's spliced
CDS, `_ORFs.gtf` its exon blocks (a naive genomic slice would return intronic sequence
for a spliced ORF), and `.nuc` the assembled transcript — the only thing that makes the
upstream −1..−6 context available at all, since the ORF sequence itself starts at +1.

Only 12 nt of flank either side of each ORF is kept. That is twice the reach of any
scoring position in either pipeline, so scoring against the derived table is *identical*
to scoring against the full transcript, not an approximation.

## The two things to know before using the outputs

**1. Which Kozak score.** `kozak_context.csv` carries two, and neither is called plain
`kozak_fraction` on purpose:

| | `full_kozak_*` | `downstream_kozak_*` |
|---|---|---|
| positions | −4, −3, +4, +5, +6 | +4, +5, +6 only |
| compare ATG with non-ATG | **never** | **yes** |
| describe or sort within one class | yes | yes |

A non-ATG smORF is defined stop-to-stop, so the three nucleotides at −3, −2, −1 *are*
the upstream in-frame stop codon — TAA, TAG or TGA, all of which begin with T. Position
−3 is therefore unfavourable in 100% of non-ATG smORFs by construction, at weight 2 of 7.
The full score measures the ORF definition rather than biology whenever the classes are
compared. Full detail in [KOZAK_SCORE_RULE.md](KOZAK_SCORE_RULE.md) and § 4 of the
generated report.

**2. Two column families in `initiation_summary.csv`.** The `scanning_*` block is a
sequence-only prediction that uses no peptide evidence; the `initiation_*` site lists
use peptide evidence and nothing else. They are adjacent, they share the word *tier*,
and they disagree by design — which is what lets prediction and MS evidence be
cross-tabulated instead of confused.

The generated `COLUMNS.md` in each output directory describes every column of every
output and is built from the output frames themselves, so it cannot drift.

## Codon tiers

Both pipelines classify start codons with the tier tables in
[non_aug_initiation_codons.md](non_aug_initiation_codons.md), compiled from Diaz de Arce
2018 (FACS-seq), Hecht 2017, Na 2018 (human N-terminomics) and Andreev 2022. **Lower is
stronger**: tier 1 is a well-established near-cognate initiator; tier 3c has never been
directly reported as an initiation site and is inferred only. The tables live in
`codon_context_common.py`, which is the file to edit if the reference changes — update
both together.

Tier 5 (E. coli evidence only) is deliberately not scanned: with no mammalian
confirmation it would make nearly every second codon a candidate.

## Quality control

Every input id appears in `qc_audit.csv` with a disposition — `analysed`,
`qc_excluded`, `untestable_no_peptide` or `excluded_not_smorf` — and the counts are
asserted to close, so nothing is silently discarded.

The Kozak step splits QC by what a failure actually threatens: *fatal* flags (no
sequence, translation mismatch, an ORF occurring twice in its transcript) drop the ORF,
while *advisory* flags affect only the reported coordinates and are carried per row in
`coordinate_flags`.

## Cohort size depends on the input

Run against the shipped master table, the cohort is the released atlas. Run against the
full pre-release annotation CSV, it is the complete discovery set, which is substantially
larger. Counts in the generated reports move accordingly; the scoring does not.

## Dependencies

Python: `pandas`, `numpy`, `scipy`, `matplotlib`, `biopython` (`Bio.Seq`). All are in
`environment.yml` / `requirements.txt` at the repo root.
