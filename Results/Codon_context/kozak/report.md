# Kozak Context of Every smORF Start Codon

Generated 2026-07-30 15:36 by `kozak_pipeline.py`
from `microprotein_master.csv`.
Positions follow the usual convention: **+1..+3 are the start codon**, +4 is the first
nucleotide after it, −1 the one before it, and there is no position 0.

## 1. Question

`initiation_pipeline.py` asks where translation could start in the 1,230
smORFs that do not begin with ATG. This analysis asks a narrower question of the whole atlas:
**how good is each smORF's annotated initiation context, and how do the non-ATG starts compare
with the ATG starts around them?**

smORFs here are called with GTFtoFASTA (Martinez et al., 2019): an ORF containing an in-frame ATG
takes that ATG as its start, and only an ORF without one is reported stop-codon-to-stop-codon.
The 2,456 ATG smORFs are therefore the internal reference — same transcripts, same ORF
caller, same sequence files, but a start that is not a boundary of the calling method.

## 2. Cohort

| quantity | value |
|---|---:|
| rows in the master table | 31,218 |
| unique ids | 31,218 |
| ... not smORF ORF ids (UniProt/reference; no ORF model) | 27,532 |
| smORFs | 3,686 |
| ... excluded by QC | 0 |
| **scored** | **3,686** |
| ... ATG | 2,456 |
| ... non-ATG | 1,230 |
| with the full 6 nt upstream window | 3,425 (92.9%) |

The cohort is whatever the input table contains, so these counts move with `--master`: the
master table shipped with this repository is the released atlas, while the upstream annotation
CSV holds the full pre-release discovery set.

## 3. Three independent readouts

**Weighted score** — five positions, weighted by reported effect size (+4 G and −3 purine at
weight 2; −4 C/A, +5 A, +6 T at weight 1), divided by the weight actually available. Positions
the transcript does not reach drop out of the denominator, so ORFs near a transcript end stay
comparable. `strong` ≥ 0.65, `adequate` ≥ 0.35, else `weak`.

**Downstream-only score** — the same scheme restricted to +4, +5 and +6. Section 4 is the reason
it exists: the upstream positions are not free to vary in a stop-to-stop ORF, so only the 3′
positions support a comparison between ATG and non-ATG starts.

**Canonical Kozak** — Kozak's own two-position rule: `strong` = purine at −3 **and** G at +4;
`adequate` = exactly one; `weak` = neither. It depends on no weighting choice of ours, which is
why it is reported alongside rather than instead.

| context | ATG % | non-ATG % |
|---|---:|---:|
| strong | 20.9 | 1.7 |
| adequate | 41 | 24.1 |
| weak | 38.1 | 74.2 |

| canonical Kozak | ATG % | non-ATG % |
|---|---:|---:|
| strong | 25.7 | 0 |
| adequate | 48.3 | 30.7 |
| weak | 26 | 69.3 |

The non-ATG column of both tables looks alarming and should not be read yet — **section 4 shows
that most of it is an artefact of how these ORFs were defined.**

## 4. ATG vs non-ATG — and why −3 must be thrown out

| | ATG | non-ATG |
|---|---:|---:|
| mean **full** weighted score | 0.431 | 0.235 |
| mean **downstream-only** score (+4/+5/+6) | 0.33 | 0.273 |
| −3 purine | 58.7% | 0.0% |
| +4 G | 40.9% | 30.7% |
| −3..−1 is an in-frame stop codon¹ | 4.8% | **100.0%** |

¹ Of the ORFs that have those three nucleotides at all: 2,389 ATG and
1,050 non-ATG. The rest sit at a transcript 5′ end.

**The −3 column is not a measurement.** These non-ATG smORFs are defined stop-codon-to-stop-codon,
so the three nucleotides immediately 5′ of the start *are* the upstream in-frame stop — TAA, TAG
or TGA. Every one of them therefore begins with T at −3, which is never a purine. The
0% above is arithmetic, not biology, and −3 carries weight 2 of 7 in the weighted score.

Two things follow, and they point the same way:

- On the **downstream-only** score, which no definition constrains, the two classes are
  still separated:
  0.33 vs 0.273,
  rank-biserial r = 0.103
  (compared with r = 0.469 on the full score).
- Rescoring the **ATG** starts as if their −3 were also non-purine gives
  **0.263** — against the non-ATG group's observed
  0.235. The gap essentially closes.

So the apparent context deficit of non-ATG smORFs is largely an artefact of how the ORFs were
defined. Their 3′ context is ordinary. This does not rescue the annotated N-termini — the point
of the initiation step stands, that a stop-to-stop start is a definitional boundary rather
than a measured initiation site — but it does mean **weak weighted context is not independent
evidence for that claim**, because the definition produces it directly.

Mann–Whitney p-values: 1.67e-122 (full) and
1.15e-07 (downstream-only); χ² p =
2.04e-259 (−3 purine) and 2.12e-09 (+4 G).
With 3,686 ORFs almost any difference reaches significance, so read the effect size, not the
p-value.

## 5. By start codon

| start_codon | initiation_tier | initiator_class | n | mean_downstream_kozak_fraction | mean_full_kozak_fraction | pct_plus4_G | pct_minus3_purine |
|---|---|---|---:|---:|---:|---:|---:|
| ATG | 0 | ATG | 2,456 | 0.33 | 0.431 | 40.9 | 58.7 |
| AGA |  | not-a-documented-initiator | 41 | 0.305 | 0.268 | 29.3 | 0 |
| GGA | 3c | non-near-cognate | 41 | 0.348 | 0.281 | 51.2 | 0 |
| GCA |  | not-a-documented-initiator | 39 | 0.34 | 0.297 | 48.7 | 0 |
| AAA |  | not-a-documented-initiator | 38 | 0.237 | 0.215 | 23.7 | 0 |
| GGG | 3c | non-near-cognate | 37 | 0.311 | 0.24 | 37.8 | 0 |
| GAG |  | not-a-documented-initiator | 33 | 0.303 | 0.219 | 36.4 | 0 |
| GCT | 3c | non-near-cognate | 32 | 0.359 | 0.303 | 56.2 | 0 |
| GGC | 3c | non-near-cognate | 32 | 0.227 | 0.188 | 25 | 0 |
| AGC | 3b | non-near-cognate | 30 | 0.217 | 0.222 | 20 | 0 |
| GAT |  | not-a-documented-initiator | 30 | 0.308 | 0.242 | 33.3 | 0 |
| GAA |  | not-a-documented-initiator | 29 | 0.397 | 0.294 | 48.3 | 0 |
| AGG | 2 | near-cognate | 29 | 0.284 | 0.248 | 31 | 0 |
| ATT | 1 | near-cognate | 29 | 0.319 | 0.255 | 37.9 | 0 |
| ACA |  | not-a-documented-initiator | 27 | 0.278 | 0.226 | 29.6 | 0 |
| AGT | 3c | non-near-cognate | 27 | 0.287 | 0.267 | 44.4 | 0 |
| ATA | 1 | near-cognate | 27 | 0.259 | 0.228 | 3.7 | 0 |
| GCC |  | not-a-documented-initiator | 26 | 0.269 | 0.253 | 30.8 | 0 |
| CAG |  | not-a-documented-initiator | 26 | 0.25 | 0.208 | 34.6 | 0 |
| CTT |  | not-a-documented-initiator | 25 | 0.32 | 0.256 | 32 | 0 |
| GTC | 3c | non-near-cognate | 24 | 0.156 | 0.17 | 4.2 | 0 |
| AAT |  | not-a-documented-initiator | 24 | 0.396 | 0.28 | 45.8 | 0 |
| ACT |  | not-a-documented-initiator | 24 | 0.177 | 0.191 | 16.7 | 0 |
| AAG | 2 | near-cognate | 23 | 0.283 | 0.244 | 30.4 | 0 |
| GGT | 3c | non-near-cognate | 23 | 0.283 | 0.267 | 43.5 | 0 |

Every non-ATG codon here is a stop-to-stop start, so read `mean_downstream_kozak_fraction`;
`mean_full_kozak_fraction` is shown beside it only to make the size of the −3 penalty visible.
Full table, including the strength breakdowns, in `kozak_by_start_codon.csv`.

## 6. By smORF type

| smorf_type | ATG score | ATG n | non-ATG score | non-ATG n |
|---|---:|---:|---:|---:|
| Iso | 0.366 | 539 | 0.276 | 48 |
| dORF | 0.326 | 475 | 0.274 | 318 |
| psORF | 0.357 | 413 | 0.267 | 45 |
| lncRNA | 0.319 | 335 | 0.253 | 158 |
| iORF | 0.244 | 154 | 0.279 | 243 |
| daoORF | 0.332 | 146 | 0.287 | 94 |
| uORF | 0.297 | 118 | 0.277 | 94 |
| eORF | 0.358 | 72 | 0.245 | 55 |
| daORF | 0.384 | 69 | 0.283 | 46 |
| uoORF | 0.176 | 64 | 0.254 | 70 |
| uaoORF | 0.348 | 28 | 0.312 | 20 |
| uaORF | 0.224 | 19 | 0.295 | 22 |
| D-Iso | 0.275 | 10 | – | – |
| N-Iso | 0.325 | 10 | 0.375 | 14 |
| doORF | 0.083 | 3 | 0.333 | 3 |
| udORF | 0 | 1 | – | – |

Downstream-only score, so the two classes are comparable within each type. Counts and every
other statistic are in `kozak_by_smorf_type.csv`.

## 7. Quality control

Context scores agree exactly across all 1288 ORFs also scored by the initiation step.

Checks are split by what they actually threaten. **Fatal** — no sequence, a sequence that does
not translate to the annotated protein, a length that is not a multiple of three, or an ORF
occurring more than once in its transcript (where the string search cannot tell which copy the
upstream window came from, so −1..−6 is unreliable). 0 smORFs are excluded
on these grounds, listed with their flags in `qc_audit.csv`.

**Advisory** — the GTF exon blocks disagree with the sequence length, or the strand in the ORF
id disagrees with the GTF. Exon blocks feed `genomic_start` and `n_exons` and nothing else, so
these ORFs keep every context measure and carry the flag in `coordinate_flags` instead:
1 smORFs, typically a spurious short leading exon block in the GTF.

## 8. Reading the primary CSV: two scores, one of them a trap

Neither score is called plain `kozak_fraction`, on purpose — the prefix forces the choice at the
point of use, which is where the mistake would otherwise be made.

| | `full_kozak_*` | `downstream_kozak_*` |
|---|---|---|
| positions | −4, −3, +4, +5, +6 | +4, +5, +6 only |
| compare ATG with non-ATG | **no** — section 4 | **yes** |
| compare within one class | yes | yes |

`which_score_to_use` restates this in every row, `upstream_is_inframe_stop` flags the affected
ORFs individually, and `COLUMNS.md` describes every column of every output.

## 9. Files

| file | rows | contents |
|---|---:|---|
| `kozak_context.csv` | 3,686 | one row per smORF — the primary table |
| `kozak_by_start_codon.csv` | 61 | one row per start codon |
| `kozak_by_smorf_type.csv` | 30 | one row per (smORF type × ATG/non-ATG) |
| `position_base_frequencies.csv` | 2,496 | base composition at −6..+6, the PWM data |
| `qc_audit.csv` | 3,686 | every smORF id and its disposition |
| `COLUMNS.md` | – | data dictionary: every column of every output, grouped by family |
| `figures/*.png` | 5 | strength, the −3 artefact, per-codon, base composition, per-type |
| `logs/run_2026-07-30_153654.log` | – | full run log |

Sequence input: `codon_context_sequences.csv.gz`, built by `build_sequence_context.py`
from the ESPRESSO long-read files (not distributed — see `DATA_AVAILABILITY.md`).

## 10. Caveats

- The weights and thresholds are constructions from reported effect sizes, not values fitted to
  data. Per-position bases are in `kozak_context.csv` if you prefer to weight differently.
- **Do not compare ATG and non-ATG smORFs on the full weighted score.** Section 4 is the reason:
  −3 is determined by the stop-to-stop definition. Use `downstream_kozak_fraction`, or compare within
  a class.
- The canonical Kozak class can never be `strong` for a stop-to-stop non-ATG smORF, for the same
  reason. Its −3 arm is always unfavourable, so `plus4_G` is the only informative half.
- Context predicts initiation **efficiency**, not whether an ORF is translated. A weak context
  is evidence about rate, not about existence.
- The upstream sequence comes from the assembled transcript. 261
  smORFs begin within 6 nt of a transcript 5′ end, so their window is truncated and the missing
  positions drop out of the denominator — a limit of the assembly rather than of the ORF.
- ATG smORFs are a reference, not a gold standard: their start codons are called by the same ORF
  finder and have not been confirmed by N-terminomics either.
