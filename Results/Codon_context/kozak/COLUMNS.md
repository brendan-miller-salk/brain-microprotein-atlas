# Data Dictionary

Generated 2026-07-30 15:36 by `kozak_pipeline.py`, from the output frames themselves.

**The one thing to know:** there are two Kozak scores, and neither is called plain
`kozak_fraction` on purpose - the prefix forces the choice at the point of use.

| | `full_kozak_*` | `downstream_kozak_*` |
|---|---|---|
| positions | -4, -3, +4, +5, +6 | +4, +5, +6 only |
| compare ATG with non-ATG | **no** | **yes** |
| compare within one class | yes | yes |

A non-ATG smORF here is defined stop-codon-to-stop-codon, so the three nucleotides at -3, -2, -1
*are* the upstream in-frame stop - TAA, TAG or TGA. All three begin with T, so -3 is
unfavourable in 100% of them by construction, at weight 2 of 7. The full score therefore
measures the ORF definition rather than biology whenever the classes are compared. The
downstream score is free of that.

Tier codes come from `non_aug_initiation_codons.md`. **Lower is stronger:**

| tier | codons | meaning |
|---|---|---|
| 0 | ATG | ATG (cognate) |
| 1 | CTG ACG GTG ATT ATA ATC TTG | Well-Established Near-Cognate |
| 2 | AAG AGG | Weak Near-Cognate |
| 3a | GCG | Non-Near-Cognate, Strongest Proteomics Evidence |
| 3b | TTT GTA AGC CGA | Non-Near-Cognate, Single-Instance Proteomics Evidence |
| 3c | GCT GGT GGG GGA GGC TTC GTC TCT TCA AGT CGT CGC CGG | Non-Near-Cognate, Inferred (not directly reported as a TIS) |

A blank tier means the codon has no documented initiation competence in mammals.


---

## `kozak_context.csv`

3,686 rows, 49 columns. One row per scored smORF - the primary table. **It carries two scores and only one of them is safe for comparing ATG with non-ATG**; see the two score families below.

### Identity

One row per scored smORF.

| column | meaning |
|---|---|
| `orf_id` | ORF identifier, as `gene_id` in the master annotation |
| `gene_name` | host gene symbol; includes the placeholders `Unnamed`, `Intergenic`, `Unknown` |
| `smorf_type` | smORF class (uORF, dORF, Iso, lncRNA, ...) |
| `discovery_origin` | how the ORF entered the atlas; empty when the input table does not carry this column |
| `protein_length_aa` | annotated protein length |
| `nt_length` | spliced CDS length; always 3x the protein length |

### Coordinates

Spliced-aware, from the GTF exon blocks. Advisory flags affect only these.

| column | meaning |
|---|---|
| `chromosome` | chromosome |
| `strand` | genomic strand |
| `genomic_coordinates` | the ORF's genomic min-max; spans introns when spliced |
| `genomic_start` | genomic coordinate of the start codon's first base |
| `n_exons` | exon blocks the ORF spans |
| `is_spliced` | `n_exons` > 1 |
| `CLICK_UCSC` | UCSC browser link, carried through from the master annotation |
| `coordinate_flags` | **list** - advisory QC. Affects `genomic_start` / `n_exons` and nothing else; the context is untouched. Empty for all but a handful of smORFs, typically a spurious short leading exon block in the GTF |

### Start codon

What is there, and how well documented it is as an initiator.

| column | meaning |
|---|---|
| `start_codon` | the ORF's first codon |
| `start_codon_class` | `ATG` or `non-ATG` - the grouping variable for every comparison |
| `initiation_tier` | tier `0`-`3c` from `non_aug_initiation_codons.md`; blank when the codon is in no tier. **Lower is stronger** |
| `initiation_tier_name` | that tier written out |
| `initiator_class` | `ATG` / `near-cognate` / `non-near-cognate` / `not-a-documented-initiator` |

### Raw context

The sequence itself. Everything below is computed from these, so reweight freely if the schemes here do not suit.

| column | meaning |
|---|---|
| `kozak_window` | the -6..+6 window as a string, `gccRccATGG`-style: upstream lowercase, codon and +4..+6 uppercase |
| `upstream_nt_available` | 0-6; how much 5' window the transcript actually provides. Below 6 the ORF starts near a transcript end |
| `pos_minus6` | base at -6 |
| `pos_minus5` | base at -5 |
| `pos_minus4` | base at -4 |
| `pos_minus3` | base at -3 |
| `pos_minus2` | base at -2 |
| `pos_minus1` | base at -1 |
| `pos_plus1` | base at +1 (start codon) |
| `pos_plus2` | base at +2 |
| `pos_plus3` | base at +3 |
| `pos_plus4` | base at +4 |
| `pos_plus5` | base at +5 |
| `pos_plus6` | base at +6. Any position the transcript does not reach is `-` |

### SCORE 1 of 2 - full

All five weighted positions, INCLUDING -3. **Valid only within one start-codon class.** A non-ATG smORF is defined stop-to-stop, so its -3 is the first base of the upstream in-frame stop and is therefore always T - unfavourable by construction, at weight 2 of 7. Comparing classes on this score measures the ORF definition, not biology.

| column | meaning |
|---|---|
| `which_score_to_use` | constant; restates the choice below in every row, so it cannot be missed by someone opening the CSV without these docs |
| `full_kozak_score` | weight achieved: +4 G and -3 purine 2 each; -4 C/A, +5 A, +6 T 1 each |
| `full_kozak_max` | weight available; positions outside the transcript drop out |
| `full_kozak_fraction` | `full_kozak_score` / `full_kozak_max` |
| `full_kozak_strength` | `strong` >= 0.65, `adequate` >= 0.35, else `weak` |

### SCORE 2 of 2 - downstream

The same scheme restricted to +4, +5 and +6, where the stop-to-stop definition constrains nothing. **This is the ATG-vs-non-ATG readout**, and the rank-biserial effect sizes for both scores are in the report's section 4.

| column | meaning |
|---|---|
| `downstream_kozak_score` | weight achieved over +4/+5/+6 |
| `downstream_kozak_max` | weight available |
| `downstream_kozak_fraction` | the comparable score |
| `downstream_kozak_strength` | same thresholds as above |

### Other readouts

Independent of the weighting scheme above.

| column | meaning |
|---|---|
| `minus3_purine` | -3 is A or G |
| `plus4_G` | +4 is G |
| `kozak_class_canonical` | Kozak's own rule (1986, 1987): `strong` = both of the above, `adequate` = exactly one, `weak` = neither. Can never be `strong` for a stop-to-stop non-ATG smORF, for the same reason as above |
| `upstream_is_inframe_stop` | -3..-1 is TAA/TAG/TGA - the per-row flag for the artefact. True for every non-ATG smORF that has the window, and for a small percentage of ATG ones |
| `consensus_matches` | positions agreeing with the full `gccRccATGG` consensus |
| `consensus_positions_testable` | positions the transcript provides |
| `consensus_match_fraction` | the ratio of the two |


---

## `kozak_by_start_codon.csv`

61 rows, 16 columns. One row per start codon.

### Group

One row per distinct start codon.

| column | meaning |
|---|---|
| `start_codon` | the codon |
| `initiation_tier` | its tier; blank when in no tier |
| `initiator_class` | `ATG` / `near-cognate` / `non-near-cognate` / `not-a-documented-initiator` |

### Statistics

Every non-ATG codon here is a stop-to-stop start, so read `mean_downstream_kozak_fraction`. The full score is present only to make the size of the -3 penalty visible.

| column | meaning |
|---|---|
| `n` | smORFs in the group |
| `mean_full_kozak_fraction` | mean full score - **within-class only**, see the primary table |
| `median_full_kozak_fraction` | its median |
| `mean_downstream_kozak_fraction` | mean downstream score - **the comparable one** |
| `pct_weighted_strong` | % with `full_kozak_strength` == strong (>= 0.65) |
| `pct_weighted_adequate` | %% adequate |
| `pct_weighted_weak` | %% weak |
| `pct_upstream_is_inframe_stop` | %% whose -3..-1 is an in-frame stop codon |
| `pct_minus3_purine` | %% with a purine at -3 |
| `pct_plus4_G` | %% with G at +4 |
| `pct_canonical_strong` | %% strong under Kozak's own two-position rule |
| `pct_canonical_weak` | %% weak under it |
| `pct_full_upstream_window` | %% with all 6 upstream positions available |


---

## `kozak_by_smorf_type.csv`

30 rows, 15 columns. One row per (smORF type x start-codon class).

### Group

One row per (smORF type x start-codon class). Because the class is a grouping column here, comparing the two rows of a type is a cross-class comparison - use the downstream score.

| column | meaning |
|---|---|
| `smorf_type` | smORF class |
| `start_codon_class` | `ATG` or `non-ATG` |

### Statistics

As above.

| column | meaning |
|---|---|
| `n` | smORFs in the group |
| `mean_full_kozak_fraction` | mean full score - **within-class only**, see the primary table |
| `median_full_kozak_fraction` | its median |
| `mean_downstream_kozak_fraction` | mean downstream score - **the comparable one** |
| `pct_weighted_strong` | % with `full_kozak_strength` == strong (>= 0.65) |
| `pct_weighted_adequate` | %% adequate |
| `pct_weighted_weak` | %% weak |
| `pct_upstream_is_inframe_stop` | %% whose -3..-1 is an in-frame stop codon |
| `pct_minus3_purine` | %% with a purine at -3 |
| `pct_plus4_G` | %% with G at +4 |
| `pct_canonical_strong` | %% strong under Kozak's own two-position rule |
| `pct_canonical_weak` | %% weak under it |
| `pct_full_upstream_window` | %% with all 6 upstream positions available |


---

## `position_base_frequencies.csv`

2,496 rows, 8 columns. Base composition at -6..+6, the PWM / sequence-logo data.

### PWM / sequence-logo data

Base composition at each context position, per group. Positions the transcript does not reach are excluded from that position's denominator, so `n_with_position` is the true n for the row - not `n_orfs`.

| column | meaning |
|---|---|
| `group_kind` | `start_codon_class`, `start_codon` (n >= 20 only) or `smorf_type` |
| `group` | the group's value |
| `n_orfs` | smORFs in the group overall |
| `position` | -6..+6; there is no position 0 |
| `base` | A, C, G or T |
| `count` | smORFs with that base at that position |
| `frequency` | `count` / `n_with_position` |
| `n_with_position` | smORFs whose transcript reaches this position - the denominator |


---

## `qc_audit.csv`

3,686 rows, 7 columns. The audit trail.

### Disposition

Every id matching the smORF ORF-id grammar and what happened to it. Checks are split by what they threaten: fatal flags mean the context itself cannot be trusted and the ORF is dropped; coordinate-only problems are advisory and appear in `coordinate_flags` in the primary table instead.

| column | meaning |
|---|---|
| `orf_id` | ORF identifier |
| `smorf_type` | smORF class |
| `start_codon_class` | `ATG` or `non-ATG` |
| `qc_flags` | **list** - fatal flags; empty when clean |
| `disposition` | `analysed` or `qc_excluded` |
| `upstream_nt_available` | as in the primary table; blank when excluded |
| `full_kozak_strength` | as in the primary table; blank when excluded |


---

## Worked example

```python
import pandas as pd
d = pd.read_csv("kozak_context.csv", low_memory=False)

# RIGHT - the comparable score
d.groupby("start_codon_class").downstream_kozak_fraction.mean()   # 0.33 vs 0.273

# WRONG - this measures the stop-to-stop definition, not initiation context
d.groupby("start_codon_class").full_kozak_fraction.mean()         # 0.431 vs 0.235

# fine: the full score within one class
d[d.start_codon_class == "ATG"].groupby("smorf_type").full_kozak_fraction.mean()
```
