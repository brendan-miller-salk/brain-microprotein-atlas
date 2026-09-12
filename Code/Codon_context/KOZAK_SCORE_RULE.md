# Rule: which Kozak score column to use

`Results/Codon_context/kozak/kozak_context.csv` has two independent scoring methods, each reported for every smORF
(both `start_codon_class` values, `ATG` and `non-ATG`):

| | `full_kozak_*` | `downstream_kozak_*` |
|---|---|---|
| positions used | −4, −3, +4, +5, +6 | +4, +5, +6 only |
| compare ATG vs non-ATG | **never** | **yes** |
| describe/sort within one class | yes | yes |

## The rule

- **Comparing or grouping by `start_codon_class`** (ATG vs non-ATG) → use `downstream_kozak_fraction`
  / `downstream_kozak_strength` only.
- **Everything else** (describing one smORF's context, sorting/filtering within a single class,
  e.g. all ATG smORFs) → `full_kozak_strength` / `full_kozak_fraction` is fine, and uses more
  information (5 positions vs. 3).

## Why

These smORFs are defined stop-codon-to-stop-codon, so the three nucleotides immediately 5′ of a
non-ATG start *are* the upstream in-frame stop codon (TAA/TAG/TGA). All three start with T, so
position −3 is a T in 100% of non-ATG smORFs — not because of biology, but because of the
definition. Position −3 is heavily weighted in `full_kozak_*`, so `full_kozak_strength` makes
non-ATG smORFs look weaker than they are. `downstream_kozak_*` drops −3 (and −4) entirely and is
unbiased between the two classes.

Full detail: [`README.md`](README.md) § "The headline result" and § "Outputs".
