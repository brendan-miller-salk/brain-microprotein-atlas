# Non-AUG Translation Initiation Codons — Comprehensive Tiered Reference

## Overview

All 64 codons ranked by evidence for translation initiation in eukaryotic/mammalian systems. Compiled from four systematic studies:

- **Diaz de Arce 2017** [1] — FACS-seq in mouse pre-B cells; tested all 9 single-mismatch codons across every -4 to +4 context
- **Hecht 2017** [2] — All 64 codons in E. coli; GFP and nanoluciferase reporters
- **Na 2017** [3] — Human N-terminomics (TAILS); mass spectrometry of N-terminal peptides from HEK293T, HUVEC, colon, substantia nigra
- **Andreev 2022** [4] — Review of non-AUG initiation in mammals

---

## Tier 0 — Canonical Start Codon

| Codon | AA | Mismatches from AUG | Efficiency | Notes |
|-------|-----|---------------------|------------|-------|
| AUG | Met | 0 | 100% (reference) | Canonical; Met-tRNAi initiation |

---

## Tier 1 — Well-Established Near-Cognates (1 mismatch, strong evidence)

| Codon | AA | Mismatch position | Max efficiency (% AUG) | Key evidence |
|-------|-----|-------------------|------------------------|--------------|
| CUG | Leu | Position 1 (A→C) | ~50% | Strongest non-AUG; 2–3× any other; c-Myc, PTEN, BAG1, FGF2, VEGFA |
| ACG | Thr | Position 2 (U→C) | ~40% | 2nd most efficient; Na 2017 (multiple genes); one instance decoded as Thr (CLK2) |
| GUG | Val | Position 1 (A→G) | ~20% | EIF4G2 uses GUG exclusively; well-established |
| AUU | Ile | Position 3 (G→U) | 5–10% | PTEN uses AUU as most efficient non-AUG TIS |
| AUA | Ile | Position 3 (G→A) | 5–10% | Well-established |
| AUC | Ile | Position 3 (G→C) | 5–10% | Well-established; annotated as start in bacterial genomes |
| UUG | Leu | Position 1 (A→U) | 5–10% | Well-established; annotated as start in bacterial genomes |

**Notes:**
- CUG is the most efficient non-AUG codon by a wide margin (2–3× any other) [4]
- ACG (position 2 mismatch) is 2nd most efficient, despite AAG/AGG (also position 2) being nearly undetectable — position 2 mismatch efficiency is codon-specific, not uniform [1]
- Almost all initiation at these codons uses Met-tRNAi; Met is incorporated at position 1 [3]
- CUG and ACG are dual-outcome: alongside the standard Met-tRNAi start, CUG can additionally use
  Leu-tRNA directly (Starck 2012) and ACG can additionally use Thr-tRNA directly (one instance,
  Na 2017) — both a Met-initiated and a cognate-initiated protein are modeled for these two
  codons, not one or the other

---

## Tier 2 — Weak Near-Cognates (1 mismatch, very low efficiency)

| Codon | AA | Mismatch position | FACS-seq efficiency | Proteomics evidence |
|-------|-----|-------------------|---------------------|---------------------|
| AAG | Lys | Position 2 (U→A) | Undetectable ("too close to background") [1] | NOT found as TIS in Na 2017 [3] |
| AGG | Arg | Position 2 (U→A) | Undetectable ("too close to background") [1] | Found as TIS in Na 2017 (TUBB3, CASE 2 — cognate Arg incorporated) [3] |

**Notes:**
- These are near-cognates (1 mismatch) but the least efficient of all single-mismatch codons [1]
- FACS-seq could not measure sequence-context dependence because expression was too close to background [1]
- **AAG was NOT found as a TIS in Na 2017.** It appears only as a P1' codon (first column) in Table 1, not as a TIS (-1 column). No proteomics evidence supports AAG as a TIS.
- **AGG was found as a TIS in Na 2017** (TUBB3), but as a CASE 2 entry — the cognate amino acid (Arg) was incorporated directly, not Met. This is one of only two documented cases of cognate-aa incorporation at a near-cognate codon (the other being ACG→Thr in CLK2). Evidence for AGG's cognate-incorporation is a single instance, no stronger than ACG's — so AGG is modeled the same way as CUG/ACG: both a Met-initiated and an Arg-initiated protein are plausible, not Arg-only.
- Met-tRNAi initiates at most near-cognate codons; Met is NOT cleaved by MetAP before K or R (large side chains)
- For K and R starts specifically: M-K and M-R are 2-aa peptides — too short for MS detection. Trypsin cleaves after K/R, so the detectable tryptic peptide starts at the residue AFTER K/R.

---

## Tier 3 — Non-Near-Cognates with Human Proteomics Evidence (2+ mismatches)

**CORRECTED** — Verified against actual Table 1 data from Na 2017 [3]. The original version of this table contained multiple errors (see correction notes below).

### How to read Na 2017 Table 1

Na 2017 Table 1 has two codon columns: "Codon at -1 position" and "Codon at first position." The TIS is in the **-1 column** for most entries (CASE 1): Met-tRNAi initiates at the -1 codon, Met is incorporated, MetAP cleaves Met when P1' is permissive, and the observed peptide starts at the P1' residue encoded by the **first column**. This was verified by two independent checks:
- peptide[0] matches the cognate aa of the first-column codon for all 68 entries
- The -1 column gives 91.9% near-cognate for N-terminal extensions, matching the paper's stated ~92%

Three TIS assignment cases exist:
- **CASE 1** (most common): TIS = -1 codon, Met incorporated, MetAP cleaves Met. Requires P1' (first codon aa) to be MetAP-permissive (A, C, G, P, S, T, V).
- **CASE 2** (rare): TIS = first codon, cognate aa incorporated (not Met). Occurs when P1' is NOT permissive and peptide starts with non-M. Examples: PPIA (UUU→Phe), TUBB3 (AGG→Arg), CLK2 (ACG→Thr).
- **CASE 3** (Met retained): TIS = first codon (AUG), Met retained. Example: METTL9.

### Tier 3a — Strongest non-near-cognate evidence (multiple instances)

| Codon | AA | Mismatches from AUG | Instances in Na 2017 | TIS case | MetAP cleaves Met? |
|-------|-----|---------------------|----------------------|----------|---------------------|
| GCG | Ala | 2 (A→G pos1, U→C pos2) | GDI1 (×2), KANSL1L (×2) — 4 instances | CASE 1 (Met cleaved) | Yes (A permissive) |

**GCG is the only non-near-cognate TIS with multiple instances in Na 2017 Table 1.** It appears 4 times: twice for GDI1 (N-terminal extension) and twice for KANSL1L (uORF), across multiple samples.

### Tier 3b — Single-instance non-near-cognate TIS codons

| Codon | AA | Mismatches from AUG | Instance in Na 2017 | TIS case | MetAP cleaves Met? |
|-------|-----|---------------------|----------------------|----------|---------------------|
| UUU | Phe | 2 (A→U pos1, U→U pos2, G→U pos3) | PPIA (N-terminal extension) | CASE 2 (cognate Phe incorporated) | N/A (no Met) |
| GUA | Val | 2 (A→G pos1, U→U pos2, G→A pos3) | FAM193B (uORF) | CASE 1 (Met cleaved) | Yes (A permissive) |
| AGC | Ser | 2 (A→A pos1, U→G pos2, G→C pos3) | NDUFB9 (uORF) | CASE 1 (Met cleaved) | Yes (A permissive) |
| CGA | Arg | 3 (A→C pos1, U→G pos2, G→A pos3) | HDAC2 (uORF) | CASE 1 (Met cleaved) | Yes (A permissive) |

**Total: 8 instances across 5 unique non-near-cognate TIS codons.**

### Correction notes (errors in the original version of this table)

The original Tier 3 table contained serious errors. Of 15 claimed non-NC TIS codons:
- **13 were fabricated** — not found as TIS in Na 2017 Table 1: ACA, CAC, CAG, CCA, CGG, GAA, GAC, GAG, GAU, GCA, GCC, GGA, GGC
- **1 had wrong gene assignments**: GCG (originally listed HDAC2; actual genes are GDI1 and KANSL1L only)
- **1 was correct**: AGC (NDUFB9)
- **3 were missed** (found in Table 1 but not in the original table): CGA (HDAC2), GUA (FAM193B), UUU (PPIA)

Root causes of the errors:
1. **Column confusion**: GCC, GCA, GGC, GGA, GGG, GGU, AGU, UCU, UCA appear in the **first column** (P1' codons), not the -1 column (TIS). They were incorrectly listed as TIS codons.
2. **Upstream codon confusion**: CGG (METTL9), CAC (TUBB3), GAA (PPIA) appear in the -1 column but are NOT the TIS — the TIS is in the first column for these CASE 2/3 entries.
3. **Entirely fabricated entries**: ACA (MED15), CAG (FCRLB), CCA (VGLL4), GAC (NDUFA2, TP53), GAG (SEPT9, MFSD4B), GAU (C11orf96) do not appear anywhere in Table 1.
4. **Mismatch count errors**: Several codons originally listed as "2 mismatches" actually have 3 (GAA, CAC, CGA, GCC, GGA, GAU, GGC, CCA).

### Note on the paper's Discussion

The Discussion of Na 2017 states: "TISs that start... at non-cognate AUG codons such as GCG, GCC, and GGC." However, in Table 1, GCC appears only as a P1' codon (first column, 10 instances) and GGC does not appear at all. This likely refers to supplemental data (Table 1 contains only the 41 peptides found in ≥2 of 4 samples; the full dataset has 90 peptides) or uses "TIS" loosely to include P1' codons. The ~92% near-cognate statistic confirms the -1 column is the TIS column.

### Note on Table 2 (Met-retained search)

Na 2017 Table 2 searched for Met-retained peptides by substituting near-cognate codons (ACG, CUG, GUG, UUG) plus four non-near-cognate control codons (**CUA, GCC, GCG, GUC**) to AUG. All 33 peptides found came from near-cognate codons; zero came from the four non-NC controls. This is expected: non-NC TIS with permissive P1' → MetAP cleaves Met → Met-retained form doesn't exist.

### Tier 3c — Codons with same amino acid as verified TIS but not directly reported

| Codon | AA | Mismatches | Notes |
|-------|-----|-----------|-------|
| GCU | Ala | 3 | Same aa as GCG; not reported in Table 1 |
| GGU | Gly | 3 | Not reported as TIS; GGC/GGA/GGG are P1' codons only |
| GGG | Gly | 2 | Appears as P1' codon (first column) in Table 1, not as TIS |
| GGA | Gly | 3 | Appears as P1' codon (first column) in Table 1, not as TIS |
| GGC | Gly | 3 | Does not appear in Table 1 |
| UUC | Phe | 3 | Same aa as UUU; not reported |
| GUC | Val | 2 | Tested as control in Table 2 (substituted to AUG); 0 hits |
| GUG | Val | 1 | Near-cognate (Tier 1), not non-NC |
| UCU | Ser | 3 | Appears as P1' codon (first column) in Table 1, not as TIS |
| UCA | Ser | 3 | Appears as P1' codon (first column) in Table 1, not as TIS |
| AGU | Ser | 2 | Appears as P1' codon (first column) in Table 1, not as TIS |
| CGU | Arg | 3 | Same aa as CGA; not reported |
| CGC | Arg | 3 | Same aa as CGA; not reported |
| CGG | Arg | 2 | Appears as -1 codon for METTL9, but TIS is AUG (CASE 3); not a TIS |

---

## Tier 4 — Stop Codon Readthrough as TIS

| Codon | Type | Mismatches from AUG | Evidence |
|-------|------|---------------------|----------|
| UGA | Stop | 3 | NOT found in Na 2017 Table 1 [3]; Hecht 2017: measurable in E. coli [2] |
| UAG | Stop | 3 | Hecht 2017: 0.007% of AUG in E. coli [2] |
| UAA | Stop | 3 | Not reported as TIS |

**Note:** The original version of this table claimed UGA was found in Na 2017 (FCRLB). This was fabricated — neither UGA nor FCRLB appear anywhere in Na 2017 Table 1. No stop codon has been confirmed as a TIS in mammalian proteomics.

---

## Tier 5 — E. coli Evidence Only (not confirmed in mammals)

Hecht 2017 [2] found that **47 of 64 codons** can initiate translation in E. coli at levels significantly above background:

| Efficiency tier | Codons | % of AUG |
|----------------|--------|----------|
| Canonical | AUG, GUG, UUG | 10–100% |
| Near-cognate | AUA, AUC, AUU, CUG | 0.1–3% |
| Non-near-cognate | 40 codons | 0.007–0.1% |
| Undetectable | 17 codons | Below background |

The 40 non-near-cognate codons with measurable E. coli initiation have not been confirmed in mammalian systems but suggest the ribosome can initiate at essentially any codon at some low rate. Mammalian confirmation is lacking for most of these.

---

## Context Factors That Modulate Non-AUG Initiation Efficiency

Non-AUG initiation is more context-dependent than AUG initiation [1, 4]. The following positions have the strongest effects:

| Position | Optimal nucleotide | Effect size | Notes |
|----------|-------------------|-------------|-------|
| **+4** | **G** | Up to 10× for CUG | Most critical for non-AUG; G→A can decrease 10× |
| **-3** | A or G (purine) | Up to 6× enhancement | Also important for AUG but less so |
| **-4** | C or A | ~70% increase | Largely specific to non-AUG; minimal effect on AUG |
| **+5** | A | Moderate | More important for non-AUG than AUG |
| **+6** | U | Moderate | More important for non-AUG than AUG |

**RNA secondary structure:** A stem-loop 14–16 nt downstream of the start codon enhances non-AUG initiation by pausing the scanning ribosome, increasing dwell time at the codon [4]. POLGARF CUG initiation is stimulated 3-fold by a downstream stem-loop, reaching 60–70% of AUG efficiency [4].

**Optimal non-AUG context (predicted):** (C/A)(A/G)(C/A)(A/C/G)CUG GAU [4]

**Leaky scanning:** Weak Kozak at upstream AUG sites (especially -3 ≠ purine, +4 ≠ G) increases the fraction of scanning ribosomes that bypass to downstream non-AUG sites [4].

---

## Met-tRNAi and MetAP Cleavage Rules

### Initiating amino acid
- **Almost all non-AUG initiation uses Met-tRNAi** — Met is incorporated at position 1 regardless of the codon [2, 3]
- **Dual-outcome codons (both Met-tRNAi and cognate-tRNA initiation are documented, so both are modeled):**
  - CUG → Met, or Leu directly (Leu-tRNA, Starck 2012)
  - ACG → Met, or Thr directly (one instance, CLK2, Na 2017)
  - AGG → Met, or Arg directly (one instance, TUBB3, Na 2017 — CASE 2); evidence is a single instance, no stronger than ACG's, so it is modeled the same way rather than as cognate-only
- **Cognate-only exception (no Met-retained counterpart reported, so only the cognate outcome is modeled):**
  - UUU → Phe (one instance, PPIA, Na 2017 — CASE 2)
- In stop-to-stop proteogenomic databases, Met is NOT included at ORF position 1 — the database entry starts with the cognate amino acid

### MetAP cleavage (Frottin 2006)
MetAP cleaves the initiator Met when the P1' residue (second amino acid) has a small side chain (radius of gyration ≤ 1.29 Å):

| P1' residue | MetAP cleaves? | Exposed N-terminus |
|-------------|----------------|---------------------|
| A (Ala) | Yes (efficient) | A |
| C (Cys) | Yes (efficient) | C |
| G (Gly) | Yes (efficient) | G |
| P (Pro) | Yes (efficient) | P |
| S (Ser) | Yes (efficient) | S |
| T (Thr) | Yes (less efficient) | T |
| V (Val) | Yes (less efficient) | V |
| K (Lys) | **No** | M-K (Met retained) |
| R (Arg) | **No** | M-R (Met retained) |
| Q (Gln) | **No** | M-Q (Met retained) |
| E (Glu) | **No** | M-E (Met retained) |
| D (Asp) | **No** | M-D (Met retained) |
| N (Asn) | **No** | M-N (Met retained) |
| H (His) | **No** | M-H (Met retained) |
| I (Ile) | **No** | M-I (Met retained) |
| L (Leu) | **No** | M-L (Met retained) |
| M (Met) | **No** | M-M (Met retained) |
| F (Phe) | **No** | M-F (Met retained) |
| W (Trp) | **No** | M-W (Met retained) |
| Y (Tyr) | **No** | M-Y (Met retained) |

### Implications for proteogenomics
- **P, G, A, S, C, T, V starts:** Consistent with Met-tRNAi initiation + MetAP cleavage. The exposed N-terminus matches the cognate amino acid.
- **K, R, Q, E, D, N, H, I, L, F, W, Y starts:** Met should be retained (M-aa). If the detected peptide starts with these amino acids without Met, either: (a) the peptide is internal (not the true N-terminus), (b) a non-Met-tRNAi mechanism is used (rare, only CUG/ACG/AGG documented), or (c) another aminopeptidase removes the M-aa dipeptide.
- **K and R starts specifically:** M-K and M-R are 2-aa peptides — too short for MS detection. Trypsin cleaves after K/R, so the detectable tryptic peptide starts at the residue AFTER K/R. This means a "K start" in a stop-to-stop database is indistinguishable from an internal K cleavage site.

---

## N-terminal Acetylation

NatA acetylates N-termini after Met removal when the exposed residue is: A, G, S, C, V, T, P [3]. This covers the P, G, A starts observed in proteogenomics. Acetylated N-terminal peptides are enriched by N-terminomics methods (TAILS, COFRADIC) but may be missed in standard shotgun searches depending on whether acetylation is included as a variable modification.

---

## Summary Table — All Codons by Tier

| Tier | Codons | Mismatches | Mammalian evidence | Confidence |
|------|--------|-----------|-------------------|------------|
| 0 | AUG | 0 | Canonical | Highest |
| 1 | CUG, ACG, GUG, AUU, AUA, AUC, UUG | 1 | FACS-seq + N-terminomics + functional studies | High |
| 2 | AAG, AGG | 1 | AAG: FACS-seq undetectable, no proteomics evidence; AGG: FACS-seq undetectable, 1 proteomics instance (TUBB3, CASE 2) | Low-moderate |
| 3a | GCG | 2 | N-terminomics (4 instances: GDI1 ×2, KANSL1L ×2) | Moderate |
| 3b | UUU, GUA, AGC, CGA | 2–3 | N-terminomics (1 instance each) | Low-moderate |
| 3c | GCU, GGU, GGG, GGA, GGC, UUC, GUC, UCU, UCA, AGU, CGU, CGC, CGG | 2–3 | Same aa as Tier 3a/3b codons; not directly reported as TIS | Low |
| 4 | UGA, UAG, UAA | 3 | No mammalian proteomics evidence; UAG in E. coli only | Very low |
| 5 | ~40 codons | 2–3 | E. coli only (Hecht 2017) | Unconfirmed in mammals |

---

## References

1. Diaz de Arce AJ, Noderer WL, Wang CL. Complete motif analysis of sequence requirements for translation initiation at non-AUG start codons. *Nucleic Acids Res.* 2018;46(2):985-994. doi:10.1093/nar/gkx1114
2. Hecht A, Glasgow JE, Jaschke PR, et al. Measurements of translation initiation from all 64 codons in E. coli. *Nucleic Acids Res.* 2017;45(6):3614-3623. doi:10.1093/nar/gkx070
3. Na CH, Barbhuiya MA, Kim MS, et al. Discovery of noncanonical translation initiation sites through mass spectrometric analysis of protein N termini. *Genome Res.* 2018;28(1):25-36. doi:10.1101/gr.226050.117
4. Andreev DE, Loughran G, Fedorova AD, et al. Non-AUG translation initiation in mammals. *Genome Biol.* 2022;23(1):111. doi:10.1186/s13059-022-02674-2
5. Ichihara K, Matsumoto A, Nishida H, et al. Combinatorial analysis of translation dynamics reveals eIF2 dependence of translation initiation at near-cognate codons. *Nucleic Acids Res.* 2021;49(12):e68. doi:10.1093/nar/gkab549
6. Frottin F, Martinez A, Peynot P, et al. The proteomics of N-terminal methionine cleavage. *Mol Cell Proteomics.* 2006;5(12):2336-2349. doi:10.1074/mcp.M600225-MCP200
7. Starck SR, Jiang V, Pavon-Eternod M, et al. Leucine-tRNA initiates at CUG start codons for protein synthesis and presentation by MHC class I. *Science.* 2012;336(6089):1719-1723. doi:10.1126/science.1220270
