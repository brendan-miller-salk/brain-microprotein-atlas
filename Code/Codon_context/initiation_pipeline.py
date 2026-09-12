#!/usr/bin/env python
"""
Initiation-site pipeline for non-ATG smORFs.

smORFs here are called with GTFtoFASTA (Martinez et al., 2019). An ORF containing an in-frame
ATG takes that ATG as its start; an ORF that does not is reported stop-codon-to-stop-codon, so
its annotated N-terminus is where the calling method had to stop rather than a measured
initiation site. For that non-ATG subset the real question is not "is this ORF real?" but:

    Where could translation actually start, what protein would that give, and is the observed
    tryptic peptide compatible with it?

The pipeline answers that for every in-frame candidate initiation site, INCLUDING the annotated
first codon, scored on identical terms.

Stages
------
    1. cohort      non-ATG smORFs from the master annotation table
    2. sequences   ORF nucleotides, exon blocks and flanking context (see
                   `build_sequence_context.py`)
    3. candidates  enumerate and score every in-frame initiation site
    4. summarise   collapse to one row per ORF
    5. write       three CSVs and a report

Outputs (under `--outdir`, by default `Results/Codon_context/initiation/`)
-------
    initiation_summary.csv     one row per testable ORF; site fields are bracketed,
                               index-aligned lists so an ORF with two plausible sites reads
                               ['annotated_start', 'downstream_alternative']
    initiation_candidates.csv  one row per (ORF, candidate site) - every score, so any
                               alternative ranking can be recomputed without re-running
    qc_audit.csv               every non-ATG entry in the input and its disposition; nothing is
                               ever silently discarded

All genomic and peptide coordinates are 1-based inclusive.

Run:
    python initiation_pipeline.py
    python initiation_pipeline.py --master /path/to/full_annotation.csv
"""

from __future__ import annotations

import argparse
import logging
import os
import re
from collections import Counter
from datetime import datetime

import pandas as pd
from Bio.Seq import Seq

import codon_context_common as common
from codon_context_common import (ADEQUATE_CONTEXT, CODON_TIER, COGNATE_STATUS,
                                 CONTEXT_POSITIONS, STOP_CODONS, STRONG_CONTEXT, TIER_NAME,
                                 TIER_ORDER, TIER_RANK, TIERS, as_list, parse_list)

# ======================================================================================
# Configuration
# ======================================================================================
#
# Paths are resolved at run time from the repository root, with `--master` / `--sequences` /
# `--outdir` (and their environment variables) overriding. Nothing here is machine-specific:
# the same code runs against the master table shipped with this repository and against the full
# annotation CSV that produced it.

# --- translation and N-terminal processing --------------------------------------------
# Non-AUG initiation almost always uses Met-tRNAi, so the initiation codon is decoded as Met,
# NOT as its cognate amino acid. These four codons are documented cognate-incorporation
# exceptions. CTG, ACG, and AGG are dual-outcome: the literature documents both a Met-tRNAi
# start and a direct cognate-tRNA start for each, so both are modeled as separate outcomes.
# TTT has only a single reported cognate-incorporation instance with no Met-retained
# counterpart, so it stays cognate-only.
COGNATE_INCORPORATION = {"CTG": "L", "ACG": "T", "AGG": "R", "TTT": "F"}
DUAL_OUTCOME_CODONS = {"CTG", "ACG", "AGG"}   # both Met-tRNAi and cognate-tRNA outcomes modeled
METAP_PERMISSIVE = set("ACGPSTV")   # MetAP excises the initiator Met (Frottin 2006)
NATA_SUBSTRATES = set("AGSCVTP")    # NatA acetylates the exposed N-terminus (Na 2018)

# A scanning ribosome is taken to recognise a site when its codon is a tier 1-3c initiator and
# its context clears this fraction of the available weight. (Positions, weights and the
# strong/adequate thresholds are shared with the Kozak pipeline - see `codon_context_common`.)
SCANNING_THRESHOLD = ADEQUATE_CONTEXT
# Written verbatim into every row of the summary, so the `scanning_*` block cannot be read out
# of context by someone who opens the CSV without the docs. That was the actual failure mode.
SCANNING_SOURCE = ("sequence-only scanning model; uses no peptide evidence - "
                   "first 5'->3' tier 1-3c codon with context fraction >= %.2f" % SCANNING_THRESHOLD)

# Canonical column names. The shipped master and the upstream annotation CSV do not carry the
# same set; `common.load_cohort` supplies whatever is missing as an empty column and says so in
# the log, rather than failing on a table that is merely smaller.
COHORT_COLS = ["gene_id", "sequence", "genomic_coordinates", "smorf_type", "gene_name",
               "discovery_origin", "Database", "peptide_sequence", "tryptic_peptide",
               "is_unique_peptide", "total_spectral_count"]

log = logging.getLogger("initiation")


# ======================================================================================
# Small helpers
# ======================================================================================

def is_tryptic_in(peptide, protein):
    """
    Is `peptide` a fully-tryptic peptide of `protein` (no missed cleavages)?

    This single test carries all of the N-terminal logic. A peptide sitting at the mature
    protein's first residue is valid as its N-terminal peptide; one sitting further in must
    follow a K or R. So when Met-tRNAi puts Met at position 1, a peptide that would have begun
    with the codon's cognate residue simply is not present, and a peptide beginning at residue
    2 is only valid if MetAP removed the Met - no special cases needed.
    """
    i = protein.find(peptide)
    if i < 0:
        return False
    return ((i == 0) or (protein[i - 1] in "KR")) and \
           ((i + len(peptide) == len(protein)) or (peptide[-1] in "KR"))


def mature_protein(orf_protein, aa_offset, codon):
    """
    List of (initiator_aa, metap_cleaves, mature_sequence) outcomes predicted for a start at
    `aa_offset`. Most codons yield exactly one outcome: Met via Met-tRNAi, cleaved by MetAP when
    residue 2 is small. CTG, ACG, and AGG yield two outcomes, since the literature documents
    both a Met-tRNAi start and a direct cognate-tRNA start for each. TTT has only a single
    reported cognate-incorporation instance with no Met-retained counterpart, so it stays
    cognate-only.
    """
    body = orf_protein[aa_offset + 1:]
    cleaves = bool(body) and body[0] in METAP_PERMISSIVE
    met_outcome = ("M", cleaves, (body if cleaves else "M" + body))

    cognate = COGNATE_INCORPORATION.get(codon)
    if cognate is None:
        return [met_outcome]
    cognate_outcome = (cognate, False, cognate + body)
    return [cognate_outcome, met_outcome] if codon in DUAL_OUTCOME_CODONS else [cognate_outcome]


def initiation_context(sequence, offset):
    """Context bases and weighted score at `offset` in `sequence`. Positions outside the
    sequence are '-' and drop out of the denominator, keeping scores comparable."""
    bases, score, available = {}, 0, 0
    for name, rel, favourable, weight in CONTEXT_POSITIONS:
        i = offset + rel
        base = sequence[i] if 0 <= i < len(sequence) else "-"
        bases[name] = base
        if base != "-":
            available += weight
            score += weight if base in favourable else 0
    fraction = score / available if available else 0.0
    bases.update(score=score, max_available=available, fraction=round(fraction, 3),
                 strength=("strong" if fraction >= STRONG_CONTEXT else
                           "adequate" if fraction >= ADEQUATE_CONTEXT else "weak"),
                 window="%s|%s|%s" % (sequence[max(0, offset - 6):offset].rjust(6, "-"),
                                      sequence[offset:offset + 3],
                                      sequence[offset + 3:offset + 6].ljust(3, "-")))
    return bases


def exons_in_transcription_order(exons, strand):
    return sorted(exons) if strand == "+" else sorted(exons, reverse=True)


def offset_to_genomic(exons, strand, nt_offset):
    """Genomic position of nucleotide `nt_offset` of the ORF, walking spliced exon blocks."""
    remaining = nt_offset
    for lo, hi in exons_in_transcription_order(exons, strand):
        length = hi - lo + 1
        if remaining < length:
            return lo + remaining if strand == "+" else hi - remaining
        remaining -= length
    raise IndexError("offset %d beyond ORF" % nt_offset)


def retained_blocks(exons, strand, nt_offset):
    """Exon blocks covering the ORF from `nt_offset` to its 3' end."""
    remaining, out = nt_offset, []
    for lo, hi in exons_in_transcription_order(exons, strand):
        length = hi - lo + 1
        if remaining >= length:
            remaining -= length
            continue
        if remaining:
            lo, hi = (lo + remaining, hi) if strand == "+" else (lo, hi - remaining)
            remaining = 0
        out.append((lo, hi))
    return out


# ======================================================================================
# Stage 1 - cohort
# ======================================================================================

def load_cohort(master_csv):
    """Unique non-ATG smORFs, plus the non-smORF non-ATG rows kept for the QC trail."""
    log.info("[1/5] reading %s", os.path.basename(master_csv))
    cohort, other, totals = common.load_cohort(
        master_csv, COHORT_COLS, log,
        keep=lambda chunk: ~chunk.sequence.astype(str).str.startswith("M"))
    other = other.drop_duplicates("sequence").reset_index(drop=True)
    log.info("      %d rows, %d unique ORFs -> %d non-ATG smORFs (+%d non-smORF, QC only)",
             totals["input_rows"], totals["input_ids"], len(cohort), len(other))
    return cohort, other, {"input_rows": totals["input_rows"],
                           "input_orfs": totals["input_ids"]}


# ======================================================================================
# Stage 2 - sequences
# ======================================================================================

def validate(orf, protein, seq):
    """Every check that must hold before an ORF is analysed. Returns a list of flags."""
    flags = []
    nt, exons = seq["nt"], seq["exons"]
    if nt is None:
        flags.append("nt_sequence_not_found")
    if exons is None:
        flags.append("exon_structure_not_found")
    if nt is None or exons is None:
        return flags
    if len(nt) % 3:
        flags.append("nt_length_not_multiple_of_3")
    if sum(hi - lo + 1 for lo, hi in exons) != len(nt):
        flags.append("exon_length_sum_!=_nt_length")
    translated = str(Seq(nt).translate())
    if translated != protein:
        flags.append("translation_mismatch")
    if "*" in translated:
        flags.append("internal_stop_codon")
    if seq["tx_copies"] != 1:
        flags.append("orf_not_unique_in_transcript")
    if re.search(r"([+-])chr", orf).group(1) != seq["strand"]:
        flags.append("strand_mismatch_id_vs_gtf")
    return flags


# ======================================================================================
# Stage 3 - candidate initiation sites
# ======================================================================================

def peptides_for(row, protein):
    """Annotated peptides present in the ORF. `peptide_sequence` (a list) is primary;
    `tryptic_peptide` is the fallback when it is absent."""
    peps = [str(p).strip().upper() for p in parse_list(row.get("peptide_sequence"))
            if str(p).strip()]
    source = "peptide_sequence"
    if not peps:
        fallback = row.get("tryptic_peptide")
        if isinstance(fallback, str) and fallback.strip():
            peps, source = [fallback.strip().upper()], "tryptic_peptide"
        else:
            return [], "none"
    return [p for p in peps if p in protein], source


def score_candidates(orf, protein, seq, peptides):
    """
    Every in-frame candidate initiation site for one ORF.

    Position 1 is always enumerated, whatever its codon - excluding it would be exactly the
    asymmetry this analysis exists to remove. Downstream sites must be tier 1-3c initiators.
    """
    nt, exons, strand = seq["nt"], seq["exons"], seq["strand"]
    # The ORF plus its retained flanks: the part of the assembled transcript any context
    # position can reach. Identical to indexing into the whole transcript - see FLANK_KEPT.
    context_seq, orf_start = common.context_sequence(seq)
    out = []
    for nt_offset in range(0, len(nt), 3):
        codon = nt[nt_offset:nt_offset + 3]
        if len(codon) < 3:
            continue
        tier = CODON_TIER.get(codon)
        if nt_offset and not tier:
            continue

        aa_offset = nt_offset // 3
        ctx = initiation_context(context_seq, orf_start + nt_offset)
        blocks = retained_blocks(exons, strand, nt_offset)
        genomic_start = offset_to_genomic(exons, strand, nt_offset)

        # A dual-outcome codon (CTG/ACG/AGG) yields two outcomes here, so this position
        # contributes two candidate dicts - everything else contributes exactly one.
        for initiator, metap, mature in mature_protein(protein, aa_offset, codon):
            # One test does all the work: every peptide must be fully tryptic in the mature
            # protein. Evaluated per outcome, since the two outcomes have different sequences.
            compatible = bool(peptides) and all(is_tryptic_in(p, mature) for p in peptides)

            out.append({
                "orf_id": orf, "aa_position": aa_offset + 1, "is_annotated_start": nt_offset == 0,
                "codon": codon, "tier": tier or "", "tier_name": TIER_NAME.get(tier, ""),
                "cognate_status": COGNATE_STATUS.get(tier, ""),
                "genomic_start": genomic_start,
                "new_start": min(lo for lo, _ in blocks), "new_end": max(hi for _, hi in blocks),
                "exon_blocks": as_list("%d-%d" % b for b in blocks),
                "context_minus4": ctx["minus4"], "context_minus3": ctx["minus3"],
                "context_plus4": ctx["plus4"], "context_plus5": ctx["plus5"],
                "context_plus6": ctx["plus6"], "context_window": ctx["window"],
                "context_score": ctx["score"], "context_max": ctx["max_available"],
                "context_fraction": ctx["fraction"], "context_strength": ctx["strength"],
                "initiator_aa": initiator, "metap_cleaves_met": metap,
                "predicted_sequence": initiator + protein[aa_offset + 1:],
                "predicted_sequence_met_cleaved": mature,
                "predicted_length_aa": len(mature),
                "nata_acetylated": bool(mature) and mature[0] in NATA_SUBSTRATES,
                "n_aa_removed": aa_offset,
                "peptide_compatible": compatible,
                "scanning_recognised": bool(tier) and ctx["fraction"] >= SCANNING_THRESHOLD,
            })
    return out


# ======================================================================================
# Stage 4 - per-ORF summary
# ======================================================================================

SITE_COLS = ["site_type", "aa_position", "codon", "cognate_status", "tier", "tier_name",
             "initiator_aa", "predicted_sequence", "predicted_sequence_met_cleaved",
             "predicted_length_aa", "metap_cleaves_met", "genomic_start"]
SITE_RENAME = {"aa_position": "initiation_codon_position", "codon": "initiation_codon",
               "tier": "initiation_tier", "tier_name": "initiation_tier_name"}


def summarise(row, protein, peptides, pep_source, candidates):
    """
    One row per ORF.

    A candidate becomes a *site* only if it is peptide-compatible AND its codon can actually
    initiate. An annotated start whose codon is in no tier describes a protein with no known
    way to begin, so it is not a site - such ORFs fall through to `no_site_reason`.
    """
    # A dual-outcome codon (CTG/ACG/AGG) at position 1 produces two candidate dicts here, both
    # sharing codon/tier/context (only initiator identity differs) - `annotated` is any one of
    # them for those shared fields; `annotated_compatible` checks across all of them.
    annotated_variants = [c for c in candidates if c["is_annotated_start"]]
    annotated = annotated_variants[0]
    annotated_compatible = any(c["peptide_compatible"] for c in annotated_variants)

    compatible = [c for c in candidates if c["peptide_compatible"]]
    sites = []
    sites.extend(dict(c, site_type="annotated_start")
                 for c in annotated_variants if c["peptide_compatible"] and c["tier"])
    downstream = [c for c in compatible if not c["is_annotated_start"]]
    if downstream:
        best_key = min((TIER_RANK[c["tier"]], c["aa_position"]) for c in downstream)
        sites.extend(dict(c, site_type="downstream_alternative") for c in downstream
                     if (TIER_RANK[c["tier"]], c["aa_position"]) == best_key)

    # Peptide compatibility describes the PEPTIDE only - whether some protein starting there
    # could have produced it. Whether that codon can actually initiate is a separate question,
    # carried by `has_supported_initiation_site`. Keeping them apart is the point.
    if annotated_compatible and downstream:
        compat = "both"
    elif annotated_compatible:
        compat = "first_codon_only"
    elif downstream:
        compat = "downstream_only"
    else:
        compat = "neither"

    if sites:
        reason = ""
    elif annotated_compatible and not annotated["tier"]:
        # The peptide fits, but the codon has no documented way to initiate - so this is a
        # non-answer, indistinguishable in practice from having no site at all.
        reason = "annotated_start_codon_not_an_initiator"
    elif any(all(p in protein[c["aa_position"] - 1:] for p in peptides) for c in candidates):
        reason = "no_site_gives_a_tryptic_match"
    else:
        reason = "peptide_not_retained_downstream"

    # Sequence-only prediction: the scanning model, which mirrors the mechanism.
    # Deliberately named `scanning_*` in the output rather than `predicted_*`: these columns
    # use no peptide evidence at all, and sit next to the `initiation_*` site lists, which use
    # nothing but. Keeping the two families distinguishable by prefix is the point.
    recognised = [c for c in candidates if c["scanning_recognised"]]
    scanned = min(recognised, key=lambda c: c["aa_position"]) if recognised else None

    rec = {
        "orf_id": row["gene_id"], "gene_name": row.get("gene_name"),
        "smorf_type": row.get("smorf_type"), "discovery_origin": row.get("discovery_origin"),
        "original_sequence": protein,
        "original_length_aa": len(protein),
        "genomic_coordinates": row.get("genomic_coordinates"),
        "annotated_start_codon": annotated["codon"],
        "annotated_start_is_initiator": bool(annotated["tier"]),
        "annotated_context_strength": annotated["context_strength"],
        "peptide_sequence": as_list(peptides), "peptide_source": pep_source,
        "n_peptides": len(peptides),
        "peptide_initiation_compatibility": compat,
        "has_supported_initiation_site": bool(sites),
        "n_initiation_sites": len(sites),
        "no_site_reason": reason,
        "scanning_source": SCANNING_SOURCE,
        "scanning_aa_position": scanned["aa_position"] if scanned else "",
        "scanning_codon": scanned["codon"] if scanned else "",
        "scanning_tier": scanned["tier"] if scanned else "",
        "scanning_tier_name": scanned["tier_name"] if scanned else "",
        "scanning_is_annotated_start": scanned["is_annotated_start"] if scanned else "",
        "scanning_peptide_compatible": scanned["peptide_compatible"] if scanned else "",
        "n_candidate_sites": len(candidates),
        "n_compatible_candidates": len(compatible),
        "is_unique_peptide": row.get("is_unique_peptide"),
        "total_spectral_count": row.get("total_spectral_count"),
    }
    for col in SITE_COLS:
        rec[SITE_RENAME.get(col, col)] = as_list(s[col] for s in sites)
    return rec


# ======================================================================================
# Stage 5 - report
# ======================================================================================

def write_report(paths, summary, candidates, qc, totals, log_path):
    n = len(summary)

    def tbl(counter, total, header):
        return "\n".join(["| %s | n | %% |" % header, "|---|---:|---:|"] +
                         ["| %s | %d | %.1f%% |" % (k, v, 100.0 * v / total)
                          for k, v in counter.most_common()])
    tiers = Counter(t for row in summary.initiation_tier for t in parse_list(row))
    text = f"""# Initiation-Site Analysis of non-ATG smORFs

Generated by `initiation_pipeline.py` from `{os.path.basename(paths['master'])}`.
All genomic and peptide coordinates are **1-based inclusive**.

## 1. Question

smORFs here are called with GTFtoFASTA (Martinez et al., 2019). An ORF with no in-frame ATG is
reported stop-codon-to-stop-codon, so its annotated N-terminus follows from the calling method
rather than from a measured initiation site, and this pipeline exists to infer where initiation
could actually occur. For every in-frame
candidate site - **including the annotated first codon** - it scores the codon's
evidence tier, its initiation context, the mature protein predicted by Met-tRNAi/MetAP, and
whether the observed tryptic peptides are compatible with it.

## 2. Cohort

| quantity | value |
|---|---:|
| rows in the master table | {totals['input_rows']:,} |
| unique ORFs in the master table | {totals['input_orfs']:,} |
| non-ATG smORFs | {totals['cohort']:,} |
| ... excluded by QC | {totals['qc_failed']:,} |
| ... with no annotated peptide (untestable) | {totals['untestable']:,} |
| **testable ORFs analysed** | **{n:,}** |
| candidate initiation sites scored | {len(candidates):,} |

The cohort is whatever the input table contains, so these counts move with `--master`: the
master table shipped with this repository is the released atlas, while the upstream annotation
CSV holds the full pre-release discovery set.

## 3. Quality control

Every analysed ORF passed all checks: the ORF nucleotide sequence translates exactly to the
annotated protein, exon block lengths sum to the nucleotide length, the strand in the ORF
identifier matches the GTF, and each ORF occurs exactly once in its transcript. Anything that
failed is listed with its flags in `qc_audit.csv`. Two tiers are structurally empty and asserted
so: **tier 0 (ATG)** - no cohort ORF contains a downstream in-frame ATG, because the ORF caller
emits ATG-initiated forms as separate `_M` records - and **tier 4 (stop readthrough)** - a
stop-to-stop ORF has no internal in-frame stop.

## 4. The annotated starts sit in weak context

{tbl(Counter(summary.annotated_context_strength), n, "annotated-start context")}

This uses no peptide evidence at all. It is the clearest single piece of support for treating
these N-termini as artefacts of the stop-to-stop definition rather than as initiation sites.

## 5. Where can initiation occur?

An ORF gets one *site* per plausible initiation point: the annotated start when it is
peptide-compatible and its codon can initiate, and the best downstream alternative (highest
tier, earliest position) when one is compatible.

{tbl(Counter(summary.peptide_initiation_compatibility), n, "peptide compatible with")}

| | ORFs |
|---|---:|
| **has a supported initiation site** | **{int(summary.has_supported_initiation_site.sum()):,}** |
| no supported site | {int((~summary.has_supported_initiation_site).sum()):,} |

{tbl(Counter(summary.no_site_reason[summary.no_site_reason != ""]), max(1, int((~summary.has_supported_initiation_site).sum())), "reason when there is no site")}

Because a shortened protein is a suffix of the longer one, a peptide is rarely discriminating -
most ORFs have several compatible candidates, which is why compatibility is reported alongside
the sequence-based prediction rather than used to select one.

### Tier of the sites found

{tbl(tiers, max(1, sum(tiers.values())), "initiation tier")}

## 6. Reading the summary CSV: two column families

`initiation_summary.csv` answers two different questions side by side, and they must not be
conflated. Every column belongs to one family or the other:

| | `scanning_*` | `initiation_*` and the other site columns |
|---|---|---|
| question | where would a scanning ribosome start? | which starts does the observed peptide permit? |
| uses peptide evidence | **no** | **yes** |
| shape | one scalar per ORF | bracketed list, one entry per site |
| example | `scanning_tier = 3a` | `initiation_tier = ['3a']` |
| may be blank because | no codon clears tier + context ({int((summary.scanning_aa_position == "").sum())} ORFs) | no site is peptide-supported ({int((~summary.has_supported_initiation_site).sum()):,} ORFs) |

They disagree by design - that is what lets sequence prediction and MS evidence be
cross-tabulated instead of confused. `scanning_source` restates the first column in full in
every row. A per-column data dictionary is in `COLUMNS.md`.

The tier codes are the same in both families, from `non_aug_initiation_codons.md`, and **lower
is stronger**: tier 1 is a well-established near-cognate initiator, tier 3c has never been
directly reported as a translation initiation site and is inferred only.

## 7. Files

| file | rows | contents |
|---|---:|---|
| `initiation_summary.csv` | {n:,} | one row per testable ORF; site fields are bracketed, index-aligned lists |
| `initiation_candidates.csv` | {len(candidates):,} | one row per (ORF, candidate site) - every score, so any alternative ranking can be recomputed |
| `qc_audit.csv` | {len(qc):,} | every non-ATG entry and its disposition; nothing silently discarded |
| `COLUMNS.md` | - | data dictionary: every column of every output, grouped by family |
| `logs/` | - | full run log |

Sequence input: `{os.path.basename(paths['sequences'])}`, built by `build_sequence_context.py`
from the ESPRESSO long-read files (not distributed - see `DATA_AVAILABILITY.md`).

## 8. Caveats

- The context weights and the scanning threshold are constructions from reported effect sizes,
  not values fitted to data. Raw per-position bases are in the candidates file.
- Predictions are sequence-based and unvalidated. N-terminomics is what would confirm any
  individual site.
- "Not an initiator" means no documented competence in mammals, not proven incapacity.
"""
    with open(paths["report"], "w") as fh:
        fh.write(text)


# ======================================================================================
# Stage 5b - data dictionary
# ======================================================================================
#
# Grouped by *family*, because the families are the thing that was confusing: the summary
# answers two different questions side by side, one using peptide evidence and one not.
# `write_columns_md` asserts that these tables describe every column actually written, so the
# dictionary cannot silently drift from the outputs.

SUMMARY_FAMILIES = [
    ("Identity", "One row per testable ORF. These identify it.", [
        ("orf_id", "ORF identifier, as `gene_id` in the master annotation"),
        ("gene_name", "host gene symbol; includes the placeholders `Unnamed`, `Intergenic`, `Unknown`"),
        ("smorf_type", "smORF class (uORF, dORF, Iso, lncRNA, ...)"),
        ("discovery_origin", "how the ORF entered the atlas; empty when the input table does not carry this column"),
        ("genomic_coordinates", "the ORF's genomic min-max; spans introns when spliced"),
        ("chromosome", "chromosome"),
        ("strand", "genomic strand"),
        ("original_sequence", "the annotated stop-to-stop protein, exactly as in the master annotation"),
        ("original_length_aa", "length of `original_sequence`"),
    ]),
    ("Annotated start", "The stop-to-stop first codon, scored on the same terms as every "
                        "downstream alternative. Uses no peptide evidence.", [
        ("annotated_start_codon", "the ORF's first codon"),
        ("annotated_start_is_initiator", "True when that codon is in tier 1-3c, i.e. has documented initiation competence in mammals"),
        ("annotated_context_strength", "`strong` / `adequate` / `weak` initiation context at position 1"),
    ]),
    ("Peptide evidence", "What was observed by MS, and what it permits. Uses peptides only.", [
        ("peptide_sequence", "**list** - the tryptic peptides observed for this ORF"),
        ("peptide_source", "which annotation column the peptides came from"),
        ("n_peptides", "number of peptides"),
        ("is_unique_peptide", "carried through from the master annotation"),
        ("total_spectral_count", "carried through from the master annotation"),
        ("peptide_initiation_compatibility", "`both` / `first_codon_only` / `downstream_only` / `neither` - which starts the peptides can come from. About the peptide alone; says nothing about codon competence"),
    ]),
    ("Supported sites", "The conclusion: starts that are BOTH peptide-compatible AND have a "
                        "codon that can initiate. Site columns are bracketed, index-aligned "
                        "lists - element i of each describes the same site. Parse with "
                        "`ast.literal_eval`. Empty lists `[]` mean no supported site.", [
        ("has_supported_initiation_site", "**the main filter.** True when there is >=1 such site"),
        ("n_initiation_sites", "how many"),
        ("no_site_reason", "why there is none; blank when there is one"),
        ("n_candidate_sites", "in-frame candidates enumerated for this ORF, before any filter. A CTG/ACG/AGG position counts twice here (Met-initiated and cognate-initiated are both enumerated)"),
        ("n_compatible_candidates", "how many of those the peptides permit"),
        ("site_type", "**list** - `annotated_start` or `downstream_alternative`"),
        ("initiation_codon_position", "**list** - 1-based amino-acid position of the site in `original_sequence`"),
        ("initiation_codon", "**list** - the codon at that site"),
        ("initiation_tier", "**list** - its evidence tier, `1` / `2` / `3a` / `3b` / `3c`. **Lower is stronger**"),
        ("initiation_tier_name", "**list** - that tier written out"),
        ("cognate_status", "**list** - `near-cognate` (tiers 1-2) or `non-near-cognate` (tiers 3a-3c)"),
        ("initiator_aa", "**list** - the residue actually placed at position 1. Met for most codons; CTG/ACG/AGG each contribute two site entries here (Met AND the cognate residue - L/T/R respectively - since both mechanisms are documented); TTT is cognate-only (Phe), with no Met-retained counterpart reported"),
        ("predicted_sequence", "**list** - the protein as synthesised from that site, initiator at position 1"),
        ("predicted_sequence_met_cleaved", "**list** - the mature protein after MetAP; equals the above when MetAP does not act"),
        ("predicted_length_aa", "**list** - length of the mature protein"),
        ("metap_cleaves_met", "**list** - whether MetAP removes the initiator Met (residue 2 is A/C/G/P/S/T/V)"),
        ("genomic_start", "**list** - genomic coordinate of the site, spliced-aware"),
    ]),
    ("Scanning model", "A SEPARATE, SEQUENCE-ONLY prediction: where a scanning ribosome would "
                       "start. Uses NO peptide evidence. One scalar per ORF, not a list. These "
                       "are expected to disagree with the site columns - keeping the two apart "
                       "is what lets prediction and MS evidence be cross-tabulated.", [
        ("scanning_source", "constant; restates the rule in every row so the block cannot be read out of context"),
        ("scanning_aa_position", "position of the first codon that is a tier 1-3c initiator AND clears context fraction %.2f, walking 5'->3'. Blank when no codon does" % SCANNING_THRESHOLD),
        ("scanning_codon", "the codon there"),
        ("scanning_tier", "its evidence tier. Same scale as `initiation_tier`, **lower is stronger** - a `3c` here means the model's best guess rests on a codon never directly reported as an initiation site"),
        ("scanning_tier_name", "that tier written out"),
        ("scanning_is_annotated_start", "True when the model picks position 1, i.e. sequence alone supports the annotated N-terminus"),
        ("scanning_peptide_compatible", "whether the peptides happen to permit the model's pick. The one place the two families meet - a cross-tabulation, not an input to either"),
    ]),
]

CANDIDATE_FAMILIES = [
    ("Identity", "One row per (ORF, candidate site, outcome). Position 1 is always present, "
                 "whatever its codon; downstream rows are tier 1-3c codons only. A CTG/ACG/AGG "
                 "candidate contributes two rows at the same position (Met-initiated and "
                 "cognate-initiated); every other codon contributes exactly one.", [
        ("orf_id", "ORF identifier"), ("chromosome", "chromosome"), ("strand", "genomic strand"),
        ("aa_position", "1-based amino-acid position of this candidate"),
        ("is_annotated_start", "True for position 1"),
        ("n_aa_removed", "residues lost relative to the annotated start"),
    ]),
    ("Codon", "What is there and how well documented it is.", [
        ("codon", "the codon"),
        ("tier", "evidence tier `1`-`3c`, blank when the codon is in no tier (possible only at position 1). **Lower is stronger**"),
        ("tier_name", "that tier written out"),
        ("cognate_status", "`near-cognate` (tiers 1-2) or `non-near-cognate` (tiers 3a-3c)"),
    ]),
    ("Initiation context", "Scored on the assembled transcript, so upstream positions are real "
                           "sequence rather than missing. Weights: +4 G and -3 purine 2 each; "
                           "-4 C/A, +5 A, +6 T 1 each.", [
        ("context_minus4", "base at -4"), ("context_minus3", "base at -3"),
        ("context_plus4", "base at +4"), ("context_plus5", "base at +5"),
        ("context_plus6", "base at +6"),
        ("context_window", "the window as a string, `gccRccATGG`-style"),
        ("context_score", "weight achieved"),
        ("context_max", "weight available; positions outside the transcript drop out"),
        ("context_fraction", "`context_score` / `context_max`"),
        ("context_strength", "`strong` >= %.2f, `adequate` >= %.2f, else `weak`" % (STRONG_CONTEXT, ADEQUATE_CONTEXT)),
    ]),
    ("Predicted protein", "What initiating here would actually produce.", [
        ("initiator_aa", "residue placed at position 1"),
        ("predicted_sequence", "protein as synthesised"),
        ("metap_cleaves_met", "whether MetAP removes the initiator Met"),
        ("predicted_sequence_met_cleaved", "mature protein after MetAP"),
        ("predicted_length_aa", "its length"),
        ("nata_acetylated", "whether NatA would acetylate the exposed N-terminus"),
    ]),
    ("Evidence and model", "The two independent readouts, per candidate.", [
        ("peptide_compatible", "peptide evidence: every observed peptide is fully tryptic in this candidate's mature protein"),
        ("scanning_recognised", "sequence only: tier 1-3c codon AND context fraction >= %.2f. The scanning model picks the earliest row where this is True" % SCANNING_THRESHOLD),
    ]),
    ("Coordinates", "Spliced-aware; derived from the GTF exon blocks.", [
        ("genomic_start", "genomic coordinate of this candidate's first base"),
        ("new_start", "genomic min of the retained blocks"),
        ("new_end", "genomic max of the retained blocks"),
        ("exon_blocks", "**list** - `start-end` blocks retained by a start here"),
    ]),
]

QC_FAMILIES = [
    ("Disposition", "Every non-ATG entry in the input and what happened to it. Nothing is "
                    "silently discarded; the pipeline asserts the counts close.", [
        ("orf_id", "ORF identifier"), ("smorf_type", "smORF class"),
        ("Database", "source database from the master annotation"),
        ("in_cohort", "True when the id parses as a smORF ORF id"),
        ("qc_flags", "**list** - reasons for exclusion; empty when clean"),
        ("disposition", "`analysed` / `untestable_no_peptide` / `qc_excluded` / `excluded_not_smorf`"),
        ("n_peptides", "peptides found; blank when not applicable"),
        ("has_site", "`has_supported_initiation_site` for analysed ORFs; blank otherwise"),
    ]),
]


def write_columns_md(paths, summary, candidates, qc):
    """Data dictionary, generated from the frames so it cannot drift from the actual columns."""
    files = [("initiation_summary.csv", summary, SUMMARY_FAMILIES,
              "One row per testable ORF. **It answers two different questions side by side** - "
              "see the two families below - so check which family a column belongs to before "
              "using it."),
             ("initiation_candidates.csv", candidates, CANDIDATE_FAMILIES,
              "One row per (ORF, candidate site), unranked and unfiltered. Every score is here, "
              "so any alternative ranking can be recomputed without re-running the pipeline."),
             ("qc_audit.csv", qc, QC_FAMILIES, "The audit trail.")]

    parts = ["""# Data Dictionary

Generated by `initiation_pipeline.py`, from the output frames themselves.

**The one thing to know:** `initiation_summary.csv` carries two parallel families of columns.
The `scanning_*` block is a sequence-only prediction that uses no peptide evidence. The
`initiation_*` site lists use peptide evidence and nothing else. They are adjacent, they share
the word *tier*, and they disagree by design.

Tier codes are shared by both families and come from `non_aug_initiation_codons.md`.
**Lower is stronger:**

| tier | codons | meaning |
|---|---|---|
%s

A blank tier means the codon has no documented initiation competence in mammals - possible only
at the annotated start, since downstream candidates are enumerated by tier.
""" % ("\n".join("| %s | %s | %s |" % (t, " ".join(TIERS[t]), TIER_NAME[t]) for t in TIER_ORDER))]

    for name, frame, families, note in files:
        described = [c for _, _, cols in families for c, _ in cols]
        missing = [c for c in frame.columns if c not in described]
        extra = [c for c in described if c not in frame.columns]
        assert not missing, "%s: columns with no description: %s" % (name, missing)
        assert not extra, "%s: described columns that do not exist: %s" % (name, extra)

        parts.append("\n---\n\n## `%s`\n\n%d rows, %d columns. %s\n"
                     % (name, len(frame), len(frame.columns), note))
        for family, fnote, cols in families:
            parts.append("### %s\n\n%s\n\n| column | meaning |\n|---|---|\n%s\n"
                         % (family, fnote,
                            "\n".join("| `%s` | %s |" % (c, d) for c, d in cols)))

    parts.append("""
---

## Worked example

```python
import pandas as pd, ast
d = pd.read_csv("initiation_summary.csv")
r = d[d.orf_id == "%s"].iloc[0]

r.annotated_start_codon          # %s   - the stop-to-stop first codon
r.scanning_tier                  # %s   - sequence only: tier of where a ribosome would start
ast.literal_eval(r.initiation_tier)   # %s - peptide evidence: tiers of the supported sites
```

The second and third lines answer different questions. Neither is derived from the other.
""" % tuple(_worked_example(summary)))

    with open(paths["columns"], "w") as fh:
        fh.write("\n".join(parts))


def _worked_example(summary):
    """A real row where the scanning model and the peptide-supported sites both exist."""
    ok = summary[(summary.scanning_tier != "") & summary.has_supported_initiation_site]
    r = (ok if len(ok) else summary).iloc[0]
    return (r.orf_id, r.annotated_start_codon, r.scanning_tier, r.initiation_tier)


# ======================================================================================
# Main
# ======================================================================================

def resolve_paths(args):
    outdir = args.outdir or os.path.join(common.DEFAULT_RESULTS_ROOT, "initiation")
    return {
        "master": common.resolve(args.master, common.ENV_MASTER, common.DEFAULT_MASTER),
        "sequences": common.resolve(args.sequences, common.ENV_SEQUENCES,
                                    common.DEFAULT_SEQUENCE_CONTEXT),
        "outdir": outdir,
        "logs": os.path.join(outdir, "logs"),
        "summary": os.path.join(outdir, "initiation_summary.csv"),
        "candidates": os.path.join(outdir, "initiation_candidates.csv"),
        "qc": os.path.join(outdir, "qc_audit.csv"),
        "columns": os.path.join(outdir, "COLUMNS.md"),
        "report": os.path.join(outdir, "report.md"),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    common.add_common_args(parser)
    args = parser.parse_args(argv)
    paths = resolve_paths(args)

    log_path = common.setup_logging("initiation", paths["logs"], [paths["outdir"]])
    log.info("=== initiation pipeline ===")
    log.info("      master:    %s", paths["master"])
    log.info("      sequences: %s", paths["sequences"])
    log.info("      outdir:    %s", paths["outdir"])

    cohort, non_smorf, totals = load_cohort(paths["master"])
    log.info("[2/5] sequence context")
    seqs = common.load_sequence_context(paths["sequences"], set(cohort.gene_id), log)

    log.info("[3/5] scoring candidate initiation sites")
    summary_rows, candidate_rows, qc_rows = [], [], []
    n_atg = n_stop = 0

    for row in cohort.itertuples():
        orf, protein = str(row.gene_id), str(row.sequence).upper()
        seq = seqs[orf]
        flags = validate(orf, protein, seq)
        qc = {"orf_id": orf, "smorf_type": row.smorf_type, "Database": row.Database,
              "in_cohort": True, "qc_flags": as_list(flags)}
        if flags:
            qc.update(disposition="qc_excluded", n_peptides="", has_site="")
            qc_rows.append(qc)
            log.warning("      QC exclusion %s: %s", orf, flags)
            continue

        n_atg += sum(1 for i in range(3, len(seq["nt"]), 3) if seq["nt"][i:i + 3] == "ATG")
        n_stop += sum(1 for i in range(3, len(seq["nt"]), 3) if seq["nt"][i:i + 3] in STOP_CODONS)

        rowd = row._asdict()
        peptides, pep_source = peptides_for(rowd, protein)
        cands = score_candidates(orf, protein, seq, peptides)
        for c in cands:
            c["chromosome"], c["strand"] = seq["chrom"], seq["strand"]
        candidate_rows.extend(cands)

        if not peptides:
            qc.update(disposition="untestable_no_peptide", n_peptides=0, has_site="")
            qc_rows.append(qc)
            continue

        rec = summarise(rowd, protein, peptides, pep_source, cands)
        rec["chromosome"], rec["strand"] = seq["chrom"], seq["strand"]
        summary_rows.append(rec)
        qc.update(disposition="analysed", n_peptides=len(peptides),
                  has_site=rec["has_supported_initiation_site"])
        qc_rows.append(qc)

    for row in non_smorf.itertuples():
        qc_rows.append({"orf_id": row.gene_id, "smorf_type": row.smorf_type,
                        "Database": row.Database, "in_cohort": False,
                        "qc_flags": as_list(["not_a_smorf_orf_id"]),
                        "disposition": "excluded_not_smorf", "n_peptides": "", "has_site": ""})

    # Structural invariants of a stop-to-stop non-ATG cohort. A failure here means the input
    # ORFs were not called that way, in which case the tier tables no longer describe them.
    assert n_atg == 0, "tier 0 not empty: %d downstream in-frame ATG" % n_atg
    assert n_stop == 0, "tier 4 not empty: %d internal in-frame stops" % n_stop
    flat = [c for cs in TIERS.values() for c in cs]
    assert len(flat) == len(set(flat)) == 27, "tier tables overlap or are the wrong size"

    summary = pd.DataFrame(summary_rows)
    candidates = pd.DataFrame(candidate_rows)
    qc = pd.DataFrame(qc_rows)

    log.info("[4/5] %d testable ORFs | %d candidate sites | %d QC rows",
             len(summary), len(candidates), len(qc))
    counts = Counter(qc.disposition)
    totals.update(cohort=len(cohort), qc_failed=counts["qc_excluded"],
                  untestable=counts["untestable_no_peptide"])
    assert counts["analysed"] + counts["qc_excluded"] + counts["untestable_no_peptide"] \
        == len(cohort), "ORF accounting does not close"

    log.info("[5/5] writing outputs")
    summary.to_csv(paths["summary"], index=False)
    candidates.to_csv(paths["candidates"], index=False)
    qc.to_csv(paths["qc"], index=False)
    write_report(paths, summary, candidates, qc, totals, log_path)
    write_columns_md(paths, summary, candidates, qc)

    log.info("      supported initiation site: %d / %d",
             int(summary.has_supported_initiation_site.sum()), len(summary))
    log.info("      compatibility: %s", dict(Counter(summary.peptide_initiation_compatibility)))
    log.info("      annotated-start context: %s", dict(Counter(summary.annotated_context_strength)))
    log.info("done. log: %s", log_path)
    return summary, candidates, qc


if __name__ == "__main__":
    main()
