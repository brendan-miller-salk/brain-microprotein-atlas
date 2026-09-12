#!/usr/bin/env python
"""
Kozak / initiation-context analysis for EVERY smORF in the brain microprotein atlas.

The other step of this pipeline, `initiation_pipeline.py`, asks where translation could
start in the smORFs that do NOT begin with ATG. This one asks a narrower question of the whole
atlas:

    How good is the initiation context of each smORF's annotated start codon, and how do the
    non-ATG starts compare with the ATG starts that surround them?

The ATG smORFs are the point. They are the internal reference: same transcripts, same ORF caller
(GTFtoFASTA, Martinez et al., 2019), same sequence pipeline - but their start is an in-frame ATG,
not the stop-to-stop boundary GTFtoFASTA falls back to when an ORF has no ATG. Scoring both on
identical terms is what makes the non-ATG numbers interpretable.

Stages
------
    1. cohort      every smORF in the master annotation table (ATG and non-ATG)
    2. sequences   ORF nucleotides, exon blocks and the upstream window (see
                   `build_sequence_context.py`)
    3. score       Kozak context of each annotated start codon
    4. aggregate   per-codon, per-smORF-type and positional base-frequency tables
    5. write       CSVs, figures and a report

Outputs (under `--outdir`, by default `Results/Codon_context/kozak/`)
-------
    kozak_context.csv              one row per smORF - the primary table
    kozak_by_start_codon.csv       one row per start codon
    kozak_by_smorf_type.csv        one row per (smORF type x ATG/non-ATG)
    position_base_frequencies.csv  base composition at -6..+6, the PWM/logo data
    qc_audit.csv                   every smORF-pattern id and its disposition

All genomic coordinates are 1-based inclusive. Context positions follow the usual convention:
+1..+3 are the start codon, +4 is the first nucleotide after it, -1 the one before it; there is
no position 0.

Run:
    python kozak_pipeline.py
    python kozak_pipeline.py --master /path/to/full_annotation.csv
"""

from __future__ import annotations

import argparse
import ast
import logging
import os
import re
from collections import Counter
from datetime import datetime

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio.Seq import Seq
from scipy import stats

import codon_context_common as common
from codon_context_common import (ADEQUATE_CONTEXT, CONTEXT_POSITIONS, COGNATE_STATUS,
                                 STOP_CODONS, STRONG_CONTEXT, fmt_p, md_table)

# ======================================================================================
# Configuration
# ======================================================================================
#
# Paths are resolved at run time from the repository root, with `--master` / `--sequences` /
# `--outdir` (and their environment variables) overriding. Nothing here is machine-specific.

# --- initiation codon tiers -----------------------------------------------------------
# The tier tables live in `codon_context_common` so the two steps cannot drift apart; tier
# "0" is ATG itself, which the initiation pipeline never had to represent because ATG ORFs were
# out of its scope.
TIERS = {"0": ("ATG",), **common.TIERS}
TIER_ORDER = ["0"] + common.TIER_ORDER
CODON_TIER = {c: t for t, cs in TIERS.items() for c in cs}
TIER_NAME = {"0": "ATG (cognate)", **common.TIER_NAME}
CODON_CLASS = {"0": "ATG", **COGNATE_STATUS}

# --- initiation context ---------------------------------------------------------------
# WEIGHTED SCORE. Positions and weights are shared with the initiation pipeline
# (`CONTEXT_POSITIONS` in `codon_context_common`): +4 G is the strongest single determinant (up
# to 10x for CUG), -3 purine next (up to 6x), then -4 C/A (~70%), +5 A and +6 U (moderate).

# DOWNSTREAM-ONLY SCORE. The non-ATG smORFs are defined stop-codon-to-stop-codon, so the three
# nucleotides at -3, -2, -1 ARE the upstream in-frame stop codon - which makes -3 a T in every
# one of them, by construction rather than by biology. Any score using -3 therefore penalises
# every non-ATG ORF automatically and cannot be compared with an ATG start. This second score
# uses only positions 3' of the codon, where no such constraint exists, and is the fair
# ATG-vs-non-ATG comparison. `upstream_is_inframe_stop` flags the affected ORFs per row.
DOWNSTREAM_POSITIONS = (("plus4", "G", 2), ("plus5", "A", 1), ("plus6", "T", 1))
# Neither score is named plain `kozak_*` in the output: the prefix forces a choice at the point
# of use, which is exactly where the mistake would otherwise be made. This string is written
# into every row of kozak_context.csv so the choice cannot be missed by someone opening the CSV
# without the docs.
SCORE_GUIDANCE = ("downstream_kozak_fraction to compare ATG with non-ATG; "
                  "full_kozak_fraction only within one class "
                  "(it includes -3, which the stop-to-stop definition fixes to T "
                  "in every non-ATG smORF)")

# CANONICAL KOZAK. Kozak's own two-position rule (Kozak 1986, 1987), reported alongside the
# weighted score because it is what most readers mean by "strong Kozak" and it depends on no
# weighting choices of ours: strong = purine at -3 AND G at +4; adequate = exactly one; weak =
# neither. The full consensus is gccRccAUGG.
KOZAK_CONSENSUS = {"minus6": "G", "minus5": "C", "minus4": "C", "minus3": "AG",
                   "minus2": "C", "minus1": "C", "plus4": "G"}
WINDOW_POSITIONS = [-6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6]

# Canonical column names. The shipped master and the upstream annotation CSV do not carry the
# same set; `common.load_cohort` supplies whatever is missing as an empty column and says so in
# the log, rather than failing on a table that is merely smaller.
COHORT_COLS = ["gene_id", "sequence", "genomic_coordinates", "smorf_type", "gene_name",
               "discovery_origin", "Database", "protein_length", "CLICK_UCSC"]

log = logging.getLogger("kozak")


# ======================================================================================
# Small helpers
# ======================================================================================

def position_base(upstream, orf_nt, pos):
    """
    Base at Kozak position `pos`, or '-' when the transcript does not reach that far.

    Positive positions are 1-based into the ORF (+1..+3 are the start codon); negative ones
    count back from the nucleotide immediately 5' of it. There is no position 0.
    """
    if pos > 0:
        return orf_nt[pos - 1] if pos - 1 < len(orf_nt) else "-"
    return upstream[pos] if len(upstream) >= -pos else "-"


def kozak_context(upstream, orf_nt):
    """
    Every context measure for one start codon.

    Three readouts are produced and none is derived from another:
      * the WEIGHTED score, which uses the reported effect sizes across five positions;
      * the DOWNSTREAM-ONLY score, which uses only +4/+5/+6 and so is immune to the
        stop-to-stop artefact at -3;
      * the CANONICAL Kozak class, which is Kozak's own -3/+4 rule.
    Positions the transcript does not reach are '-' and drop out of the denominators, so scores
    stay comparable between ORFs with and without a full 5' window.
    """
    bases = {("plus%d" if p > 0 else "minus%d") % abs(p): position_base(upstream, orf_nt, p)
             for p in WINDOW_POSITIONS}

    score = available = 0
    for name, _, favourable, weight in CONTEXT_POSITIONS:
        if bases[name] != "-":
            available += weight
            score += weight if bases[name] in favourable else 0
    fraction = score / available if available else 0.0

    d_score = d_available = 0
    for name, favourable, weight in DOWNSTREAM_POSITIONS:
        if bases[name] != "-":
            d_available += weight
            d_score += weight if bases[name] in favourable else 0
    d_fraction = d_score / d_available if d_available else 0.0

    purine_m3 = bases["minus3"] in "AG"
    g_p4 = bases["plus4"] == "G"
    matched = sum(1 for k, v in KOZAK_CONSENSUS.items() if bases[k] != "-" and bases[k] in v)
    testable = sum(1 for k in KOZAK_CONSENSUS if bases[k] != "-")

    upstream_shown = upstream[-6:] if upstream else ""
    bases.update(
        full_kozak_score=score, full_kozak_max=available, full_kozak_fraction=round(fraction, 3),
        full_kozak_strength=("strong" if fraction >= STRONG_CONTEXT else
                        "adequate" if fraction >= ADEQUATE_CONTEXT else "weak"),
        downstream_kozak_score=d_score, downstream_kozak_max=d_available,
        downstream_kozak_fraction=round(d_fraction, 3),
        downstream_kozak_strength=("strong" if d_fraction >= STRONG_CONTEXT else
                             "adequate" if d_fraction >= ADEQUATE_CONTEXT else "weak"),
        upstream_is_inframe_stop=(upstream[-3:] in STOP_CODONS),
        minus3_purine=purine_m3, plus4_G=g_p4,
        kozak_class_canonical=("strong" if purine_m3 and g_p4 else
                               "weak" if not purine_m3 and not g_p4 else "adequate"),
        consensus_matches=matched, consensus_positions_testable=testable,
        consensus_match_fraction=round(matched / testable, 3) if testable else 0.0,
        upstream_nt_available=len(upstream_shown),
        # gccRccATGG-style display: upstream lowercase, codon and downstream uppercase.
        kozak_window="%s%s%s" % (upstream_shown.lower(),
                                 orf_nt[:3], orf_nt[3:6]))
    return bases


# ======================================================================================
# Stage 1 - cohort
# ======================================================================================

def load_cohort(master_csv):
    """
    Every unique smORF in the master annotation table, ATG and non-ATG alike.

    smORFs are identified by their ORF-id grammar (`...(+|-)chrN:start-end_F:f_P:p[_M]`); the
    remaining ids are UniProt/reference proteins that have no ORF model here and cannot be
    given a transcript context.
    """
    log.info("[1/5] reading %s", os.path.basename(master_csv))
    cohort, other, totals = common.load_cohort(master_csv, COHORT_COLS, log)
    cohort["sequence"] = cohort.sequence.astype(str).str.upper()
    cohort["starts_with_M"] = cohort.sequence.str.startswith("M")

    n_other = other.gene_id.astype(str).nunique()
    log.info("      %d rows, %d unique ids -> %d smORFs (%d ATG, %d non-ATG); "
             "%d non-smORF ids ignored",
             totals["input_rows"], totals["input_ids"], len(cohort),
             int(cohort.starts_with_M.sum()), int((~cohort.starts_with_M).sum()), n_other)
    return cohort, {"input_rows": totals["input_rows"], "input_ids": totals["input_ids"],
                    "non_smorf_ids": n_other}


# ======================================================================================
# Stage 2 - sequences
# ======================================================================================

def validate(orf, protein, seq):
    """
    (fatal, advisory) flags for one ORF.

    FATAL means the Kozak context itself cannot be trusted: no sequence, a sequence that does
    not translate to the annotated protein, or an ORF that occurs more than once in its
    transcript - in which case the string search cannot tell which copy the upstream window
    belongs to, so the -1..-6 bases are unreliable.

    ADVISORY means only the reported genomic coordinates are affected. The exon blocks are used
    for `genomic_start` / `n_exons` and for nothing else, so a discrepancy there leaves every
    context measure intact and is no reason to drop the ORF.
    """
    fatal, advisory = [], []
    nt, exons = seq["nt"], seq["exons"]
    if nt is None:
        fatal.append("nt_sequence_not_found")
    if exons is None:
        fatal.append("exon_structure_not_found")
    if nt is None or exons is None:
        return fatal, advisory
    if len(nt) < 6:
        fatal.append("orf_shorter_than_context_window")
    if len(nt) % 3:
        fatal.append("nt_length_not_multiple_of_3")
    if str(Seq(nt).translate()) != protein:
        fatal.append("translation_mismatch")
    if seq["tx_copies"] != 1:
        fatal.append("orf_not_unique_in_transcript")
    if sum(hi - lo + 1 for lo, hi in exons) != len(nt):
        advisory.append("exon_length_sum_!=_nt_length")
    if re.search(r"([+-])chr", orf).group(1) != seq["strand"]:
        advisory.append("strand_mismatch_id_vs_gtf")
    return fatal, advisory


# ======================================================================================
# Stage 3 - score every annotated start
# ======================================================================================

def score_rows(cohort, seqs):
    """One scored row per smORF, plus a QC row for every id including the excluded ones."""
    rows, qc_rows = [], []
    for row in cohort.itertuples():
        orf, protein = row.gene_id, row.sequence
        seq = seqs[orf]
        fatal, advisory = validate(orf, protein, seq)
        qc = {"orf_id": orf, "smorf_type": row.smorf_type,
              "start_codon_class": "ATG" if row.starts_with_M else "non-ATG",
              "qc_flags": repr(fatal + advisory)}
        if fatal:
            qc.update(disposition="qc_excluded", upstream_nt_available="", full_kozak_strength="")
            qc_rows.append(qc)
            continue

        nt = seq["nt"]
        codon = nt[:3]
        tier = CODON_TIER.get(codon, "")
        ctx = kozak_context(seq["upstream"], nt)
        exons = seq["exons"]
        rec = {
            "orf_id": orf, "gene_name": row.gene_name, "smorf_type": row.smorf_type,
            "discovery_origin": row.discovery_origin,
            "chromosome": seq["chrom"], "strand": seq["strand"],
            "genomic_coordinates": row.genomic_coordinates,
            "genomic_start": (min(lo for lo, _ in exons) if seq["strand"] == "+"
                              else max(hi for _, hi in exons)),
            "n_exons": len(exons), "is_spliced": len(exons) > 1,
            "protein_length_aa": len(protein), "nt_length": len(nt),
            "start_codon": codon,
            "start_codon_class": "ATG" if codon == "ATG" else "non-ATG",
            "initiator_class": CODON_CLASS.get(tier, "not-a-documented-initiator"),
            "initiation_tier": tier, "initiation_tier_name": TIER_NAME.get(tier, ""),
            "kozak_window": ctx["kozak_window"],
            "upstream_nt_available": ctx["upstream_nt_available"],
        }
        for p in WINDOW_POSITIONS:
            name = ("plus%d" if p > 0 else "minus%d") % abs(p)
            rec["pos_" + name] = ctx[name]
        rec.update({
            # Immediately before the two scores, so the choice is unmissable in column order.
            "which_score_to_use": SCORE_GUIDANCE,
            "full_kozak_score": ctx["full_kozak_score"], "full_kozak_max": ctx["full_kozak_max"],
            "full_kozak_fraction": ctx["full_kozak_fraction"], "full_kozak_strength": ctx["full_kozak_strength"],
            "downstream_kozak_score": ctx["downstream_kozak_score"], "downstream_kozak_max": ctx["downstream_kozak_max"],
            "downstream_kozak_fraction": ctx["downstream_kozak_fraction"],
            "downstream_kozak_strength": ctx["downstream_kozak_strength"],
            "upstream_is_inframe_stop": ctx["upstream_is_inframe_stop"],
            "minus3_purine": ctx["minus3_purine"], "plus4_G": ctx["plus4_G"],
            "kozak_class_canonical": ctx["kozak_class_canonical"],
            "consensus_matches": ctx["consensus_matches"],
            "consensus_positions_testable": ctx["consensus_positions_testable"],
            "consensus_match_fraction": ctx["consensus_match_fraction"],
            # Advisory only - coordinates, never context. Empty for all but a handful of ORFs.
            "coordinate_flags": repr(advisory),
            "CLICK_UCSC": row.CLICK_UCSC,
        })
        rows.append(rec)
        qc.update(disposition="analysed", upstream_nt_available=ctx["upstream_nt_available"],
                  full_kozak_strength=ctx["full_kozak_strength"])
        qc_rows.append(qc)
    return pd.DataFrame(rows), pd.DataFrame(qc_rows)


# ======================================================================================
# Stage 4 - aggregate
# ======================================================================================

def summarise_group(g):
    n = len(g)
    return pd.Series({
        "n": n,
        "mean_full_kozak_fraction": round(g.full_kozak_fraction.mean(), 3),
        "median_full_kozak_fraction": round(g.full_kozak_fraction.median(), 3),
        "mean_downstream_kozak_fraction": round(g.downstream_kozak_fraction.mean(), 3),
        "pct_weighted_strong": round(100.0 * (g.full_kozak_strength == "strong").mean(), 1),
        "pct_weighted_adequate": round(100.0 * (g.full_kozak_strength == "adequate").mean(), 1),
        "pct_weighted_weak": round(100.0 * (g.full_kozak_strength == "weak").mean(), 1),
        "pct_upstream_is_inframe_stop": round(100.0 * g.upstream_is_inframe_stop.mean(), 1),
        "pct_minus3_purine": round(100.0 * g.minus3_purine.mean(), 1),
        "pct_plus4_G": round(100.0 * g.plus4_G.mean(), 1),
        "pct_canonical_strong": round(100.0 * (g.kozak_class_canonical == "strong").mean(), 1),
        "pct_canonical_weak": round(100.0 * (g.kozak_class_canonical == "weak").mean(), 1),
        "pct_full_upstream_window": round(100.0 * (g.upstream_nt_available >= 6).mean(), 1),
    })


def aggregate(d):
    by_codon = (d.groupby("start_codon").apply(summarise_group, include_groups=False)
                .reset_index())
    by_codon["n"] = by_codon.n.astype(int)
    by_codon["initiation_tier"] = by_codon.start_codon.map(lambda c: CODON_TIER.get(c, ""))
    by_codon["initiator_class"] = by_codon.initiation_tier.map(
        lambda t: CODON_CLASS.get(t, "not-a-documented-initiator"))
    by_codon = by_codon.sort_values("n", ascending=False).reset_index(drop=True)

    by_type = (d.groupby(["smorf_type", "start_codon_class"])
               .apply(summarise_group, include_groups=False).reset_index()
               .sort_values(["smorf_type", "start_codon_class"]).reset_index(drop=True))
    by_type["n"] = by_type.n.astype(int)
    return by_codon, by_type


def base_frequencies(d, min_n=20):
    """
    Base composition at each context position, per group - the PWM / sequence-logo data.

    Positions the transcript does not reach ('-') are excluded from that position's
    denominator, so each position's frequencies are over the ORFs that actually have it.

    `min_n` is deliberately low: the commonest non-ATG start has only a few dozen smORFs, so a
    higher bar would silently leave this table covering ATG alone.
    """
    groups = [("start_codon_class", v) for v in ("ATG", "non-ATG")]
    groups += [("start_codon", c) for c, n in d.start_codon.value_counts().items() if n >= min_n]
    groups += [("smorf_type", t) for t in sorted(d.smorf_type.dropna().unique())]

    out = []
    for kind, value in groups:
        g = d[d[kind] == value]
        for p in WINDOW_POSITIONS:
            col = "pos_" + (("plus%d" if p > 0 else "minus%d") % abs(p))
            counts = Counter(b for b in g[col] if b in "ACGT")
            total = sum(counts.values())
            for base in "ACGT":
                out.append({"group_kind": kind, "group": value, "n_orfs": len(g),
                            "position": p, "base": base, "count": counts.get(base, 0),
                            "frequency": round(counts.get(base, 0) / total, 4) if total else 0.0,
                            "n_with_position": total})
    return pd.DataFrame(out)


# ======================================================================================
# Stage 5 - figures
# ======================================================================================

# The start codon itself (+1..+3) is definitional, not context, and plotting it would only
# restate `start_codon`. The CSV keeps all twelve positions; the heatmap shows the nine that
# carry information.
CONTEXT_ONLY = [p for p in WINDOW_POSITIONS if p < 0 or p >= 4]


def figure_strength(d, figure_dir):
    """
    Context-strength composition by start-codon class, on all three readouts.

    The third panel is not optional. The first two both include -3, which the stop-to-stop
    definition fixes to T in every non-ATG ORF, so on their own they would read as a
    biological deficit that is not there.
    """
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    order = ["strong", "adequate", "weak"]
    colours = {"strong": "#1b7837", "adequate": "#c2a5cf", "weak": "#b2182b"}
    groups = ["ATG", "non-ATG"]

    for ax, (key, title) in zip(axes, [
            ("full_kozak_strength", "Weighted score\n(includes −3)"),
            ("kozak_class_canonical", "Canonical Kozak\n(−3 purine, +4 G)"),
            ("downstream_kozak_strength", "Downstream only (+4, +5, +6)\nthe comparable readout")]):
        bottom = np.zeros(len(groups))
        for level in order:
            vals = np.array([100.0 * (d[d.start_codon_class == g][key] == level).mean()
                             for g in groups])
            ax.bar(groups, vals, bottom=bottom, color=colours[level], label=level)
            for x, (v, b) in enumerate(zip(vals, bottom)):
                if v > 4:
                    ax.text(x, b + v / 2, "%.0f%%" % v, ha="center", va="center",
                            fontsize=9, color="white", fontweight="bold")
            bottom += vals
        ax.set_ylabel("% of smORFs")
        ax.set_title(title, fontsize=11)
        ax.set_ylim(0, 100)
    axes[0].legend(frameon=False, loc="lower right")
    axes[1].text(1, 3, "cannot be 'strong':\n−3 is the upstream stop", ha="center",
                 fontsize=8, color="white", style="italic")
    fig.suptitle("Initiation context of smORF start codons — "
                 "the first two panels are confounded, the third is not", fontweight="bold")
    fig.tight_layout()
    fig.savefig(os.path.join(figure_dir, "kozak_strength_by_class.png"), dpi=200)
    plt.close(fig)


def figure_upstream_artefact(d, comp, figure_dir):
    """
    The central result: the ATG/non-ATG gap lives entirely at -3, which the stop-to-stop
    definition fixes to T. Left, the two scores side by side; right, the counterfactual.
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    groups, colours = ["ATG", "non-ATG"], ["#2166ac", "#b2182b"]

    x = np.arange(2)
    for i, (g, c) in enumerate(zip(groups, colours)):
        sub = d[d.start_codon_class == g]
        axes[0].bar(x + (i - 0.5) * 0.36,
                    [sub.full_kozak_fraction.mean(), sub.downstream_kozak_fraction.mean()],
                    0.36, color=c, label=g)
    axes[0].set_xticks(x, ["full weighted\n(-4,-3,+4,+5,+6)", "downstream only\n(+4,+5,+6)"])
    axes[0].set_ylabel("mean context score")
    axes[0].set_title("The gap is confined to the upstream positions")
    axes[0].legend(frameon=False)

    labels = ["ATG\nas observed", "ATG rescored,\n−3 forced non-purine", "non-ATG\nas observed"]
    vals = [comp["mean_full_kozak_fraction_atg"], comp["counterfactual_atg_fraction"],
            comp["mean_full_kozak_fraction_non_atg"]]
    axes[1].bar(labels, vals, color=["#2166ac", "#7fb3d5", "#b2182b"])
    for i, v in enumerate(vals):
        axes[1].text(i, v + 0.008, "%.3f" % v, ha="center", fontsize=10, fontweight="bold")
    axes[1].set_ylabel("mean weighted context score")
    axes[1].set_title("−3 is fixed by the stop-to-stop definition")
    axes[1].tick_params(axis="x", labelsize=8)

    fig.suptitle("Why non-ATG smORFs appear to have weak Kozak context", fontweight="bold")
    fig.tight_layout()
    fig.savefig(os.path.join(figure_dir, "upstream_stop_artefact.png"), dpi=200)
    plt.close(fig)


def figure_by_codon(by_codon, figure_dir, min_n=20):
    """
    Per-codon context, ordered by evidence tier then abundance.

    Plotted on the DOWNSTREAM-ONLY score. Every non-ATG codon here belongs to a stop-to-stop
    ORF, so the full weighted score would rank them all against a -3 that is fixed to T; the
    3' positions are the only ones on which a codon can be compared with ATG.
    """
    d = by_codon[by_codon.n >= min_n].copy()
    if d.empty:
        log.warning("      no start codon reaches n >= %d; skipping per-codon figure", min_n)
        return
    d["rank"] = d.initiation_tier.map(lambda t: TIER_ORDER.index(t) if t in TIER_ORDER else 99)
    d = d.sort_values(["rank", "n"], ascending=[True, False])
    atg = d[d.start_codon == "ATG"]

    fig, ax = plt.subplots(figsize=(max(9, 0.34 * len(d)), 5))
    colour = {"ATG": "#2166ac", "near-cognate": "#4393c3",
              "non-near-cognate": "#92c5de", "not-a-documented-initiator": "#bdbdbd"}
    ax.bar(d.start_codon, d.mean_downstream_kozak_fraction,
           color=[colour[c] for c in d.initiator_class])
    for x, (v, n) in enumerate(zip(d.mean_downstream_kozak_fraction, d.n)):
        ax.text(x, v + 0.006, "%d" % n, ha="center", fontsize=6.5, rotation=90)
    if len(atg):
        atg_mean = float(atg.mean_downstream_kozak_fraction.iloc[0])
        ax.axhline(atg_mean, ls="--", lw=1.2, color="#2166ac")
        ax.text(len(d) - 0.5, atg_mean + 0.006, "ATG mean", ha="right", fontsize=8,
                color="#2166ac")
    ax.set_ylabel("mean downstream-only context score (+4, +5, +6)")
    ax.set_xlabel("start codon (bar label = number of smORFs)")
    ax.set_title("Downstream initiation context by start codon (n ≥ %d), "
                 "tier order" % min_n)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in colour.values()]
    ax.legend(handles, colour.keys(), frameon=False, fontsize=8, ncol=2)
    ax.tick_params(axis="x", labelrotation=90, labelsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(figure_dir, "kozak_by_start_codon.png"), dpi=200)
    plt.close(fig)


def figure_base_frequencies(freq, figure_dir):
    """Positional base composition, ATG vs non-ATG, with the difference between them."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    mats = {}
    for name in ("ATG", "non-ATG"):
        g = freq[(freq.group_kind == "start_codon_class") & (freq.group == name)]
        mats[name] = np.array([[g[(g.position == p) & (g.base == b)].frequency.iloc[0]
                                for p in CONTEXT_ONLY] for b in "ACGT"])
    # Full 0-1 scale: the non-ATG -3 row reaches 1.0 and must not be clipped, since that cell
    # is the result.
    for ax, name in zip(axes[:2], ("ATG", "non-ATG")):
        im = ax.imshow(mats[name], cmap="Blues", vmin=0, vmax=1.0, aspect="auto")
        ax.set_title("%s starts" % name)
        fig.colorbar(im, ax=ax, fraction=0.046, label="frequency")
    diff = mats["non-ATG"] - mats["ATG"]
    im = axes[2].imshow(diff, cmap="RdBu_r", vmin=-1.0, vmax=1.0, aspect="auto")
    axes[2].set_title("non-ATG − ATG")
    fig.colorbar(im, ax=axes[2], fraction=0.046, label="Δ frequency")
    gap = CONTEXT_ONLY.index(4) - 0.5
    for ax in axes:
        ax.set_xticks(range(len(CONTEXT_ONLY)), ["%+d" % p for p in CONTEXT_ONLY], fontsize=8)
        ax.set_yticks(range(4), list("ACGT"))
        ax.set_xlabel("position relative to start codon")
        ax.axvline(gap, color="black", lw=2)   # the start codon (+1..+3) sits here, not shown
    fig.suptitle("Base composition flanking the start codon "
                 "(the codon itself, +1..+3, is omitted)", fontweight="bold")
    fig.tight_layout()
    fig.savefig(os.path.join(figure_dir, "position_base_frequencies.png"), dpi=200)
    plt.close(fig)


def figure_by_smorf_type(by_type, figure_dir):
    """
    Mean context per smORF type, ATG and non-ATG side by side.

    Downstream-only again: this figure compares the two classes within each type, which the
    full weighted score cannot do.
    """
    piv = by_type.pivot(index="smorf_type", columns="start_codon_class",
                        values="mean_downstream_kozak_fraction")
    cnt = by_type.pivot(index="smorf_type", columns="start_codon_class", values="n").fillna(0)
    order = cnt.sum(axis=1).sort_values(ascending=False).index
    piv, cnt = piv.loc[order], cnt.loc[order]   # both, or the bar labels desynchronise
    x = np.arange(len(piv))
    fig, ax = plt.subplots(figsize=(max(9, 0.62 * len(piv)), 5))
    for i, (col, colour) in enumerate([("ATG", "#2166ac"), ("non-ATG", "#b2182b")]):
        if col not in piv:
            continue
        pos = x + (i - 0.5) * 0.4
        ax.bar(pos, piv[col].fillna(0), 0.4, label=col, color=colour)
        # n on every bar: several types have only a handful of non-ATG members, and without
        # the count their bars read as confidently as the largest type's.
        for xi, (v, k) in zip(pos, zip(piv[col].fillna(0), cnt[col])):
            ax.text(xi, v + 0.004, "%d" % k, ha="center", fontsize=6.5, rotation=90,
                    color=colour)
    ax.set_xticks(x, piv.index, rotation=45, ha="right")
    ax.set_ylabel("mean downstream-only context score (+4, +5, +6)")
    ax.set_title("Downstream initiation context by smORF type")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(figure_dir, "kozak_by_smorf_type.png"), dpi=200)
    plt.close(fig)


# ======================================================================================
# Statistics and internal consistency check
# ======================================================================================

def compare_classes(d):
    """
    ATG vs non-ATG on each readout. Sample sizes are large, so effect size is the point.

    Both scores are compared, and the difference between the two comparisons is the result:
    the full weighted score includes -3, which the stop-to-stop definition forces to T in every
    non-ATG ORF, while the downstream-only score does not. `counterfactual_atg_fraction`
    isolates exactly that - the ATG starts rescored as if their -3 were also non-purine. If it
    lands on the observed non-ATG mean, the apparent context deficit is the definition, not
    biology.
    """
    atg = d[d.start_codon_class == "ATG"]
    non = d[d.start_codon_class == "non-ATG"]
    out = {"n_atg": len(atg), "n_non_atg": len(non)}
    # Keys carry the same prefixes as the columns, so a statistic can never be reported against
    # the wrong score by accident.
    for col in ("full_kozak_fraction", "downstream_kozak_fraction"):
        u, p_u = stats.mannwhitneyu(atg[col], non[col], alternative="two-sided")
        stem = col.rsplit("_", 1)[0]          # full_kozak / downstream_kozak
        out.update({
            "mean_%s_atg" % col: round(atg[col].mean(), 3),
            "mean_%s_non_atg" % col: round(non[col].mean(), 3),
            "%s_mannwhitney_p" % stem: float(p_u),
            # Rank-biserial correlation: the readable effect size for Mann-Whitney.
            "%s_rank_biserial_r" % stem: round(2 * u / (len(atg) * len(non)) - 1, 3)})

    counter_score = (2 * (atg.pos_plus4 == "G") + (atg.pos_minus4.isin(["C", "A"]))
                     + (atg.pos_plus5 == "A") + (atg.pos_plus6 == "T"))
    counter_max = (2 * (atg.pos_plus4 != "-") + 2 * (atg.pos_minus3 != "-")
                   + (atg.pos_minus4 != "-") + (atg.pos_plus5 != "-") + (atg.pos_plus6 != "-"))
    out["counterfactual_atg_fraction"] = round((counter_score / counter_max).mean(), 3)
    # Conditional on the window existing at all: an ORF at a transcript 5' end has no -3..-1 to
    # inspect, and counting it as "not a stop" would understate a rate that is otherwise 100%.
    for label, g in (("atg", atg), ("non_atg", non)):
        have = g[g.upstream_nt_available >= 3]
        out["pct_%s_upstream_is_stop" % label] = round(100.0 * have.upstream_is_inframe_stop.mean(), 1)
        out["n_%s_with_upstream" % label] = len(have)

    for label, col in [("minus3_purine", "minus3_purine"), ("plus4_G", "plus4_G")]:
        table = np.array([[int(atg[col].sum()), int((~atg[col]).sum())],
                          [int(non[col].sum()), int((~non[col]).sum())]])
        chi2, p, _, _ = stats.chi2_contingency(table)
        out["pct_%s_atg" % label] = round(100.0 * atg[col].mean(), 1)
        out["pct_%s_non_atg" % label] = round(100.0 * non[col].mean(), 1)
        out["chi2_%s_p" % label] = float(p)
    return out


def cross_check_against_initiation(d, candidates_csv):
    """
    Internal consistency check: the non-ATG smORFs are scored at both steps of this pipeline,
    by separate code, and the two must agree exactly. Skipped when the earlier step has not
    been run.
    """
    if not os.path.exists(candidates_csv):
        log.warning("      initiation-step output not found (%s); cross-check skipped",
                    os.path.basename(candidates_csv))
        return None
    old = pd.read_csv(candidates_csv)
    old = old[old.is_annotated_start][["orf_id", "codon", "context_fraction",
                                       "context_strength", "context_score", "context_max"]]
    m = old.merge(d[["orf_id", "start_codon", "full_kozak_fraction", "full_kozak_strength",
                     "full_kozak_score", "full_kozak_max"]], on="orf_id", how="inner")
    mismatched = m[(m.codon != m.start_codon) | (m.context_fraction != m.full_kozak_fraction) |
                   (m.context_strength != m.full_kozak_strength) | (m.context_score != m.full_kozak_score) |
                   (m.context_max != m.full_kozak_max)]
    log.info("      cross-check vs initiation step: %d ORFs compared, %d mismatched",
             len(m), len(mismatched))
    assert len(mismatched) == 0, "context scores diverge from the initiation pipeline"
    return {"n_compared": len(m), "n_mismatched": len(mismatched)}


# ======================================================================================
# Report
# ======================================================================================

def write_report(paths, d, by_codon, by_type, freq, qc, comp, reg, totals, log_path):
    n = len(d)
    atg, non = d[d.start_codon_class == "ATG"], d[d.start_codon_class == "non-ATG"]
    top_codons = by_codon[by_codon.n >= 20].head(25)
    strength = pd.DataFrame({
        "context": ["strong", "adequate", "weak"],
        "ATG %": [round(100.0 * (atg.full_kozak_strength == s).mean(), 1)
                  for s in ("strong", "adequate", "weak")],
        "non-ATG %": [round(100.0 * (non.full_kozak_strength == s).mean(), 1)
                      for s in ("strong", "adequate", "weak")]})
    canon = pd.DataFrame({
        "canonical Kozak": ["strong", "adequate", "weak"],
        "ATG %": [round(100.0 * (atg.kozak_class_canonical == s).mean(), 1)
                  for s in ("strong", "adequate", "weak")],
        "non-ATG %": [round(100.0 * (non.kozak_class_canonical == s).mean(), 1)
                      for s in ("strong", "adequate", "weak")]})
    # Counts belong next to the means: several types have few non-ATG members, and a difference
    # in the third decimal there is noise.
    piv = by_type.pivot(index="smorf_type", columns="start_codon_class",
                        values=["mean_downstream_kozak_fraction", "n"])
    type_tbl = pd.DataFrame({
        "smorf_type": piv.index,
        "ATG score": piv.get(("mean_downstream_kozak_fraction", "ATG")),
        "ATG n": piv.get(("n", "ATG")),
        "non-ATG score": piv.get(("mean_downstream_kozak_fraction", "non-ATG")),
        "non-ATG n": piv.get(("n", "non-ATG")),
    }).reset_index(drop=True)
    type_tbl = type_tbl.sort_values("ATG n", ascending=False).reset_index(drop=True)
    reg_line = ("Context scores agree exactly across all %d ORFs also scored by the initiation "
                "step." % reg["n_compared"]) if reg else "Initiation step not run; not checked."
    n_advisory = int((d.coordinate_flags != "[]").sum())

    text = f"""# Kozak Context of Every smORF Start Codon

Generated {datetime.now().strftime("%Y-%m-%d %H:%M")} by `kozak_pipeline.py`
from `{os.path.basename(paths['master'])}`.
Positions follow the usual convention: **+1..+3 are the start codon**, +4 is the first
nucleotide after it, −1 the one before it, and there is no position 0.

## 1. Question

`initiation_pipeline.py` asks where translation could start in the {len(non):,}
smORFs that do not begin with ATG. This analysis asks a narrower question of the whole atlas:
**how good is each smORF's annotated initiation context, and how do the non-ATG starts compare
with the ATG starts around them?**

smORFs here are called with GTFtoFASTA (Martinez et al., 2019): an ORF containing an in-frame ATG
takes that ATG as its start, and only an ORF without one is reported stop-codon-to-stop-codon.
The {len(atg):,} ATG smORFs are therefore the internal reference — same transcripts, same ORF
caller, same sequence files, but a start that is not a boundary of the calling method.

## 2. Cohort

| quantity | value |
|---|---:|
| rows in the master table | {totals['input_rows']:,} |
| unique ids | {totals['input_ids']:,} |
| ... not smORF ORF ids (UniProt/reference; no ORF model) | {totals['non_smorf_ids']:,} |
| smORFs | {totals['cohort']:,} |
| ... excluded by QC | {totals['qc_failed']:,} |
| **scored** | **{n:,}** |
| ... ATG | {len(atg):,} |
| ... non-ATG | {len(non):,} |
| with the full 6 nt upstream window | {int((d.upstream_nt_available >= 6).sum()):,} ({100.0 * (d.upstream_nt_available >= 6).mean():.1f}%) |

The cohort is whatever the input table contains, so these counts move with `--master`: the
master table shipped with this repository is the released atlas, while the upstream annotation
CSV holds the full pre-release discovery set.

## 3. Three independent readouts

**Weighted score** — five positions, weighted by reported effect size (+4 G and −3 purine at
weight 2; −4 C/A, +5 A, +6 T at weight 1), divided by the weight actually available. Positions
the transcript does not reach drop out of the denominator, so ORFs near a transcript end stay
comparable. `strong` ≥ {STRONG_CONTEXT}, `adequate` ≥ {ADEQUATE_CONTEXT}, else `weak`.

**Downstream-only score** — the same scheme restricted to +4, +5 and +6. Section 4 is the reason
it exists: the upstream positions are not free to vary in a stop-to-stop ORF, so only the 3′
positions support a comparison between ATG and non-ATG starts.

**Canonical Kozak** — Kozak's own two-position rule: `strong` = purine at −3 **and** G at +4;
`adequate` = exactly one; `weak` = neither. It depends on no weighting choice of ours, which is
why it is reported alongside rather than instead.

{md_table(strength)}

{md_table(canon)}

The non-ATG column of both tables looks alarming and should not be read yet — **section 4 shows
that most of it is an artefact of how these ORFs were defined.**

## 4. ATG vs non-ATG — and why −3 must be thrown out

| | ATG | non-ATG |
|---|---:|---:|
| mean **full** weighted score | {comp['mean_full_kozak_fraction_atg']} | {comp['mean_full_kozak_fraction_non_atg']} |
| mean **downstream-only** score (+4/+5/+6) | {comp['mean_downstream_kozak_fraction_atg']} | {comp['mean_downstream_kozak_fraction_non_atg']} |
| −3 purine | {comp['pct_minus3_purine_atg']}% | {comp['pct_minus3_purine_non_atg']}% |
| +4 G | {comp['pct_plus4_G_atg']}% | {comp['pct_plus4_G_non_atg']}% |
| −3..−1 is an in-frame stop codon¹ | {comp['pct_atg_upstream_is_stop']}% | **{comp['pct_non_atg_upstream_is_stop']}%** |

¹ Of the ORFs that have those three nucleotides at all: {comp['n_atg_with_upstream']:,} ATG and
{comp['n_non_atg_with_upstream']:,} non-ATG. The rest sit at a transcript 5′ end.

**The −3 column is not a measurement.** These non-ATG smORFs are defined stop-codon-to-stop-codon,
so the three nucleotides immediately 5′ of the start *are* the upstream in-frame stop — TAA, TAG
or TGA. Every one of them therefore begins with T at −3, which is never a purine. The
0% above is arithmetic, not biology, and −3 carries weight 2 of 7 in the weighted score.

Two things follow, and they point the same way:

- On the **downstream-only** score, which no definition constrains, the two classes are
  {'indistinguishable' if abs(comp['downstream_kozak_rank_biserial_r']) < 0.1 else 'still separated'}:
  {comp['mean_downstream_kozak_fraction_atg']} vs {comp['mean_downstream_kozak_fraction_non_atg']},
  rank-biserial r = {comp['downstream_kozak_rank_biserial_r']}
  (compared with r = {comp['full_kozak_rank_biserial_r']} on the full score).
- Rescoring the **ATG** starts as if their −3 were also non-purine gives
  **{comp['counterfactual_atg_fraction']}** — against the non-ATG group's observed
  {comp['mean_full_kozak_fraction_non_atg']}. The gap essentially closes.

So the apparent context deficit of non-ATG smORFs is largely an artefact of how the ORFs were
defined. Their 3′ context is ordinary. This does not rescue the annotated N-termini — the point
of the initiation step stands, that a stop-to-stop start is a definitional boundary rather
than a measured initiation site — but it does mean **weak weighted context is not independent
evidence for that claim**, because the definition produces it directly.

Mann–Whitney p-values: {fmt_p(comp['full_kozak_mannwhitney_p'])} (full) and
{fmt_p(comp['downstream_kozak_mannwhitney_p'])} (downstream-only); χ² p =
{fmt_p(comp['chi2_minus3_purine_p'])} (−3 purine) and {fmt_p(comp['chi2_plus4_G_p'])} (+4 G).
With {n:,} ORFs almost any difference reaches significance, so read the effect size, not the
p-value.

## 5. By start codon

{md_table(top_codons, ['start_codon', 'initiation_tier', 'initiator_class', 'n',
                       'mean_downstream_kozak_fraction', 'mean_full_kozak_fraction',
                       'pct_plus4_G', 'pct_minus3_purine'])}

Every non-ATG codon here is a stop-to-stop start, so read `mean_downstream_kozak_fraction`;
`mean_full_kozak_fraction` is shown beside it only to make the size of the −3 penalty visible.
Full table, including the strength breakdowns, in `kozak_by_start_codon.csv`.

## 6. By smORF type

{md_table(type_tbl)}

Downstream-only score, so the two classes are comparable within each type. Counts and every
other statistic are in `kozak_by_smorf_type.csv`.

## 7. Quality control

{reg_line}

Checks are split by what they actually threaten. **Fatal** — no sequence, a sequence that does
not translate to the annotated protein, a length that is not a multiple of three, or an ORF
occurring more than once in its transcript (where the string search cannot tell which copy the
upstream window came from, so −1..−6 is unreliable). {totals['qc_failed']:,} smORFs are excluded
on these grounds, listed with their flags in `qc_audit.csv`.

**Advisory** — the GTF exon blocks disagree with the sequence length, or the strand in the ORF
id disagrees with the GTF. Exon blocks feed `genomic_start` and `n_exons` and nothing else, so
these ORFs keep every context measure and carry the flag in `coordinate_flags` instead:
{n_advisory:,} smORFs, typically a spurious short leading exon block in the GTF.

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
| `kozak_context.csv` | {n:,} | one row per smORF — the primary table |
| `kozak_by_start_codon.csv` | {len(by_codon):,} | one row per start codon |
| `kozak_by_smorf_type.csv` | {len(by_type):,} | one row per (smORF type × ATG/non-ATG) |
| `position_base_frequencies.csv` | {len(freq):,} | base composition at −6..+6, the PWM data |
| `qc_audit.csv` | {len(qc):,} | every smORF id and its disposition |
| `COLUMNS.md` | – | data dictionary: every column of every output, grouped by family |
| `figures/*.png` | 5 | strength, the −3 artefact, per-codon, base composition, per-type |
| `logs/{os.path.basename(log_path)}` | – | full run log |

Sequence input: `{os.path.basename(paths['sequences'])}`, built by `build_sequence_context.py`
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
- The upstream sequence comes from the assembled transcript. {int((d.upstream_nt_available < 6).sum()):,}
  smORFs begin within 6 nt of a transcript 5′ end, so their window is truncated and the missing
  positions drop out of the denominator — a limit of the assembly rather than of the ORF.
- ATG smORFs are a reference, not a gold standard: their start codons are called by the same ORF
  finder and have not been confirmed by N-terminomics either.
"""
    with open(paths["report"], "w") as fh:
        fh.write(text)


# ======================================================================================
# Stage 5b - data dictionary
# ======================================================================================
#
# Grouped by family, because the families are what gets confused: two scores sit side by side
# and only one of them is safe for an ATG-vs-non-ATG comparison. `write_columns_md` asserts the
# tables describe every column actually written, so the dictionary cannot drift.

CONTEXT_FAMILIES = [
    ("Identity", "One row per scored smORF.", [
        ("orf_id", "ORF identifier, as `gene_id` in the master annotation"),
        ("gene_name", "host gene symbol; includes the placeholders `Unnamed`, `Intergenic`, `Unknown`"),
        ("smorf_type", "smORF class (uORF, dORF, Iso, lncRNA, ...)"),
        ("discovery_origin", "how the ORF entered the atlas; empty when the input table does not carry this column"),
        ("protein_length_aa", "annotated protein length"),
        ("nt_length", "spliced CDS length; always 3x the protein length"),
    ]),
    ("Coordinates", "Spliced-aware, from the GTF exon blocks. Advisory flags affect only these.", [
        ("chromosome", "chromosome"), ("strand", "genomic strand"),
        ("genomic_coordinates", "the ORF's genomic min-max; spans introns when spliced"),
        ("genomic_start", "genomic coordinate of the start codon's first base"),
        ("n_exons", "exon blocks the ORF spans"), ("is_spliced", "`n_exons` > 1"),
        ("CLICK_UCSC", "UCSC browser link, carried through from the master annotation"),
        ("coordinate_flags", "**list** - advisory QC. Affects `genomic_start` / `n_exons` and nothing else; the context is untouched. Empty for all but a handful of smORFs, typically a spurious short leading exon block in the GTF"),
    ]),
    ("Start codon", "What is there, and how well documented it is as an initiator.", [
        ("start_codon", "the ORF's first codon"),
        ("start_codon_class", "`ATG` or `non-ATG` - the grouping variable for every comparison"),
        ("initiation_tier", "tier `0`-`3c` from `non_aug_initiation_codons.md`; blank when the codon is in no tier. **Lower is stronger**"),
        ("initiation_tier_name", "that tier written out"),
        ("initiator_class", "`ATG` / `near-cognate` / `non-near-cognate` / `not-a-documented-initiator`"),
    ]),
    ("Raw context", "The sequence itself. Everything below is computed from these, so reweight "
                    "freely if the schemes here do not suit.", [
        ("kozak_window", "the -6..+6 window as a string, `gccRccATGG`-style: upstream lowercase, codon and +4..+6 uppercase"),
        ("upstream_nt_available", "0-6; how much 5' window the transcript actually provides. Below 6 the ORF starts near a transcript end"),
        ("pos_minus6", "base at -6"), ("pos_minus5", "base at -5"), ("pos_minus4", "base at -4"),
        ("pos_minus3", "base at -3"), ("pos_minus2", "base at -2"), ("pos_minus1", "base at -1"),
        ("pos_plus1", "base at +1 (start codon)"), ("pos_plus2", "base at +2"),
        ("pos_plus3", "base at +3"), ("pos_plus4", "base at +4"), ("pos_plus5", "base at +5"),
        ("pos_plus6", "base at +6. Any position the transcript does not reach is `-`"),
    ]),
    ("SCORE 1 of 2 - full", "All five weighted positions, INCLUDING -3. **Valid only within one "
                            "start-codon class.** A non-ATG smORF is defined stop-to-stop, so its "
                            "-3 is the first base of the upstream in-frame stop and is therefore "
                            "always T - unfavourable by construction, at weight 2 of 7. Comparing "
                            "classes on this score measures the ORF definition, not biology.", [
        ("which_score_to_use", "constant; restates the choice below in every row, so it cannot be missed by someone opening the CSV without these docs"),
        ("full_kozak_score", "weight achieved: +4 G and -3 purine 2 each; -4 C/A, +5 A, +6 T 1 each"),
        ("full_kozak_max", "weight available; positions outside the transcript drop out"),
        ("full_kozak_fraction", "`full_kozak_score` / `full_kozak_max`"),
        ("full_kozak_strength", "`strong` >= %.2f, `adequate` >= %.2f, else `weak`" % (STRONG_CONTEXT, ADEQUATE_CONTEXT)),
    ]),
    ("SCORE 2 of 2 - downstream", "The same scheme restricted to +4, +5 and +6, where the "
                                  "stop-to-stop definition constrains nothing. **This is the "
                                  "ATG-vs-non-ATG readout**, and the rank-biserial effect sizes "
                                  "for both scores are in the report's section 4.", [
        ("downstream_kozak_score", "weight achieved over +4/+5/+6"),
        ("downstream_kozak_max", "weight available"),
        ("downstream_kozak_fraction", "the comparable score"),
        ("downstream_kozak_strength", "same thresholds as above"),
    ]),
    ("Other readouts", "Independent of the weighting scheme above.", [
        ("minus3_purine", "-3 is A or G"),
        ("plus4_G", "+4 is G"),
        ("kozak_class_canonical", "Kozak's own rule (1986, 1987): `strong` = both of the above, `adequate` = exactly one, `weak` = neither. Can never be `strong` for a stop-to-stop non-ATG smORF, for the same reason as above"),
        ("upstream_is_inframe_stop", "-3..-1 is TAA/TAG/TGA - the per-row flag for the artefact. True for every non-ATG smORF that has the window, and for a small percentage of ATG ones"),
        ("consensus_matches", "positions agreeing with the full `gccRccATGG` consensus"),
        ("consensus_positions_testable", "positions the transcript provides"),
        ("consensus_match_fraction", "the ratio of the two"),
    ]),
]

_AGG_COLS = [
    ("n", "smORFs in the group"),
    ("mean_full_kozak_fraction", "mean full score - **within-class only**, see the primary table"),
    ("median_full_kozak_fraction", "its median"),
    ("mean_downstream_kozak_fraction", "mean downstream score - **the comparable one**"),
    ("pct_weighted_strong", "%% with `full_kozak_strength` == strong (>= %.2f)" % STRONG_CONTEXT),
    ("pct_weighted_adequate", "%% adequate"), ("pct_weighted_weak", "%% weak"),
    ("pct_upstream_is_inframe_stop", "%% whose -3..-1 is an in-frame stop codon"),
    ("pct_minus3_purine", "%% with a purine at -3"), ("pct_plus4_G", "%% with G at +4"),
    ("pct_canonical_strong", "%% strong under Kozak's own two-position rule"),
    ("pct_canonical_weak", "%% weak under it"),
    ("pct_full_upstream_window", "%% with all 6 upstream positions available"),
]

BY_CODON_FAMILIES = [
    ("Group", "One row per distinct start codon.", [
        ("start_codon", "the codon"),
        ("initiation_tier", "its tier; blank when in no tier"),
        ("initiator_class", "`ATG` / `near-cognate` / `non-near-cognate` / `not-a-documented-initiator`"),
    ]),
    ("Statistics", "Every non-ATG codon here is a stop-to-stop start, so read "
                   "`mean_downstream_kozak_fraction`. The full score is present only to make the "
                   "size of the -3 penalty visible.", _AGG_COLS),
]
BY_TYPE_FAMILIES = [
    ("Group", "One row per (smORF type x start-codon class). Because the class is a grouping "
              "column here, comparing the two rows of a type is a cross-class comparison - use "
              "the downstream score.", [
        ("smorf_type", "smORF class"), ("start_codon_class", "`ATG` or `non-ATG`"),
    ]),
    ("Statistics", "As above.", _AGG_COLS),
]
FREQ_FAMILIES = [
    ("PWM / sequence-logo data", "Base composition at each context position, per group. "
                                 "Positions the transcript does not reach are excluded from that "
                                 "position's denominator, so `n_with_position` is the true n for "
                                 "the row - not `n_orfs`.", [
        ("group_kind", "`start_codon_class`, `start_codon` (n >= 20 only) or `smorf_type`"),
        ("group", "the group's value"),
        ("n_orfs", "smORFs in the group overall"),
        ("position", "-6..+6; there is no position 0"),
        ("base", "A, C, G or T"),
        ("count", "smORFs with that base at that position"),
        ("frequency", "`count` / `n_with_position`"),
        ("n_with_position", "smORFs whose transcript reaches this position - the denominator"),
    ]),
]
QC_FAMILIES = [
    ("Disposition", "Every id matching the smORF ORF-id grammar and what happened to it. "
                    "Checks are split by what they threaten: fatal flags mean the context itself "
                    "cannot be trusted and the ORF is dropped; coordinate-only problems are "
                    "advisory and appear in `coordinate_flags` in the primary table instead.", [
        ("orf_id", "ORF identifier"), ("smorf_type", "smORF class"),
        ("start_codon_class", "`ATG` or `non-ATG`"),
        ("qc_flags", "**list** - fatal flags; empty when clean"),
        ("disposition", "`analysed` or `qc_excluded`"),
        ("upstream_nt_available", "as in the primary table; blank when excluded"),
        ("full_kozak_strength", "as in the primary table; blank when excluded"),
    ]),
]


def write_columns_md(paths, d, by_codon, by_type, freq, qc):
    """Data dictionary, generated from the frames so it cannot drift from the actual columns."""
    files = [("kozak_context.csv", d, CONTEXT_FAMILIES,
              "One row per scored smORF - the primary table. **It carries two scores and only "
              "one of them is safe for comparing ATG with non-ATG**; see the two score families "
              "below."),
             ("kozak_by_start_codon.csv", by_codon, BY_CODON_FAMILIES, "One row per start codon."),
             ("kozak_by_smorf_type.csv", by_type, BY_TYPE_FAMILIES,
              "One row per (smORF type x start-codon class)."),
             ("position_base_frequencies.csv", freq, FREQ_FAMILIES,
              "Base composition at -6..+6, the PWM / sequence-logo data."),
             ("qc_audit.csv", qc, QC_FAMILIES, "The audit trail.")]

    parts = ["""# Data Dictionary

Generated %s by `kozak_pipeline.py`, from the output frames themselves.

**The one thing to know:** there are two Kozak scores, and neither is called plain
`kozak_fraction` on purpose - the prefix forces the choice at the point of use.

| | `full_kozak_*` | `downstream_kozak_*` |
|---|---|---|
| positions | -4, -3, +4, +5, +6 | +4, +5, +6 only |
| compare ATG with non-ATG | **no** | **yes** |
| compare within one class | yes | yes |

A non-ATG smORF here is defined stop-codon-to-stop-codon, so the three nucleotides at -3, -2, -1
*are* the upstream in-frame stop - TAA, TAG or TGA. All three begin with T, so -3 is
unfavourable in 100%% of them by construction, at weight 2 of 7. The full score therefore
measures the ORF definition rather than biology whenever the classes are compared. The
downstream score is free of that.

Tier codes come from `non_aug_initiation_codons.md`. **Lower is stronger:**

| tier | codons | meaning |
|---|---|---|
%s

A blank tier means the codon has no documented initiation competence in mammals.
""" % (datetime.now().strftime("%Y-%m-%d %H:%M"),
       "\n".join("| %s | %s | %s |" % (t, " ".join(TIERS[t]), TIER_NAME[t]) for t in TIER_ORDER))]

    for name, frame, families, note in files:
        described = [c for _, _, cols in families for c, _ in cols]
        missing = [c for c in frame.columns if c not in described]
        extra = [c for c in described if c not in frame.columns]
        assert not missing, "%s: columns with no description: %s" % (name, missing)
        assert not extra, "%s: described columns that do not exist: %s" % (name, extra)

        parts.append("\n---\n\n## `%s`\n\n%s rows, %d columns. %s\n"
                     % (name, format(len(frame), ","), len(frame.columns), note))
        for family, fnote, cols in families:
            parts.append("### %s\n\n%s\n\n| column | meaning |\n|---|---|\n%s\n"
                         % (family, fnote,
                            "\n".join("| `%s` | %s |" % (c, dsc) for c, dsc in cols)))

    parts.append("""
---

## Worked example

```python
import pandas as pd
d = pd.read_csv("kozak_context.csv", low_memory=False)

# RIGHT - the comparable score
d.groupby("start_codon_class").downstream_kozak_fraction.mean()   # %s vs %s

# WRONG - this measures the stop-to-stop definition, not initiation context
d.groupby("start_codon_class").full_kozak_fraction.mean()         # %s vs %s

# fine: the full score within one class
d[d.start_codon_class == "ATG"].groupby("smorf_type").full_kozak_fraction.mean()
```
""" % tuple(round(d[d.start_codon_class == c][col].mean(), 3)
            for col in ("downstream_kozak_fraction", "full_kozak_fraction")
            for c in ("ATG", "non-ATG")))

    with open(paths["columns"], "w") as fh:
        fh.write("\n".join(parts))


# ======================================================================================
# Main
# ======================================================================================

def resolve_paths(args):
    outdir = args.outdir or os.path.join(common.DEFAULT_RESULTS_ROOT, "kozak")
    return {
        "master": common.resolve(args.master, common.ENV_MASTER, common.DEFAULT_MASTER),
        "sequences": common.resolve(args.sequences, common.ENV_SEQUENCES,
                                    common.DEFAULT_SEQUENCE_CONTEXT),
        "outdir": outdir,
        "figures": os.path.join(outdir, "figures"),
        "logs": os.path.join(outdir, "logs"),
        "context": os.path.join(outdir, "kozak_context.csv"),
        "by_codon": os.path.join(outdir, "kozak_by_start_codon.csv"),
        "by_type": os.path.join(outdir, "kozak_by_smorf_type.csv"),
        "freq": os.path.join(outdir, "position_base_frequencies.csv"),
        "qc": os.path.join(outdir, "qc_audit.csv"),
        "columns": os.path.join(outdir, "COLUMNS.md"),
        "report": os.path.join(outdir, "report.md"),
        # The annotated-start rows of the initiation step, which this
        # step must reproduce exactly for the non-ATG ORFs the two cohorts share.
        "initiation_candidates": (args.initiation_candidates or
                                  os.path.join(common.DEFAULT_RESULTS_ROOT, "initiation",
                                               "initiation_candidates.csv")),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    common.add_common_args(parser)
    parser.add_argument("--initiation-candidates", default=None,
                        help="initiation_candidates.csv from initiation_pipeline.py, used as "
                             "cross-checked against; skipped when absent")
    args = parser.parse_args(argv)
    paths = resolve_paths(args)

    log_path = common.setup_logging("kozak", paths["logs"], [paths["outdir"], paths["figures"]])
    log.info("=== kozak pipeline ===")
    log.info("      master:    %s", paths["master"])
    log.info("      sequences: %s", paths["sequences"])
    log.info("      outdir:    %s", paths["outdir"])

    cohort, totals = load_cohort(paths["master"])
    log.info("[2/5] sequence context")
    seqs = common.load_sequence_context(paths["sequences"], set(cohort.gene_id), log)

    log.info("[3/5] scoring annotated start codons")
    d, qc = score_rows(cohort, seqs)
    counts = Counter(qc.disposition)
    totals.update(cohort=len(cohort), qc_failed=counts["qc_excluded"])
    assert counts["analysed"] + counts["qc_excluded"] == len(cohort), "accounting does not close"
    assert (d.start_codon == "ATG").sum() == (d.start_codon_class == "ATG").sum()
    log.info("      scored %d smORFs (%d ATG, %d non-ATG); %d excluded by QC",
             len(d), int((d.start_codon_class == "ATG").sum()),
             int((d.start_codon_class == "non-ATG").sum()), counts["qc_excluded"])
    if counts["qc_excluded"]:
        log.warning("      fatal QC flags: %s", dict(Counter(
            f for fl in qc[qc.disposition == "qc_excluded"].qc_flags for f in ast.literal_eval(fl))))
    log.info("      advisory (coordinates only): %s", dict(Counter(
        f for fl in d.coordinate_flags for f in ast.literal_eval(fl))))

    log.info("[4/5] aggregating")
    by_codon, by_type = aggregate(d)
    freq = base_frequencies(d)
    comp = compare_classes(d)
    reg = cross_check_against_initiation(d, paths["initiation_candidates"])

    log.info("[5/5] writing outputs")
    d.to_csv(paths["context"], index=False)
    by_codon.to_csv(paths["by_codon"], index=False)
    by_type.to_csv(paths["by_type"], index=False)
    freq.to_csv(paths["freq"], index=False)
    qc.to_csv(paths["qc"], index=False)
    write_columns_md(paths, d, by_codon, by_type, freq, qc)
    figure_strength(d, paths["figures"])
    figure_upstream_artefact(d, comp, paths["figures"])
    figure_by_codon(by_codon, paths["figures"])
    figure_base_frequencies(freq, paths["figures"])
    figure_by_smorf_type(by_type, paths["figures"])
    write_report(paths, d, by_codon, by_type, freq, qc, comp, reg, totals, log_path)

    log.info("      weighted strength: %s", dict(Counter(d.full_kozak_strength)))
    log.info("      canonical Kozak:   %s", dict(Counter(d.kozak_class_canonical)))
    log.info("      full score   ATG %.3f vs non-ATG %.3f (r = %s)",
             comp["mean_full_kozak_fraction_atg"], comp["mean_full_kozak_fraction_non_atg"], comp["full_kozak_rank_biserial_r"])
    log.info("      downstream   ATG %.3f vs non-ATG %.3f (r = %s)",
             comp["mean_downstream_kozak_fraction_atg"], comp["mean_downstream_kozak_fraction_non_atg"],
             comp["downstream_kozak_rank_biserial_r"])
    log.info("      ATG rescored with -3 forced non-purine: %.3f  (non-ATG observed %.3f)",
             comp["counterfactual_atg_fraction"], comp["mean_full_kozak_fraction_non_atg"])
    log.info("done. log: %s", log_path)
    return d, by_codon, by_type, freq, qc


if __name__ == "__main__":
    main()
