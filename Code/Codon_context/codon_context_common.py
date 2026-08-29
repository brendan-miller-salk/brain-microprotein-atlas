#!/usr/bin/env python
"""
Shared infrastructure for the two analysis steps of this module's pipeline.

`kozak_pipeline.py` and `initiation_pipeline.py` were written against one machine's absolute
paths and one specific annotation CSV. Everything that made them machine-specific lives here
instead, so both scripts run unchanged against either input:

    * the master table shipped with this repository (`Code/data/microprotein_master.csv`), or
    * the full upstream annotation CSV, via `--master` / `$MICROPROTEIN_MASTER`.

Three things are generalised:

`resolve_*`      every path is derived from the repository root, and every one of them can be
                 overridden on the command line or through an environment variable.

`load_cohort`    the two source tables do not carry the same columns. Columns are requested by
                 canonical name; a table that lacks one either supplies it under a known alias
                 or gets an empty column, and the substitution is logged rather than silent.

`sequence context`
                 neither pipeline reads the ~1.5 GB ESPRESSO sequence files any more. Both read
                 a small derived table (`build_sequence_context.py` writes it) holding, per ORF,
                 exactly the sequence each pipeline actually uses: the ORF nucleotides, the exon
                 blocks, and a short flanking window either side. Scoring is unchanged - see
                 `FLANK_KEPT` for why the window is wide enough to be exact.

The initiation-codon tier tables also live here. Each step used to carry its own copy, described
in both as "identical to the other"; a single copy is what actually makes that true. The two
*scoring* implementations were left as they were - `kozak_pipeline.py` cross-checks its scores
against the initiation step's, which only catches drift while the two remain separate code.
"""

from __future__ import annotations

import ast
import logging
import os
import sys

import pandas as pd

# ======================================================================================
# Paths
# ======================================================================================

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(MODULE_DIR, "..", ".."))

DEFAULT_MASTER = os.path.join(REPO_ROOT, "Code", "data", "microprotein_master.csv")
DEFAULT_SEQUENCE_CONTEXT = os.path.join(REPO_ROOT, "Code", "data",
                                        "codon_context_sequences.csv.gz")
DEFAULT_RESULTS_ROOT = os.path.join(REPO_ROOT, "Results", "Codon_context")

# Environment variables, so a whole session can be pointed at the unreleased inputs once
# instead of on every command line.
ENV_MASTER = "MICROPROTEIN_MASTER"
ENV_SEQUENCES = "CODON_CONTEXT_SEQUENCES"
ENV_ESPRESSO = "ESPRESSO_PREFIX"


def resolve(cli_value, env_var, default):
    """CLI argument, then environment variable, then the repository default."""
    return cli_value or os.environ.get(env_var) or default


def add_common_args(parser):
    """The arguments both analysis steps share. Defaults are shown in `--help`."""
    parser.add_argument("--master", default=None,
                        help="master annotation CSV (default: %s; env $%s)"
                             % (os.path.relpath(DEFAULT_MASTER, REPO_ROOT), ENV_MASTER))
    parser.add_argument("--sequences", default=None,
                        help="sequence-context table from build_sequence_context.py "
                             "(default: %s; env $%s)"
                             % (os.path.relpath(DEFAULT_SEQUENCE_CONTEXT, REPO_ROOT),
                                ENV_SEQUENCES))
    parser.add_argument("--outdir", default=None,
                        help="output directory (default: %s/<pipeline>)"
                             % os.path.relpath(DEFAULT_RESULTS_ROOT, REPO_ROOT))
    return parser


# ======================================================================================
# Logging
# ======================================================================================

def setup_logging(name, log_dir, make_dirs=()):
    """Console + file logging, and the output directory tree. Returns the log-file path."""
    from datetime import datetime

    for d in list(make_dirs) + [log_dir]:
        os.makedirs(d, exist_ok=True)
    path = os.path.join(log_dir, "run_%s.log" % datetime.now().strftime("%Y-%m-%d_%H%M%S"))
    log = logging.getLogger(name)
    log.setLevel(logging.INFO)
    log.handlers.clear()
    fmt = logging.Formatter("%(asctime)s  %(levelname)-7s  %(message)s", "%H:%M:%S")
    for h in (logging.FileHandler(path, mode="w"), logging.StreamHandler(sys.stdout)):
        h.setFormatter(fmt)
        log.addHandler(h)
    return path


# ======================================================================================
# Initiation-codon tiers, transcribed (U -> T) from non_aug_initiation_codons.md
# ======================================================================================
#
# Compiled from Diaz de Arce 2018 (FACS-seq), Hecht 2017 (E. coli), Na 2018 (human
# N-terminomics) and Andreev 2022 (review). Tier order is the confidence hierarchy and LOWER IS
# STRONGER. Tier 5 (E. coli evidence only) is deliberately absent: with no mammalian
# confirmation it would make nearly every second codon a candidate.

TIERS = {
    "1":  ("CTG", "ACG", "GTG", "ATT", "ATA", "ATC", "TTG"),
    "2":  ("AAG", "AGG"),
    "3a": ("GCG",),
    "3b": ("TTT", "GTA", "AGC", "CGA"),
    "3c": ("GCT", "GGT", "GGG", "GGA", "GGC", "TTC", "GTC", "TCT", "TCA", "AGT",
           "CGT", "CGC", "CGG"),
}
TIER_ORDER = ["1", "2", "3a", "3b", "3c"]
TIER_RANK = {t: i for i, t in enumerate(TIER_ORDER)}
CODON_TIER = {c: t for t, cs in TIERS.items() for c in cs}
TIER_NAME = {
    "1":  "Well-Established Near-Cognate",
    "2":  "Weak Near-Cognate",
    "3a": "Non-Near-Cognate, Strongest Proteomics Evidence",
    "3b": "Non-Near-Cognate, Single-Instance Proteomics Evidence",
    "3c": "Non-Near-Cognate, Inferred (not directly reported as a TIS)",
}
COGNATE_STATUS = {"1": "near-cognate", "2": "near-cognate",
                  "3a": "non-near-cognate", "3b": "non-near-cognate", "3c": "non-near-cognate"}

STOP_CODONS = ("TAA", "TAG", "TGA")

# --- initiation context ---------------------------------------------------------------
# Positions and weights from the reported effect sizes in Diaz de Arce 2018 / Andreev 2022:
# +4 G is the strongest single determinant (up to 10x for CUG), -3 purine next (up to 6x),
# then -4 C/A (~70%), +5 A and +6 U (moderate). Shared so the two steps cannot drift on the
# *definition* of the score; each still computes it with its own code.
CONTEXT_POSITIONS = (("plus4", 3, "G", 2), ("minus3", -3, "AG", 2), ("minus4", -4, "CA", 1),
                     ("plus5", 4, "A", 1), ("plus6", 5, "T", 1))
STRONG_CONTEXT, ADEQUATE_CONTEXT = 0.65, 0.35


# ======================================================================================
# Cohort
# ======================================================================================

import re  # noqa: E402  (kept beside the pattern it exists for)

# A smORF is identified by its ORF-id grammar. Ids that do not match are UniProt / reference
# proteins, which have no ORF model here and so cannot be given a transcript context.
ORF_ID_RE = re.compile(r"^.+?[+-]chr[\w.]+:\d+-\d+_F:\d+_P:\d+(_M)?$")

# Canonical column name -> other names that carry the same quantity in a different table.
# Only unambiguous equivalences belong here; anything doubtful should be left missing so the
# output column is visibly empty rather than quietly wrong.
COLUMN_ALIASES = {
    "total_spectral_count": ("total_unique_spectral_counts",),
}


def read_header(path):
    return list(pd.read_csv(path, nrows=0).columns)


def resolve_columns(path, wanted, log):
    """
    (columns to read, {canonical: source}, [canonical names that are unavailable]).

    Requesting a column that a table does not have is normal here - the shipped master carries a
    different, smaller column set than the upstream annotation CSV - so it is reported, not
    fatal. Columns that stay missing are created empty downstream.
    """
    available = set(read_header(path))
    usecols, source, missing = [], {}, []
    for name in wanted:
        if name in available:
            usecols.append(name)
            source[name] = name
            continue
        alias = next((a for a in COLUMN_ALIASES.get(name, ()) if a in available), None)
        if alias:
            usecols.append(alias)
            source[name] = alias
            log.info("      column '%s' not present; using '%s'", name, alias)
        else:
            missing.append(name)
    if missing:
        log.warning("      columns absent from the input and left empty: %s", ", ".join(missing))
    return usecols, source, missing


def load_cohort(path, wanted, log, keep=None, dedup="gene_id"):
    """
    Read the master table in chunks and return (kept rows, other rows, tallies).

    `keep` is a predicate on a chunk returning a boolean Series selecting the rows of interest
    (e.g. non-ATG only); rows that match the smORF id grammar but fail it are still returned
    separately, because both steps account for every input row in their QC trail.

    Column naming is normalised here: whatever the source called a column, the frame this
    returns uses the canonical names in `wanted`.
    """
    usecols, source, missing = resolve_columns(path, wanted, log)
    rename = {src: canon for canon, src in source.items() if src != canon}

    kept, other, n_rows, n_ids = [], [], 0, set()
    for chunk in pd.read_csv(path, usecols=usecols, chunksize=200_000, low_memory=False):
        chunk = chunk.rename(columns=rename)
        n_rows += len(chunk)
        ids = chunk.gene_id.astype(str)
        n_ids.update(ids)
        selected = chunk if keep is None else chunk[keep(chunk)]
        is_smorf = selected.gene_id.astype(str).map(lambda s: ORF_ID_RE.match(s) is not None)
        kept.append(selected[is_smorf])
        other.append(selected[~is_smorf])

    cohort = pd.concat(kept, ignore_index=True)
    other = pd.concat(other, ignore_index=True)
    for name in missing:
        cohort[name] = pd.NA
        other[name] = pd.NA
    cohort = cohort.drop_duplicates(dedup).reset_index(drop=True)
    cohort["gene_id"] = cohort.gene_id.astype(str)
    return cohort, other, {"input_rows": n_rows, "input_ids": len(n_ids)}


# ======================================================================================
# Sequence context
# ======================================================================================
#
# `build_sequence_context.py` writes this table; both analysis steps read it and neither touches the
# ESPRESSO source files. One row per ORF that could be located in them.
#
# FLANK_KEPT is what makes the substitution exact rather than approximate. The widest reach of
# any scoring position is 6 nt beyond the ORF at either end (-6 for the Kozak window, +6 for a
# candidate sited at the ORF's last codon), so 12 nt of retained flank is twice what any
# calculation can consult. Where the transcript itself provides less, the stored flank is
# correspondingly short and the missing positions score as '-' - exactly as they did when the
# full transcript was in memory.

FLANK_KEPT = 12

SEQUENCE_CONTEXT_COLUMNS = ["orf_id", "chromosome", "strand", "exon_blocks",
                            "orf_nt", "upstream_nt", "downstream_nt", "tx_offset", "tx_copies"]

# Returned for an ORF absent from the sequence table. Both pipelines' `validate()` turn `nt is
# None` into a QC exclusion, so a missing ORF is dropped with a reason rather than crashing.
MISSING_SEQUENCE = {"nt": None, "exons": None, "strand": None, "chrom": None,
                    "upstream": "", "downstream": "", "tx_offset": -1, "tx_copies": 0}


def _parse_blocks(text):
    if not isinstance(text, str) or not text:
        return None
    return [tuple(int(v) for v in block.split("-")) for block in text.split(";")]


def load_sequence_context(path, wanted, log):
    """
    {orf_id: sequence record} for every ORF in `wanted`, from the derived table.

    ORFs the table does not cover get `MISSING_SEQUENCE`, so the caller's QC stage reports them
    instead of failing.
    """
    if not os.path.exists(path):
        raise SystemExit(
            "sequence-context table not found: %s\n"
            "Build it once from the ESPRESSO sequence files:\n"
            "    python build_sequence_context.py --espresso-prefix <prefix>\n"
            "See README.md; the source files are not distributed with this repository." % path)

    table = pd.read_csv(path, low_memory=False)
    missing_cols = [c for c in SEQUENCE_CONTEXT_COLUMNS if c not in table.columns]
    if missing_cols:
        raise SystemExit("%s is missing columns %s; rebuild it with build_sequence_context.py"
                         % (path, missing_cols))

    seqs = {}
    for r in table.itertuples():
        seqs[str(r.orf_id)] = {
            "nt": str(r.orf_nt).upper(),
            "exons": _parse_blocks(r.exon_blocks),
            "strand": r.strand, "chrom": r.chromosome,
            "upstream": "" if pd.isna(r.upstream_nt) else str(r.upstream_nt).upper(),
            "downstream": "" if pd.isna(r.downstream_nt) else str(r.downstream_nt).upper(),
            "tx_offset": int(r.tx_offset), "tx_copies": int(r.tx_copies)}

    found = {orf: seqs[orf] for orf in wanted if orf in seqs}
    if len(found) < len(wanted):
        log.warning("      %d of %d ORFs are absent from %s and will fail QC",
                    len(wanted) - len(found), len(wanted), os.path.basename(path))
    return {orf: found.get(orf, dict(MISSING_SEQUENCE)) for orf in wanted}


def context_sequence(seq):
    """
    (sequence to score against, offset of the ORF's first base within it).

    The pipelines used to index into the whole assembled transcript. This reconstructs the only
    part of it they could ever reach. When the ORF was never located in a transcript
    (`tx_offset < 0`) there are no flanks and scoring falls back to the ORF alone, which is what
    the transcript-based code did in that case too.
    """
    if seq["tx_offset"] < 0:
        return seq["nt"], 0
    return seq["upstream"] + seq["nt"] + seq["downstream"], len(seq["upstream"])


# ======================================================================================
# Formatting helpers used by both reports
# ======================================================================================

def as_list(values):
    """Bracketed Python-style list, matching the input CSV's own convention."""
    return repr(list(values))


def parse_list(value):
    """
    Inverse of `as_list`; tolerant of nulls and bare strings.

    Missing values arrive as `None`, `nan` or `pd.NA` depending on whether the column was absent
    from the input table or merely empty in a row, so all three are checked - `pd.NA` in
    particular stringifies to `'<NA>'` and would otherwise be parsed as a one-element list.
    """
    if value is None or (not isinstance(value, (list, tuple, str)) and pd.isna(value)):
        return []
    if isinstance(value, list):
        return value
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return []
    try:
        parsed = ast.literal_eval(text)
    except (ValueError, SyntaxError):
        return [text]
    return list(parsed) if isinstance(parsed, (list, tuple)) else [parsed]


def fmt_p(p):
    """Tiny p-values underflow to exactly 0.0, which would print as a false certainty."""
    return "< 1e-300" if p == 0 else "%.3g" % p


def md_table(df, cols=None):
    """Markdown table. Counts render as integers, everything else to three decimals."""
    df = df[cols] if cols else df
    is_count = {c: c == "n" or c.endswith(" n") or c.startswith("n_") for c in df.columns}

    def cell(name, v):
        if isinstance(v, float) and pd.isna(v):
            return "–"
        if isinstance(v, (int, float)) and not isinstance(v, bool):
            if is_count[name]:
                return "{:,}".format(int(v))
            return ("%.3f" % v).rstrip("0").rstrip(".") if isinstance(v, float) else str(v)
        return str(v)

    head = "| " + " | ".join(df.columns) + " |"
    rule = "|" + "|".join("---:" if pd.api.types.is_numeric_dtype(df[c]) else "---"
                          for c in df.columns) + "|"
    body = ["| " + " | ".join(cell(c, v) for c, v in zip(df.columns, row)) + " |"
            for row in df.itertuples(index=False)]
    return "\n".join([head, rule] + body)
