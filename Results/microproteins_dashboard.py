import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from pathlib import Path
from urllib.parse import quote, unquote
import numpy as np
import os
import ast
import re
import html
import hashlib
import struct

# scRNA Enrichment view is always enabled.
INCLUDE_SCRNA = True


# =============================================================================
# REMOTE ASSET HOSTING (Hugging Face Datasets fallback)
# =============================================================================
# When the large figure directories (mirror_plots/, expression_profiles/,
# smorf_cartoon_figures/) are not present locally (e.g. on Streamlit Cloud),
# the dashboard streams the images from the public Hugging Face dataset.
# Override via .streamlit/secrets.toml: assets_base_url = "https://..."
DEFAULT_ASSETS_BASE_URL = "https://huggingface.co/datasets/brmiller/brain-microprotein-atlas/resolve/main"
HF_REPO_ID = "brmiller/brain-microprotein-atlas"

def _get_assets_base_url():
    try:
        return st.secrets.get("assets_base_url", DEFAULT_ASSETS_BASE_URL).rstrip("/")
    except Exception:
        return DEFAULT_ASSETS_BASE_URL


@st.cache_data(show_spinner=False)
def _list_remote_files(_repo=HF_REPO_ID):
    """Return the full file listing of the remote HF dataset (cached)."""
    try:
        from huggingface_hub import HfApi
        return HfApi().list_repo_files(repo_id=_repo, repo_type="dataset")
    except Exception as e:
        st.warning(f"Could not list remote assets ({e}); figures may be missing.")
        return []


def _remote_url(rel_path):
    """Build a public URL for a path within the HF dataset."""
    return f"{_get_assets_base_url()}/{quote(rel_path)}"


@st.cache_data(show_spinner=False)
def _localize_hf_asset(rel_path, _repo=HF_REPO_ID):
    """Download a dataset asset and return a cached local path (or None).

    The HF dataset is stored on the Xet backend, so the plain
    resolve/main/<file> URLs 403 on a direct browser GET (what st.image would
    issue). hf_hub_download reconstructs the file via the xet client
    (anonymous, requires the hf_xet package) and caches it on disk, so we serve
    that local copy to st.image instead of the un-fetchable URL.
    """
    try:
        from huggingface_hub import hf_hub_download
        return hf_hub_download(repo_id=_repo, repo_type="dataset", filename=rel_path)
    except Exception:
        return None


def _resolve_image_src(src):
    """Resolve an image/PDF source for st.image / iframe rendering.

    Local paths pass through unchanged. HF dataset URLs are Xet-backed and 403
    on a direct browser fetch, so download them server-side (cached) and return
    the local path; fall back to the original URL if the download fails.
    """
    if not src or not isinstance(src, str):
        return src
    if not src.startswith(("http://", "https://")):
        return src
    marker = "/resolve/main/"
    if "huggingface.co" in src and marker in src:
        rel = unquote(src.split(marker, 1)[1])
        local = _localize_hf_asset(rel)
        if local:
            return local
    return src


# Reference render width for the entry page's full-width figures. The body
# column is ~1130px on a 1600px window; images carry max-width:100% (see the
# stylesheet) so they shrink on narrower windows rather than overflowing.
ENTRY_FIGURE_WIDTH = 1130


@st.cache_data(show_spinner=False)
def _png_size(path):
    """(width, height) of a PNG, read from its 24-byte header. None if not a PNG.

    Cheap enough to call per render — no image library, no full decode.
    """
    try:
        with open(path, "rb") as fh:
            head = fh.read(24)
    except OSError:
        return None
    if len(head) < 24 or head[:8] != b"\x89PNG\r\n\x1a\n":
        return None
    w, h = struct.unpack(">II", head[16:24])
    return (int(w), int(h)) if w and h else None


_MEDIABOX_RE = re.compile(
    rb"/MediaBox\s*\[\s*([\d.+-]+)\s+([\d.+-]+)\s+([\d.+-]+)\s+([\d.+-]+)\s*\]"
)


def _mediabox(blob):
    """(width, height) from the first /MediaBox in a byte blob, or None."""
    m = _MEDIABOX_RE.search(blob)
    if not m:
        return None
    try:
        x0, y0, x1, y1 = (float(v) for v in m.groups())
    except ValueError:
        return None
    w, h = abs(x1 - x0), abs(y1 - y0)
    return (w, h) if w and h else None


@st.cache_data(show_spinner=False)
def _pdf_size(path):
    """(width, height) in points of a PDF's first page, from its /MediaBox.

    The PNG path gets its aspect ratio from a 24-byte header; PDFs need the
    same thing so their figures can reserve height too. These are matplotlib
    exports that use compressed object streams, so /MediaBox is usually *not*
    findable in the raw bytes — inflating each Flate stream and searching the
    result is what actually finds it. No PDF library, no full parse; the files
    are ~60 KB, so brute force is cheap and the result is cached per path.
    """
    try:
        with open(path, "rb") as fh:
            data = fh.read()
    except OSError:
        return None
    if not data.startswith(b"%PDF"):
        return None

    # Uncompressed case first — cheapest, and covers PDFs from other tools.
    size = _mediabox(data)
    if size:
        return size

    import zlib
    for m in re.finditer(rb"stream\r?\n", data):
        start = m.end()
        end = data.find(b"endstream", start)
        if end < 0:
            continue
        try:
            inflated = zlib.decompress(data[start:end])
        except zlib.error:
            continue  # not Flate, or not a self-contained stream — skip it
        size = _mediabox(inflated)
        if size:
            return size
    return None


def _pdf_reserved(src, width=ENTRY_FIGURE_WIDTH, fallback_height=420):
    """Inline PDF in a height-reserved container, mirroring _image_reserved.

    An <iframe> does reserve its own box, but only at whatever height the call
    site guessed — so a figure toggled between a PNG and a PDF rendition jumped
    as the box resized. Deriving the height from the page's own aspect ratio
    keeps both renditions the same shape and keeps the anchor targets still.
    """
    resolved = _resolve_image_src(src)
    size = _pdf_size(resolved) if isinstance(resolved, str) else None
    height = int(round(width * size[1] / size[0])) if size else fallback_height
    with st.container(height=height + 8, border=False):
        _display_pdf_inline(resolved, height=height)


def _image_reserved(src, width=ENTRY_FIGURE_WIDTH):
    """st.image, with the image's height reserved *before* it loads.

    An <img> with no intrinsic size known to the browser occupies zero height
    until it decodes, so everything below it jumps down when it arrives. On the
    Hugging Face Space the figures stream from the dataset rather than sitting on
    local disk, so that arrival is seconds late — long enough that clicking a
    section in the entry nav scrolls to a target that then slides out from under
    you (measured: 584px of drift; the browser will not compensate, scroll
    anchoring is off app-wide and re-enabling it made no difference).

    Reading the PNG header gives the aspect ratio server-side, so the height can
    be reserved with a fixed-size container and the page stops moving.
    """
    resolved = _resolve_image_src(src)
    size = _png_size(resolved) if isinstance(resolved, str) else None
    if not size:
        # Not a local PNG (e.g. an un-fetchable remote URL) — nothing to reserve.
        st.image(resolved, width=width)
        return
    # +8px so rounding can never leave the image taller than its own box, which
    # would give the container an inner scrollbar.
    height = int(round(width * size[1] / size[0])) + 8
    with st.container(height=height, border=False):
        st.image(resolved, width=width)


# =============================================================================
# PROSIT MIRROR PLOT CONFIGURATION
# =============================================================================
MIRROR_PLOT_BASE = Path(__file__).parent / "mirror_plots"
QUALITY_LEVELS = ['Strong', 'Moderate', 'Weak', 'Insufficient', 'No PROSIT']
QUALITY_EMOJI = {
    'Strong': 'Strong',
    'Moderate': 'Moderate',
    'Weak': 'Weak',
    'Insufficient': 'Insufficient',
    'No PROSIT': 'No PROSIT',
}

# ── TMT significance tiers (best/most-stringent tier a microprotein reaches) ──
# Tier 1 = strongest signal (Green) … Tier 4 = weakest (Red).
#   Tier 1: TMT q < 0.05 in ≥50% samples/condition
#   Tier 2: TMT q < 0.2  in ≥50% samples/condition
#   Tier 3: TMT q < 0.05 in ≥1 sample/condition
#   Tier 4: TMT q < 0.2  in ≥1 sample/condition
TMT_TIER_EMOJI = {
    'Tier 1': 'Tier 1',
    'Tier 2': 'Tier 2',
    'Tier 3': 'Tier 3',
    'Tier 4': 'Tier 4',
}
# ── ROSMAP RNA-seq significance (Yellow = Exploratory, Green = Significant) ──
RNA_SIG_EMOJI = {
    'Significant': 'Significant',
    'Exploratory': 'Exploratory',
}

# ── smORF decay class (NMD/NSD status × main-ORF disruption) ──────────────────
# From the external NMD analysis shipped as
# Code/data/smorf_trx_priority_color_class_assignments.tsv, whose `category`
# column encodes the NMD × main-ORF-disruption 2×2 plus two non-stop-decay
# (NSD) classes, as genome-browser track colors. NSD is decided per smORF over
# the whole transcript set — one confident no-stop transcript and none with a
# real stop — and overrides the 2×2, so it is a third decay state, not a fifth
# color. Dots below match those track colors so the dashboard reads the same
# as the BED in a browser.
# Covers Salk/TrEMBL smORFs only — Swiss-Prot microproteins were never assessed
# upstream, so they legitimately carry 'NA'.
DECAY_CLASS_TSV = (Path(__file__).resolve().parent.parent / "Code" / "data" /
                   "smorf_trx_priority_color_class_assignments.tsv")
DECAY_CLASS_LEVELS = ['No NMD · mORF Shift', 'No NMD · mORF Intact',
                      'NMD · mORF Shift', 'NMD · mORF Intact',
                      'NSD · mORF Shift', 'NSD · mORF Intact', 'NA']
DECAY_CLASS_FROM_CATEGORY = {
    'dark_green':    'No NMD · mORF Shift',
    'light_green':   'No NMD · mORF Intact',
    'gold':          'NMD · mORF Shift',
    'firebrick':     'NMD · mORF Intact',
    'nsd_disrupted': 'NSD · mORF Shift',
    'nsd_intact':    'NSD · mORF Intact',
}
DECAY_CLASS_EMOJI = {
    'No NMD · mORF Shift':  'No NMD · mORF Shift',
    'No NMD · mORF Intact': 'No NMD · mORF Intact',
    'NMD · mORF Shift':     'NMD · mORF Shift',
    'NMD · mORF Intact':    'NMD · mORF Intact',
    'NSD · mORF Shift':     'NSD · mORF Shift',
    'NSD · mORF Intact':    'NSD · mORF Intact',
    'NA':                   'NA',
}
DECAY_CLASS_HELP = (
    'Predicted decay fate of the host transcript × whether the smORF shifts '
    'the main ORF. NMD = stop trips the 50-nt rule; NSD = non-stop decay, an '
    'ORF with no in-frame stop. NSD is judged over the whole transcript set — '
    'it needs a confident no-stop transcript and no transcript carrying a real '
    'stop — so it is scored ahead of, and overrides, the NMD × disruption 2×2. '
    'NA = not assessed (Swiss-Prot microproteins were not covered by this '
    'analysis).'
)


# =============================================================================
# MAIN-TABLE COLUMN DESCRIPTIONS
# =============================================================================
# One line per display column name (as used in display_df / _column_groups).
# Single source of truth for the "Explain the column/variables to me" expander —
# shows only the descriptions for whatever columns are visible in the active
# column-view tab.
COLUMN_DESCRIPTIONS = {
    'DB': 'Database source: 🔵 reviewed Swiss-Prot vs 🟠 unreviewed Salk/TrEMBL discovery.',
    'sequence': 'The annotated microprotein amino-acid sequence, as originally called (from its annotated start through the stop codon).',
    'NonATG Alternative Sequence(s)': "Predicted sequence(s) from alternative (non-AUG) initiation candidates (nonATG_predicted_sequence), raw list format; blank for ATG-start microproteins or non-ATG starts with no compatible alternative site.",
    'NonATG Candidates': 'Alternative initiation candidates found by scanning, as codon@position (e.g. "CTG@23"); blank if the start is ATG or no alternative site was found.',
    'UCSC': "Link to view this microprotein's genomic locus in the UCSC Genome Browser.",
    'Parent Gene': 'Host gene symbol/name the smORF resides within or near.',
    'General smORF Type': 'Broad smORF category (e.g. Upstream, Downstream, Short-Isoform, lncRNA-encoded, Swiss-Prot).',
    'smORF Subtype': 'Specific smORF subtype/classification within its general type.',
    'ORF Rules': 'Predicted decay fate of the host transcript × whether the smORF shifts the main ORF (no NMD/mORF shift, no NMD/mORF intact, NMD/mORF shift, NMD/mORF intact, NSD/mORF shift, NSD/mORF intact). NMD = stop trips the 50-nt rule; NSD = non-stop decay, an ORF with no in-frame stop — scored per smORF over the whole transcript set and overriding the NMD call. NA = not assessed; Swiss-Prot microproteins were not covered by this analysis.',
    'Protein Length': 'Length of the microprotein in amino acids.',
    'Start Codon': 'Coarse call: ATG (canonical) vs nonATG (annotated with a near-cognate/non-standard start).',
    'Kozak Strength': 'Full Kozak consensus context strength around the start codon (weak/adequate/strong); Salk smORFs only.',
    'Kozak Downstream Strength': 'Kozak context strength considering only the +4 downstream position.',
    'Kozak Class': 'Canonical combined Kozak strength classification.',
    'Kozak Window': 'Local sequence context around the start codon (lowercase = flanking, uppercase = start + downstream).',
    'PhyloCSF Score': 'PhyloCSF evolutionary conservation score.',
    'ShortStop': 'ShortStop ribosome-profiling classification label (e.g. SAM-Secreted, SAM-Intracellular).',
    'ShortStop Score': 'ShortStop ML confidence score (0–1).',
    'Annotation Method': 'How the microprotein was annotated: MS (mass-spec detected) or RiboCode-ShortStop (translation evidence without MS detection).',
    'Spectra Quality': 'Best PROSIT spectral-match confidence tier (from the master Confidence column).',
    'Nt-Acetylated': 'At least one tryptic peptide carries an N-terminal acetyl mark — direct evidence that this is a genuine protein N-terminus.',
    'Nt-Acetyl Peptides': 'Tryptic peptide sequence(s) observed with an N-terminal acetyl mark.',
    'Nt-Acetyl PSMs': 'Total Nt-acetylated PSMs summed across the Nt-acetylated peptides.',
    'Nt-Acetyl PSM Fraction': 'Highest per-peptide fraction of that peptide’s PSMs carrying the Nt-acetyl mark (1.0 = every PSM acetylated).',
    'TMT Tier': 'TMT significance tier: Tier 1 (strongest) … Tier 4.',
    'RNA Significance': 'ROSMAP RNA-seq significance: Significant, Exploratory.',
    'UniProt Annotation': 'UniProt curation annotation score (1–5★) for reviewed/TrEMBL entries; blank for novel microproteins not in UniProt.',
    'NonATG Codon': 'The originally annotated (non-ATG) start codon.',
    'NonATG Valid Initiator': 'Whether the annotated codon is itself a plausible translation initiator.',
    'NonATG Has Supported Site': 'Whether current non-canonical initiation evidence/models support this annotated start actually being used (holistic verdict, distinct from Valid Initiator); Salk smORFs only.',
    'NonATG Context Strength': 'Kozak context strength around the annotated non-ATG codon.',
    'NonATG Alt Sites Found': 'Number of peptide-compatible alternative initiation sites found by scanning downstream of the annotated start.',
    'Optimal Codon Tier': 'Strongest evidence tier among the alternative initiation candidates (Well-Established Near-Cognate is best).',
    'Unique Spectral Counts (DDA)': 'Number of unique mass-spectrometry spectral counts (DDA).',
    'Razor Counts (DDA)': 'Total razor spectral counts (shared peptides assigned to this protein).',
    'TMT log2FC (50%)': 'TMT log2 fold-change (AD vs Control) — 50% missing-value threshold.',
    'TMT t-stat (50%)': 'TMT t-statistic (AD vs Control) — 50% missing-value threshold.',
    'TMT df (50%)': 'TMT t-test degrees of freedom — 50% missing-value threshold.',
    'TMT CI low (50%)': 'TMT 95% CI lower bound on log2 fold-change — 50% missing-value threshold.',
    'TMT CI high (50%)': 'TMT 95% CI upper bound on log2 fold-change — 50% missing-value threshold.',
    "TMT Cohen's d (50%)": "TMT Cohen's d effect size (AD vs Control) — 50% missing-value threshold.",
    'TMT p-val (50%)': 'TMT raw p-value — 50% missing-value threshold.',
    'TMT q-val (50%)': 'TMT q-value (BH-adjusted) — 50% missing-value threshold.',
    'TMT log2FC (0%)': 'TMT log2 fold-change (AD vs Control) — 0% missing-value threshold (stringent).',
    'TMT t-stat (0%)': 'TMT t-statistic (AD vs Control) — 0% missing-value threshold (stringent).',
    'TMT df (0%)': 'TMT t-test degrees of freedom — 0% missing-value threshold (stringent).',
    'TMT CI low (0%)': 'TMT 95% CI lower bound on log2 fold-change — 0% missing-value threshold (stringent).',
    'TMT CI high (0%)': 'TMT 95% CI upper bound on log2 fold-change — 0% missing-value threshold (stringent).',
    "TMT Cohen's d (0%)": "TMT Cohen's d effect size (AD vs Control) — 0% missing-value threshold (stringent).",
    'TMT p-val (0%)': 'TMT raw p-value — 0% missing-value threshold (stringent).',
    'TMT q-val (0%)': 'TMT q-value (BH-adjusted) — 0% missing-value threshold (stringent).',
    'MS Detect Control': 'MS detection rate in Control donors.',
    'MS Detect AD': 'MS detection rate in AD donors.',
    'Tryptic Peptides': 'Observed tryptic peptide sequence(s) identified by mass spec.',
    'Tryptic Protein ID': 'Protein ID associated with the tryptic peptides.',
    'Tryptic Start Positions': 'Start position(s) of the tryptic peptide(s) within the annotated sequence.',
    'Tryptic End Positions': 'End position(s) of the tryptic peptide(s) within the annotated sequence.',
    'ROSMAP log2FC': 'ROSMAP short-read RNA-seq log2 fold-change (AD vs Control).',
    'ROSMAP padj': 'ROSMAP short-read RNA-seq adjusted p-value.',
    'MSBB log2FC': 'MSBB short-read RNA-seq log2 fold-change (AD vs Control).',
    'MSBB padj': 'MSBB short-read RNA-seq adjusted p-value.',
    'ROSMAP Corr NonAD': 'ROSMAP: correlation between Main ORF and smORF transcripts (non-AD donors).',
    'ROSMAP Corr AD': 'ROSMAP: correlation between Main ORF and smORF transcripts (AD donors).',
    'MSBB Corr NonAD': 'MSBB: correlation between Main ORF and smORF transcripts (non-AD donors).',
    'MSBB Corr AD': 'MSBB: correlation between Main ORF and smORF transcripts (AD donors).',
    'RNA_LRT_Add_P': 'ROSMAP RNA-seq LRT additive-model p-value.',
    'RNA_LRT_Int_P': 'ROSMAP RNA-seq LRT interaction-model p-value.',
    'RP3 Default': 'RP3 default pipeline ribosome-profiling result.',
    'RP3 MM+Amb': 'RP3 multi-mapping + ambiguous-reads result.',
    'RP3 Amb': 'RP3 ambiguous-reads result.',
    'RP3 MM': 'RP3 multi-mapping result.',
    'RiboCode': 'RiboCode ORF-detection result.',
    'BLAST UniProt Match': 'UniProt accession of the best BLASTp hit.',
    'BLAST % Match': 'BLASTp percent identity to the best UniProt hit.',
    'BLAST Aln Length': 'BLASTp alignment length (residues).',
    'BLAST E-value': 'BLASTp expectation value of the best UniProt hit.',
    'BLAST Bit Score': 'BLASTp bit score of the best UniProt hit.',
}


# ── UniProt annotation score (1–5 curation completeness, rendered as gold stars) ──
def _annotation_stars(score):
    """Render a 1–5 UniProt annotation score as gold stars ('★★★☆☆'); '—' if missing."""
    if score is None or pd.isna(score):
        return '\u2014'
    try:
        n = int(round(float(score)))
    except (ValueError, TypeError):
        return '\u2014'
    n = max(0, min(5, n))
    return '\u2605' * n + '\u2606' * (5 - n)

# =============================================================================
# smORF CARTOON FIGURE CONFIGURATION
# =============================================================================
SMORF_CARTOON_DIR = Path(__file__).parent / "smorf_cartoon_figures"
SEQ_TO_COORDS_FILE = Path(__file__).parent.parent / "GTF_and_BED_files" / "Unreviewed_Brain_Microproteins_mapping_coordinates_to_sequences.tsv"

# =============================================================================
# EXPRESSION PROFILE CONFIGURATION
# =============================================================================
EXPRESSION_PROFILE_DIR = Path(__file__).parent / "expression_profiles"
GENOME_FILES_DIR = Path(__file__).parent.parent / "GTF_and_BED_files"

# =============================================================================
# smORF TYPE GROUPINGS
# =============================================================================
# smORF-type -> parent group. Kept in sync with the figure generators'
# map_smorf_label (revision_figures/_stats_figure{1..4}): udORF and isoORF are
# grouped as TrEMBL/altORF and Short-Isoform respectively.
SMORF_PARENT_GROUPS = {
    'Upstream': ['uORF', 'uoORF', 'uaORF', 'uaoORF'],
    'Downstream': ['dORF', 'doORF', 'daORF', 'daoORF'],
    'Internal ORF': ['iORF', 'oORF'],
    'TrEMBL/AltORF': ['eORF', 'udORF', 'TrEMBL'],
    'Short-Isoform': ['Iso', 'D-Iso', 'N-Iso', 'isoORF'],
    'lncRNA': ['lncRNA'],
    'psORF': ['psORF'],
    'Swiss-Prot': ['Swiss-Prot-MP'],
}
SMORF_CHILD_TO_PARENT = {
    child: parent
    for parent, children in SMORF_PARENT_GROUPS.items()
    for child in children
}
SMORF_DISPLAY_LABEL = {}

# =============================================================================
# scRNA CELL-TYPE ENRICHMENT CONFIGURATION
# =============================================================================
# Per-microprotein cell-type panel in the ID card, the single-gene analogue of
# the cell-type heatmaps in Code/scRNAseq_summary_merging_analysis/
# scRNAseq_summary.R: same diverging log2FC fill and the same
# cell_type_general grouping, but one microprotein per plot instead of an
# aggregate over smORF types.
SCRNA_SUMMARY_FILE = Path(__file__).parent / "scRNA_Enrichment" / "scRNA_Enrichment_summary.csv"
# Every tested (microprotein, cell type) pair, including the non-significant
# ones the summary above drops. Written by the same R script; the panel needs it
# to distinguish "tested, not significant" from "not tested".
SCRNA_ALL_CELLTYPES_FILE = Path(__file__).parent / "scRNA_Enrichment" / "scRNA_Enrichment_all_celltypes.csv.gz"

# Facet order from the R script's `ordered_cell_types`, with the catch-all last.
SCRNA_GENERAL_ORDER = ['Astr', 'Exc Neur', 'Immune', 'Inh Neur', 'Oli', 'Vas', 'Other']

# Okabe-Ito (colorblind-safe) — one hue per general cell type. Carries the
# grouping that facet_grid(rows = vars(cell_type_general)) carries in the R
# figure, which a single-column plotly heatmap has no strips for.
SCRNA_GENERAL_COLORS = {
    'Astr':     '#0072B2',
    'Exc Neur': '#D55E00',
    'Immune':   '#009E73',
    'Inh Neur': '#CC79A7',
    'Oli':      '#E69F00',
    'Vas':      '#56B4E9',
    'Other':    '#999999',
}

# scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027") from
# the R heatmaps, as a plotly colorscale.
SCRNA_FC_COLORSCALE = [[0.0, '#4575b4'], [0.5, '#ffffff'], [1.0, '#d73027']]

# The R heatmaps squish mean log2FC to ±0.3. Single-gene effects run larger, so
# the panel expands past this floor when the microprotein needs it (see
# _scrna_fc_limit) rather than saturating every tile.
SCRNA_FC_FLOOR = 0.3

# =============================================================================
# UCSC CUSTOM SESSION CONFIGURATION
# =============================================================================
CUSTOM_UCSC_SESSION = "AD Dark Microproteome"
CUSTOM_UCSC_USERNAME = "brmiller"
CUSTOM_UCSC_SESSION_URL = "https://genome.ucsc.edu/s/brmiller/AD%20Dark%20Microproteome"

COLORS = {
    'swiss_prot': '#74a2b7',
    'unreviewed': '#ed8651',
}


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# ── Entry-page addressing ────────────────────────────────────────────────────
# `sequence` is the de-duplication key of the merged master table (see the
# drop_duplicates in load_and_merge_all_data), so it is the microprotein's
# identity. Sequences are far too long for a URL, so the entry page is
# addressed by a short, stable digest of the sequence: ?mp=<entry id>.
ENTRY_QUERY_PARAM = "mp"


def _entry_id(seq):
    """Stable short accession-like id for a microprotein sequence."""
    s = str(seq or '')
    if not s:
        return ''
    return hashlib.sha1(s.encode('utf-8')).hexdigest()[:10]


@st.cache_data
def build_entry_id_index(_df, cache_key):
    """Map entry id -> sequence, over the *unfiltered* set.

    Built from the post-load frame rather than the filtered view so that an
    entry URL resolves regardless of which filters happen to be active.

    `_df` is underscore-prefixed (excluded from the cache key, per the
    convention used by the other index builders) to skip re-hashing a large
    frame on every rerun; `cache_key` — the row count — is what actually
    invalidates the entry, so pass it explicitly.
    """
    if _df is None or 'sequence' not in _df.columns:
        return {}
    seqs = _df['sequence'].dropna().astype(str).unique()
    return {_entry_id(s): s for s in seqs}


@st.cache_data
def build_seq_to_coords_index(_path=str(SEQ_TO_COORDS_FILE)):
    """Build sequence-to-coordinates lookup from the mapping file."""
    p = Path(_path)
    if not p.exists():
        return {}
    df = pd.read_csv(p, sep='\t')
    return dict(zip(df['sequence'], df['genomic_coordinates']))


@st.cache_data
def build_mirror_plot_index(_base=str(MIRROR_PLOT_BASE)):
    """Build index of available mirror plot images by peptide.

    Local-first: scans the on-disk directory if present, otherwise pulls the
    file listing from the Hugging Face dataset and builds URL-based entries.
    """
    base = Path(_base)
    index = {}

    def _add_entry(peptide, quality, parts, filepath):
        charge = parts[1].replace("z", "") if len(parts) > 1 else "?"
        scan = parts[-1] if len(parts) > 2 else "?"
        if peptide not in index:
            index[peptide] = {'best_quality': quality, 'plots': []}
        index[peptide]['plots'].append({
            'quality': quality,
            'charge': charge,
            'scan': scan,
            'filepath': filepath,
            'peptide': peptide,
        })

    if base.exists():
        for quality in QUALITY_LEVELS:
            quality_dir = base / quality
            if not quality_dir.exists():
                continue
            for img_file in quality_dir.glob("*.png"):
                parts = img_file.stem.split("_")
                if len(parts) >= 3:
                    _add_entry(parts[0], quality, parts, str(img_file))
        return index

    # Remote fallback
    valid_qualities = set(QUALITY_LEVELS)
    for rel in _list_remote_files():
        if not rel.startswith("mirror_plots/") or not rel.endswith(".png"):
            continue
        segs = rel.split("/")
        if len(segs) != 3:
            continue
        _, quality, fname = segs
        if quality not in valid_qualities:
            continue
        parts = Path(fname).stem.split("_")
        if len(parts) >= 3:
            _add_entry(parts[0], quality, parts, _remote_url(rel))
    return index


@st.cache_data
def build_expression_profile_index(_base=str(EXPRESSION_PROFILE_DIR)):
    """Build index of expression profile images keyed by genomic coordinates.

    Local-first; falls back to HF dataset file listing if local dir missing.
    """
    base = Path(_base)
    index = {}

    def _add(stem, coupling, filepath):
        match = re.search(r'(chr\w+)_(\d+)-(\d+)$', stem)
        if not match:
            return
        chrom, start, end = match.group(1), match.group(2), match.group(3)
        coords_key = f"{chrom}:{start}-{end}"
        gene_name = stem[:stem.rfind(f'_{chrom}')]
        index.setdefault(coords_key, []).append({
            'filepath': filepath,
            'coupling': coupling,
            'gene': gene_name,
        })

    if base.exists():
        for subdir in ['coupled', 'non_coupled']:
            subpath = base / subdir
            if not subpath.exists():
                continue
            seen_stems = set()
            for img_file in list(subpath.glob("*.png")) + list(subpath.glob("*.pdf")):
                if img_file.stem in seen_stems:
                    continue
                seen_stems.add(img_file.stem)
                _add(img_file.stem, subdir, str(img_file))
        return index

    # Remote fallback: prefer PNG over PDF per stem
    by_stem = {}  # stem -> (rel_path, coupling, ext_priority)
    for rel in _list_remote_files():
        if not rel.startswith("expression_profiles/"):
            continue
        segs = rel.split("/")
        if len(segs) != 3:
            continue
        _, subdir, fname = segs
        if subdir not in ('coupled', 'non_coupled'):
            continue
        stem, ext = Path(fname).stem, Path(fname).suffix.lower()
        if ext not in ('.png', '.pdf'):
            continue
        priority = 0 if ext == '.png' else 1
        existing = by_stem.get(stem)
        if existing is None or priority < existing[2]:
            by_stem[stem] = (rel, subdir, priority)
    for stem, (rel, subdir, _) in by_stem.items():
        _add(stem, subdir, _remote_url(rel))
    return index


def _display_pdf_inline(filepath, height=500):
    """Display a PDF inline. Local files are base64-embedded; URLs are streamed."""
    filepath = _resolve_image_src(filepath)
    try:
        if isinstance(filepath, str) and filepath.startswith(("http://", "https://")):
            src = filepath
        else:
            import base64
            with open(filepath, "rb") as f:
                b64 = base64.b64encode(f.read()).decode()
            src = f"data:application/pdf;base64,{b64}"
        st.markdown(
            f'<iframe src="{src}" '
            f'width="100%" height="{height}px" style="border:none; border-radius:8px;"></iframe>',
            unsafe_allow_html=True,
        )
    except Exception as e:
        st.warning(f"Could not display PDF: {e}")


def get_spectra_quality(tryptic_peptides_str, mirror_index):
    """Get best spectra quality for a microprotein's tryptic peptides."""

    if not tryptic_peptides_str or pd.isna(tryptic_peptides_str):
        return 'No PROSIT'
    if not mirror_index:
        return 'No PROSIT'
    try:
        peptides = ast.literal_eval(str(tryptic_peptides_str))
        if isinstance(peptides, str):
            peptides = [peptides]
    except (ValueError, SyntaxError):
        peptides = []
    # Default to 'No PROSIT': a microprotein with no PROSIT grade and no matching
    # mirror-plot spectrum has no spectral match, so it must NOT be counted as
    # 'Insufficient' (that tier is reserved for genuine low-quality PROSIT hits).
    best = 'No PROSIT'
    for pep in peptides:
        pep = pep.strip()
        if pep in mirror_index:
            q = mirror_index[pep]['best_quality']
            if QUALITY_LEVELS.index(q) < QUALITY_LEVELS.index(best):
                best = q
    return best


def get_matching_mirror_plots(tryptic_peptides_str, mirror_index):
    """Get all mirror plots for a microprotein's tryptic peptides."""
    if not tryptic_peptides_str or pd.isna(tryptic_peptides_str) or not mirror_index:
        return []
    try:
        peptides = ast.literal_eval(str(tryptic_peptides_str))
        if isinstance(peptides, str):
            peptides = [peptides]
    except (ValueError, SyntaxError):
        return []
    matches = []
    for pep in peptides:
        pep = pep.strip()
        if pep in mirror_index:
            matches.extend(mirror_index[pep]['plots'])
    matches.sort(key=lambda x: QUALITY_LEVELS.index(x['quality']))
    return matches


def create_ucsc_link(row, custom_session_id=None):
    """Create UCSC Genome Browser link."""
    def add_session_to_url(url, session_id):
        if session_id and 'genome.ucsc.edu' in url and CUSTOM_UCSC_SESSION_URL:
            position_match = re.search(r'position=([^&]+)', url)
            if position_match:
                position = position_match.group(1)
                sn = quote(CUSTOM_UCSC_SESSION)
                return (f"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38"
                        f"&hgS_doOtherUser=submit"
                        f"&hgS_otherUserName={CUSTOM_UCSC_USERNAME}"
                        f"&hgS_otherUserSessionName={sn}"
                        f"&position={position}")
            else:
                sn = quote(CUSTOM_UCSC_SESSION)
                return (f"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38"
                        f"&hgS_doOtherUser=submit"
                        f"&hgS_otherUserName={CUSTOM_UCSC_USERNAME}"
                        f"&hgS_otherUserSessionName={sn}")
        return url

    if 'CLICK_UCSC' in row and pd.notna(row['CLICK_UCSC']):
        click_ucsc = str(row['CLICK_UCSC'])
        if click_ucsc.startswith('=HYPERLINK('):
            match = re.search(r'=HYPERLINK\("([^"]+)"', click_ucsc)
            if match:
                return add_session_to_url(match.group(1), custom_session_id)
        elif click_ucsc.startswith('http'):
            return add_session_to_url(click_ucsc, custom_session_id)

    coords = None
    if 'smORF Coordinates' in row and pd.notna(row['smORF Coordinates']):
        coords = row['smORF Coordinates']
    elif 'genomic_coordinates' in row and pd.notna(row['genomic_coordinates']):
        coords = row['genomic_coordinates']

    if coords and coords.startswith('chr'):
        try:
            chrom, pos = coords.split(':')
            start, end = pos.split('-')
            position = f"{chrom}:{start}-{end}"
            if custom_session_id and CUSTOM_UCSC_SESSION_URL:
                sn = quote(CUSTOM_UCSC_SESSION)
                return (f"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38"
                        f"&hgS_doOtherUser=submit"
                        f"&hgS_otherUserName={CUSTOM_UCSC_USERNAME}"
                        f"&hgS_otherUserSessionName={sn}"
                        f"&position={position}")
            else:
                return f"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position={position}"
        except Exception:
            pass
    return None


# =============================================================================
# PAGE CONFIG (must be first Streamlit command)
# =============================================================================
st.set_page_config(
    page_title="Brain Microprotein Dashboard",
    page_icon="\U0001f9ec",
    layout="wide",
    # Filters are visible on arrival so the sidebar is discoverable rather than
    # hidden behind the toggle. Streamlit applies this on first load only; once
    # a visitor opens or closes the sidebar, their choice sticks for the session.
    # (The entry page puts nothing in the sidebar, so none is rendered there.)
    initial_sidebar_state="expanded",
)

# =============================================================================
# GLASSMORPHISM CSS
# =============================================================================
st.markdown(f"""
<style>
/* ── App gradient background ── */
.stApp, .stApp > header, [data-testid="stHeader"] {{
    background: linear-gradient(135deg,
        #1a202c 0%, #2d3748 25%, #1a202c 50%, #2d3748 75%, #1a202c 100%) !important;
}}

/* ── Always reserve the vertical scrollbar gutter. Prevents the layout
   "shake"/jitter loop when tall content (e.g. PROSIT mirror plots) toggles
   the scrollbar on/off, which resizes width-stretched images and re-toggles
   it endlessly. ── */
html, body {{
    overflow-y: scroll !important;
}}
/* Streamlit's real scroll container — keep the scrollbar gutter reserved and
   disable scroll anchoring so image loads don't nudge the viewport. */
[data-testid="stAppViewContainer"], [data-testid="stMain"], section.main {{
    scrollbar-gutter: stable !important;
    overflow-anchor: none !important;
}}
/* NOTE: scroll anchoring is NOT the lever for the entry-page nav. Re-enabling it
   here (on the scroller and every ancestor) was measured and did not compensate
   for late-loading figures — 584px of drift either way. The entry page instead
   reserves each figure's height up front so the layout never shifts; see
   _image_reserved. */

/* Cap any figure image so a fixed-width plot never forces a horizontal scrollbar. */
[data-testid="stImage"] img {{
    max-width: 100% !important;
    height: auto !important;
}}

/* ── Main container glass ── */
.main .block-container {{
    background: rgba(26, 32, 44, 0.82) !important;
    backdrop-filter: blur(12px) !important;
    border-radius: 16px !important;
    border: 1px solid rgba(116, 162, 183, 0.25) !important;
    box-shadow: 0 8px 32px rgba(0,0,0,0.35) !important;
    padding: 2rem !important;
    margin-top: 1rem !important;
}}

/* ── Sidebar gradient ── */
section[data-testid="stSidebar"] > div:first-child {{
    background: linear-gradient(180deg,
        rgba(26,32,44,0.96) 0%, rgba(45,55,72,0.96) 50%, rgba(26,32,44,0.96) 100%) !important;
    backdrop-filter: blur(10px) !important;
    border-right: 1px solid rgba(116,162,183,0.35) !important;
}}

/* ── Glass card (reusable) ── */
.glass-card {{
    background: rgba(255,255,255,0.06) !important;
    backdrop-filter: blur(14px) !important;
    border: 1px solid rgba(116,162,183,0.25) !important;
    border-radius: 14px !important;
    padding: 1.1rem 1.3rem !important;
    box-shadow: 0 4px 18px rgba(0,0,0,0.22) !important;
    margin-bottom: 0.6rem !important;
}}

.glass-card-section {{
    background: rgba(255,255,255,0.04) !important;
    backdrop-filter: blur(10px) !important;
    border: 1px solid rgba(116,162,183,0.18) !important;
    border-radius: 12px !important;
    padding: 1rem 1.2rem !important;
    margin-bottom: 0.75rem !important;
}}

/* ── ID card section header ── */
.id-section-header {{
    font-size: 0.95rem;
    font-weight: 600;
    margin-bottom: 0.6rem;
    padding-bottom: 0.35rem;
    border-bottom: 1px solid rgba(116,162,183,0.22);
    color: rgba(255,255,255,0.92) !important;
}}

/* ── ID card field ── */
.id-field-label {{
    font-size: 0.72rem;
    text-transform: uppercase;
    letter-spacing: 0.06em;
    color: rgba(255,255,255,0.50) !important;
    margin-bottom: 0.15rem;
}}
.id-field-value {{
    font-size: 1.05rem;
    font-weight: 500;
    color: #ffffff !important;
    margin-bottom: 0.7rem;
}}
/* Sequence-like values (coordinates, Kozak windows) — narrower and monospaced
   so they fit the grid cell without wrapping. */
.id-field-value.id-field-mono {{
    font-family: monospace;
    font-size: 0.82rem;
}}
.id-field-value a {{
    color: {COLORS['swiss_prot']} !important;
    text-decoration: none !important;
    border-bottom: 1px dotted rgba(116,162,183,0.6);
}}
.id-field-value a:hover {{
    color: #ffffff !important;
}}

/* ── Command-center header ── */
.cmd-header {{
    background: linear-gradient(135deg, {COLORS['swiss_prot']} 0%, {COLORS['unreviewed']} 100%);
    border-radius: 14px;
    padding: 1.2rem 1.6rem 1rem 1.6rem;
    margin-bottom: 1rem;
    box-shadow: 0 6px 24px rgba(0,0,0,0.25);
    position: relative;
    overflow: hidden;
}}
.cmd-header::before {{
    content: '';
    position: absolute;
    inset: 0;
    background: radial-gradient(ellipse at 20% 50%, rgba(255,255,255,0.12) 0%, transparent 70%);
    pointer-events: none;
}}
.cmd-top {{
    display: flex;
    justify-content: space-between;
    align-items: baseline;
    margin-bottom: 0.8rem;
    position: relative;
}}
.cmd-title {{
    font-size: 1.45rem;
    font-weight: 700;
    color: #ffffff !important;
    margin: 0;
    letter-spacing: -0.01em;
}}
.cmd-subtitle {{
    font-size: 0.78rem;
    color: rgba(255,255,255,0.65) !important;
    text-align: right;
    line-height: 1.4;
}}
.cmd-metrics {{
    display: grid;
    grid-template-columns: repeat(6, 1fr);
    gap: 0.55rem;
    position: relative;
}}
.cmd-stat {{
    background: rgba(0,0,0,0.18);
    backdrop-filter: blur(12px);
    border: 1px solid rgba(255,255,255,0.12);
    border-radius: 10px;
    padding: 0.6rem 0.5rem;
    text-align: center;
    transition: background 0.2s, transform 0.2s, box-shadow 0.2s;
}}
.cmd-stat:hover {{
    background: rgba(0,0,0,0.28);
    transform: translateY(-2px);
    box-shadow: 0 0 14px rgba(116,162,183,0.30), 0 4px 12px rgba(0,0,0,0.3);
}}
.cmd-stat .stat-label {{
    font-size: 0.65rem;
    text-transform: uppercase;
    letter-spacing: 0.06em;
    color: rgba(255,255,255,0.6);
    margin-bottom: 0.2rem;
}}
.cmd-stat .stat-val {{
    font-size: 1.45rem;
    font-weight: 700;
    color: #ffffff;
    line-height: 1.1;
}}
.cmd-stat .stat-sub {{
    font-size: 0.62rem;
    color: rgba(255,255,255,0.45);
    margin-top: 0.1rem;
}}
.cmd-stat.st-total .stat-val {{
    background: linear-gradient(135deg, #b8d8e8, #f5c4a1);
    -webkit-background-clip: text;
    -webkit-text-fill-color: transparent;
    background-clip: text;
}}
.cmd-stat.st-swiss {{
    border-color: rgba(116,162,183,0.45);
}}
.cmd-stat.st-swiss .stat-val {{
    color: #c5dce7;
}}
.cmd-stat.st-noncan {{
    border-color: rgba(237,134,81,0.45);
}}
.cmd-stat.st-noncan .stat-val {{
    color: #f5c4a1;
}}
.cmd-stat.st-swiss:hover {{
    box-shadow: 0 0 20px rgba(116,162,183,0.55), 0 4px 12px rgba(0,0,0,0.3) !important;
}}
.cmd-stat.st-noncan:hover {{
    box-shadow: 0 0 20px rgba(237,134,81,0.55), 0 4px 12px rgba(0,0,0,0.3) !important;
}}

/* ── Significance color helpers ── */
.sig-up   {{ color: #48bb78 !important; font-weight: 600; }}
.sig-down {{ color: #fc8181 !important; font-weight: 600; }}
.sig-ns   {{ color: rgba(255,255,255,0.40) !important; }}
.sig-yes  {{ color: #68d391 !important; font-weight: 600; }}
.sig-na   {{ color: rgba(255,255,255,0.30) !important; }}

/* ── Swiss-Prot / Unreviewed badges ── */
.badge-swiss {{
    display: inline-block;
    background: {COLORS['swiss_prot']};
    color: #fff !important;
    padding: 2px 10px;
    border-radius: 12px;
    font-size: 0.78rem;
    font-weight: 600;
}}
.badge-unreviewed {{
    display: inline-block;
    background: {COLORS['unreviewed']};
    color: #fff !important;
    padding: 2px 10px;
    border-radius: 12px;
    font-size: 0.78rem;
    font-weight: 600;
}}

/* ═══ Entry page (UniProt-style detail view) ═══════════════════════════════ */

/* Section picker. This used to be a list of <a href="#sec-…"> anchors, which
   worked locally but was dead on the Hugging Face Space: the huggingface.co
   wrapper embeds the app in a cross-origin iframe sized to the app's FULL
   content height, so the element the user actually scrolls is the *parent*
   huggingface.co document. A fragment link inside the iframe can only scroll
   scrollers inside that iframe, and same-origin policy forbids it touching the
   parent — measured on the live Space: parent scrollTop stayed 0 on every
   click, and the target landed at 443px instead of 72px.
   Bounding the app's own scroller does fix the landing (verified: exactly
   72px), but HF then grew the iframe to 7509px — it sizes from a postMessage,
   not from the document, so CSS cannot settle it. Selecting a section and
   re-rendering needs no scrolling at all, so it behaves identically in the
   wrapper, on the direct *.hf.space URL, on Streamlit Cloud and locally.

   Stick the element container that *holds* the picker, not the column or its
   vertical block: those stretch to the full height of the row (i.e. of the
   entry body), and a sticky box as tall as its scroll parent has nowhere to
   travel, so it never appears to stick. */
[data-testid="stElementContainer"]:has(.entry-nav-anchor) + [data-testid="stElementContainer"] {{
    position: sticky;
    top: 3.5rem;
    z-index: 5;
    background: rgba(255,255,255,0.04);
    backdrop-filter: blur(10px);
    border: 1px solid rgba(116,162,183,0.18);
    border-radius: 12px;
    padding: 0.6rem 0.5rem;
}}
/* Radio options styled as the old nav links. */
[data-testid="stElementContainer"]:has(.entry-nav-anchor) + [data-testid="stElementContainer"]
    [role="radiogroup"] > label {{
    display: block;
    padding: 0.30rem 0.55rem;
    margin-bottom: 0.1rem;
    border-left: 2px solid transparent;
    border-radius: 0 6px 6px 0;
    transition: background 0.12s ease, border-color 0.12s ease;
}}
[data-testid="stElementContainer"]:has(.entry-nav-anchor) + [data-testid="stElementContainer"]
    [role="radiogroup"] > label:hover {{
    background: rgba(116,162,183,0.16);
    border-left-color: {COLORS['swiss_prot']};
}}
[data-testid="stElementContainer"]:has(.entry-nav-anchor) + [data-testid="stElementContainer"]
    [role="radiogroup"] > label:has(input:checked) {{
    background: rgba(116,162,183,0.20);
    border-left-color: {COLORS['unreviewed']};
}}
[data-testid="stElementContainer"]:has(.entry-nav-anchor) + [data-testid="stElementContainer"]
    [role="radiogroup"] > label p {{
    font-size: 0.86rem !important;
    font-weight: 500 !important;
    color: rgba(255,255,255,0.82) !important;
}}

/* Anchor targets. scroll-margin-top keeps the section header clear of
   Streamlit's fixed toolbar when a nav link jumps to it. */
.entry-section {{
    scroll-margin-top: 4.5rem;
}}

/* Trailing runway for the last section (Sequence). An anchor can only scroll
   to the top of the viewport if enough content follows it; without this the
   final nav link lands short and the section never reaches the top, which
   reads as a broken link. 75vh clears it for the shortest realistic Sequence
   block without leaving a full blank screen below every entry. */
.entry-tail {{
    min-height: 75vh;
    pointer-events: none;
}}

/* Entry hero: gene / accession line + key-fact strip */
.entry-hero {{
    background: rgba(255,255,255,0.05);
    backdrop-filter: blur(14px);
    border: 1px solid rgba(116,162,183,0.25);
    border-left: 3px solid {COLORS['unreviewed']};
    border-radius: 12px;
    padding: 0.9rem 1.2rem 1rem 1.2rem;
    margin-bottom: 0.75rem;
}}
.entry-hero.is-swiss {{
    border-left-color: {COLORS['swiss_prot']};
}}
.entry-hero-title {{
    display: flex;
    align-items: center;
    gap: 0.7rem;
    flex-wrap: wrap;
    margin-bottom: 0.15rem;
}}
.entry-hero-title .entry-gene {{
    font-size: 1.6rem;
    font-weight: 700;
    color: #ffffff !important;
    line-height: 1.2;
}}
.entry-hero-sub {{
    font-family: monospace;
    font-size: 0.82rem;
    color: rgba(255,255,255,0.48) !important;
    margin-bottom: 0.8rem;
}}
.entry-facts {{
    display: flex;
    flex-wrap: wrap;
    gap: 0.4rem 2.2rem;
}}
.entry-facts > div {{ min-width: 7rem; }}

/* ── Password screen ── */
.login-card {{
    max-width: 420px;
    margin: 8vh auto;
    background: rgba(255,255,255,0.06);
    backdrop-filter: blur(16px);
    border: 1px solid rgba(116,162,183,0.3);
    border-radius: 18px;
    padding: 2.5rem 2rem;
    box-shadow: 0 8px 40px rgba(0,0,0,0.35);
    text-align: center;
}}
.login-card h2 {{
    background: linear-gradient(135deg, {COLORS['swiss_prot']}, {COLORS['unreviewed']});
    -webkit-background-clip: text;
    -webkit-text-fill-color: transparent;
    background-clip: text;
    font-size: 1.5rem;
    margin-bottom: 0.3rem;
}}
.login-card p {{
    color: rgba(255,255,255,0.6) !important;
    font-size: 0.85rem;
    margin-bottom: 1.2rem;
}}

/* ── Download button gradient ── */
.stDownloadButton > button {{
    background: linear-gradient(135deg, {COLORS['swiss_prot']}, {COLORS['unreviewed']}) !important;
    color: #fff !important;
    border: none !important;
    border-radius: 10px !important;
    font-weight: 600 !important;
    padding: 0.55rem 1.2rem !important;
    transition: all 0.3s ease !important;
}}
.stDownloadButton > button:hover {{
    background: linear-gradient(135deg, {COLORS['unreviewed']}, {COLORS['swiss_prot']}) !important;
    transform: translateY(-2px) !important;
    box-shadow: 0 6px 18px rgba(237,134,81,0.35) !important;
}}

/* ── Expander glass ── */
.streamlit-expanderHeader {{
    background: rgba(255,255,255,0.08) !important;
    backdrop-filter: blur(10px) !important;
    border: 1px solid rgba(116,162,183,0.2) !important;
    border-radius: 10px !important;
}}
.streamlit-expanderContent {{
    background: rgba(255,255,255,0.04) !important;
    backdrop-filter: blur(10px) !important;
    border: 1px solid rgba(116,162,183,0.15) !important;
    border-radius: 0 0 10px 10px !important;
}}

/* ── Metric containers (Streamlit native) ── */
[data-testid="metric-container"] {{
    background: rgba(255,255,255,0.06) !important;
    backdrop-filter: blur(8px) !important;
    border: 1px solid rgba(116,162,183,0.18) !important;
    border-radius: 10px !important;
}}

/* ── Dataframe glass wrapper ── */
.stDataFrame {{
    background: rgba(255,255,255,0.04) !important;
    backdrop-filter: blur(10px) !important;
    border-radius: 12px 12px 12px 12px !important;
    border: 1px solid rgba(116,162,183,0.18) !important;
}}
/* ── "How do I open an entry?" signpost above the table. The grid's row-select
      box is the only in-place control it offers, and unlabelled it reads as
      "select", so this points at it. ── */
.open-hint {{
    display: flex;
    align-items: center;
    gap: 0.6rem;
    background: rgba(237,134,81,0.10);
    border: 1px solid rgba(237,134,81,0.40);
    border-radius: 10px;
    padding: 0.5rem 0.9rem;
    margin: 0.2rem 0 0.5rem 0;
    font-size: 0.88rem;
    color: rgba(255,255,255,0.88) !important;
}}
.open-hint-arrow {{
    color: {COLORS['unreviewed']};
    font-size: 1rem;
    line-height: 1;
}}
.open-hint strong {{ color: #ffffff !important; }}

/* ── Row-click prompt below the table (the detail card it used to dock here
      now lives on its own entry page — see _render_entry_page) ── */
.detail-panel {{
    background: rgba(255,255,255,0.04);
    backdrop-filter: blur(14px);
    border: 1px solid rgba(116,162,183,0.18);
    border-top: 1px solid rgba(116,162,183,0.35);
    border-radius: 0 0 12px 12px;
    padding: 1rem 1.2rem 0.8rem 1.2rem;
    margin-top: -1rem;
    margin-bottom: 1rem;
    position: relative;
}}
.detail-panel::before {{
    content: '';
    position: absolute;
    top: 0; left: 5%; right: 5%; height: 1px;
    background: linear-gradient(90deg, transparent, rgba(116,162,183,0.5), transparent);
}}
.detail-header {{
    display: flex;
    align-items: center;
    gap: 0.7rem;
    margin-bottom: 0.6rem;
    padding-bottom: 0.5rem;
}}
.detail-header .detail-label {{
    font-size: 1.15rem;
    font-weight: 700;
    color: #ffffff;
}}
.detail-header .detail-badge {{
    margin-left: auto;
}}
.detail-prompt {{
    text-align: center;
    padding: 0.6rem 0;
    font-size: 0.85rem;
    color: rgba(255,255,255,0.4);
    font-style: italic;
}}

/* ── Inputs dark navy ── */
.stTextInput > div > div > input,
.stNumberInput input,
input[type="number"],
input[type="text"] {{
    background: #1a212d !important;
    color: #ffffff !important;
    border: 1px solid rgba(116,162,183,0.45) !important;
    border-radius: 8px !important;
}}
.stMultiSelect > div > div {{
    background: #1a212d !important;
    border: 1px solid rgba(116,162,183,0.45) !important;
    border-radius: 8px !important;
}}
.stButton > button {{
    background: #1a212d !important;
    color: #ffffff !important;
    border: 1px solid rgba(116,162,183,0.45) !important;
    border-radius: 8px !important;
    transition: all 0.3s ease !important;
}}
.stButton > button:hover {{
    background: #394254 !important;
    border: 1px solid {COLORS['swiss_prot']} !important;
    transform: translateY(-1px) !important;
    box-shadow: 0 4px 12px rgba(116,162,183,0.25) !important;
}}

/* ── Transparent blocks ── */
.element-container, div[data-testid="stVerticalBlock"],
div[data-testid="stHorizontalBlock"], section[data-testid="stSidebar"] > div {{
    background: transparent !important;
}}

/* ── Legend container ── */
.legend-container {{
    background: rgba(26,32,44,0.75) !important;
    border-radius: 12px !important;
    padding: 0.9rem !important;
    margin-top: 0.8rem !important;
    border: 1px solid rgba(116,162,183,0.3) !important;
    backdrop-filter: blur(8px) !important;
}}

/* ── Detail panel tabs ── */
[data-baseweb="tab-list"] {{
    background: rgba(255,255,255,0.06) !important;
    border-radius: 10px !important;
    padding: 4px !important;
    gap: 4px !important;
    border: 1px solid rgba(116,162,183,0.25) !important;
}}
[data-baseweb="tab"] {{
    background: transparent !important;
    color: rgba(255,255,255,0.65) !important;
    border-radius: 8px !important;
    font-size: 0.92rem !important;
    font-weight: 600 !important;
    letter-spacing: 0.02em !important;
    padding: 0.5rem 1.6rem !important;
    border: 1px solid transparent !important;
    transition: all 0.2s ease !important;
}}
[data-baseweb="tab"]:hover {{
    background: rgba(116,162,183,0.15) !important;
    color: rgba(255,255,255,0.9) !important;
    border-color: rgba(116,162,183,0.3) !important;
}}
[aria-selected="true"][data-baseweb="tab"] {{
    background: linear-gradient(135deg, rgba(116,162,183,0.30), rgba(237,134,81,0.20)) !important;
    color: #ffffff !important;
    border-color: rgba(116,162,183,0.5) !important;
    box-shadow: 0 2px 8px rgba(0,0,0,0.25) !important;
}}
[data-baseweb="tab-highlight"] {{
    display: none !important;
}}
[data-baseweb="tab-border"] {{
    display: none !important;
}}

/* ── Sidebar section headers ── */
.sidebar-section-header {{
    font-size: 0.68rem;
    font-weight: 700;
    text-transform: uppercase;
    letter-spacing: 0.10em;
    color: rgba(116,162,183,0.9) !important;
    border-left: 2px solid rgba(116,162,183,0.55);
    padding: 0.1rem 0 0.1rem 0.55rem;
    margin: 0.9rem 0 0.45rem 0;
    display: block;
}}

/* ── Table row hover & selection ── */
[data-testid="stDataFrame"] tr:hover td,
[data-testid="stDataFrame"] tr:hover th {{
    background: rgba(116,162,183,0.08) !important;
    transition: background 0.15s;
    cursor: pointer;
}}
[data-testid="stDataFrame"] tr[aria-selected="true"] td {{
    background: rgba(116,162,183,0.15) !important;
    border-left: 2px solid rgba(116,162,183,0.65) !important;
}}

/* ── Active-filter summary bar ── */
.filter-bar {{
    display: flex;
    flex-wrap: wrap;
    align-items: center;
    gap: 0.4rem;
    background: rgba(26,32,44,0.55);
    border: 1px solid rgba(116,162,183,0.22);
    border-radius: 12px;
    padding: 0.55rem 0.8rem;
    margin: -0.3rem 0 0.55rem 0;
}}
.filter-bar-title {{
    font-size: 0.7rem;
    text-transform: uppercase;
    letter-spacing: 0.07em;
    font-weight: 700;
    color: rgba(116,162,183,0.95);
    margin-right: 0.25rem;
}}
.filter-chip {{
    display: inline-flex;
    align-items: stretch;
    overflow: hidden;
    border-radius: 999px;
    border: 1px solid rgba(116,162,183,0.4);
    font-size: 0.76rem;
    line-height: 1.6;
}}
.filter-chip-key {{
    background: rgba(116,162,183,0.28);
    color: #d5e6ef !important;
    padding: 0 8px;
    font-weight: 600;
}}
.filter-chip-val {{
    background: rgba(237,134,81,0.16);
    color: #f5d3bd !important;
    padding: 0 9px;
    font-weight: 500;
}}
.filter-bar-none {{
    font-size: 0.8rem;
    color: rgba(255,255,255,0.55);
    font-style: italic;
}}
</style>
""", unsafe_allow_html=True)


# =============================================================================
# DATA LOADING — single unified source (the MASTER dataset)
# =============================================================================
# The dashboard derives every analysis view from one file,
# Code/data/microprotein_master.csv, by applying the canonical gold-standard
# filter and reproducing the column selections used by the
# Results/*/*_summary.csv generator scripts. The only inputs NOT contained in
# the master are the scRNA enrichment table (per gene x cluster) and the
# per-peptide tryptic / PROSIT-tier files, which are still read from disk.
MASTER_CSV = Path(__file__).resolve().parent.parent / "Code" / "data" / "microprotein_master.csv"


def _resolve_master_source(csv_path=MASTER_CSV):
    """Return a path pandas can read for the master table.

    The 88 MB `microprotein_master.csv` is git-ignored, so it is absent when
    the repo is checked out on Streamlit Community Cloud. The compressed
    `microprotein_master.zip` (a single-member archive) IS committed, and
    pandas reads it directly via compression inference. Prefer the plain CSV
    when present (local dev); otherwise fall back to the shipped zip.
    """
    p = Path(csv_path)
    if p.exists():
        return p
    # Second name is the legacy (misspelled) filename, kept for backward compat.
    for _zip_name in ("microprotein_master.zip", "micrprotein_master.zip"):
        z = p.with_name(_zip_name)
        if z.exists():
            return z
    return p  # let the caller surface a clear FileNotFoundError


# UniProt curation annotation score (1–5) keyed by accession. Sourced from the
# UniProt proteome export (columns Entry, Annotation). Attached only to
# microproteins that ARE genuine UniProt entries (reviewed Swiss-Prot-MP or TrEMBL).
UNIPROT_ANNOT_TSV = Path(__file__).resolve().parent.parent / "Code" / "data" / "uniprotkb_proteome_UP000005640_2026_07_13.tsv"


@st.cache_data(show_spinner=False)
def load_uniprot_annotation_scores(_tsv_path=str(UNIPROT_ANNOT_TSV)):
    """Map UniProt accession -> annotation score (1–5) from the proteome export.

    Returns an empty dict if the TSV is unavailable so the dashboard degrades
    gracefully (the annotation-score column simply stays blank).
    """
    p = Path(_tsv_path)
    if not p.exists():
        return {}
    try:
        ann = pd.read_csv(p, sep='\t', usecols=['Entry', 'Annotation'])
    except Exception:
        return {}
    scores = pd.to_numeric(ann['Annotation'], errors='coerce')
    return {str(e): float(s) for e, s in zip(ann['Entry'], scores) if pd.notna(s)}


# Per-microprotein rescue for the N-terminus filter: does at least one of its
# N-terminal (aa<=2) tryptic peptides carry 1-2 amino-acid substitutions vs.
# its matched UniProt isoform? See
# Code/Microprotein_annotation_summary/compute_nterm_peptide_substitutions.py.
NTERM_SUBSTITUTIONS_CSV = Path(__file__).resolve().parent.parent / "Code" / "data" / "blast_nterm_peptide_substitutions.csv"


@st.cache_data(show_spinner=False)
def load_nterm_substitutions(_csv_path=str(NTERM_SUBSTITUTIONS_CSV)):
    """Map microprotein sequence -> (has_rescue, detail rows).

    `has_rescue` is True when any N-terminal peptide has a comparable
    alignment (coverage>=80%) to its matched isoform with 1-2 substitutions --
    direct evidence that peptide isn't identical to the canonical protein, so
    it can't be a bare Met-excision fragment of it. `detail` carries the
    per-peptide rows for display in the microprotein detail view.

    Returns an empty dict if the file is unavailable so the N-terminus filter
    degrades gracefully to its unadjusted (stricter) behavior.
    """
    p = Path(_csv_path)
    if not p.exists():
        return {}
    try:
        sub = pd.read_csv(p)
    except Exception:
        return {}
    result = {}
    for seq, group in sub.groupby('sequence'):
        rescues = group[group['comparable'].astype(bool)
                         & group['mismatch_count'].between(1, 2)]
        result[seq] = {
            'has_rescue': not rescues.empty,
            'detail': group.to_dict('records'),
        }
    return result


@st.cache_data(show_spinner=False)
def load_decay_classes(_tsv_path=str(DECAY_CLASS_TSV)):
    """Map master gene_id -> ORF-rules label ('NMD · mORF Intact', ...).

    Keyed on the assignments TSV's `smorf_id`, which is the same identifier as
    the master's `gene_id` (the master has no `orf_id` column). The join is
    exact and 1:1; the ~867 `_alt_initiation_N` proteoforms in the TSV have no
    master counterpart and are deliberately left unmatched -- collapsing them
    onto their parent ORF would assign conflicting classes for 9 parents.

    The TSV's `category` spans six values: the NMD x main-ORF-disruption 2x2
    (dark_green/light_green/gold/firebrick) plus the two non-stop-decay
    classes (nsd_disrupted/nsd_intact), which the generator scores per smORF
    ahead of that 2x2 and which therefore override it. Categories absent from
    DECAY_CLASS_FROM_CATEGORY are dropped and fall through to 'NA', so a new
    upstream category must be added there or it silently reads as unassessed.

    Returns an empty dict if the file is unavailable so the column degrades
    gracefully to all-'NA' (e.g. on Streamlit Cloud / the HF Space).
    """
    p = Path(_tsv_path)
    if not p.exists():
        return {}
    try:
        a = pd.read_csv(p, sep='\t', usecols=['smorf_id', 'category'])
    except Exception:
        return {}
    return {r.smorf_id: DECAY_CLASS_FROM_CATEGORY[r.category]
            for r in a.itertuples()
            if r.category in DECAY_CLASS_FROM_CATEGORY}


# Best-tier ranking used to collapse the bracketed PROSIT `Confidence` list.
_CONFIDENCE_RANK = {'Strong': 0, 'Moderate': 1, 'Weak': 2, 'Insufficient': 3}


def _best_confidence(val):
    """Pick the best (lowest-rank) PROSIT tier. Accepts either a plain scalar tier
    string (current master) or a bracketed list literal of tiers (legacy master)."""
    if pd.isna(val):
        return None
    text = str(val).strip()
    if text in _CONFIDENCE_RANK:        # plain scalar tier (current master format)
        return text
    try:
        items = ast.literal_eval(text)
    except (ValueError, SyntaxError):
        return None
    if not isinstance(items, list):
        return None
    valid = [x for x in items if x in _CONFIDENCE_RANK]
    if not valid:
        return None
    return min(valid, key=lambda x: _CONFIDENCE_RANK[x])


@st.cache_data(show_spinner=False)
def load_and_filter_master(_csv_path=str(MASTER_CSV)):
    """Apply the canonical gold-standard filter to the MASTER dataset.

    This is an inlined copy of
    Code/gold_standard_filtering_criteria.load_and_filter_master — the
    authoritative filtering pipeline shared by all Results-summary generator
    scripts. Keep the two in sync if the criteria ever change.
    """
    df = pd.read_csv(_resolve_master_source(_csv_path), low_memory=False)

    # Step 1: drop contaminant / non-human entries
    if "gene_symbol" in df.columns:
        df = df[~df["gene_symbol"].str.contains(r"SHEEP|Peptide|Horse|BOVIN",
                                                case=False, na=False)]

    # Step 2: keep relevant databases
    df = df[df["Database"].str.contains("TrEMBL|Salk|Swiss-Prot", na=False)].copy()

    # Step 3: remap Database labels (TrEMBL is Salk-derived; Swiss-Prot -> Swiss-Prot-MP)
    df.loc[df["Database"] == "TrEMBL", "smorf_type"] = "TrEMBL"
    df["Database"] = df["Database"].replace({
        "TrEMBL": "Salk",
        "Swiss-Prot": "Swiss-Prot-MP",
    })

    # Step 4: evidence flags
    ribocode_pass = df["RiboCode"].astype(str).str.strip().str.upper() == "TRUE"  # noqa: E712
    df["has_RiboSAM"] = (
        ribocode_pass
        & df["shortstop_label"].isin(["SAM-Secreted", "SAM-Intracellular"])
        & (~df["smorf_type"].isin(["iORF"]))
        & (~df["smorf_type"].str.contains("iso", case=False, na=False))
    )
    df["has_LooseRiboSAM"] = (
        ribocode_pass
        & df["shortstop_label"].isin(["SAM-Secreted", "SAM-Intracellular"])
        & (~df["smorf_type"].isin(["Iso"]))
    )
    df["DDA_evidence"] = df["total_unique_spectral_counts"] > 0
    df["DIA_evidence"] = (
        df["has_LooseRiboSAM"]
        & df["Global.PG.Q.Value"].notna()
        & (df["Global.PG.Q.Value"] <= 0.01)
        & (df["Proteotypic"] == 1)
    )
    df["DIA_expression"] = (
        df["Global.PG.Q.Value"].notna()
        & (df["Global.PG.Q.Value"] <= 0.01)
        & (df["Proteotypic"] == 1)
    )
    # DIA proteomics is intentionally excluded from this dashboard: only DDA
    # mass-spec counts as MS evidence (DIA_evidence/DIA_expression columns are
    # retained for reference but are NOT used for inclusion or labeling).
    df["has_MS"] = df["DDA_evidence"]

    # Coerce numeric-evidence columns that may contain stray non-numeric
    # tokens (e.g. "FALSE") in the upstream master
    for _col in ("RP3_Default", "total_razor_spectral_counts"):
        if _col in df.columns:
            df[_col] = pd.to_numeric(df[_col], errors="coerce").fillna(0)

    # Step 5: microproteins only
    df = df[df["protein_class_length"] == "Microprotein"].copy()

    # Step 6: split by database, apply evidence filters, dedup by sequence
    mp_swiss = (
        df[df["Database"] == "Swiss-Prot-MP"]
        .loc[lambda d:
            (d["RP3_Default"] > 0)
            | (d["total_razor_spectral_counts"] > 0)]
        .drop_duplicates(subset="sequence")
    )
    mp_salk = (
        df[(df["Database"] == "Salk")
           & (df["has_RiboSAM"] | df["DDA_evidence"])]
        .sort_values("has_MS", ascending=False)
        .drop_duplicates(subset="sequence")
    )
    mp = pd.concat([mp_swiss, mp_salk], ignore_index=True)

    # Step 7: attach the smORF decay class. Must happen HERE, while gene_id is
    # still around -- the per-analysis merge downstream selects an explicit
    # column allowlist keyed on `sequence` and drops gene_id entirely.
    _decay_map = load_decay_classes()
    if _decay_map and "gene_id" in mp.columns:
        mp["NMD_Decay_Class"] = mp["gene_id"].map(_decay_map).fillna("NA")
    else:
        mp["NMD_Decay_Class"] = "NA"
    return mp


# --- Per-analysis derivations (mirror the Results-summary generator scripts) ---
def _derive_annotation(mp):
    """Salk discovery summary (Annotation Status, MS Evidence Type, DDA Grade)."""
    df = mp[mp["Database"] == "Salk"].copy()
    # Evidence Type matches Figure 7 Panel J: RiboCode-ShortStop = RiboSAM
    # translation support without DDA mass-spec detection (has_RiboSAM & ~DDA_evidence);
    # everything else (DDA-detected) is MS. DIA is not counted in this dashboard.
    df["Annotation Status"] = (df["has_RiboSAM"] & ~df["DDA_evidence"]).map(
        {True: "RiboCode-ShortStop", False: "MS"})
    df["MS_type"] = df["DDA_evidence"].map({True: "DDA", False: "No MS"})
    df["DDA Grade"] = df["Confidence"].apply(_best_confidence)
    df.loc[df["DDA Grade"].isna() & (df["Annotation Status"] == "RiboCode-ShortStop"), "DDA Grade"] = "No MS"
    df.loc[df["DDA Grade"].isna() & (df["Annotation Status"] == "MS"), "DDA Grade"] = "No PROSIT"
    # Kozak context strength + non-AUG start codon assessment — Salk-only columns
    # (Code/data/microprotein_master.csv nonATG_*/kozak_* fields), passed through
    # unrenamed so extract_unified_fields can pick them up by suffix.
    _nonatg_kozak_cols = [c for c in (
        "kozak_full_kozak_strength", "kozak_downstream_kozak_strength",
        "kozak_kozak_class_canonical", "kozak_kozak_window",
        "nonATG_annotated_start_codon", "nonATG_annotated_start_is_initiator",
        "nonATG_has_supported_initiation_site",
        "nonATG_annotated_context_strength", "nonATG_n_initiation_sites",
        "nonATG_no_site_reason", "nonATG_site_type",
        "nonATG_initiation_codon_position", "nonATG_initiation_codon",
        "nonATG_cognate_status", "nonATG_initiation_tier_name",
        "nonATG_initiator_aa", "nonATG_predicted_sequence",
    ) if c in df.columns]
    summary = df[[
        "sequence", "CLICK_UCSC", "genomic_coordinates", "smorf_type", "gene_name",
        "protein_length", "Annotation Status", "MS_type", "DDA Grade", "mean_phylocsf",
        *_nonatg_kozak_cols,
    ]].rename(columns={
        "genomic_coordinates": "smORF Coordinates",
        "smorf_type": "smORF Class",
        "gene_name": "Parent Gene",
        "protein_length": "Microprotein Length",
        "sequence": "Microprotein Sequence",
        "MS_type": "MS Evidence Type",
        "mean_phylocsf": "PhyloCSF Score",
    })
    summary["_starts_M"] = summary["Microprotein Sequence"].str.startswith("M")
    summary = summary.sort_values(by="_starts_M", ascending=False).drop(columns="_starts_M")
    return summary


def _derive_proteomics(mp):
    return mp[[
        "sequence", "CLICK_UCSC",
        "TMT_log2fc_50pct_missing", "TMT_t_statistic_50pct_missing", "TMT_df_50pct_missing",
        "TMT_conf_low_50pct_missing", "TMT_conf_high_50pct_missing", "TMT_cohens_d_50pct_missing",
        "TMT_pvalue_50pct_missing", "TMT_qvalue_50pct_missing",
        "TMT_log2fc_0pct_missing", "TMT_t_statistic_0pct_missing", "TMT_df_0pct_missing",
        "TMT_conf_low_0pct_missing", "TMT_conf_high_0pct_missing", "TMT_cohens_d_0pct_missing",
        "TMT_pvalue_0pct_missing", "TMT_qvalue_0pct_missing",
        "rate_control", "rate_ad", "Database",
        "protein_class_length", "gene_symbol", "gene_name", "protein_length",
        "start_codon", "smorf_type", "total_razor_spectral_counts",
        "total_unique_spectral_counts", "mean_phylocsf",
        "blastp_uniprot_accession_match", "microprotein_percentage_match",
        "blastp_alignment_length", "evalue", "blastp_bit",
    ]].copy()


def _derive_rp3(mp):
    return mp[mp["Database"].isin(["Salk", "Swiss-Prot-MP"])][[
        "sequence", "CLICK_UCSC", "RP3_Default", "RP3_MM_Amb", "RP3_Amb", "RP3_MM",
        "RiboCode", "Database", "gene_name", "gene_symbol", "protein_length",
        "start_codon", "smorf_type", "total_unique_spectral_counts",
        "total_razor_spectral_counts", "mean_phylocsf",
    ]].copy()


def _derive_shortread(mp):
    return mp[[
        "sequence", "CLICK_UCSC", "gene_symbol", "genomic_coordinates",
        "rosmapRNA_baseMean", "rosmapRNA_log2FoldChange", "rosmapRNA_pvalue",
        "rosmapRNA_padj", "rosmapRNA_non_smorf_hit", "rosmapRNA_body",
        "ROSMAP_BulkRNAseq_CPM", "correlation_mainORF_nonAD_rosmap",
        "correlation_mainORF_AD_rosmap", "rosmap_lrt_additive_p",
        "rosmap_lrt_interaction_p", "msbbRNA_baseMean", "msbbRNA_log2FoldChange",
        "msbbRNA_pvalue", "msbbRNA_padj", "msbbRNA_non_smorf_hit", "msbbRNA_body",
        "MSBB_BulkRNAseq_CPM", "correlation_mainORF_nonAD_msbb",
        "correlation_mainORF_AD_msbb", "Database", "gene_name", "protein_length",
        "start_codon", "smorf_type", "total_razor_spectral_counts",
        "total_unique_spectral_counts", "mean_phylocsf",
    ]].copy()


def _derive_longread(mp):
    df = mp[[
        "sequence", "CLICK_UCSC", "nanopore_baseMean", "nanopore_log2FoldChange",
        "nanopore_pvalue", "nanopore_padj", "Database", "gene_name", "gene_symbol",
        "protein_length", "start_codon", "smorf_type", "total_razor_spectral_counts",
        "total_unique_spectral_counts",
    ]].copy()
    return df[(df["Database"] == "Salk") | df["nanopore_baseMean"].notna()]


def _derive_shortstop(mp):
    df = mp[mp["Database"] == "Salk"].copy()
    df["Annotation Status"] = (df["has_RiboSAM"] & ~df["DDA_evidence"]).map(
        {True: "RiboCode-ShortStop", False: "MS"})
    return df.loc[df["shortstop_label"].notna(), [
        "sequence", "CLICK_UCSC", "gene_symbol", "gene_name", "smorf_type",
        "shortstop_label", "shortstop_score", "Annotation Status", "mean_phylocsf",
    ]].rename(columns={
        "sequence": "Microprotein Sequence",
        "gene_symbol": "smORF ID",
        "gene_name": "Gene Body (Name)",
        "smorf_type": "smORF Class",
        "shortstop_label": "ShortStop Label",
        "shortstop_score": "ShortStop Score",
        "mean_phylocsf": "PhyloCSF Score",
    })


def load_analysis_results():
    """Build every analysis view from the single master file.

    Returns an ordered dict {analysis_name: {"df": DataFrame, "type": str}}.
    All views except scRNA Enrichment are derived in-memory from the
    gold-standard-filtered master; scRNA is read from its own table because it
    is not contained in the master dataset.
    """
    mp = load_and_filter_master()

    scrna_path = SCRNA_SUMMARY_FILE
    scrna_df = None
    if INCLUDE_SCRNA and scrna_path.exists():
        try:
            scrna_df = pd.read_csv(scrna_path, low_memory=False)
        except Exception as e:
            st.warning(f"Could not load scRNA Enrichment: {e}")

    analyses = {
        "Annotation Summary": {"df": _derive_annotation(mp), "type": "unreviewed_only"},
        "Proteomics (TMT)": {"df": _derive_proteomics(mp), "type": "mixed"},
        "Proteomics + RiboSeq (RP3)": {"df": _derive_rp3(mp), "type": "mixed"},
        "Short-Read RNA in AD": {"df": _derive_shortread(mp), "type": "mixed"},
        "Long-Read RNA in AD": {"df": _derive_longread(mp), "type": "unreviewed_only"},
        "scRNA Enrichment": {"df": scrna_df, "type": "scrna"} if INCLUDE_SCRNA else None,
        "ShortStop Classification": {"df": _derive_shortstop(mp), "type": "unreviewed_only"},
    }
    # Drop conditional/None entries and any view that failed to build.
    return {k: v for k, v in analyses.items()
            if v is not None and v.get("df") is not None}


@st.cache_data(show_spinner=False)
def load_and_merge_all_data():
    """Merge every analysis view by sequence, anchored on the master dataset.

    The tryptic-peptide, PROSIT and N-terminal-acetylation columns come straight
    from the master (peptide_sequence/start/end -> Tryptic_*, Confidence for the
    per-peptide PROSIT tiers, Nt_acetyl_* for the Nt-acetylated subset), so no
    separate peptide CSV or tier-ID files are read.
    """
    analysis_files = load_analysis_results()

    # Base frame: tryptic peptides + PROSIT confidence + Nt-acetylation straight
    # from the master. The Nt_acetyl_* columns are ';'-delimited and parallel to
    # each other (one entry per Nt-acetylated tryptic peptide), but are a *subset*
    # of Tryptic_Peptides — they are not positionally aligned with it.
    mp = load_and_filter_master()
    _base_cols = [c for c in ['sequence', 'peptide_sequence', 'start', 'end',
                             'Confidence', 'NMD_Decay_Class'] if c in mp.columns]
    _acetyl_cols = [c for c in ('Nt_acetyl_tryptic_peptide', 'Nt_acetyl_N_PSMs',
                                'Nt_acetyl_PSM_fraction') if c in mp.columns]
    master_df = mp[_base_cols + _acetyl_cols].rename(
        columns={
            'peptide_sequence': 'Tryptic_Peptides',
            'start': 'Tryptic_Start_Positions',
            'end': 'Tryptic_End_Positions',
            'Nt_acetyl_tryptic_peptide': 'Nt_Acetyl_Peptides',
            'Nt_acetyl_N_PSMs': 'Nt_Acetyl_N_PSMs',
            'Nt_acetyl_PSM_fraction': 'Nt_Acetyl_PSM_Fraction',
        }
    ).copy()
    master_df['Tryptic_Peptides_present'] = master_df['Tryptic_Peptides'].notna()
    if 'Nt_Acetyl_Peptides' in master_df.columns:
        master_df['Nt_Acetylated'] = (
            master_df['Nt_Acetyl_Peptides'].notna()
            & master_df['Nt_Acetyl_Peptides'].astype(str).str.strip().ne('')
        )
    else:
        master_df['Nt_Acetylated'] = False

    for analysis_name, info in analysis_files.items():
        df = info.get('df')
        if df is not None and not df.empty:
            try:
                df = df.copy()
                if 'Microprotein Sequence' in df.columns:
                    df = df.rename(columns={'Microprotein Sequence': 'sequence'})
                elif 'sequence' not in df.columns:
                    continue
                df[f'{analysis_name}_present'] = True
                prefix = analysis_name.replace(' ', '_').replace('(', '').replace(')', '').replace('+', '_')
                cols_to_rename = {}
                for col in df.columns:
                    if col not in ['sequence', f'{analysis_name}_present']:
                        cols_to_rename[col] = f"{prefix}_{col}"
                df = df.rename(columns=cols_to_rename)
                if master_df.empty:
                    master_df = df.copy()
                else:
                    master_df = pd.merge(master_df, df, on='sequence', how='outer', suffixes=('', '_dup'))
                    dup_cols = [c for c in master_df.columns if c.endswith('_dup')]
                    master_df = master_df.drop(columns=dup_cols)
            except Exception as e:
                st.warning(f"Could not load {analysis_name}: {e}")

    presence_cols = [c for c in master_df.columns if c.endswith('_present')]
    for col in presence_cols:
        master_df[col] = master_df[col].fillna(False)

    if not master_df.empty and 'sequence' in master_df.columns:
        data_cols = [c for c in master_df.columns if c != 'sequence' and not c.endswith('_present')]
        master_df['_completeness'] = master_df[data_cols].notna().sum(axis=1)
        master_df = master_df.sort_values(['sequence', '_completeness'], ascending=[True, False])
        master_df = master_df.drop_duplicates(subset=['sequence'], keep='first')
        master_df = master_df.drop('_completeness', axis=1)

    ann_col = [c for c in master_df.columns
              if 'annotation' in c.lower() and 'summary' in c.lower() and c.endswith('_present')]
    if ann_col:
        ann_mask = master_df[ann_col[0]] == True

        # Keep original unreviewed-anchored set, OR include Swiss-Prot microproteins.
        db_cols = [c for c in master_df.columns if 'database' in c.lower()]
        if db_cols:
            swiss_mask = master_df[db_cols].apply(
                lambda row: any('swiss-prot' in str(v).lower() for v in row if pd.notna(v)),
                axis=1
            )
        else:
            swiss_mask = pd.Series(False, index=master_df.index)

        # Microprotein criterion: sequence length <= 151 aa.
        seq_len = master_df['sequence'].astype(str).str.len() if 'sequence' in master_df.columns else pd.Series(999, index=master_df.index)
        microprotein_mask = seq_len <= 151

        # Only include Swiss-Prot entries that appear in the Proteomics CSV.
        # The 752 Swiss-Prot-MP entries added after the first submission have no
        # MS evidence and should not be shown in the dashboard.
        prot_col = next((c for c in master_df.columns if c == 'Proteomics (TMT)_present'), None)
        if prot_col:
            prot_present_mask = master_df[prot_col].fillna(False).astype(bool)
        else:
            prot_present_mask = pd.Series(True, index=master_df.index)

        master_df = master_df[ann_mask | (swiss_mask & microprotein_mask & prot_present_mask)]

    return master_df


def _has_non_nterm_peptide(pos_str):
    """True if any tryptic peptide starts at aa >= 3.

    A peptide starting at aa 1 or 2 is exactly what Met excision of the *parent*
    ORF would produce, so it is not on its own evidence for the microprotein;
    a peptide starting at aa 3 or later is evidence independent of that artefact.
    `pos_str` is the master's Python-list-literal of start positions.
    """
    if not pos_str or pd.isna(pos_str):
        return False  # no tryptic data is not evidence of anything
    try:
        positions = ast.literal_eval(str(pos_str))
    except (ValueError, SyntaxError):
        return False
    if isinstance(positions, (int, float)):
        positions = [positions]
    try:
        return any(int(p) >= 3 for p in positions)
    except (TypeError, ValueError):
        return False


def _non_nterm_series(df, index=None):
    """Boolean `Has_Non_Nterm_Peptide` for `df`, derived on the fly if absent.

    Normally a plain column lookup — the column is precomputed once per data
    load. The fallback matters when a long-running process holds a dataframe
    cached before the column existed: st.cache_data keys on the decorated
    function's own bytecode, so editing a helper it calls does not invalidate
    it. Recomputing there costs a parse but keeps the filter honest; silently
    returning "no filter" would make a stale cache look like a broken widget.
    """
    if 'Has_Non_Nterm_Peptide' in df.columns:
        return df['Has_Non_Nterm_Peptide'].fillna(False).astype(bool)
    if 'Tryptic_Start_Positions' in df.columns:
        return df['Tryptic_Start_Positions'].map(_has_non_nterm_peptide).astype(bool)
    return None


def extract_unified_fields(master_df):
    """Extract and unify key fields from the merged dataset."""
    ud = master_df.copy()

    # Parent Gene
    cols = [c for c in master_df.columns
            if ('parent' in c.lower() and 'gene' in c.lower())
            or c.lower().endswith('gene_name')
            or c.lower().endswith('gene_symbol')
            or (c.lower().endswith('_gene') and 'id' not in c.lower())
            or c.lower().endswith('gene body (name)')]
    if cols:
        ud['Parent_Gene'] = master_df[cols].bfill(axis=1).iloc[:, 0]

    # smORF Class
    cols = [c for c in master_df.columns
            if ('smorf' in c.lower() and ('class' in c.lower() or 'type' in c.lower()))
            and 'id' not in c.lower().split('smorf')[1]
            and 'coordinates' not in c.lower()]
    if cols:
        ud['smORF_Class'] = master_df[cols].bfill(axis=1).iloc[:, 0]
    else:
        ud['smORF_Class'] = pd.NA

    # Swiss-Prot-MP entries have no smorf_type — assign class so they appear in group filters.
    # (After the merge, the bare 'Database' column doesn't exist; use the prefixed cols.)
    _db_cols_raw = [c for c in master_df.columns if 'database' in c.lower()]
    if _db_cols_raw:
        _is_swiss_mp = master_df[_db_cols_raw].apply(
            lambda row: any('swiss-prot' in str(v).lower() for v in row if pd.notna(v)), axis=1
        )
        ud.loc[ud['smORF_Class'].isna() & _is_swiss_mp, 'smORF_Class'] = 'Swiss-Prot-MP'

    # Display label (no remapping currently; hook for future overrides)
    ud['smORF_Display'] = ud['smORF_Class'].map(
        lambda v: SMORF_DISPLAY_LABEL.get(str(v), v) if pd.notna(v) else v
    )
    # Parent group (Upstream, Downstream, TrEMBL/AltORF, etc.)
    ud['smORF_Group'] = ud['smORF_Class'].map(
        lambda v: SMORF_CHILD_TO_PARENT.get(str(v), str(v)) if pd.notna(v) else v
    )

    # smORF decay class (NMD × main-ORF disruption). Already joined on gene_id
    # back in load_and_filter_master and carried through the merge; this only
    # normalizes rows the outer merge introduced.
    # NB: the column is named NMD_* on purpose -- the smORF_Class discovery
    # above pattern-matches any column containing 'smorf' + 'class'/'type', so
    # a name like 'smorf_decay_class' would be silently bfilled into
    # smORF_Class and shadow smorf_type.
    if 'NMD_Decay_Class' in ud.columns:
        ud['NMD_Decay_Class'] = ud['NMD_Decay_Class'].fillna('NA')
    else:
        ud['NMD_Decay_Class'] = 'NA'

    # ShortStop Label
    cols = [c for c in master_df.columns if 'shortstop' in c.lower() and 'label' in c.lower()]
    if cols:
        ud['ShortStop_Label'] = master_df[cols].bfill(axis=1).iloc[:, 0]

    # Annotation Status
    cols = [c for c in master_df.columns if 'annotation' in c.lower() and 'status' in c.lower()]
    if cols:
        ud['Annotation_Status'] = master_df[cols].bfill(axis=1).iloc[:, 0]

    # Unique Spectral Counts
    cols = [c for c in master_df.columns if 'unique_spectral_counts' in c.lower()]
    if cols:
        ud['Unique_Spectral_Counts'] = pd.to_numeric(
            master_df[cols].bfill(axis=1).iloc[:, 0], errors='coerce')

    # Razor Spectral Counts
    cols = [c for c in master_df.columns if 'razor_spectral_counts' in c.lower()]
    if cols:
        ud['Razor_Spectral_Counts'] = pd.to_numeric(
            master_df[cols].bfill(axis=1).iloc[:, 0], errors='coerce')

    # UCSC Link
    cols = [c for c in master_df.columns if 'ucsc' in c.lower() or 'CLICK_UCSC' in c]
    if cols:
        ud['UCSC_Link'] = master_df[cols].bfill(axis=1).iloc[:, 0]

    # Protein Length
    cols = [c for c in master_df.columns if 'length' in c.lower() and 'class' not in c.lower()]
    if cols:
        ud['Protein_Length'] = pd.to_numeric(
            master_df[cols].bfill(axis=1).iloc[:, 0], errors='coerce')
    # Fill missing Protein_Length from sequence length
    if 'sequence' in ud.columns:
        missing = ud['Protein_Length'].isna()
        ud.loc[missing, 'Protein_Length'] = ud.loc[missing, 'sequence'].astype(str).str.len()

    # Start Codon (the plain "ATG"/"nonATG" call). Excludes 'nonatg' columns —
    # nonATG_annotated_start_codon also matches 'start'+'codon' but holds the raw
    # triplet (e.g. "GCC"), not the coarse ATG/nonATG label.
    cols = [c for c in master_df.columns
            if 'start' in c.lower() and 'codon' in c.lower() and 'nonatg' not in c.lower()]
    if cols:
        ud['Start_Codon'] = master_df[cols].bfill(axis=1).iloc[:, 0]

    # ShortStop Score
    cols = [c for c in master_df.columns if 'shortstop' in c.lower() and 'score' in c.lower()]
    if cols:
        ud['ShortStop_Score'] = pd.to_numeric(
            master_df[cols].bfill(axis=1).iloc[:, 0], errors='coerce')

    # PhyloCSF Score
    cols = [c for c in master_df.columns if 'phylocsf' in c.lower() or 'mean_phylocsf' in c.lower()]
    if cols:
        ud['PhyloCSF_Score'] = pd.to_numeric(
            master_df[cols].bfill(axis=1).iloc[:, 0], errors='coerce')
        min_val = ud['PhyloCSF_Score'].min()
        if pd.isna(min_val):
            min_val = -1000.0
        # Fill NaN with min to push missing-data rows to the bottom of score-based sorts
        ud['PhyloCSF_Score'] = ud['PhyloCSF_Score'].fillna(min_val)

    # Tryptic peptides (+ the Nt-acetylated subset of them)
    for col in ['Tryptic_Peptides', 'Tryptic_Protein_ID', 'Tryptic_Start_Positions', 'Tryptic_End_Positions',
                'Nt_Acetyl_Peptides', 'Nt_Acetyl_N_PSMs', 'Nt_Acetyl_PSM_Fraction']:
        if col in master_df.columns:
            ud[col] = master_df[col]

    # N-terminal acetylation flag + summary numbers. The master stores ';'-delimited
    # parallel lists (one entry per Nt-acetylated peptide); collapse them to a total
    # PSM count and the best (highest) per-peptide PSM fraction so both are sortable.
    ud['Nt_Acetylated'] = (
        master_df['Nt_Acetylated'].fillna(False).astype(bool)
        if 'Nt_Acetylated' in master_df.columns
        else pd.Series(False, index=ud.index)
    )

    # N-terminal-peptide substitution rescue (see load_nterm_substitutions) --
    # False for the overwhelming majority of rows, which have no comparable
    # 1-2-substitution N-terminal peptide at all.
    _nterm_sub_map = load_nterm_substitutions()
    ud['Nterm_Substitution_Rescue'] = (
        ud['sequence'].map(lambda s: _nterm_sub_map.get(s, {}).get('has_rescue', False))
        if 'sequence' in ud.columns
        else pd.Series(False, index=ud.index)
    )
    ud['Nterm_Substitution_Detail'] = (
        ud['sequence'].map(lambda s: _nterm_sub_map.get(s, {}).get('detail', []))
        if 'sequence' in ud.columns
        else pd.Series([[]] * len(ud), index=ud.index)
    )

    # Precomputed here, inside the cached loader, so the N-terminus filter is a
    # column lookup rather than a per-rerun re-parse of the position lists once
    # per cross-filtered facet. See _has_non_nterm_peptide for the rule.
    ud['Has_Non_Nterm_Peptide'] = (
        ud['Tryptic_Start_Positions'].map(_has_non_nterm_peptide)
        if 'Tryptic_Start_Positions' in ud.columns
        else pd.Series(False, index=ud.index)
    )

    def _sum_semicolon(val):
        if val is None or (isinstance(val, float) and pd.isna(val)):
            return np.nan
        parts = pd.to_numeric(pd.Series(str(val).split(';')).str.strip(), errors='coerce')
        return parts.sum() if parts.notna().any() else np.nan

    def _max_semicolon(val):
        if val is None or (isinstance(val, float) and pd.isna(val)):
            return np.nan
        parts = pd.to_numeric(pd.Series(str(val).split(';')).str.strip(), errors='coerce')
        return parts.max() if parts.notna().any() else np.nan

    if 'Nt_Acetyl_N_PSMs' in ud.columns:
        ud['Nt_Acetyl_Total_PSMs'] = ud['Nt_Acetyl_N_PSMs'].map(_sum_semicolon)
    if 'Nt_Acetyl_PSM_Fraction' in ud.columns:
        ud['Nt_Acetyl_Max_Fraction'] = ud['Nt_Acetyl_PSM_Fraction'].map(_max_semicolon)

    # --- TMT proteomics stats ---
    for src, dest in [
        ('TMT_log2fc_50pct_missing',       'TMT_log2fc_50pct'),
        ('TMT_t_statistic_50pct_missing',  'TMT_t_statistic_50pct'),
        ('TMT_df_50pct_missing',           'TMT_df_50pct'),
        ('TMT_conf_low_50pct_missing',     'TMT_conf_low_50pct'),
        ('TMT_conf_high_50pct_missing',    'TMT_conf_high_50pct'),
        ('TMT_cohens_d_50pct_missing',     'TMT_cohens_d_50pct'),
        ('TMT_pvalue_50pct_missing',       'TMT_pvalue_50pct'),
        ('TMT_qvalue_50pct_missing',       'TMT_qvalue_50pct'),
        ('TMT_log2fc_0pct_missing',        'TMT_log2fc_0pct'),
        ('TMT_t_statistic_0pct_missing',   'TMT_t_statistic_0pct'),
        ('TMT_df_0pct_missing',            'TMT_df_0pct'),
        ('TMT_conf_low_0pct_missing',      'TMT_conf_low_0pct'),
        ('TMT_conf_high_0pct_missing',     'TMT_conf_high_0pct'),
        ('TMT_cohens_d_0pct_missing',      'TMT_cohens_d_0pct'),
        ('TMT_pvalue_0pct_missing',        'TMT_pvalue_0pct'),
        ('TMT_qvalue_0pct_missing',        'TMT_qvalue_0pct'),
        ('rate_control', 'TMT_rate_control'),
        ('rate_ad',      'TMT_rate_ad'),
    ]:
        cols = [c for c in master_df.columns if src in c]
        if cols:
            ud[dest] = pd.to_numeric(master_df[cols[0]], errors='coerce')

    # --- Short-Read RNA (ROSMAP) stats ---
    for src, dest in [
        ('rosmapRNA_log2FoldChange', 'ROSMAP_log2FC'),
        ('rosmapRNA_padj', 'ROSMAP_padj'),
        ('rosmapRNA_pvalue', 'ROSMAP_pvalue'),
        ('rosmapRNA_baseMean', 'ROSMAP_baseMean'),
        ('ROSMAP_BulkRNAseq_CPM', 'ROSMAP_CPM'),
        ('correlation_mainORF_nonAD_rosmap', 'Corr_MainORF_NonAD'),
        ('correlation_mainORF_AD_rosmap', 'Corr_MainORF_AD'),
        ('rosmap_lrt_additive_p', 'RNA_LRT_Add_P'),
        ('rosmap_lrt_interaction_p', 'RNA_LRT_Int_P'),
    ]:
        cols = [c for c in master_df.columns if src in c]
        if cols:
            ud[dest] = pd.to_numeric(master_df[cols].bfill(axis=1).iloc[:, 0], errors='coerce')

    # --- Short-Read RNA (MSBB) stats ---
    for src, dest in [
        ('msbbRNA_log2FoldChange', 'MSBB_log2FC'),
        ('msbbRNA_padj', 'MSBB_padj'),
        ('msbbRNA_pvalue', 'MSBB_pvalue'),
        ('msbbRNA_baseMean', 'MSBB_baseMean'),
        ('MSBB_BulkRNAseq_CPM', 'MSBB_CPM'),
        ('correlation_mainORF_nonAD_msbb', 'Corr_MainORF_NonAD_MSBB'),
        ('correlation_mainORF_AD_msbb', 'Corr_MainORF_AD_MSBB'),
    ]:
        cols = [c for c in master_df.columns if src in c]
        if cols:
            ud[dest] = pd.to_numeric(master_df[cols].bfill(axis=1).iloc[:, 0], errors='coerce')

    # --- Long-Read RNA (Nanopore) stats ---
    for src, dest in [
        ('nanopore_log2FoldChange', 'Nanopore_log2FC'),
        ('nanopore_padj', 'Nanopore_padj'),
        ('nanopore_pvalue', 'Nanopore_pvalue'),
        ('nanopore_baseMean', 'Nanopore_baseMean'),
    ]:
        cols = [c for c in master_df.columns if src in c]
        if cols:
            ud[dest] = pd.to_numeric(master_df[cols].bfill(axis=1).iloc[:, 0], errors='coerce')

    # --- RP3 / Ribo-Seq columns ---
    for src, dest in [
        ('RP3_Default', 'RP3_Default'),
        ('RP3_MM_Amb', 'RP3_MM_Amb'),
        ('RP3_Amb', 'RP3_Amb'),
        ('RP3_MM', 'RP3_MM'),
        ('RiboCode', 'RiboCode'),
    ]:
        cols = [c for c in master_df.columns if c.endswith(src)]
        if cols:
            ud[dest] = master_df[cols].bfill(axis=1).iloc[:, 0]

    # --- BLASTp homology vs UniProt ---
    cols = [c for c in master_df.columns if 'blastp_uniprot_accession_match' in c]
    if cols:
        ud['BLAST_UniProt_Match'] = master_df[cols].bfill(axis=1).iloc[:, 0]
    for src, dest in [
        ('microprotein_percentage_match', 'BLAST_Pct_Match'),
        ('blastp_alignment_length', 'BLAST_Aln_Length'),
        ('evalue', 'BLAST_Evalue'),
        ('blastp_bit', 'BLAST_Bit'),
    ]:
        cols = [c for c in master_df.columns if src in c]
        if cols:
            ud[dest] = pd.to_numeric(master_df[cols].bfill(axis=1).iloc[:, 0], errors='coerce')

    # --- scRNA enrichment ---
    if INCLUDE_SCRNA:
        for src, dest in [
            ('scRNA_Enrichment_logFC', 'scRNA_logFC'),
            ('scRNA_Enrichment_p_adj.glb', 'scRNA_padj'),
            ('scRNA_Enrichment_celltype', 'scRNA_celltype'),
            ('scRNA_Enrichment_cell_type_general', 'scRNA_cell_type_general'),
        ]:
            cols = [c for c in master_df.columns if c == src]
            if cols:
                if 'logFC' in src or 'padj' in src:
                    ud[dest] = pd.to_numeric(master_df[cols[0]], errors='coerce')
                else:
                    ud[dest] = master_df[cols[0]]

    # --- Additional annotation fields for ID card ---
    cols = [c for c in master_df.columns if 'smorf' in c.lower() and 'coordinates' in c.lower()]
    if cols:
        ud['smORF_Coordinates'] = master_df[cols].bfill(axis=1).iloc[:, 0]

    cols = [c for c in master_df.columns if 'ms' in c.lower() and 'evidence' in c.lower()]
    if cols:
        ud['MS_Evidence_Type'] = master_df[cols].bfill(axis=1).iloc[:, 0]

    cols = [c for c in master_df.columns if 'dda' in c.lower() and 'grade' in c.lower()]
    if cols:
        ud['DDA_Grade'] = master_df[cols].bfill(axis=1).iloc[:, 0]

    # --- Kozak context + non-AUG start codon assessment (Salk smORFs only) ---
    for src, dest in [
        ('kozak_full_kozak_strength', 'Kozak_Strength'),
        ('kozak_downstream_kozak_strength', 'Kozak_Downstream_Strength'),
        ('kozak_kozak_class_canonical', 'Kozak_Class_Canonical'),
        ('kozak_kozak_window', 'Kozak_Window'),
        ('nonATG_annotated_start_codon', 'NonATG_Annotated_Codon'),
        ('nonATG_annotated_start_is_initiator', 'NonATG_Is_Initiator'),
        ('nonATG_has_supported_initiation_site', 'NonATG_Has_Supported_Site'),
        ('nonATG_annotated_context_strength', 'NonATG_Context_Strength'),
        ('nonATG_no_site_reason', 'NonATG_No_Site_Reason'),
        # Stringified-list columns (one entry per candidate site), parsed on
        # demand in _precompute_display_columns / the detail card — same
        # convention as Tryptic_Peptides.
        ('nonATG_site_type', 'NonATG_Site_Type'),
        ('nonATG_initiation_codon_position', 'NonATG_Codon_Position'),
        ('nonATG_initiation_codon', 'NonATG_Codon'),
        ('nonATG_cognate_status', 'NonATG_Cognate_Status'),
        ('nonATG_initiation_tier_name', 'NonATG_Tier_Name'),
        ('nonATG_initiator_aa', 'NonATG_Initiator_AA'),
        ('nonATG_predicted_sequence', 'NonATG_Predicted_Sequence'),
    ]:
        cols = [c for c in master_df.columns if c.endswith(src)]
        if cols:
            ud[dest] = master_df[cols].bfill(axis=1).iloc[:, 0]

    cols = [c for c in master_df.columns if c.endswith('nonATG_n_initiation_sites')]
    if cols:
        ud['NonATG_N_Sites'] = pd.to_numeric(master_df[cols].bfill(axis=1).iloc[:, 0], errors='coerce')

    cols = [c for c in master_df.columns if 'database' in c.lower()]
    if cols:
        ud['Database'] = master_df[cols].bfill(axis=1).iloc[:, 0]
        # Keep all database labels seen across merged sources for robust filtering.
        ud['Database_All'] = master_df[cols].apply(
            lambda row: '|'.join(sorted({
                str(v).strip() for v in row
                if pd.notna(v) and str(v).strip() and str(v).strip().lower() != 'none'
            })),
            axis=1
        )

    # --- UniProt annotation score (1–5) for reviewed/TrEMBL microproteins ---
    # Attach UniProt's curation annotation score ONLY where the microprotein IS a
    # genuine UniProt entry (reviewed Swiss-Prot-MP or TrEMBL); novel Salk ORFs stay
    # blank. The microprotein's own accession = its 100%-identity BLASTp self-hit
    # (BLAST_UniProt_Match with BLAST_Pct_Match == 100), so a homolog hit never leaks
    # a neighbour's score onto a novel ORF.
    ud['UniProt_Annotation_Score'] = np.nan
    if 'BLAST_UniProt_Match' in ud.columns:
        _score_map = load_uniprot_annotation_scores()
        if _score_map:
            _is_reviewed = ud.get('Database', pd.Series('', index=ud.index)).astype(str).str.contains(
                'swiss-prot', case=False, na=False)
            _is_trembl = ud.get('smORF_Class', pd.Series('', index=ud.index)).astype(str).eq('TrEMBL')
            _pct = pd.to_numeric(ud.get('BLAST_Pct_Match', pd.Series(np.nan, index=ud.index)), errors='coerce')
            _self_acc = ud['BLAST_UniProt_Match'].where((_is_reviewed | _is_trembl) & (_pct == 100))
            # UniProt annotation scores are keyed on canonical accessions, so strip any
            # isoform suffix (e.g. Q5BLP8-2 -> Q5BLP8) before mapping.
            _self_acc = _self_acc.astype(str).str.replace(r'-\d+$', '', regex=True)
            ud['UniProt_Annotation_Score'] = pd.to_numeric(
                _self_acc.map(_score_map), errors='coerce')

    # Significance indicators (kept for sidebar filter + metrics)
    ud['TMT_Significant'] = ud.get('TMT_qvalue_50pct', pd.Series(dtype=float)).fillna(1.0) < 0.2
    ud['TMT_Highly_Significant'] = ud.get('TMT_qvalue_50pct', pd.Series(dtype=float)).fillna(1.0) < 0.05
    ud['ROSMAP_Significant'] = ud.get('ROSMAP_padj', pd.Series(dtype=float)).fillna(1.0) < 0.2
    ud['ROSMAP_Highly_Significant'] = ud.get('ROSMAP_padj', pd.Series(dtype=float)).fillna(1.0) < 0.05

    return ud


# =============================================================================
# DISK-PARQUET CACHE FOR MERGED + EXTRACTED DATAFRAME
# =============================================================================
# The full merge+extract pipeline filters the master dataset and runs many
# bfill/regex passes (~3-8 s cold). Caching the final unified_df to a single
# parquet on disk makes subsequent cold starts <0.5 s. Cache key is the mtime
# fingerprint of every source file — auto-invalidates when any input changes.

def _source_csv_paths():
    base = Path(__file__).resolve().parent
    return [
        # Single unified source for every master-derived view (microprotein set,
        # tryptic peptides, PROSIT confidence, and all per-analysis columns).
        # Falls back to the committed .zip when the git-ignored .csv is absent
        # (e.g. on Streamlit Cloud) so the fingerprint still tracks the master.
        _resolve_master_source(),
        # Not contained in the master — read directly.
        base / "scRNA_Enrichment" / "scRNA_Enrichment_summary.csv",
        # UniProt annotation-score export (Entry -> Annotation); lives in Code/data
        # so the parquet fingerprint invalidates when it changes.
        base.parent / "Code" / "data" / "uniprotkb_proteome_UP000005640_2026_07_13.tsv",
        # smORF decay-class assignments (NMD × main-ORF disruption); also in
        # Code/data so the parquet fingerprint invalidates when it changes.
        DECAY_CLASS_TSV,
    ]


@st.cache_data(show_spinner=False)
def _load_unified_df_with_disk_cache():
    """Load+merge+extract, caching the final dataframe to parquet on disk."""
    cache_dir = Path(__file__).parent / "_cache"
    cache_dir.mkdir(exist_ok=True)
    parquet_path = cache_dir / "unified_df.parquet"
    fp_path = cache_dir / "unified_df.fingerprint"

    # Include mirror_index size in fingerprint so _Spectra_Quality column
    # invalidates when figure libraries change.
    try:
        mirror_index = build_mirror_plot_index()
    except Exception:
        mirror_index = {}
    mirror_fp = f"mirror:{len(mirror_index)}"

    fp = "|".join(
        f"{p.name}:{p.stat().st_mtime_ns}:{p.stat().st_size}"
        for p in _source_csv_paths() if p.exists()
    ) + "|" + mirror_fp + "|" + f"code:{Path(__file__).stat().st_mtime_ns}"

    if parquet_path.exists() and fp_path.exists():
        try:
            if fp_path.read_text() == fp:
                return pd.read_parquet(parquet_path)
        except Exception:
            pass

    master_df = load_and_merge_all_data()
    unified_df = extract_unified_fields(master_df)
    unified_df = _precompute_display_columns(unified_df, mirror_index)
    try:
        unified_df.to_parquet(parquet_path, index=False)
        fp_path.write_text(fp)
    except Exception:
        pass  # parquet failures shouldn't block the app
    return unified_df


def _precompute_display_columns(ud, mirror_index):
    """Pre-compute heavy per-row display columns once, so row-click reruns
    can simply slice instead of running .apply() across all 6.5k rows."""
    # _Spectra_Quality — derive the best PROSIT tier directly from the master's
    # per-peptide `Confidence` list (Strong > Moderate > Weak > Insufficient).
    # This is the same authoritative source used for the annotation `DDA Grade`.
    # For the few rows without a Confidence list (e.g. some Swiss-Prot entries)
    # fall back to the mirror-plot peptide lookup, then default to 'No PROSIT'.
    if 'Confidence' in ud.columns:
        ud['_Spectra_Quality'] = ud['Confidence'].apply(_best_confidence)
        if 'Tryptic_Peptides' in ud.columns and mirror_index:
            missing_mask = ud['_Spectra_Quality'].isna()
            if missing_mask.any():
                ud.loc[missing_mask, '_Spectra_Quality'] = (
                    ud.loc[missing_mask, 'Tryptic_Peptides'].apply(
                        lambda x: get_spectra_quality(x, mirror_index)
                    )
                )
        ud['_Spectra_Quality'] = ud['_Spectra_Quality'].fillna('No PROSIT')
    elif 'Tryptic_Peptides' in ud.columns:
        ud['_Spectra_Quality'] = ud['Tryptic_Peptides'].apply(
            lambda x: get_spectra_quality(x, mirror_index)
        )
    else:
        ud['_Spectra_Quality'] = 'No PROSIT'

    # _Tryptic_Display: joined readable string from list-literal
    def _fmt_peps(val):
        if val is None or (isinstance(val, float) and pd.isna(val)) or val == '':
            return ''
        try:
            peps = ast.literal_eval(str(val))
            if isinstance(peps, str):
                peps = [peps]
            return ' \u00b7 '.join(str(p).strip() for p in peps)
        except (ValueError, SyntaxError):
            return str(val)
    if 'Tryptic_Peptides' in ud.columns:
        ud['_Tryptic_Display'] = ud['Tryptic_Peptides'].map(_fmt_peps)

    # NonATG_Best_Tier: best (lowest-rank) candidate tier for a non-AUG-annotated
    # microprotein, so it can be sorted/filtered as a scalar in the main table.
    _NONATG_TIER_RANK = {
        'Well-Established Near-Cognate': 0,
        'Weak Near-Cognate': 1,
        'Non-Near-Cognate, Strongest Proteomics Evidence': 2,
        'Non-Near-Cognate, Single-Instance Proteomics Evidence': 3,
        'Non-Near-Cognate, Inferred (not directly reported as a TIS)': 4,
    }
    def _best_tier(val):
        if val is None or (isinstance(val, float) and pd.isna(val)):
            return pd.NA
        try:
            tiers = ast.literal_eval(str(val))
        except (ValueError, SyntaxError):
            return pd.NA
        if isinstance(tiers, str):
            tiers = [tiers]
        ranked = sorted((t for t in tiers if t in _NONATG_TIER_RANK), key=_NONATG_TIER_RANK.get)
        return ranked[0] if ranked else pd.NA
    if 'NonATG_Tier_Name' in ud.columns:
        ud['NonATG_Best_Tier'] = ud['NonATG_Tier_Name'].map(_best_tier)

    # _NonATG_Candidates_Display: readable "codon@position" list for the main
    # table's Alt-Initiation view, same joined-string convention as _Tryptic_Display.
    def _fmt_candidates(codon_val, pos_val):
        if codon_val is None or (isinstance(codon_val, float) and pd.isna(codon_val)):
            return ''
        try:
            codons = ast.literal_eval(str(codon_val))
            positions = ast.literal_eval(str(pos_val)) if pd.notna(pos_val) else []
        except (ValueError, SyntaxError):
            return ''
        if isinstance(codons, str):
            codons = [codons]
        if isinstance(positions, (str, int)):
            positions = [positions]
        return ' · '.join(
            f"{c}@{p}" if p is not None else str(c)
            for c, p in zip(codons, positions or [None] * len(codons))
        )
    if 'NonATG_Codon' in ud.columns and 'NonATG_Codon_Position' in ud.columns:
        ud['_NonATG_Candidates_Display'] = [
            _fmt_candidates(c, p) for c, p in zip(ud['NonATG_Codon'], ud['NonATG_Codon_Position'])
        ]

    # _UCSC_HTML: pre-built UCSC link URL
    if 'UCSC_Link' in ud.columns:
        ud['_UCSC_HTML'] = ud['UCSC_Link'].apply(
            lambda v: create_ucsc_link({'CLICK_UCSC': v}, CUSTOM_UCSC_SESSION)
            if pd.notna(v) else None
        )

    # _DB_Emoji: vectorized
    db_col = 'Database_All' if 'Database_All' in ud.columns else ('Database' if 'Database' in ud.columns else None)
    if db_col:
        _swiss = ud[db_col].fillna('').str.lower().str.contains('swiss', na=False)
        ud['_DB_Emoji'] = np.where(_swiss, '\U0001F535 Swiss-Prot', '\U0001F7E0 Unreviewed')

    # _TMT_Tier: best (most-stringent) TMT significance tier reached per row.
    #   Tier 1: q(50%) < 0.05  ·  Tier 2: q(50%) < 0.2
    #   Tier 3: q(0%)  < 0.05  ·  Tier 4: q(0%)  < 0.2
    q50 = pd.to_numeric(ud.get('TMT_qvalue_50pct', pd.Series(index=ud.index, dtype=float)), errors='coerce')
    q0 = pd.to_numeric(ud.get('TMT_qvalue_0pct', pd.Series(index=ud.index, dtype=float)), errors='coerce')
    tmt_tier = pd.Series(pd.NA, index=ud.index, dtype='object')
    tmt_tier[q0 < 0.20] = 'Tier 4'
    tmt_tier[q0 < 0.05] = 'Tier 3'
    tmt_tier[q50 < 0.20] = 'Tier 2'
    tmt_tier[q50 < 0.05] = 'Tier 1'
    ud['_TMT_Tier'] = tmt_tier

    # _RNA_Sig: ROSMAP RNA-seq significance (Significant = padj<0.05, Exploratory = padj<0.2).
    padj = pd.to_numeric(ud.get('ROSMAP_padj', pd.Series(index=ud.index, dtype=float)), errors='coerce')
    rna_sig = pd.Series(pd.NA, index=ud.index, dtype='object')
    rna_sig[padj < 0.20] = 'Exploratory'
    rna_sig[padj < 0.05] = 'Significant'
    ud['_RNA_Sig'] = rna_sig

    # _UniProt_Annotation_Display: gold-star rendering of the 1–5 score.
    if 'UniProt_Annotation_Score' in ud.columns:
        ud['_UniProt_Annotation_Display'] = ud['UniProt_Annotation_Score'].map(_annotation_stars)

    return ud


# =============================================================================
# ACTIVE-FILTER SUMMARY
# =============================================================================
def _render_active_filter_bar(active_chips, result_count, total_count):
    """Render a compact bar listing the filters currently narrowing the view.
    ``active_chips`` is a list of (label, value) tuples for every sidebar/preset
    control that is active. (The baseline-criteria explainer that used to live
    here as "What am I looking at?" moved to _render_results_table as the
    "Explain the column/variables to me" expander, next to the column-view tabs.)
    """
    if active_chips:
        chips_html = "".join(
            "<span class='filter-chip'>"
            f"<span class='filter-chip-key'>{html.escape(str(k))}</span>"
            f"<span class='filter-chip-val'>{html.escape(str(v))}</span></span>"
            for k, v in active_chips
        )
        st.markdown(
            "<div class='filter-bar'>"
            f"<span class='filter-bar-title'>Active filters ({len(active_chips)})</span>"
            f"{chips_html}"
            "<span class='filter-bar-title' style='margin-left:auto; margin-right:0; "
            f"color:rgba(255,255,255,0.5);'>{result_count:,} of {total_count:,} shown</span>"
            "</div>",
            unsafe_allow_html=True,
        )
    else:
        st.markdown(
            "<div class='filter-bar'><span class='filter-bar-title'>Active filters</span>"
            f"<span class='filter-bar-none'>None \u2014 showing all {total_count:,} "
            "microproteins in the baseline set. Use the Quick-start presets or the sidebar "
            "to refine.</span></div>",
            unsafe_allow_html=True,
        )
    st.caption(
        f"Every view starts from a fixed gold-standard baseline of {total_count:,} brain "
        "microproteins (reviewed Swiss-Prot + unreviewed Salk/TrEMBL), applied before any "
        "sidebar filters — see Figure 1 of the manuscript for the evidence criteria."
    )


# =============================================================================
# SIDEBAR FACET HELPERS (UniProt-style checkbox facets with live counts)
# =============================================================================
def _is_swiss_mask(df):
    """Boolean mask: True for reviewed Swiss-Prot rows. Single source of truth
    for the Database_All-else-Database fallback, previously duplicated with
    slightly different fallback behavior across several call sites."""
    col = 'Database_All' if 'Database_All' in df.columns else ('Database' if 'Database' in df.columns else None)
    if col is None:
        return pd.Series(False, index=df.index)
    return df[col].fillna('').astype(str).str.lower().str.contains('swiss', na=False)


# N-terminus filter modes. A tryptic peptide starting at aa 1 or 2 is what Met
# excision of the *parent* ORF would produce, so on its own it is not evidence
# that the microprotein is real. These options select the microproteins whose
# evidence survives discounting that artefact.
NTERM_MODE_NON_NTERM = "Non-N-terminal peptides only"
NTERM_MODE_ACETYL_OR_NON_NTERM = "Nt-acetylated, substitution-distinct, or non-N-terminal"
# Facet domain. Nothing ticked means no filtering; ticked boxes are OR'd, like
# every other checkbox facet. NON_NTERM is a strict subset of
# ACETYL_OR_NON_NTERM, so ticking both is the same as ticking the second —
# exactly how the nested significance tiers (FDR<0.05 within FDR<0.2) behave.
NTERM_MODES = [NTERM_MODE_NON_NTERM, NTERM_MODE_ACETYL_OR_NON_NTERM]


def _count_checkbox(key, label, n, help=None):
    """A facet checkbox whose live count lives in the label, not the key.

    The `value=` seeding is load-bearing and the same trick as in
    _render_facet_checkboxes: Streamlit derives widget identity partly from the
    label, so a count baked into it mints a new widget on every recount, and a
    new widget falls back to its default — wiping the tick even though `key` is
    stable. Seeding the default from stored state makes that fallback a no-op.
    """
    st.checkbox(
        f"{label} ({n:,})", key=key,
        value=bool(st.session_state.get(key, False)),
        help=help,
    )


def _nterm_mask(df, mode):
    """Row mask for an N-terminus mode; None means no filtering.

    NON_NTERM keeps rows with a tryptic peptide starting at aa >= 3 — evidence
    beyond the Met-excision artefact. ACETYL_OR_NON_NTERM also readmits two
    other kinds of rows: Nt-acetylated (Nt-acetylation is co-translational and
    marks a genuine protein N-terminus, so an acetylated aa-1/2 peptide is real
    evidence), and rows whose N-terminal peptide carries 1-2 amino-acid
    substitutions vs. its matched UniProt isoform (Nterm_Substitution_Rescue,
    see compute_nterm_peptide_substitutions.py) — a peptide that differs from
    the canonical protein's sequence at that position cannot be a bare
    Met-excision fragment of it.

    Rows with no tryptic peptides at all fail both modes — absence of peptide
    data is not evidence.
    """
    if mode not in (NTERM_MODE_NON_NTERM, NTERM_MODE_ACETYL_OR_NON_NTERM):
        return None
    mask = _non_nterm_series(df)
    if mask is None:
        return None  # no peptide data at all in this frame
    if mode == NTERM_MODE_ACETYL_OR_NON_NTERM:
        if 'Nt_Acetylated' in df.columns:
            mask = mask | df['Nt_Acetylated'].fillna(False).astype(bool)
        if 'Nterm_Substitution_Rescue' in df.columns:
            mask = mask | df['Nterm_Substitution_Rescue'].fillna(False).astype(bool)
    return mask


def _render_facet_checkboxes(series, key_prefix, order=None, help=None, label_fn=None,
                             domain=None):
    """Render one st.checkbox per distinct value in `series`, each labeled with
    a live count. Sorted by count descending unless an explicit `order` list is
    given (for facets with a real ordinal/curated meaning, e.g. quality tiers).

    `series` is the CROSS-FILTERED view for this facet — the data narrowed by
    every OTHER active filter but not by this facet's own selection. That is
    what makes these counts react to boxes checked anywhere else (including in
    facets rendered further down the sidebar), while a facet's own boxes never
    zero out their own counts.

    Returns the list of currently-checked values — the same list shape
    st.multiselect used to return, so downstream `.isin(selected)` filtering
    code needs no changes, only the widget call site does.

    `label_fn`, if given, formats each value for display only (e.g. to
    capitalize "strong" -> "Strong"); the underlying value used for the
    checkbox key and the returned `selected` list is unaffected.

    `domain` (or `order`) is the facet's full value list, and every value in it
    gets a box even at count 0. A facet must never drop an option as other
    filters narrow the data: a vanished box cannot be ticked, so the only way
    back to it would be to clear the filters that hid it. Zero-count boxes stay
    clickable and simply read "(0)".
    """
    counts = series.value_counts(dropna=True)
    # Curated `order` is rendered whole; otherwise sort by count descending and
    # append any declared domain value the current view happens not to contain.
    values = list(order) if order else counts.index.tolist()
    for v in (domain or []):
        if v not in values:
            values.append(v)
    selected = []
    for i, v in enumerate(values):
        n = int(counts.get(v, 0))
        label = label_fn(v) if label_fn else v
        _key = f"{key_prefix}_{v}"
        # `value=` is seeded from the box's own stored state, and that is load-
        # bearing: Streamlit derives widget identity partly from the LABEL, so
        # the live count baked into it mints a new widget every time any other
        # facet changes the count — and a new widget silently falls back to its
        # default, wiping the tick even though `key` is stable. Seeding the
        # default with the stored value makes that fallback a no-op. Without
        # this, each click clears every other facet's selections.
        checked = st.checkbox(
            f"{label} ({n:,})", key=_key,
            value=bool(st.session_state.get(_key, False)),
            help=(help if i == 0 else None),
        )
        if checked:
            selected.append(v)
    return selected


# =============================================================================
# MAIN APP — SINGLE-PAGE LAYOUT
# =============================================================================
# ── Filter-state persistence across the entry page ──────────────────────────
# Streamlit drops widget state for any widget that did not render on a rerun,
# and the entry page returns from main() before the sidebar draws. Without a
# shadow copy, opening a microprotein silently clears every filter and the
# search box, so clicking "Back to results" lands on an unfiltered table —
# measured: a checked facet reads True during the entry render and False on the
# way back. The shadow key is a plain (non-widget) key, so nothing culls it.
#
# Not recoverable this way: the table's column sort, which st.dataframe never
# reports back to Python, and scroll position. Both are lost by design.
_PERSIST_SHADOW = '_filter_state_shadow'
_PERSIST_EXTRA = ('hero_search',)


def _persisted_filter_keys():
    """Widget keys whose values should survive an entry-page detour."""
    return [
        k for k in st.session_state
        if isinstance(k, str) and (k.startswith('f_') or k in _PERSIST_EXTRA)
    ]


def _snapshot_filter_state():
    """Mirror live filter widgets into the shadow key. Call after they render."""
    st.session_state[_PERSIST_SHADOW] = {
        k: st.session_state[k] for k in _persisted_filter_keys()
    }


def _restore_filter_state():
    """Re-seed filter widgets from the shadow key.

    Must run *before* the widgets instantiate — assigning to a widget key after
    its widget exists on the same run raises. Only fills keys Streamlit culled,
    so it never clobbers a selection the user is actively changing.
    """
    saved = st.session_state.get(_PERSIST_SHADOW)
    if not saved:
        return
    for k, v in saved.items():
        if k not in st.session_state:
            st.session_state[k] = v


def main():
    # ── Gradient banner placeholder (rendered after data loads) ──
    header_slot = st.empty()

    # Load data (parquet-cached on disk; auto-rebuilds when source CSVs change)
    with st.spinner("Loading microprotein database..."):
        try:
            unified_df = _load_unified_df_with_disk_cache()
            if 'Annotation_Status' in unified_df.columns:
                # Keep Swiss-Prot rows (which have no Annotation_Status) alongside annotated rows.
                # Only include Swiss-Prot entries that have proteomics evidence (i.e., appear in
                # the Proteomics CSV). The 752 Swiss-Prot-MP entries added after first submission
                # have no MS data and should not appear in the default view.
                has_annotation = ~(unified_df['Annotation_Status'].isna() | (unified_df['Annotation_Status'] == 'None'))
                is_swiss = unified_df.get('Database_All', unified_df.get('Database', pd.Series('', index=unified_df.index))).fillna('').str.lower().str.contains('swiss')
                _prot_present_col = next((c for c in unified_df.columns if 'proteomics' in c.lower() and c.endswith('_present')), None)
                if _prot_present_col:
                    has_proteomics = unified_df[_prot_present_col].fillna(False).astype(bool)
                    is_swiss = is_swiss & has_proteomics
                unified_df = unified_df[has_annotation | is_swiss]
            if unified_df.empty:
                st.error("No data could be loaded.")
                return
        except Exception as e:
            st.error(f"Error loading data: {e}")
            return

    # ── Pre-computed indexes (mirror plots, expression profiles, coords) ──
    mirror_index = build_mirror_plot_index()
    expression_index = build_expression_profile_index()
    seq_to_coords = build_seq_to_coords_index()
    # _Spectra_Quality is now precomputed inside the cached load (see
    # _precompute_display_columns); fall back only if the column is missing.
    if '_Spectra_Quality' not in unified_df.columns:
        if 'Tryptic_Peptides' in unified_df.columns and mirror_index:
            unified_df['_Spectra_Quality'] = unified_df['Tryptic_Peptides'].apply(
                lambda x: get_spectra_quality(x, mirror_index)
            )
        else:
            unified_df['_Spectra_Quality'] = unified_df.get(
                'Tryptic_Peptides', pd.Series(dtype=object)
            ).apply(lambda x: 'Insufficient' if (x and not pd.isna(x)) else 'No PROSIT')

    # ── Route: ?mp=<entry id> takes over the whole page with that microprotein's
    #    entry, replacing the table and the filter sidebar. Checked before any
    #    of them render, so the entry page skips the cross-filter engine and the
    #    display-table prep entirely. Resolved against the *unfiltered* frame so
    #    a shared link opens regardless of the recipient's filter state. ──
    _entry_key = st.query_params.get(ENTRY_QUERY_PARAM)
    if _entry_key:
        _render_entry_page(unified_df, _entry_key, mirror_index,
                           expression_index, seq_to_coords)
        return

    # Back on the table path: re-seed anything the entry-page detour culled.
    # Must precede the search box and the sidebar, which are the widgets it
    # restores.
    _restore_filter_state()

    # ── Search-first landing strip (a distinct "sandbox" zone to explore by example) ──
    # NOTE: filter session-state keys all use the "f_" prefix by convention.
    st.session_state.setdefault("quick_preset", None)
    st.markdown(
        "<style>div[data-testid='stVerticalBlockBorderWrapper']:has(.hero-sandbox-marker){"
        "background:linear-gradient(135deg,rgba(116,162,183,0.16),rgba(237,134,81,0.07))!important;"
        "border:1.5px dashed rgba(116,162,183,0.7)!important;"
        "border-radius:14px!important;padding:0.8rem 1.15rem 0.55rem!important;"
        "box-shadow:0 2px 16px rgba(0,0,0,0.18)!important;margin-bottom:0.8rem!important;}</style>",
        unsafe_allow_html=True,
    )
    with st.container(border=True):
        st.markdown(
            "<div class='hero-sandbox-marker' style='font-size:0.78rem; font-weight:700; "
            "letter-spacing:0.05em; text-transform:uppercase; color:#8da8b8; margin-bottom:0.25rem;'>"
            "Quick start \u00b7 example views</div>",
            unsafe_allow_html=True,
        )
        hero_search = st.text_input(
            "Search microproteins",
            key="hero_search",
            placeholder="Search a gene (e.g. BCL3) or paste an amino-acid sequence (e.g. MAASGK)\u2026",
            label_visibility="collapsed",
        )
        st.caption(
            f"Explore {len(unified_df):,} brain microproteins. These quick views are just examples to "
            "toy with \u2014 or build your own with the search above and the sidebar filters."
        )
        _quick_presets = [
            ("lncRNA Microproteins", "lncrna"),
            ("Strong Graded Microproteins", "strong_graded"),
        ]
        _active_preset = st.session_state.get("quick_preset")
        _preset_cols = st.columns(len(_quick_presets) + 1)
        for _pcol, (_plabel, _pkey) in zip(_preset_cols, _quick_presets):
            if _pcol.button(
                _plabel,
                key=f"qp_{_pkey}",
                use_container_width=True,
                type=("primary" if _active_preset == _pkey else "secondary"),
            ):
                st.session_state["quick_preset"] = None if _active_preset == _pkey else _pkey
                st.rerun()
        if _preset_cols[-1].button(
            "Show all",
            key="qp_clear",
            use_container_width=True,
            disabled=not _active_preset,
        ):
            st.session_state["quick_preset"] = None
            st.rerun()

    # ── Sidebar: ALL filters (UniProt-style checkbox facets with live counts).
    #    Filtering happens INLINE here, narrowing a running `filtered_df`
    #    sequentially as each section renders, so every facet's counts
    #    reflect everything selected in the sections above it. ──
    with st.sidebar:
        st.markdown('<div style="font-size:1.05rem; font-weight:700; color:#ffffff; margin-bottom:0.7rem; padding-bottom:0.4rem; border-bottom:1px solid rgba(116,162,183,0.25);">Filters</div>', unsafe_allow_html=True)
        st.caption("Search is at the top of the page. Use these controls to refine results.")

        base_df = unified_df.copy()

        # Unified search box (top of page): matches parent gene OR amino-acid sequence.
        if hero_search:
            _sq = hero_search.strip()
            _smask = pd.Series(False, index=base_df.index)
            if 'Parent_Gene' in base_df.columns:
                _smask = _smask | base_df['Parent_Gene'].str.contains(_sq, case=False, na=False)
            if 'sequence' in base_df.columns:
                _smask = _smask | base_df['sequence'].str.contains(_sq, case=False, na=False)
            base_df = base_df[_smask]

        # Curated quick-view preset (landing-strip buttons) — applied first so
        # it acts as a coarse pre-filter and every facet below shows counts
        # consistent with an active preset.
        _qp = st.session_state.get("quick_preset")
        if _qp == "lncrna":
            if 'smORF_Group' in base_df.columns:
                base_df = base_df[base_df['smORF_Group'] == 'lncRNA']
        elif _qp == "strong_graded":
            if '_Spectra_Quality' in base_df.columns:
                base_df = base_df[base_df['_Spectra_Quality'] == 'Strong']

        # ── Cross-filtering engine ───────────────────────────────────────────
        # Every filter is turned into a row mask over `base_df` BEFORE any widget
        # renders (selections are read straight from session_state, whose keys
        # Streamlit has already updated for this rerun). A facet is then drawn
        # against `_narrow(skip=<itself>)` — the data narrowed by all the *other*
        # filters — so checking a box updates the counts in every other facet,
        # including the ones above it, instead of only the ones below.
        _all_true = pd.Series(True, index=base_df.index)

        def _ss(key, default=None):
            return st.session_state.get(key, default)

        def _sel(prefix, domain):
            """Values in `domain` whose checkbox is currently ticked."""
            return [v for v in domain if st.session_state.get(f"{prefix}_{v}", False)]

        def _col(name):
            return base_df[name] if name in base_df.columns else None

        # Stable value domains (independent of what is currently filtered).
        _dom_status = ['Reviewed (Swiss-Prot)', 'Unreviewed (Salk/TrEMBL)']
        _dom_group = list(SMORF_PARENT_GROUPS.keys())
        _dom_dsub = list(SMORF_PARENT_GROUPS['Downstream'])
        _dom_codon = ['ATG', 'nonATG']
        _dom_kozak = ['strong', 'adequate', 'weak']
        _dom_quality = list(QUALITY_LEVELS)
        _dom_decay = list(DECAY_CLASS_LEVELS)
        _dom_annot = (sorted(unified_df['Annotation_Status'].dropna().unique())
                      if 'Annotation_Status' in unified_df.columns else [])
        _dom_shortstop = (sorted(unified_df['ShortStop_Label'].dropna().unique())
                          if 'ShortStop_Label' in unified_df.columns else [])
        _dom_tmt_sig = ['Tier 1', 'Tier 2', 'Tier 3', 'Tier 4']
        _dom_rna_sig = ['Significant', 'Exploratory']

        selected_status = _sel('f_status', _dom_status)
        selected_groups = _sel('f_smorf_group', _dom_group)
        # Sub-type only exists under Downstream; ignore stale ticks otherwise.
        selected_downstream_sub = (_sel('f_downstream_sub', _dom_dsub)
                                   if 'Downstream' in selected_groups else [])
        selected_codons = _sel('f_codons', _dom_codon)
        selected_kozak = _sel('f_kozak', _dom_kozak)
        selected_annotation = _sel('f_annotation', _dom_annot)
        selected_quality = _sel('f_quality', _dom_quality)
        selected_decay = _sel('f_decay', _dom_decay)
        selected_shortstop = _sel('f_shortstop', _dom_shortstop)
        selected_tmt_sig = _sel('f_tmt_sig', _dom_tmt_sig)
        selected_rna_sig = _sel('f_rna_sig', _dom_rna_sig)

        # Slider bounds come from the full dataset so they never move underfoot.
        def _bounds(col, cast):
            if col not in unified_df.columns:
                return None
            data = unified_df[col].dropna()
            if data.empty:
                return None
            return cast(data.min()), cast(data.max())

        _len_bounds = _bounds('Protein_Length', int)
        _spec_bounds = _bounds('Unique_Spectral_Counts', int)
        _phylo_bounds = _bounds('PhyloCSF_Score', float)
        _uni_bounds = ((1, 5) if ('UniProt_Annotation_Score' in unified_df.columns
                                  and unified_df['UniProt_Annotation_Score'].notna().any())
                       else None)

        length_range = _ss('f_length_range', _len_bounds) if _len_bounds else None
        spectral_range = _ss('f_spectral_range', _spec_bounds) if _spec_bounds else None
        phylocsf_range = _ss('f_phylocsf_range', _phylo_bounds) if _phylo_bounds else None
        uniprot_score_range = _ss('f_uniprot_range', _uni_bounds) if _uni_bounds else None

        selected_nterm = _sel('f_nterm', NTERM_MODES)

        _masks = {}

        def _add_mask(name, mask):
            if mask is not None:
                _masks[name] = mask

        # Status — one box ticked narrows; none or both means "show all".
        if len(selected_status) == 1:
            _sw = _is_swiss_mask(base_df)
            _add_mask('status', _sw if selected_status[0].startswith('Reviewed') else ~_sw)

        # smORF group + nested Downstream sub-type.
        if selected_groups and 'smORF_Class' in base_df.columns:
            _allowed = []
            for _g in selected_groups:
                _allowed.extend(SMORF_PARENT_GROUPS.get(_g, []))
            _add_mask('group', base_df['smORF_Class'].isin(_allowed))
        if selected_downstream_sub and 'smORF_Class' in base_df.columns:
            _dall = set(SMORF_PARENT_GROUPS['Downstream'])
            _add_mask('dsub', (~base_df['smORF_Class'].isin(_dall))
                              | base_df['smORF_Class'].isin(selected_downstream_sub))

        # Start codon / Kozak — rows with no value are always retained.
        if selected_codons and _col('Start_Codon') is not None:
            _add_mask('codon', base_df['Start_Codon'].isin(selected_codons)
                               | base_df['Start_Codon'].isna())
        if selected_kozak and _col('Kozak_Strength') is not None:
            _add_mask('kozak', base_df['Kozak_Strength'].isin(selected_kozak)
                               | base_df['Kozak_Strength'].isna())

        if selected_annotation and _col('Annotation_Status') is not None:
            _add_mask('annotation', base_df['Annotation_Status'].isin(selected_annotation))
        if selected_quality and _col('_Spectra_Quality') is not None:
            _add_mask('quality', base_df['_Spectra_Quality'].isin(selected_quality))
        if selected_decay and _col('NMD_Decay_Class') is not None:
            _add_mask('decay', base_df['NMD_Decay_Class'].isin(selected_decay))
        if selected_shortstop and _col('ShortStop_Label') is not None:
            _add_mask('shortstop', base_df['ShortStop_Label'].isin(selected_shortstop))

        # Significance facets. Tier 1 = strongest (q<0.05, ≥50% samples) … Tier 4
        # = weakest (q<0.2, ≥1 sample). Checked boxes are OR'd, matching every
        # other facet in the sidebar; no box ticked means no significance filter.
        # Matching on the stable tier token means a changed emoji or count in the
        # label can never silently disable a filter.
        _tmt_tiers = {
            'Tier 1': ('TMT_qvalue_50pct', 0.05),
            'Tier 2': ('TMT_qvalue_50pct', 0.2),
            'Tier 3': ('TMT_qvalue_0pct', 0.05),
            'Tier 4': ('TMT_qvalue_0pct', 0.2),
        }
        _tmt_any = None
        for _tier in selected_tmt_sig:
            _qcol, _thr = _tmt_tiers[_tier]
            if _qcol not in base_df.columns:
                continue
            _m = pd.to_numeric(base_df[_qcol], errors='coerce') < _thr
            _tmt_any = _m if _tmt_any is None else (_tmt_any | _m)
        _add_mask('tmt', _tmt_any)

        _rna_cols = {
            'Significant': 'ROSMAP_Highly_Significant',
            'Exploratory': 'ROSMAP_Significant',
        }
        _rna_any = None
        for _lvl in selected_rna_sig:
            _rcol = _rna_cols[_lvl]
            if _rcol not in base_df.columns:
                continue
            _m = base_df[_rcol].fillna(False).astype(bool)
            _rna_any = _m if _rna_any is None else (_rna_any | _m)
        _add_mask('rna', _rna_any)

        # Numeric ranges — only a mask when the slider actually narrows, and rows
        # with no value are always retained.
        def _range_mask(name, col, rng, bounds):
            if not rng or not bounds or col not in base_df.columns:
                return
            if rng[0] <= bounds[0] and rng[1] >= bounds[1]:
                return
            _add_mask(name, base_df[col].between(rng[0], rng[1]) | base_df[col].isna())

        _range_mask('length', 'Protein_Length', length_range, _len_bounds)
        _range_mask('spectral', 'Unique_Spectral_Counts', spectral_range, _spec_bounds)
        _range_mask('phylocsf', 'PhyloCSF_Score', phylocsf_range, _phylo_bounds)
        _range_mask('uniprot', 'UniProt_Annotation_Score', uniprot_score_range, _uni_bounds)

        # N-terminus: discount the Met-excision artefact of the parent ORF.
        # _add_mask ignores None, so "All microproteins" needs no branch here.
        _nt_any = None
        for _mode in selected_nterm:
            _m = _nterm_mask(base_df, _mode)
            if _m is None:
                continue
            _nt_any = _m if _nt_any is None else (_nt_any | _m)
        _add_mask('nterm', _nt_any)

        def _narrow(skip=None):
            """base_df narrowed by every active filter except the named one(s)."""
            skip = {skip} if isinstance(skip, str) else set(skip or ())
            m = _all_true
            for _name, _mk in _masks.items():
                if _name not in skip:
                    m = m & _mk
            return base_df[m]

        # The fully filtered result — every widget below is drawn against a
        # cross-filtered view, so this is computed up front rather than
        # accumulated as the sidebar renders.
        filtered_df = _narrow()

        # ── Status ──
        if 'Database_All' in base_df.columns or 'Database' in base_df.columns:
            with st.expander("Status", expanded=False):
                _v = _narrow(skip='status')
                _status_series = pd.Series(
                    np.where(_is_swiss_mask(_v), 'Reviewed (Swiss-Prot)', 'Unreviewed (Salk/TrEMBL)'),
                    index=_v.index,
                )
                _render_facet_checkboxes(
                    _status_series, key_prefix='f_status',
                    order=_dom_status,
                    help="Check one to narrow; check both (or neither) to show all.",
                )

        # ── smORF Type (+ nested Downstream Sub-type) ──
        if 'smORF_Group' in base_df.columns:
            with st.expander("smORF Type", expanded=False):
                _render_facet_checkboxes(
                    _narrow(skip={'group', 'dsub'})['smORF_Group'], key_prefix='f_smorf_group',
                    order=_dom_group,
                    help="Top-level smORF category.",
                )
                if 'Downstream' in selected_groups and 'smORF_Class' in base_df.columns:
                    _dv = _narrow(skip='dsub')
                    _dv = _dv[_dv['smORF_Class'].isin(_dom_dsub)]
                    if _dv['smORF_Class'].notna().any() or selected_downstream_sub:
                        st.markdown("**Downstream Sub-type**")
                        _render_facet_checkboxes(
                            _dv['smORF_Class'], key_prefix='f_downstream_sub',
                            order=_dom_dsub,
                            help="Refine within Downstream ORFs (doORF, daORF, etc.)",
                        )

        # ── Evidence & Quality ──
        with st.expander("Evidence & Quality", expanded=False):
            if 'Annotation_Status' in base_df.columns:
                _render_facet_checkboxes(
                    _narrow(skip='annotation')['Annotation_Status'], key_prefix='f_annotation',
                    domain=_dom_annot,
                    help="MS = detected by DDA mass spectrometry; RiboCode-ShortStop = RiboCode + "
                         "ShortStop SAM translation evidence without DDA detection (DIA proteomics "
                         "is not counted).",
                )
            st.markdown("**Spectra Quality**")
            _render_facet_checkboxes(
                _narrow(skip='quality')['_Spectra_Quality'], key_prefix='f_quality',
                order=_dom_quality,
                help="Best Prosit spectral-match confidence tier (from the master Confidence column).",
            )

        # ── Differential Expression ──
        # Checkbox facets like the rest of the sidebar: the checkbox *key* carries
        # the stable tier token and the live count lives only in the label.
        with st.expander("Differential Expression", expanded=False):
            _sig_checkbox = _count_checkbox

            has_tmt_50 = 'TMT_qvalue_50pct' in base_df.columns
            has_tmt_0  = 'TMT_qvalue_0pct'  in base_df.columns
            if has_tmt_50 or has_tmt_0:
                _tv = _narrow(skip='tmt')
                _tmt_help = (
                    "Filter by TMT proteomics q-value at the chosen stringency; "
                    "ticking several tiers keeps microproteins in any of them. "
                    "≥50% samples (stringent): peptide quantified in at least half the donors per condition. "
                    "≥1 sample/condition: peptide quantified in at least one donor per condition (lenient)."
                )
                st.markdown("**TMT-MS Significance**")
                if has_tmt_50:
                    _q50 = pd.to_numeric(_tv['TMT_qvalue_50pct'], errors='coerce')
                    _sig_checkbox("f_tmt_sig_Tier 1", "Tier 1 — FDR < 0.05, ≥50% samples",
                                  int((_q50 < 0.05).sum()), help=_tmt_help)
                    _sig_checkbox("f_tmt_sig_Tier 2", "Tier 2 — FDR < 0.2, ≥50% samples",
                                  int((_q50 < 0.2).sum()))
                if has_tmt_0:
                    _q0 = pd.to_numeric(_tv['TMT_qvalue_0pct'], errors='coerce')
                    _sig_checkbox("f_tmt_sig_Tier 3", "Tier 3 — FDR < 0.05, ≥1 sample/condition",
                                  int((_q0 < 0.05).sum()), help=(_tmt_help if not has_tmt_50 else None))
                    _sig_checkbox("f_tmt_sig_Tier 4", "Tier 4 — FDR < 0.2, ≥1 sample/condition",
                                  int((_q0 < 0.2).sum()))

            if 'ROSMAP_Significant' in base_df.columns:
                _rv = _narrow(skip='rna')
                _rna_exp_n = int(_rv['ROSMAP_Significant'].sum())
                _rna_sig_n = int(_rv['ROSMAP_Highly_Significant'].sum()) if 'ROSMAP_Highly_Significant' in _rv.columns else 0
                st.markdown("**RNA Significance**")
                _sig_checkbox("f_rna_sig_Significant", "Significant (FDR < 0.05)", _rna_sig_n,
                              help="Filter by ROSMAP bulk RNA-seq adjusted p-value")
                _sig_checkbox("f_rna_sig_Exploratory", "Exploratory (FDR < 0.2)", _rna_exp_n)

        # ── ShortStop Label ──
        if 'ShortStop_Label' in base_df.columns:
            with st.expander("ShortStop Label", expanded=False):
                _render_facet_checkboxes(
                    _narrow(skip='shortstop')['ShortStop_Label'], key_prefix='f_shortstop',
                    domain=_dom_shortstop,
                    help="ShortStop ribosome profiling classification",
                )

        # ── Score & Length Ranges ──
        # Each "in range" caption counts against the view narrowed by every other
        # filter but not by this slider, matching the facet counts above.
        with st.expander("Score & Length Ranges", expanded=False):
            def _range_slider(label, name, col, bounds, step=None, help=None):
                if not bounds or col not in base_df.columns:
                    return
                st.slider(label, min_value=bounds[0], max_value=bounds[1],
                          value=bounds, step=step, key=f"f_{name}_range", help=help)
                _rv = _narrow(skip=name)
                _rng = _ss(f"f_{name}_range", bounds)
                _n = int((_rv[col].between(_rng[0], _rng[1]) | _rv[col].isna()).sum())
                st.caption(f"→ {_n:,} microproteins in range")

            _range_slider("Protein Length (aa)", "length", 'Protein_Length', _len_bounds)
            _range_slider("Unique Spectral Counts", "spectral", 'Unique_Spectral_Counts', _spec_bounds)
            _range_slider("PhyloCSF Score", "phylocsf", 'PhyloCSF_Score', _phylo_bounds, step=0.5)
            _range_slider(
                "UniProt Annotation Score (★)", "uniprot", 'UniProt_Annotation_Score',
                _uni_bounds, step=1,
                help="UniProt curation annotation score (1–5). Reviewed / TrEMBL "
                     "entries only; microproteins without a UniProt score are always retained.",
            )

        # ── Start Codon ──
        if 'Start_Codon' in base_df.columns:
            with st.expander("Start Codon", expanded=False):
                _render_facet_checkboxes(
                    _narrow(skip='codon')['Start_Codon'], key_prefix='f_codons',
                    order=_dom_codon,
                    help="ATG (canonical) vs nonATG (annotated with a near-cognate/non-standard "
                         "start — see the Alt-Initiation column view for detail).",
                )

        # ── ORF Rules (NMD × main-ORF disruption) ──
        if 'NMD_Decay_Class' in base_df.columns:
            with st.expander("ORF Rules", expanded=False):
                _render_facet_checkboxes(
                    _narrow(skip='decay')['NMD_Decay_Class'], key_prefix='f_decay',
                    order=_dom_decay,
                    label_fn=lambda v: DECAY_CLASS_EMOJI.get(v, v),
                    help=DECAY_CLASS_HELP,
                )

        # ── Kozak Context ──
        if 'Kozak_Strength' in base_df.columns and base_df['Kozak_Strength'].notna().any():
            with st.expander("Kozak Context", expanded=False):
                _render_facet_checkboxes(
                    _narrow(skip='kozak')['Kozak_Strength'], key_prefix='f_kozak',
                    order=_dom_kozak,
                    help="Sequence-context strength around the start codon (Salk-discovered smORFs only).",
                    label_fn=lambda v: v.capitalize(),
                )

        # ── N-terminus Options ──
        with st.expander("N-terminus Options", expanded=False):
            _nv = _narrow(skip='nterm')

            def _nterm_count(mode):
                _m = _nterm_mask(_nv, mode)
                return 0 if _m is None else int(_m.sum())

            _nterm_help = (
                "A tryptic peptide starting at position 1 or 2 is exactly the peptide "
                "N-terminal Met excision of the parent ORF would produce, so on its own it "
                "is not evidence the microprotein is real.\n\n"
                "Non-N-terminal peptides only: keep microproteins with a tryptic peptide "
                "starting at position 3 or later — evidence that survives removing the "
                "M-excision peptides.\n\n"
                "Nt-acetylated, substitution-distinct, or non-N-terminal: as above, but also "
                "keep microproteins with an N-terminally acetylated peptide (co-translational, "
                "marks a genuine N-terminus) OR a position-1/2 peptide whose own sequence carries "
                "1-2 amino-acid substitutions vs. its matched UniProt isoform — a peptide that "
                "differs from the canonical protein's sequence at that position can't be a bare "
                "Met-excision fragment of it, so that's real evidence too.\n\n"
                "Ticking both is the same as ticking the second — the first is a subset. "
                "Both require tryptic evidence, so microproteins with no peptides are excluded."
            )
            for _i, _mode in enumerate(NTERM_MODES):
                _count_checkbox(f"f_nterm_{_mode}", _mode, _nterm_count(_mode),
                                help=(_nterm_help if _i == 0 else None))

        # ── Group by ──
        _GROUP_BY_COLUMN_MAP = {
            "smORF Type": "smORF_Group",
            "Annotation Method": "Annotation_Status",
            "Database Status": "_status_label",
            "Start Codon": "Start_Codon",
        }
        with st.expander("Group by", expanded=False):
            group_by_choice = st.selectbox(
                "Group by", options=["None"] + list(_GROUP_BY_COLUMN_MAP.keys()),
                key="f_group_by",
                help="Reorder the table so rows sharing this value are grouped together; "
                     "shows a per-group count breakdown above the table.",
            )
        group_by_col = None
        if group_by_choice != "None":
            group_by_col = _GROUP_BY_COLUMN_MAP[group_by_choice]
            if group_by_col == "_status_label":
                filtered_df = filtered_df.copy()
                filtered_df['_status_label'] = np.where(
                    _is_swiss_mask(filtered_df), 'Reviewed (Swiss-Prot)', 'Unreviewed (Salk/TrEMBL)'
                )
            if group_by_col not in filtered_df.columns:
                group_by_col = None

        st.markdown('<div class="sidebar-section-header">Legend</div>', unsafe_allow_html=True)
        st.markdown("""
        <div class="legend-container">
            <div style="font-weight:600; margin-bottom:0.5rem; font-size:0.85rem;">Legend</div>
            <div style="margin-bottom:0.35rem;"><span class="badge-swiss">Reviewed</span> Swiss-Prot</div>
            <div><span class="badge-unreviewed">Unreviewed</span> Salk / TrEMBL</div>
        </div>
        """, unsafe_allow_html=True)
        st.caption("Click links to view in UCSC Browser")
        st.markdown(
            "<div style='margin-top:0.6rem; font-size:0.8rem; color:#8da8b8;'>"
            "<a href='https://github.com/brendan-miller-salk/brain-microprotein-atlas' "
            "target='_blank' style='color:#74c2e1; text-decoration:none;'>"
            "GitHub Repository</a></div>",
            unsafe_allow_html=True,
        )

    # Every filter widget has now rendered — snapshot them so an entry-page
    # detour on a later rerun can be undone.
    _snapshot_filter_state()

    # ── Command-center header with live metrics ──
    total = len(filtered_df)
    genes = filtered_df['Parent_Gene'].nunique() if 'Parent_Gene' in filtered_df.columns else 0
    tmt_sig = int(filtered_df['TMT_Significant'].sum()) if 'TMT_Significant' in filtered_df.columns else 0
    rna_sig = int(filtered_df['ROSMAP_Significant'].sum()) if 'ROSMAP_Significant' in filtered_df.columns else 0

    swiss_count = int(_is_swiss_mask(filtered_df).sum())
    noncan_count = total - swiss_count

    header_slot.markdown(f"""
    <div class="cmd-header">
        <div class="cmd-top">
            <div class="cmd-title">Brain Microprotein Dashboard</div>
            <div class="cmd-subtitle">
                Saghatelian Lab &middot; Salk Institute<br>
                TMT-MS &middot; Ribo-seq &middot; Short/Long-Read RNA-seq
            </div>
        </div>
        <div class="cmd-metrics">
            <div class="cmd-stat st-total">
                <div class="stat-label">Total</div>
                <div class="stat-val">{total:,}</div>
                <div class="stat-sub">{genes:,} genes</div>
            </div>
            <div class="cmd-stat st-swiss">
                <div class="stat-label">Reviewed</div>
                <div class="stat-val">{swiss_count:,}</div>
                <div class="stat-sub">{swiss_count/total*100 if total else 0:.1f}% Swiss-Prot</div>
            </div>
            <div class="cmd-stat st-noncan">
                <div class="stat-label">Unreviewed</div>
                <div class="stat-val">{noncan_count:,}</div>
                <div class="stat-sub">{noncan_count/total*100 if total else 0:.1f}% novel</div>
            </div>
            <div class="cmd-stat">
                <div class="stat-label">smORF Groups</div>
                <div class="stat-val">{filtered_df['smORF_Group'].nunique() if 'smORF_Group' in filtered_df.columns else 0}</div>
                <div class="stat-sub">{filtered_df['smORF_Class'].nunique() if 'smORF_Class' in filtered_df.columns else 0} types</div>
            </div>
            <div class="cmd-stat">
                <div class="stat-label">TMT Significant</div>
                <div class="stat-val">{tmt_sig:,}</div>
                <div class="stat-sub">q &lt; 0.2</div>
            </div>
            <div class="cmd-stat">
                <div class="stat-label">RNA Significant</div>
                <div class="stat-val">{rna_sig:,}</div>
                <div class="stat-sub">padj &lt; 0.2</div>
            </div>
        </div>
    </div>
    """, unsafe_allow_html=True)

    # ── Active-filter summary: show exactly which sidebar controls / presets are
    #    narrowing the current view, so users always know what is going on. ──
    def _range_narrowed(col, rng):
        if not rng or col not in unified_df.columns:
            return False
        _d = unified_df[col].dropna()
        return (not _d.empty) and (rng[0] > _d.min() or rng[1] < _d.max())

    _preset_labels = {
        "lncrna": "lncRNA Microproteins",
        "strong_graded": "Strong Graded Microproteins",
    }
    _active_chips = []
    if hero_search and hero_search.strip():
        _active_chips.append(("Search", hero_search.strip()))
    if _qp in _preset_labels:
        _active_chips.append(("Quick view", _preset_labels[_qp]))
    if selected_groups:
        _active_chips.append(("smORF type", ", ".join(map(str, selected_groups))))
    if selected_downstream_sub:
        _active_chips.append(("Downstream sub-type", ", ".join(map(str, selected_downstream_sub))))
    if selected_decay:
        _active_chips.append(("ORF Rules", ", ".join(map(str, selected_decay))))
    if len(selected_status) == 1:
        _active_chips.append(("Status", selected_status[0]))
    if selected_annotation:
        _active_chips.append(("Evidence type", ", ".join(map(str, selected_annotation))))
    if selected_quality:
        _active_chips.append(("Spectra quality", ", ".join(map(str, selected_quality))))
    if selected_tmt_sig:
        _active_chips.append(("TMT-MS", ", ".join(map(str, selected_tmt_sig))))
    if selected_rna_sig:
        _active_chips.append(("RNA", ", ".join(map(str, selected_rna_sig))))
    if selected_codons:
        _active_chips.append(("Start codon", ", ".join(map(str, selected_codons))))
    if selected_kozak:
        _active_chips.append(("Kozak strength", ", ".join(map(str, selected_kozak))))
    if selected_shortstop:
        _active_chips.append(("ShortStop", ", ".join(map(str, selected_shortstop))))
    if selected_nterm:
        _active_chips.append(("N-terminus", ", ".join(map(str, selected_nterm))))
    if _range_narrowed('Protein_Length', length_range):
        _active_chips.append(("Protein length", f"{int(length_range[0])}\u2013{int(length_range[1])} aa"))
    if _range_narrowed('Unique_Spectral_Counts', spectral_range):
        _active_chips.append(("Unique spectra", f"{int(spectral_range[0])}\u2013{int(spectral_range[1])}"))
    if _range_narrowed('PhyloCSF_Score', phylocsf_range):
        _active_chips.append(("PhyloCSF", f"{phylocsf_range[0]:.1f}\u2013{phylocsf_range[1]:.1f}"))
    if _range_narrowed('UniProt_Annotation_Score', uniprot_score_range):
        _active_chips.append(("UniProt score", f"{int(uniprot_score_range[0])}\u2013{int(uniprot_score_range[1])}\u2605"))
    if group_by_col:
        _active_chips.append(("Grouped by", group_by_choice))

    _render_active_filter_bar(_active_chips, len(filtered_df), len(unified_df))

    if len(filtered_df) == 0:
        st.warning("No microproteins match your current filters. Try broadening the criteria.")
        return

    # ── Default sort: SAMs → Strong/Moderate spectra → lncRNA/uORF/dORF → high spec counts ──
    _sort_df = filtered_df.copy()
    # 1. ShortStop SAMs first (SAM-Intracellular, SAM-Secreted)
    if 'ShortStop_Label' in _sort_df.columns:
        _sort_df['_s1_sam'] = (~_sort_df['ShortStop_Label'].fillna('').str.startswith('SAM')).astype(int)
    else:
        _sort_df['_s1_sam'] = 1
    # 2. Strong or Moderate spectra quality first
    if '_Spectra_Quality' in _sort_df.columns:
        _sort_df['_s2_quality'] = _sort_df['_Spectra_Quality'].map(
            lambda q: 0 if q in ('Strong', 'Moderate') else 1
        )
    else:
        _sort_df['_s2_quality'] = 1
    # 3. lncRNA, uORF (Upstream), dORF (Downstream) smORF types first
    _lnc_u_d_types = set(
        ['lncRNA']
        + SMORF_PARENT_GROUPS.get('Upstream', [])
        + SMORF_PARENT_GROUPS.get('Downstream', [])
    )
    if 'smORF_Class' in _sort_df.columns:
        _sort_df['_s3_type'] = (~_sort_df['smORF_Class'].isin(_lnc_u_d_types)).astype(int)
    else:
        _sort_df['_s3_type'] = 1
    # 4. Highest spectral counts first
    if 'Unique_Spectral_Counts' in _sort_df.columns:
        _sort_df['_s4_spec'] = -_sort_df['Unique_Spectral_Counts'].fillna(0)
    else:
        _sort_df['_s4_spec'] = 0
    _sort_cols = ['_s1_sam', '_s2_quality', '_s3_type', '_s4_spec']
    if group_by_col and group_by_col in _sort_df.columns:
        # Group-by is the primary sort key; existing tiers act as tiebreak
        # within each group, so grouped rows still favor stronger evidence first.
        _sort_cols = [group_by_col] + _sort_cols
    filtered_df = _sort_df.sort_values(_sort_cols).drop(columns=[c for c in _sort_cols if c != group_by_col])

    # ── Prepare display table ──
    display_cols = ['sequence']


    # NonATG Alternative Sequence(s): raw nonATG_predicted_sequence (a stringified
    # Python list, e.g. "['MLWGR...', 'LLWGR...']") shown as-is, brackets included —
    # blank for ATG-start / no-alternative-site microproteins.
    if 'NonATG_Predicted_Sequence' in filtered_df.columns:
        filtered_df['NonATG Alternative Sequence(s)'] = filtered_df['NonATG_Predicted_Sequence']
        display_cols.append('NonATG Alternative Sequence(s)')

    # Non-AUG alternative-initiation candidates, joined "codon@position" summary
    # — Alt-Initiation tab only (see _render_results_table's column ordering).
    if '_NonATG_Candidates_Display' in filtered_df.columns:
        filtered_df['NonATG Candidates'] = filtered_df['_NonATG_Candidates_Display']
        display_cols.append('NonATG Candidates')

    if '_UCSC_HTML' in filtered_df.columns:
        filtered_df = filtered_df.copy()
        filtered_df['UCSC'] = filtered_df['_UCSC_HTML']
        display_cols.append('UCSC')
    elif 'UCSC_Link' in filtered_df.columns:
        filtered_df = filtered_df.copy()
        filtered_df['UCSC'] = filtered_df.apply(
            lambda row: create_ucsc_link({'CLICK_UCSC': row['UCSC_Link']}, CUSTOM_UCSC_SESSION)
            if pd.notna(row.get('UCSC_Link')) else None, axis=1
        )
        display_cols.append('UCSC')

    column_mapping = {
        'Parent_Gene': 'Parent Gene',
        'smORF_Group': 'General smORF Type',
        'smORF_Display': 'smORF Subtype',
        'Protein_Length': 'Protein Length',
        'PhyloCSF_Score': 'PhyloCSF Score',
        'ShortStop_Label': 'ShortStop',
        'ShortStop_Score': 'ShortStop Score',
        'Annotation_Status': 'Annotation Method',
        'Unique_Spectral_Counts': 'Unique Spectral Counts (DDA)',
        'Razor_Spectral_Counts': 'Razor Counts (DDA)',
        'Start_Codon': 'Start Codon',
        'TMT_log2fc_50pct': 'TMT log2FC (50%)',
        'TMT_t_statistic_50pct': 'TMT t-stat (50%)',
        'TMT_df_50pct': 'TMT df (50%)',
        'TMT_conf_low_50pct': 'TMT CI low (50%)',
        'TMT_conf_high_50pct': 'TMT CI high (50%)',
        'TMT_cohens_d_50pct': "TMT Cohen's d (50%)",
        'TMT_pvalue_50pct': 'TMT p-val (50%)',
        'TMT_qvalue_50pct': 'TMT q-val (50%)',
        'TMT_log2fc_0pct':  'TMT log2FC (0%)',
        'TMT_t_statistic_0pct': 'TMT t-stat (0%)',
        'TMT_df_0pct': 'TMT df (0%)',
        'TMT_conf_low_0pct': 'TMT CI low (0%)',
        'TMT_conf_high_0pct': 'TMT CI high (0%)',
        'TMT_cohens_d_0pct': "TMT Cohen's d (0%)",
        'TMT_pvalue_0pct':  'TMT p-val (0%)',
        'TMT_qvalue_0pct':  'TMT q-val (0%)',
        'TMT_rate_control': 'MS Detect Control',
        'TMT_rate_ad':      'MS Detect AD',
        'ROSMAP_log2FC': 'ROSMAP log2FC',
        'ROSMAP_padj': 'ROSMAP padj',
        'MSBB_log2FC': 'MSBB log2FC',
        'MSBB_padj': 'MSBB padj',
        'Corr_MainORF_NonAD': 'ROSMAP Corr NonAD',
        'Corr_MainORF_AD': 'ROSMAP Corr AD',
        'Corr_MainORF_NonAD_MSBB': 'MSBB Corr NonAD',
        'Corr_MainORF_AD_MSBB': 'MSBB Corr AD',
        'RNA_LRT_Add_P': 'RNA_LRT_Add_P',
        'RNA_LRT_Int_P': 'RNA_LRT_Int_P',
        'RP3_Default': 'RP3 Default',
        'RP3_MM_Amb': 'RP3 MM+Amb',
        'RP3_Amb': 'RP3 Amb',
        'RP3_MM': 'RP3 MM',
        'RiboCode': 'RiboCode',
        'Tryptic_Peptides': 'Tryptic Peptides',
        'Tryptic_Protein_ID': 'Tryptic Protein ID',
        'Tryptic_Start_Positions': 'Tryptic Start Positions',
        'Tryptic_End_Positions': 'Tryptic End Positions',
        'Nt_Acetylated': 'Nt-Acetylated',
        'Nt_Acetyl_Peptides': 'Nt-Acetyl Peptides',
        'Nt_Acetyl_Total_PSMs': 'Nt-Acetyl PSMs',
        'Nt_Acetyl_Max_Fraction': 'Nt-Acetyl PSM Fraction',
        'BLAST_UniProt_Match': 'BLAST UniProt Match',
        'BLAST_Pct_Match': 'BLAST % Match',
        'BLAST_Aln_Length': 'BLAST Aln Length',
        'BLAST_Evalue': 'BLAST E-value',
        'BLAST_Bit': 'BLAST Bit Score',
        'UniProt_Annotation_Score': 'UniProt Annotation Score',
        'Kozak_Strength': 'Kozak Strength',
        'Kozak_Downstream_Strength': 'Kozak Downstream Strength',
        'Kozak_Class_Canonical': 'Kozak Class',
        'Kozak_Window': 'Kozak Window',
        'NonATG_Annotated_Codon': 'NonATG Codon',
        'NonATG_Is_Initiator': 'NonATG Valid Initiator',
        'NonATG_Has_Supported_Site': 'NonATG Has Supported Site',
        'NonATG_Context_Strength': 'NonATG Context Strength',
        'NonATG_N_Sites': 'NonATG Alt Sites Found',
        'NonATG_Best_Tier': 'Optimal Codon Tier',
    }

    for orig, display_name in column_mapping.items():
        if orig in filtered_df.columns:
            filtered_df[display_name] = filtered_df[orig]
            display_cols.append(display_name)

    # UniProt annotation score (gold stars; reviewed/TrEMBL only, blank otherwise)
    if '_UniProt_Annotation_Display' in filtered_df.columns:
        filtered_df['UniProt Annotation'] = filtered_df['_UniProt_Annotation_Display']
        display_cols.append('UniProt Annotation')

    # Mirror plot quality (use pre-computed column)
    if '_Spectra_Quality' in filtered_df.columns:
        filtered_df['Spectra Quality'] = filtered_df['_Spectra_Quality'].map(
            lambda q: QUALITY_EMOJI.get(q, '\u2014')
        )
        display_cols.append('Spectra Quality')

    # TMT significance tier (Tier 1 = strongest / Dark red … Tier 4 = Green)
    if '_TMT_Tier' in filtered_df.columns:
        filtered_df['TMT Tier'] = filtered_df['_TMT_Tier'].map(
            lambda t: TMT_TIER_EMOJI.get(t, '\u2014')
        )
        display_cols.append('TMT Tier')

    # RNA significance (Yellow = Exploratory, Green = Significant)
    if '_RNA_Sig' in filtered_df.columns:
        filtered_df['RNA Significance'] = filtered_df['_RNA_Sig'].map(
            lambda t: RNA_SIG_EMOJI.get(t, '\u2014')
        )
        display_cols.append('RNA Significance')

    # smORF decay class (NMD × main-ORF disruption); dots match the BED track
    # colors. Raw NMD_Decay_Class stays untouched for the sidebar facet.
    if 'NMD_Decay_Class' in filtered_df.columns:
        filtered_df['ORF Rules'] = filtered_df['NMD_Decay_Class'].map(
            lambda v: DECAY_CLASS_EMOJI.get(v, DECAY_CLASS_EMOJI['NA'])
        )
        display_cols.append('ORF Rules')

    display_df = filtered_df[display_cols].copy()

    # ── Format tryptic peptides as readable string (precomputed if available) ──
    if 'Tryptic Peptides' in display_df.columns:
        if '_Tryptic_Display' in filtered_df.columns:
            display_df['Tryptic Peptides'] = filtered_df['_Tryptic_Display'].values
        else:
            def _fmt_peptides(val):
                if not _not_na(val):
                    return ''
                try:
                    peps = ast.literal_eval(str(val))
                    if isinstance(peps, str):
                        peps = [peps]
                    return ' \u00b7 '.join(str(p).strip() for p in peps)
                except (ValueError, SyntaxError):
                    return str(val)
            display_df['Tryptic Peptides'] = display_df['Tryptic Peptides'].map(_fmt_peptides)

    # ── Mark Swiss-Prot rows: type columns get "Swiss-Prot" label, other
    #    annotation-only fields are non-applicable (em-dash). RP3 metrics
    #    DO apply to Swiss-Prot microproteins and should display normally. ──
    _swiss_mask = None
    if 'Database_All' in filtered_df.columns or 'Database' in filtered_df.columns:
        _swiss_mask = _is_swiss_mask(filtered_df).values
        # These columns are *not applicable* to Swiss-Prot reviewed entries.
        _na_cols = ['ShortStop', 'Annotation Method']
        for _c in _na_cols:
            if _c in display_df.columns:
                display_df[_c] = display_df[_c].fillna('').astype(str)
                display_df.loc[_swiss_mask, _c] = '\u2014'
        # smORF type columns get the "Swiss-Prot" badge for reviewed entries.
        for _c in ('General smORF Type', 'smORF Subtype'):
            if _c in display_df.columns:
                display_df[_c] = display_df[_c].fillna('').astype(str)
                display_df.loc[_swiss_mask, _c] = 'Swiss-Prot'
        if 'ShortStop Score' in display_df.columns:
            _ss = pd.to_numeric(display_df['ShortStop Score'], errors='coerce')
            display_df['ShortStop Score'] = _ss.map(lambda v: f'{v:.3f}' if pd.notna(v) else '')
            display_df.loc[_swiss_mask, 'ShortStop Score'] = '\u2014'

    # ── Per-row background tint: muted teal (Swiss-Prot) vs muted orange (Unreviewed) ──
    # NOTE: pandas Styler is too slow on ~6.5k rows; we use an emoji column
    # ("DB") below for at-a-glance differentiation instead.
    if _swiss_mask is not None:
        display_df.insert(
            0,
            'DB',
            np.where(_swiss_mask, '🔵 Swiss-Prot', '🟠 Unreviewed'),
        )

    # ── Group-by count breakdown (shown above the table when grouping is active) ──
    if group_by_col and group_by_col in filtered_df.columns:
        _group_counts = filtered_df[group_by_col].value_counts()
        st.caption(
            f"Grouped by {group_by_choice}: "
            + " · ".join(f"{k}: {v:,}" for k, v in _group_counts.items())
        )

    # ── Data table (row click routes to the entry page) ──
    _render_results_table(filtered_df, display_df)

    # ── Downloads section ──
    _render_downloads_section(filtered_df, display_df)


@st.fragment
def _render_results_table(filtered_df, display_df):
    """The results table. Clicking a row navigates to that microprotein's entry
    page (see _render_entry_page). Wrapped in @st.fragment so switching the
    column view or opening the legend reruns only this block."""
    # ── Column-group presets: limit horizontal scrolling by showing a curated
    #    subset of columns. Identity columns are always shown (and pinned left).
    #    Sequence columns are pushed to the far right everywhere EXCEPT
    #    Alt-Initiation, where they're the point and stay up front. NonATG
    #    Candidates is Alt-Initiation-only — hidden from every other view. ──
    _sequence_cluster = ['sequence', 'NonATG Alternative Sequence(s)']
    _alt_init_only = ['NonATG Candidates']
    _core_identity = ['DB', 'UCSC', 'Parent Gene', 'General smORF Type', 'smORF Subtype',
                      'ORF Rules']
    # Alt-Initiation-only ordering: "NonATG Has Supported Site" sits as the 3rd
    # column (right after "sequence"), matching its position in the detail
    # card; "Optimal Codon Tier" follows, ahead of the alternative-sequence/
    # candidates columns.
    _identity_cols = (
        ['DB', 'sequence', 'NonATG Has Supported Site', 'Optimal Codon Tier',
         'NonATG Alternative Sequence(s)']
        + _alt_init_only
        + [c for c in _core_identity if c != 'DB']
    )
    # Overview column order, authored end-to-end (identity columns included, and
    # the sequence columns sitting mid-table rather than pushed right). Anything
    # here that isn't in display_df is dropped, so this is safe to over-specify.
    _OVERVIEW_ORDER = [
        'DB', 'UCSC', 'Parent Gene', 'General smORF Type', 'smORF Subtype',
        'Protein Length', 'Spectra Quality', 'PhyloCSF Score',
        'ShortStop', 'ShortStop Score', 'sequence', 'TMT Tier',
        'RNA Significance', 'Start Codon', 'ORF Rules', 'Kozak Strength',
        'NonATG Alternative Sequence(s)', 'UniProt Annotation',
    ]
    _column_groups = {
        # Overview is the one view with a hand-authored full order (see
        # _OVERVIEW_ORDER below): it interleaves the identity and sequence
        # columns rather than taking the generic identity-first /
        # sequence-last treatment, so it bypasses that assembly entirely.
        'Overview': [
            'General smORF Type', 'smORF Subtype', 'Protein Length',
            'Spectra Quality', 'PhyloCSF Score', 'ShortStop', 'ShortStop Score',
            'sequence', 'TMT Tier', 'RNA Significance', 'Start Codon',
            'ORF Rules', 'Kozak Strength', 'NonATG Alternative Sequence(s)',
            'UniProt Annotation',
        ],
        'Proteomics': [
            'Unique Spectral Counts (DDA)', 'Razor Counts (DDA)',
            'TMT log2FC (50%)', 'TMT t-stat (50%)', 'TMT df (50%)',
            'TMT CI low (50%)', 'TMT CI high (50%)', "TMT Cohen's d (50%)",
            'TMT p-val (50%)', 'TMT q-val (50%)',
            'TMT log2FC (0%)', 'TMT t-stat (0%)', 'TMT df (0%)',
            'TMT CI low (0%)', 'TMT CI high (0%)', "TMT Cohen's d (0%)",
            'TMT p-val (0%)', 'TMT q-val (0%)',
            'MS Detect Control', 'MS Detect AD', 'TMT Tier', 'Spectra Quality',
            'Tryptic Peptides', 'Tryptic Protein ID',
            'Tryptic Start Positions', 'Tryptic End Positions',
            'Nt-Acetylated', 'Nt-Acetyl Peptides', 'Nt-Acetyl PSMs',
            'Nt-Acetyl PSM Fraction',
        ],
        'Transcriptomics': [
            'ROSMAP log2FC', 'MSBB log2FC',
            'ROSMAP padj', 'MSBB padj',
            'ROSMAP Corr NonAD', 'MSBB Corr NonAD',
            'ROSMAP Corr AD', 'MSBB Corr AD',
            'RNA_LRT_Add_P', 'RNA_LRT_Int_P', 'RNA Significance',
        ],
        'Ribo-seq': [
            'RP3 Default', 'RP3 MM+Amb', 'RP3 Amb', 'RP3 MM', 'RiboCode',
        ],
        'Homology (BLAST)': [
            'BLAST UniProt Match', 'BLAST % Match', 'BLAST Aln Length',
            'BLAST E-value', 'BLAST Bit Score',
        ],
        'Alt-Initiation': [
            'General smORF Type', 'smORF Subtype', 'ORF Rules', 'Start Codon',
            'Kozak Strength', 'Kozak Downstream Strength', 'Kozak Class', 'Kozak Window',
            'NonATG Codon', 'NonATG Valid Initiator', 'NonATG Has Supported Site',
            'NonATG Context Strength', 'NonATG Alt Sites Found',
        ],
        'All columns': None,
    }
    _preset = st.segmented_control(
        'Column view',
        options=list(_column_groups.keys()),
        default='Overview',
        help='Switch between curated column sets to avoid scrolling. '
             'Identity columns (DB, UCSC, Parent Gene, smORF type/subtype) are shown in every '
             'view; sequence columns (Sequence, NonATG Alternative Sequence(s)) are pushed to '
             'the far right, except in Alt-Initiation where they stay up front. NonATG '
             'Candidates only appears in the Alt-Initiation view.',
    ) or 'Overview'

    _group = _column_groups.get(_preset)
    if _preset == 'Overview':
        _column_order = [c for c in _OVERVIEW_ORDER if c in display_df.columns]
    elif _preset == 'Alt-Initiation':
        # Alt-Initiation is only meaningful for nonATG-start smORFs — ATG-start
        # rows have no alternative-initiation candidates to show.
        if 'Start Codon' in display_df.columns:
            display_df = display_df[display_df['Start Codon'] == 'nonATG']
        # Sequence info (and NonATG Candidates) stays up front — it's the whole
        # point of this view.
        _ordered = _identity_cols + [c for c in _group if c not in _identity_cols]
        _column_order = [c for c in _ordered if c in display_df.columns]
    elif _group is None:
        # 'All columns' — core identity first, everything else in its natural
        # order, sequence cluster pushed to the very end. NonATG Candidates is
        # excluded entirely outside Alt-Initiation.
        _rest = [c for c in display_df.columns
                 if c not in _core_identity and c not in _sequence_cluster
                 and c not in _alt_init_only]
        _ordered = _core_identity + _rest + _sequence_cluster
        _column_order = [c for c in _ordered if c in display_df.columns]
    else:
        _ordered = (_core_identity + [c for c in _group if c not in _core_identity]
                    + _sequence_cluster)
        _column_order = [c for c in _ordered if c in display_df.columns]

    # ── Column legend: describes whatever columns are visible in the *current*
    #    tab, so it updates each time a different Column view is selected. ──
    with st.expander("Explain the column/variables to me", expanded=False):
        _legend_cols = _column_order if _column_order is not None else list(display_df.columns)
        st.markdown("\n".join(
            f"- **{c}** — {COLUMN_DESCRIPTIONS[c]}"
            for c in _legend_cols if c in COLUMN_DESCRIPTIONS
        ) or "No column descriptions available for this view.")

    # The row-select box in the leftmost column is the only in-place control the
    # grid offers — data-cell clicks do nothing, and st.column_config.LinkColumn
    # always opens a new tab (it calls window.open(url, "_blank") internally), so
    # a link column can't be used for same-tab navigation. Hence the signpost.
    st.markdown(
        '<div class="open-hint">'
        '<span class="open-hint-arrow">&#11013;</span>'
        '<span>Tick the <strong>select box</strong> at the start of any row to open that '
        'microprotein&rsquo;s full entry page &mdash; annotation, evidence by platform, '
        'figures and sequence.</span>'
        '</div>',
        unsafe_allow_html=True,
    )

    selection = st.dataframe(
        display_df,
        width='stretch',
        hide_index=True,
        on_select="rerun",
        selection_mode="single-row",
        # Nonce-suffixed key: returning from an entry page bumps the nonce so the
        # widget is recreated without its stale row selection (see
        # _back_to_results). Without this the table would immediately navigate
        # back to the entry we just left.
        key=f"results_table_{st.session_state.get('_table_nonce', 0)}",
        column_order=_column_order,
        column_config={
            'DB': st.column_config.TextColumn(
                'DB', help='🔵 Swiss-Prot (reviewed)  ·  🟠 Unreviewed', width='small'
            ),
            'UCSC': st.column_config.LinkColumn(
                'UCSC', help='View in UCSC Genome Browser', display_text='View'
            ),
            'sequence': st.column_config.TextColumn('Sequence', max_chars=50),
            'NonATG Alternative Sequence(s)': st.column_config.TextColumn(
                'NonATG Alternative Sequence(s)', max_chars=50,
                help='Predicted sequence(s) from alternative (non-AUG) initiation candidates, '
                     'raw list format (e.g. [\'MLWGR...\', \'LLWGR...\']); blank for ATG-start '
                     'microproteins or non-ATG starts with no compatible alternative site.'
            ),
            'Parent Gene': st.column_config.TextColumn('Parent Gene', max_chars=15),
            'General smORF Type': st.column_config.TextColumn('General smORF Type', max_chars=20),
            'smORF Subtype': st.column_config.TextColumn('smORF Subtype', max_chars=15),
            'ORF Rules': st.column_config.TextColumn(
                'ORF Rules', max_chars=27, help=DECAY_CLASS_HELP),
            'ShortStop Score': st.column_config.TextColumn('ShortStop Score', help='ShortStop ML confidence score (0-1)', max_chars=10),
            'PhyloCSF Score': st.column_config.NumberColumn('PhyloCSF Score', help='PhyloCSF evolutionary conservation score', format='%.2f'),
            'Unique Spectral Counts (DDA)': st.column_config.NumberColumn('Unique Spectral Counts (DDA)', help='Number of unique mass spectrometry spectral counts', format='%d'),
            'Razor Counts (DDA)': st.column_config.NumberColumn('Razor Counts (DDA)', help='Total razor spectral counts (shared peptides assigned to this protein)', format='%d'),
            'Protein Length': st.column_config.NumberColumn('Protein Length', help='Length of protein in amino acids', format='%d'),
            'Annotation Method': st.column_config.TextColumn('Annotation Method', help='Method used for annotation (MS = Mass Spectrometry, etc.)', max_chars=10),
            'TMT log2FC (50%)': st.column_config.NumberColumn('TMT log2FC (50%)', help='TMT log2 fold-change (AD vs Control) — 50% missing threshold', format='%.3f'),
            'TMT t-stat (50%)': st.column_config.NumberColumn('TMT t-stat (50%)', help='TMT t-statistic (AD vs Control) — 50% missing threshold', format='%.3f'),
            'TMT df (50%)': st.column_config.NumberColumn('TMT df (50%)', help='TMT t-test degrees of freedom — 50% missing threshold', format='%.1f'),
            'TMT CI low (50%)': st.column_config.NumberColumn('TMT CI low (50%)', help='TMT 95% CI lower bound on log2 fold-change — 50% missing threshold', format='%.3f'),
            'TMT CI high (50%)': st.column_config.NumberColumn('TMT CI high (50%)', help='TMT 95% CI upper bound on log2 fold-change — 50% missing threshold', format='%.3f'),
            "TMT Cohen's d (50%)": st.column_config.NumberColumn("TMT Cohen's d (50%)", help='TMT Cohen\'s d effect size (AD vs Control) — 50% missing threshold', format='%.3f'),
            'TMT p-val (50%)':  st.column_config.NumberColumn('TMT p-val (50%)',  help='TMT raw p-value — 50% missing threshold', format='%.4f'),
            'TMT q-val (50%)':  st.column_config.NumberColumn('TMT q-val (50%)',  help='TMT q-value (BH-adjusted) — 50% missing threshold', format='%.4f'),
            'TMT log2FC (0%)':  st.column_config.NumberColumn('TMT log2FC (0%)',  help='TMT log2 fold-change (AD vs Control) — 0% missing threshold (stringent)', format='%.3f'),
            'TMT t-stat (0%)': st.column_config.NumberColumn('TMT t-stat (0%)', help='TMT t-statistic (AD vs Control) — 0% missing threshold (stringent)', format='%.3f'),
            'TMT df (0%)': st.column_config.NumberColumn('TMT df (0%)', help='TMT t-test degrees of freedom — 0% missing threshold (stringent)', format='%.1f'),
            'TMT CI low (0%)': st.column_config.NumberColumn('TMT CI low (0%)', help='TMT 95% CI lower bound on log2 fold-change — 0% missing threshold (stringent)', format='%.3f'),
            'TMT CI high (0%)': st.column_config.NumberColumn('TMT CI high (0%)', help='TMT 95% CI upper bound on log2 fold-change — 0% missing threshold (stringent)', format='%.3f'),
            "TMT Cohen's d (0%)": st.column_config.NumberColumn("TMT Cohen's d (0%)", help='TMT Cohen\'s d effect size (AD vs Control) — 0% missing threshold (stringent)', format='%.3f'),
            'TMT p-val (0%)':   st.column_config.NumberColumn('TMT p-val (0%)',   help='TMT raw p-value — 0% missing threshold (stringent)', format='%.4f'),
            'TMT q-val (0%)':   st.column_config.NumberColumn('TMT q-val (0%)',   help='TMT q-value (BH-adjusted) — 0% missing threshold (stringent)', format='%.4f'),
            'MS Detect Control': st.column_config.NumberColumn('MS Detect Control', help='MS detection rate in Control donors', format='%.3f'),
            'MS Detect AD':      st.column_config.NumberColumn('MS Detect AD',      help='MS detection rate in AD donors', format='%.3f'),
            'ROSMAP log2FC': st.column_config.NumberColumn('ROSMAP log2FC', help='ROSMAP short-read RNA-seq log2 fold-change (AD vs Control)', format='%.3f'),
            'ROSMAP padj': st.column_config.NumberColumn('ROSMAP padj', help='ROSMAP short-read RNA-seq adjusted p-value', format='%.4f'),
            'MSBB log2FC': st.column_config.NumberColumn('MSBB log2FC', help='MSBB short-read RNA-seq log2 fold-change (AD vs Control)', format='%.3f'),
            'MSBB padj': st.column_config.NumberColumn('MSBB padj', help='MSBB short-read RNA-seq adjusted p-value', format='%.4f'),
            'ROSMAP Corr NonAD': st.column_config.NumberColumn('ROSMAP Corr NonAD', help='ROSMAP: correlation between Main ORF and smORF transcripts (non-AD donors)', format='%.3f'),
            'ROSMAP Corr AD': st.column_config.NumberColumn('ROSMAP Corr AD', help='ROSMAP: correlation between Main ORF and smORF transcripts (AD donors)', format='%.3f'),
            'MSBB Corr NonAD': st.column_config.NumberColumn('MSBB Corr NonAD', help='MSBB: correlation between Main ORF and smORF transcripts (non-AD donors)', format='%.3f'),
            'MSBB Corr AD': st.column_config.NumberColumn('MSBB Corr AD', help='MSBB: correlation between Main ORF and smORF transcripts (AD donors)', format='%.3f'),
            'RNA_LRT_Add_P': st.column_config.NumberColumn('RNA_LRT_Add_P', help='ROSMAP RNA-seq LRT additive model p-value', format='%.4f'),
            'RNA_LRT_Int_P': st.column_config.NumberColumn('RNA_LRT_Int_P', help='ROSMAP RNA-seq LRT interaction model p-value', format='%.4f'),
            'RP3 Default': st.column_config.TextColumn('RP3 Default', help='RP3 default pipeline result', max_chars=10),
            'RP3 MM+Amb': st.column_config.TextColumn('RP3 MM+Amb', help='RP3 multi-mapping + ambiguous result', max_chars=10),
            'RP3 Amb': st.column_config.TextColumn('RP3 Amb', help='RP3 ambiguous reads result', max_chars=10),
            'RP3 MM': st.column_config.TextColumn('RP3 MM', help='RP3 multi-mapping result', max_chars=10),
            'RiboCode': st.column_config.TextColumn('RiboCode', help='RiboCode ORF detection result', max_chars=10),
            'Tryptic Peptides': st.column_config.TextColumn('Tryptic Peptides', help='Tryptic peptide sequences', max_chars=60),
            'Tryptic Protein ID': st.column_config.TextColumn('Protein ID', help='Protein ID associated with the tryptic peptides', max_chars=15),
            'Tryptic Start Positions': st.column_config.TextColumn('Start Positions', help='Start positions of tryptic peptides', max_chars=30),
            'Tryptic End Positions': st.column_config.TextColumn('End Positions', help='End positions of tryptic peptides', max_chars=30),
            'Nt-Acetylated': st.column_config.CheckboxColumn('Nt-Acetylated', help='At least one tryptic peptide carries an N-terminal acetyl mark — direct evidence that this is a genuine protein N-terminus'),
            'Nt-Acetyl Peptides': st.column_config.TextColumn('Nt-Acetyl Peptides', help='Tryptic peptide sequence(s) observed with an N-terminal acetyl mark', max_chars=60),
            'Nt-Acetyl PSMs': st.column_config.NumberColumn('Nt-Acetyl PSMs', help='Total Nt-acetylated PSMs summed across the Nt-acetylated peptides', format='%d'),
            'Nt-Acetyl PSM Fraction': st.column_config.NumberColumn('Nt-Acetyl PSM Fraction', help='Highest per-peptide fraction of that peptide’s PSMs that carry the Nt-acetyl mark (1.0 = every PSM acetylated)', format='%.3f'),
            'Spectra Quality': st.column_config.TextColumn('Spectra Quality', help='Best Prosit spectral-match confidence tier (from the master Confidence column)', max_chars=20),
            'BLAST UniProt Match': st.column_config.TextColumn('BLAST UniProt Match', help='UniProt accession of the best BLASTp hit', max_chars=15),
            'BLAST % Match': st.column_config.NumberColumn('BLAST % Match', help='BLASTp percent identity to the best UniProt hit', format='%.1f'),
            'BLAST Aln Length': st.column_config.NumberColumn('BLAST Aln Length', help='BLASTp alignment length (residues)', format='%d'),
            'BLAST E-value': st.column_config.NumberColumn('BLAST E-value', help='BLASTp expectation value of the best UniProt hit', format='%.2e'),
            'BLAST Bit Score': st.column_config.NumberColumn('BLAST Bit Score', help='BLASTp bit score of the best UniProt hit', format='%.1f'),
            'UniProt Annotation': st.column_config.TextColumn('UniProt Annotation', help='UniProt curation annotation score (1\u20135\u2605) for reviewed / TrEMBL entries; blank for novel microproteins not in UniProt', max_chars=8),
            'UniProt Annotation Score': st.column_config.NumberColumn('UniProt Annotation Score', help='Numeric UniProt annotation score (1\u20135); use to sort by curation depth', format='%d'),
            'Kozak Strength': st.column_config.TextColumn('Kozak Strength', help='Full Kozak consensus context strength around the start codon (Salk smORFs only)', max_chars=10),
            'Kozak Downstream Strength': st.column_config.TextColumn('Kozak Downstream Strength', help='Kozak context strength considering only the +4 downstream position', max_chars=10),
            'Kozak Class': st.column_config.TextColumn('Kozak Class', help='Canonical Kozak strength classification', max_chars=10),
            'Kozak Window': st.column_config.TextColumn('Kozak Window', help='Local sequence context around the start codon (lowercase = flanking, uppercase = start + downstream)', max_chars=20),
            'NonATG Codon': st.column_config.TextColumn('NonATG Codon', help='The originally annotated (non-ATG) start codon', max_chars=6),
            'NonATG Valid Initiator': st.column_config.TextColumn('NonATG Valid Initiator', help='Whether the annotated codon is itself a plausible translation initiator', max_chars=6),
            'NonATG Has Supported Site': st.column_config.TextColumn('NonATG Has Supported Site', help='Whether current non-canonical initiation evidence/models support this annotated start actually being used', max_chars=6),
            'NonATG Context Strength': st.column_config.TextColumn('NonATG Context Strength', help='Kozak context strength around the annotated non-ATG codon', max_chars=10),
            'NonATG Alt Sites Found': st.column_config.NumberColumn('NonATG Alt Sites Found', help='Number of peptide-compatible alternative initiation sites found by scanning', format='%d'),
            'Optimal Codon Tier': st.column_config.TextColumn('Optimal Codon Tier', help='Strongest evidence tier among the alternative initiation candidates', max_chars=40),
            'NonATG Candidates': st.column_config.TextColumn('NonATG Candidates', help='Alternative initiation candidates as codon@position', max_chars=40),
        }
    )

    # ── Row click: navigate to that microprotein's entry page ──
    # display_df.index[...] maps the widget's *positional* selection back to the
    # label shared with filtered_df — the two frames diverge whenever a column
    # view sub-filters display_df (e.g. Alt-Initiation), so never index
    # filtered_df positionally here.
    selected_rows = selection.selection.rows if selection else []
    if selected_rows:
        _row = filtered_df.loc[display_df.index[selected_rows[0]]]
        st.query_params[ENTRY_QUERY_PARAM] = _entry_id(_row.get('sequence', ''))
        # scope="app": this function is a fragment, so the default fragment-scoped
        # rerun would never reach the routing check in main().
        st.rerun(scope="app")
    else:
        st.markdown('<div class="detail-panel"><div class="detail-prompt">Tick a row\'s select box above to open that microprotein\'s full entry page</div></div>', unsafe_allow_html=True)



@st.cache_data
def _read_static_file(path_str):
    """Read a static file's bytes for download (cached to avoid re-reads)."""
    p = Path(path_str)
    if p.exists():
        return p.read_bytes()
    return None


def _fmt_size(path):
    """Return human-readable file size string, or empty string if missing."""
    p = Path(path)
    if not p.exists():
        return ""
    n = p.stat().st_size
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024:
            return f"{n:.0f} {unit}"
        n /= 1024
    return f"{n:.1f} GB"


def _filtered_fasta(df):
    """Generate a FASTA string from the filtered dataframe's sequences."""
    lines = []
    for _, row in df.iterrows():
        seq = row.get("sequence", "")
        if not seq or pd.isna(seq):
            continue
        uid = row.get("UniProt_ID", "") or row.get("smORF_Coordinates", "") or "smORF"
        gene = row.get("Parent_Gene", "")
        header = f">{uid}"
        if gene and pd.notna(gene):
            header += f"|{gene}"
        lines.append(header)
        seq = str(seq).strip()
        for i in range(0, len(seq), 60):
            lines.append(seq[i:i + 60])
    return "\n".join(lines)



def _render_downloads_section(filtered_df, display_df):
    """Full downloads section rendered below the main table."""
    with st.expander("Downloads", expanded=False):

        # ── 1. Current View ──────────────────────────────────────────────────
        st.markdown(
            "<div style='font-size:0.8rem; font-weight:600; color:#8da8b8; "
            "text-transform:uppercase; letter-spacing:0.05em; margin-bottom:0.4rem;'>"
            "Current View</div>",
            unsafe_allow_html=True,
        )
        cv1, cv2, _ = st.columns([2, 2, 3])
        with cv1:
            st.download_button(
                label=f"Table (.csv) — {len(display_df):,} rows",
                data=display_df.to_csv(index=False),
                file_name=f"brain_microproteins_filtered_{len(display_df)}.csv",
                mime="text/csv",
                use_container_width=True,
                key="dl_filtered_csv",
            )
        with cv2:
            fasta_data = _filtered_fasta(filtered_df)
            fasta_count = fasta_data.count(">")
            st.download_button(
                label=f"Sequences (.fasta) — {fasta_count:,} entries",
                data=fasta_data,
                file_name=f"brain_microproteins_filtered_{fasta_count}.fasta",
                mime="text/plain",
                use_container_width=True,
                disabled=fasta_count == 0,
                key="dl_filtered_fasta",
            )

        st.divider()

        # ── 2. Genome Annotation ─────────────────────────────────────────────
        st.markdown(
            "<div style='font-size:0.8rem; font-weight:600; color:#8da8b8; "
            "text-transform:uppercase; letter-spacing:0.05em; margin-bottom:0.4rem;'>"
            "Genome Annotation — full dataset</div>",
            unsafe_allow_html=True,
        )

        _GENOME_FILES = [
            (
                "Unreviewed_Brain_Microproteins.fasta",
                "FASTA",
                "Amino acid sequences for all unreviewed brain microproteins",
                "text/plain",
            ),
            (
                "Unreviewed_Brain_Microproteins_CDS_Absent_from_UniProt.bed",
                "BED",
                "CDS genomic coordinates (hg38) — absent from UniProt",
                "text/plain",
            ),
            (
                "Unreviewed_Brain_Microproteins_Absent_from_UniProt.gtf",
                "GTF",
                "Gene transfer format annotations — absent from UniProt",
                "text/plain",
            ),
            (
                "Unreviewed_Brain_Microproteins_IDs.txt",
                "TXT",
                "Protein accession ID list",
                "text/plain",
            ),
            (
                "Unreviewed_Brain_Microproteins_genomic_coordinates.txt",
                "TXT",
                "Genomic coordinate list (chr:start-end)",
                "text/plain",
            ),
            (
                "Unreviewed_Brain_Microproteins_mapping_coordinates_to_sequences.tsv",
                "TSV",
                "Sequence ↔ coordinate mapping table",
                "text/tab-separated-values",
            ),
        ]

        for filename, fmt_tag, description, mime in _GENOME_FILES:
            fpath = GENOME_FILES_DIR / filename
            size_str = _fmt_size(fpath)
            data = _read_static_file(str(fpath))
            col_tag, col_name, col_size, col_btn = st.columns([1, 5, 1, 1])
            with col_tag:
                st.markdown(
                    f"<span style='background:#1e3a4f; color:#74c2e1; font-size:0.7rem; "
                    f"font-weight:700; padding:2px 6px; border-radius:4px; "
                    f"letter-spacing:0.04em;'>{fmt_tag}</span>",
                    unsafe_allow_html=True,
                )
            with col_name:
                st.markdown(
                    f"<span style='font-size:0.85rem;'>{filename}</span>  \n"
                    f"<span style='font-size:0.75rem; color:#8da8b8;'>{description}</span>",
                    unsafe_allow_html=True,
                )
            with col_size:
                st.markdown(
                    f"<span style='font-size:0.75rem; color:#8da8b8;'>{size_str}</span>",
                    unsafe_allow_html=True,
                )
            with col_btn:
                st.download_button(
                    label="Download",
                    data=data if data else b"",
                    file_name=filename,
                    mime=mime,
                    use_container_width=True,
                    disabled=data is None,
                    key=f"dl_genome_{filename}",
                )

        # Large combined GTF — repository only
        large_gtf = GENOME_FILES_DIR / "Ensembl_and_Unreviewed_Brain_Microproteins.gtf"
        large_size = _fmt_size(large_gtf)
        col_tag, col_name, col_size, col_btn = st.columns([1, 5, 1, 1])
        with col_tag:
            st.markdown(
                "<span style='background:#1e3a4f; color:#74c2e1; font-size:0.7rem; "
                "font-weight:700; padding:2px 6px; border-radius:4px; "
                "letter-spacing:0.04em;'>GTF</span>",
                unsafe_allow_html=True,
            )
        with col_name:
            st.markdown(
                "<span style='font-size:0.85rem;'>Ensembl_and_Unreviewed_Brain_Microproteins.gtf</span>  \n"
                "<span style='font-size:0.75rem; color:#8da8b8;'>Combined Ensembl + unreviewed annotations "
                "— available via project repository</span>",
                unsafe_allow_html=True,
            )
        with col_size:
            st.markdown(
                f"<span style='font-size:0.75rem; color:#8da8b8;'>{large_size}</span>",
                unsafe_allow_html=True,
            )
        with col_btn:
            st.markdown(
                "<span style='font-size:0.75rem; color:#8da8b8;'>too large</span>",
                unsafe_allow_html=True,
            )



def _fmt(val, fmt_str=None):
    """Format a value for display, returning '—' for NaN/None."""
    if val is None or (isinstance(val, float) and pd.isna(val)):
        return '—'
    if fmt_str:
        try:
            return fmt_str % val
        except (TypeError, ValueError):
            return str(val)
    return str(val)


def _not_na(val):
    """Return True if val is not None/NaN."""
    if val is None:
        return False
    if isinstance(val, float) and pd.isna(val):
        return False
    return True


def _sig_class_fc(val, sig_val=None, threshold=0.2):
    """CSS class for fold-change direction coloring (muted when not significant)."""
    if not _not_na(val):
        return 'sig-na'
    if sig_val is not None and _not_na(sig_val) and sig_val >= threshold:
        return 'sig-ns'
    return 'sig-up' if val > 0 else 'sig-down' if val < 0 else ''


def _sig_class_pval(val, threshold):
    """CSS class for p-value significance coloring."""
    if not _not_na(val):
        return 'sig-na'
    return 'sig-yes' if val < threshold else 'sig-ns'


# =============================================================================
# scRNA CELL-TYPE ENRICHMENT PANEL (single-microprotein heatmap + sig bars)
# =============================================================================
@st.cache_data(show_spinner=False)
def load_scrna_celltype_table(_all_path=str(SCRNA_ALL_CELLTYPES_FILE),
                              _sig_path=str(SCRNA_SUMMARY_FILE)):
    """Per-(microprotein, cell type) enrichment rows, trimmed to the plotted
    columns and indexed by sequence.

    Returns `(df, complete)`. `complete` is True when the all-cell-types export
    is present, i.e. the frame holds every cell type that was tested and the
    `significant` column separates the hits. When only the significant-only
    summary is shipped, everything in the frame is significant and `complete` is
    False, so the panel can say that cell types are missing rather than draw a
    partial picture as if it were the whole one.

    `sequence` is the join key the whole dashboard uses and it matches this
    table exactly (every sequence in it is in the master).
    """
    cols = ['sequence', 'gene', 'celltype', 'cell_type_general',
            'logFC', 'logCPM', 'p_adj.glb', 'p_adj.loc', 'significant']

    df, complete = None, False
    p_all = Path(_all_path)
    if p_all.exists():
        try:
            df = pd.read_csv(p_all, usecols=lambda c: c in cols)
            complete = True
        except Exception:
            df = None
    if df is None:
        p_sig = Path(_sig_path)
        if not p_sig.exists():
            return None, False
        try:
            df = pd.read_csv(p_sig, low_memory=False, usecols=lambda c: c in cols)
        except Exception:
            return None, False
        # This file is the p_adj.glb < 0.05 subset by construction.
        df['significant'] = True
    if 'sequence' not in df.columns:
        return None, False

    for c in ('logFC', 'logCPM', 'p_adj.glb', 'p_adj.loc'):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors='coerce')
    if 'significant' in df.columns:
        df['significant'] = df['significant'].astype(str).str.lower().isin(
            ('true', '1', 't', 'yes')) if df['significant'].dtype == object \
            else df['significant'].fillna(False).astype(bool)
    df['cell_type_general'] = df['cell_type_general'].fillna('Other').astype(str).str.strip()
    df['celltype'] = df['celltype'].astype(str).str.strip()
    return df.set_index('sequence').sort_index(), complete


def get_scrna_celltype_profile(seq):
    """Rows for one microprotein, ordered for plotting: grouped by
    cell_type_general in SCRNA_GENERAL_ORDER, alphabetical within a group.

    Returns `(prof, complete)`; `prof` is empty when the microprotein was never
    tested (no cell type reaches the model's expression floor)."""
    tbl, complete = load_scrna_celltype_table()
    if tbl is None or not seq:
        return None, complete
    try:
        prof = tbl.loc[[seq]]
    except KeyError:
        return pd.DataFrame(columns=tbl.columns), complete
    prof = prof.reset_index()
    _rank = {g: i for i, g in enumerate(SCRNA_GENERAL_ORDER)}
    prof['_g'] = prof['cell_type_general'].map(lambda g: _rank.get(g, len(_rank)))
    prof = prof.sort_values(['_g', 'celltype']).drop(columns='_g').reset_index(drop=True)
    return prof, complete


def _scrna_fc_limit(values):
    """Symmetric color limit for the log2FC tiles: the R heatmaps' ±0.3, grown
    to the next 0.1 above this microprotein's largest effect when it exceeds
    that. Symmetric so the white midpoint stays at logFC = 0."""
    values = np.asarray(values, dtype=float)
    vmax = float(np.nanmax(np.abs(values))) if values.size else 0.0
    if not np.isfinite(vmax):
        vmax = 0.0
    return max(SCRNA_FC_FLOOR, np.ceil(vmax * 10.0) / 10.0)


def _rgba(hex_color, alpha):
    """'#rrggbb' -> 'rgba(r,g,b,alpha)'. Per-bar alpha is how non-significant
    cell types are muted; go.Bar takes an array of colors but only a single
    scalar opacity."""
    h = hex_color.lstrip('#')
    r, g, b = (int(h[i:i + 2], 16) for i in (0, 2, 4))
    return f'rgba({r},{g},{b},{alpha})'


def _build_scrna_celltype_figure(prof, gene_label):
    """Two-panel figure: log2FC tiles (left) + -log10 FDR bars (right), one row
    per tested cell type. Height scales with the row count, so a microprotein
    tested in 3 cell types and one tested in 56 both get a readable panel.

    Significance is drawn, not filtered: non-significant cell types keep their
    tile and their bar (muted, no dot), so a microprotein that is up everywhere
    but only clears FDR in two cell types looks different from one that is
    genuinely flat."""
    n = len(prof)
    # Plotly stacks the first category at the bottom; reverse so the group
    # order above reads top-to-bottom, matching fct_rev() in the R heatmaps.
    prof = prof.iloc[::-1].reset_index(drop=True)
    ylabels = prof['celltype'].tolist()
    fc = prof['logFC'].to_numpy(dtype=float)
    padj = prof['p_adj.glb'].to_numpy(dtype=float)
    sig = (prof['significant'].to_numpy(dtype=bool)
           if 'significant' in prof.columns
           else np.ones(n, dtype=bool))
    # Scale the color to the significant effects when there are any: a
    # non-significant cell type can carry a large, noisy logFC that would
    # otherwise set the limit and flatten everything real into pale tiles.
    lim = _scrna_fc_limit(fc[sig] if sig.any() else fc)

    # -log10 of the global FDR — the value the R script gates on at 0.05.
    with np.errstate(divide='ignore'):
        neglog = -np.log10(np.clip(padj, 1e-300, None))
    neglog = np.where(np.isfinite(neglog), neglog, 0.0)

    hover = [
        (f"<b>{ct}</b> ({gen})<br>log2FC: {f:.3f}<br>"
         f"FDR (global): {pg:.3g}{'' if s else '  (n.s.)'}"
         f"<br>FDR (within cell type): {pl:.3g}"
         f"<br>logCPM: {cpm:.2f}")
        for ct, gen, f, pg, pl, cpm, s in zip(
            prof['celltype'], prof['cell_type_general'], fc, padj,
            prof['p_adj.loc'].to_numpy(dtype=float),
            prof['logCPM'].to_numpy(dtype=float), sig)
    ]

    height = int(min(1400, max(220, 22 * n + 130)))

    fig = make_subplots(rows=1, cols=2, shared_yaxes=True,
                        column_widths=[0.34, 0.66], horizontal_spacing=0.03,
                        subplot_titles=('log2FC', '−log10 FDR'))

    fig.add_trace(go.Heatmap(
        z=fc.reshape(-1, 1), x=['log2FC'], y=ylabels,
        zmin=-lim, zmax=lim, colorscale=SCRNA_FC_COLORSCALE, zmid=0,
        xgap=1, ygap=1,
        # Hover text must match z's (n, 1) shape, not z's flattened order.
        text=np.array(hover, dtype=object).reshape(-1, 1),
        hovertemplate='%{text}<extra></extra>',
        # Colorbar parked outside the plotting area — inside it would sit on
        # top of the significance bars. Sized in pixels, not axis fractions, so
        # a 1-cell-type panel doesn't squash the ticks into the title.
        colorbar=dict(title=dict(text='log2FC', side='right'), thickness=11,
                      lenmode='pixels', len=int(max(96, min(240, height * 0.5))),
                      x=1.015, xanchor='left', y=0.5, yanchor='middle',
                      tickfont=dict(size=9)),
    ), row=1, col=1)

    # A dot marks the cell types clearing FDR < 0.05 — the significance stars
    # convention, and the only mark on the tiles, since a heatmap cannot vary
    # opacity per cell.
    if sig.any():
        fig.add_trace(go.Scatter(
            x=['log2FC'] * int(sig.sum()),
            y=[y for y, s in zip(ylabels, sig) if s],
            mode='markers', showlegend=False,
            marker=dict(size=5, color='#111827', symbol='circle',
                        line=dict(width=1, color='rgba(255,255,255,0.85)')),
            customdata=[h for h, s in zip(hover, sig) if s],
            hovertemplate='%{customdata}<extra></extra>',
        ), row=1, col=1)

    # Bars carry the cell_type_general grouping as color (the R figure's facet
    # strips), one legend entry per group present; non-significant bars are
    # drawn at low alpha rather than dropped.
    _order = list(SCRNA_GENERAL_ORDER) + [
        g for g in prof['cell_type_general'].unique()
        if g not in SCRNA_GENERAL_ORDER  # new cluster labels upstream
    ]
    for gen in _order:
        m = (prof['cell_type_general'] == gen).to_numpy()
        if not m.any():
            continue
        base = SCRNA_GENERAL_COLORS.get(gen, '#999999')
        fig.add_trace(go.Bar(
            x=neglog[m], y=[y for y, keep in zip(ylabels, m) if keep],
            orientation='h', name=str(gen),
            marker=dict(color=[_rgba(base, 0.95 if s else 0.28)
                               for s, keep in zip(sig, m) if keep],
                        line=dict(width=0)),
            customdata=[h for h, keep in zip(hover, m) if keep],
            hovertemplate='%{customdata}<extra></extra>',
        ), row=1, col=2)

    # FDR = 0.05, the gate the manuscript's significant set is defined by. Bars
    # reaching past it are the cell types carrying a dot on the left.
    fig.add_vline(x=-np.log10(0.05), line=dict(color='rgba(255,255,255,0.5)',
                                               width=1, dash='dash'),
                  row=1, col=2)

    # Separator between cell_type_general blocks — the panel_spacing between
    # facets in the R heatmaps.
    _gen = prof['cell_type_general'].tolist()
    for i in range(len(_gen) - 1):
        if _gen[i] != _gen[i + 1]:
            for _c in (1, 2):
                fig.add_hline(y=i + 0.5, line=dict(color='rgba(255,255,255,0.22)',
                                                   width=1),
                              row=1, col=_c)

    fig.update_layout(
        height=height, barmode='overlay', bargap=0.25,
        # Left margin is left to automargin on the y axis: cluster labels run
        # from "Oli" to "InhPVALBCA8(Chandelier)" and a fixed margin either
        # clips the long ones or wastes half the panel on the short ones.
        margin=dict(l=8, r=76, t=48, b=44),
        paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)',
        font=dict(color='#e2e8f0', size=11),
        legend=dict(orientation='h', y=-0.02 - (18.0 / height), x=0,
                    font=dict(size=10), title=None,
                    bgcolor='rgba(0,0,0,0)'),
        showlegend=True,
    )
    fig.update_yaxes(tickfont=dict(size=9), showgrid=False, ticks='',
                     automargin=True, row=1, col=1)
    fig.update_xaxes(showticklabels=False, ticks='', showgrid=False, row=1, col=1)
    # Always keep the FDR = 0.05 line inside the frame, with headroom: for a
    # microprotein significant nowhere the bars stop short of it, and an axis
    # ending at the tallest bar would push the line onto the edge or off-frame.
    _xmax = float(np.nanmax(neglog)) if n else 0.0
    fig.update_xaxes(title=None, showgrid=True, gridcolor='rgba(255,255,255,0.08)',
                     zeroline=False, tickfont=dict(size=9),
                     range=[0, max(_xmax, -np.log10(0.05)) * 1.08],
                     row=1, col=2)
    for ann in fig.layout.annotations:
        ann.font.size = 11
    return fig


def _scrna_available(seq, row):
    """Whether the Single-Cell section has anything to show.

    Mirrors the early-return in _render_scrna_celltype_panel so the entry
    nav and the entry body agree on whether the section exists.
    """
    prof, _ = get_scrna_celltype_profile(seq)
    n = 0 if prof is None else len(prof)
    return bool(n) or _not_na(row.get('scRNA_logFC'))


def _render_scrna_celltype_panel(seq, row, boxed=True):
    """The Single-Cell RNA Enrichment section: headline numbers for the
    top cell type, then the per-cell-type panel."""
    prof, complete = get_scrna_celltype_profile(seq)
    scrna_fc = row.get('scRNA_logFC')

    if prof is None:
        # Table missing (not shipped / unreadable) — fall back to the master's
        # single best-cell-type summary rather than dropping the section.
        prof = pd.DataFrame()

    n = len(prof)
    n_sig = int(prof['significant'].sum()) if n and 'significant' in prof.columns else n
    if n == 0 and not _not_na(scrna_fc):
        return

    if not n:
        label = "Single-Cell RNA Enrichment"
    elif complete:
        label = (f"Single-Cell RNA Enrichment — significant in {n_sig} of "
                 f"{n} cell types tested")
    else:
        label = (f"Single-Cell RNA Enrichment — significant in {n_sig} cell type"
                 f"{'s' if n_sig != 1 else ''}")
    # On the entry page the section is already introduced by its own nav item,
    # so it renders open (boxed=False) with a plain header instead of nesting a
    # collapsed expander inside an expanded section.
    if boxed:
        _panel = st.expander(label, expanded=False)
    else:
        st.markdown(
            f'<div class="id-section-header">{html.escape(label)}</div>',
            unsafe_allow_html=True)
        _panel = st.container()

    with _panel:
        if _not_na(scrna_fc):
            scrna_pa = row.get('scRNA_padj')
            ct = _fmt(row.get('scRNA_celltype'))
            ct_gen = _fmt(row.get('scRNA_cell_type_general'))
            ct_display = f"{ct}" if ct_gen == '—' else f"{ct} ({ct_gen})"
            fc_cls = _sig_class_fc(scrna_fc, scrna_pa)
            st.markdown(f"""
        <div class="glass-card-section">
            <div class="id-section-header">Top Cell Type</div>
            <div style="display:grid; grid-template-columns: repeat(3, 1fr); gap:0.4rem 1.5rem;">
                <div><div class="id-field-label">logFC</div><div class="id-field-value {fc_cls}">{_fmt(scrna_fc, '%.4f')}</div></div>
                <div><div class="id-field-label">padj (global)</div><div class="id-field-value">{_fmt(scrna_pa, '%.4f')}</div></div>
                <div><div class="id-field-label">Cell Type</div><div class="id-field-value">{ct_display}</div></div>
            </div>
        </div>
        """, unsafe_allow_html=True)

        if n == 0:
            st.caption("Per-cell-type table not available for this microprotein.")
            return

        # Checkboxes rather than nested expanders — st.expander inside
        # st.expander is not a supported nesting in Streamlit.
        _plot_df = prof
        if complete and n_sig and n_sig < n:
            if st.checkbox("Significant cell types only",
                           key=f"scrna_ct_sig_{seq}"):
                _plot_df = prof[prof['significant']].reset_index(drop=True)

        gene = _fmt(row.get('Parent_Gene'))
        fig = _build_scrna_celltype_figure(_plot_df, gene)
        # Key off the sequence itself — hash() is salted per interpreter, so a
        # hashed key would change identity across restarts for no reason.
        st.plotly_chart(fig, use_container_width=True, key=f"scrna_ct_{seq}",
                        config={'displaylogo': False})

        _sig_mask = (_plot_df['significant'].to_numpy(dtype=bool)
                     if 'significant' in _plot_df.columns else None)
        _lim_vals = (_plot_df.loc[_sig_mask, 'logFC'] if _sig_mask is not None
                     and _sig_mask.any() else _plot_df['logFC'])
        _lim = _scrna_fc_limit(_lim_vals.to_numpy(dtype=float))
        st.caption(
            "Mathys et al. 2024 snRNA-seq, gpath association. Tiles are log2FC "
            f"(red = up with pathology, blue = down; scale clipped at ±{_lim:.1f}), "
            "bars are −log10 global FDR with the dashed line at FDR = 0.05. "
            + ("Every cell type the microprotein was tested in is shown; a dot "
               "and a full-opacity bar mark the ones passing FDR < 0.05, faded "
               "bars are tested but not significant."
               if complete else
               "Only cell types passing FDR < 0.05 are available in the shipped "
               "summary table, so cell types absent here were tested but not "
               "significant.")
        )
        if st.checkbox("Show values table", key=f"scrna_ct_tbl_{seq}"):
            _cols = ['cell_type_general', 'celltype', 'logFC', 'logCPM',
                     'p_adj.glb', 'p_adj.loc']
            if 'significant' in _plot_df.columns:
                _cols.append('significant')
            _t = _plot_df[_cols].rename(columns={
                'cell_type_general': 'General Cell Type', 'celltype': 'Cell Type',
                'logFC': 'log2FC', 'logCPM': 'logCPM',
                'p_adj.glb': 'FDR (global)', 'p_adj.loc': 'FDR (cell type)',
                'significant': 'Significant'})
            st.dataframe(_t, use_container_width=True, hide_index=True)


# =============================================================================
# ENTRY PAGE (UniProt-style detail view)
# =============================================================================
# Clicking a table row sets ?mp=<entry id> and reruns; main() then routes to
# _render_entry_page instead of the table. The page is a sticky left-hand
# section nav beside a single scrolling column of sections.
#
# Nav and body are both generated from ENTRY_SECTIONS, so a section that is
# hidden for lack of data cannot leave a dead link behind in the nav.


def _kv_section(title, fields, ncols=4):
    """Render one label/value grid — the shape shared by every entry section.

    fields: iterable of (label, value_html) or (label, value_html, extra_class),
    where extra_class is appended to `id-field-value` (e.g. a `sig-*` colour
    class, or `id-field-mono` for sequence-like values).
    """
    cells = []
    for f in fields:
        label, value = f[0], f[1]
        extra = f" {f[2]}" if len(f) > 2 and f[2] else ''
        cells.append(
            f'<div><div class="id-field-label">{label}</div>'
            f'<div class="id-field-value{extra}">{value}</div></div>'
        )
    st.markdown(
        f'<div class="glass-card-section">'
        f'<div class="id-section-header">{title}</div>'
        f'<div style="display:grid; grid-template-columns: repeat({ncols}, 1fr); '
        f'gap:0.4rem 1.5rem;">{"".join(cells)}</div>'
        f'</div>',
        unsafe_allow_html=True,
    )


def _entry_coords(row, ctx):
    """Genomic coordinates for this microprotein, or None.

    Prefers the smORF_Coordinates column (already on unified_df) and falls back
    to the sequence -> coordinates mapping file.
    """
    coords = row.get('smORF_Coordinates')
    if not _not_na(coords):
        seq = row.get('sequence')
        if seq and _not_na(seq):
            coords = ctx['seq_to_coords'].get(str(seq).strip())
    if not _not_na(coords):
        return None
    return str(coords).strip()


def _cartoon_src(coords):
    """Local path or remote URL for a smORF cartoon, or None if there isn't one."""
    if not coords or ':' not in coords:
        return None
    # chr10:73799111-73799263 -> chr10_73799111-73799263.png
    filename = coords.replace(':', '_') + '.png'
    if SMORF_CARTOON_DIR.exists():
        img_path = SMORF_CARTOON_DIR / filename
        if img_path.exists():
            return str(img_path)
    if filename in _build_cartoon_remote_set():
        return _remote_url(f"smorf_cartoon_figures/{filename}")
    return None


# ── Section availability predicates ──────────────────────────────────────────
# Each is cheap (column lookups plus cached index hits) and is called once for
# the nav and once for the body.

def _has_codon(row, ctx):
    return _not_na(row.get('Kozak_Strength')) or _not_na(row.get('NonATG_Annotated_Codon'))


def _has_tmt(row):
    return any(_not_na(row.get(c)) for c in ('TMT_log2fc', 'TMT_pvalue', 'TMT_qvalue'))


def _has_nt_acetyl(row):
    v = row.get('Nt_Acetyl_Peptides')
    return _not_na(v) and bool(str(v).strip())


def _nterm_sub_detail(row):
    """The row's N-terminal substitution records, always as a plain list.

    The column holds a list of dicts, but a parquet round-trip turns each one
    into a numpy object array — and `bool(arr)` raises for length > 1, so this
    must never be truth-tested directly (that crashed every entry with two or
    more N-terminal peptides). Length-checking is safe for list and ndarray
    alike; anything else (None, a stray NaN) degrades to empty.
    """
    detail = row.get('Nterm_Substitution_Detail')
    if detail is None:
        return []
    try:
        return list(detail)
    except TypeError:
        return []


def _has_nterm_substitution(row):
    """True only when some comparable N-terminal peptide actually differs.

    A comparable peptide with zero mismatches is the *expected* case (it is
    identical to the canonical protein, i.e. consistent with Met excision), so
    surfacing the section for those would put an all-'—' table on hundreds of
    entries. Only show it where there is something to see.
    """
    for d in _nterm_sub_detail(row):
        if not d.get('comparable'):
            continue
        n = pd.to_numeric(d.get('mismatch_count'), errors='coerce')
        if pd.notna(n) and n >= 1:
            return True
    return False


def _has_mirror(row, ctx):
    return bool(get_matching_mirror_plots(row.get('Tryptic_Peptides'), ctx['mirror_index']))


def _has_proteomics(row, ctx):
    return _has_tmt(row) or _has_nt_acetyl(row) or _has_nterm_substitution(row) or _has_mirror(row, ctx)


def _has_rna_stats(row):
    return any(_not_na(row.get(c)) for c in (
        'ROSMAP_log2FC', 'ROSMAP_pvalue', 'ROSMAP_padj',
        'MSBB_log2FC', 'MSBB_pvalue', 'MSBB_padj',
        'Nanopore_log2FC', 'Nanopore_pvalue', 'Nanopore_padj'))


def _has_expr_profile(row, ctx):
    coords = _entry_coords(row, ctx)
    return bool(coords and ctx['expression_index'].get(coords))


def _has_transcriptomics(row, ctx):
    return _has_rna_stats(row) or _has_expr_profile(row, ctx)


def _has_riboseq(row, ctx):
    return any(_not_na(row.get(c)) for c in
               ('RP3_Default', 'RP3_MM_Amb', 'RP3_Amb', 'RP3_MM', 'RiboCode'))


def _has_singlecell(row, ctx):
    return INCLUDE_SCRNA and _scrna_available(str(row.get('sequence', '')), row)


def _has_homology(row, ctx):
    return any(_not_na(row.get(c)) for c in (
        'BLAST_UniProt_Match', 'BLAST_Pct_Match', 'BLAST_Aln_Length',
        'BLAST_Evalue', 'BLAST_Bit'))


def _has_genome(row, ctx):
    if _not_na(row.get('UCSC_Link')):
        return True
    return _cartoon_src(_entry_coords(row, ctx)) is not None


def _always(row, ctx):
    return True


# ── Section renderers ────────────────────────────────────────────────────────
# The dispatcher in _render_entry_page emits each section's anchor <div>, so
# these only render content.

def _sec_overview(row, ctx):
    """Hero: identity line plus the facts you need before reading anything else."""
    gene = _fmt(row.get('Parent_Gene'))
    seq = str(row.get('sequence', ''))
    label = gene if gene != '—' else (seq[:30] + '…' if len(seq) > 30 else seq)
    badge = ('<span class="badge-swiss">Reviewed</span>' if ctx['is_swiss']
             else '<span class="badge-unreviewed">Unreviewed</span>')

    smorf_raw = _fmt(row.get('smORF_Class'))
    smorf = SMORF_DISPLAY_LABEL.get(smorf_raw, smorf_raw) if smorf_raw != '—' else '—'
    coords = _entry_coords(row, ctx) or '—'

    facts = [
        ('smORF Type', _fmt(row.get('smORF_Group'))),
        ('Subtype', smorf),
        ('Length', f"{_fmt(row.get('Protein_Length'), '%d')} aa"),
        ('Start Codon', _fmt(row.get('Start_Codon'))),
        ('MS Evidence', _fmt(row.get('MS_Evidence_Type'))),
        ('UniProt Annotation', _annotation_stars(row.get('UniProt_Annotation_Score'))),
    ]
    fact_html = ''.join(
        f'<div><div class="id-field-label">{lab}</div>'
        f'<div class="id-field-value">{val}</div></div>'
        for lab, val in facts
    )
    st.markdown(
        f'<div class="entry-hero{" is-swiss" if ctx["is_swiss"] else ""}">'
        f'<div class="entry-hero-title">'
        f'<span class="entry-gene">{html.escape(label)}</span>{badge}</div>'
        f'<div class="entry-hero-sub">{html.escape(coords)}</div>'
        f'<div class="entry-facts">{fact_html}</div>'
        f'</div>',
        unsafe_allow_html=True,
    )


def _sec_identity(row, ctx):
    smorf_raw = _fmt(row.get('smORF_Class'))
    smorf = SMORF_DISPLAY_LABEL.get(smorf_raw, smorf_raw) if smorf_raw != '—' else '—'
    _kv_section('Identity &amp; Annotation', [
        ('Parent Gene', _fmt(row.get('Parent_Gene'))),
        ('General smORF Type', _fmt(row.get('smORF_Group'))),
        ('smORF Subtype', smorf),
        ('Protein Length', f"{_fmt(row.get('Protein_Length'), '%d')} aa"),
        ('Start Codon', _fmt(row.get('Start_Codon'))),
        ('Annotation Method', _fmt(row.get('Annotation_Status'))),
        ('MS Evidence', _fmt(row.get('MS_Evidence_Type'))),
        ('DDA Grade', _fmt(row.get('DDA_Grade'))),
        ('Coordinates', _fmt(row.get('smORF_Coordinates')), 'id-field-mono'),
    ], ncols=4)


def _sec_codon(row, ctx):
    seq = str(row.get('sequence', ''))

    # ── Kozak context around the annotated start ──
    if _not_na(row.get('Kozak_Strength')):
        def _cap(col):
            return _fmt(row.get(col)).capitalize() if _not_na(row.get(col)) else '—'
        _kv_section('Kozak Context', [
            ('Full Kozak Strength', _cap('Kozak_Strength')),
            ('Downstream Strength', _cap('Kozak_Downstream_Strength')),
            ('Canonical Class', _cap('Kozak_Class_Canonical')),
            ('Kozak Window', _fmt(row.get('Kozak_Window')), 'id-field-mono'),
        ], ncols=4)

    # ── Non-AUG assessment + the alternative-initiation candidate table ──
    nonatg_codon = row.get('NonATG_Annotated_Codon')
    if not _not_na(nonatg_codon):
        return

    is_init = row.get('NonATG_Is_Initiator')
    has_support = row.get('NonATG_Has_Supported_Site')
    n_sites = row.get('NonATG_N_Sites')
    no_site_reason = _fmt(row.get('NonATG_No_Site_Reason'))
    _yn = {True: 'Yes', False: 'No'}

    _kv_section('Non-AUG Start Codon Assessment', [
        ('Annotated Codon', _fmt(nonatg_codon)),
        ('Valid Initiator?', _yn.get(is_init, _fmt(is_init))),
        ('Has Supported Site?', _yn.get(has_support, _fmt(has_support))),
        ('Context Strength', _fmt(row.get('NonATG_Context_Strength'))),
        ('Alternative Sites Found', _fmt(n_sites, '%d')),
    ], ncols=5)

    def _parse_list(val):
        try:
            parsed = ast.literal_eval(str(val))
            return parsed if isinstance(parsed, list) else [parsed]
        except (ValueError, SyntaxError):
            return []

    # Tryptic peptides observed by MS for this microprotein (as annotated) —
    # used to show which candidate sequence(s) reproduce the actual evidence.
    observed_peptides = [p for p in _parse_list(row.get('Tryptic_Peptides')) if p]

    def _supporting_peptide(pred_sequence):
        if not pred_sequence:
            return '—'
        for pep in observed_peptides:
            pos = pred_sequence.find(pep)
            if pos != -1:
                return f"{pep} (aa {pos + 1}–{pos + len(pep)})"
        return 'not detected' if observed_peptides else '—'

    # Row 1 is always the original, as-annotated translation (naive readthrough
    # from the called start, whatever its codon) — a fixed reference point
    # against which the scanned alternative candidates below can be compared.
    rows = [{
        'Candidate': 'Original (annotated)',
        'Position (scan)': 1,
        'Codon': _fmt(nonatg_codon),
        'Cognate Status': '—',
        'Tier': '—',
        'Initiator AA': seq[0] if seq else '—',
        'Predicted Sequence': seq,
        'Supporting Peptide': _supporting_peptide(seq),
    }]

    if _not_na(n_sites) and n_sites > 0:
        positions = _parse_list(row.get('NonATG_Codon_Position'))
        codons = _parse_list(row.get('NonATG_Codon'))
        cognate = _parse_list(row.get('NonATG_Cognate_Status'))
        tiers = _parse_list(row.get('NonATG_Tier_Name'))
        init_aa = _parse_list(row.get('NonATG_Initiator_AA'))
        pred_seq = _parse_list(row.get('NonATG_Predicted_Sequence'))
        n = max(len(positions), len(codons), len(cognate), len(tiers), len(init_aa), len(pred_seq))

        def _at(lst, i):
            return lst[i] if i < len(lst) else None
        for i in range(n):
            p_seq = _at(pred_seq, i) or ''
            rows.append({
                'Candidate': f'Alt {i + 1}',
                'Position (scan)': _at(positions, i),
                'Codon': _at(codons, i),
                'Cognate Status': _at(cognate, i),
                'Tier': _at(tiers, i),
                'Initiator AA': _at(init_aa, i),
                'Predicted Sequence': p_seq,
                'Supporting Peptide': _supporting_peptide(p_seq),
            })
    elif no_site_reason != '—':
        st.caption(f"No peptide-compatible alternative initiation site: {no_site_reason}")

    with st.expander(f"Alternative initiation site candidates ({len(rows) - 1})",
                     expanded=False):
        st.dataframe(pd.DataFrame(rows), width='stretch', hide_index=True)


def _sec_confidence(row, ctx):
    _kv_section('Confidence &amp; Conservation', [
        ('ShortStop Label', _fmt(row.get('ShortStop_Label'))),
        ('ShortStop Score', _fmt(row.get('ShortStop_Score'), '%.3f')),
        ('PhyloCSF Score', _fmt(row.get('PhyloCSF_Score'), '%.2f')),
        ('Unique Spectral Counts', _fmt(row.get('Unique_Spectral_Counts'), '%d')),
        ('UniProt Annotation', _annotation_stars(row.get('UniProt_Annotation_Score'))),
    ], ncols=4)


def _sec_proteomics(row, ctx):
    # ── TMT differential abundance ──
    if _has_tmt(row):
        tmt_fc, tmt_pv, tmt_qv = (row.get('TMT_log2fc'), row.get('TMT_pvalue'),
                                  row.get('TMT_qvalue'))
        _kv_section('TMT Proteomics (AD vs Control)', [
            ('log2 Fold Change', _fmt(tmt_fc, '%.4f'), _sig_class_fc(tmt_fc, tmt_qv)),
            ('p-value', _fmt(tmt_pv, '%.2e')),
            ('q-value (BH)', _fmt(tmt_qv, '%.4f'), _sig_class_pval(tmt_qv, 0.2)),
        ], ncols=3)

    # ── N-terminal acetylation (only if any acetylated peptide was observed) ──
    if _has_nt_acetyl(row):
        def _semis(val):
            if not _not_na(val):
                return []
            return [p.strip() for p in str(val).split(';')]

        def _plist(val):
            try:
                parsed = ast.literal_eval(str(val))
                return parsed if isinstance(parsed, list) else [parsed]
            except (ValueError, SyntaxError):
                return []

        _peps = _semis(row.get('Nt_Acetyl_Peptides'))
        _psms = _semis(row.get('Nt_Acetyl_N_PSMs'))
        _fracs = _semis(row.get('Nt_Acetyl_PSM_Fraction'))
        # Map each acetylated peptide back to its start position in the microprotein
        # (Tryptic_Peptides / Tryptic_Start_Positions are parallel list literals).
        _pos_map = dict(zip(_plist(row.get('Tryptic_Peptides')),
                            _plist(row.get('Tryptic_Start_Positions'))))

        def _at(lst, i):
            return lst[i] if i < len(lst) else None

        nt_rows = []
        for i, pep in enumerate(_peps):
            _f = pd.to_numeric(pd.Series([_at(_fracs, i)]), errors='coerce').iloc[0]
            _n = pd.to_numeric(pd.Series([_at(_psms, i)]), errors='coerce').iloc[0]
            nt_rows.append({
                'Acetylated Peptide': pep,
                'Start (aa)': _pos_map.get(pep),
                'Nt-Acetyl PSMs': None if pd.isna(_n) else int(_n),
                'PSM Fraction': '—' if pd.isna(_f) else f"{_f:.3f}",
            })

        _kv_section('N-terminal Acetylation (DDA)', [
            ('Acetylated Peptides', len(_peps)),
            ('Total Nt-Acetyl PSMs', _fmt(row.get('Nt_Acetyl_Total_PSMs'), '%d')),
            ('Best PSM Fraction', _fmt(row.get('Nt_Acetyl_Max_Fraction'), '%.3f')),
        ], ncols=3)
        st.dataframe(pd.DataFrame(nt_rows), width='stretch', hide_index=True)
        st.caption(
            "N-terminal acetylation is co-translational and marks a genuine protein "
            "N-terminus — an acetylated peptide starting at aa 1–2 is evidence for a real "
            "start site rather than initiator-Met excision of the parent ORF."
        )

    # ── N-terminal peptide vs. matched isoform (only if a comparable BLAST
    # re-alignment of an aa<=2 peptide was actually computed) ──
    if _has_nterm_substitution(row):
        detail = [d for d in _nterm_sub_detail(row) if d.get('comparable')]

        def _readable_subs(s):
            # substitutions is "pos:from>to" pairs, e.g. "5:R>G;9:S>T" -- render
            # as "R5G, S9T" so the position within the peptide stays visible.
            # A zero-mismatch peptide stores an empty value that comes back from
            # the CSV as NaN (which is truthy), so guard on _not_na and skip any
            # part that isn't a well-formed triple rather than unpacking blindly.
            if not _not_na(s):
                return '—'
            out = []
            for part in str(s).split(';'):
                part = part.strip()
                if ':' not in part or '>' not in part:
                    continue
                pos, change = part.split(':', 1)
                frm, to = change.split('>', 1)
                out.append(f"{frm}{pos}{to}")
            return ', '.join(out) if out else '—'

        sub_rows = [{
            'N-terminal Peptide': d.get('peptide_sequence'),
            'Start (aa)': d.get('peptide_start'),
            'Matched Isoform': d.get('uniprot_accession'),
            '% Identity': _fmt(d.get('pident'), '%.1f'),
            'Substitutions': _readable_subs(d.get('substitutions')),
        } for d in detail]

        n_rescue = sum(1 for d in detail if d.get('mismatch_count') in (1, 2))
        _kv_section('N-terminal Peptide vs. Matched Isoform', [
            ('Peptides Compared', len(detail)),
            ('With 1–2 Substitutions', n_rescue),
        ], ncols=2)
        st.dataframe(pd.DataFrame(sub_rows), width='stretch', hide_index=True)
        st.caption(
            "A peptide that differs from the matched UniProt isoform's sequence at this "
            "position can't be a bare Met-excision fragment of it — a substitution is direct "
            "evidence the peptide is genuinely distinct, even if it starts at aa 1–2."
        )

    # ── PROSIT spectral validation ──
    if _has_mirror(row, ctx):
        st.markdown('<div class="id-section-header">PROSIT Mirror Plots</div>',
                    unsafe_allow_html=True)
        _show_mirror_plots(row, ctx['mirror_index'], ctx['entry_id'])


def _sec_transcriptomics(row, ctx):
    def _de_section(title, prefix, fc_col, pv_col, pa_col, cpm_col=None):
        fc, pv, pa = row.get(fc_col), row.get(pv_col), row.get(pa_col)
        if not any(_not_na(v) for v in (fc, pv, pa)):
            return
        fields = [
            ('log2 Fold Change', _fmt(fc, '%.4f'), _sig_class_fc(fc, pa)),
            ('p-value', _fmt(pv, '%.2e')),
            ('padj', _fmt(pa, '%.4f'), _sig_class_pval(pa, 0.2)),
        ]
        if cpm_col:
            fields.append(('CPM', _fmt(row.get(cpm_col), '%.2f')))
        _kv_section(title, fields, ncols=len(fields))

    _de_section('Short-Read RNA-seq (ROSMAP)', 'ROSMAP',
                'ROSMAP_log2FC', 'ROSMAP_pvalue', 'ROSMAP_padj', 'ROSMAP_CPM')
    _de_section('Short-Read RNA-seq (MSBB)', 'MSBB',
                'MSBB_log2FC', 'MSBB_pvalue', 'MSBB_padj', 'MSBB_CPM')
    _de_section('Long-Read RNA-seq (Nanopore)', 'Nanopore',
                'Nanopore_log2FC', 'Nanopore_pvalue', 'Nanopore_padj')

    if _has_expr_profile(row, ctx):
        st.markdown('<div class="id-section-header">RNA Co-expression Profile</div>',
                    unsafe_allow_html=True)
        _show_expression_profile(row, ctx['expression_index'], ctx['seq_to_coords'],
                                 ctx['entry_id'])


def _sec_riboseq(row, ctx):
    _kv_section('Ribosome Profiling (RP3)', [
        ('RP3 Default', _fmt(row.get('RP3_Default'))),
        ('RP3 MM+Amb', _fmt(row.get('RP3_MM_Amb'))),
        ('RP3 Amb', _fmt(row.get('RP3_Amb'))),
        ('RP3 MM', _fmt(row.get('RP3_MM'))),
        ('RiboCode', _fmt(row.get('RiboCode'))),
    ], ncols=5)


def _sec_singlecell(row, ctx):
    _render_scrna_celltype_panel(str(row.get('sequence', '')), row, boxed=False)


def _sec_homology(row, ctx):
    acc = row.get('BLAST_UniProt_Match')
    if _not_na(acc):
        # Link out to the matched entry rather than printing a bare accession.
        _acc = html.escape(str(acc).strip())
        acc_html = (f'<a href="https://www.uniprot.org/uniprotkb/{quote(_acc)}" '
                    f'target="_blank" rel="noopener">{_acc}</a>')
    else:
        acc_html = '—'
    _kv_section('BLASTp Homology (vs UniProt)', [
        ('UniProt Match', acc_html),
        ('% Identity', _fmt(row.get('BLAST_Pct_Match'), '%.1f')),
        ('Aln Length', _fmt(row.get('BLAST_Aln_Length'), '%d')),
        ('E-value', _fmt(row.get('BLAST_Evalue'), '%.2e')),
        ('Bit Score', _fmt(row.get('BLAST_Bit'), '%.1f')),
    ], ncols=5)


def _sec_genome(row, ctx):
    coords = _entry_coords(row, ctx)
    src = _cartoon_src(coords)
    if src is not None:
        st.markdown('<div class="id-section-header">smORF Cartoon</div>',
                    unsafe_allow_html=True)
        st.caption(coords)
        _image_reserved(src)

    ucsc_link = row.get('UCSC_Link')
    if _not_na(ucsc_link):
        link = create_ucsc_link({'CLICK_UCSC': ucsc_link}, CUSTOM_UCSC_SESSION)
        if link:
            st.markdown(f"[View in UCSC Genome Browser]({link})")


def _sec_sequence(row, ctx):
    seq = str(row.get('sequence', ''))
    _kv_section('Sequence', [
        ('Length', f"{len(seq)} aa"),
        ('Start Codon', _fmt(row.get('Start_Codon'))),
        ('Initiator Residue', seq[0] if seq else '—'),
    ], ncols=3)
    st.code(seq, language=None)
    st.download_button(
        "Download FASTA",
        data=_filtered_fasta(row.to_frame().T),
        file_name=f"microprotein_{ctx['entry_id']}.fasta",
        mime="text/plain",
        key=f"dl_entry_fasta_{ctx['entry_id']}",
    )


# (anchor, nav title, renderer, availability predicate). The nav and the body
# are both built from this list, so they cannot disagree about which sections
# exist for a given microprotein.
# First option in the section picker, and therefore the default an entry opens
# on: the whole card is shown, and the sections below it narrow to one. A
# literal, not an f-string with a count in it — Streamlit derives widget
# identity from option text, so a changing label silently resets the selection.
ENTRY_SHOW_ALL = "All sections"

ENTRY_SECTIONS = [
    ('overview',        'Overview',                  _sec_overview,        _always),
    ('identity',        'Identity & Annotation',     _sec_identity,        _always),
    ('codon',           'Start Codon Context',       _sec_codon,           _has_codon),
    ('confidence',      'Confidence & Conservation', _sec_confidence,      _always),
    ('proteomics',      'Proteomics',                _sec_proteomics,      _has_proteomics),
    ('transcriptomics', 'Transcriptomics',           _sec_transcriptomics, _has_transcriptomics),
    ('riboseq',         'Ribosome Profiling',        _sec_riboseq,         _has_riboseq),
    ('singlecell',      'Single-Cell RNA',           _sec_singlecell,      _has_singlecell),
    ('homology',        'Homology',                  _sec_homology,        _has_homology),
    ('genome',          'Genome Context',            _sec_genome,          _has_genome),
    ('sequence',        'Sequence',                  _sec_sequence,        _always),
]


def _back_to_results():
    """Return to the table, clearing the stale row selection.

    st.dataframe keeps its selection in widget state, so a plain back-navigation
    would land on a table that still reports the row we just came from and
    immediately bounce forward again. Bumping the nonce gives the widget a fresh
    key, which resets the selection (and, unavoidably, the column sort — but
    only here, not during normal filtering).
    """
    st.session_state['_table_nonce'] = st.session_state.get('_table_nonce', 0) + 1
    st.query_params.clear()


def _render_entry_page(unified_df, entry_key, mirror_index, expression_index, seq_to_coords):
    """Full-page microprotein entry: sticky section nav beside the sections."""
    id_index = build_entry_id_index(unified_df, len(unified_df))
    seq = id_index.get(str(entry_key))
    matches = (unified_df[unified_df['sequence'].astype(str) == seq]
               if seq is not None else unified_df.iloc[0:0])

    if matches.empty:
        st.warning(
            f"No microprotein matches `{entry_key}`. The link may be stale — "
            "entry ids change if the underlying sequence changes."
        )
        if st.button("← Back to results", key="entry_back_missing"):
            _back_to_results()
            st.rerun()
        return

    row = matches.iloc[0]
    db_val = str(row.get('Database_All', row.get('Database', '')) or '').lower()
    ctx = {
        'entry_id': str(entry_key),
        'mirror_index': mirror_index,
        'expression_index': expression_index,
        'seq_to_coords': seq_to_coords,
        'is_swiss': 'swiss' in db_val,
    }

    if st.button("← Back to results", key="entry_back"):
        _back_to_results()
        st.rerun()

    sections = [s for s in ENTRY_SECTIONS if s[3](row, ctx)]

    nav_col, body_col = st.columns([1, 4], gap="medium")

    titles = [title for _a, title, _f, _av in sections]
    # "All sections" leads, so an entry opens on the full card (st.radio takes
    # its default from index 0) and the per-section options narrow from there.
    options = [ENTRY_SHOW_ALL] + titles

    with nav_col:
        # Marker div only — the CSS sticks and styles the element container
        # that follows it (the radio). Selecting re-renders the body instead of
        # scrolling to an anchor; see the stylesheet for why anchors cannot
        # work inside the Hugging Face wrapper's iframe.
        st.markdown('<div class="entry-nav-anchor"></div>', unsafe_allow_html=True)
        chosen = st.radio(
            "Section",
            options=options,
            key=f"entry_nav_{entry_key}",
            label_visibility="collapsed",
        )

    with body_col:
        show_all = chosen == ENTRY_SHOW_ALL
        for anchor, title, render_fn, _avail in sections:
            if not show_all and title != chosen:
                continue
            # Anchor target is emitted here, not by the renderers, so a new
            # section cannot be added without one. Kept for deep links even
            # though the picker no longer relies on scrolling.
            st.markdown(f'<div class="entry-section" id="sec-{anchor}"></div>',
                        unsafe_allow_html=True)
            render_fn(row, ctx)

        if show_all:
            # Runway so a deep link to the final section can still reach the
            # top where scrolling does work (direct URL / Streamlit Cloud).
            st.markdown('<div class="entry-tail"></div>', unsafe_allow_html=True)


def _show_mirror_plots(row, mirror_index, key_prefix):
    """Show Prosit mirror plot spectra for one microprotein."""
    tryptic_str = row.get('Tryptic_Peptides')

    gene = row.get('Parent_Gene', '')
    if gene and pd.notna(gene):
        st.caption(f"Gene: **{gene}**")

    try:
        all_peptides = ast.literal_eval(str(tryptic_str)) if tryptic_str and not pd.isna(tryptic_str) else []
        if isinstance(all_peptides, str):
            all_peptides = [all_peptides]
    except (ValueError, SyntaxError):
        all_peptides = []

    if all_peptides:
        badge_parts = []
        for pep in all_peptides:
            pep = pep.strip()
            if pep in mirror_index:
                q = mirror_index[pep]['best_quality']
                badge_parts.append(f"{QUALITY_EMOJI[q]} `{pep}`")
            else:
                badge_parts.append(f"— `{pep}`")
        st.markdown(" · ".join(badge_parts))

    all_plots = get_matching_mirror_plots(tryptic_str, mirror_index)
    if not all_plots:
        st.info("No Prosit mirror plot spectra available for this microprotein’s tryptic peptides.")
        return

    from collections import OrderedDict
    plots_by_peptide = OrderedDict()
    for p in all_plots:
        plots_by_peptide.setdefault(p['peptide'], []).append(p)

    peptide_list = list(plots_by_peptide.keys())
    if len(peptide_list) > 1:
        selected_peptide = st.selectbox(
            "Tryptic peptide:",
            options=peptide_list,
            format_func=lambda pep: f"{pep}  ({QUALITY_EMOJI[mirror_index[pep]['best_quality']]})",
            key=f"spectra_peptide_{key_prefix}"
        )
    else:
        selected_peptide = peptide_list[0]

    peptide_plots = plots_by_peptide[selected_peptide]
    if len(peptide_plots) > 1:
        plot_idx = st.selectbox(
            "Charge / scan:",
            options=range(len(peptide_plots)),
            format_func=lambda i: (
                f"z={peptide_plots[i]['charge']}+  "
                f"scan={peptide_plots[i]['scan']}  "
                f"{QUALITY_EMOJI[peptide_plots[i]['quality']]}"),
            key=f"spectra_plot_{key_prefix}"
        )
    else:
        plot_idx = 0

    plot_info = peptide_plots[plot_idx]
    st.markdown(
        f"**{QUALITY_EMOJI[plot_info['quality']]}** · "
        f"Peptide: `{plot_info['peptide']}` · "
        f"Charge: {plot_info['charge']}+ · "
        f"Scan: {plot_info['scan']}"
    )
    # Fixed pixel width (NOT width='stretch'): decouples the image size from the
    # container width so a vertical-scrollbar toggle can't resize the image and
    # re-trigger an infinite reflow loop (the "screen shake" on open).
    _image_reserved(plot_info['filepath'], width=760)


def _show_expression_profile(row, expression_index, seq_to_coords, key_prefix):
    """Show RNA co-expression profile for one microprotein."""
    if not expression_index:
        st.info("Expression profile index not available.")
        return

    coords = _entry_coords(row, {'seq_to_coords': seq_to_coords})
    if not coords:
        st.info("Genomic coordinates not found for this smORF.")
        return

    profiles = expression_index.get(coords, [])
    if not profiles:
        st.info("No expression profile available for this smORF.")
        return

    # If both coupled and non_coupled exist, use a radio toggle
    if len(profiles) > 1:
        coupling_options = {p['coupling']: i for i, p in enumerate(profiles)}
        chosen = st.radio(
            "Co-expression type:",
            options=list(coupling_options.keys()),
            format_func=lambda c: 'Coupled' if c == 'coupled' else 'Non-coupled',
            horizontal=True,
            key=f"expr_profile_{key_prefix}"
        )
        selected_idx = coupling_options[chosen]
    else:
        selected_idx = 0

    profile = profiles[selected_idx]
    coupling_label = "Coupled" if profile['coupling'] == 'coupled' else "Non-coupled"
    st.caption(f"{coupling_label} · {profile['gene']} · {coords}")
    fp = profile['filepath']
    # Both renditions reserve their height before the figure arrives, so the
    # coupled/non-coupled toggle and the entry-nav anchors stay put; st.image
    # for PNGs (supports expand), iframe for PDFs.
    if fp.endswith('.png'):
        _image_reserved(fp)
    else:
        _pdf_reserved(fp)


@st.cache_data(show_spinner=False)
def _build_cartoon_remote_set():
    """Set of cartoon filenames available remotely on HF (fallback only)."""
    files = set()
    for rel in _list_remote_files():
        if rel.startswith("smorf_cartoon_figures/") and rel.endswith(".png"):
            files.add(rel.split("/", 1)[1])
    return files


if __name__ == "__main__":
    main()
