import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path
from urllib.parse import quote
import numpy as np
import hashlib
import os
import ast
import re

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


# =============================================================================
# PROSIT MIRROR PLOT CONFIGURATION
# =============================================================================
MIRROR_PLOT_BASE = Path(__file__).parent / "mirror_plots"
QUALITY_LEVELS = ['Strong', 'Moderate', 'Weak', 'Insufficient', 'No PSM']
QUALITY_EMOJI = {
    'Strong': '🟢 Strong',
    'Moderate': '🟡 Moderate',
    'Weak': '🟠 Weak',
    'Insufficient': '🔴 Insufficient',
    'No PSM': '⚫ No PSM',
}

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
SMORF_PARENT_GROUPS = {
    'Upstream': ['uORF', 'uoORF', 'uaORF', 'uaoORF', 'udORF'],
    'Downstream': ['dORF', 'doORF', 'daORF', 'daoORF'],
    'Internal ORF': ['iORF', 'oORF'],
    'TrEMBL/AltORF': ['eORF', 'TrEMBL'],
    'Short-Isoform': ['Iso', 'D-Iso', 'N-Iso'],
    'lncRNA': ['lncRNA'],
    'psORF': ['psORF'],
}
SMORF_CHILD_TO_PARENT = {
    child: parent
    for parent, children in SMORF_PARENT_GROUPS.items()
    for child in children
}
SMORF_DISPLAY_LABEL = {}

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
        return 'No PSM'
    if not mirror_index:
        return 'Insufficient'
    try:
        peptides = ast.literal_eval(str(tryptic_peptides_str))
        if isinstance(peptides, str):
            peptides = [peptides]
    except (ValueError, SyntaxError):
        peptides = []
    best = 'Insufficient'
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
    layout="wide"
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
/* When detail panel is open, square the table bottom */
.table-has-detail .stDataFrame {{
    border-radius: 12px 12px 0 0 !important;
    border-bottom: none !important;
    margin-bottom: 0 !important;
}}

/* ── Docked detail panel ── */
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
</style>
""", unsafe_allow_html=True)


# =============================================================================
# PASSWORD PROTECTION
# =============================================================================
def check_password():
    """Returns True if the user entered the correct password."""
    def password_entered():
        if hashlib.sha256(st.session_state["password"].encode()).hexdigest() == os.environ.get(
            "DASHBOARD_PASSWORD_HASH",
            hashlib.sha256("revision".encode()).hexdigest()
        ):
            st.session_state["password_correct"] = True
            del st.session_state["password"]
        else:
            st.session_state["password_correct"] = False

    if st.session_state.get("password_correct", False):
        return True

    st.markdown("""
    <div class="login-card">
        <div style="font-size: 2.5rem; margin-bottom: 0.5rem;"></div>
        <h2>Brain Microproteins Dashboard</h2>
        <p>Saghatelian Lab &middot; Salk Institute for Biological Studies</p>
    </div>
    """, unsafe_allow_html=True)
    st.text_input("Password", type="password", on_change=password_entered, key="password")
    if "password_correct" in st.session_state and not st.session_state["password_correct"]:
        st.error("Password incorrect.")
    return False


if not check_password():
    st.stop()


# =============================================================================
# DATA LOADING
# =============================================================================
def load_analysis_results():
    """Load all analysis result files."""
    base_path = Path(__file__).parent
    analyses = {
        "Annotation Summary": {
            "path": base_path / "Annotations" / "Brain_Microproteins_Discovery_summary.csv",
            "type": "unreviewed_only"
        },
        "Proteomics (TMT)": {
            "path": base_path / "Proteomics" / "Proteomics_Results_summary.csv",
            "type": "mixed"
        },
        "Proteomics + RiboSeq (RP3)": {
            "path": base_path / "RP3" / "RP3_Results_summary.csv",
            "type": "unreviewed_only"
        },
        "Short-Read RNA in AD": {
            "path": base_path / "Transcriptomics" / "Short-Read_Transcriptomics_Results_summary.csv",
            "type": "mixed"
        },
        "Long-Read RNA in AD": {
            "path": base_path / "Transcriptomics" / "Long-Read_Transcriptomics_Results_summary.csv",
            "type": "unreviewed_only"
        },
        "scRNA Enrichment": {
            "path": base_path / "scRNA_Enrichment" / "scRNA_Enrichment_summary.csv",
            "type": "scrna"
        } if INCLUDE_SCRNA else None,
        "ShortStop Classification": {
            "path": base_path / "Annotations" / "ShortStop_Microproteins_summary.csv",
            "type": "unreviewed_only"
        },
    }
    # Remove any entries that were conditionally set to None
    return {k: v for k, v in analyses.items() if v is not None}


@st.cache_data(show_spinner=False)
def load_and_merge_all_data():
    """Load and merge all CSV files by sequence."""
    analysis_files = load_analysis_results()
    master_df = pd.DataFrame()

    current_dir = Path(__file__).resolve().parent
    tryptic_path = current_dir.parent / "Code" / "data" / "cleaned_tryptic_peptides_under_151aa.csv"

    if tryptic_path.exists():
        try:
            tryptic_raw = pd.read_csv(tryptic_path, low_memory=False)
            tryptic_df = tryptic_raw.rename(columns={
                'peptide_sequence': 'Tryptic_Peptides',
                'start': 'Tryptic_Start_Positions',
                'end': 'Tryptic_End_Positions',
            })
            if 'protein_id' in tryptic_df.columns:
                tryptic_df = tryptic_df.drop('protein_id', axis=1)
            tryptic_df['Tryptic_Peptides_present'] = True
            master_df = tryptic_df.copy()
        except Exception as e:
            st.error(f"Could not load tryptic peptides data: {e}")

    for analysis_name, info in analysis_files.items():
        if info['path'].exists():
            try:
                df = pd.read_csv(info['path'], low_memory=False)
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

    # TrEMBL entries have no smorf_type — assign synthetic class so they appear in the group filter
    if 'Database' in master_df.columns:
        _is_trembl = master_df['Database'].fillna('').str.contains('TrEMBL', case=False, na=False)
        ud.loc[ud['smORF_Class'].isna() & _is_trembl, 'smORF_Class'] = 'TrEMBL'

    # Display label (no remapping currently; hook for future overrides)
    ud['smORF_Display'] = ud['smORF_Class'].map(
        lambda v: SMORF_DISPLAY_LABEL.get(str(v), v) if pd.notna(v) else v
    )
    # Parent group (Upstream, Downstream, TrEMBL/AltORF, etc.)
    ud['smORF_Group'] = ud['smORF_Class'].map(
        lambda v: SMORF_CHILD_TO_PARENT.get(str(v), str(v)) if pd.notna(v) else v
    )

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

    # Start Codon
    cols = [c for c in master_df.columns if 'start' in c.lower() and 'codon' in c.lower()]
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

    # Tryptic peptides
    for col in ['Tryptic_Peptides', 'Tryptic_Protein_ID', 'Tryptic_Start_Positions', 'Tryptic_End_Positions']:
        if col in master_df.columns:
            ud[col] = master_df[col]

    # --- TMT proteomics stats ---
    for src, dest in [
        ('TMT_log2fc_50pct_missing', 'TMT_log2fc_50pct'),
        ('TMT_pvalue_50pct_missing', 'TMT_pvalue_50pct'),
        ('TMT_qvalue_50pct_missing', 'TMT_qvalue_50pct'),
        ('TMT_log2fc_0pct_missing',  'TMT_log2fc_0pct'),
        ('TMT_pvalue_0pct_missing',  'TMT_pvalue_0pct'),
        ('TMT_qvalue_0pct_missing',  'TMT_qvalue_0pct'),
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

    # Significance indicators (kept for sidebar filter + metrics)
    ud['TMT_Significant'] = ud.get('TMT_qvalue_50pct', pd.Series(dtype=float)).fillna(1.0) < 0.2
    ud['TMT_Highly_Significant'] = ud.get('TMT_qvalue_50pct', pd.Series(dtype=float)).fillna(1.0) < 0.05
    ud['ROSMAP_Significant'] = ud.get('ROSMAP_padj', pd.Series(dtype=float)).fillna(1.0) < 0.2
    ud['ROSMAP_Highly_Significant'] = ud.get('ROSMAP_padj', pd.Series(dtype=float)).fillna(1.0) < 0.05

    return ud


# =============================================================================
# DISK-PARQUET CACHE FOR MERGED + EXTRACTED DATAFRAME
# =============================================================================
# The full merge+extract pipeline reads ~16 MB across 8 CSVs and runs many
# bfill/regex passes (~3-8 s cold). Caching the final unified_df to a single
# parquet on disk makes subsequent cold starts <0.5 s. Cache key is the mtime
# fingerprint of every source CSV — auto-invalidates when any input changes.

def _source_csv_paths():
    base = Path(__file__).parent
    return [
        base.parent / "Code" / "data" / "cleaned_tryptic_peptides_under_151aa.csv",
        base / "Annotations" / "Brain_Microproteins_Discovery_summary.csv",
        base / "Proteomics" / "Proteomics_Results_summary.csv",
        base / "RP3" / "RP3_Results_summary.csv",
        base / "Transcriptomics" / "Short-Read_Transcriptomics_Results_summary.csv",
        base / "Transcriptomics" / "Long-Read_Transcriptomics_Results_summary.csv",
        base / "scRNA_Enrichment" / "scRNA_Enrichment_summary.csv",
        base / "Annotations" / "ShortStop_Microproteins_summary.csv",
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
    ) + "|" + mirror_fp

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
    # _Spectra_Quality (depends on mirror_index)
    if 'Tryptic_Peptides' in ud.columns:
        ud['_Spectra_Quality'] = ud['Tryptic_Peptides'].apply(
            lambda x: get_spectra_quality(x, mirror_index)
        )
    else:
        ud['_Spectra_Quality'] = 'No PSM'

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

    return ud


# =============================================================================
# MAIN APP — SINGLE-PAGE LAYOUT
# =============================================================================
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
            ).apply(lambda x: 'Insufficient' if (x and not pd.isna(x)) else 'No PSM')

    # ── Sidebar: ALL filters ──
    with st.sidebar:
        st.markdown('<div style="font-size:1.05rem; font-weight:700; color:#ffffff; margin-bottom:0.7rem; padding-bottom:0.4rem; border-bottom:1px solid rgba(116,162,183,0.25);">Filters</div>', unsafe_allow_html=True)
        st.markdown('<div class="sidebar-section-header">Search</div>', unsafe_allow_html=True)

        gene_search = st.text_input(
            "Gene Search",
            placeholder="e.g., BCL3, RAB3C\u2026",
            help="Search by parent gene name (partial match, case-insensitive)"
        )

        seq_search = st.text_input(
            "Sequence Search",
            placeholder="e.g., MAASGK\u2026",
            help="Search by amino acid sequence (substring match, case-insensitive)"
        )

        st.markdown('<div class="sidebar-section-header">smORF Type</div>', unsafe_allow_html=True)

        if 'smORF_Group' in unified_df.columns:
            group_options = sorted(
                [g for g in SMORF_PARENT_GROUPS.keys()
                 if unified_df['smORF_Group'].eq(g).any()],
                key=lambda g: list(SMORF_PARENT_GROUPS.keys()).index(g)
            )
            selected_groups = st.multiselect(
                "General smORF Type", options=group_options, default=[],
                help="Top-level smORF category. Leave empty for all."
            )
            # ── Downstream sub-type filter (only when Downstream group is selected) ──
            selected_downstream_sub = []
            if 'Downstream' in selected_groups and 'smORF_Class' in unified_df.columns:
                downstream_in_data = sorted([
                    t for t in SMORF_PARENT_GROUPS['Downstream']
                    if unified_df['smORF_Class'].eq(t).any()
                ])
                if downstream_in_data:
                    selected_downstream_sub = st.multiselect(
                        "Downstream Sub-type",
                        options=downstream_in_data,
                        default=[],
                        help="Refine within Downstream ORFs (doORF, daORF, etc.)"
                    )
        else:
            selected_groups = []
            selected_downstream_sub = []

        st.markdown('<div class="sidebar-section-header">Evidence &amp; Quality</div>', unsafe_allow_html=True)

        if 'Database_All' in unified_df.columns or 'Database' in unified_df.columns:
            unreviewed_only = st.checkbox(
                "Unreviewed Only",
                value=False,
                help="Show only unreviewed entries (Salk + TrEMBL); unchecked shows reviewed + unreviewed"
            )
        else:
            unreviewed_only = False

        if 'Annotation_Status' in unified_df.columns:
            ann_options = sorted(unified_df['Annotation_Status'].dropna().unique().tolist())
            selected_annotation = st.multiselect(
                "Evidence Type", options=ann_options, default=[],
                help="MS = detected by mass spectrometry; Ribo-Seq/ShortStop = ribosome profiling only; MS+Ribo = both"
            )
        else:
            selected_annotation = []

        quality_options = [q for q in QUALITY_LEVELS if q in unified_df['_Spectra_Quality'].values]
        selected_quality = st.multiselect(
            "Spectra Quality", options=quality_options, default=[],
            help="Filter by best Prosit mirror-plot spectral match quality"
        )

        st.markdown('<div class="sidebar-section-header">Significance</div>', unsafe_allow_html=True)

        tmt_sig_choice = "Any"
        has_tmt_50 = 'TMT_qvalue_50pct' in unified_df.columns
        has_tmt_0  = 'TMT_qvalue_0pct'  in unified_df.columns
        if has_tmt_50 or has_tmt_0:
            tmt_options = ["Any"]
            if has_tmt_50:
                tmt_options += [
                    "Exploratory Tier (FDR < 0.2) — ≥50% samples (stringent)",
                    "Significant Tier (FDR < 0.05) — ≥50% samples (stringent)",
                ]
            if has_tmt_0:
                tmt_options += [
                    "Exploratory Tier (FDR < 0.2) — ≥1 sample/condition",
                    "Significant Tier (FDR < 0.05) — ≥1 sample/condition",
                ]
            tmt_sig_choice = st.radio(
                "TMT-MS Significance",
                options=tmt_options,
                index=0,
                help=(
                    "Filter by TMT proteomics q-value at the chosen stringency. "
                    "≥50% samples (stringent): peptide quantified in at least half the donors per condition. "
                    "≥1 sample/condition: peptide quantified in at least one donor per condition (lenient)."
                ),
            )

        rna_sig_choice = "Any"
        if 'ROSMAP_Significant' in unified_df.columns:
            rna_sig_choice = st.radio(
                "RNA Significance",
                options=["Any", "Exploratory Tier (FDR < 0.2)", "Significant Tier (FDR < 0.05)"],
                index=0,
                help="Filter by ROSMAP bulk RNA-seq adjusted p-value",
            )

        only_tmt_sig = tmt_sig_choice == "Exploratory Tier (FDR < 0.2) — ≥50% samples (stringent)"
        only_tmt_highly_sig = tmt_sig_choice == "Significant Tier (FDR < 0.05) — ≥50% samples (stringent)"
        only_tmt_sig_strict = tmt_sig_choice == "Exploratory Tier (FDR < 0.2) — ≥1 sample/condition"
        only_tmt_highly_sig_strict = tmt_sig_choice == "Significant Tier (FDR < 0.05) — ≥1 sample/condition"
        only_rna_sig = rna_sig_choice == "Exploratory Tier (FDR < 0.2)"
        only_rna_highly_sig = rna_sig_choice == "Significant Tier (FDR < 0.05)"

        st.markdown('<div class="sidebar-section-header">Advanced</div>', unsafe_allow_html=True)

        with st.expander("Advanced Filters", expanded=False):
            length_range = None
            if 'Protein_Length' in unified_df.columns:
                len_data = unified_df['Protein_Length'].dropna()
                if len(len_data) > 0:
                    l_min, l_max = int(len_data.min()), int(len_data.max())
                    length_range = st.slider("Protein Length (aa)",
                                             min_value=l_min, max_value=l_max,
                                             value=(l_min, l_max))

            spectral_range = None
            if 'Unique_Spectral_Counts' in unified_df.columns:
                spec_data = unified_df['Unique_Spectral_Counts'].dropna()
                if len(spec_data) > 0:
                    s_min, s_max = int(spec_data.min()), int(spec_data.max())
                    spectral_range = st.slider("Unique Spectral Counts",
                                               min_value=s_min, max_value=s_max,
                                               value=(s_min, s_max))

            phylocsf_range = None
            if 'PhyloCSF_Score' in unified_df.columns:
                phylo_data = unified_df['PhyloCSF_Score'].dropna()
                if len(phylo_data) > 0:
                    p_min, p_max = float(phylo_data.min()), float(phylo_data.max())
                    phylocsf_range = st.slider("PhyloCSF Score",
                                               min_value=p_min, max_value=p_max,
                                               value=(p_min, p_max), step=0.5)

            if 'Start_Codon' in unified_df.columns:
                codon_options = sorted(unified_df['Start_Codon'].dropna().unique().tolist())
                selected_codons = st.multiselect(
                    "Start Codon", options=codon_options, default=[],
                    help="Leave empty to include all codons (ATG and near-cognate)"
                )
            else:
                selected_codons = []

            if 'ShortStop_Label' in unified_df.columns:
                ss_options = sorted(unified_df['ShortStop_Label'].dropna().unique().tolist())
                selected_shortstop = st.multiselect(
                    "ShortStop Label", options=ss_options, default=[],
                    help="ShortStop ribosome profiling classification"
                )
            else:
                selected_shortstop = []

            exclude_nterm = st.checkbox(
                "Exclude N-term M excision",
                help="Hide microproteins whose only tryptic peptides start at position 1 or 2"
            )

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
            "<a href='https://github.com/brendan-miller-salk/ad-dark-microproteome' "
            "target='_blank' style='color:#74c2e1; text-decoration:none;'>"
            "GitHub Repository</a></div>",
            unsafe_allow_html=True,
        )

    # ── Apply filters ──
    filtered_df = unified_df.copy()

    if gene_search and 'Parent_Gene' in filtered_df.columns:
        filtered_df = filtered_df[
            filtered_df['Parent_Gene'].str.contains(gene_search, case=False, na=False)
        ]

    if seq_search and 'sequence' in filtered_df.columns:
        filtered_df = filtered_df[
            filtered_df['sequence'].str.contains(seq_search.upper(), case=False, na=False)
        ]

    # smORF group + downstream sub-type filter
    if selected_groups and 'smORF_Class' in filtered_df.columns:
        allowed_types = []
        for group in selected_groups:
            allowed_types.extend(SMORF_PARENT_GROUPS.get(group, []))
        # Narrow downstream to sub-type selection if specified
        if selected_downstream_sub:
            downstream_all = set(SMORF_PARENT_GROUPS['Downstream'])
            allowed_types = [
                t for t in allowed_types
                if t not in downstream_all or t in selected_downstream_sub
            ]
        filtered_df = filtered_df[filtered_df['smORF_Class'].isin(allowed_types)]
    elif selected_downstream_sub and 'smORF_Class' in filtered_df.columns:
        # No group filter but downstream sub-type selected — filter downstream only
        filtered_df = filtered_df[filtered_df['smORF_Class'].isin(selected_downstream_sub)]

    if selected_annotation and 'Annotation_Status' in filtered_df.columns:
        filtered_df = filtered_df[filtered_df['Annotation_Status'].isin(selected_annotation)]

    if selected_codons and 'Start_Codon' in filtered_df.columns:
        filtered_df = filtered_df[
            filtered_df['Start_Codon'].isin(selected_codons) | filtered_df['Start_Codon'].isna()
        ]

    if selected_shortstop and 'ShortStop_Label' in filtered_df.columns:
        filtered_df = filtered_df[filtered_df['ShortStop_Label'].isin(selected_shortstop)]

    if spectral_range and 'Unique_Spectral_Counts' in filtered_df.columns:
        filtered_df = filtered_df[
            filtered_df['Unique_Spectral_Counts'].between(spectral_range[0], spectral_range[1]) | filtered_df['Unique_Spectral_Counts'].isna()
        ]

    if length_range and 'Protein_Length' in filtered_df.columns:
        filtered_df = filtered_df[
            filtered_df['Protein_Length'].between(length_range[0], length_range[1]) | filtered_df['Protein_Length'].isna()
        ]

    if phylocsf_range and 'PhyloCSF_Score' in filtered_df.columns:
        filtered_df = filtered_df[
            filtered_df['PhyloCSF_Score'].between(phylocsf_range[0], phylocsf_range[1]) | filtered_df['PhyloCSF_Score'].isna()
        ]

    if selected_quality:
        filtered_df = filtered_df[filtered_df['_Spectra_Quality'].isin(selected_quality)]

    if unreviewed_only:
        _db_col = 'Database_All' if 'Database_All' in filtered_df.columns else ('Database' if 'Database' in filtered_df.columns else None)
        if _db_col:
            _db_lower = filtered_df[_db_col].fillna('').astype(str).str.lower()
            _is_swiss = _db_lower.str.contains('swiss', na=False)
            filtered_df = filtered_df[~_is_swiss]

    if exclude_nterm and 'Tryptic_Start_Positions' in filtered_df.columns:
        # N-term M excision filter only applies to Isoform smORF types
        _iso_types = set(SMORF_PARENT_GROUPS.get('Short-Isoform', []))
        def _has_non_nterm(row):
            # Only apply N-term excision filter to Iso smORF types
            smorf = row.get('smORF_Class')
            if pd.isna(smorf) or str(smorf) not in _iso_types:
                return True
            pos_str = row.get('Tryptic_Start_Positions')
            if not pos_str or pd.isna(pos_str):
                return True  # keep rows with no tryptic data
            try:
                positions = ast.literal_eval(str(pos_str))
                if isinstance(positions, (int, float)):
                    positions = [positions]
                return any(int(p) > 2 for p in positions)
            except (ValueError, SyntaxError):
                return True
        filtered_df = filtered_df[filtered_df.apply(_has_non_nterm, axis=1)]

    if only_tmt_sig and 'TMT_qvalue_50pct' in filtered_df.columns:
        _q = pd.to_numeric(filtered_df['TMT_qvalue_50pct'], errors='coerce')
        filtered_df = filtered_df[_q < 0.2]

    if only_tmt_highly_sig and 'TMT_qvalue_50pct' in filtered_df.columns:
        _q = pd.to_numeric(filtered_df['TMT_qvalue_50pct'], errors='coerce')
        filtered_df = filtered_df[_q < 0.05]

    if only_tmt_sig_strict and 'TMT_qvalue_0pct' in filtered_df.columns:
        _q = pd.to_numeric(filtered_df['TMT_qvalue_0pct'], errors='coerce')
        filtered_df = filtered_df[_q < 0.2]

    if only_tmt_highly_sig_strict and 'TMT_qvalue_0pct' in filtered_df.columns:
        _q = pd.to_numeric(filtered_df['TMT_qvalue_0pct'], errors='coerce')
        filtered_df = filtered_df[_q < 0.05]

    if only_rna_sig and 'ROSMAP_Significant' in filtered_df.columns:
        filtered_df = filtered_df[filtered_df['ROSMAP_Significant']]

    if only_rna_highly_sig and 'ROSMAP_Highly_Significant' in filtered_df.columns:
        filtered_df = filtered_df[filtered_df['ROSMAP_Highly_Significant']]


    # ── Command-center header with live metrics ──
    total = len(filtered_df)
    genes = filtered_df['Parent_Gene'].nunique() if 'Parent_Gene' in filtered_df.columns else 0
    tmt_sig = int(filtered_df['TMT_Significant'].sum()) if 'TMT_Significant' in filtered_df.columns else 0
    rna_sig = int(filtered_df['ROSMAP_Significant'].sum()) if 'ROSMAP_Significant' in filtered_df.columns else 0

    swiss_count = 0
    noncan_count = 0
    db_col = 'Database_All' if 'Database_All' in filtered_df.columns else ('Database' if 'Database' in filtered_df.columns else None)
    if db_col:
        db_vals = filtered_df[db_col].fillna('').str.lower()
        swiss_count = int(db_vals.str.contains('swiss', na=False).sum())
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

    if len(filtered_df) == 0:
        st.warning("No microproteins match your current filters. Try broadening the criteria.")
        return

    # ── Prepare display table ──
    display_cols = ['sequence']

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
        'Start_Codon': 'Start Codon',
        'TMT_log2fc_50pct': 'TMT log2FC (50%)',
        'TMT_qvalue_50pct': 'TMT q-val (50%)',
        'TMT_log2fc_0pct':  'TMT log2FC (0%)',
        'TMT_qvalue_0pct':  'TMT q-val (0%)',
        'TMT_rate_control': 'MS Detect Control',
        'TMT_rate_ad':      'MS Detect AD',
        'ROSMAP_log2FC': 'RNA log2FC',
        'ROSMAP_padj': 'RNA padj',
        'Corr_MainORF_NonAD': 'Corr MainORF NonAD',
        'Corr_MainORF_AD': 'Corr MainORF AD',
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
    }

    for orig, display_name in column_mapping.items():
        if orig in filtered_df.columns:
            filtered_df[display_name] = filtered_df[orig]
            display_cols.append(display_name)

    # Mirror plot quality (use pre-computed column)
    if '_Spectra_Quality' in filtered_df.columns:
        filtered_df['Spectra Quality'] = filtered_df['_Spectra_Quality'].map(
            lambda q: QUALITY_EMOJI.get(q, '\u2014')
        )
        display_cols.append('Spectra Quality')

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
    #    smORF/RP3 fields are non-applicable (em-dash). ──
    _swiss_mask = None
    if 'Database_All' in filtered_df.columns:
        _swiss_mask = filtered_df['Database_All'].str.contains('swiss', case=False, na=False).values
        # These columns are *not applicable* to Swiss-Prot reviewed entries.
        _na_cols = ['ShortStop', 'Annotation Method',
                    'RP3 Default', 'RP3 MM+Amb', 'RP3 Amb', 'RP3 MM', 'RiboCode']
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

    # ── Data table + detail panel (wrapped in a fragment so row clicks
    #    don't restart the entire script — only this block reruns) ──
    _render_table_and_detail(
        filtered_df, display_df,
        mirror_index, expression_index, seq_to_coords,
    )

    # ── Downloads section ──
    _render_downloads_section(filtered_df, display_df)


@st.fragment
def _render_table_and_detail(filtered_df, display_df, mirror_index, expression_index, seq_to_coords):
    """Table + row-click detail panel. Wrapped in @st.fragment so a row click
    only reruns this block instead of restarting the whole script."""
    selection = st.dataframe(
        display_df,
        width='stretch',
        hide_index=True,
        on_select="rerun",
        selection_mode="single-row",
        column_config={
            'DB': st.column_config.TextColumn(
                'DB', help='🔵 Swiss-Prot (reviewed)  ·  🟠 Unreviewed', width='small'
            ),
            'UCSC': st.column_config.LinkColumn(
                'UCSC', help='View in UCSC Genome Browser', display_text='View'
            ),
            'sequence': st.column_config.TextColumn('Sequence', max_chars=50),
            'Parent Gene': st.column_config.TextColumn('Parent Gene', max_chars=15),
            'General smORF Type': st.column_config.TextColumn('General smORF Type', max_chars=20),
            'smORF Subtype': st.column_config.TextColumn('smORF Subtype', max_chars=15),
            'ShortStop Score': st.column_config.TextColumn('ShortStop Score', help='ShortStop ML confidence score (0-1)', max_chars=10),
            'PhyloCSF Score': st.column_config.NumberColumn('PhyloCSF Score', help='PhyloCSF evolutionary conservation score', format='%.2f'),
            'Unique Spectral Counts (DDA)': st.column_config.NumberColumn('Unique Spectral Counts (DDA)', help='Number of unique mass spectrometry spectral counts', format='%d'),
            'Protein Length': st.column_config.NumberColumn('Protein Length', help='Length of protein in amino acids', format='%d'),
            'Annotation Method': st.column_config.TextColumn('Annotation Method', help='Method used for annotation (MS = Mass Spectrometry, etc.)', max_chars=10),
            'TMT log2FC (50%)': st.column_config.NumberColumn('TMT log2FC (50%)', help='TMT log2 fold-change (AD vs Control) — 50% missing threshold', format='%.3f'),
            'TMT q-val (50%)':  st.column_config.NumberColumn('TMT q-val (50%)',  help='TMT q-value (BH-adjusted) — 50% missing threshold', format='%.4f'),
            'TMT log2FC (0%)':  st.column_config.NumberColumn('TMT log2FC (0%)',  help='TMT log2 fold-change (AD vs Control) — 0% missing threshold (stringent)', format='%.3f'),
            'TMT q-val (0%)':   st.column_config.NumberColumn('TMT q-val (0%)',   help='TMT q-value (BH-adjusted) — 0% missing threshold (stringent)', format='%.4f'),
            'MS Detect Control': st.column_config.NumberColumn('MS Detect Control', help='MS detection rate in Control donors', format='%.3f'),
            'MS Detect AD':      st.column_config.NumberColumn('MS Detect AD',      help='MS detection rate in AD donors', format='%.3f'),
            'RNA log2FC': st.column_config.NumberColumn('RNA log2FC', help='ROSMAP short-read RNA-seq log2 fold-change', format='%.3f'),
            'RNA padj': st.column_config.NumberColumn('RNA padj', help='ROSMAP short-read RNA-seq adjusted p-value', format='%.4f'),
            'Corr MainORF NonAD': st.column_config.NumberColumn('Corr MainORF NonAD', help='Correlation between Main ORF and smORF transcripts (non-AD, ROSMAP)', format='%.3f'),
            'Corr MainORF AD': st.column_config.NumberColumn('Corr MainORF AD', help='Correlation between Main ORF and smORF transcripts (AD, ROSMAP)', format='%.3f'),
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
            'Spectra Quality': st.column_config.TextColumn('Spectra Quality', help='Best Prosit mirror plot spectral match quality', max_chars=20),
        }
    )

    # ── Row click: Detail panel docked below table ──
    selected_rows = selection.selection.rows if selection else []
    if selected_rows:
        # CSS class to square the table bottom corners
        st.markdown('<style>.table-has-detail .stDataFrame { border-radius: 12px 12px 0 0 !important; border-bottom: none !important; }</style>', unsafe_allow_html=True)
        st.markdown('<div class="table-has-detail" style="display:none;"></div>', unsafe_allow_html=True)

        row_idx = selected_rows[0]
        _show_protein_id_card(filtered_df, display_df, row_idx)

        # ── Tabbed panel: PROSIT | Co-expression | Cartoon ──
        # RNA co-expression and smORF cartoon are only available for Unreviewed MPs
        orig_idx = selected_rows[0]
        _selected_row = filtered_df.iloc[orig_idx] if orig_idx < len(filtered_df) else None
        _db_val = str(_selected_row.get('Database_All', _selected_row.get('Database', '')) or '').lower() if _selected_row is not None else ''
        _is_swiss = 'swiss' in _db_val

        st.markdown('<div class="glass-card-section" style="margin-top:0.5rem;">', unsafe_allow_html=True)
        if _is_swiss:
            tab_prosit, = st.tabs(["PROSIT Mirror Plots"])
            with tab_prosit:
                if mirror_index:
                    _show_mirror_plots(display_df, filtered_df, orig_idx, mirror_index)
                else:
                    st.info("No PROSIT mirror plots available.")
        else:
            tab_prosit, tab_expr, tab_cartoon = st.tabs(["PROSIT Mirror Plots", "RNA Co-expression", "smORF Cartoon"])
            with tab_prosit:
                if mirror_index:
                    _show_mirror_plots(display_df, filtered_df, orig_idx, mirror_index)
                else:
                    st.info("No PROSIT mirror plots available.")
            with tab_expr:
                _show_expression_profile(filtered_df, display_df, orig_idx, expression_index, seq_to_coords)
            with tab_cartoon:
                _show_cartoon_figure(filtered_df, display_df, orig_idx, seq_to_coords)
        st.markdown('</div>', unsafe_allow_html=True)
    else:
        st.markdown('<div class="detail-panel"><div class="detail-prompt">Click a row to view protein details and MS/MS spectra</div></div>', unsafe_allow_html=True)


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
                "Unreviewed_Brain_Microproteins_CDS_absent_from_UniProt.bed",
                "BED",
                "CDS genomic coordinates (hg38) — absent from UniProt",
                "text/plain",
            ),
            (
                "Unreviewed_Brain_Microproteins_absent_from_UniProt.gtf",
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


def _show_protein_id_card(filtered_df, display_df, row_idx):
    """Show comprehensive protein ID card for selected row with glassmorphism styling."""
    orig_idx = display_df.index[row_idx]
    row = filtered_df.loc[orig_idx]

    gene = _fmt(row.get('Parent_Gene'))
    seq = str(row.get('sequence', ''))
    label = gene if gene != '—' else seq[:30]

    # Database badge
    db_val = str(row.get('Database', '')).lower()
    if 'swiss' in db_val:
        badge = '<span class="badge-swiss">Reviewed</span>'
    else:
        badge = '<span class="badge-unreviewed">Unreviewed</span>'

    # ── Card header (docked panel) ──
    st.markdown(f"""
    <div class="detail-panel">
        <div class="detail-header">
            <span style="font-size:1.3rem;"></span>
            <span class="detail-label">{label}</span>
            <span class="detail-badge">{badge}</span>
        </div>
    </div>
    """, unsafe_allow_html=True)

    # ── Core Identity ──
    smorf_raw = _fmt(row.get('smORF_Class'))
    smorf = SMORF_DISPLAY_LABEL.get(smorf_raw, smorf_raw) if smorf_raw != '—' else '—'
    smorf_group = _fmt(row.get('smORF_Group'))
    plen = _fmt(row.get('Protein_Length'), '%d')
    codon = _fmt(row.get('Start_Codon'))
    annot = _fmt(row.get('Annotation_Status'))
    ms_ev = _fmt(row.get('MS_Evidence_Type'))
    dda = _fmt(row.get('DDA_Grade'))
    ss_label = _fmt(row.get('ShortStop_Label'))
    ss_score = _fmt(row.get('ShortStop_Score'), '%.3f')
    phylo = _fmt(row.get('PhyloCSF_Score'), '%.2f')
    spec = _fmt(row.get('Unique_Spectral_Counts'), '%d')
    coords = _fmt(row.get('smORF_Coordinates'))

    st.markdown(f"""
    <div class="glass-card-section">
        <div class="id-section-header">Core Identity &amp; Annotation</div>
        <div style="display:grid; grid-template-columns: repeat(4, 1fr); gap:0.4rem 1.5rem;">
            <div><div class="id-field-label">Parent Gene</div><div class="id-field-value">{gene}</div></div>
            <div><div class="id-field-label">General smORF Type</div><div class="id-field-value">{smorf_group}</div></div>
            <div><div class="id-field-label">smORF Subtype</div><div class="id-field-value">{smorf}</div></div>
            <div><div class="id-field-label">Protein Length</div><div class="id-field-value">{plen} aa</div></div>
            <div><div class="id-field-label">Start Codon</div><div class="id-field-value">{codon}</div></div>
            <div><div class="id-field-label">Annotation Method</div><div class="id-field-value">{annot}</div></div>
            <div><div class="id-field-label">MS Evidence</div><div class="id-field-value">{ms_ev}</div></div>
            <div><div class="id-field-label">DDA Grade</div><div class="id-field-value">{dda}</div></div>
            <div><div class="id-field-label">Coordinates</div><div class="id-field-value" style="font-size:0.82rem;">{coords}</div></div>
        </div>
    </div>
    """, unsafe_allow_html=True)

    # ── Scores panel ──
    st.markdown(f"""
    <div class="glass-card-section">
        <div class="id-section-header">Confidence Scores</div>
        <div style="display:grid; grid-template-columns: repeat(4, 1fr); gap:0.4rem 1.5rem;">
            <div><div class="id-field-label">ShortStop Label</div><div class="id-field-value">{ss_label}</div></div>
            <div><div class="id-field-label">ShortStop Score</div><div class="id-field-value">{ss_score}</div></div>
            <div><div class="id-field-label">PhyloCSF Score</div><div class="id-field-value">{phylo}</div></div>
            <div><div class="id-field-label">Unique Spectral Counts</div><div class="id-field-value">{spec}</div></div>
        </div>
    </div>
    """, unsafe_allow_html=True)

    # ── TMT Proteomics (only if data exists) ──
    tmt_fc = row.get('TMT_log2fc')
    tmt_pv = row.get('TMT_pvalue')
    tmt_qv = row.get('TMT_qvalue')
    _has_tmt = any(_not_na(v) for v in [tmt_fc, tmt_pv, tmt_qv])

    if _has_tmt:
        fc_cls = _sig_class_fc(tmt_fc, tmt_qv)
        q_cls = _sig_class_pval(tmt_qv, 0.2)
        tmt_check = ' \u2705' if _not_na(tmt_qv) and tmt_qv < 0.2 else ''
        st.markdown(f"""
        <div class="glass-card-section">
            <div class="id-section-header">TMT Proteomics (AD vs Control)</div>
            <div style="display:grid; grid-template-columns: repeat(3, 1fr); gap:0.4rem 1.5rem;">
                <div><div class="id-field-label">log2 Fold Change</div><div class="id-field-value {fc_cls}">{_fmt(tmt_fc, '%.4f')}</div></div>
                <div><div class="id-field-label">p-value</div><div class="id-field-value">{_fmt(tmt_pv, '%.2e')}</div></div>
                <div><div class="id-field-label">q-value (BH)</div><div class="id-field-value {q_cls}">{_fmt(tmt_qv, '%.4f')}{tmt_check}</div></div>
            </div>
        </div>
        """, unsafe_allow_html=True)

    # ── Short-Read RNA ──
    ros_fc = row.get('ROSMAP_log2FC')
    ros_pv = row.get('ROSMAP_pvalue')
    ros_pa = row.get('ROSMAP_padj')
    ros_cpm = row.get('ROSMAP_CPM')
    _has_rna = any(_not_na(v) for v in [ros_fc, ros_pv, ros_pa])

    if _has_rna:
        fc_cls = _sig_class_fc(ros_fc, ros_pa)
        pa_cls = _sig_class_pval(ros_pa, 0.2)
        rna_check = ' \u2705' if _not_na(ros_pa) and ros_pa < 0.2 else ''
        st.markdown(f"""
        <div class="glass-card-section">
            <div class="id-section-header">Short-Read RNA-seq (ROSMAP)</div>
            <div style="display:grid; grid-template-columns: repeat(4, 1fr); gap:0.4rem 1.5rem;">
                <div><div class="id-field-label">log2 Fold Change</div><div class="id-field-value {fc_cls}">{_fmt(ros_fc, '%.4f')}</div></div>
                <div><div class="id-field-label">p-value</div><div class="id-field-value">{_fmt(ros_pv, '%.2e')}</div></div>
                <div><div class="id-field-label">padj</div><div class="id-field-value {pa_cls}">{_fmt(ros_pa, '%.4f')}{rna_check}</div></div>
                <div><div class="id-field-label">CPM</div><div class="id-field-value">{_fmt(ros_cpm, '%.2f')}</div></div>
            </div>
        </div>
        """, unsafe_allow_html=True)

    # ── Long-Read RNA ──
    nano_fc = row.get('Nanopore_log2FC')
    nano_pv = row.get('Nanopore_pvalue')
    nano_pa = row.get('Nanopore_padj')
    _has_nano = any(_not_na(v) for v in [nano_fc, nano_pv, nano_pa])

    if _has_nano:
        fc_cls = _sig_class_fc(nano_fc, nano_pa)
        pa_cls = _sig_class_pval(nano_pa, 0.2)
        nano_check = ' \u2705' if _not_na(nano_pa) and nano_pa < 0.2 else ''
        st.markdown(f"""
        <div class="glass-card-section">
            <div class="id-section-header">Long-Read RNA-seq (Nanopore)</div>
            <div style="display:grid; grid-template-columns: repeat(3, 1fr); gap:0.4rem 1.5rem;">
                <div><div class="id-field-label">log2 Fold Change</div><div class="id-field-value {fc_cls}">{_fmt(nano_fc, '%.4f')}</div></div>
                <div><div class="id-field-label">p-value</div><div class="id-field-value">{_fmt(nano_pv, '%.2e')}</div></div>
                <div><div class="id-field-label">padj</div><div class="id-field-value {pa_cls}">{_fmt(nano_pa, '%.4f')}{nano_check}</div></div>
            </div>
        </div>
        """, unsafe_allow_html=True)

    # ── RP3 / Ribo-Seq ──
    rp3_def = row.get('RP3_Default')
    rp3_mma = row.get('RP3_MM_Amb')
    rp3_amb = row.get('RP3_Amb')
    rp3_mm  = row.get('RP3_MM')
    ribo    = row.get('RiboCode')
    _has_rp3 = any(_not_na(v) for v in [rp3_def, rp3_mma, rp3_amb, rp3_mm, ribo])

    if _has_rp3:
        st.markdown(f"""
        <div class="glass-card-section">
            <div class="id-section-header">Ribosome Profiling (RP3)</div>
            <div style="display:grid; grid-template-columns: repeat(5, 1fr); gap:0.4rem 1.5rem;">
                <div><div class="id-field-label">RP3 Default</div><div class="id-field-value">{_fmt(rp3_def)}</div></div>
                <div><div class="id-field-label">RP3 MM+Amb</div><div class="id-field-value">{_fmt(rp3_mma)}</div></div>
                <div><div class="id-field-label">RP3 Amb</div><div class="id-field-value">{_fmt(rp3_amb)}</div></div>
                <div><div class="id-field-label">RP3 MM</div><div class="id-field-value">{_fmt(rp3_mm)}</div></div>
                <div><div class="id-field-label">RiboCode</div><div class="id-field-value">{_fmt(ribo)}</div></div>
            </div>
        </div>
        """, unsafe_allow_html=True)

    # ── scRNA Enrichment ──
    if INCLUDE_SCRNA:
        scrna_fc = row.get('scRNA_logFC')
        if _not_na(scrna_fc):
            scrna_pa = row.get('scRNA_padj')
            ct = _fmt(row.get('scRNA_celltype'))
            ct_gen = _fmt(row.get('scRNA_cell_type_general'))
            ct_display = f"{ct}" if ct_gen == '—' else f"{ct} ({ct_gen})"
            fc_cls = _sig_class_fc(scrna_fc, scrna_pa)
            st.markdown(f"""
        <div class="glass-card-section">
            <div class="id-section-header">Single-Cell RNA Enrichment</div>
            <div style="display:grid; grid-template-columns: repeat(3, 1fr); gap:0.4rem 1.5rem;">
                <div><div class="id-field-label">logFC</div><div class="id-field-value {fc_cls}">{_fmt(scrna_fc, '%.4f')}</div></div>
                <div><div class="id-field-label">padj (global)</div><div class="id-field-value">{_fmt(scrna_pa, '%.4f')}</div></div>
                <div><div class="id-field-label">Cell Type</div><div class="id-field-value">{ct_display}</div></div>
            </div>
        </div>
        """, unsafe_allow_html=True)

    # ── UCSC link ──
    ucsc_link = row.get('UCSC_Link')
    if ucsc_link and pd.notna(ucsc_link):
        link = create_ucsc_link({'CLICK_UCSC': ucsc_link}, CUSTOM_UCSC_SESSION)
        if link:
            st.markdown(f"[View in UCSC Genome Browser]({link})")

    # ── Full sequence ──
    with st.expander("Full amino acid sequence"):
        st.code(seq, language=None)


def _show_mirror_plots(display_df, filtered_df, row_idx, mirror_index):
    """Show Prosit mirror plot spectra for selected row (column-friendly)."""
    row = display_df.iloc[row_idx]
    orig_idx = display_df.index[row_idx]
    # Always use raw Tryptic_Peptides from filtered_df — the display_df version
    # has already been formatted as "PEP1 · PEP2" and can't be parsed by ast.literal_eval
    tryptic_str = (
        filtered_df.loc[orig_idx, 'Tryptic_Peptides']
        if 'Tryptic_Peptides' in filtered_df.columns
        else None
    )

    gene = row.get('Parent Gene', '')
    gene_suffix = f" \u2014 {gene}" if gene and pd.notna(gene) else ''
    if gene_suffix:
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
                badge_parts.append(f"\u2014 `{pep}`")
        st.markdown(" \u00b7 ".join(badge_parts))

    all_plots = get_matching_mirror_plots(tryptic_str, mirror_index)
    if not all_plots:
        st.info("No Prosit mirror plot spectra available for this microprotein\u2019s tryptic peptides.")
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
            key="spectra_peptide"
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
            key="spectra_plot"
        )
    else:
        plot_idx = 0

    plot_info = peptide_plots[plot_idx]
    st.markdown(
        f"**{QUALITY_EMOJI[plot_info['quality']]}** \u00b7 "
        f"Peptide: `{plot_info['peptide']}` \u00b7 "
        f"Charge: {plot_info['charge']}+ \u00b7 "
        f"Scan: {plot_info['scan']}"
    )
    st.image(plot_info['filepath'], width='stretch')


def _show_expression_profile(filtered_df, display_df, row_idx, expression_index, seq_to_coords):
    """Show RNA co-expression profile for selected row."""
    if not expression_index:
        st.info("Expression profile index not available.")
        return

    orig_idx = display_df.index[row_idx]
    row = filtered_df.loc[orig_idx]

    # Prefer smORF_Coordinates column (already in unified_df); fall back to seq lookup
    coords = row.get('smORF_Coordinates')
    if not _not_na(coords):
        seq = row.get('sequence')
        if seq and _not_na(seq):
            coords = seq_to_coords.get(str(seq).strip())
    if not _not_na(coords):
        st.info("Genomic coordinates not found for this smORF.")
        return

    coords = str(coords).strip()
    profiles = expression_index.get(coords, [])
    if not profiles:
        st.info("No expression profile available for this smORF.")
        return

    # If both coupled and non_coupled exist, use a radio toggle
    if len(profiles) > 1:
        coupling_options = {p['coupling']: i for i, p in enumerate(profiles)}
        labels = ['Coupled' if c == 'coupled' else 'Non-coupled' for c in coupling_options]
        chosen = st.radio(
            "Co-expression type:",
            options=list(coupling_options.keys()),
            format_func=lambda c: 'Coupled' if c == 'coupled' else 'Non-coupled',
            horizontal=True,
            key=f"expr_profile_{row_idx}"
        )
        selected_idx = coupling_options[chosen]
    else:
        selected_idx = 0

    profile = profiles[selected_idx]
    coupling_label = "Coupled" if profile['coupling'] == 'coupled' else "Non-coupled"
    st.caption(f"{coupling_label} · {profile['gene']} · {coords}")
    fp = profile['filepath']
    # Use st.image for PNGs (supports expand); fall back to iframe for PDFs
    if fp.endswith('.png'):
        st.image(fp, use_container_width=True)
    else:
        _display_pdf_inline(fp, height=420)


@st.cache_data(show_spinner=False)
def _build_cartoon_remote_set():
    """Set of cartoon filenames available remotely on HF (fallback only)."""
    files = set()
    for rel in _list_remote_files():
        if rel.startswith("smorf_cartoon_figures/") and rel.endswith(".png"):
            files.add(rel.split("/", 1)[1])
    return files


def _show_cartoon_figure(filtered_df, display_df, row_idx, seq_to_coords):
    """Show smORF cartoon figure using coordinate → image lookup."""
    orig_idx = display_df.index[row_idx]
    row = filtered_df.loc[orig_idx]

    # Use smORF_Coordinates directly; fall back to seq_to_coords TSV lookup
    coords = row.get('smORF_Coordinates')
    if not _not_na(coords):
        seq = row.get('sequence')
        if seq and _not_na(seq):
            coords = seq_to_coords.get(str(seq).strip())
    if not _not_na(coords):
        return
    coords = str(coords).strip()
    if ':' not in coords:
        return
    # Convert chr10:73799111-73799263 -> chr10_73799111-73799263.png
    filename = coords.replace(':', '_') + '.png'

    src = None
    if SMORF_CARTOON_DIR.exists():
        img_path = SMORF_CARTOON_DIR / filename
        if img_path.exists():
            src = str(img_path)
    if src is None:
        # Remote fallback
        if filename in _build_cartoon_remote_set():
            src = _remote_url(f"smorf_cartoon_figures/{filename}")

    if src is not None:
        gene = row.get('Parent_Gene', '')
        if gene and pd.notna(gene):
            st.caption(coords)
        st.image(src, use_container_width=True)
    else:
        st.info("No cartoon figure available for this smORF.")


if __name__ == "__main__":
    main()
