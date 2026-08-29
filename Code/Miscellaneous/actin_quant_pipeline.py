"""
actin_quant_pipeline.py
=======================
Generalized pipeline for quantifying subcellular F-actin distribution
from multi-channel fluorescence microscopy images (ND2 or TIFF).

Usage
-----
1. Install dependencies:
       pip install numpy scikit-image scipy nd2 matplotlib tifffile

2. Run, supplying the image folder (and optionally an output folder):
       python actin_quant_pipeline.py --image-folder /path/to/images \
           --output-folder /path/to/output

   All paths and key parameters can also be overridden on the command line.
   Run `python actin_quant_pipeline.py --help` for the full list.

Output
------
  actin_data.csv        — per-bin intensity profiles + summary metrics
  diagnostic_plots/     — overlay PNG for every image (for QC review)

Measurement outputs per cell
-----------------------------
  Bin 0 … Bin N-1      : mean actin intensity projected along the cell minor axis
                          (Bin 0 / Bin N-1 = cortical edges; central bins = nuclear region)
  AverageNuclearActin   : mean actin intensity within the nucleus mask
  AverageCellActin      : mean actin intensity within the full cell mask
  NucToCellRatio        : AverageNuclearActin / AverageCellActin  (derived in stats)
"""

# ──────────────────────────────────────────────────────────────────────────────
# DEFAULT CONFIGURATION — overridable via command-line arguments
# ──────────────────────────────────────────────────────────────────────────────

# Channel indices (0-based) in the multi-channel image
NUCLEUS_CH = 0   # nuclear stain (e.g. DAPI)
ACTIN_CH   = 2   # F-actin label (e.g. phalloidin-Cy5)
# Set to None if you have no marker channel
MARKER_CH  = 1   # optional transfection/expression marker (e.g. FITC/GFP)

# Number of spatial bins along the cell minor axis
NUM_BINS = 10

# Cell selection mode per image group.
# Map filename prefixes (lowercase) to one of:
#   "all"     — quantify every valid nucleus-matched cell in the image
#   "central" — select a single cell: most central + largest nucleus
#
# Override on the CLI with --group-mode prefix=mode (repeatable), e.g.
#   --group-mode control=all --group-mode treatment=central
# Leave empty ({}) to use "all" for every image.
GROUP_SELECTION_MODE = {}

# Nuclear segmentation parameters
NUC_MIN_AREA        = 300     # minimum nucleus area in pixels
NUC_MAX_AREA_FRAC   = 0.15   # maximum nucleus area as fraction of image area
NUC_MIN_CIRCULARITY = 0.35   # 0 = line, 1 = perfect circle
NUC_DIST_SIGMA      = 25     # Gaussian sigma on distance map (larger = fewer splits)
NUC_PEAK_DIST       = 80     # minimum separation between nucleus seed peaks (pixels)

# Cell segmentation parameters
CELL_MIN_AREA          = 500   # minimum cell area in pixels
CELL_FG_THRESH_FRAC    = 0.20  # foreground threshold as fraction of Otsu value
                                # (lower = include more dim peripheral actin)

# Nucleus-to-cell area ratio limits (sanity filter for bad matches)
NUC_CELL_RATIO_MIN = 0.02
NUC_CELL_RATIO_MAX = 0.70

# ──────────────────────────────────────────────────────────────────────────────
# Dependencies
# ──────────────────────────────────────────────────────────────────────────────

import argparse
import os
import sys
import math
import glob
import csv
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import nd2
    HAS_ND2 = True
except ImportError:
    HAS_ND2 = False

try:
    import tifffile
    HAS_TIFF = True
except ImportError:
    HAS_TIFF = False

from skimage import filters, morphology, measure, segmentation, feature
from scipy import ndimage

# ──────────────────────────────────────────────────────────────────────────────
# Image I/O
# ──────────────────────────────────────────────────────────────────────────────

def read_image(path):
    """Read an ND2 or TIFF file and return a (C, H, W) float64 array.

    Z-stacks (Z, C, H, W) are max-projected over Z before returning.
    """
    ext = os.path.splitext(path)[1].lower()
    if ext == ".nd2":
        if not HAS_ND2:
            raise ImportError("nd2 package required: pip install nd2")
        with nd2.ND2File(path) as f:
            img = f.asarray()
        img = np.squeeze(img)
    else:
        if not HAS_TIFF:
            raise ImportError("tifffile package required: pip install tifffile")
        img = tifffile.imread(path)

    if img.ndim == 2:
        img = img[np.newaxis]
    elif img.ndim == 3:
        if img.shape[0] > 10:          # (H, W, C) → (C, H, W)
            img = np.moveaxis(img, -1, 0)
    elif img.ndim == 4:                # (Z, C, H, W) → max project over Z
        img = np.max(img, axis=0)

    return img.astype(np.float64)


def list_images(folder):
    """Return sorted list of ND2/TIFF file paths in *folder*."""
    patterns = ["*.nd2", "*.ND2", "*.tif", "*.tiff", "*.TIF", "*.TIFF"]
    files = []
    for pat in patterns:
        files.extend(glob.glob(os.path.join(folder, pat)))
    return sorted(set(files))

# ──────────────────────────────────────────────────────────────────────────────
# Nuclear segmentation
# ──────────────────────────────────────────────────────────────────────────────

def segment_nuclei(dapi):
    """Segment individual nuclei from a DAPI channel image.

    Steps:
      1. Gaussian smooth → Otsu threshold → morphological open → fill holes
      2. Distance transform, smoothed to avoid over-splitting elongated nuclei
      3. Watershed seeded from distance-map peaks
      4. Merge over-split segments that share the same binary connected component
      5. Filter by area and circularity

    Returns
    -------
    labels : (H, W) int array  — labelled nucleus mask (0 = background)
    props  : list of RegionProperties
    """
    h, w = dapi.shape
    max_area = NUC_MAX_AREA_FRAC * h * w

    smoothed = ndimage.gaussian_filter(dapi, sigma=2)
    thresh   = filters.threshold_otsu(smoothed)
    binary   = smoothed > thresh
    binary   = morphology.binary_opening(binary, morphology.disk(3))
    binary   = ndimage.binary_fill_holes(binary)

    distance        = ndimage.distance_transform_edt(binary)
    distance_smooth = ndimage.gaussian_filter(distance, sigma=NUC_DIST_SIGMA)

    peaks = feature.peak_local_max(
        distance_smooth,
        min_distance=NUC_PEAK_DIST,
        labels=binary.astype(int),
        exclude_border=False,
    )
    markers = np.zeros_like(binary, dtype=int)
    for idx, (r, c) in enumerate(peaks, start=1):
        markers[r, c] = idx
    markers = ndimage.label(morphology.dilation(markers > 0, morphology.disk(3)))[0]

    if markers.max() > 1:
        labels = segmentation.watershed(-distance, markers, mask=binary)
    else:
        labels = measure.label(binary)

    # Merge watershed segments that belong to the same binary connected component
    # (prevents curved/kidney-shaped nuclei from being split into two labels)
    cc_labels = measure.label(binary)
    ws_to_cc  = {}
    for rp in measure.regionprops(labels):
        cc_vals = cc_labels[labels == rp.label]
        cc_vals = cc_vals[cc_vals > 0]
        if len(cc_vals):
            ws_to_cc[rp.label] = int(np.bincount(cc_vals).argmax())

    cc_to_ws = defaultdict(list)
    for ws_lbl, cc_lbl in ws_to_cc.items():
        cc_to_ws[cc_lbl].append(ws_lbl)

    merged = np.zeros_like(labels)
    for new_lbl, (_, ws_lbls) in enumerate(cc_to_ws.items(), start=1):
        for wl in ws_lbls:
            merged[labels == wl] = new_lbl
    labels = merged

    # Filter by area and circularity
    good = []
    for rp in measure.regionprops(labels):
        if rp.area < NUC_MIN_AREA or rp.area > max_area:
            continue
        circ = (4 * math.pi * rp.area / rp.perimeter ** 2
                if rp.perimeter > 0 else 0)
        if circ < NUC_MIN_CIRCULARITY:
            continue
        good.append(rp)

    clean = np.zeros_like(labels)
    for i, rp in enumerate(good, start=1):
        clean[labels == rp.label] = i

    return clean, measure.regionprops(clean)

# ──────────────────────────────────────────────────────────────────────────────
# Cell segmentation
# ──────────────────────────────────────────────────────────────────────────────

def segment_cells(actin, nucleus_labels):
    """Segment cell bodies from the actin channel using seeded watershed.

    Foreground is defined at CELL_FG_THRESH_FRAC × Otsu to capture dim
    peripheral actin without including background.  Nucleus labels seed
    the watershed; border pixels seed the background.

    Returns
    -------
    labels : (H, W) int array
    props  : list of RegionProperties
    """
    h, w     = actin.shape
    smoothed = ndimage.gaussian_filter(actin, sigma=2)
    thresh   = filters.threshold_otsu(smoothed)

    fg_mask = smoothed > thresh * CELL_FG_THRESH_FRAC
    fg_mask = ndimage.binary_fill_holes(fg_mask)
    fg_mask = morphology.binary_closing(fg_mask, morphology.disk(5))
    fg_mask = ndimage.binary_fill_holes(fg_mask)

    markers      = nucleus_labels.copy()
    bg_val       = markers.max() + 1
    markers[~fg_mask & (markers == 0)] = bg_val

    # Seed image border as background so neighbour cells are not absorbed
    edge = 10
    border = np.zeros_like(fg_mask)
    border[:edge, :] = border[-edge:, :] = True
    border[:, :edge] = border[:, -edge:] = True
    markers[border & fg_mask & (nucleus_labels == 0)] = bg_val

    gradient  = filters.sobel(smoothed)
    ws_labels = segmentation.watershed(gradient, markers, mask=fg_mask)
    ws_labels[ws_labels == bg_val] = 0

    clean = np.zeros_like(ws_labels)
    for i, rp in enumerate(measure.regionprops(ws_labels), start=1):
        if rp.area >= CELL_MIN_AREA:
            clean[ws_labels == rp.label] = i

    return clean, measure.regionprops(clean)

# ──────────────────────────────────────────────────────────────────────────────
# Nucleus–cell matching
# ──────────────────────────────────────────────────────────────────────────────

def match_nuclei_to_cells(nuc_labels, cell_labels):
    """Map each cell label to the nucleus label with greatest overlap inside it.

    Returns dict: {cell_label: nucleus_label or None}
    """
    matches = {}
    for cp in measure.regionprops(cell_labels):
        cell_mask  = cell_labels == cp.label
        nuc_inside = np.unique(nuc_labels[cell_mask])
        nuc_inside = nuc_inside[nuc_inside > 0]
        if len(nuc_inside) == 0:
            matches[cp.label] = None
        elif len(nuc_inside) == 1:
            matches[cp.label] = int(nuc_inside[0])
        else:
            best = max(nuc_inside,
                       key=lambda n: np.sum(cell_mask & (nuc_labels == n)))
            matches[cp.label] = int(best)
    return matches


def _valid_candidates(nuc_props, cell_labels, cell_props, matches):
    """Return list of (cell_prop, nuc_prop) passing the area ratio filter."""
    nuc_to_cell = {nl: cl for cl, nl in matches.items() if nl is not None}
    valid = []
    for np_ in nuc_props:
        cl = nuc_to_cell.get(np_.label)
        if cl is None:
            continue
        cell_area = np.sum(cell_labels == cl)
        nuc_frac  = np_.area / cell_area if cell_area > 0 else 0
        if not (NUC_CELL_RATIO_MIN <= nuc_frac <= NUC_CELL_RATIO_MAX):
            continue
        cp = next((c for c in cell_props if c.label == cl), None)
        if cp is not None:
            valid.append((cp, np_))
    return valid

# ──────────────────────────────────────────────────────────────────────────────
# Cell selection strategies
# ──────────────────────────────────────────────────────────────────────────────

def select_all_cells(nuc_props, cell_labels, cell_props, matches):
    """Return every nucleus-matched cell that passes the area ratio filter."""
    return _valid_candidates(nuc_props, cell_labels, cell_props, matches)


def select_central_cell(nuc_props, cell_labels, cell_props, matches, h, w):
    """Return the single most central, largest-nucleus cell as (cell_prop, nuc_prop).

    Score = norm_distance − 0.5 × norm_nucleus_area  (lower = better)
    """
    candidates = _valid_candidates(nuc_props, cell_labels, cell_props, matches)
    if not candidates:
        return None, None

    cy_img, cx_img = h / 2.0, w / 2.0
    scored = []
    for cp, np_ in candidates:
        ncy, ncx = np_.centroid
        dist = math.hypot(ncy - cy_img, ncx - cx_img)
        scored.append((cp, np_, dist))

    max_dist = max(s[2] for s in scored) or 1.0
    max_area = max(s[1].area for s in scored) or 1.0

    best = min(scored,
               key=lambda s: s[2] / max_dist - 0.5 * s[1].area / max_area)
    return best[0], best[1]

# ──────────────────────────────────────────────────────────────────────────────
# Actin distribution profiling
# ──────────────────────────────────────────────────────────────────────────────

def projected_profile(actin, cell_mask, cell_region, nuc_region=None,
                      num_bins=NUM_BINS):
    """Project all cell pixels onto the cell minor axis and bin by position.

    If nuc_region is provided the axis is re-centred on the nucleus centroid
    and the extent is made symmetric, so the nucleus always falls in the
    central bins (bins N//2−1 and N//2).

    Empty bins at the edges (asymmetric cells) are filled by nearest-neighbour
    interpolation from the closest populated bin.

    Returns 1-D array of length num_bins, or empty array on failure.
    """
    cy, cx    = cell_region.centroid
    angle_rad = -cell_region.orientation + math.pi / 2  # perpendicular to major axis
    ax_x      = math.cos(angle_rad)
    ax_y      = math.sin(angle_rad)

    rows, cols   = np.where(cell_mask)
    intensities  = actin[rows, cols]
    proj         = (cols - cx) * ax_x + (rows - cy) * ax_y

    if proj.max() - proj.min() < 1:
        return np.array([])

    if nuc_region is not None:
        ncy, ncx  = nuc_region.centroid
        nuc_proj  = (ncx - cx) * ax_x + (ncy - cy) * ax_y
        proj      = proj - nuc_proj
        max_ext   = max(abs(proj.min()), abs(proj.max()))
        p_min, p_max = -max_ext, max_ext
    else:
        p_min, p_max = proj.min(), proj.max()

    edges     = np.linspace(p_min, p_max, num_bins + 1)
    bin_means = np.full(num_bins, np.nan)
    for i in range(num_bins):
        mask = ((proj >= edges[i]) & (proj < edges[i + 1])
                if i < num_bins - 1
                else (proj >= edges[i]) & (proj <= edges[i + 1]))
        if np.any(mask):
            bin_means[i] = float(np.mean(intensities[mask]))

    # Fill empty bins by nearest-neighbour interpolation
    for i in range(num_bins):
        if np.isnan(bin_means[i]):
            for d in range(1, num_bins):
                if i + d < num_bins and not np.isnan(bin_means[i + d]):
                    bin_means[i] = bin_means[i + d]; break
                if i - d >= 0     and not np.isnan(bin_means[i - d]):
                    bin_means[i] = bin_means[i - d]; break

    return bin_means

# ──────────────────────────────────────────────────────────────────────────────
# CSV output
# ──────────────────────────────────────────────────────────────────────────────

def write_row(csv_path, group, cell_id, bin_means, actin_nuc, actin_cell):
    """Append one cell's data (10 bin rows) to actin_data.csv."""
    header = ["DataGroup", "ImageID", "Bin", "Average",
              "AverageNuclearActin", "AverageCellActin"]
    write_header = not os.path.exists(csv_path)
    with open(csv_path, "a", newline="") as fh:
        writer = csv.writer(fh, delimiter=";")
        if write_header:
            writer.writerow(header)
        for i, val in enumerate(bin_means):
            writer.writerow([group, cell_id, i,
                             f"{val:.6f}", f"{actin_nuc:.6f}", f"{actin_cell:.6f}"])

# ──────────────────────────────────────────────────────────────────────────────
# Diagnostic plot
# ──────────────────────────────────────────────────────────────────────────────

def save_diagnostic(path, nucleus, actin, marker,
                    nuc_labels, nuc_props,
                    cell_labels, cell_props,
                    selected_cells):
    """Save a 2×3 diagnostic overlay PNG for one image."""
    h, w = actin.shape
    selected_cell_labels = {cp.label for cp, _ in selected_cells}
    selected_nuc_labels  = {np_.label for _, np_ in selected_cells}

    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle(os.path.basename(path), fontsize=12, fontweight="bold")

    axes[0, 0].imshow(nucleus, cmap="Blues");  axes[0, 0].set_title("Nuclear channel")
    axes[0, 1].imshow(actin,   cmap="Reds");   axes[0, 1].set_title("Actin channel")
    if marker is not None:
        axes[0, 2].imshow(marker, cmap="Greens"); axes[0, 2].set_title("Marker channel")
    else:
        axes[0, 2].set_visible(False)

    # Nucleus outlines
    axes[1, 0].imshow(nucleus, cmap="gray")
    for rp in nuc_props:
        col = "lime" if rp.label in selected_nuc_labels else "yellow"
        lw  = 2.5   if rp.label in selected_nuc_labels else 1.0
        for c in measure.find_contours((nuc_labels == rp.label).astype(float), 0.5):
            axes[1, 0].plot(c[:, 1], c[:, 0], color=col, linewidth=lw)
    axes[1, 0].set_title(f"Nuclei — {len(selected_nuc_labels)} selected (lime)")

    # Cell outlines
    axes[1, 1].imshow(actin, cmap="gray")
    cmap_set = plt.cm.Set2(np.linspace(0, 1, max(len(cell_props), 1)))
    for ci, cp in enumerate(cell_props):
        col = "cyan"              if cp.label in selected_cell_labels else cmap_set[ci % len(cmap_set)]
        lw  = 2.5                 if cp.label in selected_cell_labels else 1.0
        for c in measure.find_contours((cell_labels == cp.label).astype(float), 0.5):
            axes[1, 1].plot(c[:, 1], c[:, 0], color=col, linewidth=lw)
    axes[1, 1].set_title(f"Cells — {len(selected_cell_labels)} selected (cyan)")

    # Composite overlay
    composite = np.zeros((h, w, 3), dtype=float)
    composite[:, :, 0] = actin   / (actin.max()   or 1)
    composite[:, :, 2] = nucleus / (nucleus.max() or 1)
    if marker is not None:
        composite[:, :, 1] = marker / (marker.max() or 1)
    axes[1, 2].imshow(np.clip(composite, 0, 1))
    for cp, np_ in selected_cells:
        for c in measure.find_contours((cell_labels == cp.label).astype(float), 0.5):
            axes[1, 2].plot(c[:, 1], c[:, 0], "white", linewidth=2)
        for c in measure.find_contours((nuc_labels == np_.label).astype(float), 0.5):
            axes[1, 2].plot(c[:, 1], c[:, 0], "yellow", linewidth=1.5)
    axes[1, 2].set_title("Composite (white=cell, yellow=nucleus)")

    for ax in axes.flat:
        ax.axis("off")
    plt.tight_layout()
    plt.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)

# ──────────────────────────────────────────────────────────────────────────────
# Group assignment
# ──────────────────────────────────────────────────────────────────────────────

def get_group(basename):
    """Assign an image to an experimental group based on filename prefix.

    The group name drives cell selection mode (via GROUP_SELECTION_MODE).
    Override this function to implement custom naming conventions.
    """
    name = basename.lower()
    for prefix in sorted(GROUP_SELECTION_MODE.keys(), key=len, reverse=True):
        if name.startswith(prefix.lower()):
            return prefix
    return name  # default: group = full basename

# ──────────────────────────────────────────────────────────────────────────────
# Per-file processing
# ──────────────────────────────────────────────────────────────────────────────

def process_file(fpath, csv_path, diag_folder):
    """Process one image file, write results to CSV, save diagnostic plot.

    Returns a dict summarising the result for the pipeline summary table.
    """
    basename = os.path.splitext(os.path.basename(fpath))[0]
    group    = get_group(basename)
    mode     = GROUP_SELECTION_MODE.get(group, globals().get("_DEFAULT_MODE", "all"))

    print(f"  {basename}  group={group}  mode={mode}")

    img    = read_image(fpath)
    actin  = img[ACTIN_CH]
    nucleus = img[NUCLEUS_CH]
    marker  = img[MARKER_CH] if MARKER_CH is not None and MARKER_CH < img.shape[0] else None
    h, w   = actin.shape

    # Segmentation
    nuc_labels, nuc_props = segment_nuclei(nucleus)
    if not nuc_props:
        print("    → no nuclei found, skipping")
        return {"file": basename, "group": group, "status": "NO_NUCLEI", "n_cells": 0}

    cell_labels, cell_props = segment_cells(actin, nuc_labels)
    matches = match_nuclei_to_cells(nuc_labels, cell_labels)

    # Cell selection
    if mode == "central":
        best_cell, best_nuc = select_central_cell(
            nuc_props, cell_labels, cell_props, matches, h, w)
        selected = [(best_cell, best_nuc)] if best_cell is not None else []
    else:
        selected = select_all_cells(nuc_props, cell_labels, cell_props, matches)

    if not selected:
        print("    → no valid cells found, skipping")
        save_diagnostic(
            os.path.join(diag_folder, f"{basename}.png"),
            nucleus, actin, marker, nuc_labels, nuc_props,
            cell_labels, cell_props, [])
        return {"file": basename, "group": group, "status": "NO_CELLS", "n_cells": 0}

    # Profile and write
    cells_written = 0
    for ci, (cp, np_) in enumerate(selected, start=1):
        cell_mask = cell_labels == cp.label
        nuc_mask  = nuc_labels  == np_.label

        actin_cell = float(np.mean(actin[cell_mask]))
        actin_nuc  = float(np.mean(actin[nuc_mask]))
        bins       = projected_profile(actin, cell_mask, cp, nuc_region=np_)

        if len(bins) < NUM_BINS or np.any(np.isnan(bins)):
            print(f"    cell {ci}: invalid profile, skipping")
            continue

        cell_id = f"{basename}_cell{ci}" if len(selected) > 1 else basename
        write_row(csv_path, group, cell_id, bins, actin_nuc, actin_cell)
        cells_written += 1
        print(f"    cell {ci}: nuc={actin_nuc:.0f}  cell={actin_cell:.0f}  "
              f"ratio={actin_nuc/actin_cell:.3f}")

    save_diagnostic(
        os.path.join(diag_folder, f"{basename}.png"),
        nucleus, actin, marker, nuc_labels, nuc_props,
        cell_labels, cell_props, selected)

    status = "OK" if cells_written > 0 else "NO_PROFILE"
    print(f"    → {cells_written}/{len(selected)} cells written  [{status}]")
    return {"file": basename, "group": group, "status": status, "n_cells": cells_written}

# ──────────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────────

def _parse_group_mode(items):
    """Parse repeated --group-mode prefix=mode arguments into a dict."""
    out = {}
    for it in items or []:
        if "=" not in it:
            raise argparse.ArgumentTypeError(
                f"--group-mode expects prefix=mode, got: {it!r}")
        prefix, mode = it.split("=", 1)
        mode = mode.strip().lower()
        if mode not in ("all", "central"):
            raise argparse.ArgumentTypeError(
                f"--group-mode mode must be 'all' or 'central', got: {mode!r}")
        out[prefix.strip().lower()] = mode
    return out


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Quantify subcellular F-actin distribution from "
                    "multi-channel ND2/TIFF microscopy images.")
    p.add_argument("--image-folder", "-i", required=True,
                   help="Folder containing ND2 or TIFF input images.")
    p.add_argument("--output-folder", "-o", default=None,
                   help="Folder to write outputs into. Default: "
                        "<image-folder>/actin_quant_output")
    p.add_argument("--csv-name", default="actin_data.csv",
                   help="Name of the per-cell CSV file written into the "
                        "output folder (default: actin_data.csv).")
    p.add_argument("--diag-subfolder", default="diagnostic_plots",
                   help="Sub-folder (inside --output-folder) for diagnostic "
                        "PNG overlays (default: diagnostic_plots).")
    p.add_argument("--nucleus-ch", type=int, default=NUCLEUS_CH,
                   help=f"Nuclear channel index, 0-based (default: {NUCLEUS_CH}).")
    p.add_argument("--actin-ch", type=int, default=ACTIN_CH,
                   help=f"Actin channel index, 0-based (default: {ACTIN_CH}).")
    p.add_argument("--marker-ch", type=int, default=MARKER_CH,
                   help=f"Marker channel index, 0-based; pass -1 for none "
                        f"(default: {MARKER_CH}).")
    p.add_argument("--num-bins", type=int, default=NUM_BINS,
                   help=f"Number of spatial bins (default: {NUM_BINS}).")
    p.add_argument("--group-mode", action="append", default=None,
                   metavar="PREFIX=MODE",
                   help="Map a filename prefix to a cell-selection mode "
                        "('all' or 'central'). May be passed multiple times. "
                        "If omitted, all images use 'all'.")
    p.add_argument("--default-mode", choices=("all", "central"), default="all",
                   help="Cell-selection mode for images whose prefix is not "
                        "matched by any --group-mode (default: all).")
    p.add_argument("--overwrite-csv", action="store_true",
                   help="Delete an existing output CSV before writing.")
    return p.parse_args(argv)


def main(argv=None):
    global NUCLEUS_CH, ACTIN_CH, MARKER_CH, NUM_BINS, GROUP_SELECTION_MODE
    global _DEFAULT_MODE

    args = parse_args(argv)

    image_folder = os.path.abspath(os.path.expanduser(args.image_folder))
    if not os.path.isdir(image_folder):
        print(f"Image folder not found: {image_folder}")
        sys.exit(1)

    output_folder = (os.path.abspath(os.path.expanduser(args.output_folder))
                     if args.output_folder
                     else os.path.join(image_folder, "actin_quant_output"))
    diag_folder = os.path.join(output_folder, args.diag_subfolder)
    csv_path    = os.path.join(output_folder, args.csv_name)

    NUCLEUS_CH = args.nucleus_ch
    ACTIN_CH   = args.actin_ch
    MARKER_CH  = None if args.marker_ch is not None and args.marker_ch < 0 else args.marker_ch
    NUM_BINS   = args.num_bins
    GROUP_SELECTION_MODE = _parse_group_mode(args.group_mode)
    _DEFAULT_MODE = args.default_mode

    image_files = list_images(image_folder)
    if not image_files:
        print(f"No images found in: {image_folder}")
        sys.exit(1)

    os.makedirs(output_folder, exist_ok=True)
    os.makedirs(diag_folder, exist_ok=True)

    if args.overwrite_csv and os.path.exists(csv_path):
        os.remove(csv_path)

    print(f"Found {len(image_files)} image(s) in {image_folder}\n")
    results = []
    for i, fpath in enumerate(image_files, start=1):
        print(f"[{i}/{len(image_files)}] {os.path.basename(fpath)}")
        try:
            r = process_file(fpath, csv_path, diag_folder)
        except Exception as e:
            print(f"    ERROR: {e}")
            r = {"file": os.path.splitext(os.path.basename(fpath))[0],
                 "group": "?", "status": "ERROR", "n_cells": 0}
        results.append(r)

    # Summary
    print(f"\n{'='*70}")
    print("PIPELINE SUMMARY")
    print(f"{'='*70}")
    total_cells = sum(r["n_cells"] for r in results)
    n_ok        = sum(1 for r in results if r["status"] == "OK")
    print(f"  Images processed : {len(results)}")
    print(f"  Images OK        : {n_ok}")
    print(f"  Total cells      : {total_cells}")
    print(f"  Output CSV       : {csv_path}")
    print(f"  Diagnostic plots : {diag_folder}/")
    print()

    groups = sorted({r["group"] for r in results})
    for g in groups:
        g_results = [r for r in results if r["group"] == g]
        g_cells   = sum(r["n_cells"] for r in g_results)
        print(f"  {g:20s}  {g_cells} cells  from  {len(g_results)} images")

    print(f"\n{'─'*70}")
    for r in results:
        flag = "  " if r["status"] == "OK" else ">>"
        print(f"  {flag} {r['file']:30s}  group={r['group']:15s}  "
              f"cells={r['n_cells']:2d}  {r['status']}")


if __name__ == "__main__":
    main()
