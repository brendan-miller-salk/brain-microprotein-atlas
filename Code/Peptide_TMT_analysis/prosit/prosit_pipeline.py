#!/usr/bin/env python3
"""
prosit_pipeline.py — Unified PROSIT TMT Spectral Prediction + SA Pipeline

Processes ALL tryptic peptides (not just one-hit-wonders) across round 1
(b1-b50) and round 2 (round2/b1-b14).  For each peptide, collects up to 5
candidate PSMs (filtered by Purity >= 0.5, ranked by Hyperscore), predicts
MS/MS with Prosit_2020_intensity_TMT via Koina, then computes spectral angle
(SA) for each candidate vs. the predicted spectrum and selects the best
representative by evidence score (Strong > Moderate > Weak > Insufficient,
SA as tiebreaker) for the mirror plot.

Output:
  - all_peptides_prosit.csv          (Phase 1: peptides + best PSM info)
  - predicted_spectra/prosit_tmt_predicted_all.msp   (Phase 1: MSP library)
  - all_peptides_with_SA.csv         (Phase 2: augmented with SA scores)
  - mirror_plots/                    (Phase 2: 3-panel annotated PNGs)

Usage:
  conda activate prosit
  python prosit_pipeline.py              # full run (Phase 1 + 2)
  python prosit_pipeline.py --test       # test mode: 20 peptides
  python prosit_pipeline.py --phase1     # Phase 1 only (predict)
  python prosit_pipeline.py --phase2     # Phase 2 only (SA; requires Phase 1 output)
  python prosit_pipeline.py --workers 32 # control parallelism
"""

import os
import sys
import re
import glob
import time
import argparse
import numpy as np
import pandas as pd
from collections import Counter
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
import warnings

# =============================================================================
# Configuration
# =============================================================================
WORKSPACE = "/project2/hassy_472/brendajm/scratch1_transfer/tmt_rosmap"
PROSIT_DIR = os.path.join(WORKSPACE, "prosit")
INPUT_CSV = os.path.join(PROSIT_DIR, "cleaned_tryptic_peptides_detailed_under_151aa.csv")

# Output paths (full mode)
OUTPUT_CSV_P1 = os.path.join(PROSIT_DIR, "all_peptides_prosit.csv")
OUTPUT_SPECTRA_DIR = os.path.join(PROSIT_DIR, "predicted_spectra")
MSP_OUTPUT = os.path.join(OUTPUT_SPECTRA_DIR, "prosit_tmt_predicted_all.msp")
OUTPUT_CSV_P2 = os.path.join(PROSIT_DIR, "all_peptides_with_SA.csv")
OUTPUT_CSV_ANNOTATED = os.path.join(PROSIT_DIR, "cleaned_tryptic_peptides_detailed_under_151aa_with_SA.csv")
PLOT_DIR = os.path.join(PROSIT_DIR, "mirror_plots")
CANDIDATE_CSV = os.path.join(PROSIT_DIR, "candidate_psms.csv")

# Output paths (test mode)
OUTPUT_CSV_P1_TEST = os.path.join(PROSIT_DIR, "all_peptides_prosit_test.csv")
MSP_OUTPUT_TEST = os.path.join(OUTPUT_SPECTRA_DIR,
                                "prosit_tmt_predicted_all_test.msp")
OUTPUT_CSV_P2_TEST = os.path.join(PROSIT_DIR, "all_peptides_with_SA_test.csv")
OUTPUT_CSV_ANNOTATED_TEST = os.path.join(PROSIT_DIR, "cleaned_tryptic_peptides_detailed_under_151aa_with_SA_test.csv")
PLOT_DIR_TEST = os.path.join(PROSIT_DIR, "mirror_plots_test")
CANDIDATE_CSV_TEST = os.path.join(PROSIT_DIR, "candidate_psms_test.csv")

# Koina / PROSIT settings
MODEL_NAME = "Prosit_2020_intensity_TMT"
SERVER_URL = "koina.wilhelmlab.org:443"
NCE = 38.0
FRAGMENTATION = "HCD"
DEFAULT_CHARGE = 2
BATCH_SIZE = 500
MAX_PEPTIDE_LEN = 30
MIN_INTENSITY = 1e-4

# Only scan the appended shortstop search results
RESULT_SUBDIR = "shortstop_proteogenomics_appended_results_cpm05/DDA"

# Fragment matching
PPM_TOLERANCE = 20.0

# TMT reporter ion exclusion window (matches MSFragger clear-mz settings)
CLEAR_MZ_LO = 125.5
CLEAR_MZ_HI = 131.5

# Representative selection
MIN_PURITY = 0.5
MAX_CANDIDATES = 5

# SA classification thresholds (degrees)
SA_HIGH = 20.0
SA_MODERATE = 45.0
SA_POOR = 70.0


# =============================================================================
# Batch enumeration — round 1 + round 2
# =============================================================================
def enumerate_batches(workspace):
    """
    Return list of (batch_label, batch_dir) tuples for all batches.
    batch_label is the relative path from workspace (e.g., 'b1' or 'round2/b3').
    batch_dir is the absolute path.
    """
    batches = []

    # Round 1: b1-b50 at workspace root
    for d in sorted(glob.glob(os.path.join(workspace, "b[0-9]*"))):
        if os.path.isdir(d):
            label = os.path.basename(d)
            batches.append((label, d))

    # Round 2: round2/b1-b14
    r2_root = os.path.join(workspace, "round2")
    if os.path.isdir(r2_root):
        for d in sorted(glob.glob(os.path.join(r2_root, "b[0-9]*"))):
            if os.path.isdir(d):
                label = f"round2/{os.path.basename(d)}"
                batches.append((label, d))

    return batches


# =============================================================================
# Phase 1 — Step 1: Read tryptic peptides
# =============================================================================
def read_tryptic_peptides(csv_path):
    import ast as _ast
    master_df = pd.read_csv(csv_path)
    print(f"  Loaded {len(master_df)} protein rows from "
          f"{os.path.basename(csv_path)}")

    # Explode bracketed peptide lists into one row per peptide
    rows = []
    for _, row in master_df.iterrows():
        protein_id = row["protein_id"]
        try:
            peps = _ast.literal_eval(row["peptide_sequence"])
        except (ValueError, SyntaxError):
            continue
        if not isinstance(peps, list):
            peps = [peps]
        for pep in peps:
            pep = str(pep).strip()
            if pep and pep.lower() != "nan":
                rows.append({"protein_id": protein_id,
                             "tryptic_peptide": pep})

    df = pd.DataFrame(rows)
    print(f"  Expanded to {len(df)} individual tryptic peptides "
          f"(from {len(master_df)} protein rows)")
    return df, master_df


# =============================================================================
# Phase 1 — Step 2: Find ALL peptides + compute global spectral counts
# =============================================================================
def find_all_peptides(workspace, peptide_set):
    """
    Scan peptide.tsv across all batches (round1 + round2).
    For each tryptic peptide found, collect spectral count per batch.

    Returns:
        all_pep_df: DataFrame with per-batch peptide observations
        global_counts: dict { peptide -> total spectral count across batches }
    """
    batches = enumerate_batches(workspace)
    n_r1 = sum(1 for b, _ in batches if "/" not in b)
    n_r2 = sum(1 for b, _ in batches if "/" in b)
    print(f"  Found {len(batches)} batch directories "
          f"(round 1: {n_r1}, round 2: {n_r2})")

    records = []

    for batch_label, batch_dir in tqdm(batches, desc="  Scanning peptide.tsv",
                                        unit="batch"):
        pep_file = os.path.join(batch_dir, RESULT_SUBDIR, "peptide.tsv")
        if not os.path.isfile(pep_file):
            continue
        try:
            df = pd.read_csv(pep_file, sep="\t", usecols=[
                "Peptide", "Charges", "Spectral Count", "Assigned Modifications"
            ], low_memory=False)
            mask = df["Peptide"].isin(peptide_set)
            matched = df.loc[mask].copy()
            if len(matched) == 0:
                continue
            matched["Batch"] = batch_label
            records.append(matched)
        except Exception as e:
            tqdm.write(f"  Warning: {batch_label}: {str(e)[:60]}")

    if not records:
        print("  WARNING: No matching peptides found in any peptide.tsv!")
        return pd.DataFrame(), {}

    all_pep = pd.concat(records, ignore_index=True)
    print(f"  Total peptide-batch observations: {len(all_pep)}")

    # Global spectral counts
    global_counts = all_pep.groupby("Peptide")["Spectral Count"].sum().to_dict()
    n_found = len(global_counts)
    n_ohw = sum(1 for c in global_counts.values() if c == 1)
    print(f"  Unique peptides found:   {n_found}")
    print(f"  One-hit-wonders (SC=1):  {n_ohw}")
    print(f"  Multi-PSM (SC>1):        {n_found - n_ohw}")

    return all_pep, global_counts


# =============================================================================
# Phase 1 — Step 3: Get candidate PSMs per peptide (Purity + Hyperscore)
# =============================================================================
def get_candidate_psms(workspace, all_pep_df):
    """
    For each peptide, scan psm.tsv across all batches it was found in.
    Filter by Purity >= MIN_PURITY, then keep up to MAX_CANDIDATES PSMs
    by Hyperscore (all sharing the same charge & modified peptide as the
    rank-1 hit so they can be compared against the same PROSIT prediction).

    Returns DataFrame with up to MAX_CANDIDATES rows per peptide.
    Includes a Candidate_Rank column (1 = best Hyperscore).
    """
    if all_pep_df.empty:
        return pd.DataFrame()

    batches_needed = all_pep_df["Batch"].unique()
    all_peptides = set(all_pep_df["Peptide"])

    psm_records = []

    for batch_label in tqdm(sorted(batches_needed),
                            desc="  Reading PSM details", unit="batch"):
        batch_dir = os.path.join(workspace, batch_label)
        psm_file = os.path.join(batch_dir, RESULT_SUBDIR, "psm.tsv")
        if not os.path.isfile(psm_file):
            continue
        try:
            for chunk in pd.read_csv(
                psm_file, sep="\t",
                usecols=["Spectrum", "Peptide", "Modified Peptide", "Charge",
                         "Assigned Modifications", "Hyperscore", "Purity"],
                chunksize=50000, low_memory=False,
            ):
                mask = chunk["Peptide"].isin(all_peptides)
                matched = chunk.loc[mask].copy()
                if len(matched) > 0:
                    matched["Batch"] = batch_label
                    psm_records.append(matched)
        except Exception as e:
            tqdm.write(f"  Warning: {batch_label}: {str(e)[:60]}")

    if not psm_records:
        print("  WARNING: No PSM-level details found!")
        return pd.DataFrame()

    psm_df = pd.concat(psm_records, ignore_index=True)
    print(f"  Total PSM rows collected: {len(psm_df)}")

    psm_df["Hyperscore"] = pd.to_numeric(psm_df["Hyperscore"], errors="coerce")
    psm_df["Purity"] = pd.to_numeric(psm_df["Purity"], errors="coerce")

    # For each peptide: prefer non-acetylated PSMs (PROSIT-compatible),
    # filter by purity, lock charge & mods to rank-1, keep top MAX_CANDIDATES
    candidates = []
    n_acetyl_rescued = 0
    n_acetyl_only = 0
    for peptide, group in psm_df.groupby("Peptide"):
        # Prefer non-N-term-acetylated PSMs (n[43] is unsupported by PROSIT)
        is_acetyl = group["Modified Peptide"].str.contains(
            r'^n\[43\]', na=False)
        non_acetyl = group[~is_acetyl]
        if len(non_acetyl) > 0:
            pool = non_acetyl
            rescued = is_acetyl.any()
            if rescued:
                n_acetyl_rescued += 1
        else:
            # All PSMs are N-term acetylated — peptide will be excluded later
            pool = group
            rescued = False
            n_acetyl_only += 1

        pure = pool[pool["Purity"] >= MIN_PURITY]
        if len(pure) == 0:
            # Fallback: keep the single best Hyperscore regardless of purity
            pure = pool.nlargest(1, "Hyperscore")
        # Rank-1 determines the target charge & modifications
        rank1 = pure.nlargest(1, "Hyperscore").iloc[0]
        target_charge = rank1["Charge"]
        target_mod = rank1["Modified Peptide"]
        matching = pure[
            (pure["Charge"] == target_charge)
            & (pure["Modified Peptide"] == target_mod)
        ]
        top = matching.nlargest(MAX_CANDIDATES, "Hyperscore").copy()
        top["Candidate_Rank"] = range(1, len(top) + 1)
        top["Acetyl_Rescued"] = rescued
        candidates.append(top)

    if n_acetyl_rescued > 0:
        print(f"  Rescued {n_acetyl_rescued} peptides by preferring "
              f"non-acetylated PSMs over n[43] rank-1")
    if n_acetyl_only > 0:
        print(f"  {n_acetyl_only} peptides have ONLY n[43] PSMs "
              f"(will be excluded downstream)")

    cand_df = pd.concat(candidates, ignore_index=True)
    n_peptides = cand_df["Peptide"].nunique()
    n_multi = (cand_df.groupby("Peptide").size() > 1).sum()
    print(f"  Candidates: {n_peptides} peptides, "
          f"{n_multi} with multiple candidates, "
          f"{len(cand_df)} total PSMs")

    return cand_df


# =============================================================================
# Phase 1 — Step 4: Convert FragPipe mods to PROSIT UNIMOD format
# =============================================================================
def parse_fragpipe_mods(modified_peptide, assigned_mods):
    """
    Convert FragPipe Modified Peptide to PROSIT-formatted sequence.
    Returns PROSIT-ready sequence string, or None if unparseable/unsupported.
    """
    if pd.isna(modified_peptide):
        return None

    mod_pep = str(modified_peptide).strip()

    # --- N-terminal modification ---
    nterm_tag = "[UNIMOD:737]-"  # default: TMT
    bare_seq = mod_pep

    nterm_match = re.match(r'^n\[(\d+)\](.+)$', mod_pep)
    if nterm_match:
        nterm_mass = int(nterm_match.group(1))
        bare_seq = nterm_match.group(2)
        if nterm_mass == 43:
            # N-term acetylation — not supported by PROSIT TMT model
            return None
        elif nterm_mass == 230:
            nterm_tag = "[UNIMOD:737]-"
        else:
            nterm_tag = "[UNIMOD:737]-"

    # --- Build residue-by-residue ---
    result = []
    i = 0
    while i < len(bare_seq):
        aa = bare_seq[i]
        inline_mod = None
        if i + 1 < len(bare_seq) and bare_seq[i + 1] == '[':
            bracket_end = bare_seq.index(']', i + 1)
            mass_str = bare_seq[i + 2:bracket_end]
            inline_mod = int(float(mass_str))
            i = bracket_end + 1
        else:
            i += 1

        if aa == 'C':
            result.append("C[UNIMOD:4]")
        elif aa == 'K':
            result.append("K[UNIMOD:737]")
        elif aa == 'M' and inline_mod == 147:
            result.append("M[UNIMOD:35]")
        else:
            result.append(aa)

    return nterm_tag + "".join(result)


# =============================================================================
# Phase 1 — Step 5: Koina API predictions
# =============================================================================
def predict_with_koina(peptide_seqs, charges, nce=38.0, batch_size=500):
    """
    Batched gRPC predictions via Koina.
    Returns dict with 'intensities', 'mz', 'annotation' arrays.
    """
    from koinapy import Koina
    print(f"  Connecting to Koina ({SERVER_URL})...")
    model = Koina(MODEL_NAME, server_url=SERVER_URL, ssl=True)

    n = len(peptide_seqs)
    all_intensities, all_mz, all_annotations = [], [], []
    failed_indices = []

    for start in tqdm(range(0, n, batch_size), desc="  PROSIT predictions",
                      unit="batch"):
        end = min(start + batch_size, n)
        batch_peps = peptide_seqs[start:end]
        batch_charges = charges[start:end]

        inputs = {
            "peptide_sequences": np.array(batch_peps, dtype=object
                                          ).reshape(-1, 1),
            "precursor_charges": np.array(batch_charges, dtype=np.int32
                                          ).reshape(-1, 1),
            "collision_energies": np.full((end - start, 1), nce,
                                         dtype=np.float32),
            "fragmentation_types": np.array(
                [FRAGMENTATION] * (end - start), dtype=object
            ).reshape(-1, 1),
        }

        max_retries = 3
        for attempt in range(max_retries):
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    result = model.predict(inputs)
                all_intensities.append(result["intensities"])
                all_mz.append(result["mz"])
                all_annotations.append(result["annotation"])
                break
            except Exception as e:
                err_msg = str(e)
                if attempt < max_retries - 1 and "not supported" not in err_msg:
                    wait = 5 * (attempt + 1)
                    print(f"\n  Batch {start}-{end} failed "
                          f"(attempt {attempt+1}): {err_msg[:80]}")
                    print(f"  Retrying in {wait}s...")
                    time.sleep(wait)
                else:
                    # Per-peptide fallback
                    print(f"\n  Batch {start}-{end} failed. "
                          f"Falling back to per-peptide prediction...")
                    for idx in range(len(batch_peps)):
                        single_in = {
                            "peptide_sequences": np.array(
                                [batch_peps[idx]], dtype=object
                            ).reshape(1, 1),
                            "precursor_charges": np.array(
                                [batch_charges[idx]], dtype=np.int32
                            ).reshape(1, 1),
                            "collision_energies": np.full(
                                (1, 1), nce, dtype=np.float32),
                            "fragmentation_types": np.array(
                                [FRAGMENTATION], dtype=object
                            ).reshape(1, 1),
                        }
                        try:
                            with warnings.catch_warnings():
                                warnings.simplefilter("ignore")
                                r = model.predict(single_in)
                            all_intensities.append(r["intensities"])
                            all_mz.append(r["mz"])
                            all_annotations.append(r["annotation"])
                        except Exception as e2:
                            print(f"    SKIP: {batch_peps[idx][:50]} "
                                  f"— {str(e2)[:60]}")
                            all_intensities.append(np.full(
                                (1, 174), np.nan, dtype=np.float32))
                            all_mz.append(np.full(
                                (1, 174), np.nan, dtype=np.float32))
                            all_annotations.append(np.full(
                                (1, 174), b"", dtype=object))
                            failed_indices.append(start + idx)
                    break

    predictions = {
        "intensities": np.vstack(all_intensities),
        "mz": np.vstack(all_mz),
        "annotation": np.vstack(all_annotations),
    }

    n_ok = n - len(failed_indices)
    print(f"  Predictions received: {n_ok}/{n} succeeded")
    if failed_indices:
        print(f"  WARNING: {len(failed_indices)} peptides failed prediction")

    return predictions


# =============================================================================
# Phase 1 — Step 6: Write MSP spectral library
# =============================================================================
def write_msp(output_path, peptides, charges, prosit_seqs, predictions):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    intensities = predictions["intensities"]
    mz_arr = predictions["mz"]
    ann_arr = predictions["annotation"]

    count = 0
    with open(output_path, "w") as f:
        for i in range(len(peptides)):
            spec_int = intensities[i]
            spec_mz = mz_arr[i]
            spec_ann = ann_arr[i]

            if np.all(np.isnan(spec_int)):
                continue

            valid = ((spec_mz > 0) & (spec_int > MIN_INTENSITY)
                     & (~np.isnan(spec_int)))
            if not np.any(valid):
                continue

            filt_mz = spec_mz[valid]
            filt_int = spec_int[valid]
            filt_ann = spec_ann[valid]

            max_int = filt_int.max()
            if max_int > 0:
                filt_int = filt_int / max_int * 10000.0

            order = np.argsort(filt_mz)
            filt_mz = filt_mz[order]
            filt_int = filt_int[order]
            filt_ann = filt_ann[order]

            f.write(f"Name: {peptides[i]}/{charges[i]}\n")
            f.write(f"Comment: PROSIT_seq={prosit_seqs[i]} NCE={NCE} "
                    f"Model={MODEL_NAME} Frag={FRAGMENTATION}\n")
            f.write(f"Num Peaks: {len(filt_mz)}\n")

            for j in range(len(filt_mz)):
                ann_str = filt_ann[j]
                if isinstance(ann_str, bytes):
                    ann_str = ann_str.decode("utf-8", errors="replace")
                f.write(f"{filt_mz[j]:.4f}\t{filt_int[j]:.1f}"
                        f"\t\"{ann_str}\"\n")

            f.write("\n")
            count += 1

    print(f"  Wrote {count} spectra to {os.path.basename(output_path)}")
    return count


# =============================================================================
# Phase 2 — MSP parser
# =============================================================================
def parse_msp(msp_path):
    """Parse MSP file -> dict keyed by 'PEPTIDE/charge'."""
    spectra = {}
    current_name = None
    current_mz, current_int, current_ann = [], [], []

    with open(msp_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                if current_name and current_mz:
                    spectra[current_name] = {
                        "mz": np.array(current_mz, dtype=np.float64),
                        "intensity": np.array(current_int, dtype=np.float64),
                        "annotation": current_ann,
                    }
                current_name = None
                current_mz, current_int, current_ann = [], [], []
                continue
            if line.startswith("Name:"):
                current_name = line.split("Name:")[1].strip()
            elif line.startswith("Comment:") or line.startswith("Num Peaks:"):
                continue
            else:
                parts = line.split("\t")
                if len(parts) >= 2:
                    try:
                        current_mz.append(float(parts[0]))
                        current_int.append(float(parts[1]))
                        ann = parts[2].strip('"') if len(parts) > 2 else ""
                        current_ann.append(ann)
                    except ValueError:
                        pass

    if current_name and current_mz:
        spectra[current_name] = {
            "mz": np.array(current_mz, dtype=np.float64),
            "intensity": np.array(current_int, dtype=np.float64),
            "annotation": current_ann,
        }
    return spectra


# =============================================================================
# Phase 2 — Observed spectrum extraction
# =============================================================================
def extract_observed_spectrum_fast(mzml_path, scan_number):
    """Extract a single MS2 spectrum from mzML by scan number."""
    from pyteomics import mzml
    if not os.path.isfile(mzml_path):
        return None, None
    try:
        with mzml.MzML(mzml_path) as reader:
            for target_id in [
                f"controllerType=0 controllerNumber=1 scan={scan_number}",
                f"scan={scan_number}",
            ]:
                try:
                    spectrum = reader.get_by_id(target_id)
                    mz = spectrum.get("m/z array", np.array([]))
                    intensity = spectrum.get("intensity array", np.array([]))
                    if len(mz) > 0:
                        # Exclude TMT reporter ion window
                        mask = ~((mz >= CLEAR_MZ_LO) & (mz <= CLEAR_MZ_HI))
                        return (mz[mask].astype(np.float64),
                                intensity[mask].astype(np.float64))
                except KeyError:
                    continue
    except Exception:
        pass

    # Sequential fallback
    try:
        with mzml.MzML(mzml_path) as reader:
            for spectrum in reader:
                if f"scan={scan_number}" in spectrum.get("id", ""):
                    mz = spectrum.get("m/z array", np.array([]))
                    intensity = spectrum.get("intensity array", np.array([]))
                    if len(mz) > 0:
                        mask = ~((mz >= CLEAR_MZ_LO) & (mz <= CLEAR_MZ_HI))
                        return (mz[mask].astype(np.float64),
                                intensity[mask].astype(np.float64))
                    return None, None
    except Exception:
        pass
    return None, None


# =============================================================================
# Phase 2 — Peak matching & spectral angle
# =============================================================================
def match_peaks(pred_mz, pred_int, obs_mz, obs_int, ppm_tol=20.0):
    matched_pred, matched_obs, matched_idxs = [], [], []
    used_obs = set()  # prevent one observed peak matching multiple predicted ions
    for i, p_mz in enumerate(pred_mz):
        tol = p_mz * ppm_tol / 1e6
        diffs = np.abs(obs_mz - p_mz)
        min_idx = np.argmin(diffs)
        if diffs[min_idx] <= tol and min_idx not in used_obs:
            matched_pred.append(pred_int[i])
            matched_obs.append(obs_int[min_idx])
            matched_idxs.append(i)
            used_obs.add(min_idx)
    return (np.array(matched_pred), np.array(matched_obs),
            len(matched_pred), len(pred_mz), matched_idxs)


def compute_spectral_angle(pred_int, obs_int):
    if len(pred_int) == 0 or len(obs_int) == 0:
        return 90.0, 0.0
    pred_norm = np.sqrt(np.sum(pred_int ** 2))
    obs_norm = np.sqrt(np.sum(obs_int ** 2))
    if pred_norm == 0 or obs_norm == 0:
        return 90.0, 0.0
    cos_sim = np.clip(
        np.dot(pred_int, obs_int) / (pred_norm * obs_norm), -1, 1)
    sa_rad = np.arccos(cos_sim)
    return np.degrees(sa_rad), 1.0 - (2.0 / np.pi) * sa_rad


# =============================================================================
# Phase 2 — Classification helpers
# =============================================================================
def classify_sa(sa_degrees):
    if sa_degrees <= SA_HIGH:     return "High"
    if sa_degrees <= SA_MODERATE: return "Moderate"
    if sa_degrees <= SA_POOR:     return "Poor"
    return "Terrible"


def classify_match_coverage(n_cleavages, total_sites):
    """Coverage based on unique cleavage sites, not fragment count."""
    if total_sites == 0:
        return "Sparse", 0.0
    pct = 100.0 * n_cleavages / total_sites
    if pct <= 25:
        return "Sparse", pct
    if pct <= 50:
        return "Partial", pct
    return "Majority", pct


def overall_confidence(sa_rating, match_rating):
    if sa_rating == "High":
        return "Strong" if match_rating != "Sparse" else "Moderate"
    if sa_rating == "Moderate":
        return "Weak" if match_rating == "Sparse" else "Moderate"
    if sa_rating == "Poor":
        return "Weak" if match_rating == "Majority" else "Insufficient"
    return "Insufficient"


def parse_ion_annotation(ann_str):
    m = re.match(r"([by])(\d+)\+(\d+)", ann_str)
    if m:
        return m.group(1), int(m.group(2)), int(m.group(3))
    return None, None, None


def compute_unique_cleavages(matched_idxs, pred_ann, peptide_len):
    """Count unique backbone cleavage sites evidenced by matched ions.
    b_i → site i;  y_i → site (peptide_len - i).
    Returns (n_unique, total_sites).
    """
    sites = set()
    for idx in matched_idxs:
        ion_t, ion_ord, _ = parse_ion_annotation(pred_ann[idx])
        if ion_t == "b":
            sites.add(ion_ord)
        elif ion_t == "y":
            sites.add(peptide_len - ion_ord)
    return len(sites), max(peptide_len - 1, 1)


def parse_spectrum_id(spectrum_id, batch, workspace):
    """
    Parse FragPipe Spectrum column: 'rushtmt_b1_21.59377.59377.3'
    or round2: 'rosmaptmt_r2_b1_01.12345.12345.2'.
    batch can be 'b1' or 'round2/b1'.
    Returns: (mzml_path, scan_number)
    """
    parts = spectrum_id.split(".")
    if len(parts) < 4:
        return None, None
    mzml_basename = parts[0]
    scan_number = int(parts[1])
    mzml_path = os.path.join(workspace, batch, f"{mzml_basename}.mzML")
    return mzml_path, scan_number


# =============================================================================
# Phase 2 — Worker function for multiprocessing
# =============================================================================
def _process_one_peptide(args):
    """
    Worker function for multiprocessing.Pool.
    Evaluates up to MAX_CANDIDATES observed spectra per peptide against the
    PROSIT prediction and picks the one with the best evidence score
    (Strong > Moderate > Weak > Insufficient), using SA as tiebreaker.
    Fully self-contained — no shared state.
    """
    (row_dict, candidate_dicts, pred_mz, pred_int, pred_ann,
     workspace, ppm_tol, plot_dir) = args

    # Higher = better evidence rank
    CONF_RANK = {"Strong": 4, "Moderate": 3, "Weak": 2, "Insufficient": 1,
                 "N/A": 0}

    peptide = row_dict["tryptic_peptide"]
    charge = int(row_dict["Charge"])
    protein_id = str(row_dict.get("protein_id", ""))

    result = {
        "index": row_dict["_index"],
        "SA_degrees": np.nan,
        "SA_normalized": np.nan,
        "SA_rating": "Error",
        "Match_coverage": "N/A",
        "Match_coverage_pct": 0.0,
        "Confidence": "N/A",
        "Matched_fragments": 0,
        "Total_predicted_fragments": 0,
        "Unique_cleavages": 0,
        "Total_cleavage_sites": 0,
        "Candidates_evaluated": 0,
        "Selected_Spectrum": "",
        "Selected_Batch": "",
        "status": "",
    }

    if pred_mz is None:
        result["SA_rating"] = "No prediction"
        result["status"] = "no_pred"
        return result

    # Filter out fragment ions with charge >= precursor charge (physically impossible)
    frag_mask = []
    for ann in pred_ann:
        _, _, ion_chg = parse_ion_annotation(ann)
        frag_mask.append(ion_chg is None or ion_chg < charge)
    frag_mask = np.array(frag_mask, dtype=bool)
    if not np.all(frag_mask):
        pred_mz  = pred_mz[frag_mask]
        pred_int = pred_int[frag_mask]
        pred_ann = [a for a, keep in zip(pred_ann, frag_mask) if keep]

    result["Total_predicted_fragments"] = len(pred_mz)

    # Evaluate each candidate and pick by most matched fragments,
    # then best evidence tier, then lowest SA as tiebreakers
    best_score = (-1, 91.0, -1)  # (n_clv, -sa_deg, n_match)
    best_candidate = None
    best_obs_data = None
    n_evaluated = 0

    for cand in candidate_dicts:
        cand_spectrum = cand["Spectrum"]
        cand_batch = cand["Batch"]

        mzml_path, scan_number = parse_spectrum_id(
            cand_spectrum, cand_batch, workspace)
        if mzml_path is None:
            continue

        obs_mz, obs_int = extract_observed_spectrum_fast(
            mzml_path, scan_number)
        if obs_mz is None or len(obs_mz) == 0:
            continue

        matched_pred, matched_obs, n_match, n_pred, m_idxs = match_peaks(
            pred_mz, pred_int, obs_mz, obs_int, ppm_tol=ppm_tol)
        sa_deg, sa_norm = compute_spectral_angle(matched_pred, matched_obs)
        n_clv, total_sites = compute_unique_cleavages(
            m_idxs, pred_ann, len(peptide))
        sa_rat = classify_sa(sa_deg)
        match_cov, _ = classify_match_coverage(n_clv, total_sites)
        conf = overall_confidence(sa_rat, match_cov)
        n_evaluated += 1

        # Unique cleavages first, then lowest SA, then fragment count
        cand_score = (n_clv, -sa_deg, n_match)
        if cand_score > best_score:
            best_score = cand_score
            best_candidate = cand
            best_obs_data = (obs_mz, obs_int, n_match, n_pred,
                             sa_deg, sa_norm, scan_number,
                             n_clv, total_sites)

    result["Candidates_evaluated"] = n_evaluated

    if best_candidate is None:
        result["SA_rating"] = "No observed"
        result["status"] = "no_obs_any_candidate"
        return result

    obs_mz, obs_int, n_match, n_pred, sa_deg, sa_norm, scan_number, \
        n_clv, total_sites = best_obs_data
    sa_rat = classify_sa(sa_deg)
    match_cov, match_pct = classify_match_coverage(n_clv, total_sites)
    conf = overall_confidence(sa_rat, match_cov)

    result.update({
        "SA_degrees": sa_deg,
        "SA_normalized": sa_norm,
        "SA_rating": sa_rat,
        "Match_coverage": match_cov,
        "Match_coverage_pct": round(match_pct, 1),
        "Confidence": conf,
        "Matched_fragments": n_match,
        "Total_predicted_fragments": n_pred,
        "Unique_cleavages": n_clv,
        "Total_cleavage_sites": total_sites,
        "Selected_Spectrum": best_candidate["Spectrum"],
        "Selected_Batch": best_candidate["Batch"],
        "status": "ok",
    })

    # Mirror plot — only for the best representative
    acetyl_rescued = bool(row_dict.get("Acetyl_Rescued", False))
    if plot_dir:
        try:
            _create_mirror_plot(
                peptide, charge, pred_mz, pred_int, pred_ann,
                obs_mz, obs_int, sa_deg, sa_norm, n_match,
                plot_dir, scan_number, ppm_tol,
                protein_id=protein_id,
                acetyl_rescued=acetyl_rescued,
                n_cleavages=n_clv, total_sites=total_sites)
        except Exception as e:
            result["status"] += f";plot_err:{str(e)[:40]}"

    return result


# =============================================================================
# Phase 2 — 3-panel mirror plot (called inside worker)
# =============================================================================
def _create_mirror_plot(peptide, charge, pred_mz, pred_int, pred_ann,
                        obs_mz, obs_int, sa_deg, sa_norm, n_matched,
                        plot_dir, scan_number, ppm_tol,
                        protein_id="", acetyl_rescued=False,
                        n_cleavages=0, total_sites=0):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    # -- Aesthetics --
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Helvetica", "Arial", "DejaVu Sans"]
    plt.rcParams["pdf.fonttype"] = 42   # TrueType — editable in Illustrator
    plt.rcParams["ps.fonttype"] = 42

    B_CLR, Y_CLR = "#3366CC", "#CC3333"
    UNMATCH_CLR = "#B0B0B0"
    EVIDENCE_CLR = {
        "Strong": "#2E8B57", "Moderate": "#D4880F",
        "Weak": "#C44E52", "Insufficient": "#8C8C8C",
    }

    fig, (ax_seq, ax_spec, ax_ppm) = plt.subplots(
        3, 1, figsize=(10, 10),
        gridspec_kw={"height_ratios": [0.6, 4, 1.0], "hspace": 0.12},
    )
    fig.subplots_adjust(left=0.10, right=0.95, top=0.85, bottom=0.08)

    # Parse annotations
    ion_types, ion_ords = [], []
    for ann in pred_ann:
        t, o, _ = parse_ion_annotation(ann)
        ion_types.append(t)
        ion_ords.append(o)

    # Build match info
    match_info = []
    used_obs_plot = set()
    for i, p_mz in enumerate(pred_mz):
        tol = p_mz * ppm_tol / 1e6
        diffs = np.abs(obs_mz - p_mz)
        min_idx = np.argmin(diffs)
        if diffs[min_idx] <= tol and min_idx not in used_obs_plot:
            signed_ppm = (obs_mz[min_idx] - p_mz) / p_mz * 1e6
            match_info.append({
                "pred_idx": i, "obs_idx": min_idx,
                "pred_mz": p_mz, "obs_mz": obs_mz[min_idx],
                "ann": pred_ann[i] if i < len(pred_ann) else "",
                "ion_type": ion_types[i], "ion_ord": ion_ords[i],
                "signed_ppm": signed_ppm,
            })
            used_obs_plot.add(min_idx)

    matched_b = {m["ion_ord"] for m in match_info if m["ion_type"] == "b"}
    matched_y = {m["ion_ord"] for m in match_info if m["ion_type"] == "y"}

    # -- Panel 1: Peptide Sequence --
    n_aa = len(peptide)
    # Adaptive sizing: shrink for long peptides
    if n_aa <= 15:
        aa_fs, aa_pad = 20, 0.3
    elif n_aa <= 25:
        aa_fs, aa_pad = 16, 0.2
    else:
        aa_fs, aa_pad = 12, 0.15
    for i, aa in enumerate(peptide):
        # Tint matched residues, white for unmatched
        has_b = (i + 1) in matched_b
        has_y = (n_aa - i) in matched_y
        if has_b and has_y:
            box_fc = "#E8E0F0"   # lavender (both)
        elif has_b:
            box_fc = "#DCE6F1"   # light blue
        elif has_y:
            box_fc = "#F5E0E0"   # light rose
        else:
            box_fc = "white"
        ax_seq.text(
            i + 0.5, 0.5, aa, ha="center", va="center", fontsize=aa_fs,
            fontweight="bold",
            bbox=dict(boxstyle=f"round,pad={aa_pad}", facecolor=box_fc,
                      edgecolor="#BBBBBB", linewidth=0.8, alpha=0.95))

    ion_fs = 14 if n_aa <= 15 else (11 if n_aa <= 25 else 8)
    ion_lw = 2.5 if n_aa <= 15 else (1.8 if n_aa <= 25 else 1.2)
    for i in range(1, n_aa):
        c = B_CLR if i in matched_b else "#E0E0E0"
        ax_seq.plot([i - 0.1] * 2, [0.1, 0.9], color=c, linewidth=ion_lw)
        ax_seq.text(i - 0.1, 0.02, f"b{i}", ha="center", va="top",
                    fontsize=ion_fs, color=c, fontweight="bold")
    for i in range(1, n_aa):
        c = Y_CLR if i in matched_y else "#E0E0E0"
        ax_seq.plot([n_aa - i + 0.1] * 2, [0.1, 0.9], color=c, linewidth=ion_lw)
        ax_seq.text(n_aa - i + 0.1, 0.98, f"y{i}", ha="center", va="bottom",
                    fontsize=ion_fs, color=c, fontweight="bold")

    ax_seq.set_xlim(0, n_aa)
    ax_seq.set_ylim(-0.15, 1.15)
    ax_seq.axis("off")

    match_cov, match_pct = classify_match_coverage(n_cleavages, total_sites)
    sa_rat = classify_sa(sa_deg)
    conf = overall_confidence(sa_rat, match_cov)

    # Two-tier title: bold peptide (black) + colored metrics
    title_top = f"{peptide}/{charge}+"
    if protein_id and protein_id != "nan":
        title_top += f"  |  {protein_id}"
    fig.suptitle(title_top, fontsize=22, fontweight="bold", y=0.95)

    title_line2 = (
        f"SA = {sa_deg:.1f}\u00b0 ({sa_rat})  \u2022  "
        f"Cleavages: {n_cleavages}/{total_sites} ({match_pct:.0f}%, {match_cov})  \u2022  "
        f"Evidence: {conf}")
    if acetyl_rescued:
        title_line2 += "  |  \u26a0 Acetyl-rescued PSM"
    ecol = EVIDENCE_CLR.get(conf, "#333333")

    # -- Panel 2: Mirror Spectrum --
    obs_max = obs_int.max() if len(obs_int) > 0 else 1.0
    pred_max = pred_int.max() if len(pred_int) > 0 else 1.0
    obs_norm = obs_int / obs_max * 100 if obs_max > 0 else obs_int
    pred_norm = pred_int / pred_max * 100 if pred_max > 0 else pred_int

    # Gold mirror shading - lighter observed, deeper predicted
    ax_spec.axhspan(0, 145, color="#FFE49C", alpha=0.5, zorder=0)
    ax_spec.axhspan(-145, 0, color="#FFF3D4", alpha=0.4, zorder=0)

    # Unmatched peaks
    ax_spec.vlines(obs_mz, 0, obs_norm, linewidth=0.5, color=UNMATCH_CLR)
    ax_spec.vlines(pred_mz, 0, -pred_norm, linewidth=0.5, color=UNMATCH_CLR)

    # Matched peaks
    for m in match_info:
        c = (B_CLR if m["ion_type"] == "b"
             else (Y_CLR if m["ion_type"] == "y" else "gray"))
        ax_spec.vlines([m["obs_mz"]], [0], [obs_norm[m["obs_idx"]]],
                       linewidth=1.5, color=c)
        ax_spec.vlines([m["pred_mz"]], [0], [-pred_norm[m["pred_idx"]]],
                       linewidth=1.5, color=c)

    # Annotations: ALL matched peaks, greedy overlap avoidance
    if match_info:
        x_range = ax_spec.get_xlim()[1] - ax_spec.get_xlim()[0]
        if x_range == 0:
            x_range = 1.0
        min_dx = x_range * 0.018
        min_dy = 5.0
        sorted_matches = sorted(match_info,
                                key=lambda m: obs_norm[m["obs_idx"]],
                                reverse=True)
        placed = []
        for m in sorted_matches:
            x = m["obs_mz"]
            y = obs_norm[m["obs_idx"]]
            too_close = any(abs(x - px) < min_dx and abs(y - py) < min_dy
                           for px, py in placed)
            if too_close:
                continue
            placed.append((x, y))
            c = (B_CLR if m["ion_type"] == "b"
                 else (Y_CLR if m["ion_type"] == "y" else "gray"))
            ax_spec.text(
                x, y, f" {m['ann']}",
                rotation=45, va="bottom", ha="left",
                fontsize=14, color=c, clip_on=True)

    ax_spec.axhline(y=0, color="black", linewidth=0.8)
    ax_spec.set_ylabel("Relative Intensity (%)", fontsize=20)
    ax_spec.yaxis.set_major_formatter(
        plt.FuncFormatter(lambda y, _: f"{abs(y):.0f}"))

    # Subtle horizontal gridlines
    ax_spec.yaxis.grid(True, alpha=0.12, linestyle="--", color="gray", zorder=0)
    ax_spec.set_axisbelow(True)

    # Remove right spine only (keep top)
    ax_spec.spines["right"].set_visible(False)

    legend_h = [
        Line2D([0], [0], color=B_CLR, linewidth=1.5, label="b-ions"),
        Line2D([0], [0], color=Y_CLR, linewidth=1.5, label="y-ions"),
        Line2D([0], [0], color=UNMATCH_CLR, linewidth=0.8, label="Unmatched"),
    ]
    ax_spec.legend(handles=legend_h, loc="upper right", fontsize=17,
                   framealpha=0.9, edgecolor="#CCCCCC")
    ax_spec.set_ylim(-140, 140)
    ax_spec.tick_params(labelsize=17)
    ax_spec.text(0.01, 0.97, "Observed", transform=ax_spec.transAxes,
                 fontsize=14, va="top", ha="left", color="0.35",
                 fontstyle="italic")
    ax_spec.text(0.01, 0.12, title_line2, transform=ax_spec.transAxes,
                 fontsize=14, va="bottom", ha="left", color=ecol,
                 fontstyle="italic")
    ax_spec.text(0.01, 0.03, "Predicted (PROSIT)",
                 transform=ax_spec.transAxes,
                 fontsize=14, va="bottom", ha="left", color="0.35",
                 fontstyle="italic")

    # -- Panel 3: PPM Error (lollipop) --
    ax_ppm.set_facecolor("#E8F0FE")
    ax_ppm.patch.set_alpha(0.5)
    ax_ppm.spines["right"].set_visible(False)

    if match_info:
        for ion_t, clr in [("b", B_CLR), ("y", Y_CLR)]:
            mzs = [m["pred_mz"] for m in match_info
                    if m["ion_type"] == ion_t]
            ppms = [m["signed_ppm"] for m in match_info
                    if m["ion_type"] == ion_t]
            if mzs:
                # Lollipop: thin stem from zero + dot cap
                ax_ppm.vlines(mzs, [0] * len(mzs), ppms,
                              colors=clr, linewidth=0.6, alpha=0.45, zorder=2)
                ax_ppm.scatter(mzs, ppms, c=clr, s=80, alpha=0.8,
                               label=f"{ion_t}-ions", edgecolors="white",
                               linewidths=0.5, zorder=3)
        o_mz = [m["pred_mz"] for m in match_info
                if m["ion_type"] not in ("b", "y")]
        o_ppm = [m["signed_ppm"] for m in match_info
                 if m["ion_type"] not in ("b", "y")]
        if o_mz:
            ax_ppm.vlines(o_mz, [0] * len(o_mz), o_ppm,
                          colors="gray", linewidth=0.5, alpha=0.35, zorder=2)
            ax_ppm.scatter(o_mz, o_ppm, c="gray", s=50, alpha=0.5,
                           label="other", edgecolors="none", zorder=3)
        ax_ppm.axhline(0, color="black", linewidth=0.5)
        ax_ppm.axhline(ppm_tol, color="gray", linewidth=0.4,
                       linestyle="--", alpha=0.6)
        ax_ppm.axhline(-ppm_tol, color="gray", linewidth=0.4,
                       linestyle="--", alpha=0.6)
        ax_ppm.set_ylim(-25, 25)
        ax_ppm.set_yticks([-20, 0, 20])
        ax_ppm.set_ylabel("ppm error", fontsize=20)
        ax_ppm.set_xlabel("m/z", fontsize=20)
        ax_ppm.legend(fontsize=17, loc="upper right", edgecolor="#CCCCCC")
        ax_ppm.tick_params(labelsize=17)
        ax_ppm.set_xlim(ax_spec.get_xlim())
    else:
        ax_ppm.text(0.5, 0.5, "No matched peaks", ha="center", va="center",
                    fontsize=17, color="gray", transform=ax_ppm.transAxes)
        ax_ppm.set_xlabel("m/z", fontsize=20)
        ax_ppm.set_ylabel("ppm error", fontsize=20)

    # Save into confidence subdirectory
    conf_subdir = os.path.join(plot_dir, conf)
    os.makedirs(conf_subdir, exist_ok=True)
    safe_name = re.sub(r'[^\w]', '_', peptide)
    out_stem = f"{safe_name}_{charge}_{scan_number}"
    plt.savefig(os.path.join(conf_subdir, out_stem + ".pdf"),
               dpi=300, bbox_inches="tight")
    plt.savefig(os.path.join(conf_subdir, out_stem + ".png"),
               dpi=300, bbox_inches="tight")
    plt.close()


# =============================================================================
# Phase 2 — Main SA computation with multiprocessing
# =============================================================================
def run_phase2(csv_path, msp_path, out_csv, plot_dir, n_workers,
               no_plots=False, cand_csv=None):
    print("\n" + "=" * 65)
    print(" Phase 2: Spectral Angle Computation (multiprocessing)")
    print("  Representative selection: best SA from up to "
          f"{MAX_CANDIDATES} candidates (Purity >= {MIN_PURITY})")
    print("=" * 65)

    # Step 1: Read Phase 1 CSV
    print("\n[1/5] Reading Phase 1 results...")
    df = pd.read_csv(csv_path)
    print(f"  {len(df)} peptides to process")

    # Load candidate PSMs for representative selection
    cand_lookup = {}  # peptide -> list of candidate row dicts
    if cand_csv and os.path.isfile(cand_csv):
        cand_df = pd.read_csv(cand_csv)
        for peptide, grp in cand_df.groupby("Peptide"):
            cand_lookup[peptide] = grp.sort_values(
                "Candidate_Rank").to_dict("records")
        n_multi = sum(1 for v in cand_lookup.values() if len(v) > 1)
        print(f"  Loaded {len(cand_df)} candidate PSMs for "
              f"{len(cand_lookup)} peptides ({n_multi} with alternatives)")
    else:
        print("  No candidate PSMs file found — "
              "using primary PSM only (original behavior)")

    # Step 2: Parse MSP
    print("\n[2/5] Parsing predicted spectra from MSP...")
    pred_spectra = parse_msp(msp_path)
    print(f"  {len(pred_spectra)} predicted spectra loaded")

    # Step 3: Prepare worker tasks
    print(f"\n[3/5] Computing spectral angles ({n_workers} workers)...")
    if not no_plots:
        os.makedirs(plot_dir, exist_ok=True)
        print(f"  Mirror plots -> {plot_dir}")

    tasks = []
    for idx, row in df.iterrows():
        peptide = row["tryptic_peptide"]
        charge = int(row["Charge"])
        pred_key = f"{peptide}/{charge}"

        pred = pred_spectra.get(pred_key)
        if pred:
            p_mz = pred["mz"]
            p_int = pred["intensity"]
            p_ann = pred["annotation"]
        else:
            p_mz, p_int, p_ann = None, None, None

        row_dict = row.to_dict()
        row_dict["_index"] = idx

        # Build candidate list for this peptide
        candidates = cand_lookup.get(peptide, [])
        if not candidates:
            # Fallback: use the primary row itself as the sole candidate
            candidates = [{"Spectrum": row["Spectrum"],
                           "Batch": row["Batch"]}]

        tasks.append((
            row_dict, candidates, p_mz, p_int, p_ann,
            WORKSPACE, PPM_TOLERANCE,
            plot_dir if not no_plots else None,
        ))

    # Run with multiprocessing
    results = []
    with Pool(processes=n_workers) as pool:
        for res in tqdm(pool.imap_unordered(_process_one_peptide, tasks),
                        total=len(tasks), desc="  Processing", unit="PSM"):
            results.append(res)

    # Sort back by index
    results.sort(key=lambda r: r["index"])

    # Step 4: Augment CSV
    print("\n[4/5] Saving augmented CSV...")
    df["SA_degrees"] = [r["SA_degrees"] for r in results]
    df["SA_normalized"] = [r["SA_normalized"] for r in results]
    df["SA_rating"] = [r["SA_rating"] for r in results]
    df["Matched_fragments"] = [r["Matched_fragments"] for r in results]
    df["Total_predicted_fragments"] = [r["Total_predicted_fragments"]
                                       for r in results]
    df["Unique_cleavages"] = [r["Unique_cleavages"] for r in results]
    df["Total_cleavage_sites"] = [r["Total_cleavage_sites"]
                                  for r in results]
    df["Match_coverage"] = [r["Match_coverage"] for r in results]
    df["Match_coverage_pct"] = [r["Match_coverage_pct"] for r in results]
    df["Confidence"] = [r["Confidence"] for r in results]
    df["Candidates_evaluated"] = [r["Candidates_evaluated"] for r in results]
    df["Selected_Spectrum"] = [r["Selected_Spectrum"] for r in results]
    df["Selected_Batch"] = [r["Selected_Batch"] for r in results]

    df.to_csv(out_csv, index=False)
    print(f"  Saved -> {os.path.basename(out_csv)}")

    # Representative selection stats
    n_improved = sum(1 for r in results
                     if r["Candidates_evaluated"] > 1
                     and r["Selected_Spectrum"]
                     and r["Selected_Spectrum"] != "")
    n_reselected = sum(
        1 for i, r in enumerate(results)
        if r["Candidates_evaluated"] > 1
        and r["Selected_Spectrum"]
        and r["Selected_Spectrum"] != df.iloc[i]["Spectrum"]
    )
    if n_improved > 0:
        print(f"\n  Representative selection:")
        print(f"    Peptides with multiple candidates evaluated: {n_improved}")
        print(f"    Re-selected (different from rank-1):        {n_reselected}")

    # Log issues
    errors = [r for r in results if r["status"] and r["status"] != "ok"]
    if errors:
        print(f"\n  Notes: {len(errors)} peptides had issues:")
        for status, count in Counter(
            r["status"].split(";")[0] for r in errors
        ).most_common():
            print(f"    {status}: {count}")

    # Step 5: Summary
    print("\n[5/5] Summary")
    valid = df.dropna(subset=["SA_degrees"])
    if len(valid) > 0:
        print(f"  Peptides with SA scores: {len(valid)}/{len(df)}")
        print(f"  SA (degrees): mean={valid['SA_degrees'].mean():.1f}\u00b0, "
              f"median={valid['SA_degrees'].median():.1f}\u00b0")
        print(f"  SA (normalized): "
              f"mean={valid['SA_normalized'].mean():.3f}, "
              f"median={valid['SA_normalized'].median():.3f}")

        # OHW vs multi-PSM breakdown
        if "Is_OHW" in valid.columns:
            for label, mask in [("One-hit-wonders", valid["Is_OHW"]),
                                ("Multi-PSM", ~valid["Is_OHW"])]:
                sub = valid[mask]
                if len(sub) > 0:
                    print(f"\n  {label} ({len(sub)} peptides):")
                    print(f"    SA mean={sub['SA_degrees'].mean():.1f}\u00b0, "
                          f"median={sub['SA_degrees'].median():.1f}\u00b0")

        print(f"\n  SA rating distribution:")
        for rating in ["High", "Moderate", "Poor", "Terrible"]:
            n = (valid["SA_rating"] == rating).sum()
            pct = 100 * n / len(valid) if len(valid) > 0 else 0
            print(f"    {rating:10s}: {n:4d} ({pct:5.1f}%)")
        print(f"\n  Match coverage distribution:")
        for cov in ["Sparse", "Partial", "Majority"]:
            n = (valid["Match_coverage"] == cov).sum()
            pct = 100 * n / len(valid) if len(valid) > 0 else 0
            print(f"    {cov:10s}: {n:4d} ({pct:5.1f}%)")
        print(f"\n  Overall evidence distribution:")
        for conf in ["Strong", "Moderate", "Weak", "Insufficient"]:
            n = (valid["Confidence"] == conf).sum()
            pct = 100 * n / len(valid) if len(valid) > 0 else 0
            print(f"    {conf:10s}: {n:4d} ({pct:5.1f}%)")
    else:
        print("  No valid SA scores computed.")

    if not no_plots:
        n_plots = sum(
            len([f for f in files if f.endswith((".pdf", ".png"))])
            for _, _, files in os.walk(plot_dir)
        )
        print(f"\n  Mirror plots saved: {n_plots} in {plot_dir}")

    print("\n" + "=" * 65)
    print(" Phase 2 complete!")
    print(f"   Output CSV:    {out_csv}")
    if not no_plots:
        print(f"   Mirror plots:  {plot_dir}")
    print("=" * 65)


# =============================================================================
# Main
# =============================================================================
def main():
    parser = argparse.ArgumentParser(
        description="PROSIT TMT Pipeline — predict spectra + compute SA")
    parser.add_argument("--test", action="store_true",
                        help="Test mode: 20 peptides, 4 workers")
    parser.add_argument("--phase1", action="store_true",
                        help="Run Phase 1 only (prediction)")
    parser.add_argument("--phase2", action="store_true",
                        help="Run Phase 2 only (SA computation)")
    parser.add_argument("--workers", type=int, default=None,
                        help="Number of parallel workers for Phase 2 "
                             "(default: min(32, cpu_count))")
    parser.add_argument("--no-plots", action="store_true",
                        help="Skip mirror plot generation")
    parser.add_argument("--annotate-only", action="store_true",
                        help="Only append bracketed metrics to the master "
                             "CSV from existing all_peptides_with_SA CSV")
    args = parser.parse_args()

    test_mode = args.test
    n_workers = args.workers or min(32, cpu_count())

    # Determine which phases to run
    run_p1 = not args.phase2  # Phase 1 unless --phase2 only
    run_p2 = not args.phase1  # Phase 2 unless --phase1 only
    if args.annotate_only:
        run_p1 = False
        run_p2 = False

    # Paths based on mode
    if test_mode:
        csv_p1 = OUTPUT_CSV_P1_TEST
        msp_out = MSP_OUTPUT_TEST
        csv_p2 = OUTPUT_CSV_P2_TEST
        csv_ann = OUTPUT_CSV_ANNOTATED_TEST
        plot_d = PLOT_DIR_TEST
        cand_csv = CANDIDATE_CSV_TEST
        test_n = 20
        n_workers = min(4, n_workers)
    else:
        csv_p1 = OUTPUT_CSV_P1
        msp_out = MSP_OUTPUT
        csv_p2 = OUTPUT_CSV_P2
        csv_ann = OUTPUT_CSV_ANNOTATED
        plot_d = PLOT_DIR
        cand_csv = CANDIDATE_CSV

    # =========================================================================
    # Phase 1
    # =========================================================================
    if run_p1:
        print("=" * 65)
        print(" Phase 1: PROSIT TMT Spectral Prediction (all peptides)")
        if test_mode:
            print(f" *** TEST MODE: {test_n} peptides, "
                  f"{n_workers} workers ***")
        print("=" * 65)

        # Step 1: Read CSV
        print("\n[1/6] Reading tryptic peptide CSV...")
        tryptic_df, master_df = read_tryptic_peptides(INPUT_CSV)
        peptide_set = set(tryptic_df["tryptic_peptide"].tolist())

        # Step 2: Find ALL peptides across batches
        print("\n[2/6] Scanning all batches for peptides (round 1 + 2)...")
        all_pep_df, global_counts = find_all_peptides(WORKSPACE, peptide_set)
        if all_pep_df.empty:
            print("  No peptides found. Exiting.")
            return

        # Step 3: Get candidate PSMs per peptide (purity-filtered, top-5)
        print(f"\n[3/6] Retrieving candidate PSMs per peptide "
              f"(Purity >= {MIN_PURITY}, top {MAX_CANDIDATES})...")
        all_candidates_df = get_candidate_psms(WORKSPACE, all_pep_df)
        if all_candidates_df.empty:
            print("  No PSM details retrieved. Exiting.")
            return

        # Save full candidate list for Phase 2 representative selection
        all_candidates_df.to_csv(cand_csv, index=False)
        print(f"  Saved candidate PSMs -> {os.path.basename(cand_csv)}")

        # Use rank-1 candidates as primary for PROSIT prediction
        best_psm_df = all_candidates_df[
            all_candidates_df["Candidate_Rank"] == 1
        ].copy().reset_index(drop=True)

        # Add metadata columns
        best_psm_df["Global_Spectral_Count"] = (
            best_psm_df["Peptide"].map(global_counts))
        best_psm_df["Is_OHW"] = best_psm_df["Global_Spectral_Count"] == 1
        best_psm_df["Round"] = best_psm_df["Batch"].apply(
            lambda b: 2 if b.startswith("round2/") else 1)

        # Exclude peptides whose rank-1 is still N-term acetylated
        # (these are peptides with NO non-acetylated PSMs anywhere)
        has_ac = best_psm_df["Modified Peptide"].str.contains(
            r'^n\[43\]', na=False)
        n_acetyl = has_ac.sum()
        if n_acetyl > 0:
            best_psm_df = best_psm_df[~has_ac].copy()
            print(f"  Excluded {n_acetyl} peptides with only n[43] PSMs "
                  f"(no non-acetylated alternative found)")

        # Test mode sampling
        if test_mode:
            has_ox = best_psm_df["Modified Peptide"].str.contains(
                r'\[147\]', na=False)
            has_cys = best_psm_df["Peptide"].str.contains('C', na=False)
            is_ohw = best_psm_df["Is_OHW"]
            is_r2 = best_psm_df["Round"] == 2

            samples = []
            for mask, label, n_take in [
                (is_ohw & ~is_r2, "OHW_r1", 5),
                (~is_ohw & ~is_r2, "multi_r1", 5),
                (is_ohw & is_r2, "OHW_r2", 3),
                (~is_ohw & is_r2, "multi_r2", 3),
                (has_ox, "oxidized", 2),
                (has_cys, "cys", 2),
            ]:
                subset = best_psm_df[mask]
                subset = subset[
                    subset["Peptide"].str.len() <= MAX_PEPTIDE_LEN]
                if len(subset) > 0:
                    samples.append(subset.head(n_take))
            best_psm_df = pd.concat(samples).drop_duplicates(
                subset="Peptide").head(test_n)
            print(f"  Test mode: selected {len(best_psm_df)} diverse "
                  f"peptides")
            for _, row in best_psm_df.iterrows():
                r = "R2" if row["Round"] == 2 else "R1"
                ohw = "OHW" if row["Is_OHW"] else f"SC={row['Global_Spectral_Count']}"
                print(f"    [{r}] [{ohw}] {row['Modified Peptide']}")

        # Step 4: Convert mods
        print("\n[4/6] Converting modifications to PROSIT format...")
        prosit_seqs = []
        valid_mask = []
        for _, row in best_psm_df.iterrows():
            seq = parse_fragpipe_mods(row["Modified Peptide"],
                                       row["Assigned Modifications"])
            if seq is None:
                prosit_seqs.append(None)
                valid_mask.append(False)
                continue
            bare = re.sub(r'\[UNIMOD:\d+\]|-', '', seq)
            if len(bare) > MAX_PEPTIDE_LEN:
                prosit_seqs.append(seq)
                valid_mask.append(False)
            else:
                prosit_seqs.append(seq)
                valid_mask.append(True)

        best_psm_df = best_psm_df.copy()
        best_psm_df["PROSIT_Sequence"] = prosit_seqs
        best_psm_df["Valid"] = valid_mask

        valid_df = best_psm_df[best_psm_df["Valid"]].copy()
        n_skip = len(best_psm_df) - len(valid_df)
        print(f"  Valid for prediction: {len(valid_df)}/{len(best_psm_df)}"
              f" ({n_skip} skipped)")

        if len(valid_df) == 0:
            print("  No valid peptides. Exiting.")
            return

        pep_list = valid_df["Peptide"].tolist()
        charge_list = valid_df["Charge"].astype(int).tolist()
        prosit_list = valid_df["PROSIT_Sequence"].tolist()

        # Step 5: Koina
        print(f"\n[5/6] Calling Koina API ({MODEL_NAME}, NCE={NCE})...")
        predictions = predict_with_koina(prosit_list, charge_list,
                                          nce=NCE, batch_size=BATCH_SIZE)

        # Step 6: Write outputs
        print(f"\n[6/6] Writing Phase 1 outputs...")
        n_spectra = write_msp(msp_out, pep_list, charge_list,
                               prosit_list, predictions)

        # Build output CSV — deduplicate tryptic_df for merge so
        # one PSM row per unique peptide (protein_id kept for reference)
        tryptic_unique = tryptic_df.drop_duplicates(
            subset="tryptic_peptide", keep="first")
        out_df = valid_df[[
            "Peptide", "Charge", "Modified Peptide",
            "Assigned Modifications", "PROSIT_Sequence",
            "Spectrum", "Batch", "Hyperscore",
            "Global_Spectral_Count", "Is_OHW", "Round",
        ]].copy()
        out_df = out_df.rename(columns={"Peptide": "tryptic_peptide"})
        out_df = out_df.merge(tryptic_unique, on="tryptic_peptide",
                              how="left")
        out_df.to_csv(csv_p1, index=False)
        print(f"  Saved CSV -> {os.path.basename(csv_p1)}")

        n_ohw = out_df["Is_OHW"].sum()
        n_multi = (~out_df["Is_OHW"]).sum()
        n_r1 = (out_df["Round"] == 1).sum()
        n_r2 = (out_df["Round"] == 2).sum()

        print("\n" + "=" * 65)
        print(" Phase 1 complete!")
        print(f"   Input protein rows:     {len(master_df)}")
        print(f"   Unique input peptides:  {len(peptide_set)}")
        print(f"   Peptides in PSMs:       {len(global_counts)}")
        print(f"   Valid for prediction:   {len(valid_df)}")
        print(f"   One-hit-wonders:        {n_ohw}")
        print(f"   Multi-PSM:              {n_multi}")
        print(f"   Round 1 (best PSM):     {n_r1}")
        print(f"   Round 2 (best PSM):     {n_r2}")
        print(f"   Spectra in MSP:         {n_spectra}")
        print(f"   MSP file:  {msp_out}")
        print(f"   CSV file:  {csv_p1}")
        print("=" * 65)

    # =========================================================================
    # Phase 2
    # =========================================================================
    if run_p2:
        run_phase2(csv_p1, msp_out, csv_p2, plot_d, n_workers,
                   no_plots=args.no_plots, cand_csv=cand_csv)

    # =========================================================================
    # Bracketed annotation output (from existing/all-new Phase 2 CSV)
    # =========================================================================
    if run_p2 or args.annotate_only:
        if not os.path.exists(csv_p2):
            print(f"\nERROR: Cannot annotate because Phase 2 CSV is missing: "
                  f"{csv_p2}")
            print("Run Phase 2 first, or run without --annotate-only.")
            return

        print("\n" + "=" * 65)
        print(" Appending bracketed metrics to master CSV")
        print("=" * 65)
        import ast as _ast
        sa_df = pd.read_csv(csv_p2)

        # Build peptide -> metrics lookup from Phase 2 output
        lookup = {}
        for _, r in sa_df.iterrows():
            pep = r["tryptic_peptide"]
            lookup[pep] = (
                r["SA_degrees"],
                r["SA_normalized"],
                r.get("Match_coverage", None),
                r.get("Match_coverage_pct", None),
                r.get("Confidence", None),
            )

        # Read original master CSV to preserve all columns
        if run_p1:
            ann_master = master_df.copy()
        else:
            ann_master = pd.read_csv(INPUT_CSV)

        sa_deg_col, sa_norm_col = [], []
        match_cov_col, match_cov_pct_col, conf_col = [], [], []
        for _, row in ann_master.iterrows():
            try:
                peps = _ast.literal_eval(row["peptide_sequence"])
            except (ValueError, SyntaxError):
                peps = []
            if not isinstance(peps, list):
                peps = [peps]

            degs, norms = [], []
            covs, cov_pcts, confs = [], [], []
            for pep in peps:
                pep = str(pep).strip()
                if pep in lookup:
                    d, n, cov, cov_pct, conf = lookup[pep]
                    degs.append(round(d, 2) if pd.notna(d) else None)
                    norms.append(round(n, 4) if pd.notna(n) else None)
                    covs.append(str(cov) if pd.notna(cov) else None)
                    cov_pcts.append(round(cov_pct, 1)
                                    if pd.notna(cov_pct) else None)
                    confs.append(str(conf) if pd.notna(conf) else None)
                else:
                    degs.append(None)
                    norms.append(None)
                    covs.append(None)
                    cov_pcts.append(None)
                    confs.append(None)
            sa_deg_col.append(str(degs))
            sa_norm_col.append(str(norms))
            match_cov_col.append(str(covs))
            match_cov_pct_col.append(str(cov_pcts))
            conf_col.append(str(confs))

        ann_master["SA_degrees"] = sa_deg_col
        ann_master["SA_normalized"] = sa_norm_col
        ann_master["Match_coverage"] = match_cov_col
        ann_master["Match_coverage_pct"] = match_cov_pct_col
        ann_master["Confidence"] = conf_col
        ann_master.to_csv(csv_ann, index=False)

        n_with_sa = sum(1 for v in sa_deg_col if "None" not in v)
        print(f"  Saved -> {os.path.basename(csv_ann)}")
        print(f"  {len(ann_master)} protein rows "
              f"({n_with_sa} with full SA coverage)")
        print(f"  Metrics matched for {len(lookup)} unique peptides")


if __name__ == "__main__":
    main()
