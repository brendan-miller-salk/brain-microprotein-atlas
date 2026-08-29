#!/usr/bin/env python3
"""
Gold Standard Filtering Criteria for Brain Microprotein Discovery Pipeline
==========================================================================

Reference script documenting the canonical filtering strategy applied to the
MASTER dataset (2025-08-28_Brain_Microproteins_Discovery_MASTER_dataset.csv).

Use this as the authoritative template whenever loading and filtering the
master dataset for any downstream analysis (figures, tables, statistics, etc.).

Author: Brendan Miller
"""

import pandas as pd


def load_and_filter_master(csv_path: str) -> pd.DataFrame:
    """
    Load the MASTER dataset and apply the gold-standard filtering pipeline.

    Parameters
    ----------
    csv_path : str
        Path to the MASTER CSV file.

    Returns
    -------
    pd.DataFrame
        Filtered and deduplicated microprotein dataframe.
    """

    df = pd.read_csv(csv_path, low_memory=False)

    # ── Step 1: Remove contaminant / non-human entries ──────────────────────
    if "gene_symbol" in df.columns:
        pattern = r"SHEEP|Peptide|Horse|BOVIN"
        df = df[~df["gene_symbol"].str.contains(pattern, case=False, na=False)]

    # ── Step 2: Keep only relevant databases ────────────────────────────────
    df = df[df["Database"].str.contains("TrEMBL|Salk|Swiss-Prot", na=False)].copy()

    # ── Step 3: Remap Database labels ───────────────────────────────────────
    # TrEMBL entries are Salk-derived; Swiss-Prot → Swiss-Prot-MP
    df.loc[df["Database"] == "TrEMBL", "smorf_type"] = "TrEMBL"
    df["Database"] = df["Database"].replace({
        "TrEMBL":     "Salk",
        "Swiss-Prot": "Swiss-Prot-MP",
    })

    # ── Step 4: Build evidence flags ────────────────────────────────────────

    # RiboCode passes — robust to either boolean or string representation
    ribocode_pass = df["RiboCode"].astype(str).str.strip().str.upper() == "TRUE"

    # Strict RiboSAM: RiboCode + ShortStop SAM, no iORF, no iso-containing types
    df["has_RiboSAM"] = (
        ribocode_pass
        & df["shortstop_label"].isin(["SAM-Secreted", "SAM-Intracellular"])
        & (~df["smorf_type"].isin(["iORF"]))
        & (~df["smorf_type"].str.contains("iso", case=False, na=False))
    )

    # Loose RiboSAM: RiboCode + ShortStop SAM, exclude only exact 'Iso'
    df["has_LooseRiboSAM"] = (
        ribocode_pass
        & df["shortstop_label"].isin(["SAM-Secreted", "SAM-Intracellular"])
        & (~df["smorf_type"].isin(["Iso"]))
    )

    # DDA evidence: any unique spectral counts > 0
    df["DDA_evidence"] = df["total_unique_spectral_counts"] > 0

    # DIA evidence (novel ORFs): requires RiboSAM translation support + DIA detection
    # Used for Salk entries only — novel ORFs need both translation and proteomic evidence
    df["DIA_evidence"] = (
        df["has_LooseRiboSAM"]
        & df["Global.PG.Q.Value"].notna()
        & (df["Global.PG.Q.Value"] <= 0.01)
        & (df["Proteotypic"] == 1)
    )

    # DIA expression (known proteins): DIA detection alone, no translation gate
    # Used for Swiss-Prot entries — already curated, expression evidence is sufficient
    df["DIA_expression"] = (
        df["Global.PG.Q.Value"].notna()
        & (df["Global.PG.Q.Value"] <= 0.01)
        & (df["Proteotypic"] == 1)
    )

    # Combined MS flag
    df["has_MS"] = df["DDA_evidence"] | df["DIA_evidence"]

    # Coerce numeric-evidence columns that may contain stray non-numeric
    # tokens (e.g. "FALSE") in the upstream master
    for _col in ("RP3_Default", "total_razor_spectral_counts"):
        if _col in df.columns:
            df[_col] = pd.to_numeric(df[_col], errors="coerce").fillna(0)

    # ── Step 5: Filter to microproteins ─────────────────────────────────────
    df = df[df["protein_class_length"] == "Microprotein"].copy()

    # ── Step 6: Split by database, apply Salk-specific filtering, dedup ─────

    # Swiss-Prot-MP: require at least one line of evidence, then deduplicate by sequence
    #   - RP3_Default > 0         : Ribo-seq expression
    #   - total_razor_spectral_counts > 0 : any MS detection
    #   - DIA_expression           : DIA proteomics (no translation gate needed)
    mp_swiss = (
        df[df["Database"] == "Swiss-Prot-MP"]
        .loc[lambda d:
            (d["RP3_Default"] > 0)
            | (d["total_razor_spectral_counts"] > 0)
            | d["DIA_expression"]
        ]
        .drop_duplicates(subset="sequence")
    )

    # Salk (novel): must have at least one line of evidence, then
    # MS-prioritized deduplication (rows with MS evidence kept first)
    mp_salk = (
        df[(df["Database"] == "Salk")
           & (df["has_RiboSAM"] | df["DDA_evidence"] | df["DIA_evidence"])]
        .sort_values("has_MS", ascending=False)
        .drop_duplicates(subset="sequence")
    )

    mp = pd.concat([mp_swiss, mp_salk], ignore_index=True)

    return mp


# ── R equivalent (copy-paste into any R script) ────────────────────────────
R_TEMPLATE = r"""
library(dplyr)
library(stringr)

df <- read.csv(DATA_CSV, stringsAsFactors = FALSE) %>%
  filter(!grepl("SHEEP|Peptide|Horse|BOVIN", gene_symbol, ignore.case = TRUE)) %>%
  filter(str_detect(Database, "TrEMBL|Salk|Swiss-Prot")) %>%
  mutate(
    smorf_type = ifelse(Database == "TrEMBL", "TrEMBL", smorf_type),
    Database   = case_when(
      Database == "TrEMBL"     ~ "Salk",
      Database == "Swiss-Prot" ~ "Swiss-Prot-MP",
      TRUE                     ~ Database
    )
  )

df <- df %>%
  mutate(
    # RiboCode passes — robust to either logical TRUE or character "True"
    ribocode_pass = (RiboCode == TRUE) | (RiboCode == "True"),

    has_RiboSAM = ribocode_pass &
      (shortstop_label %in% c("SAM-Secreted", "SAM-Intracellular")) &
      (!smorf_type %in% c("iORF")) &
      (!grepl("iso", smorf_type, ignore.case = TRUE)),

    has_LooseRiboSAM = ribocode_pass &
      (shortstop_label %in% c("SAM-Secreted", "SAM-Intracellular")) &
      (!smorf_type %in% c("Iso")),

    DDA_evidence = total_unique_spectral_counts > 0,

    # DIA_evidence: for novel Salk ORFs — requires RiboSAM + DIA detection
    DIA_evidence = has_LooseRiboSAM &
      !is.na(Global.PG.Q.Value) & (Global.PG.Q.Value <= 0.01) &
      (Proteotypic == 1),

    # DIA_expression: for Swiss-Prot — DIA detection alone (no translation gate)
    DIA_expression = !is.na(Global.PG.Q.Value) & (Global.PG.Q.Value <= 0.01) &
      (Proteotypic == 1),

    has_MS = DDA_evidence | DIA_evidence
  )

# Swiss-Prot-MP: require at least one line of evidence (Ribo-seq, MS, or DIA)
mp_swiss <- df %>%
  filter(protein_class_length == "Microprotein", Database == "Swiss-Prot-MP") %>%
  filter(
    (!is.na(RP3_Default) & RP3_Default > 0) |
    (!is.na(total_razor_spectral_counts) & total_razor_spectral_counts > 0) |
    DIA_expression
  ) %>%
  filter(!duplicated(sequence))

mp_salk <- df %>%
  filter(protein_class_length == "Microprotein", Database == "Salk") %>%
  filter(has_RiboSAM | DDA_evidence | DIA_evidence) %>%
  arrange(desc(has_MS)) %>%
  filter(!duplicated(sequence))

mp <- bind_rows(mp_swiss, mp_salk)
"""


if __name__ == "__main__":
    # Example usage — update path as needed
    CSV_PATH = (
        "/Users/brendanmiller/Library/CloudStorage/Box-Box/brain_smorfs/annotation/"
        "2025-08-28_annotation_for_ribocode_plaqueomics_proteogenomics_DIA_shortstop/"
        "final_files/2025-08-28_Brain_Microproteins_Discovery_MASTER_dataset.csv"
    )

    mp = load_and_filter_master(CSV_PATH)

    print(f"Total microproteins after filtering: {len(mp)}")
    print(f"  Swiss-Prot-MP: {(mp['Database'] == 'Swiss-Prot-MP').sum()}")
    print(f"  Salk (novel):  {(mp['Database'] == 'Salk').sum()}")
    print(f"  With MS:       {mp['has_MS'].sum()}")
    print(f"  RNA only:      {(~mp['has_MS']).sum()}")
