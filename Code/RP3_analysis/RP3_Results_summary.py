import numpy as np
import pandas as pd
import os
import sys

# === 1. Load and filter MASTER file using gold standard criteria ===
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from gold_standard_filtering_criteria import load_and_filter_master

mp = load_and_filter_master(os.path.join(os.path.dirname(__file__), '..', 'data', 'microprotein_master.csv'))

# === 2. Keep Salk (novel) plus Reviewed (Swiss-Prot-MP) microproteins ===
rp3_microproteins = mp[mp['Database'].isin(['Salk', 'Swiss-Prot-MP'])].copy()

print(f"Salk + Reviewed (Swiss-Prot-MP) microproteins with evidence: {len(rp3_microproteins):,}")

# === Extract relevant columns ===
summary_salk = rp3_microproteins[[
    'sequence',
    'CLICK_UCSC',
    'RP3_Default',
    'RP3_MM_Amb',
    'RP3_Amb',
    'RP3_MM',
    'RiboCode',
    'Database',
    'gene_name',
    'gene_symbol',
    'protein_length',
    'start_codon',
    'smorf_type',
    'total_unique_spectral_counts',
    'total_razor_spectral_counts',
    'mean_phylocsf'
]]

outdir = '../../Results/RP3'
os.makedirs(outdir, exist_ok=True)

# === Optional: save combined output ===
summary_salk.to_csv(os.path.join(outdir, 'RP3_Results_summary.csv'), index=False)

# === Preview ===
print(summary_salk.head())

# === Print Number of psORFs with  total_unique_spectral_counts > 0 and RP3_MM_Amb > 0 using smorf_type = "psORF"===
num_psorfs_rp3_mm_amb = summary_salk[(summary_salk['total_unique_spectral_counts'] > 0) & (summary_salk['smorf_type'] == 'psORF')].shape[0]
print(f"Number of psORFs with total_unique_spectral_counts > 0 and RP3_MM_Amb > 0: {num_psorfs_rp3_mm_amb}")

# == Save psORFs with RP3_MM_Amb > 0 to separate file ===
psorf_rp3_mm_amb = summary_salk[(summary_salk['total_unique_spectral_counts'] > 0) & (summary_salk['smorf_type'] == 'psORF')]
psorf_rp3_mm_amb.to_csv(os.path.join(outdir, 'RP3_psORFs.csv'), index=False)