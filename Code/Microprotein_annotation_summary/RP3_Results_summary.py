import numpy as np
import pandas as pd
import os
import sys

# === 1. Load and filter MASTER file using gold standard criteria ===
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from gold_standard_filtering_criteria import load_and_filter_master

mp = load_and_filter_master(os.path.join(os.path.dirname(__file__), '..', 'data', 'microprotein_master.csv'))

# === 2. Keep Salk (novel) plus Reviewed (Swiss-Prot-MP) microproteins ===
df = mp[mp['Database'].isin(['Salk', 'Swiss-Prot-MP'])].copy()

# === Extract relevant columns ===
summary_salk = df[[
    'sequence',
    'CLICK_UCSC',
    'RP3_Default',
    'RP3_MM_Amb',
    'RP3_Amb',
    'RP3_MM',
    'RiboCode',
    'Database',
    'gene_symbol',
    'protein_length',
    'start_codon',
    'smorf_type',
    'mean_phylocsf'
]]

outdir = '/Users/brendanmiller/Library/CloudStorage/Box-Box/brendan_alan_shared_folder/AD_paper/Github/Results/RP3'

# === Optional: save combined output ===
summary_salk.to_csv(os.path.join(outdir, 'RP3_Results_summary.csv'), index=False)

# === Preview ===
print(summary_salk.head())

outdir_supplement = '/Users/brendanmiller/Library/CloudStorage/Box-Box/brendan_alan_shared_folder/AD_paper/supplementary'
summary_salk.to_csv(os.path.join(outdir_supplement, 'Table5_RP3_Results_summary.csv'), index=False)