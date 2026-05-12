import pandas as pd
import os
import sys

# === 1. Load and filter MASTER file using gold standard criteria ===
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from gold_standard_filtering_criteria import load_and_filter_master

mp = load_and_filter_master(os.path.join(os.path.dirname(__file__), '..', 'data', 'microprotein_master.csv'))

# === 2. Set annotation label for Salk entries ===
mp['Annotation Status'] = mp['has_MS'].map({True: 'MS', False: 'RiboCode_SAM'})

print(f"Counts by Annotation Status (Salk):")
print(mp.loc[mp['Database'] == 'Salk', 'Annotation Status'].value_counts())

# === 3. Extract relevant columns (both Salk and Swiss-Prot-MP) ===
cols = [
    'sequence',
    'CLICK_UCSC',
    'TMT_log2fc_50pct_missing',
    'TMT_pvalue_50pct_missing',
    'TMT_qvalue_50pct_missing',
    'TMT_log2fc_0pct_missing',
    'TMT_pvalue_0pct_missing',
    'TMT_qvalue_0pct_missing',
    'rate_control',
    'rate_ad',
    'Database',
    'protein_class_length',
    'gene_symbol',
    'gene_name',
    'protein_length',
    'start_codon',
    'smorf_type',
    'total_razor_spectral_counts',
    'total_unique_spectral_counts',
    'mean_phylocsf'
]

summary_all = mp[cols].copy()

# Swiss-Prot-MP filtering is handled centrally by the gold standard
# (load_and_filter_master): require Ribo-seq, MS, or DIA evidence.

print(f"Swiss-Prot-MP entries: {summary_all[summary_all['Database']=='Swiss-Prot-MP'].shape[0]}")
print(f"Total rows in output: {summary_all.shape[0]}")

outdir = '/Users/brendanmiller/Library/CloudStorage/Box-Box/brendan_alan_shared_folder/AD_paper/Github/Results/Proteomics'

# === Optional: save combined output ===
summary_all.to_csv(os.path.join(outdir, 'Proteomics_Results_summary.csv'), index=False)