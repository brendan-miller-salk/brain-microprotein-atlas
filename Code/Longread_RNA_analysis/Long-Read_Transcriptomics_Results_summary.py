import pandas as pd
import os
import sys

# === 1. Load and filter MASTER file using gold standard criteria ===
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from gold_standard_filtering_criteria import load_and_filter_master

mp = load_and_filter_master(os.path.join(os.path.dirname(__file__), '..', 'data', 'microprotein_master.csv'))

# === 2. Extract relevant columns (both Salk and Swiss-Prot-MP) ===
cols = [
    'sequence',
    'CLICK_UCSC',
    'nanopore_baseMean',
    'nanopore_log2FoldChange',
    'nanopore_pvalue',
    'nanopore_padj',
    'Database',
    'gene_name',
    'gene_symbol',
    'protein_length',
    'start_codon',
    'smorf_type',
    'total_razor_spectral_counts',
    'total_unique_spectral_counts'
]

summary_all = mp[cols].copy()

# Remove rows with missing values in 'nanopore_baseMean' for Swiss-Prot-MP entries
summary_all = summary_all[
    (summary_all['Database'] == 'Salk') |
    summary_all['nanopore_baseMean'].notna()
]

outdir = '/Users/brendanmiller/Library/CloudStorage/Box-Box/brendan_alan_shared_folder/AD_paper/Github/Results/Transcriptomics'

# === Optional: save combined output ===
summary_all.to_csv(os.path.join(outdir, 'Long-Read_Transcriptomics_Results_summary.csv'), index=False)
