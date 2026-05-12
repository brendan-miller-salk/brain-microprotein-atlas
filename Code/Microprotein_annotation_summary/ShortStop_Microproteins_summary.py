import pandas as pd
import os
import sys

# === 1. Load and filter MASTER file using gold standard criteria ===
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from gold_standard_filtering_criteria import load_and_filter_master

mp = load_and_filter_master(os.path.join(os.path.dirname(__file__), '..', 'data', 'microprotein_master.csv'))

# === 2. Keep only Salk (novel) microproteins ===
df = mp[mp['Database'] == 'Salk'].copy()

# === 3. Set annotation label based on which evidence is present ===
df['Annotation Status'] = df['has_MS'].map({True: 'MS', False: 'RiboCode_SAM'})

# === Keep only rows with valid shortstop annotations ===
summary = df.loc[df['shortstop_label'].notna(), [
    'sequence',
    'CLICK_UCSC',
    'gene_symbol',
    'gene_name',
    'smorf_type',
    'shortstop_label',
    'shortstop_score',
    'Annotation Status',
    'mean_phylocsf'
]].rename(columns={
    'sequence': 'Microprotein Sequence',
    'gene_symbol': 'smORF ID',
    'gene_name': 'Gene Body (Name)',
    'smorf_type': 'smORF Class',
    'shortstop_label': 'ShortStop Label',
    'shortstop_score': 'ShortStop Score',
    'mean_phylocsf': 'PhyloCSF Score'
    
})

# === Print ShortStop label counts ===
print(summary['ShortStop Label'].value_counts())
print("\nCounts by Annotation Status:")
print(summary['Annotation Status'].value_counts())

# === Save summary ===
outdir = '/Users/brendanmiller/Library/CloudStorage/Box-Box/brendan_alan_shared_folder/AD_paper/Github/Results/Annotations'
summary.to_csv(os.path.join(outdir, 'ShortStop_Microproteins_summary.csv'), index=False)
