import pandas as pd
import os
import ast
import sys

# === 1. Load and filter MASTER file using gold standard criteria ===
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from gold_standard_filtering_criteria import load_and_filter_master

mp = load_and_filter_master(os.path.join(os.path.dirname(__file__), '..', 'data', 'microprotein_master.csv'))

# === 2. Keep only Salk (novel) microproteins ===
df = mp[mp['Database'] == 'Salk'].copy()

# === 3. Set annotation label based on which evidence is present ===
df['Annotation Status'] = df['has_MS'].map({True: 'MS', False: 'RiboCode_SAM'})

# Additional summary showing DDA vs DIA breakdown
print("\nBreakdown of MS evidence types:")
df['MS_type'] = df.apply(lambda row:
    'DDA+DIA' if row['DDA_evidence'] and row['DIA_evidence'] else
    'DDA only' if row['DDA_evidence'] else
    'DIA only' if row['DIA_evidence'] else
    'No MS', axis=1)

print(df[df['has_MS']]['MS_type'].value_counts())

# === 7. Derive DDA Grade from Confidence (pick best tier from bracketed list) ===
_confidence_rank = {'Strong': 0, 'Moderate': 1, 'Weak': 2, 'Insufficient': 3}

def best_confidence(val):
    if pd.isna(val):
        return None
    try:
        items = ast.literal_eval(val)
    except (ValueError, SyntaxError):
        return None
    if not isinstance(items, list):
        return None
    valid = [x for x in items if x in _confidence_rank]
    if not valid:
        return None
    return min(valid, key=lambda x: _confidence_rank[x])

df['DDA Grade'] = df['Confidence'].apply(best_confidence)

# Fill blank DDA Grade: RiboCode_SAM → "No MS", MS without grade → "No PROSIT"
df.loc[df['DDA Grade'].isna() & (df['Annotation Status'] == 'RiboCode_SAM'), 'DDA Grade'] = 'No MS'
df.loc[df['DDA Grade'].isna() & (df['Annotation Status'] == 'MS'), 'DDA Grade'] = 'No PROSIT'

# 8) select & rename summary columns
summary = df[[
    'sequence',
    'CLICK_UCSC',
    'genomic_coordinates',
    'smorf_type',
    'gene_name',
    'protein_length',
    'Annotation Status',
    'MS_type',
    'DDA Grade',
    'mean_phylocsf'
]].rename(columns={
    'genomic_coordinates': 'smORF Coordinates',
    'smorf_type':           'smORF Class',
    'gene_name':          'Parent Gene',
    'protein_length':       'Microprotein Length',
    'sequence':             'Microprotein Sequence',
    'MS_type':              'MS Evidence Type',
    'mean_phylocsf':       'PhyloCSF Score'
})

# Show sequences that start with M first
summary['start_codon'] = summary['Microprotein Sequence'].str.startswith('M')
summary = summary.sort_values(by='start_codon', ascending=False)
summary.drop(columns='start_codon', inplace=True)

# 7) write out
outdir = '../../Results/Annotations'
os.makedirs(outdir, exist_ok=True)
summary.to_csv(os.path.join(outdir, 'Brain_Microproteins_Discovery_summary.csv'),
               index=False)

# GitHub version - supplementary output commented out for relative path compatibility
# outdir_supplement = '/Users/brendanmiller/Library/CloudStorage/Box-Box/brendan_alan_shared_folder/AD_paper/supplementary'
# summary.to_csv(os.path.join(outdir_supplement, 'Table1_Brain_Microproteins_Discovery_summary.csv'),
#                index=False)

# preview
print(summary.head())

# Print dimensions
print(f"Number of rows: {summary.shape[0]}")

# print smORF Class 
print(f"Unique smORF Classes: {summary['smORF Class'].unique()}")

print("\nCounts of Annotation Status:")
print(df['Annotation Status'].value_counts())

