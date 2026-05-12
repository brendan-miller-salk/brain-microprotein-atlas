import pandas as pd
import os
import sys

# === 1. Load and filter MASTER file using gold standard criteria ===
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from gold_standard_filtering_criteria import load_and_filter_master

mp = load_and_filter_master(os.path.join(os.path.dirname(__file__), '..', 'data', 'microprotein_master.csv'))

# === 2. Extract relevant columns (both Salk and Swiss-Prot-MP) ===
summary_all = mp[[
    'sequence',
    'CLICK_UCSC',
    'gene_symbol',
    'genomic_coordinates',
    'rosmapRNA_baseMean',
    'rosmapRNA_log2FoldChange',
    'rosmapRNA_pvalue',
    'rosmapRNA_padj',
    'rosmapRNA_non_smorf_hit',
    'rosmapRNA_body',
    'ROSMAP_BulkRNAseq_CPM',
    'correlation_mainORF_nonAD_rosmap',
    'correlation_mainORF_AD_rosmap',
    'rosmap_lrt_additive_p',
    'rosmap_lrt_interaction_p',
    'msbbRNA_baseMean',
    'msbbRNA_log2FoldChange',
    'msbbRNA_pvalue',
    'msbbRNA_padj',
    'msbbRNA_non_smorf_hit',
    'msbbRNA_body',
    'MSBB_BulkRNAseq_CPM',
    'Database',
    'gene_name',
    'protein_length',
    'start_codon',
    'smorf_type',
    'total_razor_spectral_counts',
    'total_unique_spectral_counts',
    'mean_phylocsf'
]].copy()

# Swiss-Prot-MP filtering is handled centrally by the gold standard
# (load_and_filter_master): require Ribo-seq, MS, or DIA evidence.

# === Create output directories ===
outdir = '../../Results/Transcriptomics'
os.makedirs(outdir, exist_ok=True)

# === Save output files ===
summary_all.to_csv(os.path.join(outdir, 'Short-Read_Transcriptomics_Results_summary.csv'), index=False)

# === Preview ===
print("Short-Read Transcriptomics Results Summary:")
print(f"Unreviewed (Salk) microproteins: {(summary_all['Database'] == 'Salk').sum()}")
print(f"Swiss-Prot-MP microproteins: {(summary_all['Database'] == 'Swiss-Prot-MP').sum()}")
print(f"Total microproteins: {len(summary_all)}")
print("\nFirst few rows of combined data:")
print(summary_all.head())