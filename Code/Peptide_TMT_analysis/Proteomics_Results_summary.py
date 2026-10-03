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
    'TMT_t_statistic_50pct_missing',
    'TMT_df_50pct_missing',
    'TMT_conf_low_50pct_missing',
    'TMT_conf_high_50pct_missing',
    'TMT_cohens_d_50pct_missing',
    'TMT_pvalue_50pct_missing',
    'TMT_qvalue_50pct_missing',
    'TMT_log2fc_0pct_missing',
    'TMT_t_statistic_0pct_missing',
    'TMT_df_0pct_missing',
    'TMT_conf_low_0pct_missing',
    'TMT_conf_high_0pct_missing',
    'TMT_cohens_d_0pct_missing',
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

# === 4. Add b/y fragment-ion coverage of each tryptic peptide ===
# Bracketed per-peptide lists aligned with peptide_sequence (theoretical b/y
# ladder of the PROSIT-selected PSM; prosit/prosit_pipeline.py --annotate-only),
# plus per-microprotein rollups over its peptides.
import ast

ladder_cols = ['ladder_b_coverage_pct', 'ladder_y_coverage_pct',
               'ladder_by_union_coverage_pct', 'ladder_longest_by_run',
               'ladder_consec5']
pep = pd.read_csv(os.path.join(os.path.dirname(__file__), '..', 'data',
                               'cleaned_tryptic_peptides_detailed_under_151aa_with_SA.csv'),
                  usecols=['sequence', 'peptide_sequence'] + ladder_cols)


def _vals(cell):
    if pd.isna(cell):
        return []
    return [v for v in ast.literal_eval(cell) if v is not None]


pep['max_ladder_by_union_coverage_pct'] = pep['ladder_by_union_coverage_pct'].map(
    lambda c: max(_vals(c), default=None))
pep['max_ladder_longest_by_run'] = pep['ladder_longest_by_run'].map(
    lambda c: max(_vals(c), default=None))
pep['n_peptides_ladder_consec5'] = pep['ladder_consec5'].map(
    lambda c: sum(_vals(c)) if _vals(c) else None)

summary_all = summary_all.merge(pep, on='sequence', how='left', validate='one_to_one')
for c in ['max_ladder_longest_by_run', 'n_peptides_ladder_consec5']:
    summary_all[c] = summary_all[c].astype('Int64')
print(f"Microproteins with b/y ladder coverage: "
      f"{summary_all['max_ladder_by_union_coverage_pct'].notna().sum()}")

# Swiss-Prot-MP filtering is handled centrally by the gold standard
# (load_and_filter_master): require Ribo-seq, MS, or DIA evidence.

print(f"Swiss-Prot-MP entries: {summary_all[summary_all['Database']=='Swiss-Prot-MP'].shape[0]}")
print(f"Total rows in output: {summary_all.shape[0]}")

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
outdir = os.path.join(REPO_ROOT, 'Results', 'Proteomics')

# === Optional: save combined output ===
summary_all.to_csv(os.path.join(outdir, 'Proteomics_Results_summary.csv'), index=False)