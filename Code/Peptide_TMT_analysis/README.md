# TMT Proteomics Analysis

This module processes TMT (Tandem Mass Tag) proteomics data from FragPipe and performs quantitative analysis. Scripts in fragpipe_results_processing_scripts/ are the same onsed used to generate the summary data outputs here.

## Overview
- **Input**: FragPipe TMT results from multiple rounds
- **Output**: Quantified protein abundances and statistical comparisons
- **Main Script**: `Proteomics_Results_summary.py`

## Analysis Steps

### 1. Data Processing (`fragpipe_results_processing_scripts/`)

#### Protein ID Processing:
1. **`process_proteinID_from_TMT_round1.py`** - Process protein identifications from TMT round 1
2. **`process_proteinID_frmo_TMT_round2.py`** - Process protein identifications from TMT round 2
3. **`find_unique_tryptic_peptides.py`** - Identify unique tryptic peptides

#### Matrix Generation and Correction:
4. **`generate_raw_TMT_intensity_matrix_across_batches.R`** - Generate raw intensity matrices
5. **`TMT_regressed_corrrected_matrix.R`** - Apply batch correction and normalization

#### Statistical Analysis (TAMPOR):
6. **`TAMPOR_round1.R`** - TMT analysis for round 1
7. **`TAMPOR_round2.R`** - TMT analysis for round 2  
8. **`TAMPOR_combined_rounds.R`** - Combined analysis across rounds
9. **`TMT_ANOVA.R`** - ANOVA testing across conditions

### 2. PROSIT Spectral Validation (`prosit/`)

Generalized version of the PROSIT spectral-prediction and spectral-angle (SA)
validation pipeline used to vet PSMs for noncanonical microprotein peptides
across the TMT-ROSMAP FragPipe searches (round 1 batches `b1`-`b50` and
round 2 batches `round2/b1`-`round2/b14`). The original project-specific
pipeline lives at `brain_smorfs/prosit/` on the analysis HPC; the scripts
here are a pruned, parameterized re-distribution so the workflow is
reproducible from any FragPipe project.

#### Pipeline overview

For each tryptic peptide assigned to a microprotein:

1. **Phase 1 - prediction.** Scan FragPipe `peptide.tsv` / `psm.tsv` across
   every batch, collect up to 5 candidate PSMs per peptide
   (`Purity >= MIN_PURITY = 0.5`, ranked by Hyperscore, locked to the rank-1
   PSM's charge state and modified peptide). Convert FragPipe modifications
   to PROSIT-compatible UNIMOD form (`C[UNIMOD:4]`, `K[UNIMOD:737]`,
   `M[UNIMOD:35]`, N-term TMT `[UNIMOD:737]-`). PSMs whose only evidence is
   N-terminal acetylation (`n[43]`, unsupported by `Prosit_2020_intensity_TMT`)
   are dropped; peptides with at least one non-acetylated alternative are
   flagged `Acetyl_Rescued = True`. Predictions are obtained from the
   `Prosit_2020_intensity_TMT` model via the
   [Koina](https://koina.wilhelmlab.org) gRPC API (`NCE = 38`, `HCD`) and
   written as an MSP spectral library.
2. **Phase 2 - spectral angle.** For each candidate PSM, the corresponding
   MS2 scan is pulled from the source mzML with `pyteomics`, the TMT reporter
   window (125.5-131.5 m/z) is masked, predicted/observed peaks are matched
   at 20 ppm with each observed peak claimable by only one predicted ion, and
   predicted fragments with `fragment_charge >= precursor_charge` are dropped
   from both numerator and denominator. Match coverage is then computed as
   the fraction of unique backbone cleavage sites evidenced by matched
   ions - a `b_i` ion supports cleavage site `i` and a `y_i` ion supports
   site `peptide_len - i`, so duplicate ion-series hits at the same site are
   not double-counted (denominator = `peptide_len - 1`). The candidate with
   the most matched fragments (then evidence tier, then lowest SA) is
   selected as the representative spectrum and a 3-panel mirror plot
   (sequence ladder + mirror spectrum + ppm-error lollipop) is rendered.

#### Scoring rubric

| SA (degrees) | Rating   |
|--------------|----------|
| <= 20        | High     |
| 20-45        | Moderate |
| 45-70        | Poor     |
| > 70         | Terrible |

| Match coverage (% backbone cleavage sites evidenced) | Rating   |
|------------------------------------------------------|----------|
| <= 25%                                               | Sparse   |
| > 25 and <= 50%                                      | Partial  |
| > 50%                                                | Majority |

| SA \ Coverage | Sparse       | Partial      | Majority     |
|---------------|--------------|--------------|--------------|
| High          | Moderate     | **Strong**   | **Strong**   |
| Moderate      | Weak         | Moderate     | Moderate     |
| Poor          | Insufficient | Insufficient | Weak         |
| Terrible      | Insufficient | Insufficient | Insufficient |

The composite tier is written to the `Confidence` column of
`all_peptides_with_SA.csv`.

#### Files

- [prosit/prosit_pipeline.py](prosit/prosit_pipeline.py) - two-phase PROSIT
  prediction + SA pipeline (multiprocessing).
- [prosit/run_prosit_pipeline.sh](prosit/run_prosit_pipeline.sh) - SLURM
  wrapper that activates the `prosit` conda env and dispatches the pipeline.

#### PROSIT inputs

- A protein-level CSV with columns `protein_id` and a list-like
  `peptide_sequence` column (e.g. `['AAAPQAPAAR', 'RAAAPQAPAAR']`). The
  microprotein default is `cleaned_tryptic_peptides_detailed_under_151aa.csv`.
- FragPipe results laid out as
  `{batch}/shortstop_proteogenomics_appended_results_cpm05/DDA/{peptide,psm}.tsv`
  (the `RESULT_SUBDIR` constant in `prosit_pipeline.py` may be edited for
  other search configurations).
- Source `*.mzML` files in each batch directory, named to match the FragPipe
  `Spectrum` field.

#### PROSIT outputs

- `all_peptides_prosit.csv` - Phase 1, one row per peptide that produced a
  PROSIT prediction.
- `predicted_spectra/prosit_tmt_predicted_all.msp` - PROSIT MSP library.
- `candidate_psms.csv` - up to 5 candidate PSMs per peptide.
- `all_peptides_with_SA.csv` - Phase 2 augmentation (`SA_degrees`,
   `SA_normalized`, `SA_rating`, `Match_coverage`, `Match_coverage_pct`
   (% of unique backbone cleavage sites evidenced by matched b/y ions),
  `Confidence`, `Candidates_evaluated`, `Selected_Spectrum`,
  `Selected_Batch`, `Acetyl_Rescued`).
- `cleaned_tryptic_peptides_detailed_under_151aa_with_SA.csv` -
  bracket-aligned per-protein metrics merged back onto the input CSV.
- `mirror_plots/{Strong,Moderate,Weak,Insufficient}/*.{png,pdf}` -
  per-peptide 3-panel diagnostic figures.

### 3. Results Summary
- **`Proteomics_Results_summary.py`** - Generates final summary tables and statistics

## Usage

### Generate final summary:
```bash
cd ..
python Proteomics_Results_summary.py
```

### Run PROSIT validation:
```bash
cd prosit
python prosit_pipeline.py \
    --workspace  /path/to/fragpipe_root \
    --input-csv  cleaned_tryptic_peptides_detailed_under_151aa.csv \
    --output-dir prosit_out
```

## Dependencies
- Python: pandas, numpy, matplotlib, tqdm, pyteomics, koinapy, scipy
- R: Various statistical packages (dplyr, ggplot2, etc.)
- FragPipe output files (peptide.tsv, psm.tsv, source mzML)
- Network access to `koina.wilhelmlab.org:443` for PROSIT predictions

## Outputs
- Processed protein abundance matrices
- Statistical test results
- Microprotein quantification summaries
- Quality control metrics
- PROSIT MSP spectral library + per-peptide SA scores, Confidence tiers,
  and mirror-plot diagnostics
