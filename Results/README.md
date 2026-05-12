# Brain Microproteins Dashboard & Results

This directory contains:

1. The **interactive Streamlit dashboard**
   ([microproteins_dashboard.py](microproteins_dashboard.py)) used to explore
   the brain microprotein atlas and its differential expression in
   Alzheimer's disease.
2. The **summary CSVs** that the dashboard (and the manuscript) consumes,
   produced by `bash run_all_analyses.sh --mode=run` from the repo root.
3. The **figure libraries** referenced by the dashboard (mirror plots,
   expression-profile triptychs, smORF cartoons).
4. A **supplemental-table builder**
   ([generate_supplemental_tables.py](generate_supplemental_tables.py)) that
   consolidates everything into `supplementary/Supplemental_Tables.xlsx`.

## Live dashboard

A hosted version is available at
**<https://brain-microprotein-atlas.streamlit.app/>**.

## Local launch

```bash
# from repo root
cd Results
pip install -r requirements.txt
bash launch_dashboard.sh                # http://localhost:8505
```

`launch_dashboard.sh` activates the `github` conda env (edit if your env is
named differently), runs `streamlit run microproteins_dashboard.py
--server.port 8505 --server.headless true --server.runOnSave true`.

The **single-cell RNA enrichment** view is enabled by default.

## Figure assets: hosted vs. local

The dashboard renders three large figure libraries:

| Directory | Files | Approx. size |
|-----------|------:|-------------:|
| `mirror_plots/{Strong,Moderate,Weak,Insufficient}/` | 5,687 | 2.3 GB |
| `expression_profiles/{coupled,non_coupled}/`        | 8,652 | 1.9 GB |
| `smorf_cartoon_figures/`                            | 8,688 | 623 MB |

These are too large to ship in the GitHub repo, so they are **mirrored as a
public Hugging Face dataset**:
**[brmiller/brain-microprotein-atlas](https://huggingface.co/datasets/brmiller/brain-microprotein-atlas)**.

The dashboard is **local-first**: at startup it scans
`Results/{mirror_plots,expression_profiles,smorf_cartoon_figures}` and uses
on-disk files when present. If a directory is missing (e.g. on Streamlit
Cloud), it falls back to streaming images directly from the Hugging Face
dataset over HTTPS — no token required, no download step. The fallback is
indexed once via `huggingface_hub.HfApi().list_repo_files()` and cached for
the session, so per-image rendering only costs the actual file fetch.

To override the host (e.g. private mirror, S3, CloudFront) copy
`.streamlit/secrets.toml.example` to `.streamlit/secrets.toml` and set
`assets_base_url`. To force-download the assets locally, see
[../DATA_AVAILABILITY.md](../DATA_AVAILABILITY.md).

## Dashboard features

| Feature | Description |
|---------|-------------|
| Single-page multi-view UI | Sidebar selector switches between analyses; all views share the same filtered table. |
| **Annotation Summary** | Discovery metadata, smORF type, sequence, length, classification. |
| **Proteomics (TMT)** | Per-microprotein TMT evidence (PSMs, q-values, fold-changes) for Swiss-Prot vs. unreviewed. |
| **Proteomics + RiboSeq (RP3)** | Joint MS + RiboCode evidence view. |
| **Short-Read RNA in AD** | ROSMAP DLPFC + MSBB DESeq2 statistics (AD vs. control). |
| **Long-Read RNA in AD** | Nanopore ESPRESSO differential expression. |
| **scRNA Enrichment** | Cell-type-specific stats from Mathys et al. 2024. |
| **ShortStop Classification** | ML-based smORF classification results. |
| Filters | Gene/sequence search; smORF type (general + downstream sub-types); database (Reviewed vs Unreviewed); evidence type; PSM `Confidence` tier (Strong/Moderate/Weak/Insufficient/No PSM); statistical thresholds (TMT FDR, RNA-seq p-value); sequence completeness. |
| Mirror-plot gallery | PROSIT 3-panel diagnostics (sequence ladder + mirror spectrum + ppm-error lollipop) embedded inline per peptide. |
| Expression-profile viewer | PDF/PNG main-ORF / smORF triptychs keyed by genomic coordinates. |
| smORF cartoons | Per-locus cartoons indexed by `chrX_start-end`. |
| UCSC Genome Browser links | One-click jump to a custom UCSC session per microprotein. |
| Row-selection detail panel | ID-card + tabbed detail panels (Mirror Plots, Expression Profiles, Annotations). |
| CSV export | Download the currently filtered table. |
| Glassmorphism theme | Color-coded Swiss-Prot (`#74a2b7`) vs. unreviewed (`#ed8651`). |
| Optional password gate | Set `DASHBOARD_PASSWORD_HASH` env var to require login. |

## Files in this directory

### Scripts

| File | Purpose |
|------|---------|
| [microproteins_dashboard.py](microproteins_dashboard.py) | Main Streamlit app. |
| [launch_dashboard.sh](launch_dashboard.sh) | Local launcher (port 8505). |
| [generate_supplemental_tables.py](generate_supplemental_tables.py) | Builds `supplementary/Supplemental_Tables.xlsx` (formatted, one sheet per S-table). Run with `--include-scrna` to add the scRNA-seq tab. |
| [convert_expr_profiles.py](convert_expr_profiles.py) | Utility for converting expression-profile figures (placeholder). |
| [convert_hires.py](convert_hires.py) | Utility for high-resolution figure export (placeholder). |
| [requirements.txt](requirements.txt) | Dashboard-only Python deps (`streamlit`, `pandas`, `plotly`, `pathlib2`, `huggingface_hub`). |
| `microproteins_dashboard_backup.py`, `microproteins_dashboard_broken.py`, `microproteins_dashboard_pre_restore_*.py` | Historical snapshots — not used at runtime. |

### Data directories

| Directory | Contents |
|-----------|----------|
| `Annotations/` | `Brain_Microproteins_Discovery_summary.csv`, `ShortStop_Microproteins_summary.csv`, `smORF_type_definitions.csv`. |
| `Proteomics/` | `Proteomics_Results_summary.csv` (TMT evidence table). |
| `RP3/` | `RP3_Results_summary.csv`, `RP3_psORFs.csv`, RiboCode BED/GTF/TXT, `mapping_groups_rpkm*.txt`, `ribocode_results_ORFs_category.pdf`. |
| `ShortStop/` | `ShortStop_Microproteins_summary.csv`. |
| `Transcriptomics/` | `Short-Read_Transcriptomics_Results_summary.csv`, `Long-Read_Transcriptomics_Results_summary.csv`. |
| `scRNA_Enrichment/` | `scRNA_Enrichment_summary.csv`, UpSet plots + input matrices, cell-type heatmaps (`heatmap_PSM.pdf`, `heatmap_log2FC.pdf`, `cell_type_smorf_type_heatmap.pdf`), `volcano_all_celltypes.pdf`. |
| `expression_profiles/coupled/`, `expression_profiles/non_coupled/` | Per-pair triptych figures (PNG/PDF) routed by main-ORF/smORF coupling (`|Δr| > 0.1`). Files named `GENE_chrX_start-end.{png,pdf}`. |
| `mirror_plots/{Strong,Moderate,Weak,Insufficient}/` | PROSIT 3-panel mirror plots stratified by `Confidence` tier. |
| `smorf_cartoon_figures/` | Per-locus smORF cartoons keyed by genomic coordinate. |

### Datasets the dashboard reads

The app loads CSVs from each of the directories above plus
`../Code/data/cleaned_tryptic_peptides_under_151aa.csv` for the peptide
view. Files are merged by protein sequence; missing files degrade gracefully
(the corresponding view is hidden / disabled).

## Dependencies

```
streamlit >= 1.28.0
pandas    >= 1.5.0
plotly    >= 5.17.0
pathlib2  >= 2.3.7
huggingface_hub >= 0.24.0
```

Install with `pip install -r requirements.txt`. The dashboard does **not**
require the R / bioinformatics dependencies of the main analysis pipeline —
it only reads the summary CSVs.

## Regenerating the input CSVs

From the repository root:

```bash
bash run_all_analyses.sh --mode=run
```

This populates every subdirectory listed above from
`Code/data/microprotein_master.csv` and the per-module summary scripts. See
the top-level [README.md](../README.md) for details.
