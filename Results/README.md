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

Hosted versions (both online 24/7, same app):

- Streamlit Community Cloud: **<https://brain-microprotein-atlas.streamlit.app/>**
- Hugging Face Space (public, no password): **<https://huggingface.co/spaces/brmiller/brain-microprotein-atlas-app>**

The Hugging Face Space is deployed with
[push_to_hf_space.py](push_to_hf_space.py), which uploads the git-tracked file
set and streams figure images from the
[brmiller/brain-microprotein-atlas](https://huggingface.co/datasets/brmiller/brain-microprotein-atlas)
dataset. It runs password-free via the `DASHBOARD_PUBLIC=1` Space variable.

## Local launch

```bash
# from repo root
cd Results
pip install -r requirements.txt
bash launch_dashboard.sh                # http://localhost:8505
```

`launch_dashboard.sh` runs `streamlit run microproteins_dashboard.py
--server.port 8505 --server.headless true --server.runOnSave true`. It uses
whatever Python is on your `PATH`; set `DASHBOARD_CONDA_ENV=<env>` first if you
want it to activate a conda environment.

The **single-cell RNA enrichment** view is enabled by default.

## Figure assets: hosted vs. local

The dashboard renders three large figure libraries:

| Directory | Files | Approx. size |
|-----------|------:|-------------:|
| `mirror_plots/{Strong,Moderate,Weak,Insufficient}/` | 11,374 | 2.6 GB |
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
| Filters | Gene/sequence search above the table, plus cross-filtering sidebar facets in this top-to-bottom order: **Status** (Reviewed Swiss-Prot vs. Unreviewed Salk/TrEMBL) → **smORF Type** (general + nested Downstream sub-types) → **Evidence & Quality** (annotation method; PSM `Confidence` tier — Strong/Moderate/Weak/Insufficient/No PROSIT) → **Differential Expression** (TMT-MS FDR tiers 1–4, ROSMAP RNA FDR tiers) → **ShortStop Label** → **Score & Length Ranges** (protein length, unique spectral counts, PhyloCSF, UniProt annotation score) → **Start Codon** (ATG vs. nonATG) → **ORF Rules** (NMD/NSD × main-ORF disruption — see below) → **Kozak Context** (weak/adequate/strong; Salk smORFs only) → **N-terminus Options** (two checkbox options — see below). A final **Group by** section reorders rather than filters. Every facet is cross-filtered: its counts reflect all *other* active filters, ticked boxes within one facet are OR'd, and facets are AND'd together. |
| N-terminus options | A tryptic peptide starting at aa 1–2 is what Met excision of the parent ORF would produce, so it is not on its own evidence for the microprotein. *Non-N-terminal peptides only* keeps microproteins with a peptide starting at aa ≥ 3; *Nt-acetylated, substitution-distinct, or non-N-terminal* additionally readmits two other kinds of row — those carrying an Nt-acetylated peptide (Nt-acetylation is co-translational and marks a genuine N-terminus), and those whose N-terminal peptide carries 1–2 amino-acid substitutions vs. its matched UniProt isoform (`Nterm_Substitution_Rescue`), which a bare Met-excision fragment could not. Microproteins with no tryptic peptides at all fail both modes — absence of peptide data is not evidence. Ticked boxes are OR'd like every other facet, and the first option is a subset of the second, so ticking both equals ticking the second. Row filters only — peptide lists and spectral counts are unchanged. |
| N-terminal acetylation | `Nt_acetyl_*` columns from the master surface in the Proteomics column view (peptide, total PSMs, best PSM fraction) and as a per-peptide table in the detail card; Nt-acetylation marks a genuine protein N-terminus. |
| ORF Rules | Predicted decay fate of the host transcript × whether the smORF shifts the main ORF: 🟢 no NMD/mORF shift, 🟩 no NMD/mORF intact, 🟡 NMD/mORF shift, 🔴 NMD/mORF intact, 🟣 NSD/mORF shift, 🟪 NSD/mORF intact, ⚪ NA. NMD = the stop trips the 50-nt rule. NSD = non-stop decay, an ORF with no in-frame stop; unlike the NMD calls it is decided **per smORF over the whole transcript set** — it needs at least one confident no-stop transcript and no transcript carrying a real stop — and it overrides the NMD × disruption call, so it is a third decay state rather than a fifth colour. Read from `Code/data/smorf_trx_priority_color_class_assignments.tsv` and joined 1:1 on `gene_id`; the dots match the smORF Rules genome-browser track so the two read alike. Each NMD call is the class of one *representative* compatible transcript — a smORF whose transcripts disagree still shows a single colour — and all classes are predictions from transcript architecture, not measured decay. NA = not assessed: Swiss-Prot microproteins were never covered upstream, and alternative proteoforms are deliberately left unmatched. |
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
| [launch_hf_dashboard.sh](launch_hf_dashboard.sh) | Deploys the dashboard to both the Hugging Face Space and GitHub/Streamlit Cloud. Needs `HF_TOKEN` or a prior `hf auth login`. |
| [push_to_hf_space.py](push_to_hf_space.py) | Uploads the app to the Hugging Face Space (called by `launch_hf_dashboard.sh`). |
| [upload_to_hf.py](upload_to_hf.py) | Uploads the figure libraries to the Hugging Face dataset. |
| [requirements.txt](requirements.txt) | Dashboard-only Python deps (`streamlit`, `pandas`, `plotly`, `pathlib2`, `huggingface_hub`, `hf_xet`, `pyarrow`). |

### Data directories

| Directory | Contents |
|-----------|----------|
| `Annotations/` | `Brain_Microproteins_Discovery_summary.csv`, `ShortStop_Microproteins_summary.csv`, `smORF_type_definitions.csv`. |
| `Codon_context/` | `kozak/kozak_context.csv` + figures and `initiation/initiation_summary.csv`, `initiation_candidates.csv`, with per-run `report.md`, `COLUMNS.md`, and `qc_audit.csv`. |
| `Proteomics/` | `Proteomics_Results_summary.csv` (TMT evidence table). |
| `RP3/` | `RP3_Results_summary.csv`, `RP3_psORFs.csv`, RiboCode BED/GTF/TXT, `mapping_groups_rpkm*.txt`, `ribocode_results_ORFs_category.pdf`. |
| `ShortStop/` | `ShortStop_Microproteins_summary.csv`. |
| `Transcriptomics/` | `Short-Read_Transcriptomics_Results_summary.csv`, `Long-Read_Transcriptomics_Results_summary.csv`. |
| `scRNA_Enrichment/` | `scRNA_Enrichment_summary.csv` (significant pairs only), `scRNA_Enrichment_all_celltypes.csv.gz` (every tested microprotein x cell-type pair, read by the ID card's cell-type panel), UpSet plots + input matrices, cell-type heatmaps (`heatmap_PSM.pdf`, `heatmap_log2FC.pdf`, `cell_type_smorf_type_heatmap.pdf`), `volcano_all_celltypes.pdf`. |
| `expression_profiles/coupled/`, `expression_profiles/non_coupled/` | Per-pair triptych figures (PNG/PDF) routed by main-ORF/smORF coupling (`|Δr| > 0.1`). Files named `GENE_chrX_start-end.{png,pdf}`. |
| `mirror_plots/{Strong,Moderate,Weak,Insufficient}/` | PROSIT 3-panel mirror plots stratified by `Confidence` tier. |
| `smorf_cartoon_figures/` | Per-locus smORF cartoons keyed by genomic coordinate. |

Not everything above ships in a clone. The summary CSVs do, but the three figure
libraries and the large RiboCode `ribocode_results*.{bed,gtf,txt}` files in `RP3/`
are excluded by `.gitignore` — the figures stream from Hugging Face at runtime,
and the RiboCode files are on Zenodo. See
[../DATA_AVAILABILITY.md](../DATA_AVAILABILITY.md).

### Datasets the dashboard reads

The app loads CSVs from each of the directories above, plus two files from
`../Code/data/`: `microprotein_master.csv` (or `microprotein_master.zip`), which
it filters with the same gold-standard criteria as the analysis pipeline, and
`uniprotkb_proteome_UP000005640_2026_07_13.tsv` for UniProt annotation scores.
Files are merged by protein sequence; missing files degrade gracefully (the
corresponding view is hidden / disabled).

## Dependencies

```
streamlit       >= 1.37.0
pandas          >= 1.5.0
plotly          >= 5.17.0
pathlib2        >= 2.3.7
huggingface_hub >= 0.32.0
hf_xet          >= 1.0.0
pyarrow         >= 14.0
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
