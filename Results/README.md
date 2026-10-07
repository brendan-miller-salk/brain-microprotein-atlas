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

## Transcript-level assignment rules

A PROSIT tier says how well a spectrum matches its peptide and to which smORF it
maps. It does not say **which proteoform** it came from, **where translation
started**, or **which transcriptoform** encodes it. Four rule sets help with
transcript-level assignment, each exposed as a sidebar facet and detailed in the
table below:

| Rules | Question | Most useful for | Detail |
|---|---|---|---|
| **Start-codon rules** | With no in-frame ATG, where *could* initiation occur? | non-ATG smORFs, whose annotated start is a calling artifact | [Codon_context](../Code/Codon_context/README.md) · [report](Codon_context/initiation/report.md) |
| **ORF rules** | Stop codon trips the 50-nt NMD rule? Host's main ORF disrupted? | `iORF`, `Iso` — smORFs inside a host gene, where displacing the main ORF is the real question | [Code/data](../Code/data/README.md) |
| **Peptide specificity** | Do the entry's MS peptides also occur (identical, I=L, same mass, or as a non-tryptic fragment) in another human UniProt protein? | TrEMBL and N-truncated `Iso` entries, whose peptides may come from the full-length protein; any unreviewed entry before follow-up | [Peptide_specificity](../Code/Peptide_specificity/README.md) |
| **N-terminus options** | Is the identifying peptide a genuine microprotein N-terminus, or just what Met excision of the host ORF would produce? | smORFs whose only peptide starts at aa 1–2 and so could be shared with the parent protein | see `N-terminus options` below |

These layers stack: an `iORF` that disrupts the host's main ORF, with high
transcript coverage and a strong PROSIT grade, for example.

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
| Filters | Gene/sequence search above the table, plus cross-filtering sidebar facets in this top-to-bottom order: **Status** (Reviewed Swiss-Prot vs. Unreviewed Salk/TrEMBL) → **smORF Type** (general + nested Downstream sub-types) → **Evidence & Quality** (annotation method; PSM `Confidence` tier — Strong/Moderate/Weak/Insufficient/No PROSIT; Ribosome Coverage (RiboCode-SAM, Coverage, Flanking); Flanking Ribo-seq — see below; Ribo-seq Density vs Host CDS — see below; Peptide Evidence; Fragment-Ion Coverage; Peptide Specificity — see below) → **Differential Expression** (TMT-MS FDR tiers 1–4, ROSMAP RNA FDR tiers) → **ShortStop Label** → **Score & Length Ranges** (protein length, unique spectral counts, PhyloCSF, UniProt annotation score) → **Start Codon** (ATG vs. nonATG) → **ORF Rules** (NMD/NSD × main-ORF disruption — see below) → **Kozak Context** (weak/adequate/strong; Salk smORFs only) → **N-terminus Options** (two checkbox options — see below). A final **Group by** section reorders rather than filters. Every facet is cross-filtered: its counts reflect all *other* active filters, ticked boxes within one facet are OR'd, and facets are AND'd together. |
| N-terminus options | A tryptic peptide starting at aa 1–2 is what Met excision of the parent ORF would produce, so it is not on its own evidence for the microprotein. *Non-N-terminal peptides only* keeps microproteins with a peptide starting at aa ≥ 3; *Nt-acetylated, substitution-distinct, or non-N-terminal* additionally readmits two other kinds of row — those carrying an Nt-acetylated peptide (Nt-acetylation is co-translational and marks a genuine N-terminus), and those whose N-terminal peptide carries 1–2 amino-acid substitutions vs. its matched UniProt isoform (`Nterm_Substitution_Rescue`), which a bare Met-excision fragment could not. Microproteins with no tryptic peptides at all fail both modes — absence of peptide data is not evidence. Ticked boxes are OR'd like every other facet, and the first option is a subset of the second, so ticking both equals ticking the second. Row filters only — peptide lists and spectral counts are unchanged. |
| Peptide specificity | For the unreviewed (Salk + TrEMBL) entries: whether each MS peptide also occurs in another human UniProt protein (UP000005640 incl. TrEMBL and isoforms). Levels: **Unique**, **Deamidation-only match**, **Same-mass peptide elsewhere** (b/y ions may separate them), **Shared tryptic peptide**, **Possible fragment · same gene / · other protein** (the peptide sits inside that protein at a site trypsin does not cut, so the fully tryptic search never considered it there), **No MS peptides**, and NA for Swiss-Prot. It appears as a sidebar facet (the six peptide-match levels only; **No MS peptides** and NA have no box but still show in the column), a Proteomics-view column (`(partial)` = only some peptides match), and a per-peptide entry/twin table in the detail card. Non-unique means ambiguous by peptide specificity alone, not disproved. Read from `Code/data/unreviewed_microproteins_peptide_flags.csv`, joined 1:1 on `gene_id`. |
| Flanking Ribo-seq | For Salk and TrEMBL smORFs: adult-brain Ribo-seq P-sites within 50/100/250/1000 nt up- and downstream of the CDS, walked along the spliced host transcript. In the sidebar, two rows of boxes (**Upstream**, **Downstream**) keep smORFs with at least one P-site within the ticked window. Windows nest, so several ticks on one side equal the largest; the two sides are AND'd. The **Flanking** box under Ribosome Coverage is a shortcut for at least one P-site within 1 kb on either side. The Ribo-seq column view adds P-site counts per window and the frame-0 % at 250 nt. The entry card adds a per-window table (P-sites, nt available, density relative to the CDS, overlap with annotated CDS, frame 0/1/2 %) and a P-site plot coloured by frame relative to the smORF start. The plot is read live from the `adult_psite` bigWigs: local `Code/data/browser_tracks/` if present, else downloaded once from the Hugging Face dataset; it needs `pyBigWig`. Frame % excludes annotated-CDS positions and the 6 nt before the start codon. Read from `Code/data/smorf_flank_psites.csv` and `smorf_flank_structure.csv`, joined 1:1 on `gene_id`. |
| Ribo-seq Density vs Host CDS | For Salk and TrEMBL smORFs: Ribo-seq density on the smORF relative to the main CDS of its host gene, as a ratio of RPKMs (no pseudocount). uORFs and dORFs that do not overlap the host CDS are compared on all P-sites (all frames); smORFs that overlap or sit inside the host CDS on frame-0 P-sites (the P-site Frame 0 RPKM column), so the host's ribosomes on shared codons are not counted for the smORF. The host CDS is the smORF transcript's own annotated CDS when it has one, else the host gene's MANE Select (or Ensembl canonical) CDS. A pair counts only when the smORF and the host start codon lie on one transcript model (GENCODE, ENCODE4 or ESPRESSO long-read), the host has ≥ 20 P-sites, and the smORF is not an in-frame class; `Iso`/`N-Iso`/`D-Iso` and TrEMBL entries overlapping the host CDS are left out, since their P-sites are largely the host's. The sidebar has three nested boxes (**Denser than host (>1×)**, **≥2× host density**, **≥10× host density**) that keep smORFs with ≥ 10 P-sites whose density exceeds the host's by at least that factor; ticking several equals ticking the lowest. The Ribo-seq column view adds **Host CDS Frame 0 RPKM** (next to the smORF's own P-site Frame 0 RPKM), **smORF/Host Ratio (frame 0)**, **smORF/Host Ratio (all frames)**, **Host CDS Comparison** (the call: ORF denser, Host denser, Low ORF coverage, Low host coverage, Not on shared mRNA, In-frame), **Comparison Basis** (which ratio the call used) and **Host CDS** (symbol and transcript). Swiss-Prot, lncRNA, psORF and eORF entries are not compared, nor are smORFs without an identifiable host CDS; most TrEMBL entries overlap the host CDS or lack their own P-site counts, so only a handful are compared. Density is whole-ORF ribosome occupancy, not protein output, and start/stop pile-up on short ORFs is not removed. Read from `Code/data/smorf_host_cds_density.csv`, joined 1:1 on `gene_id`. |
| Parent Gene (TrEMBL) | For TrEMBL entries, Parent Gene is the gene at the CDS locus (GENCODE v50 host genes) rather than the UniProt gene name wherever the two differ: placeholder names (LOC..., the accession itself), names UniProt has since changed, and paralogues (e.g. a CDS on CALM2 named CALM1). Unnamed loci show their Ensembl gene ID. Search still matches the old name. Read from `Code/data/trembl_locus_genes.csv`, joined 1:1 on `gene_id`. |
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
| [upload_figures_to_hf.sh](upload_figures_to_hf.sh) | Wrapper for `upload_to_hf.py`: authenticates, checks the three figure directories, and picks a staging parent outside the Box-synced tree. Not part of the dashboard deploy. |
| [push_browser_tracks.py](push_browser_tracks.py) | Publishes the UCSC browser tracks from `Code/data/browser_tracks/` to the Hugging Face dataset, where `bigDataUrl` points. Uploads only changed files, then verifies each by hash. |
| [push_browser_tracks.sh](push_browser_tracks.sh) | Wrapper for `push_browser_tracks.py`: authenticates, reports the track inventory, then pushes. |
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
It also joins `unreviewed_microproteins_peptide_flags.csv` (Peptide Specificity),
`smorf_flank_psites.csv` (Flanking Ribo-seq),
`smorf_host_cds_density.csv` (Ribo-seq Density vs Host CDS),
`trembl_locus_genes.csv` (TrEMBL Parent Gene)
and the other `Code/data/` tables listed in [Code/data/README.md](../Code/data/README.md).
Files are merged by protein sequence or master `gene_id`; missing files degrade gracefully (the
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
