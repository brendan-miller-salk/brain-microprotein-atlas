# Brain Microprotein Atlas

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.20045162-blue)](https://doi.org/10.5281/zenodo.20045162)
[![Streamlit](https://img.shields.io/badge/dashboard-live-FF4B4B?logo=streamlit&logoColor=white)](https://brain-microprotein-atlas.streamlit.app/)
[![HF Dataset](https://img.shields.io/badge/HF-figures-yellow)](https://huggingface.co/datasets/brmiller/brain-microprotein-atlas)
[![PRIDE](https://img.shields.io/badge/PRIDE-PXD071500-blue)](https://www.ebi.ac.uk/pride/archive/projects/PXD071500)

Sequences, abundance data, analysis code, and an interactive dashboard for
the manuscript, *"A microprotein atlas of the human frontal cortex in
Alzheimer's disease." (in revision)*

---

## Step 1 — Clone

```bash
git clone https://github.com/brendan-miller-salk/brain-microprotein-atlas.git
cd brain-microprotein-atlas
```

## Step 2A — Launch the dashboard locally

```bash
cd Results
pip install -r requirements.txt
bash launch_dashboard.sh           # opens at http://localhost:8505
```

Password: `revision`. Figure assets (~5 GB) stream on demand from the
[Hugging Face dataset](https://huggingface.co/datasets/brmiller/brain-microprotein-atlas);
no local download is required.

Or use the **hosted version**: <https://brain-microprotein-atlas.streamlit.app/>

## Step 2B — Reproduce the analyses

```bash
bash setup_environment.sh                # one-time: creates conda env + deps
bash run_all_analyses.sh --mode=run      # runs every analysis end-to-end
```

This executes every summary analysis against the master CSV in
`Code/data/microprotein_master.csv` and regenerates every results table
under `Results/`. Use `--mode=docs` instead to print the step list
without executing.

**[Zenodo · 10.5281/zenodo.20045162](https://doi.org/10.5281/zenodo.20045162)**
holds the bulk archive that does not fit on GitHub: the full
`microprotein_master.csv` and supplemental tables, raw RP3 / RiboCode
outputs, TMT proteomics evidence tables, processed RNA-seq count matrices,
GTF/BED/FASTA reference files, and the complete figure libraries
(mirror plots, expression-profile triptychs, smORF cartoons).
See [DATA_AVAILABILITY.md](DATA_AVAILABILITY.md) and
[DATA_REQUIREMENTS.md](DATA_REQUIREMENTS.md) for the full data manifest.

---

## Repository layout

```
brain-microprotein-atlas/
├── README.md                       <- you are here
├── DATA_AVAILABILITY.md            <- raw data sources and access
├── DATA_REQUIREMENTS.md            <- shipped input files
├── environment.yml                 <- conda env spec
├── requirements.txt                <- pip deps
├── setup_environment.sh            <- installs env + Python deps
├── run_all_analyses.sh             <- reproduces every analysis (Step 2B)
│
├── Code/                           <- analysis pipelines
│   ├── data/                              shared input CSVs (master tables, etc.)
│   ├── Microprotein_annotation_summary/   smORF discovery + annotation
│   ├── Peptide_TMT_analysis/              TMT-MS proteomics (ROSMAP / MSBB)
│   ├── RP3_analysis/                      Ribo-seq (RP3 + RiboCode)
│   ├── Shortread_RNA_analysis/            ROSMAP/MSBB short-read RNA-seq
│   ├── Longread_RNA_analysis/             ESPRESSO long-read isoforms
│   ├── scRNAseq_summary_merging_analysis/ Mathys 2024 cell-type enrichment
│   └── Miscellanous/                      figure helpers, microscopy, misc.
│
├── Results/                        <- outputs + dashboard
│   ├── microproteins_dashboard.py         Streamlit app (Step 2A)
│   ├── launch_dashboard.sh                local launcher (port 8505)
│   ├── requirements.txt                   dashboard-only Python deps
│   ├── upload_to_hf.py                    figure -> Hugging Face uploader
│   ├── generate_supplemental_tables.py    builds Supplemental_Tables.xlsx
│   │
│   ├── Annotations/                       discovery + ShortStop summaries
│   ├── Proteomics/                        TMT evidence table
│   ├── RP3/                               RiboCode results (BED/GTF/TXT)
│   ├── Transcriptomics/                   short- + long-read DE results
│   ├── scRNA_Enrichment/                  cell-type enrichment + heatmaps
│   ├── ShortStop/                         ML smORF classification results
│   ├── mirror_plots/                      PROSIT MS2 mirror plots (HF-hosted)
│   ├── expression_profiles/               main-ORF / smORF coupling triptychs (HF-hosted)
│   └── smorf_cartoon_figures/             per-locus cartoons (HF-hosted)
│
└── GTF_and_BED_files/              <- genome coordinate exports (BED, GTF, FASTA)
```

---

## Citation

Miller B. *et al.* "A microprotein atlas of the human frontal cortex in Alzheimer's disease" (manuscript in revision).

## Contact

Brendan Miller · Saghatelian Lab · Salk Institute · brmiller@salk.edu
