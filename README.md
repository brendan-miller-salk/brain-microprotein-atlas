# Brain Microprotein Atlas

[![HF Space](https://img.shields.io/badge/HF-dashboard-yellow?logo=huggingface)](https://huggingface.co/spaces/brmiller/brain-microprotein-atlas-app)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.20045161-blue)](https://doi.org/10.5281/zenodo.20045161)
[![HF Dataset](https://img.shields.io/badge/HF-figures-yellow)](https://huggingface.co/datasets/brmiller/brain-microprotein-atlas)
[![PRIDE](https://img.shields.io/badge/PRIDE-PXD071500-blue)](https://www.ebi.ac.uk/pride/archive/projects/PXD071500)

## Explore

**[Open the interactive atlas →](https://huggingface.co/spaces/brmiller/brain-microprotein-atlas-app)**

Search any microprotein and see its sequence, PROSIT spectra, and
AD differential expression.

> ### Build your own subset
> Stack any of the filters below, then **export the result as a filtered
> `.csv` or `.fasta`**. You are not restricted to the fixed confidence tiers
> further down this page — define the microprotein set your analysis actually
> needs and download exactly that.

Filters available:

| Filter | Narrows by |
|---|---|
| **Status** | Reviewed (Swiss-Prot) vs unreviewed (Salk/TrEMBL) |
| **smORF Type** | Top-level category — `iORF`, `Iso`, `uORF`, `dORF`, `lncRNA` … plus downstream sub-types |
| **Evidence & Quality** | TMT-MS ID, RiboCode-ShortStop, and PROSIT confidence tier |
| **Differential Expression** | TMT tiers 1–4 (FDR × sample coverage); ROSMAP RNA-seq FDR < 0.05 / < 0.2 |
| **ShortStop Label** | Ribosome-profiling classification |
| **Score & Length** | Sliders for protein length, unique spectral counts, PhyloCSF |
| **Start Codon** | ATG vs non-ATG (near-cognate / non-standard) |
| **ORF Rules** | Predicted decay fate (NMD / NSD / neither) × whether the smORF shifts the host's main ORF |
| **Kozak Context** | Sequence-context strength at the start codon (Non-UniProt smORFs only) |
| **N-terminus Options** | M-excision, acetylation |

Each microprotein also links straight into the **[AD Dark Microproteome UCSC session](https://genome.ucsc.edu/s/brmiller/AD%20Dark%20Microproteome)**
(hg38), which loads the Ribo-seq, long-read isoform, PROSIT peptide, and smORF
rules tracks at that locus.

*(The app sleeps after 48 h idle; the first load can take a minute to wake.)*

---

## Download

**GRCh38 / hg38**, UCSC-style `chr` names.

### Complete microprotein/smORF atlas regardless of confidence tiers (i.e., most permissive)

[![GTF](https://img.shields.io/badge/GTF-8250df?style=for-the-badge)](https://github.com/brendan-miller-salk/brain-microprotein-atlas/releases/download/reference-files/Unreviewed_Brain_Microproteins_Absent_from_UniProt.gtf)
[![BED](https://img.shields.io/badge/BED-8250df?style=for-the-badge)](https://github.com/brendan-miller-salk/brain-microprotein-atlas/releases/download/reference-files/Unreviewed_Brain_Microproteins_CDS_Absent_from_UniProt.bed)
[![FASTA](https://img.shields.io/badge/FASTA-8250df?style=for-the-badge)](https://github.com/brendan-miller-salk/brain-microprotein-atlas/releases/download/reference-files/Unreviewed_Brain_Microproteins.fasta)
[![IDs](https://img.shields.io/badge/IDs-8250df?style=for-the-badge)](https://github.com/brendan-miller-salk/brain-microprotein-atlas/releases/download/reference-files/Unreviewed_Brain_Microproteins_IDs.txt)

All evidence types in one bundle — **4,814 microprotein candidates, 5,681 FASTA records**.
Use this as the all-encompassing microprotein candidate list with further filters of interest (e.g., PROSIT, ORF rules, N-term rules)

### PROSIT-based-confidence subsets

Graded by best **PROSIT spectral-angle** peptide match in ROSMAP TMT-MS.
Ribo-ShortStop covers microproteins seen only by Ribo-seq.

| Subset | Microproteins | Records | Files |
|---|--:|--:|---|
| 🟢 **Strong** — strictest spectral angle, ≥25% ion coverage | 1,067 | 1,328 | [GTF](https://github.com/brendan-miller-salk/brain-microprotein-atlas/releases/download/reference-files/Unreviewed_Brain_Microproteins_Strong.gtf) · [BED](https://github.com/brendan-miller-salk/brain-microprotein-atlas/releases/download/reference-files/Unreviewed_Brain_Microproteins_CDS_Strong.bed) · [FASTA](https://github.com/brendan-miller-salk/brain-microprotein-atlas/releases/download/reference-files/Unreviewed_Brain_Microproteins_Strong.fasta) · [IDs](https://github.com/brendan-miller-salk/brain-microprotein-atlas/releases/download/reference-files/Unreviewed_Brain_Microproteins_IDs_Strong.txt) |
| 🟡 **Moderate** — moderate agreement | 1,433 | 1,847 | [GTF](https://github.com/brendan-miller-salk/brain-microprotein-atlas/releases/download/reference-files/Unreviewed_Brain_Microproteins_Moderate.gtf) · [BED](https://github.com/brendan-miller-salk/brain-microprotein-atlas/releases/download/reference-files/Unreviewed_Brain_Microproteins_CDS_Moderate.bed) · [FASTA](https://github.com/brendan-miller-salk/brain-microprotein-atlas/releases/download/reference-files/Unreviewed_Brain_Microproteins_Moderate.fasta) · [IDs](https://github.com/brendan-miller-salk/brain-microprotein-atlas/releases/download/reference-files/Unreviewed_Brain_Microproteins_IDs_Moderate.txt) |
| 🟠 **Weak** — low agreement, use with caution | 343 | 456 | [GTF](https://github.com/brendan-miller-salk/brain-microprotein-atlas/releases/download/reference-files/Unreviewed_Brain_Microproteins_Weak.gtf) · [BED](https://github.com/brendan-miller-salk/brain-microprotein-atlas/releases/download/reference-files/Unreviewed_Brain_Microproteins_CDS_Weak.bed) · [FASTA](https://github.com/brendan-miller-salk/brain-microprotein-atlas/releases/download/reference-files/Unreviewed_Brain_Microproteins_Weak.fasta) · [IDs](https://github.com/brendan-miller-salk/brain-microprotein-atlas/releases/download/reference-files/Unreviewed_Brain_Microproteins_IDs_Weak.txt) |
| ⚪ **Insufficient** — below reliable confidence | 158 | 203 | [GTF](https://github.com/brendan-miller-salk/brain-microprotein-atlas/releases/download/reference-files/Unreviewed_Brain_Microproteins_Insufficient.gtf) · [BED](https://github.com/brendan-miller-salk/brain-microprotein-atlas/releases/download/reference-files/Unreviewed_Brain_Microproteins_CDS_Insufficient.bed) · [FASTA](https://github.com/brendan-miller-salk/brain-microprotein-atlas/releases/download/reference-files/Unreviewed_Brain_Microproteins_Insufficient.fasta) · [IDs](https://github.com/brendan-miller-salk/brain-microprotein-atlas/releases/download/reference-files/Unreviewed_Brain_Microproteins_IDs_Insufficient.txt) |
| 🔵 **Ribo-ShortStop** — Ribo-seq translation, no DDA MS | 1,592 | 1,592 | [GTF](https://github.com/brendan-miller-salk/brain-microprotein-atlas/releases/download/reference-files/Unreviewed_Brain_Microproteins_Ribo_ShortStop.gtf) · [BED](https://github.com/brendan-miller-salk/brain-microprotein-atlas/releases/download/reference-files/Unreviewed_Brain_Microproteins_CDS_Ribo_ShortStop.bed) · [FASTA](https://github.com/brendan-miller-salk/brain-microprotein-atlas/releases/download/reference-files/Unreviewed_Brain_Microproteins_Ribo_ShortStop.fasta) · [IDs](https://github.com/brendan-miller-salk/brain-microprotein-atlas/releases/download/reference-files/Unreviewed_Brain_Microproteins_IDs_Ribo_ShortStop.txt) |


Full record counts, column definitions, and per-file details:
**[GTF_and_BED_files/README.md](GTF_and_BED_files/README.md)**

---

## Repository

```
├── Code/
│   ├── gold_standard_filtering_criteria.py  canonical master-table filter, imported everywhere
│   ├── data/                                shared input tables
│   ├── Microprotein_annotation_summary/     smORF discovery, classification, GTF/BED/FASTA export
│   ├── Codon_context/                       Kozak strength and non-ATG initiation
│   ├── RP3_analysis/                        Ribo-seq translation evidence (RiboCode)
│   ├── Peptide_TMT_analysis/                TMT-MS proteomics + PROSIT spectral validation
│   ├── Shortread_RNA_analysis/              ROSMAP/MSBB short-read DESeq2
│   ├── Longread_RNA_analysis/               Nanopore ESPRESSO isoforms
│   ├── scRNAseq_summary_merging_analysis/   Mathys 2024 cell-type enrichment
│   └── Miscellaneous/                       actin microscopy quantitation (standalone)
├── Results/               summary CSVs + the Streamlit dashboard
└── GTF_and_BED_files/     hg38 coordinate exports (GTF, BED, FASTA)
```

Module-by-module guide: **[Code/README.md](Code/README.md)** ·
Dashboard and outputs: **[Results/README.md](Results/README.md)**

### Run summary statistics yourself

```bash
git clone https://github.com/brendan-miller-salk/brain-microprotein-atlas.git
cd brain-microprotein-atlas

bash setup_environment.sh            # one-time: conda env + deps
bash run_all_analyses.sh --mode=run  # regenerate every summary table
```

Runs end to end on the shipped data; use `--mode=docs` to print the steps
without executing.

The bulk files that are too large for git — the combined Ensembl + microprotein
GTF, the full uncompressed master table, the figure libraries — are optional and
live on Zenodo:

```bash
./download_zenodo_assets.sh          # list what is available
./download_zenodo_assets.sh --all    # fetch and extract everything (~5.1 GB)
```

Each archive is checksum-verified against the Zenodo record before it is
extracted.

Upstream raw-data processing (RNA-seq, TMT-MS, Ribo-seq) is
documented but needs controlled-access data — see
[DATA_AVAILABILITY.md](DATA_AVAILABILITY.md) and
[DATA_REQUIREMENTS.md](DATA_REQUIREMENTS.md).

Run the dashboard locally:

```bash
cd Results && pip install -r requirements.txt && bash launch_dashboard.sh
```

---

## Citation

Miller B. *et al.* "A microprotein atlas of the human frontal cortex in
Alzheimer's disease." *Nature Aging* (2026).
doi:[10.1038/s43587-026-01207-x](https://doi.org/10.1038/s43587-026-01207-x)

Brendan Miller · Saghatelian Lab · Salk Institute · brmiller@salk.edu
