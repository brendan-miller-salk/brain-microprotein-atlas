# Miscellaneous

Actin fluorescence-microscopy quantitation: the image-analysis pipeline behind
the F-actin distribution measurements in the manuscript.

This module is standalone. It is *not* part of
`bash run_all_analyses.sh --mode=run` — it does not read the master table and
produces none of the summary CSVs under `Results/`.

## Files

| File | Purpose |
|------|---------|
| `actin_quant_pipeline.py` | Quantifies subcellular F-actin distribution from multi-channel fluorescence images (ND2 or TIFF). |

## Usage

```bash
pip install numpy scikit-image scipy nd2 matplotlib tifffile

python actin_quant_pipeline.py \
    --image-folder  /path/to/images \
    --output-folder /path/to/output
```

`python actin_quant_pipeline.py --help` lists every parameter. Channel indices
default to nucleus = 0 and actin = 2; override them if your acquisition differs.

## Outputs

- `actin_data.csv` — per-bin intensity profiles and summary metrics per cell.
  Bin 0 and Bin N-1 are the cortical edges; central bins cover the nuclear
  region. Also reports `AverageNuclearActin`, `AverageCellActin`, and the
  derived `NucToCellRatio`.
- `diagnostic_plots/` — one overlay PNG per input image, for QC review.

## Inputs

The downsampled PNG exports the pipeline was run on are **not in a clone** —
52 MB of images, mirrored to Zenodo instead. From the repo root:

```bash
./download_zenodo_assets.sh actin_micrographs
```

That unpacks the 25 PNGs into this directory (`Code/Miscellaneous/`), which is
what `--image-folder` then wants. To place them somewhere else, download
`actin_micrographs.tar.gz` from
[Zenodo](https://doi.org/10.5281/zenodo.20045161) by hand and extract it
anywhere:

```bash
tar -xzf actin_micrographs.tar.gz -C /path/to/images
```

It holds 25 files: `ev_01.png` … `ev_13.png` (13 images, contiguous) and
`psorf2_01.png` … `psorf2_13.png` (12 images — the numbering skips
`psorf2_09`).

The raw ND2 acquisitions are not redistributed — `.nd2` and `.tif` are excluded
by `.gitignore`. See
[../../DATA_AVAILABILITY.md](../../DATA_AVAILABILITY.md).
