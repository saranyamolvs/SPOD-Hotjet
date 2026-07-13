# SPOD-Hotjet Shadowgraph Pipeline

This repository contains a batch pipeline for preprocessing impinging-jet shadowgraph TIFF images and computing spectral proper orthogonal decomposition (SPOD) for a Mach 1.5 nozzle data campaign.

## Expected data layout

Place the raw image campaign under an input root such as `RP`:

```text
RP/
  hD2/
    NPR2p5/*.tif
    NPR3p67/*.tif
    NPR5/*.tif        # `NPT5` is also accepted and normalized in outputs
  hD4/
  hD6/
  hD8/
  hD10/
  hD12/
  hD20/
```

The script discovers all `hD*` folders and all `NPR*`/`NPT*` subfolders automatically. Outputs are written to a separate directory so raw data are not modified.

## Install

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Run the full campaign

You must provide the camera frame spacing `dt` in seconds so SPOD frequencies are reported in Hz:

```bash
python scripts/shadowgraph_spod_pipeline.py --input-root RP --output-root SPOD_results --dt 1.0e-5
```

For a quick trial on a small subset:

```bash
python scripts/shadowgraph_spod_pipeline.py --input-root RP --output-root SPOD_trial --dt 1.0e-5 --hds 2 --nprs 2.5 --max-frames 300
```

## What is saved for each h/D and NPR case

Each case is written to `SPOD_results/<hD>/<NPR>/` with:

- `mean_raw.png`: time-averaged raw shadowgraph image.
- `flat_field_background.png`: low-pass background/illumination estimate used for correction.
- `mean_enhanced.png`: mean of contrast-enhanced frames.
- `rms_fluctuation.png`: RMS image highlighting unsteady shock/shear-layer/acoustic structures.
- `enhanced_frame_*.png`: representative enhanced frames for figure panels and quality control.
- `spod_eigenvalues.csv`: frequency-by-mode SPOD energy table for plotting and reporting.
- `spod_spectrum.png`: log-scale SPOD eigenvalue spectrum.
- `spod_results.npz`: compressed arrays containing frequencies, eigenvalues, real mode images, coefficients, and block starts.
- `top_frequency_*Hz/spod_mode_*.png`: mode-shape figures at the most energetic frequencies.
- `case_summary.json`: parameters and metadata for reproducibility.

A campaign-level `run_summary.json` is also saved in the output root.

## Preprocessing method

For each case, the pipeline:

1. Loads TIFF frames in natural filename order.
2. Converts color images to grayscale if needed.
3. Estimates a background/flat-field image by Gaussian filtering the raw temporal mean.
4. Subtracts and normalizes by the background level.
5. Applies percentile contrast stretching and gamma correction for clearer publication figures.
6. Forms zero-mean fluctuation fields for SPOD.

Important tunable options:

- `--clip-low` and `--clip-high`: percentile limits for contrast stretching. Start with `1` and `99`; use `0.5`/`99.5` if weak structures are clipped.
- `--flat-sigma`: Gaussian scale for background correction. Increase for very smooth illumination gradients; decrease if the background estimate removes large flow structures.
- `--gamma`: display gamma. Values below 1 brighten weak features.

## SPOD settings

The implementation uses a Welch/block SPOD formulation with a Hann window and overlap. Key options:

- `--nperseg`: frames per block. Larger values improve frequency resolution but reduce the number of independent blocks.
- `--overlap`: fractional overlap between blocks; `0.5` is a common starting point.
- `--max-modes`: number of SPOD modes saved at each frequency.
- `--save-top-frequencies`: number of frequencies with the largest first-mode energy exported as PNG mode figures.

For journal-quality comparisons, keep preprocessing and SPOD settings identical across all NPR and h/D cases, and report `dt`, `nperseg`, overlap, image resolution, number of frames, and any cropping or masking applied outside this script.
