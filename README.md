# psp-cosmo

**Cosmology: PSP (Phase State Parameter) model vs ΛCDM.**  
Implementation, χ²/BIC tests, and out-of-sample predictions on Pantheon+, DESI, CMB data.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXX)

## Key idea

The goal of this repository is to provide a **reproducible, falsifiable comparison** between the PSP model and the standard ΛCDM cosmology.

- **PSP model:** Parameters $B = 2.591$ and $w = 0.0367$ are **fixed** from torus geometry (no fitting to data). Only $H_0$ and $\Omega_m$ are fitted.
- **Test strategy:** **Out-of-sample prediction**. The model is trained on Pantheon+ SNe data, and the $\chi^2$ is evaluated on independent datasets (DESI DR1, CMB) that were not used in the training.
- **Falsifiability:** Predictions for DESI Y5 (2027–2028) and CMB-S4 are explicitly stated in the documentation. If $w = -1$ is confirmed with precision better than 0.5%, the PSP model is falsified in this channel.

## Current Results

### 1. Out-of-sample test on DESI DR1 H(z)
The model was trained on Pantheon+ and tested on 5 DESI DR1 points without additional fitting.

| Model | $\chi^2$ (DESI DR1) | Notes |
|-------|---------------------|-------|
| **PSP** | **3.77** | Trained on Pantheon+, tested on unseen DESI data |
| ΛCDM | 8.56 | Standard model with $w=-1$ |
| **Difference** | **$\Delta\chi^2 = -4.79$** | PSP provides a significantly better fit to these out-of-sample data |

> **Note:** The improvement in $\chi^2$ is achieved **without** fitting $B$ or $w$. These parameters are fixed by the geometric derivation.

### 2. α_H Scan (SNe+BAO+CMB, zHD)
This repository reproduces the scan of the parameter $\alpha_H$ for different normalization factors (`norm`).

![Alpha_H Scan](docs/alpha_H_scan.png)

*Figure: $\Delta\chi^2$ vs $\alpha_H$ for three normalization levels (norm=0.5, 1.0, 2.5). Data: SNe+BAO+CMB, zHD.*

- **Key observation:** For `norm=0.5`, the curve shows a strong linear growth in $\Delta\chi^2$, indicating tension with the data at high $\alpha_H$.
- **Best fit:** The minimum $\Delta\chi^2$ occurs near $\alpha_H \approx 0.0$ for all normalization levels, consistent with the geometric derivation.

*The script to reproduce this plot is available in `scripts/plot_alpha_H_scan.py`.*

## How to Reproduce

### Prerequisites
- Python 3.8+
- DVC (for data management)

### Installation
```bash
git clone https://github.com/Chernoknizhny/psp-cosmo.git
cd psp-cosmo
python -m venv venv
# Windows: venv\Scripts\activate
# macOS/Linux: source venv/bin/activate
pip install -r requirements.txt
