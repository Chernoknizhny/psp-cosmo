
# psp-cosmo

**Cosmology: PSP (Phase State Parameter) model vs ΛCDM.**  
Implementation, χ²/BIC tests, and out-of-sample predictions on Pantheon+, DESI, CMB data.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXX)

## Key idea

The goal of this repository is to provide a **reproducible, falsifiable comparison** between the PSP model and the standard ΛCDM cosmology.

- **PSP model:** Parameters $B = 2.591$ and $w = 0.0367$ are **fixed** from torus geometry (no fitting to data). Only $H_0$ and $\Omega_m$ are fitted.
- **Test strategy:** **Out-of-sample prediction**. The model is trained on Pantheon+ SNe data, and the $\chi^2$ is evaluated on independent datasets (DESI DR1, CMB).
- **Falsifiability:** If $w = -1$ is confirmed with precision better than 0.5%, the PSP model is falsified in this channel.

## Current Results

### 1. Out-of-sample test on DESI DR1 H(z)

| Model | $\chi^2$ (DESI DR1) | Notes |
|-------|---------------------|-------|
| **PSP** | **3.77** | Trained on Pantheon+, tested on unseen DESI data |
| ΛCDM | 8.56 | Standard model with $w=-1$ |
| **Difference** | **$\Delta\chi^2 = -4.79$** | PSP provides a significantly better fit |

> **Note:** The improvement in $\chi^2$ is achieved **without** fitting $B$ or $w$.

### 2. α_H Scan (SNe+BAO+CMB, zHD)

![Alpha_H Scan](docs/alpha_H_scan.png)

*Figure: $\Delta\chi^2$ vs $\alpha_H$ for three normalization levels (norm=0.5, 1.0, 2.5).*

---

## How to Reproduce

### Installation
```bash
git clone https://github.com/Chernoknizhny/psp-cosmo.git
cd psp-cosmo
python -m venv venv
# Windows: venv\Scripts\activate
# macOS/Linux: source venv/bin/activate
pip install -r requirements.txt
Get Data
bash
dvc pull
Run Analyses
Compute χ²: python scripts/compute_chi2_desi.py
Reproduce Plot: python scripts/plot_alpha_H_scan.py
Full workflow: python scripts/retro_forecast.py
Structure
src/psp_model.py: PSP metric and H(z) calculation.
src/lcdm_model.py: ΛCDM H(z) for comparison.
src/stats.py: χ2, BIC, Δχ2 functions.
config/params.yaml: Fixed parameters (B,w) and priors.
scripts/: Ready-to-run scripts.
data/: Metadata; raw data via DVC.
results/: Output plots and tables.
Data Sources
Pantheon+: Scolnic et al. (2022).
DESI DR1: Aghamousa et al. (2022).
CMB: Planck Collaboration.
How to Cite
If you use this code or results, please cite:

This repository (Zenodo DOI).
Chernoknizhny, E. "PSP Metric as a Solution to Einstein's Equations." Zenodo. [DOI: 10.5281/zenodo.22952157]
Chernoknizhny, E. "Metric with Two Scale Factors." Zenodo. [DOI: 10.5281/zenodo.22946097]
License
MIT License.

text

6.  Прокрути вниз. В поле **Commit changes** напиши кратко: `update README with results and structure`.
7.  Нажми зелёную кнопку **Commit changes**.

✅ Готово. Теперь на главной странице репозитория есть нормальная аннотация с твоими числами.

---

## Шаг 2. Создаём `requirements.txt` (список библиотек)

Это нужно, чтобы любой человек мог поставить те же версии Python‑библиотек.

1.  На главной странице репозитория нажми **Add file** → **Create new file**.
2.  В поле имени файла напиши: `requirements.txt`.
3.  В большое окно редактора вставь этот текст:

```text
numpy
scipy
matplotlib
pandas
emcee
corner
dvc
Внизу нажми Commit changes, сообщение: add requirements.txt.
✅ Готово.

Шаг 3. Создаём config/params.yaml (чтобы доказать, что B и w не фитятся)
Это самый важный файл для защиты от критики «ты просто подогнал параметры».

Нажми Add file → Create new file.
В имени файла напиши: config/params.yaml.
(GitHub сам создаст папку config и положит туда файл).
В редактор вставь этот текст:
yaml
psp:
  B: 2.591
  w: 0.0367
  fixed: true
  fitted:
    - H0
    - Om

lcdm:
  w: -1.0
  fixed: true

priors:
  planck_Om:
    mean: 0.315
    sigma: 0.007
