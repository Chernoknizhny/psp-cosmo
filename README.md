
# psp-cosmo

**Cosmology: PSP (Phase State Parameter) model vs ΛCDM.**  
Implementation, χ²/BIC tests, and out-of-sample predictions on Pantheon+, DESI, CMB data.

## Key idea

- **PSP model:** parameters $B = 2.591$, $w = 0.0367$ are fixed from torus geometry (no fitting). Only $H_0$ and $\Omega_m$ are fitted.
- **Test:** out-of-sample prediction (e.g., Pantheon+ → DESI DR1) with explicit $\chi^2$ and $\Delta\text{BIC}$.
- **Goal:** reproducible, falsifiable comparisons with ΛCDM.

## Current results (example — replace with your latest numbers)

| Metric | PSP | ΛCDM | Notes |
|--------|-----|------|-------|
| $\chi^2$ (Pantheon+) | 1234.5 | 1236.7 | Example values |
| $\chi^2$ (DESI DR1, out-of-sample) | 3.75 | 8.56 | PSP trained on Pantheon+, tested on unseen DESI |
| $\Delta\text{BIC}$ | −5.2 | — | Negative ΔBIC indicates stronger support for PSP |

> **Important:** The PSP parameters $B$ and $w$ are not fitted to these data; they are derived from the torus geometry. This is a key difference from standard dynamical dark energy parametrizations.

## Quick start

```bash
git clone https://github.com/username/psp-cosmo.git
cd psp-cosmo
python -m venv venv
# Windows: venv\Scripts\activate
# macOS/Linux: source venv/bin/activate
pip install -r requirements.txt
Usage examples
Compute χ² on DESI DR1 (out-of-sample test):
bash
python scripts/compute_chi2_desi.py
Plot Δχ² vs α_H scan:
bash
python scripts/plot_alpha_H_scan.py
Run full retro-forecast workflow:
bash
python scripts/retro_forecast.py
Structure
src/ — core cosmology and PSP/ΛCDM models.
scripts/ — ready-to-run scripts (χ² calculation, plots, forecasts).
notebooks/ — exploratory analysis.
config/ — fixed parameters (
B
B, 
w
w) and priors.
data/ — metadata; large files via DVC.
results/ — plots, tables, logs.
Data
Large data files (Pantheon+ covariance matrix, DESI points) are managed with DVC.
After cloning:

bash
dvc pull
How to cite
If you use this code or results in your work, please cite:

The repository (Zenodo DOI will be available after linking with Zenodo).
The underlying model papers (Zenodo links):
License
MIT.

text

---

### Что важно сделать прямо сейчас

1. **Подставь реальные числа в таблицу «Current results».** Это самое сильное место для внешнего читателя. Если у тебя сейчас $\chi^2_{\text{PSP}} = 3.75$ на DESI DR1, а у ΛCDM = 8.56 — оставь эти цифры: они сразу показывают преимущество PSP в тесте на новых данных.
2. **Замени `username` в команде `git clone`** на свой реальный username на GitHub.
3. **В секции «How to cite» оставь ссылки на Zenodo** — это сразу связывает код с публикациями и снимает вопрос «откуда взялась модель».

---

### Про твой график (Δχ² vs α_H)

В секции «Usage examples» я добавил `plot_alpha_H_scan.py`. Если у тебя скрипт ещё не готов, я могу прямо сейчас написать его скелет: он будет брать CSV с колонками `alpha_H`, `norm`, `delta_chi2` и строить три линии (как на твоём изображении), сохраняя в `results/plots/`.

Скажи, какой вариант тебе важнее сейчас:
- «Дай скелет скрипта `plot_alpha_H_scan.py` под этот график»
- «Дай скелет `compute_chi2_desi.py` с формулой $\chi^2 = \sum \frac{(H_{\text{model}} - H_{\text{data}})^2}{\sigma^2}$»
- «Сначала закончим с README и структурой, код потом»

Что делаем?
