# psp-cosmo

**Cosmology: PSP (Phase State Parameter) model vs ΛCDM.**  
Implementation, χ²/BIC tests, and out-of-sample predictions on Pantheon+, DESI, CMB data.

## Key idea
- **PSP model:** parameters $B=2.591$, $w=0.0367$ are fixed from torus geometry (no fitting). Only $H_0$ and $\Omega_m$ are fitted.
- **Test:** out-of-sample prediction (e.g., Pantheon+ → DESI DR1) with explicit $\chi^2$ and $\Delta\text{BIC}$.
- **Goal:** reproducible, falsifiable comparisons with ΛCDM.

## Quick start
```bash
git clone https://github.com/username/psp-cosmo.git
cd psp-cosmo
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate
pip install -r requirements.txt
