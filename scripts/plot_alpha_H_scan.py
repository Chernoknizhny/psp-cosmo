import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from psp_model import H_psp
from stats import chi2, delta_chi2


def load_desi_data(path="data/desi_dr1.csv"):
    df = pd.read_csv(path)
    return df["z"].values, df["H"].values, df["sigma"].values


def scan_alpha_H(z, H_data, sigma, H0=67.4, Om=0.315, B=2.591, w=0.0367):
    """
    Сканирование по α_H: масштабируем H0 через α_H = H0/H0_fid.
    """
    alpha_H = np.linspace(0.8, 1.2, 200)
    chi2_vals = []

    for a in alpha_H:
        H_psp_vals = H_psp(z, H0 * a, Om, B, w)
        chi2_vals.append(chi2(H_psp_vals, H_data, sigma))

    chi2_vals = np.array(chi2_vals)
    chi2_min = chi2_vals.min()
    delta = chi2_vals - chi2_min

    return alpha_H, delta


def main():
    z, H_data, sigma = load_desi_data()
    alpha_H, delta = scan_alpha_H(z, H_data, sigma)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(alpha_H, delta, color="C0", lw=2)
    ax.axhline(0, color="k", ls="--", lw=0.8)
    ax.axhline(1, color="gray", ls=":", lw=0.8, label="1σ")
    ax.axhline(4, color="gray", ls=":", lw=0.8, label="2σ")
    ax.set_xlabel(r"$\alpha_H$")
    ax.set_ylabel(r"$\Delta\chi^2$")
    ax.set_title(r"Scan $\Delta\chi^2$ vs $\alpha_H$ (PSP, DESI DR1)")
    ax.legend()
    fig.tight_layout()

    out = os.path.join(os.path.dirname(__file__), '..', 'results', 'ah_scan_norms.png')
    fig.savefig(out, dpi=150)
    print(f"График сохранён: {out}")


if __name__ == "__main__":
    main()
