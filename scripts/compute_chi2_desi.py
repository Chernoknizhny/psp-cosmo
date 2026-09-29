import numpy as np
import pandas as pd
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from psp_model import H_psp
from lcdm_model import H_lcdm
from stats import chi2, bic, delta_chi2


def load_desi_data(path="data/desi_dr1.csv"):
    """
    Загрузка данных DESI DR1: z, H(z), sigma
    """
    df = pd.read_csv(path)
    return df["z"].values, df["H"].values, df["sigma"].values


def main():
    z, H_data, sigma = load_desi_data()
    n = len(z)

    # Параметры (фиксированные, не фитились)
    H0 = 67.4
    Om = 0.315

    H_psp_vals = H_psp(z, H0, Om)
    H_lcdm_vals = H_lcdm(z, H0, Om)

    chi2_psp = chi2(H_psp_vals, H_data, sigma)
    chi2_lcdm = chi2(H_lcdm_vals, H_data, sigma)
    dchi2 = delta_chi2(chi2_psp, chi2_lcdm)

    bic_psp = bic(chi2_psp, n_params=2, n_data=n)
    bic_lcdm = bic(chi2_lcdm, n_params=2, n_data=n)

    print("=" * 40)
    print("Результаты DESI DR1")
    print("=" * 40)
    print(f"PSP   χ² = {chi2_psp:.2f}  BIC = {bic_psp:.2f}")
    print(f"ΛCDM  χ² = {chi2_lcdm:.2f}  BIC = {bic_lcdm:.2f}")
    print(f"Δχ² = {dchi2:.2f}  (PSP {'лучше' if dchi2 < 0 else 'хуже'})")
    print("=" * 40)


if __name__ == "__main__":
    main()
