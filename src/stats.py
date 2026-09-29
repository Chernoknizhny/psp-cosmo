import numpy as np


def chi2(H_model, H_data, sigma):
    """
    Расчёт χ² = Σ (H_model - H_data)² / σ²
    """
    return np.sum(((H_model - H_data) / sigma) ** 2)


def bic(chi2_val, n_params, n_data):
    """
    Bayesian Information Criterion: BIC = χ² + k·ln(n)
    """
    return chi2_val + n_params * np.log(n_data)


def delta_chi2(chi2_model, chi2_ref):
    """
    Разница χ²: Δχ² = χ²_model - χ²_ref
    """
    return chi2_model - chi2_ref
