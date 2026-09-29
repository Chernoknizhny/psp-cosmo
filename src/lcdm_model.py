import numpy as np


def H_lcdm(z, H0, Om):
    """
    Стандартная ΛCDM: H(z) = H0 * sqrt(Om*(1+z)^3 + (1-Om))
    w = -1, плоская Вселенная.
    """
    return H0 * np.sqrt(Om * (1 + z) ** 3 + (1 - Om))


def dH_lcdm(z, H0, Om):
    """
    Производная dH/dz для ΛCDM.
    """
    dz = 1e-6
    return (H_lcdm(z + dz, H0, Om) - H_lcdm(z, H0, Om)) / dz
