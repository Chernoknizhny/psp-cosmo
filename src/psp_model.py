import numpy as np


def H_psp(z, H0, Om, B=2.591, w=0.0367):
    """
    Модель PSP: H(z) = H0 * sqrt(Om*(1+z)^3 + (1-Om) * f(z, B, w))

    Параметры B и w зафиксированы и не фитятся.
    """
    E_de = Om * (1 + z) ** 3 + (1 - Om) * (1 + z) ** (3 * (1 + w)) * B ** ((1 + z) ** (-w) - 1)
    return H0 * np.sqrt(E_de)


def dH_psp(z, H0, Om, B=2.591, w=0.0367):
    """
    Производная dH/dz для PSP.
    """
    dz = 1e-6
    return (H_psp(z + dz, H0, Om, B, w) - H_psp(z, H0, Om, B, w)) / dz
