import numpy as np
import emcee
import corner
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from psp_model import H_psp
from stats import chi2


def load_desi_data(path="data/desi_dr1.csv"):
    df = pd.read_csv(path)
    return df["z"].values, df["H"].values, df["sigma"].values


def log_prior(theta):
    H0, Om = theta
    if 50 < H0 < 80 and 0.1 < Om < 0.5:
        return 0.0
    return -np.inf


def log_likelihood(theta, z, H_data, sigma):
    H0, Om = theta
    H_model = H_psp(z, H0, Om)
    return -0.5 * chi2(H_model, H_data, sigma)


def log_posterior(theta, z, H_data, sigma):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, z, H_data, sigma)


def main():
    z, H_data, sigma = load_desi_data()

    ndim = 2
    nwalkers = 32
    nsteps = 5000
    burn = 1000

    p0 = np.random.uniform(low=[60, 0.25], high=[75, 0.40], size=(nwalkers, ndim))

    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=(z, H_data, sigma))
    print("Запуск MCMC...")
    sampler.run_mcmc(p0, nsteps, progress=True)

    samples = sampler.get_chain(discard=burn, flat=True)

    fig = corner.corner(
        samples,
        labels=[r"$H_0$", r"$\Omega_m$"],
        truths=[67.4, 0.315],
        show_titles=True,
        title_fmt=".3f",
    )
    out = os.path.join(os.path.dirname(__file__), '..', 'results', 'corner_psp.png')
    fig.savefig(out, dpi=150)
    print(f"Corner plot сохранён: {out}")

    H0_med, Om_med = np.median(samples, axis=0)
    H0_err = np.std(samples[:, 0])
    Om_err = np.std(samples[:, 1])
    print(f"H0 = {H0_med:.2f} ± {H0_err:.2f}")
    print(f"Om = {Om_med:.3f} ± {Om_err:.3f}")


if __name__ == "__main__":
    main()
