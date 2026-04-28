"""
H(z) analysis for φ-recursive cosmology (background model, now falsified).

Provides functions to fit the Hubble parameter H(z) against the
φ-recursive model (forward and reverse branches) and ΛCDM baseline.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import minimize, OptimizeResult
from astropy.cosmology import FlatLambdaCDM

_PHI = (1 + np.sqrt(5)) / 2
_LN_PHI = np.log(_PHI)
_PHI_CONJ = _PHI - 1
_LN_PHI_CONJ = np.log(_PHI_CONJ)


def load_hz_data(data_path: str) -> tuple[NDArray, NDArray, NDArray]:
    """Load H(z) data from a CSV file.

    Parameters
    ----------
    data_path:
        Path to CSV with columns: z, H, sigma_H.

    Returns
    -------
    z_data, h_data, sigma_h : arrays
    """
    data = np.loadtxt(data_path, delimiter=",", skiprows=1)
    return data[:, 0], data[:, 1], data[:, 2]


def h_model(
    z: NDArray,
    om: float,
    t0: float,
    sigma: float = 1.0,
    h0: float = 70.0,
) -> NDArray:
    """φ-recursive H(z) model.

    Parameters
    ----------
    z : array
        Redshift values.
    om : float
        Matter density parameter Ω_m.
    t0 : float
        Characteristic time scale [s].
    sigma : float
        +1 for forward branch, -1 for reverse (conjugate φ^{-1}) branch.
    h0 : float
        Approximate Hubble constant [km/s/Mpc] used for unit conversion baseline.

    Returns
    -------
    H(z) in km/s/Mpc
    """
    ol = 1.0 - om
    ln_r = _LN_PHI if sigma > 0 else _LN_PHI_CONJ
    hubble = (sigma * ln_r / t0) * np.sqrt(om * (1 + z) ** 3 + ol)
    return hubble * 3.08568e19 / 3.15576e16


def fit_forward_branch(
    z_data: NDArray, h_data: NDArray, sigma_h: NDArray
) -> OptimizeResult:
    """Fit the forward (φ) branch to H(z) data."""

    def chi2(params: NDArray) -> float:
        om, t0 = params
        h_pred = h_model(z_data, om, t0, sigma=1.0)
        return float(np.sum(((h_data - h_pred) / sigma_h) ** 2))

    return minimize(chi2, [0.3, 1e-17], bounds=[(0.1, 0.4), (1e-18, 1e-16)])


def fit_reverse_branch(
    z_data: NDArray, h_data: NDArray, sigma_h: NDArray
) -> OptimizeResult:
    """Fit the reverse (φ^{-1}) branch to H(z) data."""

    def chi2(params: NDArray) -> float:
        om, t0 = params
        h_pred = h_model(z_data, om, t0, sigma=-1.0)
        return float(np.sum(((h_data - np.abs(h_pred)) / sigma_h) ** 2))

    return minimize(chi2, [0.3, 1e-17], bounds=[(0.1, 0.4), (1e-18, 1e-16)])


def run_hz_analysis(
    data_path: str = "real_hz.csv",
    output_path: str = "hz_comparison_dual.png",
) -> dict:
    """Run full H(z) analysis: load data, fit both branches, plot, return results.

    Parameters
    ----------
    data_path : str
        Path to the H(z) CSV data file.
    output_path : str
        Where to save the comparison figure.

    Returns
    -------
    dict with keys: z_data, h_data, sigma_h, res_forward, res_reverse, h_lcdm
    """
    z_data, h_data, sigma_h = load_hz_data(data_path)

    res_forward = fit_forward_branch(z_data, h_data, sigma_h)
    res_reverse = fit_reverse_branch(z_data, h_data, sigma_h)

    om_f, t0_f = res_forward.x
    om_r, t0_r = res_reverse.x

    cosmo_lcdm = FlatLambdaCDM(H0=70, Om0=0.3)
    h_lcdm = cosmo_lcdm.H(z_data).value

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.figure(figsize=(8, 5))
    plt.errorbar(z_data, h_data, yerr=sigma_h, fmt="o", label="Data")
    plt.plot(z_data, h_model(z_data, om_f, t0_f, sigma=1.0), label="Fibonacci Forward")
    plt.plot(
        z_data,
        np.abs(h_model(z_data, om_r, t0_r, sigma=-1.0)),
        "--",
        label="Fibonacci Reverse (abs)",
    )
    plt.plot(z_data, h_lcdm, ":", label="ΛCDM")
    plt.xlabel("Redshift z")
    plt.ylabel("H(z) [km/s/Mpc]")
    plt.title("H(z) Comparison: φ-Recursive vs ΛCDM")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()

    print(f"Forward  fitted Ωm={om_f:.3f}, t0={t0_f:.2e} s, χ²={res_forward.fun:.2f}")
    print(f"Reverse  fitted Ωm={om_r:.3f}, t0={t0_r:.2e} s, χ²={res_reverse.fun:.2f}")
    print(f"Figure saved → {output_path}")

    return {
        "z_data": z_data,
        "h_data": h_data,
        "sigma_h": sigma_h,
        "res_forward": res_forward,
        "res_reverse": res_reverse,
        "h_lcdm": h_lcdm,
    }


if __name__ == "__main__":
    run_hz_analysis()
