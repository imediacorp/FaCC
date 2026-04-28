"""
CMB low-ℓ oscillation detection for φ-perturbation analysis.

Fits dual log-periodic oscillations (forward φ and conjugate φ^{-1} modes)
to residuals of the Planck TT low-ℓ power spectrum after subtracting a
smooth polynomial ΛCDM baseline.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import curve_fit

_PHI = (1 + np.sqrt(5)) / 2
_LN_PHI = np.log(_PHI)
_PHI_CONJ = _PHI - 1
_LN_PHI_CONJ = np.log(_PHI_CONJ)


def load_cmb_data(data_path: str) -> tuple[NDArray, NDArray, NDArray]:
    """Load CMB low-ℓ data from a CSV file.

    Parameters
    ----------
    data_path : str
        Path to CSV with columns: ell, C_ell, sigma.

    Returns
    -------
    ell, cl_data, sigma_cl : arrays
    """
    data = np.loadtxt(data_path, skiprows=1, delimiter=",")
    return data[:, 0], data[:, 1], data[:, 2]


def preprocess_cmb_data(
    ell: NDArray, cl_data: NDArray, sigma_cl: NDArray
) -> tuple[NDArray, NDArray, NDArray]:
    """Replace zero/negative uncertainties and mask invalid entries.

    Zero uncertainties are replaced with 5% of |C_ℓ| + 10.0 to avoid
    divide-by-zero in the fit.
    """
    estimated_sigma = 0.05 * np.abs(cl_data) + 10.0
    sigma_cl = np.where(sigma_cl <= 0, estimated_sigma, sigma_cl)

    valid = np.isfinite(sigma_cl) & np.isfinite(cl_data) & (ell > 0) & (sigma_cl > 0)
    return ell[valid], cl_data[valid], sigma_cl[valid]


def compute_residuals(ell: NDArray, cl_data: NDArray, poly_degree: int | None = None) -> tuple[NDArray, NDArray]:
    """Subtract a polynomial baseline (in log-ℓ space) and return residuals.

    Parameters
    ----------
    ell : array
        Multipole moments.
    cl_data : array
        Power spectrum values.
    poly_degree : int, optional
        Polynomial degree. Defaults to min(6, n-1).

    Returns
    -------
    cl_lcdm : array
        Smooth ΛCDM baseline.
    residuals : array
        cl_data - cl_lcdm.
    """
    degree = poly_degree if poly_degree is not None else min(6, len(ell) - 1)
    coeffs = np.polyfit(np.log(ell), cl_data, degree)
    cl_lcdm = np.polyval(coeffs, np.log(ell))
    return cl_lcdm, cl_data - cl_lcdm


def dual_osc_model(
    ell: NDArray,
    amp_phi: float,
    phase_phi: float,
    amp_conj: float,
    phase_conj: float,
) -> NDArray:
    """Dual log-periodic oscillation model for CMB residuals.

    f(ℓ) = A_φ × cos(2π × log(ℓ) / ln(φ) + φ_0)
           + A_{φ^{-1}} × cos(2π × log(ℓ) / ln(φ^{-1}) + φ_1)
    """
    return amp_phi * np.cos(2 * np.pi * np.log(ell) / _LN_PHI + phase_phi) + \
           amp_conj * np.cos(2 * np.pi * np.log(ell) / _LN_PHI_CONJ + phase_conj)


def fit_dual_oscillations(
    ell: NDArray, residuals: NDArray, sigma_cl: NDArray
) -> tuple[NDArray, NDArray]:
    """Fit dual log-periodic oscillations to CMB residuals.

    Returns
    -------
    popt : array
        Best-fit parameters [amp_phi, phase_phi, amp_conj, phase_conj].
    pcov : array
        Parameter covariance matrix.
    """
    amp_guess = float(np.std(residuals))
    popt, pcov = curve_fit(
        dual_osc_model,
        ell,
        residuals,
        sigma=sigma_cl,
        p0=[amp_guess, 0.0, amp_guess / 2, 0.0],
        absolute_sigma=True,
        maxfev=5000,
    )
    return popt, pcov


def run_cmb_analysis(
    data_path: str = "real_cmb_lowl.csv",
    output_path: str = "cmb_osc_dual.png",
) -> dict:
    """Run full CMB low-ℓ oscillation analysis and save figure.

    Parameters
    ----------
    data_path : str
        Path to CMB low-ℓ CSV data file.
    output_path : str
        Where to save the residual plot.

    Returns
    -------
    dict with keys: ell, residuals, sigma_cl, popt, pcov
    """
    ell_raw, cl_raw, sigma_raw = load_cmb_data(data_path)
    ell, cl_data, sigma_cl = preprocess_cmb_data(ell_raw, cl_raw, sigma_raw)
    cl_lcdm, residuals = compute_residuals(ell, cl_data)
    popt, pcov = fit_dual_oscillations(ell, residuals, sigma_cl)

    amp_phi, phase_phi, amp_conj, phase_conj = popt
    print(f"Fitted amp_φ:    {amp_phi:.2e}, phase_φ:    {phase_phi:.2f}")
    print(f"Fitted amp_φ⁻¹:  {amp_conj:.2e}, phase_φ⁻¹:  {phase_conj:.2f}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.figure(figsize=(10, 6))
    plt.errorbar(ell, residuals, yerr=sigma_cl, fmt=".", label="Residuals", alpha=0.6)
    plt.plot(
        ell,
        dual_osc_model(ell, *popt),
        "r-",
        label="Dual Log-Periodic Fit",
        linewidth=2,
    )
    plt.xscale("log")
    plt.xlabel("ℓ")
    plt.ylabel("ΔC_ℓ")
    plt.title("CMB Low-ℓ Residuals — Dual φ/φ⁻¹ Mode Fit")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()

    print(f"Figure saved → {output_path}")
    return {
        "ell": ell,
        "residuals": residuals,
        "sigma_cl": sigma_cl,
        "popt": popt,
        "pcov": pcov,
    }


if __name__ == "__main__":
    run_cmb_analysis()
