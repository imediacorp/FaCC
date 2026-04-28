"""
DESI DR2 BAO analysis for φ-oscillation search.

Fetches the public DESI DR2 (and DR1) BAO distance measurements from the
CobayaSampler/bao_data GitHub repository, computes ΛCDM predictions via
astropy, and searches the ΛCDM residuals for φ-log-periodic modulation in
redshift.

The φ-modulation model for BAO residuals:
    ΔR(z) = A_φ × cos(2π × ln(1+z) / ln(φ) + δ)

where R is any distance ratio (D_M/r_s, D_H/r_s, D_V/r_s).
"""

from __future__ import annotations

import io
import numpy as np
from numpy.typing import NDArray
from scipy.optimize import curve_fit

_PHI = (1 + np.sqrt(5)) / 2
_LN_PHI = np.log(_PHI)

# ── Data URLs ─────────────────────────────────────────────────────────────────

_BASE = "https://raw.githubusercontent.com/CobayaSampler/bao_data/master"

URLS = {
    "dr2_mean": f"{_BASE}/desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_mean.txt",
    "dr2_cov":  f"{_BASE}/desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_cov.txt",
    "dr1_mean": f"{_BASE}/desi_2024_gaussian_bao_ALL_GCcomb_mean.txt",
    "dr1_cov":  f"{_BASE}/desi_2024_gaussian_bao_ALL_GCcomb_cov.txt",
}

# Tracer labels for individual files
_DR2_TRACERS = {
    "BGS z≈0.3":   "desi_bao_dr2/desi_gaussian_bao_BGS_BRIGHT-21.35_GCcomb_mean.txt",
    "LRG z=0.4-0.6": "desi_bao_dr2/desi_gaussian_bao_LRG_GCcomb_z0.4-0.6_mean.txt",
    "LRG z=0.6-0.8": "desi_bao_dr2/desi_gaussian_bao_LRG_GCcomb_z0.6-0.8_mean.txt",
    "ELG z=1.1-1.6": "desi_bao_dr2/desi_gaussian_bao_ELG_LOPnotqso_GCcomb_z1.1-1.6_mean.txt",
    "QSO":          "desi_bao_dr2/desi_gaussian_bao_QSO_GCcomb_mean.txt",
    "Lyman-α":      "desi_bao_dr2/desi_gaussian_bao_Lya_GCcomb_mean.txt",
}


# ── Fetching ──────────────────────────────────────────────────────────────────

def _fetch_text(url: str) -> str:
    """Download a text file from `url`; raises on HTTP errors."""
    import requests
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    return resp.text


def parse_bao_mean(text: str) -> tuple[NDArray, NDArray, list[str]]:
    """Parse a CobayaSampler BAO mean file.

    Returns
    -------
    z_arr : array  — redshift of each measurement
    val_arr : array — distance ratio value
    qty_arr : list[str] — quantity type per row
    """
    rows = [
        line.split()
        for line in text.splitlines()
        if line.strip() and not line.startswith("#")
    ]
    z_arr   = np.array([float(r[0]) for r in rows])
    val_arr = np.array([float(r[1]) for r in rows])
    qty_arr = [r[2] for r in rows]
    return z_arr, val_arr, qty_arr


def parse_bao_cov(text: str) -> NDArray:
    """Parse a CobayaSampler BAO covariance file → 2-D array."""
    return np.loadtxt(io.StringIO(text))


def fetch_desi_bao(release: str = "dr2") -> dict:
    """Fetch and parse DESI BAO data for a given release ('dr2' or 'dr1').

    Returns
    -------
    dict with keys: z, val, qty, cov, sigma
    """
    mean_text = _fetch_text(URLS[f"{release}_mean"])
    cov_text  = _fetch_text(URLS[f"{release}_cov"])
    z, val, qty = parse_bao_mean(mean_text)
    cov = parse_bao_cov(cov_text)
    sigma = np.sqrt(np.diag(cov))
    return {"z": z, "val": val, "qty": qty, "cov": cov, "sigma": sigma}


# ── ΛCDM predictions ──────────────────────────────────────────────────────────

def sound_horizon_rd(
    H0: float = 67.36,
    Om0: float = 0.315,
    Ob0: float = 0.0493,
) -> float:
    """Approximate sound horizon at the baryon drag epoch [Mpc].

    Aubourg et al. (2015) fitting formula, accurate to ~0.1% for
    standard cosmologies (arXiv:1411.1074, Eq. 26).
    """
    h = H0 / 100.0
    omh2 = Om0 * h ** 2
    obh2 = Ob0 * h ** 2
    return 147.33 * (omh2 / 0.1430) ** (-0.255) * (obh2 / 0.02222) ** (-0.128)


def lcdm_distance_ratios(
    z_arr: NDArray,
    qty_arr: list[str],
    H0: float = 67.36,
    Om0: float = 0.315,
    Ob0: float = 0.0493,
) -> NDArray:
    """Compute ΛCDM predictions for the given BAO distance ratios.

    Parameters
    ----------
    z_arr, qty_arr
        Redshifts and quantity types from `fetch_desi_bao`.
    H0, Om0, Ob0
        Cosmological parameters.

    Returns
    -------
    pred : array of ΛCDM predicted values matching the order of z_arr/qty_arr.
    """
    from astropy.cosmology import FlatLambdaCDM
    cosmo = FlatLambdaCDM(H0=H0, Om0=Om0, Ob0=Ob0)
    r_d = sound_horizon_rd(H0=H0, Om0=Om0, Ob0=Ob0)

    preds = []
    for z, qty in zip(z_arr, qty_arr):
        if qty == "DM_over_rs":
            val = cosmo.comoving_transverse_distance(z).value / r_d
        elif qty == "DH_over_rs":
            val = (2.998e5 / cosmo.H(z).value) / r_d
        elif qty == "DV_over_rs":
            D_M = cosmo.comoving_transverse_distance(z).value
            D_H = 2.998e5 / cosmo.H(z).value
            D_V = (z * D_M ** 2 * D_H) ** (1.0 / 3.0)
            val = D_V / r_d
        else:
            val = np.nan
        preds.append(val)
    return np.array(preds)


# ── φ-oscillation model ───────────────────────────────────────────────────────

def phi_bao_residual_model(
    z: NDArray,
    A_phi: float,
    delta: float,
) -> NDArray:
    """φ-log-periodic modulation model for BAO residuals.

    ΔR(z) = A_φ × cos(2π × ln(1+z) / ln(φ) + δ)
    """
    return A_phi * np.cos(2 * np.pi * np.log(1.0 + z) / _LN_PHI + delta)


def fit_phi_oscillation(
    z: NDArray,
    residuals: NDArray,
    sigma: NDArray,
) -> tuple[NDArray, NDArray]:
    """Fit the φ-log-periodic model to BAO residuals.

    Returns
    -------
    popt : [A_phi, delta]
    pcov : 2×2 covariance matrix
    """
    amp_guess = float(np.std(residuals)) if len(residuals) > 1 else 0.1
    try:
        popt, pcov = curve_fit(
            phi_bao_residual_model,
            z,
            residuals,
            sigma=sigma,
            p0=[amp_guess, 0.0],
            absolute_sigma=True,
            maxfev=2000,
        )
    except RuntimeError:
        popt = np.array([0.0, 0.0])
        pcov = np.diag([np.inf, np.inf])
    return popt, pcov


def phi_scale_redshifts(n_range: tuple[int, int] = (-3, 5)) -> NDArray:
    """Redshifts at which ln(1+z) = n × ln(φ), i.e., z = φⁿ − 1."""
    return np.array([_PHI ** n - 1 for n in range(*n_range)])


# ── Main analysis ─────────────────────────────────────────────────────────────

def run_desi_bao_analysis(
    H0: float = 67.36,
    Om0: float = 0.315,
    Ob0: float = 0.0493,
) -> dict:
    """Fetch DR2 and DR1, compute ΛCDM predictions and φ-fit residuals.

    Returns
    -------
    dict with keys:
        dr2, dr1 — raw fetch dicts
        dr2_pred, dr1_pred — ΛCDM predictions
        dr2_resid, dr1_resid — fractional residuals (obs - pred) / pred
        dr2_popt, dr2_pcov — φ-fit parameters for DR2
        rd — sound horizon [Mpc]
        phi_z — φ-harmonic redshifts
    """
    dr2 = fetch_desi_bao("dr2")
    dr1 = fetch_desi_bao("dr1")

    dr2_pred = lcdm_distance_ratios(dr2["z"], dr2["qty"], H0=H0, Om0=Om0, Ob0=Ob0)
    dr1_pred = lcdm_distance_ratios(dr1["z"], dr1["qty"], H0=H0, Om0=Om0, Ob0=Ob0)

    dr2_resid = (dr2["val"] - dr2_pred) / dr2_pred
    dr1_resid = (dr1["val"] - dr1_pred) / dr1_pred
    dr2_sigma_frac = dr2["sigma"] / dr2_pred
    dr1_sigma_frac = dr1["sigma"] / dr1_pred

    popt, pcov = fit_phi_oscillation(dr2["z"], dr2_resid, dr2_sigma_frac)

    return {
        "dr2": dr2,
        "dr1": dr1,
        "dr2_pred": dr2_pred,
        "dr1_pred": dr1_pred,
        "dr2_resid": dr2_resid,
        "dr1_resid": dr1_resid,
        "dr2_sigma_frac": dr2_sigma_frac,
        "dr1_sigma_frac": dr1_sigma_frac,
        "dr2_popt": popt,
        "dr2_pcov": pcov,
        "rd": sound_horizon_rd(H0=H0, Om0=Om0, Ob0=Ob0),
        "phi_z": phi_scale_redshifts(),
    }


if __name__ == "__main__":
    res = run_desi_bao_analysis()
    print(f"Sound horizon r_d = {res['rd']:.2f} Mpc")
    print(f"DR2  A_φ = {res['dr2_popt'][0]:.4f}, δ = {res['dr2_popt'][1]:.3f}")
    print(f"DR2 residuals: {res['dr2_resid'].round(4)}")
    print(f"φ-harmonic redshifts: {res['phi_z'].round(3)}")
