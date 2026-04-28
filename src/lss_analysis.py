"""
Large-scale structure (LSS) P(k) analysis for φ-scale oscillation detection.

Searches for φ-spaced features in the matter power spectrum using direct
peak-finding and spline-residual analysis.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from scipy.signal import find_peaks
from scipy.interpolate import UnivariateSpline
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM  # noqa: F401 — available for callers

_PHI = (1 + np.sqrt(5)) / 2


def load_pk_data(data_path: str) -> tuple[NDArray, NDArray, NDArray]:
    """Load P(k) data from a CSV file.

    Parameters
    ----------
    data_path : str
        Path to CSV with columns: k, P(k), sigma_Pk.

    Returns
    -------
    k_data, pk_data, sigma_pk : arrays
    """
    data = np.loadtxt(data_path, delimiter=",", skiprows=1)
    return data[:, 0], data[:, 1], data[:, 2]


def compute_phi_scales(
    k_bao: float = 0.02, n_range: tuple[int, int] = (-5, 6)
) -> NDArray:
    """Compute expected φ-spaced k-scales around the BAO scale.

    Parameters
    ----------
    k_bao : float
        BAO scale [h/Mpc].
    n_range : tuple[int, int]
        Range of exponents (inclusive start, exclusive end).

    Returns
    -------
    phi_scales : array
        k values φ^n × k_BAO.
    """
    return k_bao * np.array([_PHI ** n for n in range(*n_range)])


def find_direct_peaks(
    k_data: NDArray, pk_data: NDArray, sigma_pk: NDArray, threshold_sigma: float = 3.0
) -> NDArray:
    """Find peaks in P(k) above a σ threshold.

    Parameters
    ----------
    threshold_sigma : float
        Minimum height relative to mean σ.

    Returns
    -------
    k_peaks : array
        k-values of detected peaks.
    """
    peaks, _ = find_peaks(pk_data, height=sigma_pk.mean() * threshold_sigma)
    return k_data[peaks]


def compute_spline_residuals(
    k_data: NDArray, pk_data: NDArray
) -> tuple[NDArray, NDArray]:
    """Fit a smooth spline baseline and return residuals.

    Returns
    -------
    pk_smooth : array
        Smooth spline baseline.
    residuals : array
        pk_data - pk_smooth.
    """
    log_k = np.log(k_data)
    log_pk = np.log(pk_data)
    spline = UnivariateSpline(log_k, log_pk, s=len(k_data) * 0.5)
    pk_smooth = np.exp(spline(log_k))
    return pk_smooth, pk_data - pk_smooth


def find_residual_peaks(
    k_data: NDArray, residuals: NDArray, sigma_pk: NDArray, prominence_sigma: float = 0.5
) -> NDArray:
    """Find peaks in P(k) residuals above a prominence threshold.

    Returns
    -------
    k_peaks : array
        k-values of residual peaks.
    """
    peaks, _ = find_peaks(residuals, height=0, prominence=sigma_pk.mean() * prominence_sigma)
    return k_data[peaks]


def count_phi_matches(k_peaks: NDArray, phi_scales: NDArray, tolerance: float = 0.005) -> int:
    """Count how many φ-scales have a detected peak within `tolerance`."""
    if len(k_peaks) == 0:
        return 0
    return sum(
        np.min(np.abs(k_peaks - ks)) < tolerance for ks in phi_scales
    )


def run_lss_analysis(
    data_path: str = "real_pk_lowk.csv",
    output_path: str = "pk_phi.png",
    k_bao: float = 0.02,
) -> dict:
    """Run full LSS P(k) φ-scale analysis and save figure.

    Parameters
    ----------
    data_path : str
        Path to P(k) CSV data file.
    output_path : str
        Where to save the analysis figure.
    k_bao : float
        BAO scale used to generate φ-spaced reference scales [h/Mpc].

    Returns
    -------
    dict with keys: k_data, pk_data, sigma_pk, phi_scales, k_peaks,
                    k_residual_peaks, pk_smooth, residuals
    """
    k_data, pk_data, sigma_pk = load_pk_data(data_path)
    phi_scales = compute_phi_scales(k_bao=k_bao)
    pk_smooth, residuals = compute_spline_residuals(k_data, pk_data)
    k_peaks = find_direct_peaks(k_data, pk_data, sigma_pk)
    k_residual_peaks = find_residual_peaks(k_data, residuals, sigma_pk)

    print("=" * 60)
    print("ANALYSIS RESULTS")
    print("=" * 60)

    if len(k_peaks) > 0:
        n_match = count_phi_matches(k_peaks, phi_scales)
        print(f"\nDirect peaks found: {len(k_peaks)}")
        print(f"φ-scale matches:    {n_match} / {len(phi_scales)}")
        print(f"Peak locations (k): {k_peaks}")
    else:
        print("\nNo direct peaks found. P(k) appears monotonically decreasing.")

    if len(k_residual_peaks) > 0:
        n_res_match = count_phi_matches(k_residual_peaks, phi_scales)
        print(f"\nResidual peaks found:          {len(k_residual_peaks)}")
        print(f"φ-scale matches in residuals:  {n_res_match} / {len(phi_scales)}")
        print(f"Residual peak locations (k):   {k_residual_peaks}")
    else:
        print("\nNo significant oscillations detected around smooth fit.")

    print("\nExpected φ-scales:")
    for i, ks in enumerate(phi_scales):
        print(f"  φ^{i-5} × k_BAO = {ks:.6f} h/Mpc")

    # --- Plot ---
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))

    ax1.errorbar(k_data, pk_data, yerr=sigma_pk, fmt=".", label="Data", alpha=0.7)
    ax1.plot(k_data, pk_smooth, "g-", label="Smooth fit", linewidth=2, alpha=0.7)
    for ks in phi_scales:
        ax1.axvline(ks, ls="--", color="r", alpha=0.3)
    if len(k_peaks) > 0:
        peaks_idx = np.searchsorted(k_data, k_peaks)
        ax1.scatter(
            k_peaks, pk_data[peaks_idx], color="orange", s=100, marker="*",
            label="Direct peaks", zorder=5
        )
    ax1.set_xscale("log")
    ax1.set_xlabel("k [h/Mpc]")
    ax1.set_ylabel("P(k) [(Mpc/h)³]")
    ax1.set_title("Matter Power Spectrum")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    ax2.errorbar(k_data, residuals, yerr=sigma_pk, fmt=".", label="Residuals", alpha=0.7)
    ax2.axhline(0, color="k", linewidth=0.5)
    for j, ks in enumerate(phi_scales):
        ax2.axvline(ks, ls="--", color="r", alpha=0.3, label="φ-scales" if j == 0 else "")
    if len(k_residual_peaks) > 0:
        res_idx = np.searchsorted(k_data, k_residual_peaks)
        ax2.scatter(
            k_residual_peaks, residuals[res_idx], color="blue", s=100,
            marker="*", label="Residual peaks", zorder=5
        )
    ax2.set_xscale("log")
    ax2.set_xlabel("k [h/Mpc]")
    ax2.set_ylabel("ΔP(k) [(Mpc/h)³]")
    ax2.set_title("Oscillations Around Smooth Fit")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"\nFigure saved → {output_path}")

    return {
        "k_data": k_data,
        "pk_data": pk_data,
        "sigma_pk": sigma_pk,
        "phi_scales": phi_scales,
        "k_peaks": k_peaks,
        "k_residual_peaks": k_residual_peaks,
        "pk_smooth": pk_smooth,
        "residuals": residuals,
    }


if __name__ == "__main__":
    run_lss_analysis()
