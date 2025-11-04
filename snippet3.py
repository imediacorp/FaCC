#!/usr/bin/env python
"""
extract_desi_pk.py
Extracts galaxy power spectrum data from a DESI FITS file.
"""

import argparse
import json
import numpy as np
from scipy.signal import find_peaks
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
from pathlib import Path


# ----------------------------------------------------------------------
def load_config(cfg_path: Path):
    """Loads configuration from a JSON file."""
    with open(cfg_path, "r") as f:
        return json.load(f)["lss"]

# ----------------------------------------------------------------------
def load_pk_data(csv_path: Path):
    """Loads power spectrum data from a CSV file."""
    data = np.loadtxt(csv_path, delimiter=",", skiprows=1)
    if data.shape[1] < 3:
        raise ValueError("CSV must have at least 3 columns: k, Pk, sigma_Pk")
    k, pk, sigma = data[:, 0], data[:, 1], data[:, 2]
    return k, pk, sigma

# ----------------------------------------------------------------------
def main():
   """Main function to execute the power spectrum extraction."""

    parser = argparse.ArgumentParser(description='Extract galaxy power spectrum data from FITS files.')
    parser.add_argument('--config', type=Path, default=Path('config.json'),
                        help='Path to the configuration file (JSON).')
    parser.add_argument('--data', type=Path, default=Path('pk_data.csv'),
                        help='Path to the power spectrum data CSV file.')
    parser.add_argument('--outfig', type=Path, default=Path('pk_phi.png'),
                        help='Path to save the resulting figure.')
    args = parser.parse_args()

    # Load configuration
    cfg = load_config(args.config)

    # Load power spectrum data
    k_data, pk_data, sigma_pk = load_pk_data(args.data)

    # Expected φ-scales (k_BAO)
    k_bao = cfg["k_bao"]

    # Direct peak finding
    height_threshold = sigma_pk.mean() * 3  # 3σ
    direct_peaks, _ = find_peaks(pk_data, height=height_threshold)
    k_direct = k_data[direct_peaks] if len(direct_peaks) > 0 else np.array([])

    # Residuals analysis
    log_k = np.log(k_data)
    log_pk = np.log(np.maximum(pk_data, 1e-10))  # Safe log
    spline = UnivariateSpline(log_k, log_pk, s=len(k_data) * 0.5, k=3)
    pk_smooth = np.exp(spline(log_k))
    residuals = pk_data - pk_smooth

    prom_factor = cfg.get("peak_prominence_factor", 0.5)
    res_peaks, _ = find_peaks(np.abs(residuals), prominence=sigma_pk.mean() * prom_factor)
    k_res = k_data[res_peaks] if len(res_peaks) > 0 else np.array([])

    # Matching
    tolerance = cfg["tolerance"]

    direct_matches = count_matches(k_direct, [k_bao], tolerance)
    res_matches = count_matches(k_res, [k_bao], tolerance)

    # Reporting
    print("=" * 60)
    print("LSS φ-SCALE ANALYSIS")
    print("=" * 60)
    print(f"Data points      : {len(k_data)}")
    print(f"φ-scales (k_BAO={k_bao:.4f}) : {len([k_bao])}")
    print("\n--- Direct peaks ---")
    print(f"Found {len(k_direct)} peaks → {direct_matches}/{len([k_bao])} matches")
    if len(k_direct):
        print(f"   k = {k_direct}")
    print("\n--- Residual peaks ---")
    print(f"Found {len(k_res)} peaks → {res_matches}/{len([k_bao])} matches")
    if len(k_res):
        print(f"   k = {k_res}")

    # Plotting
    fig = plt.figure(figsize=(10, 10), constrained_layout=True)
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 1], hspace=0.1)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex=ax1)

    # Top panel
    ax1.errorbar(k_data, pk_data, yerr=sigma_pk, fmt=".", label="Data", alpha=0.7, capsize=2)
    ax1.plot(k_data, pk_smooth, "g-", lw=2, label="Smooth spline", alpha=0.8)
    for i, ks in enumerate(k_bao):
        ax1.axvline(ks, ls="--", color="r", alpha=0.4)
        if i == 0 or i == len(k_bao)-1 or abs(ks - k_bao) < 1e-4:
            ax1.text(ks, ax1.get_ylim()[1]*0.95, f"{ks:.4f}", rotation=90,
                     va="top", ha="center", fontsize=8, color="r", alpha=0.7)
    if len(k_direct):
        ax1.scatter(k_direct, pk_data[direct_peaks], color="orange", s=120, marker="*", label="Direct peaks", zorder=5)
    ax1.set_xscale("log")
    ax1.set_ylabel(r"$P(k)\;[({\rm Mpc}/h)^3]$")
    ax1.set_title("Galaxy Power Spectrum")
    ax1.legend(loc="upper right")
    ax1.grid(True, alpha=0.3)

    # Bottom panel
    ax2.errorbar(k_data, residuals, yerr=sigma_pk, fmt=".", label="Residuals", alpha=0.7, capsize=2)
    ax2.axhline(0, color="k", lw=1)
    for i, ks in enumerate(k_bao):
        ax2.axvline(ks, ls="--", color="r", alpha=0.4)
    if len(k_res):
        ax2.scatter(k_res, residuals[res_peaks], color="blue", s=120, marker="*", label="Residual peaks", zorder=5)
    ax2.set_xscale("log")
    ax2.set_xlabel(r"$k\; [h\,{\rm Mpc}^{-1}]$")
    ax2.set_ylabel(r"$\Delta P(k)$")
    ax2.set_title("Oscillations around Smooth Fit")
    ax2.legend(loc="upper right")
    ax2.grid(True, alpha=0.3)

    plt.savefig(args.outfig, dpi=250, bbox_inches=None)
    print(f"\nFigure saved → {args.outfig}")
    plt.show()

def count_matches(k_vals, k_baos, tolerance):
    """Counts number of times k_vals match a k_bao within tolerance."""
    matches = 0
    for k_bao in k_baos:
        if np.any(np.abs(k_vals - k_bao) < tolerance):
            matches += 1
    return matches

if __name__ == "__main__":
    main()
