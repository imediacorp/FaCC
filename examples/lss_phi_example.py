"""
LSS P(k) φ-scale feature detection — worked example.

Loads a matter power spectrum CSV, fits a smooth spline baseline,
identifies oscillatory residuals, and reports how many peaks
coincide with φ-spaced scales around k_BAO.

Usage
-----
    python examples/lss_phi_example.py
    python examples/lss_phi_example.py --data path/to/pk.csv
    python examples/lss_phi_example.py --data real_pk_lowk.csv --output pk_phi_example.png

Input CSV format
----------------
    k,Pk,sigma_Pk
    0.001,12345.6,617.3
    ...
"""

import argparse
import os
import sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

from lss_analysis import (
    load_pk_data,
    compute_phi_scales,
    compute_spline_residuals,
    find_direct_peaks,
    find_residual_peaks,
    count_phi_matches,
)


def main(data_path: str, output_path: str) -> None:
    print(f"Loading P(k) data from: {data_path}")
    k, pk, sigma = load_pk_data(data_path)
    print(f"  {len(k)} k-bins  |  k ∈ [{k.min():.4f}, {k.max():.4f}] h/Mpc")

    phi_scales = compute_phi_scales(k_bao=0.02)
    print(f"\nφ-spaced reference scales (φⁿ × k_BAO):")
    for n, ks in zip(range(-5, 6), phi_scales):
        in_range = "✓" if k.min() <= ks <= k.max() else "—"
        print(f"  n={n:+d}  k = {ks:.5f} h/Mpc  {in_range}")

    pk_smooth, residuals = compute_spline_residuals(k, pk)
    print(f"\nSpline baseline fitted.")

    k_peaks = find_direct_peaks(k, pk, sigma, threshold_sigma=3.0)
    k_res_peaks = find_residual_peaks(k, residuals, sigma, prominence_sigma=0.5)

    if len(k_peaks) > 0:
        n_match = count_phi_matches(k_peaks, phi_scales)
        print(f"\nDirect peaks: {len(k_peaks)}  |  φ-matches: {n_match}/{len(phi_scales)}")
        print(f"  Peak k-values: {k_peaks.round(5)}")
    else:
        print("\nNo direct peaks found (P(k) likely monotonically decreasing).")

    if len(k_res_peaks) > 0:
        n_res_match = count_phi_matches(k_res_peaks, phi_scales)
        print(f"\nResidual peaks: {len(k_res_peaks)}  |  φ-matches: {n_res_match}/{len(phi_scales)}")
        print(f"  Peak k-values: {k_res_peaks.round(5)}")
    else:
        print("\nNo significant residual peaks detected.")

    # Plot
    try:
        import matplotlib
        matplotlib.use("Agg")
    except ImportError as e:
        print(f"\nSkipping plot (matplotlib unavailable: {e})")
        return
    import matplotlib.pyplot as plt

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    ax1.errorbar(k, pk, yerr=sigma, fmt=".", alpha=0.6, label="P(k) data")
    ax1.plot(k, pk_smooth, "g-", lw=2, label="Spline baseline")
    for n, ks in zip(range(-5, 6), phi_scales):
        if k.min() <= ks <= k.max():
            ax1.axvline(ks, ls="--", color="r", alpha=0.3,
                        label="φ-scales" if n == -5 else "")
    if len(k_peaks) > 0:
        idx = np.searchsorted(k, k_peaks).clip(0, len(k) - 1)
        ax1.scatter(k_peaks, pk[idx], marker="*", s=150, color="orange",
                    zorder=5, label="Direct peaks")
    ax1.set_xscale("log"); ax1.set_yscale("log")
    ax1.set_ylabel("P(k) [(Mpc/h)³]")
    ax1.set_title("Matter Power Spectrum — φ-Scale Analysis")
    ax1.legend(fontsize=8); ax1.grid(True, alpha=0.3)

    ax2.errorbar(k, residuals, yerr=sigma, fmt=".", alpha=0.6, label="Residuals")
    ax2.axhline(0, color="k", lw=0.5)
    for ks in phi_scales:
        if k.min() <= ks <= k.max():
            ax2.axvline(ks, ls="--", color="r", alpha=0.3)
    if len(k_res_peaks) > 0:
        ri = np.searchsorted(k, k_res_peaks).clip(0, len(k) - 1)
        ax2.scatter(k_res_peaks, residuals[ri], marker="*", s=150, color="blue",
                    zorder=5, label="Residual peaks")
    ax2.set_xscale("log")
    ax2.set_xlabel("k [h/Mpc]"); ax2.set_ylabel("ΔP(k) [(Mpc/h)³]")
    ax2.set_title("Oscillation Residuals Around Smooth Baseline")
    ax2.legend(fontsize=8); ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    print(f"\nFigure saved → {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="LSS P(k) φ-scale analysis example")
    parser.add_argument("--data",   default=os.path.join(ROOT, "real_pk_lowk.csv"))
    parser.add_argument("--output", default="lss_phi_example.png")
    args = parser.parse_args()
    main(args.data, args.output)
