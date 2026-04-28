"""
DESI DR2 BAO φ-oscillation search — worked example.

Fetches the public DESI DR2 and DR1 BAO distance measurements,
computes ΛCDM predictions, and fits a φ-log-periodic model to the
fractional residuals ΔR(z) = (obs − ΛCDM) / ΛCDM.

Usage
-----
    python examples/desi_bao_example.py
    python examples/desi_bao_example.py --H0 67.4 --Om0 0.31

Requires an internet connection on first run (data cached by OS).
"""

import argparse
import os
import sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

from desi_bao import (
    fetch_desi_bao,
    lcdm_distance_ratios,
    sound_horizon_rd,
    phi_bao_residual_model,
    fit_phi_oscillation,
    phi_scale_redshifts,
)

_PHI = (1 + np.sqrt(5)) / 2


def main(H0: float, Om0: float, Ob0: float) -> None:
    print(f"Cosmology: H₀={H0:.2f}  Ω_m={Om0:.3f}  Ω_b={Ob0:.4f}")
    rd = sound_horizon_rd(H0=H0, Om0=Om0, Ob0=Ob0)
    print(f"Sound horizon r_d = {rd:.3f} Mpc\n")

    print("Fetching DESI DR2 BAO data …")
    dr2 = fetch_desi_bao("dr2")
    dr2_pred = lcdm_distance_ratios(dr2["z"], dr2["qty"], H0=H0, Om0=Om0, Ob0=Ob0)
    dr2_resid = (dr2["val"] - dr2_pred) / dr2_pred
    dr2_sigma_frac = dr2["sigma"] / dr2_pred

    print("Fetching DESI DR1 BAO data …")
    dr1 = fetch_desi_bao("dr1")
    dr1_pred = lcdm_distance_ratios(dr1["z"], dr1["qty"], H0=H0, Om0=Om0, Ob0=Ob0)
    dr1_resid = (dr1["val"] - dr1_pred) / dr1_pred

    print(f"\n{'z':>6}  {'Quantity':>14}  {'Observed':>10}  {'ΛCDM':>10}  "
          f"{'Residual%':>10}  {'σ%':>6}")
    print("-" * 64)
    for z, val, qty, pred, res, sig in zip(
        dr2["z"], dr2["val"], dr2["qty"], dr2_pred, dr2_resid, dr2_sigma_frac
    ):
        print(f"{z:6.3f}  {qty:>14}  {val:10.4f}  {pred:10.4f}  "
              f"{res*100:+9.2f}%  {sig*100:5.2f}%")

    popt, pcov = fit_phi_oscillation(dr2["z"], dr2_resid, dr2_sigma_frac)
    perr = np.sqrt(np.diag(pcov))
    A_phi, delta = popt
    significance = abs(A_phi) / perr[0] if perr[0] < np.inf else 0.0

    print(f"\nφ-oscillation fit to DR2 residuals:")
    print(f"  A_φ = {A_phi:.5f} ± {perr[0]:.5f}  ({significance:.1f}σ)")
    print(f"  δ   = {delta:.3f} ± {perr[1]:.3f} rad")
    print(f"\nφ-harmonic redshifts (z = φⁿ − 1):")
    for pz in phi_scale_redshifts():
        if 0.05 < pz < 2.5:
            print(f"  z = {pz:.3f}")

    # Plot
    try:
        import matplotlib
        matplotlib.use("Agg")
    except ImportError as e:
        print(f"\nSkipping plot (matplotlib unavailable: {e})")
        return
    import matplotlib.pyplot as plt

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

    qty_color = {"DM_over_rs": "#1f77b4", "DH_over_rs": "#ff7f0e", "DV_over_rs": "#2ca02c"}
    qty_label = {"DM_over_rs": "D_M/r_d", "DH_over_rs": "D_H/r_d", "DV_over_rs": "D_V/r_d"}

    z_fine = np.linspace(0.1, 2.5, 300)
    for qty, col in qty_color.items():
        pred_fine = lcdm_distance_ratios(z_fine, [qty] * len(z_fine), H0=H0, Om0=Om0, Ob0=Ob0)
        ax1.plot(z_fine, pred_fine, "--", color=col, lw=1.5, alpha=0.7, label=f"ΛCDM {qty_label[qty]}")

    seen = set()
    for z, val, qty, sig in zip(dr2["z"], dr2["val"], dr2["qty"], dr2["sigma"]):
        col = qty_color[qty]
        lbl = f"DR2 {qty_label[qty]}" if qty not in seen else ""
        ax1.errorbar(z, val, yerr=sig, fmt="o", color=col, ms=7, label=lbl, capsize=3)
        seen.add(qty)
    for z, val, qty, sig in zip(dr1["z"], dr1["val"], dr1["qty"], dr1["sigma"]):
        col = qty_color[qty]
        ax1.errorbar(z, val, yerr=sig, fmt="s", color=col, ms=6, mfc="white",
                     mew=1.5, capsize=3, alpha=0.6)

    ax1.set_xlabel("Redshift z"); ax1.set_ylabel("Distance ratio")
    ax1.set_title("DESI BAO Distance Ratios vs ΛCDM")
    ax1.legend(fontsize=8, ncol=2); ax1.grid(True, alpha=0.3)

    z_fit = np.linspace(0.1, 2.5, 400)
    fit_curve = phi_bao_residual_model(z_fit, *popt)
    ax2.plot(z_fit, fit_curve * 100, "r-", lw=2,
             label=f"φ-fit  A_φ={A_phi:.4f} ({significance:.1f}σ)")
    ax2.axhline(0, color="k", lw=0.5)

    seen2 = set()
    for z, res, sig, qty in zip(dr2["z"], dr2_resid, dr2_sigma_frac, dr2["qty"]):
        col = qty_color[qty]
        lbl = f"DR2 {qty_label[qty]}" if qty not in seen2 else ""
        ax2.errorbar(z, res * 100, yerr=sig * 100, fmt="o", color=col,
                     ms=7, label=lbl, capsize=3)
        seen2.add(qty)

    for pz in phi_scale_redshifts():
        if 0.05 < pz < 2.5:
            ax2.axvline(pz, ls=":", color="purple", alpha=0.4)

    ax2.set_xlabel("Redshift z"); ax2.set_ylabel("(obs − ΛCDM) / ΛCDM  [%]")
    ax2.set_title("Fractional Residuals + φ-Oscillation Fit")
    ax2.legend(fontsize=8, ncol=2); ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    out = "desi_bao_example.png"
    plt.savefig(out, dpi=150)
    print(f"\nFigure saved → {out}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="DESI DR2 BAO φ-oscillation example")
    parser.add_argument("--H0",  type=float, default=67.36)
    parser.add_argument("--Om0", type=float, default=0.315)
    parser.add_argument("--Ob0", type=float, default=0.0493)
    args = parser.parse_args()
    main(args.H0, args.Om0, args.Ob0)
