"""
G–φ connection analysis — worked example.

Demonstrates the three independent angles on a potential link between
the gravitational constant G and the golden ratio φ:

  1. Varying-G model constrained by pulsar timing
  2. SGWB spectrum with φ-spaced ISCO frequencies
  3. GWTC-3 BH mass-ratio clustering test

Usage
-----
    python examples/g_phi_example.py
    python examples/g_phi_example.py --A_G_frac 0.5 --A_phi_gw 0.1
"""

import argparse
import os
import sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

from g_phi_analysis import (
    varying_g_model,
    gdot_over_g,
    pulsar_constraint_amplitude,
    f_isco,
    phi_harmonic_masses,
    omega_gw_powerlaw,
    omega_gw_phi_modulated,
    load_gwtc3_builtin,
    compute_mass_ratios,
    ks_test_phi_clustering,
    count_phi_ratio_matches,
    _G0,
    _PHI,
    _LN_PHI,
)

_T_NOW = 4.35e17  # s


def main(A_G_frac: float, A_phi_gw: float) -> None:
    # ── 1. Varying G ──────────────────────────────────────────────────────────
    A_G_max = pulsar_constraint_amplitude(_T_NOW)
    A_G = A_G_max * A_G_frac
    print("=" * 60)
    print("1. VARYING-G MODEL")
    print("=" * 60)
    print(f"  Pulsar constraint (Ġ/G < 10⁻¹² yr⁻¹):  A_G_max = {A_G_max:.3e}")
    print(f"  Using A_G = {A_G_frac:.0%} of limit:           A_G     = {A_G:.3e}")
    print(f"  G varies by ±{A_G*100:.4f}% over one φ-cycle")

    t_arr  = np.linspace(0.05 * _T_NOW, _T_NOW, 500)
    G_arr  = varying_g_model(t_arr, A_G=A_G)
    gd_arr = gdot_over_g(t_arr, A_G=A_G)
    print(f"  Max |Ġ/G| reached: {np.max(np.abs(gd_arr)):.3e} yr⁻¹  "
          f"({'within' if np.max(np.abs(gd_arr)) < 1e-12 else 'EXCEEDS'} pulsar limit)")

    # ── 2. GW spectrum ────────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("2. GW CHIRP FREQUENCY / SGWB SPECTRUM")
    print("=" * 60)
    phi_masses = phi_harmonic_masses(m_ref=30.0)
    f_isco_arr = f_isco(phi_masses)
    print(f"  {'n':>4}  {'M [M☉]':>10}  {'f_ISCO [Hz]':>12}  {'LIGO band':>10}")
    print("  " + "-" * 44)
    for i, (m, f) in enumerate(zip(phi_masses, f_isco_arr)):
        n = i - 4
        band = "✓" if 10 <= f <= 500 else ("LF" if f < 10 else "HF")
        print(f"  {n:>+4}  {m:>10.2f}  {f:>12.3f}  {band:>10}")

    f_arr = np.logspace(-2, 3, 800)
    omega_base = omega_gw_powerlaw(f_arr)
    omega_mod  = omega_gw_phi_modulated(f_arr, A_phi=A_phi_gw)
    max_mod_pct = np.max(np.abs(omega_mod / omega_base - 1)) * 100
    print(f"\n  SGWB modulation amplitude A_φ = {A_phi_gw:.2f}")
    print(f"  Max deviation from power-law:  {max_mod_pct:.1f}%")

    # ── 3. GWTC-3 mass ratios ─────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("3. GWTC-3 MASS-RATIO CLUSTERING TEST")
    print("=" * 60)
    m1, m2, events = load_gwtc3_builtin()
    q = compute_mass_ratios(m1, m2)
    log_q = np.log(q)
    ks_stat, p_val = ks_test_phi_clustering(log_q)
    n_matches, phi_lines = count_phi_ratio_matches(log_q)

    print(f"  Events in GWTC-3: {len(m1)}")
    print(f"  Mass ratio range: q = [{q.min():.2f}, {q.max():.2f}]")
    print(f"  KS statistic (vs uniform log q): {ks_stat:.3f}")
    print(f"  p-value: {p_val:.4f}  →  "
          f"{'reject' if p_val < 0.05 else 'consistent with'} uniform distribution")
    print(f"  Events within ±15% of a φⁿ ratio: {n_matches}")
    print(f"  Expected φⁿ values tested: {len(phi_lines)}")

    # ── Plot ──────────────────────────────────────────────────────────────────
    try:
        import matplotlib
        matplotlib.use("Agg")
    except ImportError as e:
        print(f"\nSkipping plot (matplotlib unavailable: {e})")
        return
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec

    fig = plt.figure(figsize=(14, 10))
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.35, wspace=0.3)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])

    t_Gyr = t_arr / 3.15576e16
    ax1.plot(t_Gyr, G_arr / _G0, color="#1f77b4", lw=2)
    ax1.axhline(1.0, ls="--", color="gray", alpha=0.6)
    ax1.set_xlabel("Cosmic time [Gyr]"); ax1.set_ylabel("G(t) / G₀")
    ax1.set_title("Varying G(t)  —  φ-Log-Periodic Model")
    ax1.grid(True, alpha=0.3)

    ax2.plot(t_Gyr, gd_arr, color="#d62728", lw=2)
    ax2.axhline(1e-12, ls=":", color="green", lw=1.5, label="Pulsar limit ±10⁻¹²")
    ax2.axhline(-1e-12, ls=":", color="green", lw=1.5)
    ax2.axhspan(-1e-12, 1e-12, alpha=0.06, color="green")
    ax2.set_xlabel("Cosmic time [Gyr]"); ax2.set_ylabel("Ġ/G  [yr⁻¹]")
    ax2.set_title("Time Derivative of G")
    ax2.legend(fontsize=8); ax2.grid(True, alpha=0.3)

    ax3.loglog(f_arr, omega_base, "--", color="gray", lw=1.5, label="Power-law baseline")
    ax3.loglog(f_arr, omega_mod, color="#1f77b4", lw=2,
               label=f"φ-modulated  (A_φ={A_phi_gw:.2f})")
    for i, (m, f) in enumerate(zip(phi_masses, f_isco_arr)):
        if 0.01 <= f <= 1000:
            ax3.axvline(f, ls="--", color="red", alpha=0.4)
            ax3.text(f * 1.05, omega_base[0] * 2, f"φ^{i-4}", fontsize=7, color="red")
    ax3.set_xlabel("f [Hz]"); ax3.set_ylabel("Ω_GW(f)")
    ax3.set_title("Stochastic GW Background + φ-ISCO Markers")
    ax3.legend(fontsize=8); ax3.grid(True, alpha=0.3, which="both")

    bins = np.arange(log_q.min() - 0.05, log_q.max() + 0.2, 0.15)
    ax4.hist(log_q, bins=bins, color="#1f77b4", alpha=0.7, edgecolor="white")
    for lp in phi_lines:
        if log_q.min() <= lp <= log_q.max():
            n = int(round(lp / _LN_PHI))
            ax4.axvline(lp, ls="--", color="red", alpha=0.5)
            ax4.text(lp + 0.02, ax4.get_ylim()[1] * 0.9 if ax4.get_ylim()[1] > 0 else 1,
                     f"φ^{n}", fontsize=8, color="red")
    ax4.set_xlabel("ln(q) = ln(m₁/m₂)"); ax4.set_ylabel("Number of events")
    ax4.set_title(f"GWTC-3 Mass Ratios  (KS p={p_val:.3f})")
    ax4.grid(True, alpha=0.3)

    out = "g_phi_example.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    print(f"\nFigure saved → {out}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="G–φ connection analysis example")
    parser.add_argument("--A_G_frac",  type=float, default=0.5,
                        help="Fraction of pulsar A_G limit to use (default: 0.5)")
    parser.add_argument("--A_phi_gw",  type=float, default=0.08,
                        help="SGWB φ-modulation amplitude (default: 0.08)")
    args = parser.parse_args()
    main(args.A_G_frac, args.A_phi_gw)
