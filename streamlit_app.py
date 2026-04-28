"""
FaCC — Fibonacci as Cosmological Constant
Interactive Analysis Dashboard

Tabs
----
1. 📊 Forecast        — DESI Year-5 sensitivity (CAMB required)
2. 🔬 Systematics     — Systematic error budget  (CAMB required)
3. 📈 Power Spectrum  — P(k) visualisation        (CAMB required)
4. 🌌 H(z) Analysis  — Falsified background model (data only)
5. 📡 CMB Residuals  — Dual φ/φ⁻¹ oscillation fit (data only)
6. 🔭 LSS P(k)       — φ-scale feature detection  (data only)
7. 📋 Summary         — Key results overview

Run:
    streamlit run streamlit_app.py
"""
from __future__ import annotations

import os
import sys

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import streamlit as st

# ── path setup ────────────────────────────────────────────────────────────────
ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(ROOT, "src"))

# ── optional heavy imports (CAMB / Fortran) ───────────────────────────────────
try:
    from phi_modulation import PhiModulationModel
    HAS_CAMB = True
except Exception:
    HAS_CAMB = False

try:
    from systematics import SystematicErrorBudget
    HAS_SYSTEMATICS = True
except Exception:
    HAS_SYSTEMATICS = False

# ── lightweight data-analysis imports (always available) ─────────────────────
from hz_analysis import (
    load_hz_data, fit_forward_branch, fit_reverse_branch, h_model,
)
from cmb_analysis import (
    load_cmb_data, preprocess_cmb_data, compute_residuals,
    fit_dual_oscillations, dual_osc_model,
)
from lss_analysis import (
    load_pk_data, compute_phi_scales, compute_spline_residuals,
    find_direct_peaks, find_residual_peaks,
)
from astropy.cosmology import FlatLambdaCDM

_PHI = (1 + np.sqrt(5)) / 2

# ── page config ───────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="FaCC — φ-Modulation Dashboard",
    page_icon="🌌",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown(
    """
    <style>
    .block-container { padding-top: 1.2rem; }
    div[data-testid="metric-container"] {
        background: #f4f6fa; border-radius: 8px; padding: 0.4rem 0.8rem;
    }
    </style>
    """,
    unsafe_allow_html=True,
)

st.title("φ-Modulation Analysis  ·  FaCC Dashboard")
st.caption(
    "Fibonacci as Cosmological Constant — "
    "falsified background · testable perturbations"
)

# ── sidebar ───────────────────────────────────────────────────────────────────
sb = st.sidebar
sb.header("⚙️  Parameters")

if HAS_CAMB:
    _m0 = PhiModulationModel()
    sb.info(f"φ = {_m0.phi:.8f}   ·   ln φ = {_m0.lnphi:.6f}")
else:
    sb.warning("CAMB not installed — forecast tabs disabled.")
sb.markdown("---")

with sb.expander("📐  Modulation & Forecast", expanded=True):
    A_phi = sb.slider(
        "A_φ  (modulation amplitude)",
        0.001, 0.05, 0.01, 0.001, format="%.3f",
    )
    z_eff = sb.slider("Redshift  z_eff", 0.0, 2.0, 0.8, 0.1)
    k_min = sb.slider(
        "k_min  [h/Mpc]", 0.001, 0.05, 0.01, 0.001, format="%.3f",
    )
    k_max = sb.slider("k_max  [h/Mpc]", 0.10, 1.0, 0.30, 0.01)
    include_sys = sb.checkbox(
        "Include systematic errors", True,
        disabled=not HAS_SYSTEMATICS,
    )

with sb.expander("🔭  Cosmology (Planck 2018)", expanded=False):
    H0    = sb.slider("H₀  [km/s/Mpc]", 60.0, 75.0, 67.36, 0.1)
    ombh2 = sb.slider("Ω_b h²", 0.020, 0.025, 0.02237, 0.0001, format="%.4f")
    omch2 = sb.slider("Ω_c h²", 0.10, 0.14, 0.1200, 0.001)
    ns    = sb.slider("n_s", 0.90, 1.00, 0.9649, 0.001)

sb.markdown("---")
sb.markdown(
    "**FaCC Framework**  \n"
    "[GitHub ↗](https://github.com/imediacorp/FaCC)  ·  "
    "Bryan David Persaud"
)


# ── cached computations ───────────────────────────────────────────────────────

@st.cache_data(show_spinner="Running DESI forecast…")
def _forecast(A_phi, k_min, k_max, H0, ombh2, omch2, ns, include_sys):
    params = dict(H0=H0, ombh2=ombh2, omch2=omch2,
                  As=2.1e-9, ns=ns, tau=0.0544)
    m = PhiModulationModel(params=params)
    if include_sys and HAS_SYSTEMATICS:
        return m.forecast_desi_sensitivity_with_systematics(
            A_phi_true=A_phi, k_min=k_min, k_max=k_max,
            n_k=100, include_systematics=True,
        )
    return m.forecast_desi_sensitivity(
        A_phi_true=A_phi, k_min=k_min, k_max=k_max, n_k=100,
    )


@st.cache_data(show_spinner="Computing power spectrum…")
def _power_spectrum(k_min, k_max, z_eff, H0, ombh2, omch2, ns):
    params = dict(H0=H0, ombh2=ombh2, omch2=omch2,
                  As=2.1e-9, ns=ns, tau=0.0544)
    m = PhiModulationModel(params=params)
    return m.get_base_power_spectrum(
        k_min=k_min * 0.5, k_max=k_max * 2, npoints=500, z=z_eff,
    )


@st.cache_data(show_spinner="Fitting H(z) data…")
def _hz():
    path = os.path.join(ROOT, "real_hz.csv")
    z, h, s = load_hz_data(path)
    res_f = fit_forward_branch(z, h, s)
    res_r = fit_reverse_branch(z, h, s)
    return z, h, s, res_f, res_r


@st.cache_data(show_spinner="Fitting CMB residuals…")
def _cmb():
    path = os.path.join(ROOT, "real_cmb_lowl.csv")
    ell_raw, cl_raw, sig_raw = load_cmb_data(path)
    ell, cl, sig = preprocess_cmb_data(ell_raw, cl_raw, sig_raw)
    cl_base, residuals = compute_residuals(ell, cl)
    popt, pcov = fit_dual_oscillations(ell, residuals, sig)
    return ell, cl, sig, cl_base, residuals, popt, pcov


@st.cache_data(show_spinner="Analysing LSS P(k)…")
def _lss():
    path = os.path.join(ROOT, "real_pk_lowk.csv")
    k, pk, s = load_pk_data(path)
    phi_scales = compute_phi_scales()
    pk_smooth, residuals = compute_spline_residuals(k, pk)
    k_peaks = find_direct_peaks(k, pk, s)
    k_res   = find_residual_peaks(k, residuals, s)
    return k, pk, s, phi_scales, pk_smooth, residuals, k_peaks, k_res


# ── shared plotly theme ────────────────────────────────────────────────────────
_LAYOUT = dict(template="plotly_white", font=dict(family="sans-serif", size=12))
_COLORS = dict(
    base="#1a1a2e", mod="#1f77b4", lcdm="#555555",
    stat="#1a1a2e", sys="#d62728", total="#1f77b4",
    fwd="#1f77b4",  rev="#ff7f0e",
    photoz="#d62728", bias="#2ca02c", geom="#1f77b4",
    phi_line="rgba(214,39,40,0.35)",
)


# ── tab layout ────────────────────────────────────────────────────────────────
(
    tab_forecast, tab_sys, tab_ps,
    tab_hz, tab_cmb, tab_lss,
    tab_summary,
) = st.tabs([
    "📊 Forecast",
    "🔬 Systematics",
    "📈 Power Spectrum",
    "🌌 H(z) Analysis",
    "📡 CMB Residuals",
    "🔭 LSS P(k)",
    "📋 Summary",
])


# ══════════════════════════════════════════════════════════════════════════════
# Tab 1 — Forecast Analysis
# ══════════════════════════════════════════════════════════════════════════════
with tab_forecast:
    st.header("DESI Year-5 Forecast")

    if not HAS_CAMB:
        st.error(
            "**CAMB is not installed.**  "
            "Install a Fortran compiler then `pip install camb`.  "
            "The H(z), CMB, and LSS tabs work without it."
        )
    else:
        try:
            fc = _forecast(A_phi, k_min, k_max, H0, ombh2, omch2, ns, include_sys)

            s_stat = fc.get("sigma_Aphi_stat", fc["sigma_Aphi"])
            s_sys  = fc.get("sigma_Aphi_sys",  0.0)
            s_tot  = fc.get("sigma_Aphi_total", fc["sigma_Aphi"])
            SNR    = fc["SNR"]
            k_fc   = fc["k"]
            Pk_b   = fc["Pk_base"]
            Pk_m   = fc["Pk_mod"]
            mf     = fc["mod_factor"]

            # ── metrics ──────────────────────────────────────────────────────
            c1, c2, c3, c4 = st.columns(4)
            c1.metric("A_φ (true)",         f"{A_phi:.4f}")
            c2.metric("σ_Aφ  statistical",  f"{s_stat:.2e}")
            c3.metric(
                "σ_Aφ  systematic",
                f"{s_sys:.2e}" if (include_sys and HAS_SYSTEMATICS) else "—",
            )
            delta_snr = f"+{SNR - 3:.2f} σ vs 3σ" if SNR > 3 else f"{SNR - 3:.2f} σ vs 3σ"
            c4.metric("SNR", f"{SNR:.2f} σ", delta=delta_snr)
            st.markdown("---")

            # ── 2 × 2 chart grid ─────────────────────────────────────────────
            fig = make_subplots(
                rows=2, cols=2,
                subplot_titles=[
                    "φ-Modulation Signal  (ΔP/P)",
                    "Power Spectrum Comparison",
                    "Error Budget vs k",
                    "Detection Significance",
                ],
                vertical_spacing=0.14,
                horizontal_spacing=0.10,
            )

            # P1 — modulation ratio
            ratio_pct = (Pk_m / Pk_b - 1) * 100
            fig.add_trace(
                go.Scatter(
                    x=k_fc, y=ratio_pct, name="ΔP/P [%]",
                    line=dict(color=_COLORS["mod"], width=2),
                    fill="tozeroy", fillcolor="rgba(31,119,180,0.10)",
                ),
                row=1, col=1,
            )
            fig.add_hline(y=0, line_dash="dash", line_color="gray",
                          opacity=0.5, row=1, col=1)

            # P2 — P(k) comparison
            fig.add_trace(
                go.Scatter(x=k_fc, y=Pk_b, name="ΛCDM baseline",
                           line=dict(color=_COLORS["lcdm"], width=2)),
                row=1, col=2,
            )
            fig.add_trace(
                go.Scatter(x=k_fc, y=Pk_m, name="φ-modulated",
                           line=dict(color=_COLORS["mod"], width=2, dash="dash")),
                row=1, col=2,
            )

            # P3 — error budget
            sigma_P = fc.get("sigma_P", Pk_b * 0.05)
            fig.add_trace(
                go.Scatter(x=k_fc, y=sigma_P / Pk_b * 100, name="Statistical",
                           line=dict(color=_COLORS["stat"], width=2)),
                row=2, col=1,
            )
            if include_sys and HAS_SYSTEMATICS and "systematic_budget" in fc:
                sbd = fc["systematic_budget"]
                fig.add_trace(
                    go.Scatter(x=k_fc, y=sbd["sigma_P_sys"] / Pk_b * 100,
                               name="Systematic",
                               line=dict(color=_COLORS["sys"], width=2, dash="dash")),
                    row=2, col=1,
                )
                fig.add_trace(
                    go.Scatter(x=k_fc, y=sbd["sigma_P_total"] / Pk_b * 100,
                               name="Total",
                               line=dict(color=_COLORS["total"], width=2, dash="dot")),
                    row=2, col=1,
                )

            # P4 — SNR bars
            if include_sys and HAS_SYSTEMATICS and s_sys > 0:
                cats   = ["Statistical", "Systematic", "Total"]
                vals   = [A_phi / s_stat, A_phi / s_sys, SNR]
                colors = [_COLORS["stat"], _COLORS["sys"], _COLORS["total"]]
            else:
                cats   = ["Statistical", "Total"]
                vals   = [SNR, SNR]
                colors = [_COLORS["stat"], _COLORS["total"]]
            fig.add_trace(
                go.Bar(
                    x=cats, y=vals,
                    marker_color=colors, opacity=0.8,
                    text=[f"{v:.2f}σ" for v in vals],
                    textposition="outside",
                    name="SNR", showlegend=False,
                ),
                row=2, col=2,
            )
            fig.add_hline(y=3, line_dash="dash", line_color="red",
                          annotation_text="3σ", annotation_font_color="red",
                          row=2, col=2)
            fig.add_hline(y=5, line_dash="dot", line_color="orange",
                          annotation_text="5σ", annotation_font_color="orange",
                          row=2, col=2)

            # Axes
            for r, c in [(1, 1), (1, 2), (2, 1)]:
                fig.update_xaxes(type="log", title_text="k [h/Mpc]", row=r, col=c)
            fig.update_yaxes(type="log", title_text="P(k) [(Mpc/h)³]", row=1, col=2)
            fig.update_yaxes(title_text="ΔP/P [%]",           row=1, col=1)
            fig.update_yaxes(title_text="Relative error [%]", row=2, col=1)
            fig.update_yaxes(title_text="SNR [σ]",            row=2, col=2)
            fig.update_layout(
                height=680,
                legend=dict(orientation="h", yanchor="bottom", y=1.01, x=0),
                **_LAYOUT,
            )
            st.plotly_chart(fig, use_container_width=True)

        except Exception as exc:
            st.error(f"Forecast failed: {exc}")
            st.exception(exc)


# ══════════════════════════════════════════════════════════════════════════════
# Tab 2 — Systematic Error Budget
# ══════════════════════════════════════════════════════════════════════════════
with tab_sys:
    st.header("Systematic Error Budget")

    if not HAS_CAMB:
        st.error("CAMB required for this tab.")
    elif not HAS_SYSTEMATICS:
        st.warning("systematics module unavailable.")
    else:
        try:
            k_ps, _, Pk_ps = _power_spectrum(k_min, k_max, z_eff, H0, ombh2, omch2, ns)

            from scipy.interpolate import interp1d
            k_tgt = np.logspace(np.log10(k_min), np.log10(k_max), 120)
            Pk_tgt = interp1d(
                k_ps, Pk_ps[0], kind="linear",
                bounds_error=False, fill_value="extrapolate",
            )(k_tgt)
            sig_stat = Pk_tgt * 0.10

            seb = SystematicErrorBudget(z_eff=z_eff)
            res = seb.compute_systematic_budget(k_tgt, Pk_tgt, sig_stat)

            def _share(arr):
                return float(np.mean(arr / (res["sigma_P_total"] + 1e-30)) * 100)

            c1, c2, c3 = st.columns(3)
            c1.metric("Photo-z share",  f"{_share(res['sigma_P_photoz']):.1f}%")
            c2.metric("Bias share",     f"{_share(res['sigma_P_bias']):.1f}%")
            c3.metric("Geometry share", f"{_share(res['sigma_P_geometry']):.1f}%")
            st.markdown("---")

            fig = make_subplots(
                rows=2, cols=2,
                subplot_titles=[
                    "Individual Systematic Components",
                    "Total Error Breakdown",
                    "Systematic Fraction vs k",
                    "Mean-Squared Contributions",
                ],
                specs=[[{}, {}], [{}, {"type": "pie"}]],
                vertical_spacing=0.14,
                horizontal_spacing=0.10,
            )

            # P1 — components
            for label, arr, color, dash in [
                ("Photo-z",   res["sigma_P_photoz"],   _COLORS["photoz"], "solid"),
                ("Bias",      res["sigma_P_bias"],      _COLORS["bias"],   "dash"),
                ("Geometry",  res["sigma_P_geometry"],  _COLORS["geom"],   "dot"),
            ]:
                fig.add_trace(
                    go.Scatter(x=k_tgt, y=arr / Pk_tgt * 100, name=label,
                               line=dict(color=color, dash=dash, width=2)),
                    row=1, col=1,
                )

            # P2 — breakdown
            for label, arr, color, dash in [
                ("Statistical", sig_stat,             _COLORS["stat"],  "solid"),
                ("Systematic",  res["sigma_P_sys"],   _COLORS["sys"],   "dash"),
                ("Total",       res["sigma_P_total"], _COLORS["total"], "dot"),
            ]:
                fig.add_trace(
                    go.Scatter(x=k_tgt, y=arr / Pk_tgt * 100, name=label,
                               line=dict(color=color, dash=dash, width=2),
                               showlegend=True),
                    row=1, col=2,
                )

            # P3 — fraction
            fig.add_trace(
                go.Scatter(
                    x=k_tgt, y=res["fraction_sys"] * 100,
                    fill="tozeroy", fillcolor="rgba(148,103,189,0.15)",
                    line=dict(color="mediumpurple", width=2),
                    name="Sys fraction", showlegend=False,
                ),
                row=2, col=1,
            )
            fig.add_hline(y=50, line_dash="dash", line_color="gray",
                          annotation_text="50%", row=2, col=1)

            # P4 — pie
            ms = {
                "Photo-z":  float(np.mean(res["sigma_P_photoz"]  ** 2)),
                "Bias":     float(np.mean(res["sigma_P_bias"]    ** 2)),
                "Geometry": float(np.mean(res["sigma_P_geometry"] ** 2)),
            }
            fig.add_trace(
                go.Pie(
                    labels=list(ms.keys()), values=list(ms.values()),
                    marker_colors=[_COLORS["photoz"], _COLORS["bias"], _COLORS["geom"]],
                    hole=0.35, showlegend=False,
                    textinfo="label+percent",
                ),
                row=2, col=2,
            )

            for r, c in [(1, 1), (1, 2), (2, 1)]:
                fig.update_xaxes(type="log", title_text="k [h/Mpc]", row=r, col=c)
            fig.update_yaxes(title_text="Relative error [%]",    row=1, col=1)
            fig.update_yaxes(title_text="Relative error [%]",    row=1, col=2)
            fig.update_yaxes(title_text="Systematic fraction [%]",
                             range=[0, 105], row=2, col=1)
            fig.update_layout(height=680, **_LAYOUT)
            st.plotly_chart(fig, use_container_width=True)

        except Exception as exc:
            st.error(f"Systematics analysis failed: {exc}")
            st.exception(exc)


# ══════════════════════════════════════════════════════════════════════════════
# Tab 3 — Power Spectrum Visualisation
# ══════════════════════════════════════════════════════════════════════════════
with tab_ps:
    st.header("Power Spectrum Visualisation")

    if not HAS_CAMB:
        st.error("CAMB required for this tab.")
    else:
        try:
            k_ps, _, Pk_ps = _power_spectrum(k_min, k_max, z_eff, H0, ombh2, omch2, ns)
            Pk_full = Pk_ps[0]

            params = dict(H0=H0, ombh2=ombh2, omch2=omch2,
                          As=2.1e-9, ns=ns, tau=0.0544)
            _m = PhiModulationModel(params=params)
            Pk_mod_full, mod_full = _m.apply_phi_modulation(
                k_ps, Pk_full, A_phi=A_phi,
            )

            mask = (k_ps >= k_min) & (k_ps <= k_max)
            k_v, Pb_v, Pm_v, mf_v = (
                k_ps[mask], Pk_full[mask], Pk_mod_full[mask], mod_full[mask],
            )

            fig = make_subplots(
                rows=3, cols=1,
                shared_xaxes=True,
                subplot_titles=[
                    f"P(k) at z = {z_eff:.1f}",
                    "φ-Modulation Factor  [1 + A_φ cos(2π ln(k/k₀)/ln φ)]",
                    "Relative Difference  [(P_mod − P_base) / P_base]",
                ],
                vertical_spacing=0.07,
                row_heights=[0.40, 0.30, 0.30],
            )

            # Row 1 — P(k)
            fig.add_trace(
                go.Scatter(x=k_v, y=Pb_v, name="ΛCDM baseline",
                           line=dict(color=_COLORS["lcdm"], width=2.5)),
                row=1, col=1,
            )
            fig.add_trace(
                go.Scatter(x=k_v, y=Pm_v, name="φ-modulated",
                           line=dict(color=_COLORS["mod"], width=2.5, dash="dash")),
                row=1, col=1,
            )

            # Row 2 — modulation factor
            fig.add_trace(
                go.Scatter(
                    x=k_v, y=mf_v, name="Mod factor",
                    line=dict(color="#2ca02c", width=2),
                    fill="tonexty",
                    fillcolor="rgba(44,160,44,0.08)",
                    showlegend=False,
                ),
                row=2, col=1,
            )
            fig.add_hline(y=1, line_dash="dash", line_color="gray",
                          opacity=0.6, row=2, col=1)
            fig.add_hrect(
                y0=1 - A_phi, y1=1 + A_phi,
                fillcolor="gray", opacity=0.08, line_width=0,
                row=2, col=1,
            )

            # Row 3 — relative difference
            ratio_pct = (Pm_v / Pb_v - 1) * 100
            fig.add_trace(
                go.Scatter(
                    x=k_v, y=ratio_pct, name="ΔP/P [%]",
                    line=dict(color=_COLORS["mod"], width=2),
                    fill="tozeroy", fillcolor="rgba(31,119,180,0.10)",
                    showlegend=False,
                ),
                row=3, col=1,
            )
            fig.add_hline(y=0, line_dash="dash", line_color="gray",
                          opacity=0.5, row=3, col=1)

            fig.update_xaxes(type="log", title_text="k [h/Mpc]", row=3, col=1)
            fig.update_yaxes(type="log", title_text="P(k) [(Mpc/h)³]", row=1, col=1)
            fig.update_yaxes(title_text="Modulation factor",  row=2, col=1)
            fig.update_yaxes(title_text="ΔP/P [%]",          row=3, col=1)
            fig.update_layout(
                height=780,
                legend=dict(orientation="h", yanchor="bottom", y=1.01, x=0),
                **_LAYOUT,
            )
            st.plotly_chart(fig, use_container_width=True)

        except Exception as exc:
            st.error(f"Power spectrum computation failed: {exc}")
            st.exception(exc)


# ══════════════════════════════════════════════════════════════════════════════
# Tab 4 — H(z) Analysis
# ══════════════════════════════════════════════════════════════════════════════
with tab_hz:
    st.header("H(z) Analysis — Falsified Background Model")
    st.markdown(
        "The **φ-recursive background** proposes "
        "H(z) ∝ ln(φ) × √(Ω_m(1+z)³ + Ω_Λ).  "
        "Both the forward (φ) and reverse (φ⁻¹) branches are fit against "
        "cosmic-chronometer data and compared with ΛCDM."
    )

    hz_path = os.path.join(ROOT, "real_hz.csv")
    if not os.path.exists(hz_path):
        st.error(f"`real_hz.csv` not found in `{ROOT}`. Run from the repo root.")
    else:
        try:
            z_d, h_d, s_d, res_f, res_r = _hz()
            om_f, t0_f = res_f.x
            om_r, t0_r = res_r.x

            # ΛCDM reference
            cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
            h_lcdm_d = cosmo.H(z_d).value
            chi2_lcdm = float(np.sum(((h_d - h_lcdm_d) / s_d) ** 2))

            c1, c2, c3 = st.columns(3)
            c1.metric("ΛCDM χ²",             f"{chi2_lcdm:.2f}")
            c2.metric("φ-Forward χ²",         f"{res_f.fun:.2f}")
            c3.metric("φ⁻¹-Reverse χ² (|H|)", f"{res_r.fun:.2f}")

            worst = min(res_f.fun, res_r.fun)
            if worst > chi2_lcdm:
                st.error(
                    f"Best φ-model is **worse** than ΛCDM by "
                    f"Δχ² = {worst - chi2_lcdm:.1f}  →  background model **falsified**."
                )
            st.markdown("---")

            z_fine = np.linspace(float(z_d.min()), float(z_d.max()), 400)
            h_lcdm_fine = cosmo.H(z_fine).value
            h_fwd_fine  = h_model(z_fine, om_f, t0_f, sigma=1.0)
            h_rev_fine  = np.abs(h_model(z_fine, om_r, t0_r, sigma=-1.0))

            fig = go.Figure()

            fig.add_trace(go.Scatter(
                x=z_d, y=h_d,
                error_y=dict(type="data", array=s_d, visible=True,
                             thickness=1.2, width=4),
                mode="markers",
                name="Cosmic chronometers",
                marker=dict(color="black", size=7, symbol="circle"),
            ))
            fig.add_trace(go.Scatter(
                x=z_fine, y=h_lcdm_fine,
                name=f"ΛCDM  (χ²={chi2_lcdm:.1f})",
                line=dict(color=_COLORS["lcdm"], width=2.5, dash="dot"),
            ))
            fig.add_trace(go.Scatter(
                x=z_fine, y=h_fwd_fine,
                name=f"φ-Forward  Ω_m={om_f:.3f}  (χ²={res_f.fun:.1f})",
                line=dict(color=_COLORS["fwd"], width=2.5),
            ))
            fig.add_trace(go.Scatter(
                x=z_fine, y=h_rev_fine,
                name=f"φ⁻¹-Reverse |…|  Ω_m={om_r:.3f}  (χ²={res_r.fun:.1f})",
                line=dict(color=_COLORS["rev"], width=2.5, dash="dash"),
            ))

            # Residual inset as annotation band
            fig.update_layout(
                xaxis_title="Redshift z",
                yaxis_title="H(z)  [km/s/Mpc]",
                title="H(z) Comparison: φ-Recursive Models vs ΛCDM",
                legend=dict(orientation="h", yanchor="bottom", y=1.02, x=0),
                height=500,
                **_LAYOUT,
            )
            st.plotly_chart(fig, use_container_width=True)

            # χ² residuals
            h_fwd_d = h_model(z_d, om_f, t0_f, sigma=1.0)
            h_rev_d = np.abs(h_model(z_d, om_r, t0_r, sigma=-1.0))

            fig2 = go.Figure()
            for label, h_model_d, color in [
                ("ΛCDM",       h_lcdm_d, _COLORS["lcdm"]),
                ("φ-Forward",  h_fwd_d,  _COLORS["fwd"]),
                ("φ⁻¹-Reverse",h_rev_d,  _COLORS["rev"]),
            ]:
                pull = (h_d - h_model_d) / s_d
                fig2.add_trace(go.Scatter(
                    x=z_d, y=pull, mode="markers+lines",
                    name=label,
                    marker=dict(color=color, size=6),
                    line=dict(color=color, width=1, dash="dot"),
                ))
            fig2.add_hline(y=0, line_dash="dash", line_color="gray", opacity=0.5)
            fig2.add_hrect(y0=-1, y1=1, fillcolor="gray", opacity=0.06, line_width=0)
            fig2.update_layout(
                xaxis_title="Redshift z",
                yaxis_title="Pull  (H_data − H_model) / σ",
                title="Normalised Residuals (Pull)",
                height=280,
                margin=dict(t=40),
                **_LAYOUT,
            )
            st.plotly_chart(fig2, use_container_width=True)

            with st.expander("Fit details"):
                st.write(
                    f"**φ-Forward:**  Ω_m = {om_f:.4f},  "
                    f"t₀ = {t0_f:.3e} s,  χ² = {res_f.fun:.2f}"
                )
                st.write(
                    f"**φ⁻¹-Reverse:**  Ω_m = {om_r:.4f},  "
                    f"t₀ = {t0_r:.3e} s,  χ² = {res_r.fun:.2f}"
                )
                st.write(f"**ΛCDM reference:**  χ² = {chi2_lcdm:.2f}")

        except Exception as exc:
            st.error(f"H(z) analysis failed: {exc}")
            st.exception(exc)


# ══════════════════════════════════════════════════════════════════════════════
# Tab 5 — CMB Residuals
# ══════════════════════════════════════════════════════════════════════════════
with tab_cmb:
    st.header("CMB Low-ℓ Residuals — Dual φ / φ⁻¹ Oscillation Fit")
    st.markdown(
        "A polynomial baseline in log ℓ approximates the ΛCDM spectrum.  "
        "Residuals Δ*C*_ℓ are fit with the **dual log-periodic model**:  \n"
        "ΔC_ℓ = A_φ cos(2π ln ℓ / ln φ + φ₀) + "
        "A_{φ⁻¹} cos(2π ln ℓ / ln(φ−1) + φ₁)"
    )

    cmb_path = os.path.join(ROOT, "real_cmb_lowl.csv")
    if not os.path.exists(cmb_path):
        st.error(f"`real_cmb_lowl.csv` not found in `{ROOT}`.")
    else:
        try:
            ell, cl, sig, cl_base, residuals, popt, pcov = _cmb()
            amp_phi, ph_phi, amp_conj, ph_conj = popt
            perr = np.sqrt(np.diag(pcov))

            c1, c2, c3, c4 = st.columns(4)
            c1.metric("A_φ",      f"{amp_phi:.2e}",  f"±{perr[0]:.1e}")
            c2.metric("φ₀  [rad]",f"{ph_phi:.3f}")
            c3.metric("A_{φ⁻¹}", f"{amp_conj:.2e}", f"±{perr[2]:.1e}")
            c4.metric("φ₁  [rad]",f"{ph_conj:.3f}")
            st.markdown("---")

            fig = make_subplots(
                rows=2, cols=1,
                shared_xaxes=True,
                subplot_titles=[
                    "Planck TT Low-ℓ  +  Polynomial ΛCDM Baseline",
                    "Residuals  +  Dual Log-Periodic Fit",
                ],
                vertical_spacing=0.10,
                row_heights=[0.45, 0.55],
            )

            # Row 1 — raw data + baseline
            fig.add_trace(go.Scatter(
                x=ell, y=cl,
                error_y=dict(type="data", array=sig, visible=True,
                             thickness=1, width=3, color="rgba(0,0,0,0.4)"),
                mode="markers",
                name="Planck TT",
                marker=dict(color="black", size=4),
            ), row=1, col=1)

            fig.add_trace(go.Scatter(
                x=ell, y=cl_base,
                name="Polynomial baseline",
                line=dict(color="firebrick", width=2, dash="dash"),
            ), row=1, col=1)

            # Row 2 — residuals + fit
            fig.add_trace(go.Scatter(
                x=ell, y=residuals,
                error_y=dict(type="data", array=sig, visible=True,
                             thickness=1, width=3,
                             color="rgba(31,119,180,0.4)"),
                mode="markers",
                name="Residuals ΔC_ℓ",
                marker=dict(color=_COLORS["mod"], size=4),
            ), row=2, col=1)

            ell_fine = np.exp(
                np.linspace(np.log(ell.min()), np.log(ell.max()), 600)
            )
            fit_fine = dual_osc_model(ell_fine, *popt)

            fig.add_trace(go.Scatter(
                x=ell_fine, y=fit_fine,
                name=f"Dual log-periodic fit  (A_φ={amp_phi:.1e}, A_{{φ⁻¹}}={amp_conj:.1e})",
                line=dict(color="crimson", width=2.5),
            ), row=2, col=1)

            # Individual mode contributions
            from cmb_analysis import _LN_PHI, _LN_PHI_CONJ  # type: ignore
            fit_phi  = amp_phi  * np.cos(2 * np.pi * np.log(ell_fine) / _LN_PHI  + ph_phi)
            fit_conj = amp_conj * np.cos(2 * np.pi * np.log(ell_fine) / _LN_PHI_CONJ + ph_conj)

            fig.add_trace(go.Scatter(
                x=ell_fine, y=fit_phi,
                name="φ mode only",
                line=dict(color=_COLORS["fwd"], width=1.5, dash="dot"),
                opacity=0.7,
            ), row=2, col=1)
            fig.add_trace(go.Scatter(
                x=ell_fine, y=fit_conj,
                name="φ⁻¹ mode only",
                line=dict(color=_COLORS["rev"], width=1.5, dash="dot"),
                opacity=0.7,
            ), row=2, col=1)

            fig.add_hline(y=0, line_dash="solid", line_color="gray",
                          line_width=0.5, opacity=0.5, row=2, col=1)

            fig.update_xaxes(type="log", title_text="Multipole ℓ", row=2, col=1)
            fig.update_yaxes(title_text="ℓ(ℓ+1)C_ℓ / 2π  [μK²]", row=1, col=1)
            fig.update_yaxes(title_text="ΔC_ℓ  [μK²]",            row=2, col=1)
            fig.update_layout(
                height=680,
                legend=dict(orientation="h", yanchor="bottom", y=1.01, x=0),
                **_LAYOUT,
            )
            st.plotly_chart(fig, use_container_width=True)

        except Exception as exc:
            st.error(f"CMB analysis failed: {exc}")
            st.exception(exc)


# ══════════════════════════════════════════════════════════════════════════════
# Tab 6 — LSS P(k) φ-Scale Detection
# ══════════════════════════════════════════════════════════════════════════════
with tab_lss:
    st.header("LSS P(k) — φ-Scale Feature Detection")
    st.markdown(
        "A smooth spline baseline is subtracted to isolate oscillatory residuals.  "
        "Dashed red lines mark the expected φ-spaced scales "
        "k_n = φⁿ × k_BAO."
    )

    pk_path = os.path.join(ROOT, "real_pk_lowk.csv")
    if not os.path.exists(pk_path):
        st.error(f"`real_pk_lowk.csv` not found in `{ROOT}`.")
    else:
        try:
            k_d, pk_d, s_d, phi_scales, pk_sm, res_pk, k_pk, k_rp = _lss()

            def _phi_matches(k_found, tol=0.005):
                if len(k_found) == 0:
                    return 0
                return sum(
                    np.min(np.abs(k_found - ks)) < tol
                    for ks in phi_scales
                )

            c1, c2, c3, c4 = st.columns(4)
            c1.metric("Direct peaks",       str(len(k_pk)))
            c2.metric("Residual peaks",     str(len(k_rp)))
            c3.metric("φ-scale matches",
                       f"{_phi_matches(k_rp)} / {len(phi_scales)}")
            c4.metric("φ  (golden ratio)",  f"{_PHI:.6f}")
            st.markdown("---")

            fig = make_subplots(
                rows=2, cols=1,
                shared_xaxes=True,
                subplot_titles=[
                    "Matter Power Spectrum  P(k)  +  Smooth Baseline",
                    "Oscillation Residuals  ΔP(k)",
                ],
                vertical_spacing=0.10,
                row_heights=[0.50, 0.50],
            )

            # Row 1 — P(k)
            fig.add_trace(go.Scatter(
                x=k_d, y=pk_d,
                error_y=dict(type="data", array=s_d, visible=True,
                             thickness=1, width=3,
                             color="rgba(31,119,180,0.35)"),
                mode="markers",
                name="P(k) data",
                marker=dict(color=_COLORS["mod"], size=5),
            ), row=1, col=1)

            fig.add_trace(go.Scatter(
                x=k_d, y=pk_sm,
                name="Spline baseline",
                line=dict(color="#2ca02c", width=2.5),
            ), row=1, col=1)

            if len(k_pk) > 0:
                idx = np.searchsorted(k_d, k_pk).clip(0, len(k_d) - 1)
                fig.add_trace(go.Scatter(
                    x=k_pk, y=pk_d[idx],
                    mode="markers", name="Direct peaks",
                    marker=dict(color="#ff7f0e", size=13, symbol="star",
                                line=dict(color="darkorange", width=1)),
                ), row=1, col=1)

            # Row 2 — residuals
            fig.add_trace(go.Scatter(
                x=k_d, y=res_pk,
                error_y=dict(type="data", array=s_d, visible=True,
                             thickness=1, width=3,
                             color="rgba(31,119,180,0.30)"),
                mode="markers",
                name="Residuals",
                marker=dict(color=_COLORS["mod"], size=4),
            ), row=2, col=1)

            fig.add_hline(y=0, line_dash="solid", line_color="black",
                          line_width=0.5, row=2, col=1)

            if len(k_rp) > 0:
                ri = np.searchsorted(k_d, k_rp).clip(0, len(k_d) - 1)
                fig.add_trace(go.Scatter(
                    x=k_rp, y=res_pk[ri],
                    mode="markers", name="Residual peaks",
                    marker=dict(color="royalblue", size=13, symbol="star",
                                line=dict(color="darkblue", width=1)),
                ), row=2, col=1)

            # φ-scale vertical lines in both rows
            for i, ks in enumerate(phi_scales):
                if k_d.min() <= ks <= k_d.max():
                    label = f"φ^{i - 5}·k_BAO"
                    for row in (1, 2):
                        fig.add_vline(
                            x=ks,
                            line_dash="dash",
                            line_color=_COLORS["phi_line"],
                            line_width=1.5,
                            annotation_text=label if row == 1 else "",
                            annotation_font_size=9,
                            annotation_font_color="firebrick",
                            row=row, col=1,
                        )

            fig.update_xaxes(type="log", title_text="k [h/Mpc]", row=2, col=1)
            fig.update_yaxes(type="log", title_text="P(k) [(Mpc/h)³]", row=1, col=1)
            fig.update_yaxes(title_text="ΔP(k) [(Mpc/h)³]",            row=2, col=1)
            fig.update_layout(
                height=720,
                legend=dict(orientation="h", yanchor="bottom", y=1.01, x=0),
                **_LAYOUT,
            )
            st.plotly_chart(fig, use_container_width=True)

            # φ-scale table
            with st.expander("Expected φ-scale reference table"):
                import pandas as pd
                df = pd.DataFrame({
                    "n":           list(range(-5, 6)),
                    "k_φⁿ [h/Mpc]": [float(f"{s:.6f}") for s in phi_scales],
                    "In data range": [
                        "✓" if k_d.min() <= s <= k_d.max() else "—"
                        for s in phi_scales
                    ],
                })
                st.dataframe(df, use_container_width=True, hide_index=True)

        except Exception as exc:
            st.error(f"LSS analysis failed: {exc}")
            st.exception(exc)


# ══════════════════════════════════════════════════════════════════════════════
# Tab 7 — Summary
# ══════════════════════════════════════════════════════════════════════════════
with tab_summary:
    st.header("Analysis Summary")

    col_l, col_r = st.columns([1, 1])

    with col_l:
        st.markdown("### Active Parameters")
        params_md = f"""
| Parameter | Value |
|-----------|-------|
| A_φ | {A_phi:.4f} |
| z_eff | {z_eff:.2f} |
| k range | {k_min:.3f} – {k_max:.2f} h/Mpc |
| H₀ | {H0:.2f} km/s/Mpc |
| Ω_b h² | {ombh2:.4f} |
| Ω_c h² | {omch2:.3f} |
| n_s | {ns:.4f} |
"""
        st.markdown(params_md)

    with col_r:
        st.markdown("### Scientific Status")
        st.error("❌  Background model **falsified**  (H(z), CMB, P(k))")
        st.success("✅  Perturbation model **testable** with near-future surveys")
        st.info("🔭  DESI Year-5 constrains A_φ to the forecast σ_Aφ level")

    if HAS_CAMB:
        st.markdown("---")
        st.markdown("### Forecast at current parameters")
        try:
            fc_s = _forecast(A_phi, k_min, k_max, H0, ombh2, omch2, ns, include_sys)
            SNR_s = fc_s["SNR"]
            s_s   = fc_s.get("sigma_Aphi_total", fc_s["sigma_Aphi"])

            col1, col2 = st.columns(2)
            col1.metric("Signal-to-noise", f"{SNR_s:.2f} σ")
            col2.metric("Total σ_Aφ",      f"{s_s:.2e}")

            if SNR_s >= 5:
                st.success(f"Very strong detection: {SNR_s:.2f}σ  (≥ 5σ)")
            elif SNR_s >= 3:
                st.success(f"Strong detection: {SNR_s:.2f}σ  (≥ 3σ)")
            elif SNR_s >= 2:
                st.warning(f"Marginal detection: {SNR_s:.2f}σ  (2–3σ)")
            else:
                st.info(f"Below detection threshold: {SNR_s:.2f}σ  (< 2σ)")
        except Exception:
            st.warning("Could not compute forecast (check CAMB installation).")

    st.markdown("---")
    st.caption(
        "**FaCC** — Fibonacci as Cosmological Constant  ·  "
        "Bryan David Persaud  ·  "
        "[GitHub](https://github.com/imediacorp/FaCC)"
    )
