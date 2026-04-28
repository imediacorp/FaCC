"""
G–φ Connection Analysis.

Three independent angles on whether the gravitational constant G
and the golden ratio φ are related:

1. Varying-G cosmology (Brans-Dicke inspired):
   G(t) = G_0 × [1 + A_G cos(2π ln(t/t_0) / ln φ + δ)]
   Constrained by binary pulsar Ġ/G ≲ 10⁻¹² yr⁻¹.

2. GW chirp-frequency φ-modulation:
   Since f_ISCO ∝ 1/(G M), a φ-spacing in BH masses produces
   φ-spaced ISCO frequencies.  We model the stochastic GW
   background (SGWB) Ω_GW(f) and overlay φ-harmonic lines.

3. GWTC-3 mass-ratio clustering:
   Test whether GW event mass ratios q = m1/m2 cluster near
   φ^n values using a KS test against uniform log(q) and
   a histogram binned at log(φ) intervals.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from scipy.stats import kstest

_PHI = (1 + np.sqrt(5)) / 2
_LN_PHI = np.log(_PHI)
_G0 = 6.674e-11          # m^3 kg^-1 s^-2
_C = 2.998e8              # m/s
_MSUN_KG = 1.989e30       # kg
_PC_M = 3.086e16          # m per parsec

# ─── 1. Varying-G model ───────────────────────────────────────────────────────

def varying_g_model(
    t: NDArray,
    A_G: float = 1e-5,
    t0: float = 4.35e17,
    delta: float = 0.0,
) -> NDArray:
    """φ-log-periodic varying-G model.

    G(t) = G₀ [1 + A_G cos(2π ln(t/t₀) / ln φ + δ)]

    Parameters
    ----------
    t : array
        Cosmic time [s].
    A_G : float
        Modulation amplitude (dimensionless).
    t0 : float
        Reference time scale [s].  Default ≈ age of universe.
    delta : float
        Phase offset [rad].

    Returns
    -------
    G_t : array
        G(t) in m³ kg⁻¹ s⁻².
    """
    phase = 2 * np.pi * np.log(t / t0) / _LN_PHI + delta
    return _G0 * (1.0 + A_G * np.cos(phase))


def gdot_over_g(
    t: NDArray,
    A_G: float = 1e-5,
    t0: float = 4.35e17,
    delta: float = 0.0,
) -> NDArray:
    """Compute Ġ/G for the varying-G model [yr⁻¹].

    Ġ/G = −A_G (2π / ln φ) sin(2π ln(t/t₀) / ln φ + δ) / t
    """
    S_PER_YR = 3.15576e7
    phase = 2 * np.pi * np.log(t / t0) / _LN_PHI + delta
    gdot_g_per_s = -A_G * (2 * np.pi / _LN_PHI) * np.sin(phase) / t
    return gdot_g_per_s * S_PER_YR


def pulsar_constraint_amplitude(t_obs: float = 4.35e17) -> float:
    """Maximum A_G allowed by pulsar Ġ/G < 10⁻¹² yr⁻¹ at present epoch.

    From |Ġ/G| ≤ 10⁻¹² yr⁻¹:
        A_G ≤ t_obs × 10⁻¹² × ln φ / (2π × 3.15576×10⁷)
    """
    S_PER_YR = 3.15576e7
    gdot_limit_per_s = 1e-12 / S_PER_YR
    return gdot_limit_per_s * t_obs * _LN_PHI / (2 * np.pi)


# ─── 2. GW chirp-frequency φ-modulation ──────────────────────────────────────

def f_isco(m_total_msun: float | NDArray, G: float = _G0) -> NDArray:
    """ISCO gravitational-wave frequency for a Schwarzschild BH [Hz].

    f_ISCO = c³ / (6^(3/2) π G M)

    Parameters
    ----------
    m_total_msun : float or array
        Total binary mass in solar masses.
    G : float
        Gravitational constant [m³ kg⁻¹ s⁻²].

    Returns
    -------
    f_ISCO : array [Hz]
    """
    m_kg = np.asarray(m_total_msun, dtype=float) * _MSUN_KG
    return _C ** 3 / (6.0 ** 1.5 * np.pi * G * m_kg)


def phi_harmonic_masses(
    m_ref: float = 30.0, n_range: tuple[int, int] = (-4, 5)
) -> NDArray:
    """Generate a series of masses spaced by φ around a reference mass.

    m_n = m_ref × φ^n   for n in range(*n_range)
    """
    return m_ref * np.array([_PHI ** n for n in range(*n_range)])


def omega_gw_powerlaw(
    f: NDArray,
    omega_ref: float = 1e-9,
    f_ref: float = 25.0,
    alpha: float = 2.0 / 3.0,
) -> NDArray:
    """SGWB energy-density spectrum (power-law background).

    Ω_GW(f) = Ω_ref × (f / f_ref)^α

    For CBC background: α = 2/3 (inspiral dominated).
    """
    return omega_ref * (f / f_ref) ** alpha


def omega_gw_phi_modulated(
    f: NDArray,
    omega_ref: float = 1e-9,
    f_ref: float = 25.0,
    alpha: float = 2.0 / 3.0,
    A_phi: float = 0.05,
    phase: float = 0.0,
) -> NDArray:
    """SGWB with superimposed φ-log-periodic modulation.

    Ω_GW(f) = Ω_ref (f/f_ref)^α [1 + A_φ cos(2π ln(f/f_ref) / ln φ + δ)]
    """
    base = omega_gw_powerlaw(f, omega_ref, f_ref, alpha)
    mod = 1.0 + A_phi * np.cos(2 * np.pi * np.log(f / f_ref) / _LN_PHI + phase)
    return base * mod


# ─── 3. GWTC mass-ratio clustering ───────────────────────────────────────────

GWTC3_MASSES: list[tuple[str, float, float]] = [
    # (event, m1_source, m2_source) — GWTC-3 median values [M_sun]
    ("GW150914",  35.6, 30.6),
    ("GW151012",  23.3, 13.6),
    ("GW151226",  13.7,  7.7),
    ("GW170104",  31.0, 20.1),
    ("GW170608",  11.0,  7.6),
    ("GW170729",  50.6, 34.3),
    ("GW170809",  35.2, 23.8),
    ("GW170814",  30.7, 25.3),
    ("GW170818",  35.5, 26.8),
    ("GW170823",  39.6, 29.4),
    ("GW190408_181802", 24.6, 18.4),
    ("GW190412",  30.1,  8.3),
    ("GW190421_213856", 42.0, 32.7),
    ("GW190503_185404", 44.0, 31.6),
    ("GW190512_180714", 23.2, 12.2),
    ("GW190513_205428", 35.6, 17.8),
    ("GW190517_055101", 37.4, 25.2),
    ("GW190519_153544", 66.0, 40.5),
    ("GW190521",        85.1, 65.4),
    ("GW190602_175927", 69.1, 47.8),
    ("GW190620_030421", 57.1, 35.4),
    ("GW190630_185205", 35.1, 23.6),
    ("GW190701_203306", 53.9, 40.0),
    ("GW190706_222641", 67.0, 38.2),
    ("GW190707_093326",  12.0,  7.6),
    ("GW190708_232457",  17.1, 11.9),
    ("GW190720_000836",  13.0,  5.8),
    ("GW190727_060333",  38.1, 26.0),
    ("GW190728_064510",  10.8,  6.1),
    ("GW190731_140936",  41.5, 28.8),
    ("GW190803_022701",  37.7, 27.0),
    ("GW190814",         23.2,  2.6),
    ("GW190828_063405",  32.1, 26.5),
    ("GW190828_065509",  24.0, 10.1),
    ("GW190910_112807",  44.8, 32.8),
    ("GW190915_235702",  35.9, 24.3),
    ("GW190924_021846",   8.9,  5.0),
    ("GW190929_012149",  80.8, 23.1),
    ("GW190930_133541",  12.3,  7.6),
    ("GW191103_012549",  11.4,  7.5),
    ("GW191105_143521",  10.7,  7.7),
    ("GW191109_010717",  65.1, 47.4),
    ("GW191113_071753",   6.7,  1.9),
    ("GW191126_115259",  12.0,  6.2),
    ("GW191127_050227",  13.6,  7.6),
    ("GW191129_134029",  10.7,  6.7),
    ("GW191204_110529",  11.8,  5.2),
    ("GW191204_171526",  11.9,  8.3),
    ("GW191215_223052",  24.0, 13.4),
    ("GW191216_213338",  12.1,  7.7),
    ("GW191222_033537",  47.6, 34.1),
    ("GW191230_180458",  64.7, 46.6),
    ("GW200105_162426",   8.9,  1.9),
    ("GW200112_155838",  47.6, 36.1),
    ("GW200115_042309",   5.7,  1.5),
    ("GW200128_022011",  40.1, 26.8),
    ("GW200129_065458",  34.5, 28.9),
    ("GW200202_154313",  10.1,  7.3),
    ("GW200208_130117",  38.6, 27.4),
    ("GW200208_222617",  56.1, 37.9),
    ("GW200209_085452",  46.3, 26.0),
    ("GW200210_092254",  24.1,  2.8),
    ("GW200216_220804",  36.2, 23.7),
    ("GW200219_094415",  40.0, 26.1),
    ("GW200220_061928",  87.0, 57.4),
    ("GW200220_124850",  28.1, 10.1),
    ("GW200224_222234",  40.2, 32.8),
    ("GW200225_060421",  19.3, 13.4),
    ("GW200302_015811",  18.8, 11.9),
    ("GW200306_093714",  16.1,  8.3),
    ("GW200308_173609",  76.5, 45.0),
    ("GW200311_115853",  43.4, 32.2),
    ("GW200316_215756",  13.1,  7.8),
    ("GW200322_091133",  49.6, 28.4),
]


def load_gwtc3_builtin() -> tuple[NDArray, NDArray, list[str]]:
    """Return built-in GWTC-3 source masses (no file required).

    Returns
    -------
    m1, m2 : arrays of float  [M_sun]
    events : list of str
    """
    events = [row[0] for row in GWTC3_MASSES]
    m1 = np.array([row[1] for row in GWTC3_MASSES])
    m2 = np.array([row[2] for row in GWTC3_MASSES])
    return m1, m2, events


def compute_mass_ratios(m1: NDArray, m2: NDArray) -> NDArray:
    """Compute q = m1/m2 (always ≥ 1)."""
    return np.maximum(m1, m2) / np.minimum(m1, m2)


def phi_expected_log_ratios(n_range: tuple[int, int] = (-2, 4)) -> NDArray:
    """Expected ln(φ^n) values for integer n in n_range."""
    return np.array([n * _LN_PHI for n in range(*n_range)])


def ks_test_phi_clustering(log_q: NDArray) -> tuple[float, float]:
    """KS test: is log(q) consistent with uniform distribution?

    Under the null hypothesis that mass ratios carry no preferred scale,
    log(q) is approximately uniform in [0, log(q_max)].

    Returns
    -------
    ks_stat, p_value
    """
    q_min, q_max = log_q.min(), log_q.max()
    cdf_uniform = lambda x: (x - q_min) / (q_max - q_min)
    ks_stat, p_value = kstest(log_q, cdf_uniform)
    return float(ks_stat), float(p_value)


def count_phi_ratio_matches(
    log_q: NDArray,
    n_range: tuple[int, int] = (-2, 4),
    tolerance: float = 0.15,
) -> tuple[int, NDArray]:
    """Count mass ratios within `tolerance` of a φ^n value (in log space).

    Parameters
    ----------
    log_q : array
        Natural log of mass ratios.
    tolerance : float
        Half-window in log-ratio units (~15% by default).

    Returns
    -------
    n_matches : int
    phi_lines : array of log(φ^n) reference positions
    """
    phi_lines = phi_expected_log_ratios(n_range)
    n_matches = int(sum(
        np.min(np.abs(log_q - lp)) < tolerance for lp in phi_lines
    ))
    return n_matches, phi_lines


def run_g_phi_analysis() -> dict:
    """Run all three G–φ analyses and return a results dict.

    Returns
    -------
    dict with keys:
        t_arr, G_arr, gdot_arr, A_G_max
        f_arr, omega_base, omega_mod, phi_masses, f_isco_arr
        m1, m2, events, q, log_q, ks_stat, p_value, n_matches, phi_lines
    """
    # 1. Varying G
    t_now = 4.35e17          # s (age of universe)
    t_arr = np.linspace(0.1 * t_now, t_now, 500)
    A_G_max = pulsar_constraint_amplitude(t_now)
    G_arr = varying_g_model(t_arr, A_G=A_G_max)
    gdot_arr = gdot_over_g(t_arr, A_G=A_G_max)

    # 2. GW frequency
    f_arr = np.logspace(-2, 3, 800)       # 0.01–1000 Hz
    omega_base = omega_gw_powerlaw(f_arr)
    omega_mod  = omega_gw_phi_modulated(f_arr)
    phi_masses = phi_harmonic_masses()
    f_isco_arr = f_isco(phi_masses)

    # 3. Mass-ratio clustering
    m1, m2, events = load_gwtc3_builtin()
    q = compute_mass_ratios(m1, m2)
    log_q = np.log(q)
    ks_stat, p_value = ks_test_phi_clustering(log_q)
    n_matches, phi_lines = count_phi_ratio_matches(log_q)

    return {
        # Varying G
        "t_arr": t_arr,
        "G_arr": G_arr,
        "gdot_arr": gdot_arr,
        "A_G_max": A_G_max,
        # GW spectrum
        "f_arr": f_arr,
        "omega_base": omega_base,
        "omega_mod": omega_mod,
        "phi_masses": phi_masses,
        "f_isco_arr": f_isco_arr,
        # Mass ratios
        "m1": m1,
        "m2": m2,
        "events": events,
        "q": q,
        "log_q": log_q,
        "ks_stat": ks_stat,
        "p_value": p_value,
        "n_matches": n_matches,
        "phi_lines": phi_lines,
    }


if __name__ == "__main__":
    res = run_g_phi_analysis()
    print(f"Max allowed A_G (pulsar constraint): {res['A_G_max']:.2e}")
    print(f"GWTC-3 KS stat: {res['ks_stat']:.3f}, p={res['p_value']:.3f}")
    print(f"φ-ratio matches (±15% log window): {res['n_matches']}")
    print(f"φ-spaced ISCO frequencies [Hz]: {res['f_isco_arr'].round(2)}")
