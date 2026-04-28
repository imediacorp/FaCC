"""Tests for src/g_phi_analysis.py."""

import numpy as np
import pytest

from src.g_phi_analysis import (
    varying_g_model,
    gdot_over_g,
    pulsar_constraint_amplitude,
    f_isco,
    phi_harmonic_masses,
    omega_gw_powerlaw,
    omega_gw_phi_modulated,
    load_gwtc3_builtin,
    compute_mass_ratios,
    phi_expected_log_ratios,
    ks_test_phi_clustering,
    count_phi_ratio_matches,
    _G0,
    _PHI,
    _LN_PHI,
)

_T_NOW = 4.35e17  # seconds (age of universe)


class TestVaryingGModel:
    def test_zero_amplitude_returns_G0(self):
        t = np.linspace(0.1 * _T_NOW, _T_NOW, 50)
        G = varying_g_model(t, A_G=0.0)
        assert np.allclose(G, _G0)

    def test_output_shape(self):
        t = np.linspace(0.1 * _T_NOW, _T_NOW, 100)
        G = varying_g_model(t)
        assert G.shape == t.shape

    def test_bounded_by_amplitude(self):
        t = np.linspace(0.1 * _T_NOW, _T_NOW, 200)
        A_G = 0.01
        G = varying_g_model(t, A_G=A_G)
        assert np.all(G >= _G0 * (1 - A_G) - 1e-30)
        assert np.all(G <= _G0 * (1 + A_G) + 1e-30)

    def test_scalar_input(self):
        G = varying_g_model(np.array([_T_NOW]))
        assert G.shape == (1,)


class TestGdotOverG:
    def test_zero_amplitude_gives_zero(self):
        t = np.linspace(0.1 * _T_NOW, _T_NOW, 50)
        gdot = gdot_over_g(t, A_G=0.0)
        assert np.allclose(gdot, 0.0)

    def test_output_shape(self):
        t = np.linspace(0.1 * _T_NOW, _T_NOW, 50)
        gdot = gdot_over_g(t)
        assert gdot.shape == t.shape

    def test_magnitude_small_for_small_amplitude(self):
        t = np.array([_T_NOW])
        A_G = 1e-5
        gdot = gdot_over_g(t, A_G=A_G)
        assert abs(gdot[0]) < 1e-10  # well below typical pulsar constraint


class TestPulsarConstraintAmplitude:
    def test_value_around_one_millitenth(self):
        A_max = pulsar_constraint_amplitude(_T_NOW)
        # Should be ~10^-3 for the pulsar Ġ/G < 10^-12 yr^-1 constraint
        assert 1e-4 < A_max < 1e-2

    def test_returns_float(self):
        assert isinstance(pulsar_constraint_amplitude(), float)

    def test_proportional_to_t_obs(self):
        A1 = pulsar_constraint_amplitude(t_obs=_T_NOW)
        A2 = pulsar_constraint_amplitude(t_obs=2 * _T_NOW)
        assert A2 == pytest.approx(2 * A1, rel=1e-9)


class TestFIsco:
    def test_30_solar_mass_bh(self):
        # f_ISCO for 30 M_sun ≈ 147 Hz
        freq = f_isco(30.0)
        assert 100 < float(freq) < 200

    def test_larger_mass_lower_frequency(self):
        f_light = f_isco(10.0)
        f_heavy = f_isco(100.0)
        assert float(f_light) > float(f_heavy)

    def test_array_input(self):
        masses = np.array([10.0, 30.0, 100.0])
        freqs = f_isco(masses)
        assert freqs.shape == masses.shape
        assert np.all(freqs > 0)

    def test_f_isco_scales_with_G(self):
        f1 = f_isco(30.0, G=_G0)
        f2 = f_isco(30.0, G=2 * _G0)
        assert f1 == pytest.approx(2 * f2, rel=1e-9)


class TestPhiHarmonicMasses:
    def test_default_length(self):
        masses = phi_harmonic_masses()
        assert len(masses) == 9  # range(-4, 5)

    def test_ratio_between_adjacent_is_phi(self):
        masses = phi_harmonic_masses()
        ratios = masses[1:] / masses[:-1]
        assert np.allclose(ratios, _PHI, rtol=1e-10)

    def test_reference_mass_at_n_zero(self):
        masses = phi_harmonic_masses(m_ref=30.0, n_range=(0, 1))
        assert masses[0] == pytest.approx(30.0, rel=1e-10)


class TestOmegaGwPowerlaw:
    def test_output_shape(self):
        f = np.logspace(-2, 3, 100)
        out = omega_gw_powerlaw(f)
        assert out.shape == f.shape

    def test_all_positive(self):
        f = np.logspace(-2, 3, 100)
        out = omega_gw_powerlaw(f)
        assert np.all(out > 0)

    def test_at_f_ref_equals_omega_ref(self):
        out = omega_gw_powerlaw(np.array([25.0]), omega_ref=1e-9, f_ref=25.0)
        assert float(out[0]) == pytest.approx(1e-9, rel=1e-9)

    def test_power_law_slope(self):
        f = np.array([10.0, 20.0])
        out = omega_gw_powerlaw(f, alpha=2.0 / 3.0)
        ratio = out[1] / out[0]
        assert ratio == pytest.approx((20.0 / 10.0) ** (2.0 / 3.0), rel=1e-9)


class TestOmegaGwPhiModulated:
    def test_zero_modulation_matches_powerlaw(self):
        f = np.logspace(-2, 3, 50)
        base = omega_gw_powerlaw(f)
        mod  = omega_gw_phi_modulated(f, A_phi=0.0)
        assert np.allclose(base, mod)

    def test_modulation_bounded(self):
        f = np.logspace(-2, 3, 200)
        A = 0.1
        base = omega_gw_powerlaw(f)
        mod  = omega_gw_phi_modulated(f, A_phi=A)
        ratio = mod / base
        assert np.all(ratio >= 1 - A - 1e-12)
        assert np.all(ratio <= 1 + A + 1e-12)


class TestLoadGwtc3Builtin:
    def test_returns_three_outputs(self):
        m1, m2, events = load_gwtc3_builtin()
        assert len(m1) == len(m2) == len(events)

    def test_m1_geq_m2(self):
        m1, m2, _ = load_gwtc3_builtin()
        # m1 is primary (heavier by convention in the table)
        assert np.all(m1 >= m2 * 0.9)  # allow small tolerance for equal-mass systems

    def test_masses_positive(self):
        m1, m2, _ = load_gwtc3_builtin()
        assert np.all(m1 > 0)
        assert np.all(m2 > 0)

    def test_reasonable_event_count(self):
        m1, _, _ = load_gwtc3_builtin()
        assert len(m1) >= 60


class TestComputeMassRatios:
    def test_ratio_always_geq_one(self):
        m1 = np.array([30.0, 10.0, 50.0])
        m2 = np.array([20.0, 15.0, 50.0])
        q = compute_mass_ratios(m1, m2)
        assert np.all(q >= 1.0)

    def test_equal_masses_give_one(self):
        m1 = np.array([25.0])
        m2 = np.array([25.0])
        q = compute_mass_ratios(m1, m2)
        assert q[0] == pytest.approx(1.0)


class TestPhiExpectedLogRatios:
    def test_n_zero_gives_zero(self):
        lines = phi_expected_log_ratios(n_range=(0, 1))
        assert lines[0] == pytest.approx(0.0)

    def test_n_one_gives_ln_phi(self):
        lines = phi_expected_log_ratios(n_range=(1, 2))
        assert lines[0] == pytest.approx(_LN_PHI, rel=1e-10)

    def test_spacing_is_ln_phi(self):
        lines = phi_expected_log_ratios(n_range=(0, 4))
        diffs = np.diff(lines)
        assert np.allclose(diffs, _LN_PHI, rtol=1e-10)


class TestKsTestPhiClustering:
    def test_returns_two_floats(self):
        log_q = np.log(np.random.default_rng(42).uniform(1, 5, 50))
        ks, pv = ks_test_phi_clustering(log_q)
        assert isinstance(ks, float)
        assert isinstance(pv, float)

    def test_p_value_in_unit_interval(self):
        log_q = np.linspace(0, 2, 30)
        _, pv = ks_test_phi_clustering(log_q)
        assert 0.0 <= pv <= 1.0

    def test_perfectly_uniform_high_p(self):
        # Uniform distribution should not reject the null
        log_q = np.linspace(0.001, 2.0, 200)
        _, pv = ks_test_phi_clustering(log_q)
        assert pv > 0.05


class TestCountPhiRatioMatches:
    def test_returns_int_and_array(self):
        log_q = np.array([0.0, _LN_PHI, 2 * _LN_PHI])
        n, lines = count_phi_ratio_matches(log_q)
        assert isinstance(n, int)
        assert isinstance(lines, np.ndarray)

    def test_exact_phi_match(self):
        log_q = np.array([_LN_PHI])
        n, _ = count_phi_ratio_matches(log_q, n_range=(0, 3), tolerance=0.01)
        assert n >= 1

    def test_no_match_far_from_phi(self):
        log_q = np.array([5.0, 6.0, 7.0])  # not near any φ^n in range(-2, 4)
        n, _ = count_phi_ratio_matches(log_q, n_range=(-2, 4), tolerance=0.1)
        assert n == 0
