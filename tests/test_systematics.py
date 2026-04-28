"""Tests for src/systematics.py."""

import numpy as np
import pytest

from src.systematics import SystematicErrorBudget


K = np.logspace(-2, -0.5, 50)
PK = 1e4 * (K / 0.05) ** (-1.5)
SIGMA_STAT = PK * 0.05


class TestSystematicErrorBudgetInit:
    def test_default_z_eff(self):
        s = SystematicErrorBudget()
        assert s.z_eff == pytest.approx(0.8)

    def test_custom_z_eff(self):
        s = SystematicErrorBudget(z_eff=1.2)
        assert s.z_eff == pytest.approx(1.2)

    def test_photo_z_error_increases_with_z(self):
        s_low = SystematicErrorBudget(z_eff=0.5)
        s_high = SystematicErrorBudget(z_eff=1.5)
        assert s_high.sigma_z_photo > s_low.sigma_z_photo


class TestPhotoZError:
    def setup_method(self):
        self.s = SystematicErrorBudget(z_eff=0.8)

    def test_output_shape(self):
        err = self.s.photo_z_error(K, PK)
        assert err.shape == K.shape

    def test_non_negative(self):
        err = self.s.photo_z_error(K, PK)
        assert np.all(err >= 0)

    def test_relative_error_capped_at_10_percent(self):
        # The implementation clips relative error at 0.1 (10%)
        err = self.s.photo_z_error(K, PK)
        rel_err = err / PK
        assert np.all(rel_err <= 0.1 + 1e-12)

    def test_zero_Pk_gives_zero_error(self):
        err = self.s.photo_z_error(K, np.zeros_like(K))
        np.testing.assert_allclose(err, 0.0)

    def test_custom_sigma_z(self):
        err_small = self.s.photo_z_error(K, PK, sigma_z=0.01)
        err_large = self.s.photo_z_error(K, PK, sigma_z=0.05)
        assert np.all(err_large >= err_small)


class TestBiasUncertainty:
    def setup_method(self):
        self.s = SystematicErrorBudget(z_eff=0.8)

    def test_output_shape(self):
        err = self.s.bias_uncertainty(K, PK)
        assert err.shape == K.shape

    def test_non_negative(self):
        err = self.s.bias_uncertainty(K, PK)
        assert np.all(err >= 0)

    def test_proportional_to_Pk(self):
        err = self.s.bias_uncertainty(K, PK)
        ratio = err / PK
        np.testing.assert_allclose(ratio, ratio[0], rtol=1e-10)

    def test_scales_with_sigma_b(self):
        err_small = self.s.bias_uncertainty(K, PK, sigma_b=0.02)
        err_large = self.s.bias_uncertainty(K, PK, sigma_b=0.10)
        np.testing.assert_allclose(err_large / err_small, 5.0, rtol=1e-10)


class TestSurveyGeometryError:
    def setup_method(self):
        self.s = SystematicErrorBudget(z_eff=0.8)

    def test_output_shape(self):
        err = self.s.survey_geometry_error(K, PK)
        assert err.shape == K.shape

    def test_non_negative(self):
        err = self.s.survey_geometry_error(K, PK)
        assert np.all(err >= 0)

    def test_larger_survey_reduces_error(self):
        err_small = self.s.survey_geometry_error(K, PK, V_survey=10.0)
        err_large = self.s.survey_geometry_error(K, PK, V_survey=200.0)
        # Larger survey → larger k_min_survey → smaller suppression at given k → less error at low k
        assert err_small[0] >= err_large[0]


class TestComputeSystematicBudget:
    def setup_method(self):
        self.s = SystematicErrorBudget(z_eff=0.8)

    def test_output_keys(self):
        result = self.s.compute_systematic_budget(K, PK, SIGMA_STAT)
        for key in ("sigma_P_total", "sigma_P_sys", "sigma_P_photoz",
                    "sigma_P_bias", "sigma_P_geometry", "fraction_sys"):
            assert key in result

    def test_total_geq_stat(self):
        result = self.s.compute_systematic_budget(K, PK, SIGMA_STAT)
        assert np.all(result["sigma_P_total"] >= SIGMA_STAT - 1e-20)

    def test_total_geq_sys(self):
        result = self.s.compute_systematic_budget(K, PK, SIGMA_STAT)
        assert np.all(result["sigma_P_total"] >= result["sigma_P_sys"] - 1e-20)

    def test_fraction_in_unit_interval(self):
        result = self.s.compute_systematic_budget(K, PK, SIGMA_STAT)
        assert np.all(result["fraction_sys"] >= 0)
        assert np.all(result["fraction_sys"] <= 1 + 1e-10)

    def test_no_systematics_returns_stat_only(self):
        result = self.s.compute_systematic_budget(
            K, PK, SIGMA_STAT,
            include_photoz=False, include_bias=False, include_geometry=False
        )
        np.testing.assert_allclose(result["sigma_P_sys"], 0.0)
        np.testing.assert_allclose(result["sigma_P_total"], SIGMA_STAT)

    def test_individual_flags(self):
        r_all = self.s.compute_systematic_budget(K, PK, SIGMA_STAT)
        r_no_photoz = self.s.compute_systematic_budget(K, PK, SIGMA_STAT, include_photoz=False)
        # Turning off a source should reduce or maintain total systematic error
        assert np.all(r_no_photoz["sigma_P_sys"] <= r_all["sigma_P_sys"] + 1e-20)


class TestPropagateToAphi:
    def setup_method(self):
        self.s = SystematicErrorBudget(z_eff=0.8)
        self.dP_dAphi = PK * 0.5

    def test_returns_positive_float(self):
        sigma = self.s.propagate_to_Aphi(K, SIGMA_STAT, self.dP_dAphi)
        assert sigma > 0
        assert np.isfinite(sigma)

    def test_larger_sigma_sys_gives_larger_sigma_Aphi(self):
        s1 = self.s.propagate_to_Aphi(K, SIGMA_STAT * 1.0, self.dP_dAphi)
        s2 = self.s.propagate_to_Aphi(K, SIGMA_STAT * 10.0, self.dP_dAphi)
        assert s2 > s1


class TestComputeTotalAphiError:
    def setup_method(self):
        self.s = SystematicErrorBudget(z_eff=0.8)

    def test_quadrature_formula(self):
        stat, sys = 0.003, 0.004
        total = self.s.compute_total_Aphi_error(stat, sys)
        assert total == pytest.approx(np.sqrt(stat**2 + sys**2), rel=1e-10)

    def test_zero_sys_gives_stat(self):
        total = self.s.compute_total_Aphi_error(0.005, 0.0)
        assert total == pytest.approx(0.005)

    def test_total_geq_either_component(self):
        stat, sys = 0.003, 0.004
        total = self.s.compute_total_Aphi_error(stat, sys)
        assert total >= stat
        assert total >= sys
