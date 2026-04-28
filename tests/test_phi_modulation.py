"""Tests for src/phi_modulation.py — pure-math functions only (no CAMB required)."""

import numpy as np
import pytest
from unittest.mock import patch, MagicMock

from src.phi_modulation import PhiModulationModel

PHI = (1 + np.sqrt(5)) / 2


class TestPhiModulationModelInit:
    def test_default_params(self):
        model = PhiModulationModel()
        assert model.params["H0"] == pytest.approx(67.36)
        assert model.params["ombh2"] == pytest.approx(0.02237)
        assert model.params["omch2"] == pytest.approx(0.1200)

    def test_custom_params(self):
        custom = {"H0": 70.0, "ombh2": 0.022, "omch2": 0.12, "As": 2.1e-9, "ns": 0.96, "tau": 0.05}
        model = PhiModulationModel(params=custom)
        assert model.params["H0"] == pytest.approx(70.0)

    def test_golden_ratio(self):
        model = PhiModulationModel()
        assert model.phi == pytest.approx(PHI, rel=1e-10)
        assert model.lnphi == pytest.approx(np.log(PHI), rel=1e-10)

    def test_phi_bounds(self):
        model = PhiModulationModel()
        assert 1.6 < model.phi < 1.7


class TestApplyPhiModulation:
    def setup_method(self):
        self.model = PhiModulationModel()
        self.k = np.logspace(-2, 0, 100)
        self.Pk = np.ones(100) * 1e4

    def test_output_shape(self):
        Pk_mod, mod = self.model.apply_phi_modulation(self.k, self.Pk)
        assert Pk_mod.shape == self.k.shape
        assert mod.shape == self.k.shape

    def test_zero_amplitude_is_identity(self):
        Pk_mod, mod = self.model.apply_phi_modulation(self.k, self.Pk, A_phi=0.0)
        np.testing.assert_allclose(Pk_mod, self.Pk)
        np.testing.assert_allclose(mod, np.ones_like(self.k))

    def test_modulation_bounded(self):
        # |A_phi * cos(...)| <= A_phi, so mod in [1 - A_phi, 1 + A_phi]
        A_phi = 0.05
        _, mod = self.model.apply_phi_modulation(self.k, self.Pk, A_phi=A_phi)
        assert np.all(mod >= 1.0 - A_phi - 1e-12)
        assert np.all(mod <= 1.0 + A_phi + 1e-12)

    def test_modulation_always_positive(self):
        _, mod = self.model.apply_phi_modulation(self.k, self.Pk, A_phi=0.01)
        assert np.all(mod > 0)

    def test_Pk_mod_positive(self):
        Pk_mod, _ = self.model.apply_phi_modulation(self.k, self.Pk, A_phi=0.01)
        assert np.all(Pk_mod > 0)

    def test_phase_shifts_signal(self):
        _, mod0 = self.model.apply_phi_modulation(self.k, self.Pk, A_phi=0.01, phi_phase=0.0)
        _, mod_pi = self.model.apply_phi_modulation(self.k, self.Pk, A_phi=0.01, phi_phase=np.pi)
        # At phi_phase offset by pi, the cosine flips — peaks become troughs
        assert not np.allclose(mod0, mod_pi)

    def test_pivot_scale_effect(self):
        _, mod1 = self.model.apply_phi_modulation(self.k, self.Pk, A_phi=0.01, k_pivot=0.01)
        _, mod2 = self.model.apply_phi_modulation(self.k, self.Pk, A_phi=0.01, k_pivot=0.1)
        assert not np.allclose(mod1, mod2)

    def test_log_periodic_structure(self):
        # At k = k_pivot * phi^n the argument 2pi * log(phi^n) / ln(phi) = 2pi*n
        # So cos(...) = cos(phi_phase), giving mod = 1 + A_phi * cos(phi_phase)
        A_phi = 0.05
        k_pivot = 0.05
        phi_phase = 0.0
        k_at_phi = np.array([k_pivot * PHI**n for n in range(-3, 4)])
        _, mod = self.model.apply_phi_modulation(k_at_phi, np.ones_like(k_at_phi), A_phi=A_phi, k_pivot=k_pivot, phi_phase=phi_phase)
        expected = 1.0 + A_phi * np.cos(phi_phase)
        np.testing.assert_allclose(mod, expected, atol=1e-10)

    def test_handles_zero_k_gracefully(self):
        k_with_zero = np.array([0.0, 0.01, 0.1])
        Pk_small = np.ones(3) * 1e4
        Pk_mod, mod = self.model.apply_phi_modulation(k_with_zero, Pk_small)
        assert np.all(np.isfinite(Pk_mod))
        assert np.all(np.isfinite(mod))

    def test_array_Pk_preserved_shape(self):
        Pk_varying = self.k**-2 * 1e4
        Pk_mod, _ = self.model.apply_phi_modulation(self.k, Pk_varying, A_phi=0.02)
        assert Pk_mod.shape == self.k.shape


class TestGetBasePowerSpectrumMocked:
    """Tests for get_base_power_spectrum using a CAMB mock to avoid the Fortran dependency."""

    def setup_method(self):
        self.model = PhiModulationModel()

    def _make_camb_mock(self, k, z=0.0):
        mock_results = MagicMock()
        Pk = np.ones((1, len(k))) * 1e4
        mock_results.get_matter_power_spectrum.return_value = (k, np.array([z]), Pk)
        return mock_results

    def test_returns_three_arrays(self):
        k = np.logspace(-4, 1, 500)
        mock_results = self._make_camb_mock(k)
        with patch("camb.get_results", return_value=mock_results), \
             patch("camb.CAMBparams"):
            kh, z_arr, pk = self.model.get_base_power_spectrum()
        assert len(kh) == 500
        assert len(z_arr) == 1
        assert pk.shape[1] == 500

    def test_custom_k_range(self):
        k = np.logspace(-3, 0, 200)
        mock_results = self._make_camb_mock(k)
        with patch("camb.get_results", return_value=mock_results), \
             patch("camb.CAMBparams"):
            kh, _, _ = self.model.get_base_power_spectrum(k_min=1e-3, k_max=1.0, npoints=200)
        assert len(kh) == 200


class TestForecastDesiSensitivityMocked:
    """Smoke-tests forecast_desi_sensitivity with CAMB mocked out."""

    def setup_method(self):
        self.model = PhiModulationModel()

    def _make_camb_mock(self, npoints=500):
        k = np.logspace(-2, np.log10(0.6), npoints)
        Pk = (k / 0.05) ** (-1.5) * 5e3
        mock_results = MagicMock()
        mock_results.get_matter_power_spectrum.return_value = (k, np.array([0.8]), Pk[np.newaxis, :])
        return mock_results

    def test_output_keys(self):
        mock_results = self._make_camb_mock()
        with patch("camb.get_results", return_value=mock_results), \
             patch("camb.CAMBparams"):
            result = self.model.forecast_desi_sensitivity(A_phi_true=0.01)
        for key in ("k", "Pk_base", "Pk_mod", "sigma_P", "sigma_Aphi", "SNR", "mod_factor"):
            assert key in result

    def test_snr_positive(self):
        mock_results = self._make_camb_mock()
        with patch("camb.get_results", return_value=mock_results), \
             patch("camb.CAMBparams"):
            result = self.model.forecast_desi_sensitivity(A_phi_true=0.01)
        assert result["SNR"] > 0

    def test_larger_amplitude_gives_higher_snr(self):
        mock_results = self._make_camb_mock()
        with patch("camb.get_results", return_value=mock_results), \
             patch("camb.CAMBparams"):
            r1 = self.model.forecast_desi_sensitivity(A_phi_true=0.005)
        mock_results2 = self._make_camb_mock()
        with patch("camb.get_results", return_value=mock_results2), \
             patch("camb.CAMBparams"):
            r2 = self.model.forecast_desi_sensitivity(A_phi_true=0.02)
        assert r2["SNR"] > r1["SNR"]
