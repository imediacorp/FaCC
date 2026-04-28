"""Tests for src/desi_bao.py — all network calls are mocked."""

import numpy as np
import pytest
from unittest.mock import patch, MagicMock

from src.desi_bao import (
    parse_bao_mean,
    parse_bao_cov,
    sound_horizon_rd,
    lcdm_distance_ratios,
    phi_bao_residual_model,
    fit_phi_oscillation,
    phi_scale_redshifts,
    fetch_desi_bao,
)

_PHI = (1 + np.sqrt(5)) / 2
_LN_PHI = np.log(_PHI)

_SAMPLE_MEAN = """# [z] [value at z] [quantity]
0.295 7.94 DV_over_rs
0.510 13.59 DM_over_rs
0.510 21.86 DH_over_rs
0.706 17.35 DM_over_rs
0.706 19.46 DH_over_rs
"""

_SAMPLE_COV = """1.0 0.0 0.0 0.0 0.0
0.0 1.0 0.0 0.0 0.0
0.0 0.0 1.0 0.0 0.0
0.0 0.0 0.0 1.0 0.0
0.0 0.0 0.0 0.0 1.0
"""


class TestParseBaoMean:
    def test_correct_lengths(self):
        z, val, qty = parse_bao_mean(_SAMPLE_MEAN)
        assert len(z) == len(val) == len(qty) == 5

    def test_z_values(self):
        z, _, _ = parse_bao_mean(_SAMPLE_MEAN)
        assert z[0] == pytest.approx(0.295)
        assert z[1] == pytest.approx(0.510)

    def test_quantity_strings(self):
        _, _, qty = parse_bao_mean(_SAMPLE_MEAN)
        assert qty[0] == "DV_over_rs"
        assert qty[1] == "DM_over_rs"
        assert qty[2] == "DH_over_rs"

    def test_ignores_comment_lines(self):
        text = "# comment\n0.3 8.0 DV_over_rs\n"
        z, val, qty = parse_bao_mean(text)
        assert len(z) == 1


class TestParseBaoCov:
    def test_shape(self):
        cov = parse_bao_cov(_SAMPLE_COV)
        assert cov.shape == (5, 5)

    def test_diagonal_values(self):
        cov = parse_bao_cov(_SAMPLE_COV)
        assert np.allclose(np.diag(cov), 1.0)


class TestSoundHorizonRd:
    def test_planck_value(self):
        # Planck 2018 gives r_d ≈ 147 Mpc
        rd = sound_horizon_rd(H0=67.36, Om0=0.315, Ob0=0.0493)
        assert 145.0 < rd < 150.0

    def test_increases_with_lower_Om0(self):
        rd_low  = sound_horizon_rd(H0=67.36, Om0=0.28)
        rd_high = sound_horizon_rd(H0=67.36, Om0=0.35)
        assert rd_low > rd_high

    def test_returns_float(self):
        assert isinstance(sound_horizon_rd(), float)


class TestLcdmDistanceRatios:
    def test_output_length(self):
        z = np.array([0.3, 0.5, 1.0])
        qty = ["DV_over_rs", "DM_over_rs", "DH_over_rs"]
        pred = lcdm_distance_ratios(z, qty)
        assert len(pred) == 3

    def test_all_positive(self):
        z = np.array([0.3, 0.5, 1.0, 2.0])
        qty = ["DM_over_rs", "DH_over_rs", "DM_over_rs", "DV_over_rs"]
        pred = lcdm_distance_ratios(z, qty)
        assert np.all(pred > 0)

    def test_dm_increases_with_z(self):
        z = np.array([0.3, 0.6, 1.0])
        qty = ["DM_over_rs"] * 3
        pred = lcdm_distance_ratios(z, qty)
        assert np.all(np.diff(pred) > 0)

    def test_nan_for_unknown_qty(self):
        z = np.array([0.5])
        pred = lcdm_distance_ratios(z, ["UNKNOWN_QUANTITY"])
        assert np.isnan(pred[0])


class TestPhiBaoResidualModel:
    def test_zero_amplitude(self):
        z = np.linspace(0.1, 2.0, 30)
        out = phi_bao_residual_model(z, A_phi=0.0, delta=0.0)
        assert np.allclose(out, 0.0)

    def test_amplitude_bounded(self):
        z = np.linspace(0.1, 2.0, 100)
        A = 0.05
        out = phi_bao_residual_model(z, A_phi=A, delta=0.0)
        assert np.all(np.abs(out) <= A + 1e-12)

    def test_output_shape(self):
        z = np.linspace(0.2, 1.5, 50)
        out = phi_bao_residual_model(z, A_phi=0.01, delta=1.0)
        assert out.shape == z.shape

    def test_phase_shift(self):
        z = np.array([0.618])  # z = φ-1, so ln(1+z)/ln(φ) = 1
        out0 = phi_bao_residual_model(z, A_phi=1.0, delta=0.0)
        out_pi = phi_bao_residual_model(z, A_phi=1.0, delta=np.pi)
        assert np.isclose(out0, -out_pi, rtol=1e-9)


class TestFitPhiOscillation:
    def test_returns_popt_pcov(self):
        z = np.linspace(0.1, 2.0, 13)
        true_A, true_d = 0.02, 0.5
        resid = phi_bao_residual_model(z, true_A, true_d)
        sigma = np.full_like(z, 0.005)
        popt, pcov = fit_phi_oscillation(z, resid, sigma)
        assert len(popt) == 2
        assert pcov.shape == (2, 2)

    def test_recovers_known_amplitude(self):
        z = np.linspace(0.1, 2.0, 50)
        true_A = 0.03
        resid = phi_bao_residual_model(z, true_A, 0.0)
        sigma = np.full_like(z, 1e-4)
        popt, _ = fit_phi_oscillation(z, resid, sigma)
        assert abs(popt[0]) == pytest.approx(true_A, rel=0.05)

    def test_graceful_on_flat_residuals(self):
        z = np.linspace(0.1, 2.0, 13)
        resid = np.zeros_like(z)
        sigma = np.ones_like(z) * 0.01
        popt, pcov = fit_phi_oscillation(z, resid, sigma)
        assert len(popt) == 2  # should not raise


class TestPhiScaleRedshifts:
    def test_default_length(self):
        zs = phi_scale_redshifts()
        assert len(zs) == 8  # range(-3, 5)

    def test_n_zero_gives_zero(self):
        zs = phi_scale_redshifts(n_range=(0, 1))
        assert zs[0] == pytest.approx(0.0, abs=1e-10)

    def test_n_one_gives_phi_minus_one(self):
        zs = phi_scale_redshifts(n_range=(1, 2))
        assert zs[0] == pytest.approx(_PHI - 1, rel=1e-10)


class TestFetchDesiBao:
    @patch("src.desi_bao._fetch_text")
    def test_returns_dict_with_expected_keys(self, mock_fetch):
        mock_fetch.side_effect = [_SAMPLE_MEAN, _SAMPLE_COV]
        result = fetch_desi_bao("dr2")
        assert set(result.keys()) >= {"z", "val", "qty", "cov", "sigma"}

    @patch("src.desi_bao._fetch_text")
    def test_sigma_is_sqrt_of_cov_diag(self, mock_fetch):
        mock_fetch.side_effect = [_SAMPLE_MEAN, _SAMPLE_COV]
        result = fetch_desi_bao("dr2")
        assert np.allclose(result["sigma"], np.sqrt(np.diag(result["cov"])))

    @patch("src.desi_bao._fetch_text")
    def test_correct_number_of_measurements(self, mock_fetch):
        mock_fetch.side_effect = [_SAMPLE_MEAN, _SAMPLE_COV]
        result = fetch_desi_bao("dr2")
        assert len(result["z"]) == 5
