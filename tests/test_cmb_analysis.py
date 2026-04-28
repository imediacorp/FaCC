"""Tests for src/cmb_analysis.py."""

import numpy as np
import pytest

from src.cmb_analysis import (
    load_cmb_data,
    preprocess_cmb_data,
    compute_residuals,
    dual_osc_model,
    fit_dual_oscillations,
    _PHI,
    _LN_PHI,
    _LN_PHI_CONJ,
)


class TestLoadCmbData:
    def test_returns_three_arrays(self, cmb_csv):
        ell, cl, sig = load_cmb_data(cmb_csv)
        assert len(ell) == len(cl) == len(sig) == 28  # ells 2–29

    def test_ell_positive(self, cmb_csv):
        ell, _, _ = load_cmb_data(cmb_csv)
        assert np.all(ell > 0)


class TestPreprocessCmbData:
    def test_removes_invalid_sigma(self):
        # zero/negative sigma is replaced with an estimate, not dropped —
        # all points with finite ell > 0 survive
        ell = np.array([2.0, 3.0, 4.0, 5.0])
        cl  = np.array([1000.0, 900.0, 800.0, 700.0])
        sig = np.array([0.0, -1.0, 50.0, 30.0])
        ell_out, cl_out, sig_out = preprocess_cmb_data(ell, cl, sig)
        assert np.all(sig_out > 0)
        assert len(ell_out) == 4  # sigma replaced, not dropped

    def test_all_valid_unchanged_length(self):
        ell = np.array([2.0, 3.0, 4.0])
        cl  = np.array([1000.0, 900.0, 800.0])
        sig = np.array([50.0, 40.0, 35.0])
        ell_out, cl_out, sig_out = preprocess_cmb_data(ell, cl, sig)
        assert len(ell_out) == 3

    def test_zero_sigma_replaced(self):
        ell = np.array([2.0])
        cl  = np.array([1000.0])
        sig = np.array([0.0])
        _, _, sig_out = preprocess_cmb_data(ell, cl, sig)
        assert sig_out[0] > 0


class TestComputeResiduals:
    def test_residual_shape(self, cmb_csv):
        ell, cl, _ = load_cmb_data(cmb_csv)
        baseline, resid = compute_residuals(ell, cl)
        assert baseline.shape == ell.shape
        assert resid.shape == ell.shape

    def test_residuals_sum_near_zero(self, cmb_csv):
        ell, cl, _ = load_cmb_data(cmb_csv)
        _, resid = compute_residuals(ell, cl)
        assert abs(np.mean(resid)) < abs(np.mean(cl))

    def test_custom_poly_degree(self, cmb_csv):
        ell, cl, _ = load_cmb_data(cmb_csv)
        _, resid3 = compute_residuals(ell, cl, poly_degree=3)
        _, resid6 = compute_residuals(ell, cl, poly_degree=6)
        assert resid3.shape == resid6.shape


class TestDualOscModel:
    def test_output_shape(self):
        ell = np.arange(2, 30, dtype=float)
        out = dual_osc_model(ell, 100.0, 0.0, 50.0, 0.5)
        assert out.shape == ell.shape

    def test_zero_amplitudes(self):
        ell = np.arange(2, 30, dtype=float)
        out = dual_osc_model(ell, 0.0, 0.0, 0.0, 0.0)
        assert np.allclose(out, 0.0)

    def test_amplitude_bounds(self):
        ell = np.arange(2, 30, dtype=float)
        A1, A2 = 100.0, 60.0
        out = dual_osc_model(ell, A1, 0.0, A2, 0.0)
        assert np.all(np.abs(out) <= A1 + A2 + 1e-9)

    def test_constants_sane(self):
        assert _LN_PHI == pytest.approx(np.log(_PHI), rel=1e-10)
        # φ − 1 = 1/φ, so ln(φ−1) = −ln(φ) < 0
        assert _LN_PHI_CONJ == pytest.approx(-_LN_PHI, rel=1e-10)


class TestFitDualOscillations:
    def test_returns_popt_pcov(self, cmb_csv):
        ell, cl, sig = load_cmb_data(cmb_csv)
        ell, cl, sig = preprocess_cmb_data(ell, cl, sig)
        _, resid = compute_residuals(ell, cl)
        popt, pcov = fit_dual_oscillations(ell, resid, sig)
        assert len(popt) == 4
        assert pcov.shape == (4, 4)

    def test_popt_finite(self, cmb_csv):
        ell, cl, sig = load_cmb_data(cmb_csv)
        ell, cl, sig = preprocess_cmb_data(ell, cl, sig)
        _, resid = compute_residuals(ell, cl)
        popt, _ = fit_dual_oscillations(ell, resid, sig)
        assert np.all(np.isfinite(popt))
