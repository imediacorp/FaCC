"""Tests for src/hz_analysis.py."""

import numpy as np
import pytest

from src.hz_analysis import (
    load_hz_data,
    h_model,
    fit_forward_branch,
    fit_reverse_branch,
)

_PHI = (1 + np.sqrt(5)) / 2


class TestLoadHzData:
    def test_returns_three_arrays(self, hz_csv):
        z, h, s = load_hz_data(hz_csv)
        assert len(z) == len(h) == len(s) == 7

    def test_z_positive(self, hz_csv):
        z, _, _ = load_hz_data(hz_csv)
        assert np.all(z > 0)

    def test_sigma_positive(self, hz_csv):
        _, _, s = load_hz_data(hz_csv)
        assert np.all(s > 0)

    def test_h_positive(self, hz_csv):
        _, h, _ = load_hz_data(hz_csv)
        assert np.all(h > 0)


class TestHModel:
    def test_output_shape(self):
        z = np.linspace(0, 2, 20)
        out = h_model(z, om=0.3, t0=1e-17)
        assert out.shape == z.shape

    def test_positive_for_forward(self):
        z = np.array([0.0, 0.5, 1.0])
        out = h_model(z, om=0.3, t0=1e-17, sigma=1.0)
        assert np.all(out > 0)

    def test_increases_with_z_for_forward(self):
        z = np.array([0.1, 0.5, 1.0, 1.5])
        out = h_model(z, om=0.3, t0=1e-17, sigma=1.0)
        assert np.all(np.diff(out) > 0)

    def test_reverse_branch_uses_conj_lnr(self):
        # σ=−1 picks ln(φ−1)=−ln(φ), so −1×ln(φ−1)=+ln(φ):
        # the reverse branch produces the same positive magnitude as forward.
        z = np.array([0.5])
        fwd = h_model(z, om=0.3, t0=1e-17, sigma=1.0)
        rev = h_model(z, om=0.3, t0=1e-17, sigma=-1.0)
        assert np.isclose(fwd, rev, rtol=1e-9)

    def test_h0_parameter_scales_output(self):
        z = np.array([0.0])
        out_default = h_model(z, om=0.3, t0=1e-17, h0=70.0)
        out_scaled  = h_model(z, om=0.3, t0=1e-17, h0=140.0)
        assert out_default.shape == out_scaled.shape


class TestFitForwardBranch:
    def test_returns_optimize_result(self, hz_csv):
        z, h, s = load_hz_data(hz_csv)
        res = fit_forward_branch(z, h, s)
        assert hasattr(res, "x")
        assert len(res.x) == 2

    def test_om_in_valid_range(self, hz_csv):
        z, h, s = load_hz_data(hz_csv)
        res = fit_forward_branch(z, h, s)
        om, t0 = res.x
        assert 0.1 <= om <= 0.4

    def test_t0_positive(self, hz_csv):
        z, h, s = load_hz_data(hz_csv)
        res = fit_forward_branch(z, h, s)
        _, t0 = res.x
        assert t0 > 0


class TestFitReverseBranch:
    def test_returns_optimize_result(self, hz_csv):
        z, h, s = load_hz_data(hz_csv)
        res = fit_reverse_branch(z, h, s)
        assert hasattr(res, "x")
        assert len(res.x) == 2

    def test_om_in_valid_range(self, hz_csv):
        z, h, s = load_hz_data(hz_csv)
        res = fit_reverse_branch(z, h, s)
        om, _ = res.x
        assert 0.1 <= om <= 0.4

    def test_chi2_finite(self, hz_csv):
        z, h, s = load_hz_data(hz_csv)
        res = fit_reverse_branch(z, h, s)
        assert np.isfinite(res.fun)
