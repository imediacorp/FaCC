"""Tests for src/lss_analysis.py."""

import numpy as np
import pytest

from src.lss_analysis import (
    load_pk_data,
    compute_phi_scales,
    compute_spline_residuals,
    find_direct_peaks,
    find_residual_peaks,
    count_phi_matches,
)

_PHI = (1 + np.sqrt(5)) / 2


class TestLoadPkData:
    def test_returns_three_arrays(self, pk_csv):
        k, pk, s = load_pk_data(pk_csv)
        assert len(k) == len(pk) == len(s) == 40

    def test_k_positive(self, pk_csv):
        k, _, _ = load_pk_data(pk_csv)
        assert np.all(k > 0)

    def test_pk_positive(self, pk_csv):
        _, pk, _ = load_pk_data(pk_csv)
        assert np.all(pk > 0)


class TestComputePhiScales:
    def test_default_length(self):
        scales = compute_phi_scales()
        assert len(scales) == 11  # range(-5, 6)

    def test_ratio_between_adjacent_is_phi(self):
        scales = compute_phi_scales()
        ratios = scales[1:] / scales[:-1]
        assert np.allclose(ratios, _PHI, rtol=1e-10)

    def test_custom_k_bao(self):
        scales = compute_phi_scales(k_bao=0.05)
        # The n=0 term should equal k_bao
        assert scales[5] == pytest.approx(0.05, rel=1e-10)  # n=0 is index 5 in range(-5,6)

    def test_custom_n_range(self):
        scales = compute_phi_scales(n_range=(0, 3))
        assert len(scales) == 3

    def test_all_positive(self):
        scales = compute_phi_scales()
        assert np.all(scales > 0)


class TestComputeSplineResiduals:
    def test_shapes_match_input(self, pk_csv):
        k, pk, _ = load_pk_data(pk_csv)
        smooth, resid = compute_spline_residuals(k, pk)
        assert smooth.shape == k.shape
        assert resid.shape == k.shape

    def test_smooth_positive(self, pk_csv):
        k, pk, _ = load_pk_data(pk_csv)
        smooth, _ = compute_spline_residuals(k, pk)
        assert np.all(smooth > 0)

    def test_residuals_zero_mean_approx(self, pk_csv):
        k, pk, _ = load_pk_data(pk_csv)
        _, resid = compute_spline_residuals(k, pk)
        # For a smooth power law, residuals should be small relative to pk
        assert np.std(resid) < np.std(pk)


class TestFindDirectPeaks:
    def test_returns_array(self, pk_csv):
        k, pk, s = load_pk_data(pk_csv)
        peaks = find_direct_peaks(k, pk, s)
        assert isinstance(peaks, np.ndarray)

    def test_peaks_within_k_range(self, pk_csv):
        k, pk, s = load_pk_data(pk_csv)
        peaks = find_direct_peaks(k, pk, s)
        if len(peaks) > 0:
            assert np.all(peaks >= k.min())
            assert np.all(peaks <= k.max())

    def test_high_threshold_returns_fewer_peaks(self, pk_csv):
        k, pk, s = load_pk_data(pk_csv)
        peaks_low  = find_direct_peaks(k, pk, s, threshold_sigma=1.0)
        peaks_high = find_direct_peaks(k, pk, s, threshold_sigma=10.0)
        assert len(peaks_low) >= len(peaks_high)


class TestFindResidualPeaks:
    def test_returns_array(self, pk_csv):
        k, pk, s = load_pk_data(pk_csv)
        _, resid = compute_spline_residuals(k, pk)
        peaks = find_residual_peaks(k, resid, s)
        assert isinstance(peaks, np.ndarray)

    def test_peaks_in_k_range(self, pk_csv):
        k, pk, s = load_pk_data(pk_csv)
        _, resid = compute_spline_residuals(k, pk)
        peaks = find_residual_peaks(k, resid, s)
        if len(peaks) > 0:
            assert np.all(peaks >= k.min())
            assert np.all(peaks <= k.max())


class TestCountPhiMatches:
    def test_zero_peaks(self):
        phi_scales = compute_phi_scales()
        assert count_phi_matches(np.array([]), phi_scales) == 0

    def test_exact_match(self):
        phi_scales = compute_phi_scales()
        # A peak exactly at one φ-scale should count as 1
        k_peaks = np.array([phi_scales[3]])
        assert count_phi_matches(k_peaks, phi_scales, tolerance=1e-9) == 1

    def test_no_match(self):
        phi_scales = compute_phi_scales()
        k_peaks = np.array([99.9])  # far outside any φ-scale
        assert count_phi_matches(k_peaks, phi_scales, tolerance=0.001) == 0

    def test_multiple_matches(self):
        phi_scales = compute_phi_scales()
        k_peaks = phi_scales[:3]  # exactly the first three scales
        assert count_phi_matches(k_peaks, phi_scales, tolerance=1e-9) == 3
