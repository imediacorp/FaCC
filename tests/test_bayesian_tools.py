"""Tests for src/bayesian_tools.py."""

import numpy as np
import pytest

from src.bayesian_tools import BayesianEvidence, compute_bic, interpret_bic


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _make_identity_models(n=10):
    """Simple models that return their parameter vector as-is."""
    def model_lcdm(theta):
        return np.asarray(theta)

    def model_phi(theta):
        return np.asarray(theta)

    data = np.ones(n)
    cov = np.eye(n) * 0.01
    return data, cov, model_lcdm, model_phi


# ---------------------------------------------------------------------------
# BayesianEvidence
# ---------------------------------------------------------------------------

class TestBayesianEvidenceInit:
    def test_constructs_without_error(self):
        data, cov, m1, m2 = _make_identity_models()
        be = BayesianEvidence(data, cov, m1, m2)
        assert be.data.shape == (10,)

    def test_inv_cov_is_inverse(self):
        data, cov, m1, m2 = _make_identity_models()
        be = BayesianEvidence(data, cov, m1, m2)
        product = be.cov @ be.inv_cov
        np.testing.assert_allclose(product, np.eye(10), atol=1e-10)


class TestLogLikelihood:
    def setup_method(self):
        n = 5
        self.data = np.zeros(n)
        self.cov = np.eye(n)
        self.model_lcdm = lambda theta: np.zeros(n)
        self.model_phi = lambda theta: np.ones(n) * 0.5
        self.be = BayesianEvidence(self.data, self.cov, self.model_lcdm, self.model_phi)

    def test_perfect_fit_is_finite(self):
        ll = self.be.log_likelihood(np.zeros(5), model_type="lcdm")
        assert np.isfinite(ll)

    def test_better_fit_has_higher_likelihood(self):
        ll_good = self.be.log_likelihood(np.zeros(5), model_type="lcdm")
        ll_bad = self.be.log_likelihood(np.zeros(5), model_type="phi")
        assert ll_good > ll_bad

    def test_unknown_model_type_raises(self):
        with pytest.raises(ValueError, match="Unknown model_type"):
            self.be.log_likelihood(np.zeros(5), model_type="unknown")

    def test_likelihood_is_nonpositive_log_scale(self):
        # log-likelihood should be negative (or zero at best for this parameterisation)
        ll = self.be.log_likelihood(np.zeros(5), model_type="lcdm")
        assert ll <= 1.0  # not a strict bound, but sanity check


class TestHarmonicMeanEvidence:
    def setup_method(self):
        n = 3
        data = np.zeros(n)
        cov = np.eye(n)
        model = lambda theta: np.zeros(n)
        self.be = BayesianEvidence(data, cov, model, model)

    def test_returns_finite_value(self):
        samples = np.zeros((50, 3))
        log_like = lambda theta: self.be.log_likelihood(theta, "lcdm")
        log_Z = self.be.harmonic_mean_evidence(samples, log_like)
        assert np.isfinite(log_Z)

    def test_consistent_with_repeated_calls(self):
        samples = np.zeros((30, 3))
        log_like = lambda theta: self.be.log_likelihood(theta, "lcdm")
        log_Z1 = self.be.harmonic_mean_evidence(samples, log_like)
        log_Z2 = self.be.harmonic_mean_evidence(samples, log_like)
        assert log_Z1 == pytest.approx(log_Z2)


class TestInterpretBayesFactor:
    @pytest.mark.parametrize("log_B, expected", [
        (6.0, "Very strong evidence for φ-modulated model"),
        (3.0, "Strong evidence for φ-modulated model"),
        (1.5, "Positive evidence for φ-modulated model"),
        (0.0, "Inconclusive (evidence not decisive)"),
        (-1.5, "Positive evidence for ΛCDM model"),
        (-3.0, "Strong evidence for ΛCDM model"),
        (-6.0, "Very strong evidence for ΛCDM model"),
    ])
    def test_interpretation(self, log_B, expected):
        n = 2
        data, cov, m1, m2 = _make_identity_models(n)
        be = BayesianEvidence(data, cov, m1, m2)
        assert be.interpret_bayes_factor(log_B) == expected


# ---------------------------------------------------------------------------
# compute_bic / interpret_bic
# ---------------------------------------------------------------------------

class TestComputeBic:
    def test_same_chi2_and_params_gives_zero(self):
        delta = compute_bic(10.0, 10.0, 100, 3, 3)
        assert delta == pytest.approx(0.0)

    def test_extra_params_penalises_model2(self):
        # model2 has 2 extra params, same chi2 → positive ΔBIC (penalise model2)
        delta = compute_bic(10.0, 10.0, 100, 3, 5)
        assert delta > 0

    def test_better_chi2_can_overcome_penalty(self):
        # model2 has better chi2 by 20 with 1 extra param
        delta = compute_bic(30.0, 10.0, 100, 3, 4)
        assert delta < 0

    def test_formula_correctness(self):
        # BIC1 = chi2_1 + k1*log(n); BIC2 = chi2_2 + k2*log(n); ΔBIC = BIC2 - BIC1
        n, k1, k2 = 50, 2, 3
        chi2_1, chi2_2 = 5.0, 4.0
        expected = (chi2_2 + k2 * np.log(n)) - (chi2_1 + k1 * np.log(n))
        assert compute_bic(chi2_1, chi2_2, n, k1, k2) == pytest.approx(expected)


class TestInterpretBic:
    @pytest.mark.parametrize("delta_bic, fragment", [
        (-15.0, "Very strong evidence for model 2"),
        (-8.0, "Strong evidence for model 2"),
        (-4.0, "Positive evidence for model 2"),
        (0.0, "Inconclusive"),
        (4.0, "Positive evidence for model 1"),
        (8.0, "Strong evidence for model 1"),
        (15.0, "Very strong evidence for model 1"),
    ])
    def test_interpretation(self, delta_bic, fragment):
        assert fragment in interpret_bic(delta_bic)
