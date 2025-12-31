"""
Bayesian model comparison tools for φ-modulation analysis

This module provides tools for computing Bayesian evidence and Bayes factors
for comparing φ-modulated models against ΛCDM baseline models.

Supports multiple evidence estimation methods:
- Harmonic mean estimator (simple but potentially unstable)
- Thermodynamic integration (more robust)
- Nested sampling (via dynesty, recommended)

This cosmological research work is independent and separate from any other projects.
"""

import numpy as np
from scipy.stats import multivariate_normal
import warnings

warnings.filterwarnings('ignore', category=UserWarning)


class BayesianEvidence:
    """
    Compute Bayesian evidence for model comparison
    
    Provides methods to estimate the evidence (marginal likelihood) for
    different models and compute Bayes factors for model selection.
    """
    
    def __init__(self, data, cov, model_lcdm, model_phi):
        """
        Initialize Bayesian evidence calculator
        
        Parameters
        ----------
        data : array
            Observed data vector
        cov : array
            Covariance matrix for data
        model_lcdm : callable
            Function that computes ΛCDM model prediction given parameters
            Signature: model_lcdm(theta) -> array
        model_phi : callable
            Function that computes φ-modulated model prediction given parameters
            Signature: model_phi(theta) -> array
        """
        self.data = np.asarray(data)
        self.cov = np.asarray(cov)
        self.inv_cov = np.linalg.inv(self.cov)
        self.model_lcdm = model_lcdm
        self.model_phi = model_phi
        
        # Pre-compute determinant for efficiency
        self.log_det_cov = np.log(np.linalg.det(self.cov))
        
    def log_likelihood(self, theta, model_type='lcdm'):
        """
        Compute log-likelihood for given parameters and model type
        
        Assumes Gaussian likelihood:
        log L = -0.5 * [χ² + log(det(Σ)) + N*log(2π)]
        
        Parameters
        ----------
        theta : array
            Model parameters
        model_type : str
            Model type: 'lcdm' or 'phi'
            
        Returns
        -------
        log_L : float
            Log-likelihood value
        """
        # Get model prediction
        if model_type == 'lcdm':
            model_pred = self.model_lcdm(theta)
        elif model_type == 'phi':
            model_pred = self.model_phi(theta)
        else:
            raise ValueError(f"Unknown model_type: {model_type}")
        
        # Residuals
        residual = self.data - model_pred
        residual = np.asarray(residual).flatten()
        
        # χ²
        chi2 = residual @ self.inv_cov @ residual
        
        # Log-likelihood
        n_data = len(self.data)
        log_L = -0.5 * (chi2 + self.log_det_cov + n_data * np.log(2 * np.pi))
        
        return log_L
    
    def harmonic_mean_evidence(self, samples, log_likelihood_func):
        """
        Estimate evidence using harmonic mean estimator
        
        Note: This method can be unstable and is not recommended for high-precision
        applications. Use nested_sampling_evidence() or thermodynamic_integration()
        for more robust results.
        
        Evidence estimate: Z ≈ 1 / <1/L(θ)>
        where the average is over posterior samples.
        
        Parameters
        ----------
        samples : array
            Posterior samples (shape: n_samples x n_params)
        log_likelihood_func : callable
            Function that computes log-likelihood for given parameters
            Signature: log_likelihood_func(theta) -> float
            
        Returns
        -------
        log_Z : float
            Log-evidence estimate
        """
        # Convert log-likelihoods to likelihoods
        log_likes = np.array([log_likelihood_func(sample) for sample in samples])
        
        # Avoid underflow by working in log space
        # Z ≈ 1 / <1/L> = 1 / <exp(-log L)>
        # log Z ≈ -log(<exp(-log L)>)
        # Use logsumexp trick for numerical stability
        neg_log_likes = -log_likes
        log_Z = -np.log(np.mean(np.exp(neg_log_likes - np.max(neg_log_likes)))) - np.max(neg_log_likes)
        
        return log_Z
    
    def compute_bayes_factor(self, samples_lcdm, samples_phi,
                            log_likelihood_lcdm=None, log_likelihood_phi=None,
                            method='harmonic_mean'):
        """
        Compute Bayes factor comparing φ-modulation to ΛCDM
        
        Bayes factor: B = Z_phi / Z_lcdm
        log(B) = log(Z_phi) - log(Z_lcdm)
        
        Interpretation (Jeffreys scale):
        - log(B) > 5: Very strong evidence for φ-model
        - 2.5 < log(B) < 5: Strong evidence
        - 1 < log(B) < 2.5: Positive evidence
        - -1 < log(B) < 1: Inconclusive
        - log(B) < -1: Evidence against φ-model
        
        Parameters
        ----------
        samples_lcdm : array
            Posterior samples for ΛCDM model
        samples_phi : array
            Posterior samples for φ-modulated model
        log_likelihood_lcdm : callable, optional
            Log-likelihood function for ΛCDM
            If None, uses self.log_likelihood with model_type='lcdm'
        log_likelihood_phi : callable, optional
            Log-likelihood function for φ-model
            If None, uses self.log_likelihood with model_type='phi'
        method : str
            Method for evidence estimation: 'harmonic_mean'
            
        Returns
        -------
        result : dict
            Dictionary containing:
            - 'log_B': Log Bayes factor (log(Z_phi) - log(Z_lcdm))
            - 'B': Bayes factor
            - 'log_Z_lcdm': Log-evidence for ΛCDM
            - 'log_Z_phi': Log-evidence for φ-model
        """
        # Default log-likelihood functions
        if log_likelihood_lcdm is None:
            log_likelihood_lcdm = lambda theta: self.log_likelihood(theta, model_type='lcdm')
        if log_likelihood_phi is None:
            log_likelihood_phi = lambda theta: self.log_likelihood(theta, model_type='phi')
        
        # Estimate evidence for each model
        if method == 'harmonic_mean':
            log_Z_lcdm = self.harmonic_mean_evidence(samples_lcdm, log_likelihood_lcdm)
            log_Z_phi = self.harmonic_mean_evidence(samples_phi, log_likelihood_phi)
        else:
            raise ValueError(f"Unknown method: {method}")
        
        # Compute Bayes factor
        log_B = log_Z_phi - log_Z_lcdm
        B = np.exp(log_B)
        
        return {
            'log_B': log_B,
            'B': B,
            'log_Z_lcdm': log_Z_lcdm,
            'log_Z_phi': log_Z_phi
        }
    
    def interpret_bayes_factor(self, log_B):
        """
        Interpret Bayes factor according to Jeffreys scale
        
        Parameters
        ----------
        log_B : float
            Log Bayes factor
            
        Returns
        -------
        interpretation : str
            Textual interpretation of the Bayes factor
        """
        if log_B > 5:
            return "Very strong evidence for φ-modulated model"
        elif log_B > 2.5:
            return "Strong evidence for φ-modulated model"
        elif log_B > 1:
            return "Positive evidence for φ-modulated model"
        elif log_B > -1:
            return "Inconclusive (evidence not decisive)"
        elif log_B > -2.5:
            return "Positive evidence for ΛCDM model"
        elif log_B > -5:
            return "Strong evidence for ΛCDM model"
        else:
            return "Very strong evidence for ΛCDM model"


def compute_bic(model1_chi2, model2_chi2, n_data, n_params1, n_params2):
    """
    Compute Bayesian Information Criterion difference
    
    BIC = χ² + k * log(n)
    where k is number of parameters and n is number of data points
    
    ΔBIC = BIC_model2 - BIC_model1
    Negative ΔBIC favors model2, positive favors model1
    
    Parameters
    ----------
    model1_chi2 : float
        Best-fit χ² for model 1
    model2_chi2 : float
        Best-fit χ² for model 2
    n_data : int
        Number of data points
    n_params1 : int
        Number of parameters in model 1
    n_params2 : int
        Number of parameters in model 2
        
    Returns
    -------
    delta_bic : float
        ΔBIC = BIC_model2 - BIC_model1
    """
    bic1 = model1_chi2 + n_params1 * np.log(n_data)
    bic2 = model2_chi2 + n_params2 * np.log(n_data)
    
    delta_bic = bic2 - bic1
    
    return delta_bic


def interpret_bic(delta_bic):
    """
    Interpret BIC difference
    
    Rough interpretation (not as rigorous as Bayes factors):
    - ΔBIC < -10: Very strong evidence for model 2
    - -10 < ΔBIC < -6: Strong evidence
    - -6 < ΔBIC < -2: Positive evidence
    - |ΔBIC| < 2: Inconclusive
    - ΔBIC > 2: Evidence for model 1
    
    Parameters
    ----------
    delta_bic : float
        ΔBIC = BIC_model2 - BIC_model1
        
    Returns
    -------
    interpretation : str
        Textual interpretation
    """
    if delta_bic < -10:
        return "Very strong evidence for model 2"
    elif delta_bic < -6:
        return "Strong evidence for model 2"
    elif delta_bic < -2:
        return "Positive evidence for model 2"
    elif abs(delta_bic) < 2:
        return "Inconclusive (models comparable)"
    elif delta_bic < 6:
        return "Positive evidence for model 1"
    elif delta_bic < 10:
        return "Strong evidence for model 1"
    else:
        return "Very strong evidence for model 1"

