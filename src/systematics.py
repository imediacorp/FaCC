"""
Systematic error analysis for φ-modulation forecasts

This module provides tools to model and propagate systematic uncertainties
that affect DESI power spectrum measurements, including:
- Photometric redshift errors
- Galaxy bias uncertainties
- Survey geometry effects
- Overall systematic error budget

This cosmological research work is independent and separate from any other projects.
"""

import numpy as np
from scipy.interpolate import interp1d
import warnings

warnings.filterwarnings('ignore', category=UserWarning)


class SystematicErrorBudget:
    """
    Compute systematic error contributions to power spectrum forecasts
    
    Provides methods to estimate systematic uncertainties from various sources
    and combine them into an overall error budget for φ-modulation parameter constraints.
    """
    
    def __init__(self, z_eff=0.8):
        """
        Initialize systematic error budget calculator
        
        Parameters
        ----------
        z_eff : float
            Effective redshift of the analysis
        """
        self.z_eff = z_eff
        
        # Default DESI systematic error parameters
        # Based on DESI Year 5 expected performance
        self.sigma_z_photo = 0.02 * (1 + z_eff)  # Photo-z error
        self.sigma_bias = 0.05  # Relative uncertainty in galaxy bias
        self.sigma_fnl = 0.1  # Uncertainty in local PNG parameter (if relevant)
        
    def photo_z_error(self, k, Pk, sigma_z=None):
        """
        Estimate power spectrum error from photometric redshift uncertainties
        
        Photo-z errors smear the radial direction, affecting P(k) at all scales,
        with stronger impact at higher k (smaller scales).
        
        Parameters
        ----------
        k : array
            Wavenumbers [h/Mpc]
        Pk : array
            Power spectrum values [(Mpc/h)^3]
        sigma_z : float, optional
            Photo-z error (defaults to self.sigma_z_photo)
            
        Returns
        -------
        sigma_P_photoz : array
            Error on P(k) from photo-z uncertainties [(Mpc/h)^3]
        """
        if sigma_z is None:
            sigma_z = self.sigma_z_photo
        
        # Convert redshift error to distance error
        # r(z) = ∫_0^z c/H(z') dz'
        # Using approximate: dr/dz ≈ c/H(z) ≈ c/(H0 * sqrt(Ωm(1+z)^3))
        # Simplified: Δr/r ≈ Δz/(1+z) at low z
        H0 = 70.0  # km/s/Mpc (approximate)
        c = 299792.458  # km/s
        H_z = H0 * np.sqrt(0.3 * (1 + self.z_eff)**3 + 0.7)  # Simplified ΛCDM
        dr_dz = c / H_z  # Mpc/h
        
        # Photo-z error in distance units
        sigma_r = sigma_z * dr_dz
        
        # Error on P(k) scales roughly as: σ_P/P ≈ k * σ_r
        # More precisely: photo-z errors damp power at high k
        # Using approximation: σ_P/P ≈ (k * σ_r)^2 / 2 for small errors
        relative_error = (k * sigma_r)**2 / 2.0
        
        # Limit relative error to reasonable values
        relative_error = np.clip(relative_error, 0, 0.1)
        
        sigma_P_photoz = Pk * relative_error
        
        return sigma_P_photoz
    
    def bias_uncertainty(self, k, Pk, sigma_b=None):
        """
        Estimate power spectrum error from galaxy bias uncertainties
        
        Galaxy bias relates galaxy power spectrum to matter power spectrum:
        P_gal(k) = b^2 * P_matter(k)
        
        Uncertainties in bias propagate as: σ_P/P = 2 * σ_b/b
        
        Parameters
        ----------
        k : array
            Wavenumbers [h/Mpc]
        Pk : array
            Power spectrum values [(Mpc/h)^3]
        sigma_b : float, optional
            Relative uncertainty in bias (defaults to self.sigma_bias)
            
        Returns
        -------
        sigma_P_bias : array
            Error on P(k) from bias uncertainties [(Mpc/h)^3]
        """
        if sigma_b is None:
            sigma_b = self.sigma_bias
        
        # Bias uncertainty propagates as: σ_P/P = 2 * σ_b/b
        # Assuming b ≈ 1-2 for typical galaxies
        relative_error = 2.0 * sigma_b
        
        sigma_P_bias = Pk * relative_error
        
        return sigma_P_bias
    
    def survey_geometry_error(self, k, Pk, V_survey=100.0):
        """
        Estimate power spectrum error from survey geometry effects
        
        Survey window function, mask, and geometry introduce mode coupling
        and additional uncertainties beyond simple cosmic variance.
        
        Parameters
        ----------
        k : array
            Wavenumbers [h/Mpc]
        Pk : array
            Power spectrum values [(Mpc/h)^3]
        V_survey : float
            Survey volume [(Gpc/h)^3]
            
        Returns
        -------
        sigma_P_geometry : array
            Error on P(k) from survey geometry [(Mpc/h)^3]
        """
        # Survey geometry effects are typically subdominant to cosmic variance
        # but become important at low k (large scales)
        # Simple approximation: adds ~10-20% to cosmic variance error at low k
        
        # Scale-dependent: larger effect at low k
        k_min_survey = 2 * np.pi / (V_survey**(1/3))  # Approximate minimum k
        suppression_factor = 1.0 / (1.0 + (k / k_min_survey)**2)
        
        # Additional 10-20% error from geometry effects
        relative_error = 0.15 * suppression_factor
        
        sigma_P_geometry = Pk * relative_error
        
        return sigma_P_geometry
    
    def compute_systematic_budget(self, k, Pk, sigma_P_stat, 
                                  include_photoz=True, include_bias=True,
                                  include_geometry=True):
        """
        Compute total systematic error budget
        
        Combines various systematic error sources and adds them in quadrature
        with statistical errors.
        
        Parameters
        ----------
        k : array
            Wavenumbers [h/Mpc]
        Pk : array
            Power spectrum values [(Mpc/h)^3]
        sigma_P_stat : array
            Statistical error on P(k) [(Mpc/h)^3]
        include_photoz : bool
            Include photo-z systematic errors
        include_bias : bool
            Include bias systematic errors
        include_geometry : bool
            Include survey geometry systematic errors
            
        Returns
        -------
        result : dict
            Dictionary containing:
            - 'sigma_P_total': Total error (stat + sys) [(Mpc/h)^3]
            - 'sigma_P_sys': Systematic error component [(Mpc/h)^3]
            - 'sigma_P_photoz': Photo-z error [(Mpc/h)^3]
            - 'sigma_P_bias': Bias error [(Mpc/h)^3]
            - 'sigma_P_geometry': Geometry error [(Mpc/h)^3]
            - 'fraction_sys': Fraction of total error from systematics
        """
        # Initialize systematic error arrays
        sigma_P_sys_squared = np.zeros_like(Pk)
        
        sigma_P_photoz = np.zeros_like(Pk)
        sigma_P_bias = np.zeros_like(Pk)
        sigma_P_geometry = np.zeros_like(Pk)
        
        # Add photo-z errors
        if include_photoz:
            sigma_P_photoz = self.photo_z_error(k, Pk)
            sigma_P_sys_squared += sigma_P_photoz**2
        
        # Add bias errors
        if include_bias:
            sigma_P_bias = self.bias_uncertainty(k, Pk)
            sigma_P_sys_squared += sigma_P_bias**2
        
        # Add geometry errors
        if include_geometry:
            sigma_P_geometry = self.survey_geometry_error(k, Pk)
            sigma_P_sys_squared += sigma_P_geometry**2
        
        # Total systematic error (quadrature sum)
        sigma_P_sys = np.sqrt(sigma_P_sys_squared)
        
        # Total error (statistical + systematic in quadrature)
        sigma_P_total_squared = sigma_P_stat**2 + sigma_P_sys**2
        sigma_P_total = np.sqrt(sigma_P_total_squared)
        
        # Fraction from systematics
        fraction_sys = sigma_P_sys / (sigma_P_total + 1e-20)
        
        return {
            'sigma_P_total': sigma_P_total,
            'sigma_P_sys': sigma_P_sys,
            'sigma_P_photoz': sigma_P_photoz,
            'sigma_P_bias': sigma_P_bias,
            'sigma_P_geometry': sigma_P_geometry,
            'fraction_sys': fraction_sys
        }
    
    def propagate_to_Aphi(self, k, sigma_P_sys, dP_dAphi):
        """
        Propagate systematic errors to A_φ parameter constraint
        
        Converts systematic errors on P(k) to additional uncertainty on A_φ
        using Fisher matrix formalism.
        
        Parameters
        ----------
        k : array
            Wavenumbers [h/Mpc]
        sigma_P_sys : array
            Systematic error on P(k) [(Mpc/h)^3]
        dP_dAphi : array
            Derivative dP/dA_φ [(Mpc/h)^3]
            
        Returns
        -------
        sigma_Aphi_sys : float
            Additional systematic uncertainty on A_φ
        """
        # Fisher information from systematic errors
        # Treat systematic errors as additional uncertainty
        # F_sys = Σ_k (dP/dA_φ)^2 / (σ_sys^2)
        sigma_P_sys_sq = sigma_P_sys**2 + 1e-20  # Avoid division by zero
        
        F_sys = np.sum(dP_dAphi**2 / sigma_P_sys_sq)
        
        # Systematic error contribution to A_φ
        sigma_Aphi_sys = 1.0 / np.sqrt(F_sys + 1e-20)
        
        return sigma_Aphi_sys
    
    def compute_total_Aphi_error(self, sigma_Aphi_stat, sigma_Aphi_sys):
        """
        Combine statistical and systematic errors on A_φ
        
        Parameters
        ----------
        sigma_Aphi_stat : float
            Statistical uncertainty on A_φ
        sigma_Aphi_sys : float
            Systematic uncertainty on A_φ
            
        Returns
        -------
        sigma_Aphi_total : float
            Total uncertainty on A_φ
        """
        # Add in quadrature
        sigma_Aphi_total = np.sqrt(sigma_Aphi_stat**2 + sigma_Aphi_sys**2)
        
        return sigma_Aphi_total

