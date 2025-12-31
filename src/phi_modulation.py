"""
φ-modulation analysis for DESI forecasts

Implements a two-parameter extension to ΛCDM testing for log-periodic
oscillations in the matter power spectrum at scales determined by the golden ratio φ.

This cosmological research work is independent and separate from any other projects.
It originated as a thought experiment and hypothesis exploring whether the Golden Ratio
might be fundamental to cosmic structure, given self-similarity patterns observed
from plants to galaxies.
"""

import numpy as np
import camb
from scipy.interpolate import interp1d
import warnings

warnings.filterwarnings('ignore', category=UserWarning)


class PhiModulationModel:
    """
    Implements φ-modulated power spectrum within ΛCDM framework
    
    This class provides methods to:
    - Generate ΛCDM baseline power spectra using CAMB
    - Apply log-periodic φ-modulation to power spectra
    - Compute BAO signatures with modulation
    - Forecast DESI sensitivity using Fisher matrix analysis
    """
    
    def __init__(self, params=None):
        """
        Initialize with cosmological parameters
        
        Parameters
        ----------
        params : dict, optional
            Cosmological parameters. Defaults to Planck 2018.
            Expected keys: 'H0', 'ombh2', 'omch2', 'As', 'ns', 'tau'
        """
        if params is None:
            self.params = {
                'H0': 67.36,
                'ombh2': 0.02237,
                'omch2': 0.1200,
                'As': 2.1e-9,
                'ns': 0.9649,
                'tau': 0.0544
            }
        else:
            self.params = params
        
        # Golden ratio
        self.phi = (1 + np.sqrt(5)) / 2
        self.lnphi = np.log(self.phi)
        
    def get_base_power_spectrum(self, k_min=1e-4, k_max=10, npoints=500, z=0.0):
        """
        Get ΛCDM power spectrum using CAMB
        
        Parameters
        ----------
        k_min : float
            Minimum wavenumber [h/Mpc]
        k_max : float
            Maximum wavenumber [h/Mpc]
        npoints : int
            Number of k points
        z : float
            Redshift (default 0.0)
            
        Returns
        -------
        k : array
            Wavenumbers [h/Mpc]
        z_arr : array
            Redshift array
        Pk : array
            Power spectrum P(k) [(Mpc/h)^3]
        """
        # Set up CAMB parameters
        pars = camb.CAMBparams()
        pars.set_cosmology(
            H0=self.params['H0'],
            ombh2=self.params['ombh2'],
            omch2=self.params['omch2'],
            tau=self.params['tau']
        )
        pars.InitPower.set_params(
            As=self.params['As'],
            ns=self.params['ns']
        )
        
        # Set redshift and k range
        pars.set_matter_power(redshifts=[z], kmax=k_max)
        
        # Compute results
        results = camb.get_results(pars)
        powers = results.get_matter_power_spectrum(
            minkh=k_min, maxkh=k_max, npoints=npoints
        )
        
        return powers.k, powers.z, powers.P_k
    
    def apply_phi_modulation(self, k, Pk, A_phi=0.01, phi_phase=0.0, k_pivot=0.05):
        """
        Apply φ-modulation to power spectrum
        
        The modulation is:
        P_mod(k) = P(k) × [1 + A_φ × cos(2π × log(k/k_pivot) / ln(φ) + φ_0)]
        
        Parameters
        ----------
        k : array
            Wavenumbers [h/Mpc]
        Pk : array
            Power spectrum values [(Mpc/h)^3]
        A_phi : float
            Amplitude of modulation (typically < 0.01)
        phi_phase : float
            Phase offset [radians]
        k_pivot : float
            Pivot scale for log-periodicity [h/Mpc]
            
        Returns
        -------
        Pk_mod : array
            Modulated power spectrum [(Mpc/h)^3]
        modulation : array
            Modulation factor (1 + A_φ * cos(...))
        """
        # Log-periodic modulation
        # Avoid log(0) by using k/k_pivot with small offset if needed
        k_ratio = np.where(k > 0, k / k_pivot, 1e-10)
        modulation = 1 + A_phi * np.cos(
            2 * np.pi * np.log(k_ratio) / self.lnphi + phi_phase
        )
        
        return Pk * modulation, modulation
    
    def compute_bao_signature(self, z=0.5, A_phi=0.01, r_min=80, r_max=120, n_r=200):
        """
        Compute BAO signature with φ-modulation
        
        Computes the correlation function ξ(r) via Fourier transform of P(k).
        
        Parameters
        ----------
        z : float
            Redshift for power spectrum evaluation
        A_phi : float
            Amplitude of φ-modulation
        r_min : float
            Minimum separation [Mpc/h]
        r_max : float
            Maximum separation [Mpc/h]
        n_r : int
            Number of r points
            
        Returns
        -------
        r : array
            Separation [Mpc/h]
        xi_base : array
            Base ΛCDM correlation function
        xi_mod : array
            Modulated correlation function
        """
        # Get power spectrum
        k_full, z_arr, Pk_full = self.get_base_power_spectrum(
            k_min=1e-4, k_max=2.0, npoints=1000, z=z
        )
        Pk_base = Pk_full[0]  # Extract z=0 slice
        
        # Interpolate P(k) to finer grid if needed
        # Ensure we have good coverage around BAO scale (k ~ 0.1 h/Mpc)
        k = k_full
        Pk = Pk_base
        
        # Define r range around BAO scale
        r = np.linspace(r_min, r_max, n_r)
        
        # Compute correlation function via Fourier transform
        # ξ(r) = ∫ P(k) sin(kr)/(kr) * k^2 dk / (2π^2)
        xi_base = np.zeros_like(r)
        xi_mod = np.zeros_like(r)
        
        # Apply modulation to get modulated P(k)
        Pk_mod, _ = self.apply_phi_modulation(k, Pk, A_phi=A_phi)
        
        for i, ri in enumerate(r):
            # Base correlation function
            integrand_base = Pk * np.sin(k * ri) / (k * ri + 1e-10) * k**2
            xi_base[i] = np.trapz(integrand_base, k)
            
            # Modulated correlation function
            integrand_mod = Pk_mod * np.sin(k * ri) / (k * ri + 1e-10) * k**2
            xi_mod[i] = np.trapz(integrand_mod, k)
        
        # Normalize
        xi_base /= (2 * np.pi**2)
        xi_mod /= (2 * np.pi**2)
        
        return r, xi_base, xi_mod
    
    def forecast_desi_sensitivity(self, A_phi_true=0.01, k_min=0.01, k_max=0.3, n_k=50):
        """
        Forecast DESI sensitivity using Fisher matrix approximation
        
        Parameters
        ----------
        A_phi_true : float
            True value of A_φ to forecast sensitivity for
        k_min : float
            Minimum k [h/Mpc] for DESI reliable range
        k_max : float
            Maximum k [h/Mpc] for DESI reliable range
        n_k : int
            Number of k bins
            
        Returns
        -------
        result : dict
            Dictionary containing:
            - 'k': k values [h/Mpc]
            - 'Pk_base': Base power spectrum
            - 'Pk_mod': Modulated power spectrum
            - 'sigma_P': Error on P(k) per bin
            - 'sigma_Aphi': Forecast error on A_φ
            - 'SNR': Signal-to-noise ratio (A_φ_true / σ_Aφ)
        """
        # DESI Year 5 specifications
        V_survey = 100  # (Gpc/h)^3 for DESI Y5
        z_eff = 0.8
        n_gal = 3e-4  # (h/Mpc)^3
        
        # Get power spectrum at effective k range
        k = np.logspace(np.log10(k_min), np.log10(k_max), n_k)
        
        # Base P(k) at z=z_eff
        _, _, Pk_full = self.get_base_power_spectrum(
            k_min=k_min*0.5, k_max=k_max*2, npoints=500, z=z_eff
        )
        Pk_base_array = Pk_full[0]
        
        # Interpolate to our k grid
        k_full = np.logspace(np.log10(k_min*0.5), np.log10(k_max*2), 500)
        Pk_interp = interp1d(k_full, Pk_base_array, kind='linear', 
                            bounds_error=False, fill_value='extrapolate')
        Pk_base = Pk_interp(k)
        
        # Modulated P(k)
        Pk_mod, mod_factor = self.apply_phi_modulation(
            k, Pk_base, A_phi=A_phi_true, k_pivot=0.05
        )
        
        # Cosmic variance error per k-bin
        # Δk = k * (log10(k_max) - log10(k_min)) / n_k
        Delta_k = k * (np.log10(k_max) - np.log10(k_min)) / n_k
        
        # Number of independent modes in each bin
        # N_modes = V_survey * k^2 * Δk / (2π^2)
        N_modes = V_survey * k**2 * Delta_k / (2 * np.pi**2)
        
        # Add shot noise term: P_shot = 1/n_gal
        P_shot = 1.0 / n_gal
        
        # Error on P(k) including cosmic variance and shot noise
        # σ_P = P * sqrt(2 / N_modes) for cosmic variance
        # Add shot noise contribution
        sigma_P_cv = Pk_base * np.sqrt(2 / (N_modes + 1e-10))
        sigma_P = np.sqrt(sigma_P_cv**2 + (P_shot / np.sqrt(N_modes + 1e-10))**2)
        
        # Signal: derivative with respect to A_phi
        # dP/dA_phi = P_base * cos(...) when A_phi is small
        # More precisely: dP/dA_phi = P_base * (mod_factor - 1) / A_phi_true
        dP_dA = Pk_base * (mod_factor - 1) / (A_phi_true + 1e-10)
        
        # Fisher matrix element for A_phi
        # F_Aphi = Σ_k (dP/dA_phi)^2 / σ_P^2
        F_Aphi = np.sum(dP_dA**2 / (sigma_P**2 + 1e-20))
        
        # Forecast uncertainty
        sigma_Aphi = 1.0 / np.sqrt(F_Aphi + 1e-20)
        
        # Signal-to-noise ratio
        SNR = A_phi_true / (sigma_Aphi + 1e-20)
        
        return {
            'k': k,
            'Pk_base': Pk_base,
            'Pk_mod': Pk_mod,
            'sigma_P': sigma_P,
            'sigma_Aphi': sigma_Aphi,
            'SNR': SNR,
            'mod_factor': mod_factor
        }
    
    def forecast_desi_sensitivity_with_systematics(self, A_phi_true=0.01, 
                                                   k_min=0.01, k_max=0.3, 
                                                   n_k=50, include_systematics=True):
        """
        Forecast DESI sensitivity including systematic error budget
        
        This method extends forecast_desi_sensitivity() by including systematic
        error contributions from photo-z errors, bias uncertainties, and survey geometry.
        
        Parameters
        ----------
        A_phi_true : float
            True value of A_φ to forecast sensitivity for
        k_min : float
            Minimum k [h/Mpc] for DESI reliable range
        k_max : float
            Maximum k [h/Mpc] for DESI reliable range
        n_k : int
            Number of k bins
        include_systematics : bool
            Whether to include systematic error contributions
            
        Returns
        -------
        result : dict
            Dictionary containing all fields from forecast_desi_sensitivity() plus:
            - 'sigma_Aphi_stat': Statistical uncertainty only
            - 'sigma_Aphi_sys': Systematic uncertainty contribution
            - 'sigma_Aphi_total': Total uncertainty (stat + sys)
            - 'systematic_budget': Detailed systematic error breakdown
        """
        # Get base forecast (statistical only)
        base_result = self.forecast_desi_sensitivity(
            A_phi_true=A_phi_true, k_min=k_min, k_max=k_max, n_k=n_k
        )
        
        if not include_systematics:
            # Return base result with consistent naming
            result = base_result.copy()
            result['sigma_Aphi_stat'] = result['sigma_Aphi']
            result['sigma_Aphi_sys'] = 0.0
            result['sigma_Aphi_total'] = result['sigma_Aphi']
            return result
        
        # Import systematics module (avoid circular import)
        try:
            from .systematics import SystematicErrorBudget
        except ImportError:
            # If systematics module not available, return base result
            warnings.warn("systematics module not available, returning statistical-only forecast")
            result = base_result.copy()
            result['sigma_Aphi_stat'] = result['sigma_Aphi']
            result['sigma_Aphi_sys'] = 0.0
            result['sigma_Aphi_total'] = result['sigma_Aphi']
            return result
        
        # Compute systematic error budget
        z_eff = 0.8  # DESI effective redshift
        sys_budget = SystematicErrorBudget(z_eff=z_eff)
        
        k = base_result['k']
        Pk_base = base_result['Pk_base']
        sigma_P_stat = base_result['sigma_P']
        dP_dAphi = Pk_base * (base_result['mod_factor'] - 1) / (A_phi_true + 1e-10)
        
        # Get systematic error breakdown
        systematic_result = sys_budget.compute_systematic_budget(
            k, Pk_base, sigma_P_stat
        )
        
        # Propagate systematic errors to A_φ
        sigma_Aphi_sys = sys_budget.propagate_to_Aphi(
            k, systematic_result['sigma_P_sys'], dP_dAphi
        )
        
        # Total uncertainty
        sigma_Aphi_stat = base_result['sigma_Aphi']
        sigma_Aphi_total = sys_budget.compute_total_Aphi_error(
            sigma_Aphi_stat, sigma_Aphi_sys
        )
        
        # Update result dictionary
        result = base_result.copy()
        result['sigma_Aphi_stat'] = sigma_Aphi_stat
        result['sigma_Aphi_sys'] = sigma_Aphi_sys
        result['sigma_Aphi_total'] = sigma_Aphi_total
        result['sigma_Aphi'] = sigma_Aphi_total  # Update main uncertainty field
        result['SNR'] = A_phi_true / (sigma_Aphi_total + 1e-20)  # Update SNR
        result['systematic_budget'] = systematic_result
        
        return result

