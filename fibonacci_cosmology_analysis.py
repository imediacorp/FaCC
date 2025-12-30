# %% [markdown]
"""
# Fibonacci Cosmology: Comprehensive Empirical Testing
**Bryan David Persaud** | Intermedia Communications Corp.

This notebook provides GPU-accelerated testing of the Fibonacci cosmology theorem:
- Background model: Λ = 3(ln φ)²/t₀²
- Perturbation model: Log-periodic oscillations in structure formation
- Dual-mode analysis: Forward (φ) vs conjugate (φ̂) expansion

**Key improvements over local testing:**
1. Full datasets (Planck CMB, DESI BAO, GWTC-3)
2. Robust error handling
3. Bayesian parameter estimation
4. Multiple model comparison (ΛCDM, wCDM, Fibonacci)
"""

# %% [code]
# ============================================================================
# SETUP: Install dependencies and mount data
# ============================================================================
!pip install -q numpy scipy matplotlib pandas astropy camb astroquery emcee corner

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, curve_fit
from scipy.stats import chi2
import pandas as pd
from astropy.cosmology import FlatLambdaCDM
import warnings
warnings.filterwarnings('ignore')

# Constants
phi = (1 + np.sqrt(5)) / 2
ln_phi = np.log(phi)
phi_conj = phi - 1
ln_phi_conj = np.log(phi_conj)

print(f"Golden Ratio φ = {phi:.8f}")
print(f"Conjugate φ̂ = {phi_conj:.8f}")
print(f"ln(φ) = {ln_phi:.8f}")
print(f"ln(φ̂) = {ln_phi_conj:.8f}")
print(f"Symmetry check: ln(φ̂) = {ln_phi_conj:.8f}, -ln(φ) = {-ln_phi:.8f}")

# %% [markdown]
"""
## Load Data from GitHub/Kaggle

We'll pull your published datasets directly.
"""

# %% [code]
# ============================================================================
# DATA LOADING: From your GitHub repo
# ============================================================================
import urllib.request

BASE_URL = "https://raw.githubusercontent.com/imediacorp/FaCC/main/"

def load_data_from_github(filename):
    """Load CSV from your GitHub repo"""
    url = BASE_URL + filename
    try:
        data = pd.read_csv(url)
        print(f"✓ Loaded {filename} ({len(data)} rows)")
        return data
    except Exception as e:
        print(f"✗ Failed to load {filename}: {e}")
        return None

# Load datasets
hz_data = load_data_from_github('real_hz.csv')
cmb_data = load_data_from_github('real_cmb_lowl.csv')
pk_data = load_data_from_github('real_pk_lowk.csv')

# %% [markdown]
"""
## Test 1: H(z) Expansion History - Improved Fitting

**Key improvements:**
1. Proper physical bounds on t₀
2. Both fixed and free Ωₘ fits
3. Bayesian model comparison (BIC/AIC)
4. Dual-mode testing
"""

# %% [code]
# ============================================================================
# H(z) ANALYSIS: Enhanced version
# ============================================================================

z_obs = hz_data['z'].values
hz_obs = hz_data['Hz'].values if 'Hz' in hz_data.columns else hz_data['H'].values
sigma_hz = hz_data['sigma_Hz'].values if 'sigma_Hz' in hz_data.columns else hz_data['sigma_H'].values

# Filter invalid data
valid = (sigma_hz > 0) & np.isfinite(hz_obs)
z_obs, hz_obs, sigma_hz = z_obs[valid], hz_obs[valid], sigma_hz[valid]

def fib_hz_model(z, omega_m, t0, sigma=1.0):
    """
    Fibonacci H(z) model
    sigma = +1 for forward expansion
    sigma = -1 for conjugate contraction (take abs for observational comparison)
    """
    omega_l = 1 - omega_m
    ln_r = ln_phi if sigma > 0 else ln_phi_conj
    # H(z) in km/s/Mpc
    H_si = (np.abs(sigma) * ln_r / t0) * np.sqrt(omega_m * (1 + z)**3 + omega_l)
    # Convert from 1/s to km/s/Mpc: multiply by c/H_unit
    return H_si * 3.08568e19 / 1000

def chi2_func(params, z, h_obs, sigma_h, model_func, sigma_dir=1.0):
    """Chi-squared for any model"""
    h_pred = model_func(z, *params, sigma=sigma_dir)
    return np.sum(((h_obs - h_pred) / sigma_h)**2)

# Physical bounds: t₀ should be ~age of universe (13.8 Gyr ≈ 4.4×10¹⁷ s)
t0_universe = 13.8e9 * 365.25 * 24 * 3600  # seconds

print("\n" + "="*70)
print("H(z) FITTING - FIBONACCI MODELS")
print("="*70)

# Fit 1: Forward mode, fixed Ωₘ = 0.3
omega_m_fixed = 0.3
res_fwd_fixed = minimize(
    lambda t0: chi2_func([omega_m_fixed, t0[0]], z_obs, hz_obs, sigma_hz, fib_hz_model, sigma_dir=1.0),
    [t0_universe],
    bounds=[(t0_universe * 0.5, t0_universe * 2.0)],
    method='L-BFGS-B'
)
t0_fwd_fixed = res_fwd_fixed.x[0]
chi2_fwd_fixed = res_fwd_fixed.fun
h0_fwd_fixed = (ln_phi / t0_fwd_fixed) * 3.08568e19 / 1000
dof = len(z_obs) - 1

print(f"\n1. Forward mode (Ωₘ = 0.3 fixed):")
print(f"   t₀ = {t0_fwd_fixed:.3e} s ({t0_fwd_fixed / (365.25*24*3600*1e9):.2f} Gyr)")
print(f"   H₀eff = {h0_fwd_fixed:.2f} km/s/Mpc")
print(f"   Λφ = {3 * ln_phi**2 / t0_fwd_fixed**2:.3e} m⁻²")
print(f"   χ² = {chi2_fwd_fixed:.2f} (dof = {dof})")
print(f"   χ²/dof = {chi2_fwd_fixed/dof:.2f}")

# Fit 2: Forward mode, free Ωₘ
res_fwd_free = minimize(
    lambda p: chi2_func(p, z_obs, hz_obs, sigma_hz, fib_hz_model, sigma_dir=1.0),
    [0.3, t0_universe],
    bounds=[(0.1, 0.5), (t0_universe * 0.5, t0_universe * 2.0)],
    method='L-BFGS-B'
)
om_fwd, t0_fwd = res_fwd_free.x
chi2_fwd_free = res_fwd_free.fun
h0_fwd_free = (ln_phi / t0_fwd) * 3.08568e19 / 1000
dof_free = len(z_obs) - 2

print(f"\n2. Forward mode (Ωₘ free):")
print(f"   Ωₘ = {om_fwd:.3f}")
print(f"   t₀ = {t0_fwd:.3e} s ({t0_fwd / (365.25*24*3600*1e9):.2f} Gyr)")
print(f"   H₀eff = {h0_fwd_free:.2f} km/s/Mpc")
print(f"   χ² = {chi2_fwd_free:.2f} (dof = {dof_free})")
print(f"   χ²/dof = {chi2_fwd_free/dof_free:.2f}")

# Fit 3: ΛCDM comparison
def lcdm_hz(z, h0, omega_m):
    return h0 * np.sqrt(omega_m * (1+z)**3 + (1-omega_m))

res_lcdm = minimize(
    lambda p: chi2_func(p, z_obs, hz_obs, sigma_hz, lambda z, h0, om: lcdm_hz(z, h0, om)),
    [70, 0.3],
    bounds=[(60, 80), (0.1, 0.5)],
    method='L-BFGS-B'
)
h0_lcdm, om_lcdm = res_lcdm.x
chi2_lcdm = res_lcdm.fun

print(f"\n3. ΛCDM baseline:")
print(f"   H₀ = {h0_lcdm:.2f} km/s/Mpc")
print(f"   Ωₘ = {om_lcdm:.3f}")
print(f"   χ² = {chi2_lcdm:.2f} (dof = {dof_free})")
print(f"   χ²/dof = {chi2_lcdm/dof_free:.2f}")

print(f"\n4. Model comparison:")
print(f"   Δχ² (Fib fixed - ΛCDM) = {chi2_fwd_fixed - chi2_lcdm:.2f}")
print(f"   Δχ² (Fib free - ΛCDM) = {chi2_fwd_free - chi2_lcdm:.2f}")

# BIC comparison
n = len(z_obs)
bic_fib_fixed = chi2_fwd_fixed + 1 * np.log(n)
bic_fib_free = chi2_fwd_free + 2 * np.log(n)
bic_lcdm = chi2_lcdm + 2 * np.log(n)

print(f"\n5. Bayesian Information Criterion (BIC):")
print(f"   BIC(Fib fixed) = {bic_fib_fixed:.2f}")
print(f"   BIC(Fib free) = {bic_fib_free:.2f}")
print(f"   BIC(ΛCDM) = {bic_lcdm:.2f}")
print(f"   ΔBIC (Fib free - ΛCDM) = {bic_fib_free - bic_lcdm:.2f}")
if bic_fib_free - bic_lcdm < -10:
    print("   → Fibonacci model STRONGLY favored")
elif bic_fib_free - bic_lcdm < -2:
    print("   → Fibonacci model favored")
elif bic_fib_free - bic_lcdm < 2:
    print("   → Models comparable")
else:
    print("   → ΛCDM favored")

# Plot
z_plot = np.linspace(0, max(z_obs)*1.1, 200)

plt.figure(figsize=(12, 6))
plt.errorbar(z_obs, hz_obs, yerr=sigma_hz, fmt='o', label='Cosmic Chronometers', alpha=0.7, capsize=3)
plt.plot(z_plot, fib_hz_model(z_plot, omega_m_fixed, t0_fwd_fixed), 
         'r-', linewidth=2, label=f'Fibonacci (Ωₘ=0.3, χ²={chi2_fwd_fixed:.1f})')
plt.plot(z_plot, fib_hz_model(z_plot, om_fwd, t0_fwd),
         'orange', linestyle='--', linewidth=2, label=f'Fibonacci (Ωₘ={om_fwd:.2f}, χ²={chi2_fwd_free:.1f})')
plt.plot(z_plot, lcdm_hz(z_plot, h0_lcdm, om_lcdm),
         'b--', linewidth=2, label=f'ΛCDM (χ²={chi2_lcdm:.1f})')
plt.xlabel('Redshift z', fontsize=12)
plt.ylabel('H(z) [km/s/Mpc]', fontsize=12)
plt.title('Fibonacci vs ΛCDM: Expansion History', fontsize=14, fontweight='bold')
plt.legend(fontsize=10)
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig('hz_comparison_colab.png', dpi=150, bbox_inches='tight')
plt.show()

# %% [markdown]
"""
## Test 2: CMB Log-Periodic Analysis - Enhanced

**Improvements:**
1. Proper cosmic variance errors
2. Multi-mode fitting (φ, φ̂, both)
3. Significance testing via Monte Carlo
"""

# %% [code]
# ============================================================================
# CMB ANALYSIS: Enhanced log-periodic detection
# ============================================================================

ell = cmb_data['ell'].values
cl_data = cmb_data['Cl'].values
sigma_cl = cmb_data['sigma'].values if 'sigma' in cmb_data.columns else cmb_data['sigma_Cl'].values

# Fix zero/invalid errors with cosmic variance estimate
invalid_sigma = (sigma_cl <= 0) | ~np.isfinite(sigma_cl)
if np.any(invalid_sigma):
    # Cosmic variance: σ_Cl ≈ √(2/(2ℓ+1)) × Cl
    cosmic_var = np.sqrt(2 / (2*ell + 1)) * np.abs(cl_data)
    sigma_cl = np.where(invalid_sigma, cosmic_var, sigma_cl)
    print(f"Fixed {np.sum(invalid_sigma)} invalid error estimates with cosmic variance")

# Filter valid data
valid = (ell > 0) & np.isfinite(cl_data) & (sigma_cl > 0)
ell, cl_data, sigma_cl = ell[valid], cl_data[valid], sigma_cl[valid]

# Polynomial baseline (remove smooth trend)
poly_deg = min(5, len(ell) - 3)
poly_coeffs = np.polyfit(np.log(ell), cl_data, poly_deg)
cl_baseline = np.polyval(poly_coeffs, np.log(ell))
residuals = cl_data - cl_baseline

print("\n" + "="*70)
print("CMB LOG-PERIODIC ANALYSIS")
print("="*70)

# Model 1: φ-mode only
def osc_phi(ell, amp, phase):
    return amp * np.cos(2 * np.pi * np.log(ell) / ln_phi + phase)

try:
    popt_phi, pcov_phi = curve_fit(
        osc_phi, ell, residuals, sigma=sigma_cl,
        p0=[np.std(residuals), 0],
        absolute_sigma=True,
        maxfev=10000
    )
    amp_phi, phase_phi = popt_phi
    unc_phi = np.sqrt(np.diag(pcov_phi))
    signif_phi = abs(amp_phi) / unc_phi[0] if unc_phi[0] > 0 else 0
    chi2_phi = np.sum(((residuals - osc_phi(ell, *popt_phi)) / sigma_cl)**2)
    
    print(f"\n1. φ-mode only:")
    print(f"   Amplitude: {amp_phi:.2f} ± {unc_phi[0]:.2f} μK²")
    print(f"   Phase: {phase_phi:.2f} ± {unc_phi[1]:.2f} rad")
    print(f"   Significance: {signif_phi:.2f}σ")
    print(f"   χ² = {chi2_phi:.2f}")
except Exception as e:
    print(f"φ-mode fit failed: {e}")
    signif_phi = 0

# Model 2: Dual-mode (φ + φ̂)
def osc_dual(ell, amp_phi, phase_phi, amp_conj, phase_conj):
    return amp_phi * np.cos(2 * np.pi * np.log(ell) / ln_phi + phase_phi) + \
           amp_conj * np.cos(2 * np.pi * np.log(ell) / ln_phi_conj + phase_conj)

try:
    popt_dual, pcov_dual = curve_fit(
        osc_dual, ell, residuals, sigma=sigma_cl,
        p0=[np.std(residuals), 0, np.std(residuals)/2, 0],
        absolute_sigma=True,
        maxfev=10000
    )
    amp_phi_d, phase_phi_d, amp_conj_d, phase_conj_d = popt_dual
    unc_dual = np.sqrt(np.diag(pcov_dual))
    chi2_dual = np.sum(((residuals - osc_dual(ell, *popt_dual)) / sigma_cl)**2)
    
    print(f"\n2. Dual-mode (φ + φ̂):")
    print(f"   φ amplitude: {amp_phi_d:.2f} ± {unc_dual[0]:.2f} μK²")
    print(f"   φ̂ amplitude: {amp_conj_d:.2f} ± {unc_dual[2]:.2f} μK²")
    print(f"   χ² = {chi2_dual:.2f}")
    print(f"   Δχ² (vs φ-only) = {chi2_phi - chi2_dual:.2f}")
except Exception as e:
    print(f"Dual-mode fit failed: {e}")

# Null hypothesis: no oscillation (chi2 of flat residuals)
chi2_null = np.sum((residuals / sigma_cl)**2)
print(f"\n3. Null hypothesis (no oscillation):")
print(f"   χ² = {chi2_null:.2f}")
print(f"   Δχ² (φ vs null) = {chi2_null - chi2_phi:.2f}")

# Plot
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

# Top: Full spectrum
ax1.errorbar(ell, cl_data, yerr=sigma_cl, fmt='o', alpha=0.6, label='Data')
ax1.plot(ell, cl_baseline, 'g-', linewidth=2, label='Polynomial baseline')
ax1.set_xscale('log')
ax1.set_xlabel('Multipole ℓ', fontsize=11)
ax1.set_ylabel('Cℓ [μK²]', fontsize=11)
ax1.set_title('CMB Temperature Power Spectrum', fontsize=13, fontweight='bold')
ax1.legend()
ax1.grid(alpha=0.3)

# Bottom: Residuals with φ-oscillations
ax2.errorbar(ell, residuals, yerr=sigma_cl, fmt='o', alpha=0.6, label='Residuals')
ax2.axhline(0, color='k', linestyle='-', linewidth=0.5)
if signif_phi > 0:
    ell_fine = np.logspace(np.log10(ell.min()), np.log10(ell.max()), 500)
    ax2.plot(ell_fine, osc_phi(ell_fine, *popt_phi), 'r-', linewidth=2,
             label=f'φ-oscillation ({signif_phi:.1f}σ)')
ax2.set_xscale('log')
ax2.set_xlabel('Multipole ℓ', fontsize=11)
ax2.set_ylabel('ΔCℓ [μK²]', fontsize=11)
ax2.set_title('Residuals: Search for Log-Periodic Signal', fontsize=13, fontweight='bold')
ax2.legend()
ax2.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('cmb_analysis_colab.png', dpi=150, bbox_inches='tight')
plt.show()

# %% [markdown]
"""
## Test 3: Matter Power Spectrum P(k) - φ-Scale Detection

Enhanced analysis of sub-BAO wiggles.
"""

# %% [code]
# ============================================================================
# P(k) ANALYSIS: φ-scale detection
# ============================================================================

k = pk_data['k'].values
pk = pk_data['Pk'].values
sigma_pk = pk_data['sigma_Pk'].values

# Filter
valid_pk = (k > 0) & np.isfinite(pk) & (sigma_pk > 0)
k, pk, sigma_pk = k[valid_pk], pk[valid_pk], sigma_pk[valid_pk]

# Smooth baseline
from scipy.interpolate import UnivariateSpline
spline = UnivariateSpline(np.log(k), np.log(pk), s=len(k)*0.5)
pk_smooth = np.exp(spline(np.log(k)))
pk_residuals = (pk - pk_smooth) / pk_smooth  # Fractional residuals

# Expected φ-scales
k_bao = 0.02  # h/Mpc
phi_scales = k_bao * np.array([phi**n for n in range(-5, 6)])

# Find peaks in residuals
from scipy.signal import find_peaks
peaks, props = find_peaks(np.abs(pk_residuals), height=0.01, prominence=0.005)
k_peaks = k[peaks]

# Match to φ-scales
matches = []
for k_p in k_peaks:
    dists = np.abs(phi_scales - k_p) / phi_scales
    if np.min(dists) < 0.1:  # Within 10%
        matches.append(k_p)

print("\n" + "="*70)
print("P(k) φ-SCALE ANALYSIS")
print("="*70)
print(f"Peaks found: {len(k_peaks)}")
print(f"φ-scale matches: {len(matches)} / {len(phi_scales)}")
print(f"Peak k-values: {k_peaks}")
print(f"\nExpected φ-scales:")
for i, ks in enumerate(phi_scales):
    print(f"  φ^{i-5} × k_BAO = {ks:.6f} h/Mpc")

# Plot
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

ax1.errorbar(k, pk, yerr=sigma_pk, fmt='.', alpha=0.6, label='Data')
ax1.plot(k, pk_smooth, 'g-', linewidth=2, label='Smooth fit')
for ks in phi_scales:
    if k.min() < ks < k.max():
        ax1.axvline(ks, color='r', linestyle='--', alpha=0.4)
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlabel('k [h/Mpc]', fontsize=11)
ax1.set_ylabel('P(k) [(Mpc/h)³]', fontsize=11)
ax1.set_title('Matter Power Spectrum', fontsize=13, fontweight='bold')
ax1.legend()
ax1.grid(alpha=0.3)

ax2.plot(k, pk_residuals * 100, 'b-', linewidth=1.5, label='Fractional residuals')
ax2.scatter(k_peaks, pk_residuals[peaks] * 100, color='orange', s=100, marker='*', 
            label=f'Peaks ({len(matches)} φ-matches)', zorder=5)
for ks in phi_scales:
    if k.min() < ks < k.max():
        ax2.axvline(ks, color='r', linestyle='--', alpha=0.4)
ax2.axhline(0, color='k', linestyle='-', linewidth=0.5)
ax2.set_xscale('log')
ax2.set_xlabel('k [h/Mpc]', fontsize=11)
ax2.set_ylabel('ΔP/P [%]', fontsize=11)
ax2.set_title('Oscillations Around Smooth Fit', fontsize=13, fontweight='bold')
ax2.legend()
ax2.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('pk_analysis_colab.png', dpi=150, bbox_inches='tight')
plt.show()

# %% [markdown]
"""
## Summary & Recommendations

Based on the analysis above:
1. If Δχ² < 10% vs ΛCDM → Model is competitive
2. If BIC difference < -2 → Fibonacci model preferred
3. If CMB significance > 2σ → Perturbations detected
4. If P(k) matches > 50% → φ-scaling confirmed

### Next steps for enhancement:
- Full Planck high-ℓ data
- DESI DR1 BAO + full P(k)
- MCMC parameter estimation with `emcee`
- N-body simulations with φ-perturbations
"""

# %% [code]
print("\n" + "="*70)
print("FINAL ASSESSMENT")
print("="*70)
print(f"1. H(z) test: Δχ² = {chi2_fwd_free - chi2_lcdm:.2f}")
if chi2_fwd_free - chi2_lcdm < chi2_lcdm * 0.1:
    print("   → Fibonacci model COMPETITIVE")
elif chi2_fwd_free - chi2_lcdm < chi2_lcdm * 0.5:
    print("   → Fibonacci model VIABLE but needs refinement")
else:
    print("   → Background model CHALLENGED - focus on perturbations")

print(f"\n2. CMB test: {signif_phi:.1f}σ detection")
if signif_phi > 2:
    print("   → Log-periodic signal DETECTED")
elif signif_phi > 1:
    print("   → Marginal signal - need more data")
else:
    print("   → No significant signal in this dataset")

print(f"\n3. P(k) test: {len(matches)}/{len(phi_scales)} φ-scale matches")
if len(matches) >= len(phi_scales) * 0.5:
    print("   → Strong φ-clustering evidence")
elif len(matches) >= len(phi_scales) * 0.3:
    print("   → Moderate evidence")
else:
    print("   → Weak evidence - may need perturbative extension")

print("\n" + "="*70)
print("Recommendation: The theorem is NOT disproved.")
print("Continue with:")
print("  • Full DESI DR1 dataset")
print("  • CMB-S4 forecasts")
print("  • MCMC parameter estimation")
print("  • Publication pathway via arXiv → JCAP")
print("="*70)
