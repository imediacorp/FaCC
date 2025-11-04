import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import os

# Constants
phi = (1 + np.sqrt(5)) / 2
ln_phi = np.log(phi)
phi_conj = phi - 1
ln_phi_conj = np.log(phi_conj)

# Load data
data = np.loadtxt('real_cmb_lowl.csv', skiprows=1, delimiter=',')
ell, cl_data, sigma_cl = data[:,0], data[:,1], data[:,2]

# Handle zero uncertainties
zero_sigma_count = np.sum(sigma_cl <= 0)
if zero_sigma_count > 0:
    estimated_sigma = 0.05 * np.abs(cl_data) + 10.0
    sigma_cl = np.where(sigma_cl <= 0, estimated_sigma, sigma_cl)

# Filter invalid
valid_mask = np.isfinite(sigma_cl) & np.isfinite(cl_data) & (ell > 0) & (sigma_cl > 0)
ell = ell[valid_mask]
cl_data = cl_data[valid_mask]
sigma_cl = sigma_cl[valid_mask]

# Polynomial baseline
poly_degree = min(6, len(ell) - 1)
baseline_coeffs = np.polyfit(np.log(ell), cl_data, poly_degree)
cl_lcdm = np.polyval(baseline_coeffs, np.log(ell))

# Residuals
residuals = cl_data - cl_lcdm

# Dual oscillation model
def osc_model_dual(ell, amp_phi, phase_phi, amp_conj, phase_conj):
    return amp_phi * np.cos(2 * np.pi * np.log(ell) / ln_phi + phase_phi) + \
           amp_conj * np.cos(2 * np.pi * np.log(ell) / ln_phi_conj + phase_conj)


# Fit
amp_guess = np.std(residuals)
popt, pcov = curve_fit(osc_model_dual, ell, residuals, sigma=sigma_cl,
                       p0=[amp_guess, 0, amp_guess/2, 0], absolute_sigma=True, maxfev=5000)
amp_phi_fit, phase_phi_fit, amp_conj_fit, phase_conj_fit = popt
print(f"Fitted amp_phi: {amp_phi_fit:.2e}, phase_phi: {phase_phi_fit:.2f}")
print(f"Fitted amp_conj: {amp_conj_fit:.2e}, phase_conj: {phase_conj_fit:.2f}")

# Plot
plt.figure(figsize=(10, 6))
plt.errorbar(ell, residuals, yerr=sigma_cl, fmt='.', label='Residuals', alpha=0.6)
plt.plot(ell, osc_model_dual(ell, *popt), 'r-', label='Dual Log-Periodic Fit', linewidth=2)
plt.xscale('log')
plt.xlabel('ℓ')
plt.ylabel('ΔC_ℓ')
plt.title('CMB Residuals - Dual Mode')
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig('cmb_osc_dual.png', dpi=150)
print("Plot saved as 'cmb_osc_dual.png'")
plt.show()
