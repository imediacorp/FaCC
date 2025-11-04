# hz_fitter.py
import numpy as np
from scipy.optimize import minimize, minimize_scalar
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM

phi = (1 + np.sqrt(5)) / 2
ln_phi = np.log(phi)

# Load data
data = np.loadtxt('real_hz.csv', delimiter=',', skiprows=1)
z_data, h_data, sigma_h = data[:, 0], data[:, 1], data[:, 2]

# VECTORIZED h_model
def h_model(z, om, t0):
    z = np.atleast_1d(z)
    ol = 1 - om
    hubble = (ln_phi / t0) * np.sqrt(om * (1 + z)**3 + ol)
    return (hubble * 3.08568e19).flatten()

# Fix Ω_m = 0.3, fit t0 only
om_fixed = 0.3

def chi2_t0(t0):
    h_pred = h_model(z_data, om_fixed, t0)
    return np.sum(((h_data - h_pred) / sigma_h)**2)

# TIGHTER BOUNDS
res = minimize_scalar(chi2_t0, bounds=(1.2e18, 1.6e18), method='bounded')
t0_fit = res.x
chi2_fit = res.fun
h0_eff = (ln_phi / t0_fit) * 3.08568e19
lambda_phi = 3 * (ln_phi)**2 / t0_fit**2

# ΛCDM
def chi2_lcdm(p):
    om, h0 = p
    cosmo = FlatLambdaCDM(H0=h0, Om0=om)
    return np.sum(((h_data - cosmo.H(z_data).value) / sigma_h)**2)

res_lcdm = minimize(chi2_lcdm, [0.3, 70], bounds=[(0.1, 0.5), (60, 80)])
om_lcdm, h0_lcdm = res_lcdm.x
chi2_lcdm = res_lcdm.fun

# PRINT
print("="*60)
print("H(z) FIT RESULTS (FINAL)")
print("="*60)
print(f"Fibonacci: Ω_m={om_fixed:.4f}, t₀={t0_fit:.2e} s")
print(f"          H₀^eff={h0_eff:.2f} km/s/Mpc, Λ_φ={lambda_phi:.2e} m⁻²")
print(f"          χ²={chi2_fit:.2f} (dof={len(z_data)-1})")
print(f"ΛCDM:     Ω_m={om_lcdm:.4f}, H₀={h0_lcdm:.2f}, χ²={chi2_lcdm:.2f}")
print(f"Δχ² = {chi2_fit - chi2_lcdm:.2f}")

# PLOT
cosmo_lcdm = FlatLambdaCDM(H0=h0_lcdm, Om0=om_lcdm)
plt.figure(figsize=(9,6))
plt.errorbar(z_data, h_data, yerr=sigma_h, fmt='o', label='Data', alpha=0.7)
plt.plot(z_data, h_model(z_data, om_fixed, t0_fit), 'r-', lw=2, label=f'Fibonacci (H₀^eff={h0_eff:.1f})')
plt.plot(z_data, cosmo_lcdm.H(z_data).value, 'k--', lw=2, label=f'ΛCDM (H₀={h0_lcdm:.1f})')
plt.xlabel('Redshift z'); plt.ylabel('H(z) [km/s/Mpc]')
plt.legend(); plt.grid(alpha=0.3)
plt.title('Fibonacci Cosmological Constant: H(z) Fit')
plt.tight_layout()
plt.savefig('hz_comparison_real.png', dpi=150)
plt.show()
