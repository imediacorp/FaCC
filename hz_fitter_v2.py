import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM  # For baseline comparison

# Constants from theorem
phi = (1 + np.sqrt(5)) / 2
ln_phi = np.log(phi)
phi_conj = phi - 1
ln_phi_conj = np.log(phi_conj)

# Load data
data = np.loadtxt('real_hz.csv', delimiter=',', skiprows=1)  # Columns: z, H, sigma_H
z_data, h_data, sigma_h = data[:, 0], data[:, 1], data[:, 2]

# Model H(z) [km/s/Mpc] - Recursive with matter (perturbative)
def h_model(z, om, t0, sigma=1.0, h0=70):  # sigma = 1 for forward, -1 for reverse
    ol = 1 - om
    ln_r = ln_phi if sigma > 0 else ln_phi_conj  # Use conj for reverse magnitude
    hubble = (sigma * ln_r / t0) * np.sqrt(om * (1 + z)**3 + ol)  # Sign for direction
    return hubble * 3.08568e19 / 3.15576e16  # Convert to km/s/Mpc

# Chi-squared for forward
def chi2_forward(params):
    om, t0 = params
    h_pred = h_model(z_data, om, t0, sigma=1.0)
    return np.sum(((h_data - h_pred) / sigma_h)**2)

# Chi-squared for reverse (using abs for fit since data is positive)
def chi2_reverse(params):
    om, t0 = params
    h_pred = h_model(z_data, om, t0, sigma=-1.0)
    return np.sum(((h_data - np.abs(h_pred)) / sigma_h)**2)

# Fit forward
res_forward = minimize(chi2_forward, [0.3, 1e-17], bounds=[(0.1, 0.4), (1e-18, 1e-16)])
om_fit_f, t0_fit_f = res_forward.x
chi2_fit_f = res_forward.fun
print(f"Forward Fitted Om: {om_fit_f:.3f}, t0: {t0_fit_f:.2e} s, chi2: {chi2_fit_f:.2f}")

# Fit reverse (exploratory)
res_reverse = minimize(chi2_reverse, [0.3, 1e-17], bounds=[(0.1, 0.4), (1e-18, 1e-16)])
om_fit_r, t0_fit_r = res_reverse.x
chi2_fit_r = res_reverse.fun
print(f"Reverse Fitted Om: {om_fit_r:.3f}, t0: {t0_fit_r:.2e} s, chi2: {chi2_fit_r:.2f}")

# Baseline ΛCDM
cosmo_lcdm = FlatLambdaCDM(H0=70, Om0=0.3)
h_lcdm = cosmo_lcdm.H(z_data).value

# Plot
plt.errorbar(z_data, h_data, yerr=sigma_h, fmt='o', label='Data')
plt.plot(z_data, h_model(z_data, om_fit_f, t0_fit_f, sigma=1.0), label='Fibonacci Forward')
plt.plot(z_data, np.abs(h_model(z_data, om_fit_r, t0_fit_r, sigma=-1.0)), '--', label='Fibonacci Reverse (abs)')
plt.plot(z_data, h_lcdm, ':', label='ΛCDM')
plt.xlabel('Redshift z')
plt.ylabel('H(z) [km/s/Mpc]')
plt.legend()
plt.savefig('hz_comparison_dual.png')
plt.show()
