# hz_forecast.py
import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM

phi = (1 + np.sqrt(5)) / 2
z = np.linspace(0, 3, 500)

# ΛCDM baseline
cosmo = FlatLambdaCDM(H0=68, Om0=0.3)
h_lcdm = cosmo.H(z).value

# Fibonacci damping: H(z) = H_LCDM * (1 + B * cos(log(1+z) / ln_phi))
B = 0.03
h_fib = h_lcdm * (1 + B * np.cos(2 * np.pi * np.log(1 + z + 1e-3) / np.log(phi)))

plt.figure(figsize=(9,6))
plt.plot(z, h_lcdm, 'k-', lw=2, label='ΛCDM')
plt.plot(z, h_fib, 'r-', lw=2, label='Fibonacci Damping (B=0.03)')
for n in range(-1, 3):
    z_phi = np.exp(n * np.log(phi)) - 1
    if 0 < z_phi < 3:
        plt.axvline(z_phi, color='r', ls='--', alpha=0.6)
plt.xlabel('Redshift z')
plt.ylabel('H(z) [km/s/Mpc]')
plt.legend()
plt.grid(alpha=0.3)
plt.title('JWST Forecast: Recursive Damping in H(z)')
plt.tight_layout()
plt.savefig('hz_forecast.png', dpi=150)
plt.show()
