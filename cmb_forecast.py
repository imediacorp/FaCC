# cmb_forecast.py
import numpy as np
import matplotlib.pyplot as plt

phi = (1 + np.sqrt(5)) / 2
ell = np.arange(2, 2500, dtype=float)  # ← float array
A = 0.015  # 1.5% modulation

# Standard acoustic peaks (simplified)
ell_term = np.power(ell * (ell + 1), -1.0)  # ← np.power avoids integer error
damping = np.exp(-ell / 200.0)

cl_standard = 1e6 * ell_term * damping * (
    1 + 5 * np.exp(-((ell - 220) / 50)**2) +
      3 * np.exp(-((ell - 550) / 60)**2) +
      2 * np.exp(-((ell - 800) / 70)**2)
)

# Add phi-oscillation
osc = 1 + A * np.cos(2 * np.pi * np.log(ell / 220) / np.log(phi))
cl_fib = cl_standard * osc

plt.figure(figsize=(10,6))
plt.plot(ell, cl_standard, 'k-', lw=1.5, label='ΛCDM (CAMB-like)')
plt.plot(ell, cl_fib, 'r-', lw=1.5, label='Fibonacci Perturbations')
for n in range(-2, 3):
    l_phi = 220 * phi**n
    if 2 < l_phi < 2500:
        plt.axvline(l_phi, color='r', ls='--', alpha=0.6)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'Multipole $\ell$')
plt.ylabel(r'$\ell(\ell+1)C_\ell / 2\pi$ [$\mu$K$^2$]')
plt.xlim(2, 2500)
plt.legend()
plt.grid(alpha=0.3, which='both')
plt.title('CMB-S4 Forecast: Log-Periodic Fibonacci Modulation')
plt.tight_layout()
plt.savefig('cmb_forecast.png', dpi=150)
plt.show()
