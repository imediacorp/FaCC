# forecast_plots.py
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

phi = (1 + np.sqrt(5)) / 2
k = np.logspace(-2, 0, 1000)
k_bao = 0.1
A = 0.02

# Standard P(k)
pk_standard = 1e5 * k**(-1.5) * (1 + 10 * np.exp(-((k - k_bao)/0.02)**2))

# Add phi-oscillations
osc = 1 + A * np.cos(2 * np.pi * np.log(k / k_bao) / np.log(phi))
pk_fib = pk_standard * osc

plt.figure(figsize=(9,5))
plt.loglog(k, pk_standard, 'k-', label='ΛCDM')
plt.loglog(k, pk_fib, 'r-', label='Fibonacci Perturbations')
for n in range(-3,4):
    plt.axvline(k_bao * phi**n, color='r', ls='--', alpha=0.5)
plt.xlabel('k [h/Mpc]'); plt.ylabel('P(k)')
plt.legend(); plt.grid(alpha=0.3)
plt.title('DESI Y6 Forecast: φ-Oscillations in P(k)')
plt.tight_layout()
plt.savefig('pk_forecast.png', dpi=150)
plt.show()
