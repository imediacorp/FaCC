# cmb_osc_detector.py
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

phi = (1 + np.sqrt(5)) / 2
ln_phi = np.log(phi)

# Load real low-ℓ CMB
data = np.loadtxt('real_cmb_lowl.csv', delimiter=',', skiprows=1)
ell, cl_data, sigma_cl = data[:, 0], data[:, 1], data[:, 2]

# Polynomial baseline
coeffs = np.polyfit(np.log(ell), cl_data, 5)
cl_fit = np.polyval(coeffs, np.log(ell))
residuals = cl_data - cl_fit

def osc_model(ell, amp, phase):
    return amp * np.cos(2 * np.pi * np.log(ell) / ln_phi + phase)

popt, pcov = curve_fit(osc_model, ell, residuals, sigma=sigma_cl,
                       p0=[np.std(residuals), 0], absolute_sigma=True, maxfev=5000)
amp, phase = popt
unc = np.sqrt(np.diag(pcov))
signif = abs(amp) / unc[0] if unc[0] > 0 else 0

print("="*60)
print("CMB LOG-PERIODIC FIT")
print("="*60)
print(f"Amplitude: {amp:.2f} ± {unc[0]:.2f} μK²  ({signif:.1f}σ)")

# Plot
plt.figure(figsize=(9,5))
plt.errorbar(ell, residuals, yerr=sigma_cl, fmt='o', alpha=0.7, label='Residuals')
plt.plot(ell, osc_model(ell, *popt), 'r-', lw=2, label=f'φ-oscillation ({signif:.1f}σ)')
plt.xscale('log')
plt.xlabel('ℓ'); plt.ylabel('ΔC_ℓ [μK²]')
plt.title('Search for Fibonacci Harmonics in CMB')
plt.legend(); plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig('cmb_osc_real.png', dpi=150)
plt.show()
