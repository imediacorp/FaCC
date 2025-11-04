# snippet.py
import numpy as np
from scipy.signal import find_peaks
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt

phi = (1 + np.sqrt(5)) / 2

# Load real P(k)
data = np.loadtxt('real_pk_lowk.csv', delimiter=',', skiprows=1)
k_data, pk_data, sigma_pk = data[:, 0], data[:, 1], data[:, 2]

k_bao = 0.1
phi_scales = k_bao * np.array([phi**n for n in range(-5, 6)])

# Detrend
log_k = np.log(k_data)
log_pk = np.log(pk_data)
spline = UnivariateSpline(log_k, log_pk, s=len(k_data)*0.5)
pk_smooth = np.exp(spline(log_k))
residuals = pk_data - pk_smooth

# Find peaks in residuals
residual_peaks, _ = find_peaks(np.abs(residuals), height=0, prominence=sigma_pk.mean()*0.6)
k_peaks = k_data[residual_peaks]

# Match to φ-scales
tolerance = 0.008
matches = [any(abs(kp - ks) < tolerance for ks in phi_scales) for kp in k_peaks]
n_matches = sum(matches)

print("="*60)
print("P(k) φ-SCALE SEARCH")
print("="*60)
print(f"Residual peaks found: {len(k_peaks)}")
print(f"φ-scale matches: {n_matches}/11")
if n_matches > 0:
    matched_ks = [k_peaks[i] for i, m in enumerate(matches) if m]
    print(f"Matched k: {matched_ks}")

# -------------------------------------------------
# PLOT (fixed & polished)
# -------------------------------------------------
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), constrained_layout=True)

# Top panel – full P(k)
ax1.errorbar(k_data, pk_data, yerr=sigma_pk, fmt='.', alpha=0.7, label='CAMB P(k)')
ax1.plot(k_data, pk_smooth, 'g-', lw=2, label='Smooth fit')
for ks in phi_scales:
    ax1.axvline(ks, color='r', ls='--', alpha=0.5)
ax1.set_xscale('log')
ax1.set_xlabel('k [h/Mpc]')
ax1.set_ylabel('P(k) [(Mpc/h)$^3$]')
ax1.set_title('Matter Power Spectrum')
ax1.legend()
ax1.grid(alpha=0.3)

# Bottom panel – residuals
ax2.errorbar(k_data, residuals, yerr=sigma_pk, fmt='.', alpha=0.7, label='Residuals')
ax2.axhline(0, color='k', lw=0.8)
for ks in phi_scales:
    ax2.axvline(ks, color='r', ls='--', alpha=0.5)
if len(k_peaks) > 0:
    ax2.scatter(k_peaks, residuals[residual_peaks],
                color='orange', s=80, zorder=5,
                label=f'{n_matches} φ-match{"es" if n_matches!=1 else ""}')
ax2.set_xscale('log')
ax2.set_xlabel('k [h/Mpc]')
ax2.set_ylabel(r'$\Delta P(k)$ [(Mpc/h)$^3$]')
ax2.set_title('Residual Oscillations')
ax2.legend()
ax2.grid(alpha=0.3)

plt.suptitle('Search for Golden-Ratio Scales in $P(k)$', fontsize=14)
plt.savefig('pk_phi_real.png', dpi=150)
plt.show()
