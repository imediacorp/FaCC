import numpy as np
from scipy.signal import find_peaks
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM

# Constants
phi = (1 + np.sqrt(5)) / 2

# Load data
data = np.loadtxt('real_pk_lowk.csv', delimiter=',', skiprows=1)  # k, Pk, sigma_Pk
k_data, pk_data, sigma_pk = data[:,0], data[:,1], data[:,2]

# Expected φ-scales (e.g., around BAO k~0.02 h/Mpc)
k_bao = 0.02  # Example
phi_scales = k_bao * np.array([phi**n for n in range(-5, 6)])  # Multiples

# METHOD 1: Direct peak finding
peaks, _ = find_peaks(pk_data, height=sigma_pk.mean()*3)  # 3σ threshold
k_peaks = k_data[peaks]

# METHOD 2: Analyze oscillations around smooth fit
# Fit smooth baseline to remove the monotonic trend
log_k = np.log(k_data)
log_pk = np.log(pk_data)
spline = UnivariateSpline(log_k, log_pk, s=len(k_data)*0.5)  # Smooth fit
pk_smooth = np.exp(spline(log_k))

# Calculate residuals (oscillations)
residuals = pk_data - pk_smooth

# Find peaks in residuals
residual_peaks, properties = find_peaks(residuals, height=0, prominence=sigma_pk.mean()*0.5)
k_residual_peaks = k_data[residual_peaks]

print("=" * 60)
print("ANALYSIS RESULTS")
print("=" * 60)

# Check direct peaks
if len(k_peaks) > 0:
    matches = [np.min(np.abs(k_peaks - ks)) < 0.005 for ks in phi_scales]
    print(f"\nDirect peaks found: {len(k_peaks)}")
    print(f"φ-scale matches: {sum(matches)} / {len(phi_scales)}")
    print(f"Peak locations (k): {k_peaks}")
else:
    print("\nNo direct peaks found. P(k) appears monotonically decreasing.")

# Check residual peaks
if len(k_residual_peaks) > 0:
    residual_matches = [np.min(np.abs(k_residual_peaks - ks)) < 0.005 for ks in phi_scales]
    print(f"\nResidual oscillation peaks found: {len(k_residual_peaks)}")
    print(f"φ-scale matches in residuals: {sum(residual_matches)} / {len(phi_scales)}")
    print(f"Residual peak locations (k): {k_residual_peaks}")
else:
    print("\nNo significant oscillations detected around smooth fit.")

print("\nExpected φ-scales:")
for i, ks in enumerate(phi_scales):
    print(f"  φ^{i-5} × k_BAO = {ks:.6f} h/Mpc")

# Create comprehensive plot
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))

# Top panel: Full power spectrum
ax1.errorbar(k_data, pk_data, yerr=sigma_pk, fmt='.', label='Data', alpha=0.7)
ax1.plot(k_data, pk_smooth, 'g-', label='Smooth fit', linewidth=2, alpha=0.7)
for ks in phi_scales:
    ax1.axvline(ks, ls='--', color='r', alpha=0.3)
if len(k_peaks) > 0:
    ax1.scatter(k_peaks, pk_data[peaks], color='orange', s=100, marker='*',
                label='Direct peaks', zorder=5)
ax1.set_xscale('log')
ax1.set_xlabel('k [h/Mpc]')
ax1.set_ylabel('P(k) [(Mpc/h)³]')
ax1.set_title('Matter Power Spectrum')
ax1.legend()
ax1.grid(True, alpha=0.3)

# Bottom panel: Residuals showing oscillations
ax2.errorbar(k_data, residuals, yerr=sigma_pk, fmt='.', label='Residuals', alpha=0.7)
ax2.axhline(0, color='k', linestyle='-', linewidth=0.5)
for ks in phi_scales:
    ax2.axvline(ks, ls='--', color='r', alpha=0.3, label='φ-scales' if ks == phi_scales[0] else '')
if len(k_residual_peaks) > 0:
    ax2.scatter(k_residual_peaks, residuals[residual_peaks], color='blue', s=100,
                marker='*', label='Residual peaks', zorder=5)
ax2.set_xscale('log')
ax2.set_xlabel('k [h/Mpc]')
ax2.set_ylabel('ΔP(k) [(Mpc/h)³]')
ax2.set_title('Oscillations Around Smooth Fit')
ax2.legend()
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('pk_phi.png', dpi=150)
print("\nPlot saved as 'pk_phi.png'")
plt.show()
