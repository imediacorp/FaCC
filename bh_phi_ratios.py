import numpy as np
import pandas as pd
from scipy.stats import chisquare

# Load GWTC-3 CSV (user download: masses.csv)
df = pd.read_csv('GWTC3_masses.csv')  # Columns: event, m1_source, m2_source
ratios = df['m1_source'] / df['m2_source']
log_ratios = np.log(ratios)

# Test φ-scaling: expected log(φ^n) for n=-2 to 2
phi = (1 + np.sqrt(5))/2
expected = np.log(phi**np.arange(-2,3))
observed_binned = np.histogram(log_ratios, bins=5)[0]  # Bin for chi2
chi2, p = chisquare(observed_binned, f_exp=np.full(5, len(ratios)/5))
print(f"Chi² φ-fit: {chi2:.2f}, p-value: {p:.3f}")

# Plot histogram
plt.hist(ratios, bins=20, alpha=0.7, label='Observed')
plt.axvline(phi, color='r', ls='--', label='φ')
plt.axvline(1/phi, color='b', ls=':', label='1/φ')
plt.xlabel('m1/m2'); plt.ylabel('Count'); plt.legend()
plt.savefig('bh_ratios.png')
