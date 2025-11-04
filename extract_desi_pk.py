#!/usr/bin/env python
"""
Extract real DESI DR2 P(k) from FITS → CSV for lss_phi_analyzer.py
"""

from astropy.io import fits
import numpy as np
import pandas as pd
from pathlib import Path

# Paths
fits_dir = Path("desi_dr2")
pk_file = fits_dir / "pk_monopole.fits"
cov_file = fits_dir / "covariance.fits"
csv_file = Path("pk_data_real.csv")

# Load
print("Loading P(k) FITS...")
hdul_pk = fits.open(pk_file)
hdul_cov = fits.open(cov_file)

# Extract (adjust column names if needed)
k = hdul_pk[1].data['k']           # h/Mpc
pk = hdul_pk[1].data['pk0']        # (Mpc/h)^3
cov = hdul_cov[1].data             # Full covariance
sigma_pk = np.sqrt(np.diag(cov))   # Diagonal errors

# Save as CSV
df = pd.DataFrame({
    'k': k,
    'Pk': pk,
    'sigma_Pk': sigma_pk
})
df.to_csv(csv_file, index=False)
print(f"Saved {len(k)} points → {csv_file}")

# Optional: Save full covariance
np.save("desi_dr2_cov.npy", cov)
print("Covariance saved → desi_dr2_cov.npy")
