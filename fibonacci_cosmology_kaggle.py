#%% md
"""
# Fibonacci Cosmology: Kaggle-Ready Notebook (Script Version)
**Author:** Bryan David Persaud  
**Affiliation:** Intermedia Communications Corp.  

This Kaggle-ready script-notebook reproduces the core exploratory analyses from the repository:
- Load H(z), low-ℓ CMB TT, and low-k P(k) datasets from GitHub (with offline fallback to local CSVs)
- Plot basic diagnostics and residuals used in the Fibonacci cosmology perturbation tests

Notes for Kaggle usage:
- If Internet is enabled, data are fetched from the GitHub repo. If not, place the CSV files in the working directory (e.g., upload as Kaggle Dataset) and the loader will fall back to local files.
- This script uses the "percent" notebook format (#%% cells). You can run it as a Kaggle "Notebook" (Python Script) or convert to a Jupyter notebook using jupytext locally if desired.
"""

#%%
# ============================================================================
# SETUP
# ============================================================================
import sys
import warnings
warnings.filterwarnings('ignore')

# Optional: lightweight installs if needed on Kaggle
try:
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from scipy.optimize import curve_fit
except Exception:
    # Fallback installs (Kaggle will cache these in sessions)
    import subprocess, sys
    subprocess.run([sys.executable, '-m', 'pip', 'install', '-q', 'numpy', 'pandas', 'matplotlib', 'scipy'])
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from scipy.optimize import curve_fit

#%%
# ============================================================================
# CONSTANTS AND HELPERS
# ============================================================================
phi = (1 + 5 ** 0.5) / 2
ln_phi = np.log(phi)
phi_conj = phi - 1
ln_phi_conj = np.log(phi_conj)

print(f"Golden Ratio φ = {phi:.8f}")
print(f"Conjugate φ̂ = {phi_conj:.8f}")
print(f"ln(φ) = {ln_phi:.8f}")
print(f"ln(φ̂) = {ln_phi_conj:.8f}")

BASE_URL = "https://raw.githubusercontent.com/imediacorp/FaCC/main/"

def load_csv(name, local_fallback=True):
    """Try GitHub first; fall back to local if disabled/no internet."""
    import io, urllib.request
    url = BASE_URL + name
    try:
        with urllib.request.urlopen(url, timeout=10) as r:
            data = r.read()
        df = pd.read_csv(io.BytesIO(data))
        print(f"✓ Loaded from GitHub: {name} ({len(df)} rows)")
        return df
    except Exception as e:
        if local_fallback:
            try:
                df = pd.read_csv(name)
                print(f"✓ Loaded locally: {name} ({len(df)} rows)")
                return df
            except Exception as e2:
                print(f"✗ Could not load {name} from GitHub or local: {e2}")
                return None
        else:
            print(f"✗ Could not load {name} from GitHub: {e}")
            return None

#%%
# ============================================================================
# LOAD DATASETS
# ============================================================================
hz = load_csv('real_hz.csv')
cmb = load_csv('real_cmb_lowl.csv')
pk  = load_csv('real_pk_lowk.csv')

# Basic integrity checks
for name, df in [('H(z)', hz), ('CMB low-ℓ', cmb), ('P(k)', pk)]:
    if df is None:
        print(f"Warning: {name} missing; downstream plots will be skipped.")

#%% md
"""
## Quick-Look Plots
These reproduce the minimal visual diagnostics used in the repository. For publication-level figures and statistical tests, see scripts like `hz_fitter_v2.py`, `cmb_osc_detector_v2.py`, and `lss_phi_analyzer.py` in the repo.
"""

#%%
# H(z) quick plot
if hz is not None:
    z = hz[hz.columns[0]].values
    H = hz[hz.columns[1]].values
    s = hz[hz.columns[2]].values

    plt.figure(figsize=(6,4))
    plt.errorbar(z, H, yerr=s, fmt='o', ms=3, alpha=0.8, label='H(z) data')
    plt.xlabel('z')
    plt.ylabel('H(z) [km/s/Mpc]')
    plt.title('Cosmic Chronometer H(z) Data')
    plt.grid(True, alpha=0.2)
    plt.legend()
    plt.tight_layout()
    plt.show()

#%%
# CMB low-ℓ quick plot
if cmb is not None:
    ell = cmb[cmb.columns[0]].values
    C   = cmb[cmb.columns[1]].values
    sC  = cmb[cmb.columns[2]].values

    plt.figure(figsize=(6,4))
    plt.errorbar(ell, C, yerr=sC, fmt='.', ms=4, alpha=0.8, label='TT low-ℓ')
    plt.xlabel('ℓ')
    plt.ylabel('C_ℓ')
    plt.title('CMB TT (low-ℓ)')
    plt.grid(True, alpha=0.2)
    plt.legend()
    plt.tight_layout()
    plt.show()

#%%
# P(k) quick plot
if pk is not None:
    kval = pk[pk.columns[0]].values
    Pk   = pk[pk.columns[1]].values
    sPk  = pk[pk.columns[2]].values

    plt.figure(figsize=(6,4))
    plt.errorbar(kval, Pk, yerr=sPk, fmt='o', ms=3, alpha=0.8, label='P(k)')
    plt.xscale('log')
    plt.xlabel('k [h/Mpc]')
    plt.ylabel('P(k)')
    plt.title('Matter Power Spectrum (low-k)')
    plt.grid(True, which='both', alpha=0.2)
    plt.legend()
    plt.tight_layout()
    plt.show()

#%% md
"""
## Optional: Minimal Log-Periodic Residual Fit (Demonstration)
This cell demonstrates a minimal residual fit using a toy ΛCDM-like smooth curve and a cosine in log-k for P(k). It is for illustration only; see `lss_phi_analyzer.py` for a robust version.
"""

#%%
if pk is not None:
    kval = pk[pk.columns[0]].values
    Pk   = pk[pk.columns[1]].values
    sPk  = pk[pk.columns[2]].values

    # Smooth baseline (toy): a power-law with a soft turnover near BAO
    def baseline(k, A, n, k0):
        return A * (k / k0) ** n / (1 + (k / (k0*5))**2)

    # Log-periodic modulation with ln(φ) period
    def modulated(k, A, n, k0, B, phase):
        base = baseline(k, A, n, k0)
        return base * (1.0 + B * np.cos(np.log(k + 1e-9)/ln_phi + phase))

    # Initial guesses
    A0 = np.median(Pk)
    n0 = -1.0
    k00 = 0.02
    B0 = 0.05
    phase0 = 0.0

    try:
        popt, pcov = curve_fit(modulated, kval, Pk, sigma=np.maximum(sPk, 1e-10),
                               p0=[A0, n0, k00, B0, phase0], maxfev=10000)
        A1, n1, k01, B1, phase1 = popt
        fit = modulated(kval, *popt)

        plt.figure(figsize=(6,4))
        plt.plot(kval, Pk, 'o', ms=3, alpha=0.6, label='Data')
        plt.plot(kval, fit, '-', lw=2, label='Toy fit (log-periodic)')
        plt.xscale('log')
        plt.xlabel('k [h/Mpc]')
        plt.ylabel('P(k)')
        plt.title('Toy Log-Periodic Fit (period ~ ln φ)')
        plt.grid(True, which='both', alpha=0.2)
        plt.legend()
        plt.tight_layout()
        plt.show()

        resid = (Pk - fit)
        plt.figure(figsize=(6,3))
        plt.plot(kval, resid, '.', ms=4)
        plt.axhline(0, color='k', lw=1)
        plt.xscale('log')
        plt.xlabel('k [h/Mpc]')
        plt.ylabel('Residual')
        plt.title('Residuals')
        plt.grid(True, which='both', alpha=0.2)
        plt.tight_layout()
        plt.show()

        print(f"Estimated modulation amplitude B ≈ {B1:.4f}")
    except Exception as e:
        print(f"Fit failed (expected if data are sparse or noisy): {e}")

#%% md
"""
## Next Steps
- For full analyses and publication figures, see the scripts in this repo:
  - `hz_fitter_v2.py` (H(z) falsification tests)
  - `cmb_osc_detector_v2.py` (CMB low-ℓ residuals)
  - `lss_phi_analyzer.py` (P(k) residual analysis)
- On Kaggle, you can attach this GitHub repo as a Dataset or enable Internet and use the loader above.
"""
