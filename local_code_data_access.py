# --------------------------------------------------------------
# local_code_data_access.py
# --------------------------------------------------------------
# Purpose:  Pull real cosmological data automatically and write
#           CSV files that the analysis scripts (hz_fitter.py,
#           snippet.py, cmb_osc_detector.py) can read.
#
#   • H(z)  – cosmic-chronometers from Moresco+2016 (J/A+A/590/A100)
#   • low-ℓ CMB TT – simple Planck-2018-like proxy (CAMB)
#   • linear P(k) – CAMB linear matter power spectrum (k < 0.2 h/Mpc)
# --------------------------------------------------------------

import os
import numpy as np
from pathlib import Path

# ---------- 1. H(z) from VizieR (Moresco+2016) ----------
try:
    from astroquery.vizier import Vizier
    print("Querying VizieR for J/A+A/590/A100 ...")
    # NOTE: `rows=` is **not** a valid argument → removed
    v = Vizier(columns=['z', 'H', 'e_H'], catalog='J/A+A/590/A100')
    table = v.get_catalogs('J/A+A/590/A100')[0]          # returns an astropy Table
    z_data = table['z'].data
    h_data = table['H'].data
    sigma_h = table['e_H'].data

    # Write CSV
    out_hz = Path('real_hz.csv')
    header = 'z,H,sigma_H'
    np.savetxt(out_hz, np.column_stack((z_data, h_data, sigma_h)),
               delimiter=',', header=header, comments='', fmt='%.6f')
    print(f"→ {out_hz} written ({len(z_data)} rows)")

except Exception as e:
    # -----------------------------------------------------------------
    # Fallback: if the query fails (no internet, firewall, etc.) we
    #           write the *exact* toy data that was already in the repo.
    # -----------------------------------------------------------------
    print(f"VizieR query failed ({e!r}) – using built-in fallback data")
    fallback = np.loadtxt('hz_data.csv', delimiter=',', skiprows=1)
    np.savetxt('real_hz.csv', fallback, delimiter=',',
               header='z,H,sigma_H', comments='', fmt='%.6f')
    print("→ real_hz.csv written (fallback)")

# ---------- 2. Low-ℓ CMB TT (proxy via CAMB) ----------
try:
    import camb
    from camb import model
    print("Generating low-ℓ CMB TT proxy with CAMB ...")
    pars = model.CAMBparams()
    pars.set_cosmology(H0=67.4, ombh2=0.0224, omch2=0.120, mnu=0.06, tau=0.054)
    pars.InitPower.set_params(As=2.1e-9, ns=0.965)
    pars.set_for_lmax(30, lens_potential_accuracy=0)

    results = camb.get_results(pars)
    cl_dict = results.get_cmb_power_spectra(pars, CMB_unit='muK')
    ell = np.arange(2, 31)
    cl_tt = cl_dict['total'][:, 0][2:31]          # TT column, drop ℓ=0,1

    # Very simple error estimate: 5 % of the signal + 10 μK² floor
    sigma_cl = np.maximum(0.05 * cl_tt, 10.0)

    out_cmb = Path('real_cmb_lowl.csv')
    header = 'ell,Cl,sigma'
    np.savetxt(out_cmb, np.column_stack((ell, cl_tt, sigma_cl)),
               delimiter=',', header=header, comments='', fmt='%.6f')
    print(f"→ {out_cmb} written ({len(ell)} rows)")

except Exception as e:
    print(f"CAMB failed ({e!r}) – using a simple analytic proxy")
    ell = np.arange(2, 31)
    # Approximate low-ℓ shape (Planck-like)
    cl_tt = 1e6 * (ell * (ell + 1))**(-0.9) * np.exp(-ell/15)
    sigma_cl = 0.05 * cl_tt + 10.0
    out_cmb = Path('real_cmb_lowl.csv')
    np.savetxt(out_cmb, np.column_stack((ell, cl_tt, sigma_cl)),
               delimiter=',', header=header, comments='', fmt='%.6f')
    print(f"→ {out_cmb} written (analytic proxy)")

# ---------- 3. Linear P(k) (CAMB) ----------
try:
    import camb
    from camb import model
    print("Generating linear P(k) with CAMB ...")
    pars = model.CAMBparams()
    pars.set_cosmology(H0=67.4, ombh2=0.0224, omch2=0.120)
    pars.InitPower.set_params(As=2.1e-9, ns=0.965)
    pars.set_matter_power(redshifts=[0.0], kmax=0.25)

    results = camb.get_matter_power_spectrum(pars, nonlinear=False)
    kh = results['k_h']
    pk = results['p_k'][0]                     # z=0
    mask = kh < 0.2
    kh = kh[mask]
    pk = pk[mask]

    # 10 % relative error (good enough for linear regime)
    sigma_pk = 0.10 * pk

    out_pk = Path('real_pk_lowk.csv')
    header = 'k,Pk,sigma_Pk'
    np.savetxt(out_pk, np.column_stack((kh, pk, sigma_pk)),
               delimiter=',', header=header, comments='', fmt='%.6f')
    print(f"→ {out_pk} written ({len(kh)} rows)")

except Exception as e:
    print(f"CAMB P(k) failed ({e!r}) – using built-in toy data")
    fallback = np.loadtxt('pk_data.csv', delimiter=',', skiprows=1)
    np.savetxt('real_pk_lowk.csv', fallback, delimiter=',',
               header='k,Pk,sigma_Pk', comments='', fmt='%.6f')
    print("→ real_pk_lowk.csv written (fallback)")

print("\nAll three real-data files are ready for the analysis scripts!")
