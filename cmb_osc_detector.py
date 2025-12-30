# cmb_osc_detector.py  –  dual‑mode CMB oscillation detector
# -------------------------------------------------------------
# 2025‑09‑05  (C)  Your Name / <your‑institution>

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import corner
import emcee
from scipy.linalg import inv, det
import os, sys

# ----------------------------------------------------------------------
# 1. Data loading -------------------------------------------------------
# ----------------------------------------------------------------------
def _maybe_download_planck_r3(base, dest_dir):
    """
    Download text spectra and covariance if `base` is an HTTP(S) URL.
    The expected text files are:
      - Cl_EE.dat, Cl_BB.dat, Cl_TE.dat (two columns: ell, Cl_ell)
      - cov_R3.dat (square matrix text, N_ell x N_ell)

    This helper enables easy use on Google Colab by pointing to an online
    dataset (e.g., a GitHub raw/Zenodo/LAMBDA link prefix). Files are cached
    in `dest_dir`.
    """
    import urllib.request

    os.makedirs(dest_dir, exist_ok=True)
    filenames = ["Cl_EE.dat", "Cl_BB.dat", "Cl_TE.dat", "cov_R3.dat"]
    downloaded = []
    for fn in filenames:
        local_path = os.path.join(dest_dir, fn)
        if os.path.exists(local_path) and os.path.getsize(local_path) > 0:
            downloaded.append(local_path)
            continue
        url = base.rstrip("/") + "/" + fn
        try:
            print(f"Downloading {url} → {local_path}")
            urllib.request.urlretrieve(url, local_path)
            downloaded.append(local_path)
        except Exception as e:
            print(f"WARNING: Could not download {url}: {e}")
    return all(os.path.exists(os.path.join(dest_dir, fn)) for fn in filenames)


def load_planck_r3(cov_dir="data/planck_r3", url_base=None):
    """
    Load the Planck R3-like angular power spectra and a covariance matrix.

    Supports three modes:
      1) Local directory with text files: `cov_dir` points to a folder
         containing `Cl_EE.dat`, `Cl_BB.dat`, `Cl_TE.dat`, `cov_R3.dat`.
      2) Remote base URL: pass `url_base` (HTTP/HTTPS). Files will be downloaded
         to `cov_dir` (created if needed) and cached.
      3) Auto online: set `cov_dir="auto"` to download from default URLs
         defined below (recommended for Google Colab).

    Returns:
        ell : np.ndarray   (ℓ values)
        Cl_dict : dict with keys {'EE','BB','TE'} each np.ndarray(ℓ)
        cov : np.ndarray   (full covariance, shape (Nℓ, Nℓ))
    """
    # Default online dataset (ASCII) compatible with this loader.
    # You can replace these with your own mirrors. The below URLs point to
    # a small prepared set hosted as plain text. If they become unavailable,
    # pass your own `url_base` that hosts the four required files.
    DEFAULT_URL_BASE = os.environ.get(
        "PLANCK_R3_URL_BASE",
        # Example placeholder: replace with your own mirror if needed
        "https://raw.githubusercontent.com/jetbrains-research/sample-cosmology-datasets/main/planck_r3_ascii"
    )

    # Decide source
    if cov_dir == "auto" and url_base is None:
        url_base = DEFAULT_URL_BASE
        cov_dir = os.path.join("data", "planck_r3")

    # If url_base is an http(s) link, try downloading
    if url_base and url_base.startswith(("http://", "https://")):
        ok = _maybe_download_planck_r3(url_base, cov_dir)
        if not ok:
            print("WARNING: Online download incomplete. Will try to proceed or fall back.")

    # Try to read local text files
    try:
        Cl_EE = np.loadtxt(os.path.join(cov_dir, "Cl_EE.dat"))
        Cl_BB = np.loadtxt(os.path.join(cov_dir, "Cl_BB.dat"))
        Cl_TE = np.loadtxt(os.path.join(cov_dir, "Cl_TE.dat"))
        ell = Cl_EE[:, 0].astype(int)
        Cl_dict = {'EE': Cl_EE[:, 1], 'BB': Cl_BB[:, 1], 'TE': Cl_TE[:, 1]}
        cov_path = os.path.join(cov_dir, "cov_R3.dat")
        if os.path.exists(cov_path):
            cov = np.loadtxt(cov_path)
        else:
            # Fallback: diagonal covariance with 10% fractional error
            print("WARNING: cov_R3.dat not found. Using diagonal covariance fallback (results are not for publication).")
            N = len(ell)
            std = 0.1 * np.maximum(1e-12, np.abs(Cl_dict['EE']))
            cov = np.diag(std[:N]**2)
        return ell, Cl_dict, cov
    except Exception as e:
        # As a last resort, try to load from FITS if provided (requires astropy)
        print(f"Text data load failed: {e}")
        print("Attempting FITS-based load (requires astropy). Set url_base to a directory containing FITS tables.")
        try:
            from astropy.io import fits
            # Expected FITS names if hosted similarly
            ee_path = os.path.join(cov_dir, "Cl_R3_EE.fits")
            bb_path = os.path.join(cov_dir, "Cl_R3_BB.fits")
            te_path = os.path.join(cov_dir, "Cl_R3_TE.fits")
            cov_path = os.path.join(cov_dir, "COV_R3.fits")
            with fits.open(ee_path) as hdul:
                data = hdul[1].data
                ell = np.array(data["ELL"]) if "ELL" in data.names else np.arange(len(data))
                EE = np.array(data.field(1))
            with fits.open(bb_path) as hdul:
                BB = np.array(hdul[1].data.field(1))
            with fits.open(te_path) as hdul:
                TE = np.array(hdul[1].data.field(1))
            with fits.open(cov_path) as hdul:
                cov = np.array(hdul[1].data)
            Cl_dict = {"EE": EE, "BB": BB, "TE": TE}
            return ell, Cl_dict, cov
        except Exception as e2:
            raise FileNotFoundError(
                "Could not load Planck R3 data. Provide local text files, set url_base to a valid HTTP directory containing the four .dat files, or place FITS tables with the expected names."
            ) from e2

# ----------------------------------------------------------------------
# 2. Duality model (ΔCℓ) ------------------------------------------------
# ----------------------------------------------------------------------
def delta_Cl(ell, Aphi=1.0, phi_phi=0.0, Ahatphi=0.5, phi_hatphi=0.0,
             use_conjugate=True, baseline=0.0):
    """
    Compute the dual‑oscillation correction ΔCℓ = Aφ cos(α) + Aĥφ
cos(-α+Δ)
    where α = 2π log ℓ / ln φ, Δ = φφ - φĥφ.
    Parameters
    ----------
    ell : array_like
        Multipole moments.
    Aphi, phi_phi : float
        Amplitude and phase of the forward mode.
    Ahatphi, phi_hatphi : float
        Amplitude and phase of the conjugate mode.
    use_conjugate : bool
        If False, the second term is omitted (single‑mode case).
    baseline : float
        Optional additive offset (nuisance) to absorb foreground/calibr.
    """
    # Golden ratio and its reciprocal
    phi = (1 + np.sqrt(5)) / 2.0
    ln_phi = np.log(phi)
    # Oscillation argument
    alpha = 2.0 * np.pi * np.log(ell) / ln_phi
    # Forward term
    dC = Aphi * np.cos(alpha + phi_phi)
    if use_conjugate:
        # Using cos(−α + φĥφ) = cos(α - φĥφ)
        dC += Ahatphi * np.cos(alpha - phi_hatphi)
    return dC + baseline

# ----------------------------------------------------------------------
# 3. Likelihood ---------------------------------------------------------
# ----------------------------------------------------------------------
def log_likelihood(theta, ell, Cl_dict, cov, spectrum='EE',
use_conjugate=True):
    """
    Gaussian log‑likelihood with full covariance.
    theta = [Aphi, phi_phi, Ahatphi, phi_hatphi, baseline]
    """
    Aphi, phi_phi, Ahatphi, phi_hatphi, baseline = theta
    # Predict full theory: Cl_theory = Cl_smooth + ΔCℓ
    # Here we treat the Planck R3 spectrum as the observed smooth baseline.
    Cl_obs = Cl_dict[spectrum]
    # Compute ΔCℓ
    dCl = delta_Cl(ell, Aphi, phi_phi, Ahatphi, phi_hatphi,
                   use_conjugate=use_conjugate, baseline=baseline)
    Cl_th = Cl_obs + dCl
    # Residual vector
    resid = Cl_th - Cl_obs
    # Compute χ² = rᵀ Σ⁻¹ r (use scipy.linalg.inv for speed on modest N_ell)
    invcov = inv(cov)
    chi2 = resid @ invcov @ resid
    # Log‑det term (constant for fixed Σ)
    logdet = np.log(det(cov))
    return -0.5 * (chi2 + logdet)

# ----------------------------------------------------------------------
# 4. Priors (flat) -------------------------------------------------------
# ----------------------------------------------------------------------
def log_prior(theta, use_conjugate=True):
    Aphi, phi_phi, Ahatphi, phi_hatphi, baseline = theta
    # Flat priors in reasonable ranges
    if -5.0 < Aphi < 5.0 and -np.pi < phi_phi < np.pi:
        if use_conjugate:
            if -5.0 < Ahatphi < 5.0 and -np.pi < phi_hatphi < np.pi:
                if -1.0 < baseline < 1.0:
                    return 0.0
        else:
            if -1.0 < baseline < 1.0:
                return 0.0
    return -np.inf

def log_posterior(theta, ell, Cl_dict, cov, spectrum='EE',
use_conjugate=True):
    lp = log_prior(theta, use_conjugate)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, ell, Cl_dict, cov,
                               spectrum=spectrum,
use_conjugate=use_conjugate)

# ----------------------------------------------------------------------
# 5. MCMC driver ---------------------------------------------------------
# ----------------------------------------------------------------------
def run_mcmc(ell, Cl_dict, cov, nwalkers=64, nsteps=5000,
             spectrum='EE', use_conjugate=True, out_dir="chains"):
    os.makedirs(out_dir, exist_ok=True)
    # Parameter order: Aphi, phi_phi, Ahatphi, phi_hatphi, baseline
    ndim = 5 if use_conjugate else 3
    # Initialize walkers in a tiny Gaussian ball around a guess
    p0 = np.array([0.5, 0.0, 0.2, 0.0, 0.0])  # <- adjust as needed
    # For single‑mode we only need Aphi, phi_phi, baseline
    p0_single = np.array([0.5, 0.0, 0.0])
    p0 = p0_single if not use_conjugate else p0

    # emcee.EnsembleSampler
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior,
                                    args=(ell, Cl_dict, cov, spectrum,
use_conjugate))

    # Burn‑in
    print("Burn‑in...")
    pos, prob, state = sampler.run_mcmc(p0 + 1e-3*np.random.randn(nwalkers, ndim), nsteps//2)
    sampler.reset()

    # Production
    print("Production...")
    sampler.run_mcmc(pos, nsteps, progress=True)

    # Save chain
    chain = sampler.get_chain(flat=False)   # (nsteps, nwalkers, ndim)
    np.save(os.path.join(out_dir, f"chain_{'dual' if use_conjugate else 'single'}.npy"), chain)
    return sampler

# ----------------------------------------------------------------------
# 6. Post‑processing -----------------------------------------------------
# ----------------------------------------------------------------------
def plot_corner(chain, labels=None, fname="corner.png"):
    if labels is None:
        labels = [r"$A_\phi$", r"$\phi_\phi$", r"$A_{\hat\phi}$",
                  r"$\phi_{\hat\phi}$", r"$b$"][:chain.shape[-1]]
    fig = corner.corner(chain, labels=labels, quantiles=[0.16,0.5,0.84],
                        show_titles=True, title_fmt=".3f")
    fig.savefig(fname, dpi=300)

def plot_residuals(ell, Cl_obs, dCl_best, spectrum='EE',
fname="residuals.png"):
    plt.figure(figsize=(8,5))
    plt.plot(ell, dCl_best, 'C0-', label="ΔCℓ (dual fit)")
    plt.plot(ell, Cl_obs, 'C1-', label="Observed")
    plt.xlabel(r"$\ell$")
    plt.ylabel(r"$C_\ell$ [$\mu$K$^2$]")
    plt.title(f"Residuals for {spectrum}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(fname, dpi=300)
    plt.show()

# ----------------------------------------------------------------------
# 7. Main ---------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__":
    # Example usage (replace paths with your own)
    # Colab-friendly: honor environment variables to auto-download from the web
    #   - Set PLANCK_R3_URL_BASE to an HTTP(S) directory with the four .dat files
    #   - Optionally set PLANCK_R3_DIR for the local cache directory
    ell, Cl_dict, cov = load_planck_r3(
        os.environ.get("PLANCK_R3_DIR", "data/planck_r3"),
        url_base=os.environ.get("PLANCK_R3_URL_BASE")
    )
    spec = "EE"  # Choose spectrum to fit: 'EE', 'BB', or 'TE'

    # Run dual‑mode MCMC
    sampler_dual = run_mcmc(
        ell, Cl_dict, cov,
        nwalkers=96, nsteps=4000,
        spectrum=spec, use_conjugate=True,
        out_dir="chains/dual",
    )

    # Run single‑mode MCMC (no conjugate term) for Bayes factor
    sampler_single = run_mcmc(
        ell, Cl_dict, cov,
        nwalkers=96, nsteps=4000,
        spectrum=spec, use_conjugate=False,
        out_dir="chains/single",
    )

    # Pull out posterior samples
    chain_dual = sampler_dual.get_chain(flat=True)
    chain_single = sampler_single.get_chain(flat=True)

    # Posterior means (or MAP)
    theta_dual = np.mean(chain_dual, axis=0)
    theta_single = np.mean(chain_single, axis=0)

    # Plot corner for the dual fit
    plot_corner(
        chain_dual,
        labels=[r"$A_\phi$", r"$\phi_\phi$", r"$A_{\hat\phi}$", r"$\phi_{\hat\phi}$", r"$b$"],
        fname="corner_dual.png",
    )

    # Compute best‑fit ΔCℓ and residuals
    dCl_dual = delta_Cl(ell, *theta_dual, use_conjugate=True, baseline=theta_dual[-1])
    plot_residuals(ell, Cl_dict[spec], dCl_dual, spectrum=spec, fname="residuals.png")

    # Rough Bayes factor estimate via a simple harmonic mean approximation (not ideal, but a start)
    def harmonic_mean_est(samples, logp):
        return np.mean(np.exp(logp(samples)))

    # Define a single‑mode log posterior wrapper (enforce Ahatphi=0)
    def log_posterior_single(theta3):
        return log_posterior(np.append(theta3, 0.0), ell, Cl_dict, cov, spectrum=spec, use_conjugate=False)

    # Using a subset of the chains to keep the estimation cheap
    idx = np.random.choice(len(chain_dual), size=min(2000, len(chain_dual)), replace=False)
    sub_chain_dual = chain_dual[idx]
    sub_chain_single = chain_single[idx]

    Z_dual = harmonic_mean_est(
        sub_chain_dual,
        lambda th: log_posterior(th, ell, Cl_dict, cov, spectrum=spec, use_conjugate=True),
    )
    Z_single = harmonic_mean_est(
        sub_chain_single,
        lambda th: log_posterior_single(th),
    )

    BF = Z_dual / Z_single if Z_single != 0 else np.inf
    print(f"Bayes factor (dual / single) ≈ {BF:.2f}  → log10(BF) ≈ {np.log10(BF):.2f}")