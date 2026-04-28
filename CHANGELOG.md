# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - 2026-04-28

### Added
- **`src/hz_analysis.py`** — H(z) analysis module: load cosmic-chronometer data,
  fit φ-recursive forward and reverse branches, compare with ΛCDM
- **`src/cmb_analysis.py`** — CMB low-ℓ oscillation module: polynomial ΛCDM
  baseline subtraction, dual log-periodic φ/φ⁻¹ oscillation fit to residuals
- **`src/lss_analysis.py`** — LSS P(k) module: spline baseline, direct and
  residual peak finding, φ-scale feature detection with match counting
- **`src/g_phi_analysis.py`** — G–φ connection module: φ-log-periodic varying-G
  model (pulsar-constrained), SGWB spectrum with φ-spaced ISCO markers,
  GWTC-3 BH mass-ratio KS clustering test (74 events built-in)
- **`src/desi_bao.py`** — Live DESI DR2 BAO integration: fetches DR2 and DR1
  distance ratios from CobayaSampler/bao_data, computes ΛCDM predictions via
  astropy (Aubourg+2015 sound horizon), fits φ-log-periodic residual model
- **Interactive Plotly dashboard** — full rewrite of `streamlit_app.py`;
  9 tabs with zoom/hover/download on all charts (replaces static matplotlib)
  - 🛰 DESI BAO tab: live DR2/DR1 data with ΛCDM overlay and φ-oscillation fit
  - 🌊 G–φ Connection tab: 2×2 figure for G(t), Ġ/G, SGWB, mass-ratio histogram
- **Comprehensive test suite** — 184 tests across 8 files (from ~66);
  new `conftest.py` with shared fixtures; all network/CAMB calls mocked
- **Collaboration infrastructure** — `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`,
  GitHub issue templates (bug, feature), PR template with checklist
- **Example scripts** — `examples/lss_phi_example.py`, `desi_bao_example.py`,
  `g_phi_example.py` with CLI argument support
- **`docs/index.md`** — navigation index for all documentation

### Changed
- Matplotlib imports made lazy (inside `run_*` functions only) — fixes
  `ImportError` on Apple Silicon Macs with x86_64 matplotlib installations
- `setup.py`: corrected `README.md` filename casing (was `Readme.md`, broke
  PyPI packaging on Linux); synced `extras_require` with `requirements-dev.txt`
- `README.md`: added CI/Python/licence badges, Troubleshooting section,
  Contributing link
- Theory docs (`CosmicPerturbationHypothesis.md`, `FibonacciCosmologyAnalysis.md`,
  `CosmologicalTestResultsAnalysis.md`) moved from root into `docs/`

### Removed
- `RELEASE_NOTES_v0.1.0.md` — content consolidated into `CHANGELOG.md`

## [0.1.1] - 2025-01-16

### Added
- `SystematicErrorBudget` class for modeling systematic uncertainties in forecasts
  - Photo-z error modeling
  - Galaxy bias uncertainty propagation
  - Survey geometry effects
  - Comprehensive error budget calculation
- `BayesianEvidence` class for Bayesian model comparison
  - Harmonic mean estimator for evidence calculation
  - Bayes factor computation (φ-modulation vs ΛCDM)
  - BIC computation and interpretation utilities
- `forecast_desi_sensitivity_with_systematics()` method in `PhiModulationModel`
  - Extended forecast method that includes systematic error contributions
  - Returns detailed systematic error breakdown
- `CITATION.cff` file for proper citation metadata
- `examples/systematics_example.py` - Example script demonstrating systematic error analysis
- Updated `src/__init__.py` to export new modules

### Changed
- Enhanced `PhiModulationModel` to integrate with systematic error analysis
- Improved module organization and exports

## [0.1.0] - 2025-01-01

### Added
- `PhiModulationModel` class implementing φ-modulated power spectrum analysis
- CAMB integration for ΛCDM baseline power spectrum generation
- DESI Year 5 forecast analysis with Fisher matrix methodology
- `notebooks/01_desi_forecasts.ipynb` - Comprehensive forecast analysis notebook
- Log-periodic modulation: `P(k) = P_ΛCDM(k) × [1 + A_φ × cos(2π × log(k/k_pivot) / ln(φ) + φ_0)]`
- BAO signature computation with φ-modulation
- Independence statement documentation (`INDEPENDENCE.md`)
- Updated `requirements.txt` with CAMB, emcee, dynesty, corner dependencies
- `.gitignore` for clean repository management
- MIT License file
- GitHub Actions CI workflow for Python testing

### Changed
- Updated all GitHub repository references to `FaCC`
- Enhanced README with new framework documentation
- Improved code organization with `src/` directory structure

### Fixed
- Repository independence clearly documented
- All URL references updated to correct repository name

