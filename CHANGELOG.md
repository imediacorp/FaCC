# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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

