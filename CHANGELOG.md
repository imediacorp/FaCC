# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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

