# Release v0.1.0 - Initial Release

**Release Date:** January 2025

## Overview

This is the initial release of the FaCC (Fibonacci Cosmology) repository, introducing a scientifically defensible framework for testing φ-modulation as an empirical pattern in cosmic structure.

## Key Features

### Core Analysis Framework

- **PhiModulationModel Class**: A comprehensive Python class implementing φ-modulated power spectrum analysis within the ΛCDM framework
  - CAMB integration for generating physically accurate ΛCDM power spectra
  - Log-periodic modulation: `P(k) = P_ΛCDM(k) × [1 + A_φ × cos(2π × log(k/k_pivot) / ln(φ) + φ_0)]`
  - BAO signature computation with φ-modulation
  - DESI Year 5 sensitivity forecasts using Fisher matrix analysis

### Analysis Tools

- **DESI Forecast Notebook** (`notebooks/01_desi_forecasts.ipynb`): Comprehensive analysis notebook featuring:
  - Power spectrum ratio analysis
  - BAO signature comparison
  - DESI SNR vs amplitude forecasts
  - 4-panel publication-quality figures

### Documentation

- Complete README with usage examples
- Independence statement (`INDEPENDENCE.md`) clarifying separation from other projects
- Comprehensive code documentation
- Release notes and changelog

### Scientific Approach

This release implements a **two-parameter extension to ΛCDM** (amplitude A_φ and phase φ_0), representing a shift from hypothesis-proving to hypothesis-testing with proper statistical rigor.

**Key Scientific Claims:**
- Tests whether the Golden Ratio (φ ≈ 1.618) appears in cosmic structure as an empirical pattern
- Forecasts suggest DESI Year 5 can detect oscillations with amplitude A_φ ≳ 0.005 at 3σ confidence
- Model is testable, falsifiable, and respects ΛCDM as the baseline

## Installation

```bash
pip install -r requirements.txt
```

**Key Dependencies:**
- CAMB (for ΛCDM power spectrum calculations)
- NumPy, SciPy, Matplotlib, Pandas
- Astropy
- Optional: emcee, dynesty, corner (for Bayesian analysis)

## Quick Start

```python
from src.phi_modulation import PhiModulationModel

# Initialize model
model = PhiModulationModel()

# Generate power spectrum
k, z, Pk = model.get_base_power_spectrum(k_min=0.01, k_max=0.3)

# Apply φ-modulation
Pk_mod, mod_factor = model.apply_phi_modulation(k, Pk[0], A_phi=0.01)

# Forecast DESI sensitivity
forecast = model.forecast_desi_sensitivity(A_phi_true=0.01)
print(f"Forecast σ_Aφ = {forecast['sigma_Aphi']:.4f}, SNR = {forecast['SNR']:.2f}σ")
```

## Repository

- **GitHub:** https://github.com/imediacorp/FaCC
- **License:** MIT
- **Documentation:** See README.md

## Citation

If you use this code in your research, please cite:

```bibtex
@misc{persaud_fibonacci_cosmology,
  author       = {Bryan David Persaud},
  title        = {Fibonacci Cosmology: Falsified Background, Testable Perturbations},
  year         = {2025},
  howpublished = {GitHub repository},
  url          = {https://github.com/imediacorp/FaCC}
}
```

## Acknowledgments

This work represents independent cosmological research exploring whether the Golden Ratio manifests in cosmic structure formation. It originated as a thought experiment and hypothesis motivated by self-similarity patterns observed from plants to galaxies.

