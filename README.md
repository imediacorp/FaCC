# Fibonacci Cosmology: Falsified Background, Testable Perturbations

[![Python CI](https://github.com/imediacorp/FaCC/actions/workflows/python-ci.yml/badge.svg)](https://github.com/imediacorp/FaCC/actions/workflows/python-ci.yml)
[![Python 3.9–3.12](https://img.shields.io/badge/python-3.9--3.12-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![arXiv](https://img.shields.io/badge/arXiv-preprint-b31b1b.svg)](https://arxiv.org/search/?searchtype=author&query=Persaud)

Author: Bryan David Persaud  
Affiliation: Intermedia Communications Corp.  
Contact: bryan@imediacorp.com

---

## Overview
This repository accompanies two complementary works investigating whether the Golden Ratio (φ ≈ 1.618) manifests in cosmology:

1. **Background model (falsified)**: tests φ-recursion in the expansion history. Result: decisively ruled out by H(z), CMB low-ℓ, and P(k) data.
2. **Perturbation model (falsifiable)**: proposes weak log-periodic modulations in structure growth, yielding φ-spaced features in the matter power spectrum and CMB. Provides concrete forecasts for near-future surveys.

The code here reproduces the core figures, basic fits, and forecast visuals, and lets you explore signals using public summary datasets.

**Note on Independence**: This cosmological research work is completely independent and separate from any other projects. It originated as a simple thought experiment and hypothesis that the Golden Ratio (Fibonacci sequence) might be fundamental to cosmic structure, given the self-similarity patterns observed from plants to galaxies.

### Theorem overview (short)
- **Golden Ratio as recursive constant**: The hypothesis explores φ as a fundamental recursion constant. In a background variant (now falsified), φ enters the expansion rate via recursive scaling; in the perturbation variant, φ induces log-periodic modulations in clustering statistics.
- **Duality for perturbations**: Two coupled modes appear naturally — a forward (φ) and a conjugate (φ−1) branch — leading to a dual log-periodic signature. In practice, this means residuals may contain cosines periodic in log(k) or log(ℓ) with periods set by ln(φ) and ln(φ−1).
- **Falsifiable prediction**: If these dual log-periodic features are absent at the sensitivity of current/forthcoming surveys, the perturbation hypothesis is ruled out.

---

## System Requirements

### Python
Python 3.9–3.12 is required.

### Fortran compiler (required for CAMB)
`camb` compiles Fortran code on installation. You must have a Fortran compiler installed **before** running `pip install -r requirements.txt`:

**macOS:**
```bash
brew install gcc   # provides gfortran
```

**Ubuntu/Debian:**
```bash
sudo apt-get install gfortran libgfortran5
```

**Conda (all platforms):**
```bash
conda install -c conda-forge gfortran
```

If `pip install camb` fails with a compiler error, a missing Fortran compiler is almost always the cause.

---

## Installation

```bash
# 1. Clone the repository
git clone https://github.com/imediacorp/FaCC.git
cd FaCC

# 2. Install Fortran compiler (see above — required for CAMB)

# 3. Create and activate a virtual environment
python -m venv .venv
source .venv/bin/activate          # macOS/Linux
# .venv\Scripts\activate           # Windows

# 4. Install Python dependencies
pip install -r requirements.txt

# 5. (Optional) Install development tools for testing and linting
pip install -r requirements-dev.txt
```

**Conda users:**
```bash
conda create -n facc python=3.11
conda activate facc
conda install -c conda-forge gfortran
pip install -r requirements.txt
```

---

## Quick Start

```python
from src.phi_modulation import PhiModulationModel

# Initialize model with Planck 2018 parameters
model = PhiModulationModel()

# Generate base ΛCDM power spectrum
k, z, Pk = model.get_base_power_spectrum(k_min=0.01, k_max=0.3)

# Apply φ-modulation: P_mod = P_ΛCDM × [1 + A_φ × cos(2π × log(k/k_pivot) / ln(φ) + φ_0)]
Pk_mod, mod_factor = model.apply_phi_modulation(k, Pk[0], A_phi=0.01)

# Forecast DESI sensitivity using Fisher matrix
result = model.forecast_desi_sensitivity(A_phi_true=0.01)
print(f"Forecast σ_Aφ = {result['sigma_Aphi']:.4f}, SNR = {result['SNR']:.2f}σ")
```

**Run the forecast notebook:**
```bash
jupyter notebook notebooks/01_desi_forecasts.ipynb
```

---

## Papers
- `fibonacci_cosmology_falsified.pdf` — Background model and its observational falsification
- `fibonacci_perturbations.pdf` — Perturbation-level hypothesis and detectability forecasts

If you use this repository, please cite the corresponding paper(s).

### Paper summaries
- **Background (falsified)**: A φ-recursive background expansion was tested against observations. H(z) fits, low-ℓ CMB residuals, and matter P(k) collectively rule it out with large Δχ². Key takeaway: φ is not a viable background constant for cosmic acceleration.
- **Perturbations (falsifiable)**: Keeping ΛCDM background, introduce dual log-periodic modulations governed by φ and φ−1 in the growth of structure. Predicts φ-spaced features in P(k) around the BAO scale and small CMB residual oscillations. Forecasts suggest near-future surveys could detect or refute these signals decisively.

---

## Repository Structure

```
.
├── src/                            # Core analysis modules
│   ├── __init__.py
│   ├── phi_modulation.py           # PhiModulationModel class (CAMB-based)
│   ├── bayesian_tools.py           # BayesianEvidence, BIC utilities
│   ├── systematics.py              # SystematicErrorBudget class
│   ├── hz_analysis.py              # H(z) fitting functions (importable)
│   ├── cmb_analysis.py             # CMB oscillation detection (importable)
│   └── lss_analysis.py             # LSS P(k) analysis (importable)
├── tests/                          # pytest test suite
│   ├── __init__.py
│   ├── test_phi_modulation.py
│   ├── test_bayesian_tools.py
│   └── test_systematics.py
├── notebooks/                      # Analysis notebooks
│   └── 01_desi_forecasts.ipynb
├── examples/
│   └── systematics_example.py
├── docs/                           # Documentation and notes
├── .github/workflows/python-ci.yml
├── fibonacci_cosmology_falsified.pdf
├── fibonacci_perturbations.pdf
├── real_hz.csv
├── real_cmb_lowl.csv
├── real_pk_lowk.csv
├── hz_fitter_v2.py                 # Legacy H(z) script (logic also in src/hz_analysis.py)
├── cmb_osc_detector_v2.py          # Legacy CMB script (logic also in src/cmb_analysis.py)
├── lss_phi_analyzer.py             # Legacy P(k) script (logic also in src/lss_analysis.py)
├── requirements.txt
└── requirements-dev.txt
```

---

## Usage Guide

### Core Framework

**Run the DESI forecast analysis notebook:**
```bash
jupyter notebook notebooks/01_desi_forecasts.ipynb
```

**Use `PhiModulationModel` in your own code:**
```python
from src.phi_modulation import PhiModulationModel

model = PhiModulationModel()
# Or with custom cosmological parameters:
model = PhiModulationModel(params={
    'H0': 67.36, 'ombh2': 0.02237, 'omch2': 0.1200,
    'As': 2.1e-9, 'ns': 0.9649, 'tau': 0.0544
})

k, z, Pk = model.get_base_power_spectrum(k_min=0.01, k_max=0.3, npoints=500)
Pk_mod, mod_factor = model.apply_phi_modulation(k, Pk[0], A_phi=0.01, phi_phase=0.0, k_pivot=0.05)
r, xi_base, xi_mod = model.compute_bao_signature(z=0.5, A_phi=0.01)

forecast = model.forecast_desi_sensitivity(A_phi_true=0.01)
print(f"σ_Aφ = {forecast['sigma_Aphi']:.4f}, SNR = {forecast['SNR']:.2f}σ")
```

### Legacy Analysis Scripts

These scripts run standalone from the repo root. Data CSVs must be present in the project root.

```bash
python hz_fitter_v2.py          # H(z) falsification → hz_comparison_dual.png
python cmb_osc_detector_v2.py   # CMB log-periodic residuals → cmb_osc_dual.png
python lss_phi_analyzer.py      # P(k) φ-scale analysis → pk_phi.png
python forecast_plots.py        # Forecast visualization → pk_forecast.png
python fetch_desi_bao_pk.py     # Reconstruct DESI P(k) → pk_data_real.csv
```

The same logic is also available as importable functions in `src/hz_analysis.py`, `src/cmb_analysis.py`, and `src/lss_analysis.py`.

---

## Running Tests

```bash
pip install -r requirements-dev.txt
pytest tests/ -v --cov=src --cov-report=term-missing
```

Tests run without CAMB, covering all pure-Python analysis logic.

---

## Data Descriptions
- `real_hz.csv`: Cosmic-chronometer H(z) catalog (columns: z, H, sigma_H).
- `real_cmb_lowl.csv`: Planck TT low-ℓ summary (columns: ell, C_ell, sigma).
- `real_pk_lowk.csv`: SDSS/BOSS P(k) for low-k scales (columns: k, P(k), sigma_Pk).
- `pk_data.csv`: Simple synthetic P(k) used by examples.
- `pk_data_real.csv`: Reconstructed P(k) produced by `fetch_desi_bao_pk.py` (DESI DR2-like).

---

## Data Sources & Attribution
- H(z) cosmic chronometers: aggregated public compilations (see file header/commit history).
- Planck low-ℓ TT: simplified summary exported to `real_cmb_lowl.csv`; refer to Planck Collaboration papers for authoritative spectra.
- Low-k P(k): `real_pk_lowk.csv` is a reduced sample for demonstrating analysis workflows.
- DESI DR2 BAO summary: fetched from `data.desi.lbl.gov` when internet is enabled; otherwise `fetch_desi_bao_pk.py` uses a small offline fallback.

Please credit the original survey teams (Planck Collaboration, SDSS/BOSS/eBOSS, DESI Collaboration).

---

## Using on Kaggle

- Kaggle profile: https://www.kaggle.com/bryanpersaud
- FaCC Notebook: https://www.kaggle.com/code/bryanpersaud/faac-notebook/

```python
!pip -q install numpy scipy matplotlib pandas astropy sympy pillow requests beautifulsoup4
```

Upload `real_hz.csv`, `real_cmb_lowl.csv`, and `real_pk_lowk.csv` as Kaggle Dataset files if internet is off.

## Run in Google Colab

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/imediacorp/FaCC/blob/main/fibonacci_cosmology_analysis.ipynb)

```python
!pip install -q numpy scipy matplotlib pandas astropy camb emcee corner
!git clone https://github.com/imediacorp/FaCC.git
%cd FaCC
```

---

## Reproducibility
- Scripts are deterministic given the fixed input CSVs; random seeds are not used.
- Baselines: Some routines approximate ΛCDM trends (e.g., polynomial baseline for low-ℓ CMB) to isolate residuals; these are analysis conveniences, not full Boltzmann codes.
- Units: P(k) is plotted in customary units; axes are labeled in each figure.

---

## License
MIT License — see [LICENSE](LICENSE).

---

## Citation

```bibtex
@misc{persaud_fibonacci_cosmology,
  author       = {Bryan David Persaud},
  title        = {Fibonacci Cosmology: Falsified Background, Testable Perturbations},
  year         = {2025},
  howpublished = {GitHub repository},
  url          = {https://github.com/imediacorp/FaCC}
}
```

---

## Project Independence

This cosmological research work is completely independent and separate from any other projects. See [INDEPENDENCE.md](INDEPENDENCE.md) for a detailed statement.

---

## Troubleshooting

**`ModuleNotFoundError: No module named 'camb'`**  
CAMB requires a Fortran compiler. Install gfortran first (see System Requirements above), then `pip install camb`. The dashboard H(z), CMB, LSS, DESI BAO, and G–φ tabs all work without CAMB.

**`ImportError: mach-o file, but is an incompatible architecture`** (Apple Silicon)  
You have an x86_64 package installed in an arm64 Python environment (or vice versa). Create a fresh virtual environment with the native arm64 Python:
```bash
/usr/bin/python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

**`ssl.SSLCertVerificationError` on macOS Python 3.13**  
Run the certificate installer bundled with Python:
```bash
open /Applications/Python\ 3.13/Install\ Certificates.command
```

**`scipy.optimize` fails to converge in CMB fit**  
The dual log-periodic fit can be sensitive to initial conditions. If `curve_fit` raises `RuntimeError`, your CMB data range may be too narrow. Ensure `real_cmb_lowl.csv` covers ℓ = 2–30.

**DESI BAO tab shows a network error**  
The tab fetches live data from GitHub. Check your internet connection. Data is cached for 1 hour per session — refreshing the page will retry.

**Tests fail with `ModuleNotFoundError: No module named 'src'`**  
Run pytest from the repo root:
```bash
cd /path/to/FaCC && pytest tests/ -v
```

---

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for the PR process, code style guide, and testing requirements.

---

## Support

Questions or issues: please open a [GitHub issue](https://github.com/imediacorp/FaCC/issues) or email bryan@imediacorp.com.
