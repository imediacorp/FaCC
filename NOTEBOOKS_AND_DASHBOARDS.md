# Notebooks and Interactive Tools

## Colab Notebooks

### Main Analysis Notebook
**Location:** `fibonacci_cosmology_analysis.ipynb` (root directory)

**Colab Link:**
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/imediacorp/FaCC/blob/main/fibonacci_cosmology_analysis.ipynb)

This notebook contains the legacy analysis including:
- H(z) cosmic chronometer fits
- CMB low-ℓ residual analysis
- Matter power spectrum P(k) analysis
- Legacy φ-modulation tests

### New DESI Forecasts Notebook
**Location:** `notebooks/01_desi_forecasts.ipynb`

**To run in Colab:**
1. Upload the notebook to Colab
2. Install dependencies:
   ```python
   !pip install camb numpy scipy matplotlib
   ```
3. Clone or upload the repository files (or run cells that import from GitHub)

This notebook includes:
- **New φ-modulation framework** using `PhiModulationModel`
- DESI Year 5 forecast analysis
- 4-panel forecast visualization
- Systematic error analysis (when systematics module is available)
- Fisher matrix sensitivity calculations

### Kaggle Notebook
**Location:** `fibonacci_cosmology_kaggle.ipynb`

**Kaggle Link:** https://www.kaggle.com/code/bryanpersaud/faac-notebook/

Optimized for Kaggle environment with robust data loading.

## Streamlit Dashboard

**Status:** ❌ Not currently implemented

There is no Streamlit dashboard in the repository at this time. However, creating one would be straightforward given the modular code structure.

### Potential Dashboard Features

If you'd like a Streamlit dashboard, it could include:

1. **Interactive Forecast Tool**
   - Sliders for A_φ, redshift, k-range
   - Real-time power spectrum visualization
   - Forecast sensitivity plots
   - SNR calculations

2. **Systematic Error Explorer**
   - Toggle different systematic error sources
   - Visualize error breakdown
   - Compare statistical vs. systematic uncertainties

3. **Parameter Space Explorer**
   - Interactive cosmological parameter exploration
   - φ-modulation parameter sweeps
   - Real-time model comparison

4. **Data Visualization**
   - Upload and visualize power spectrum data
   - Compare with φ-modulated models
   - Residual analysis

### Creating a Streamlit Dashboard

If you want to create one, I can help set it up. It would need:

```python
# Example structure
import streamlit as st
from src.phi_modulation import PhiModulationModel
from src.systematics import SystematicErrorBudget
import matplotlib.pyplot as plt
import numpy as np

st.title("φ-Modulation Analysis Dashboard")
# ... interactive widgets and plots
```

Would require adding `streamlit` to `requirements.txt`.

## Quick Access

### For Colab Users

**Main notebook:**
```python
# In Colab, you can clone the repo:
!git clone https://github.com/imediacorp/FaCC.git
cd FaCC
!pip install camb numpy scipy matplotlib pandas astropy
```

**Or use the direct link:**
- Click the Colab badge in README.md
- Or visit: https://colab.research.google.com/github/imediacorp/FaCC/blob/main/fibonacci_cosmology_analysis.ipynb

### For Local Jupyter

```bash
# Install dependencies
pip install -r requirements.txt

# Run notebook
jupyter notebook fibonacci_cosmology_analysis.ipynb
# or
jupyter notebook notebooks/01_desi_forecasts.ipynb
```

## Recommendations

1. **For new analysis:** Use `notebooks/01_desi_forecasts.ipynb` - it uses the latest framework
2. **For legacy analysis:** Use `fibonacci_cosmology_analysis.ipynb`
3. **For interactive exploration:** Consider creating a Streamlit dashboard (I can help!)

Would you like me to create a Streamlit dashboard for the φ-modulation analysis?

