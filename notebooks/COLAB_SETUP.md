# Colab Setup Instructions

## Quick Setup for Colab

### Step 1: Clone Repository

```python
!git clone https://github.com/imediacorp/FaCC.git
%cd FaCC
```

### Step 2: Install Dependencies

```python
!pip install camb numpy scipy matplotlib pandas astropy
```

### Step 3: Run Notebook

The notebook `notebooks/01_desi_forecasts.ipynb` has been updated with robust path handling that works in Colab. Just run the cells in order.

### Alternative: Direct Import in Colab

If you want to import directly in a Colab notebook, use this cell:

```python
import sys
import os

# Clone repository (if not already done)
if not os.path.exists('FaCC'):
    !git clone https://github.com/imediacorp/FaCC.git

# Change to repository directory
os.chdir('FaCC')

# Add src to path
sys.path.insert(0, 'src')

# Install dependencies
!pip install -q camb numpy scipy matplotlib pandas astropy

# Now import
from phi_modulation import PhiModulationModel
from systematics import SystematicErrorBudget  # if needed

print("✅ Setup complete!")
```

## Verifying the Setup

Run this to verify everything works:

```python
# Check imports
from phi_modulation import PhiModulationModel

# Initialize model
model = PhiModulationModel()

# Test basic functionality
k, z, Pk = model.get_base_power_spectrum(k_min=0.01, k_max=0.1, npoints=50)
print(f"✅ Success! Generated power spectrum with {len(k)} k-bins")
print(f"   Golden ratio φ = {model.phi:.8f}")
```

## Troubleshooting

### ModuleNotFoundError

If you get `ModuleNotFoundError: No module named 'phi_modulation'`:

1. Make sure you've cloned the repository: `!git clone https://github.com/imediacorp/FaCC.git`
2. Change to the directory: `%cd FaCC`
3. Check that src directory exists: `!ls -la src/`
4. Add src to path: `sys.path.insert(0, 'src')`

### CAMB Errors

If CAMB installation fails:

```python
# Try installing with system dependencies
!apt-get install -y gfortran
!pip install camb
```

### Path Issues

If path issues persist, use absolute paths:

```python
import sys
import os

# Get absolute path to src
repo_path = '/content/FaCC'  # Colab default after clone
src_path = os.path.join(repo_path, 'src')

if src_path not in sys.path:
    sys.path.insert(0, src_path)

from phi_modulation import PhiModulationModel
```

