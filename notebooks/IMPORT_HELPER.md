# Import Helper for Notebooks

## Problem
When running notebooks (especially in Colab), imports may fail with `ModuleNotFoundError: No module named 'phi_modulation'`.

## Solutions

### Option 1: Use Absolute Path (Recommended for Colab)

```python
import sys
import os

# Add src directory to path (works from notebook directory)
src_path = os.path.abspath(os.path.join(os.getcwd(), '..', 'src'))
if src_path not in sys.path:
    sys.path.insert(0, src_path)

from phi_modulation import PhiModulationModel
```

### Option 2: Clone Repository in Colab

```python
# In Colab, clone the repository first
!git clone https://github.com/imediacorp/FaCC.git
%cd FaCC

# Then import
import sys
sys.path.insert(0, 'src')

from phi_modulation import PhiModulationModel
```

### Option 3: Use Full Module Path

```python
import sys
import os

# Get the project root (two levels up from notebooks/)
project_root = os.path.abspath(os.path.join(os.getcwd(), '..'))
src_path = os.path.join(project_root, 'src')

if src_path not in sys.path:
    sys.path.insert(0, src_path)

from phi_modulation import PhiModulationModel
```

### Option 4: Install as Package (For Local Development)

```python
# In notebook directory or project root
import sys
import os

# Add project root to path
project_root = os.path.abspath(os.path.join(os.getcwd(), '..'))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

# Now import from src
from src.phi_modulation import PhiModulationModel
```

## Quick Fix for Current Notebook

Add this cell at the beginning of your notebook:

```python
import sys
import os

# Determine the project root
if 'notebooks' in os.getcwd():
    # We're in notebooks/ directory
    project_root = os.path.abspath(os.path.join(os.getcwd(), '..'))
else:
    # We're in project root
    project_root = os.getcwd()

# Add src to path
src_path = os.path.join(project_root, 'src')
if src_path not in sys.path:
    sys.path.insert(0, src_path)

# Verify
print(f"Project root: {project_root}")
print(f"Added to path: {src_path}")
print(f"Files in src: {os.listdir(src_path) if os.path.exists(src_path) else 'NOT FOUND'}")

# Now import
from phi_modulation import PhiModulationModel
from systematics import SystematicErrorBudget  # if needed
```

## For Colab Specifically

```python
# Step 1: Clone repository
!git clone https://github.com/imediacorp/FaCC.git
%cd FaCC

# Step 2: Install dependencies
!pip install camb numpy scipy matplotlib pandas astropy

# Step 3: Set up path and import
import sys
sys.path.insert(0, 'src')

from phi_modulation import PhiModulationModel
print("✅ Import successful!")
```

## Verification

After setting up the path, verify with:

```python
import phi_modulation
print(f"✅ phi_modulation module found at: {phi_modulation.__file__}")

model = PhiModulationModel()
print(f"✅ Model initialized: φ = {model.phi:.8f}")
```

