# Fix for Colab - CAMB API Issue

## Problem
You're seeing this error in Colab:
```
AttributeError: 'tuple' object has no attribute 'k'
```

This happens because the Colab repository has an older version of the code that uses the wrong CAMB API.

## Quick Fix in Colab

Run this cell in your Colab notebook to fix the code:

```python
# Fix the CAMB API issue in phi_modulation.py
import os

phi_mod_path = 'src/phi_modulation.py'

# Read the file
with open(phi_mod_path, 'r') as f:
    content = f.read()

# Check if it needs fixing
if 'return powers.k, powers.z, powers.P_k' in content:
    print("Fixing CAMB API issue...")
    
    # Replace the old code with the fixed version
    old_code = """        powers = results.get_matter_power_spectrum(
            minkh=k_min, maxkh=k_max, npoints=npoints
        )
        
        return powers.k, powers.z, powers.P_k"""
    
    new_code = """        kh, z_arr, pk = results.get_matter_power_spectrum(
            minkh=k_min, maxkh=k_max, npoints=npoints
        )
        
        return kh, z_arr, pk"""
    
    content = content.replace(old_code, new_code)
    
    # Write back
    with open(phi_mod_path, 'w') as f:
        f.write(content)
    
    print("✅ Fixed! Please restart the runtime (Runtime > Restart runtime)")
    print("   Then re-run your import cells.")
else:
    print("✅ Code is already fixed!")
```

## Alternative: Manual Edit

Or manually edit `src/phi_modulation.py` in Colab:

1. Find the `get_base_power_spectrum` method (around line 100)
2. Change this:
   ```python
   powers = results.get_matter_power_spectrum(...)
   return powers.k, powers.z, powers.P_k
   ```
   
3. To this:
   ```python
   kh, z_arr, pk = results.get_matter_power_spectrum(...)
   return kh, z_arr, pk
   ```

4. Save the file
5. Restart runtime: Runtime > Restart runtime
6. Re-run your cells

## Best Solution: Pull Latest Code

The fix is already in the GitHub repository. Pull the latest changes:

```python
# Pull latest code from GitHub
!cd FaCC && git pull origin main

# Restart runtime
# Runtime > Restart runtime

# Then re-run your cells
```

## Verify the Fix

After fixing, verify it works:

```python
from phi_modulation import PhiModulationModel

model = PhiModulationModel()
k, z, pk = model.get_base_power_spectrum(k_min=0.01, k_max=0.1, npoints=10)
print(f"✅ Success! k shape: {k.shape}, pk shape: {pk.shape}")
```

