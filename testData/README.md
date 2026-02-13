# Test Data Directory

This directory contains sample data files for testing pyIrena.

## Files

Place your test HDF5 or text data files here.

### Expected Format

#### HDF5 Files (NXcanSAS)

Files should follow the NXcanSAS standard with:
- Q values (scattering vector in 1/Å)
- Intensity values (in cm⁻¹)
- Optional: Error values, dQ values

#### Text Files

Tab or space-separated columns:
```
Q           Intensity   Error
0.001       1000.0      10.0
0.002       950.0       9.5
...
```

## Creating Test Data

You can create synthetic test data with Python:

```python
import numpy as np
import h5py

# Create synthetic SAXS data
q = np.logspace(-3, 0, 100)
intensity = 1000 * q**(-4) + 0.01  # Porod + background
error = 0.05 * intensity

# Save as HDF5 (simplified, not full NXcanSAS)
with h5py.File('testData/synthetic_data.h5', 'w') as f:
    f.create_dataset('Q', data=q)
    f.create_dataset('Intensity', data=intensity)
    f.create_dataset('Error', data=error)
```

## Sample Data

Copy your experimental data files here for testing the GUI.
