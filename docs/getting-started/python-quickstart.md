# Quick Start Guide

This guide will help you get started with pyIrena in just a few minutes.

## Installation

```bash
pip install pyirena
```

Or install with plotting support:

```bash
pip install pyirena[plotting]
```

## Your First Fit

Here's a complete example of fitting small-angle scattering data:

```python
import numpy as np
from pyirena.core.unified import UnifiedFitModel

# 1. Load your data (or create synthetic data for this example)
q = np.logspace(-3, 0, 100)  # Q values in 1/Å
intensity = np.ones_like(q)  # Replace with your measured intensity
error = 0.1 * intensity      # Replace with your error values

# 2. Create a model with one structural level
model = UnifiedFitModel(num_levels=1)

# 3. Set initial parameter guesses
model.levels[0].Rg = 50.0    # Radius of gyration in Ångströms
model.levels[0].G = 1000.0   # Guinier prefactor
model.levels[0].P = 4.0      # Power law slope (4 = Porod)
model.levels[0].B = 1e-3     # Power law prefactor
model.background = 0.01      # Flat background

# 4. Configure what to fit
model.levels[0].fit_Rg = True
model.levels[0].fit_G = True
model.levels[0].fit_P = False  # Keep P fixed at 4.0
model.levels[0].fit_B = True
model.fit_background = True

# 5. Run the fit
results = model.fit(q, intensity, error, verbose=1)

# 6. View results
print(model.get_parameter_summary())
print(f"\nFit successful: {results['success']}")
print(f"Reduced χ²: {results['reduced_chi_squared']:.4f}")
```

## Loading Real Data

### From NXcanSAS HDF5 files

```python
from pyirena.core.unified import load_data_from_nxcansas

# Load data
data = load_data_from_nxcansas('path/to/your/data.h5')

# Extract arrays
q = data['Q']
intensity = data['Intensity']
error = data['Error']

# Now fit as above
model = UnifiedFitModel(num_levels=1)
results = model.fit(q, intensity, error)
```

### From text files

```python
import numpy as np

# Load three-column data: Q, I, Error
data = np.loadtxt('data.txt', skiprows=1)  # Skip header row
q = data[:, 0]
intensity = data[:, 1]
error = data[:, 2]
```

## Plotting Results

If you have matplotlib installed:

```python
from pyirena.plotting.unified_plots import plot_fit_results

# Plot the fit
plot_fit_results(
    model,
    show_residuals=True,
    log_scale=True,
    save_path='fit_result.png'
)
```

Or make your own plot:

```python
import matplotlib.pyplot as plt

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)

# Plot data and fit
ax1.errorbar(model.q_data, model.I_data, model.error_data,
             fmt='o', label='Data', markersize=4)
ax1.plot(model.q_data, model.fit_intensity, 'r-', label='Fit', linewidth=2)
ax1.set_ylabel('Intensity (cm⁻¹)')
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.legend()
ax1.grid(True, alpha=0.3)

# Plot residuals
residuals = (model.I_data - model.fit_intensity) / model.error_data
ax2.plot(model.q_data, residuals, 'o', markersize=4)
ax2.axhline(0, color='k', linestyle='--', linewidth=1)
ax2.set_xlabel('Q (Å⁻¹)')
ax2.set_ylabel('Residuals (σ)')
ax2.set_xscale('log')
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('my_fit.png', dpi=150)
plt.show()
```

## Multi-Level Fitting

For hierarchical structures (e.g., primary particles + aggregates):

```python
# Create 2-level model
model = UnifiedFitModel(num_levels=2)

# Level 1: Small primary particles
model.levels[0].Rg = 20.0       # Small Rg
model.levels[0].G = 100.0
model.levels[0].P = 4.0         # Sharp interface
model.levels[0].B = 1e-3

# Level 2: Large aggregates
model.levels[1].Rg = 200.0      # Large Rg
model.levels[1].G = 5000.0
model.levels[1].P = 2.5         # Mass fractal
model.levels[1].link_RGCO = True  # Link to level 1 Rg

# Fit
results = model.fit(q, intensity, error)
```

## Common Use Cases

### 1. Spherical Particles (Porod Scattering)

```python
model = UnifiedFitModel(num_levels=1)
model.levels[0].P = 4.0         # Porod law
model.levels[0].fit_P = False   # Keep fixed
```

### 2. Mass Fractal Aggregates

```python
model = UnifiedFitModel(num_levels=1)
model.levels[0].mass_fractal = True
model.levels[0].P = 2.5         # Fractal dimension
# B is calculated automatically
```

### 3. With Correlation Function

```python
model = UnifiedFitModel(num_levels=1)
model.levels[0].correlations = True
model.levels[0].ETA = 60.0      # Correlation distance
model.levels[0].PACK = 0.3      # Packing factor
model.levels[0].fit_ETA = True
model.levels[0].fit_PACK = True
```

## Parameter Guidelines

| Parameter | Typical Range | Physical Meaning |
|-----------|--------------|------------------|
| Rg | 1-10,000 Å | Size of structure |
| G | 0.01-10⁶ cm⁻¹ | Contrast × volume fraction |
| P | 0-6 | Surface characteristics |
| B | 10⁻²⁰-1 cm⁻¹Å⁻ᴾ | Surface area related |

### Power Law Slope Interpretation

- **P ≈ 4**: Sharp interface (Porod law)
- **P ≈ 3**: Surface fractal
- **P = 2-3**: Mass fractal
- **P > 4**: Diffuse interface

## Tips for Better Fits

1. **Start simple**: Begin with 1 level, add more if needed
2. **Fix parameters**: Don't fit everything at once
3. **Check bounds**: Ensure parameter limits make physical sense
4. **Use good initial guesses**: Look at your data first
5. **Examine residuals**: Random scatter = good fit

## Next Steps

- Read the full [User Guide](../USAGE_GUIDE.md)
- Explore [Examples](../pyirena/examples/)
- Check [API Documentation](api.md)
- Learn the [theory](theory.md)

## Getting Help

- Check [GitHub Issues](https://github.com/jilavsky/pyirena/issues)
- Read the [FAQ](faq.md)
- Contact: ilavsky@aps.anl.gov
