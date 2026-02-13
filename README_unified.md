# Unified Fit Model for Small-Angle Scattering

A comprehensive Python implementation of the Unified fit model (Beaucage method) for analyzing hierarchical structures in small-angle scattering (SAXS/SANS/USAXS) data.

## Overview

The Unified fit model combines multiple structural levels using Guinier-Porod crossover functions to describe complex hierarchical materials. Each level represents a distinct structural feature (e.g., primary particles, aggregates, etc.).

### Mathematical Model

For each structural level *i*, the intensity is given by:

```
I_i(q) = G_i × exp(-q²R_g,i²/3) + exp(-q²R_g(i-1)²/3) × B_i × {[erf(k×q×R_g,i/√6)]³/q}^P_i
```

Where:
- **G_i**: Guinier prefactor (low-q amplitude)
- **R_g,i**: Radius of gyration
- **B_i**: Power law prefactor (Porod constant)
- **P_i**: Power law slope
- **k**: Correction factor (1.0 for P > 3, 1.06 for P ≤ 3)

Optional correlation function (Born-Green approximation):
```
I_i(q) = I_i(q) / (1 + PACK × f(q, ETA))
```

### References

- Beaucage, G. (1995). *J. Appl. Cryst.* **28**, 717-728
- Beaucage, G. (1996). *J. Appl. Cryst.* **29**, 134-146
- Beaucage, G. et al. (1997). *Macromolecules* **30**, 4158-4166

## Installation

### Requirements

```bash
pip install -r requirements.txt
```

Required packages:
- numpy >= 1.20.0
- scipy >= 1.7.0
- h5py >= 3.0.0
- matplotlib >= 3.3.0 (optional, for plotting)

### Files

- **unified.py**: Core model implementation
- **unified_utils.py**: Plotting and analysis utilities
- **unified_demo.py**: Demonstration scripts with examples
- **hdf5code.py**: HDF5/NXcanSAS data loading functions
- **requirements.txt**: Package dependencies

## Quick Start

### Basic Example

```python
import numpy as np
from unified import UnifiedFitModel

# Create model with 1 structural level
model = UnifiedFitModel(num_levels=1)

# Set initial parameters
model.levels[0].Rg = 50.0    # Radius of gyration [Å]
model.levels[0].G = 1000.0   # Guinier prefactor [cm^-1]
model.levels[0].P = 4.0      # Power law slope
model.levels[0].B = 1e-3     # Porod constant
model.background = 0.01      # Flat background

# Fit to data
results = model.fit(q_data, intensity_data, error_data)

# Display results
print(model.get_parameter_summary())
```

### Loading Data from NXcanSAS

```python
from unified import load_data_from_nxcansas

# Load data from HDF5 file
data = load_data_from_nxcansas('/path/to/file.h5')

# Access arrays
q = data['Q']
intensity = data['Intensity']
error = data['Error']
dq = data['dQ']
```

### Multi-Level Fitting

```python
# Create model with 2 levels
model = UnifiedFitModel(num_levels=2)

# Level 1: Primary particles
model.levels[0].Rg = 20.0
model.levels[0].G = 100.0
model.levels[0].P = 4.0
model.levels[0].B = 1e-3

# Level 2: Aggregates
model.levels[1].Rg = 200.0
model.levels[1].G = 10000.0
model.levels[1].P = 2.5  # Mass fractal
model.levels[1].link_RGCO = True  # Link cutoff to level 1 Rg

# Perform fit
results = model.fit(q, intensity, error)
```

### Plotting Results

```python
from unified_utils import plot_fit_results, plot_guinier_analysis, plot_porod_analysis

# Plot fit with residuals
plot_fit_results(model, show_residuals=True, log_scale=True)

# Guinier plot
plot_guinier_analysis(model, level_idx=0)

# Porod plot
plot_porod_analysis(model, level_idx=0)
```

## Features

### Model Capabilities

- **Multiple structural levels** (1-5 levels)
- **Hierarchical linking** via RgCO parameter
- **Correlation functions** (Born-Green sphere amplitude)
- **Mass fractal mode** with automatic B calculation
- **Automatic B estimation** from Hammouda relationship
- **Slit smearing** support for USAXS data
- **Parameter bounds** and constraints
- **Robust optimization** using scipy.optimize.least_squares

### Analysis Tools

- **Invariant calculation**: Q = ∫ I(q) q² dq
- **Surface/volume ratio**: For Porod scattering (P ≈ 4)
- **Size distribution moments**: Average radius, volume, surface
- **Parameter estimation**: Automatic initial guess from data
- **Export capabilities**: Text files with fit results and data tables

### Plotting Functions

- `plot_fit_results()`: Data, model, levels, and residuals
- `plot_guinier_analysis()`: ln(I) vs Q² for Guinier region
- `plot_porod_analysis()`: I×Q⁴ vs Q for Porod region
- All plots support log/linear scales and export to files

## API Reference

### UnifiedLevel

Dataclass representing one structural level.

**Parameters:**
- `Rg`: Radius of gyration [Å]
- `G`: Guinier prefactor [cm⁻¹]
- `P`: Power law slope [dimensionless]
- `B`: Porod constant [cm⁻¹ Å⁻ᴾ]
- `ETA`: Correlation distance [Å]
- `PACK`: Packing factor [dimensionless]
- `RgCO`: Cutoff radius [Å]
- `K`: Correction factor (auto-calculated)

**Flags:**
- `correlations`: Enable correlation function
- `mass_fractal`: Enable mass fractal mode
- `link_RGCO`: Link RgCO to previous level
- `link_B`: Auto-estimate B from G, Rg, P

**Fitting control:**
- `fit_Rg`, `fit_G`, `fit_P`, `fit_B`, etc.: Enable parameter fitting
- `Rg_limits`, `G_limits`, etc.: (min, max) bounds for each parameter

### UnifiedFitModel

Main class for the Unified fit model.

#### Methods

**`__init__(num_levels=1)`**
- Initialize model with specified number of levels

**`calculate_intensity(q)`**
- Calculate model intensity at given q values
- Returns: numpy array of intensities

**`fit(q, intensity, error=None, method='trf', max_iterations=1000, verbose=0)`**
- Fit model to experimental data
- Args:
  - `q`: Scattering vector [1/Å]
  - `intensity`: Measured intensity [cm⁻¹]
  - `error`: Uncertainties [cm⁻¹]
  - `method`: Optimization method ('trf', 'dogbox', 'lm')
  - `max_iterations`: Maximum iterations
  - `verbose`: Verbosity level (0, 1, 2)
- Returns: Dictionary with fit results

**`calculate_invariant(level_idx=0, q_min=0.0, q_max=None)`**
- Calculate scattering invariant for a level
- Returns: Dictionary with invariant and surface/volume ratio

**`get_parameter_summary()`**
- Generate formatted summary of parameters
- Returns: String with parameter table

### Utility Functions

**`load_data_from_nxcansas(file_path, use_slit_smeared=False)`**
- Load data from NXcanSAS HDF5 file

**`plot_fit_results(model, show_residuals=True, log_scale=True, save_path=None)`**
- Plot fit results with data and model

**`estimate_initial_parameters(q, intensity, num_levels=1)`**
- Estimate initial parameters from data

**`export_fit_results(model, filename, include_levels=True)`**
- Export results to text file

**`calculate_size_distribution_moments(model, level_idx=0)`**
- Calculate size parameters assuming spherical particles

## Examples

### Example 1: Single Level Fit

```python
from unified import UnifiedFitModel
import numpy as np

# Generate synthetic data
q = np.logspace(-3, 0, 100)
# ... (load your data here)

# Create and configure model
model = UnifiedFitModel(num_levels=1)
model.levels[0].Rg = 50.0
model.levels[0].G = 1000.0
model.levels[0].P = 4.0
model.levels[0].B = 1e-3

# Configure fitting
model.levels[0].fit_Rg = True
model.levels[0].fit_G = True
model.levels[0].fit_P = False  # Fix P = 4 (Porod)
model.levels[0].fit_B = True

# Set bounds
model.levels[0].Rg_limits = (1.0, 500.0)
model.levels[0].G_limits = (1.0, 1e6)

# Fit
results = model.fit(q, intensity, error)

print(f"Chi-squared: {results['chi_squared']:.4e}")
print(f"Reduced chi-squared: {results['reduced_chi_squared']:.4f}")
```

### Example 2: Mass Fractal

```python
model = UnifiedFitModel(num_levels=1)

# Enable mass fractal mode
model.levels[0].mass_fractal = True
model.levels[0].Rg = 30.0
model.levels[0].G = 500.0
model.levels[0].P = 2.5  # Fractal dimension

# In mass fractal mode, B is calculated automatically:
# B = (G × P / Rg^P) × exp(Γ(P/2))

results = model.fit(q, intensity, error)
```

### Example 3: Correlation Function

```python
model = UnifiedFitModel(num_levels=1)

# Enable correlations
model.levels[0].correlations = True
model.levels[0].ETA = 60.0    # Correlation distance
model.levels[0].PACK = 0.3    # Packing factor

# Fit correlation parameters
model.levels[0].fit_ETA = True
model.levels[0].fit_PACK = True

results = model.fit(q, intensity, error)
```

### Example 4: Hierarchical Structure (3 levels)

```python
model = UnifiedFitModel(num_levels=3)

# Level 1: Primary particles (smallest)
model.levels[0].Rg = 10.0
model.levels[0].G = 50.0
model.levels[0].P = 4.0

# Level 2: Aggregates
model.levels[1].Rg = 100.0
model.levels[1].G = 1000.0
model.levels[1].P = 3.0
model.levels[1].link_RGCO = True  # Link to level 1

# Level 3: Clusters (largest)
model.levels[2].Rg = 1000.0
model.levels[2].G = 50000.0
model.levels[2].P = 2.0
model.levels[2].link_RGCO = True  # Link to level 2

results = model.fit(q, intensity, error)
```

## Running Demonstrations

The package includes comprehensive demonstrations:

```bash
python unified_demo.py
```

This will run four demonstrations:
1. **Demo 1**: Single-level fit with synthetic data
2. **Demo 2**: Two-level hierarchical structure
3. **Demo 3**: Correlation function fitting
4. **Demo 4**: Loading real data from NXcanSAS file (optional)

## Parameter Guidelines

### Typical Ranges

| Parameter | Typical Range | Notes |
|-----------|--------------|-------|
| Rg | 1-10,000 Å | Structure size |
| G | 0.01-10⁶ cm⁻¹ | Depends on contrast and volume fraction |
| P | 0-6 | P=4 for sharp interface (Porod), P<3 for fractals |
| B | 10⁻²⁰-1 cm⁻¹ Å⁻ᴾ | Depends on P and surface area |
| ETA | Rg-10×Rg | Correlation distance |
| PACK | 0-0.64 | Packing factor (max ≈ 0.64 for hard spheres) |

### Physical Interpretation

- **P ≈ 4**: Sharp interface (Porod law)
- **P ≈ 3**: Surface fractal
- **P < 3**: Mass fractal
- **P > 4**: Diffuse interface or polydispersity effects

### K Factor (Auto-calculated)

- K = 1.0 for P > 3
- K = 1.06 for P ≤ 3 (mass fractals and weak decays)

## Advanced Features

### Slit Smearing

For USAXS data with slit geometry:

```python
model.use_slit_smearing = True
model.slit_length = 0.01  # [1/Å]

results = model.fit(q, intensity, error)
```

Note: Current implementation is simplified. For production, implement proper Lake (1967) integration.

### Parameter Linking

Link B to G, Rg, and P (Hammouda relationship):

```python
model.levels[0].link_B = True
# B will be calculated as: B = G × exp(-P/2) × (3P/2)^(P/2) / Rg^P
```

Link RgCO to previous level:

```python
model.levels[1].link_RGCO = True
# RgCO for level 2 = Rg of level 1
```

### Invariant Analysis

```python
inv = model.calculate_invariant(level_idx=0, q_min=0.0, q_max=1.0)

print(f"Invariant: {inv['invariant']:.4e} cm⁻¹ Å⁻³")
print(f"Invariant: {inv['invariant_cm4']:.4e} cm⁻⁴")

if inv['surface_to_volume'] is not None:
    print(f"S/V: {inv['surface_to_volume']:.4e} m²/cm³")
```

## Troubleshooting

### Common Issues

**1. Fit doesn't converge**
- Check initial parameter guesses (use `estimate_initial_parameters()`)
- Adjust parameter bounds
- Fix some parameters (e.g., P = 4 for Porod)
- Reduce number of fitted parameters initially

**2. Unphysical parameters**
- Tighten parameter bounds
- Enable parameter linking (link_B, link_RGCO)
- Check data quality and Q range

**3. High chi-squared**
- Consider adding more levels
- Enable correlation function
- Check for systematic errors in data

**4. Negative intensities**
- Check background level
- Adjust G and B bounds to positive values

## File Formats

### NXcanSAS Input

The code uses `hdf5code.readGenericNXcanSAS()` to read data from NXcanSAS-formatted HDF5 files. The function extracts:
- Q: Scattering vector [1/Å]
- Intensity: [cm⁻¹]
- Error: Uncertainties [cm⁻¹]
- dQ: Q resolution [1/Å]

### Output Files

**Text export** (`export_fit_results()`):
- Parameter summary
- Data table: Q, I_data, I_model, Residuals, Level contributions
- Tab-delimited format for easy import into Excel/Origin

## Performance Notes

- Typical fit time: 1-10 seconds for 100-1000 data points
- Multi-level fits (3+ levels) may require 10-60 seconds
- Use `verbose=1` or `verbose=2` to monitor progress
- For large datasets (>1000 points), consider rebinning data

## License

This implementation is based on the published work of Greg Beaucage and follows the mathematical formulations described in the references above.

## Contributing

Improvements welcome! Areas for contribution:
- Proper slit smearing implementation (Lake integration)
- Additional correlation functions
- GUI interface
- Parallel fitting for batch processing
- MCMC uncertainty estimation

## Contact

For questions about the Unified fit model theory, consult the references.
For issues with this Python implementation, check the code documentation.

## Acknowledgments

- Greg Beaucage for developing the Unified fit model
- Igor Pro Irena package (original implementation)
- Small-angle scattering community

---

**Version**: 1.0.0
**Date**: 2024
**Python**: 3.8+
