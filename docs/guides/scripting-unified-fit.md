# Unified Fit Model - Quick Usage Guide

## Installation

```bash
pip install -r requirements.txt
```

## Basic Workflow

### 1. Import the modules

```python
from unified import UnifiedFitModel, load_data_from_nxcansas
from unified_utils import plot_fit_results, export_fit_results
import numpy as np
```

### 2. Load your data

**Option A: From NXcanSAS HDF5 file**
```python
data = load_data_from_nxcansas('/path/to/your/data.h5')
q = data['Q']
intensity = data['Intensity']
error = data['Error']
```

**Option B: From arrays**
```python
# Your q, intensity, and error arrays
q = np.array([...])
intensity = np.array([...])
error = np.array([...])
```

### 3. Create and configure the model

```python
# Create model with desired number of levels
model = UnifiedFitModel(num_levels=1)

# Set initial parameter values
model.levels[0].Rg = 50.0      # Radius of gyration [Angstroms]
model.levels[0].G = 1000.0     # Guinier prefactor [cm^-1]
model.levels[0].P = 4.0        # Power law slope
model.levels[0].B = 1e-3       # Porod constant
model.background = 0.01        # Background [cm^-1]
```

### 4. Configure which parameters to fit

```python
# Enable/disable fitting for each parameter
model.levels[0].fit_Rg = True
model.levels[0].fit_G = True
model.levels[0].fit_P = False   # Fix P = 4 (Porod law)
model.levels[0].fit_B = True
model.fit_background = True

# Set parameter bounds
model.levels[0].Rg_limits = (1.0, 500.0)
model.levels[0].G_limits = (1.0, 1e6)
model.levels[0].B_limits = (1e-10, 1.0)
```

### 5. Perform the fit

```python
results = model.fit(q, intensity, error, verbose=1)

# Check if fit succeeded
if results['success']:
    print("Fit successful!")
    print(f"Chi-squared: {results['chi_squared']:.4e}")
    print(f"Reduced chi-squared: {results['reduced_chi_squared']:.4f}")
else:
    print(f"Fit failed: {results['message']}")
```

### 6. View results

```python
# Print parameter summary
print(model.get_parameter_summary())

# Access fitted parameters
print(f"Fitted Rg = {model.levels[0].Rg:.2f} Angstroms")
print(f"Fitted G = {model.levels[0].G:.2e} cm^-1")
```

### 7. Plot and export

```python
# Plot results
plot_fit_results(model, show_residuals=True, log_scale=True,
                 save_path='fit_results.png')

# Export to text file
export_fit_results(model, 'fit_results.txt', include_levels=True)
```

---

## Common Scenarios

### Scenario 1: Porod Scattering (Sharp Interface)

```python
model = UnifiedFitModel(num_levels=1)
model.levels[0].Rg = 50.0
model.levels[0].G = 1000.0
model.levels[0].P = 4.0        # Porod law
model.levels[0].B = 1e-3
model.levels[0].fit_P = False  # Fix P = 4

results = model.fit(q, intensity, error)

# Calculate surface/volume ratio
inv = model.calculate_invariant(0)
if inv['surface_to_volume'] is not None:
    print(f"Surface/Volume = {inv['surface_to_volume']:.4e} m²/cm³")
```

### Scenario 2: Mass Fractal

```python
model = UnifiedFitModel(num_levels=1)
model.levels[0].mass_fractal = True  # Enable mass fractal mode
model.levels[0].Rg = 30.0
model.levels[0].G = 500.0
model.levels[0].P = 2.5              # Fractal dimension
# B is calculated automatically in mass fractal mode

results = model.fit(q, intensity, error)
```

### Scenario 3: Two-Level Hierarchical Structure

```python
model = UnifiedFitModel(num_levels=2)

# Level 1: Primary particles
model.levels[0].Rg = 20.0
model.levels[0].G = 100.0
model.levels[0].P = 4.0
model.levels[0].B = 1e-3

# Level 2: Aggregates
model.levels[1].Rg = 200.0
model.levels[1].G = 10000.0
model.levels[1].P = 2.5
model.levels[1].link_RGCO = True  # Link cutoff to level 1 Rg

results = model.fit(q, intensity, error)
```

### Scenario 4: With Correlation Function

```python
model = UnifiedFitModel(num_levels=1)
model.levels[0].Rg = 30.0
model.levels[0].G = 500.0
model.levels[0].P = 4.0
model.levels[0].B = 5e-4

# Enable correlations
model.levels[0].correlations = True
model.levels[0].ETA = 60.0     # Correlation distance [Angstroms]
model.levels[0].PACK = 0.3     # Packing factor

# Fit correlation parameters
model.levels[0].fit_ETA = True
model.levels[0].fit_PACK = True

results = model.fit(q, intensity, error)
```

### Scenario 5: Automatic Parameter Estimation

```python
from unified_utils import estimate_initial_parameters, apply_parameters_from_dict

# Estimate initial parameters from data
initial_params = estimate_initial_parameters(q, intensity, num_levels=1)

# Apply to model
model = UnifiedFitModel(num_levels=1)
apply_parameters_from_dict(model, initial_params)

print("Initial guesses:")
print(model.get_parameter_summary())

# Then fit
results = model.fit(q, intensity, error)
```

---

## Tips for Successful Fitting

### 1. Start Simple
- Begin with 1 level
- Fix P = 4 initially
- Add complexity only if needed

### 2. Good Initial Guesses
- Use `estimate_initial_parameters()` for automatic estimation
- Or manually estimate:
  - **Rg**: From Guinier region (where ln(I) vs Q² is linear)
  - **G**: Intensity at lowest Q
  - **B**: From high-Q plateau in I×Q⁴ plot
  - **Background**: Intensity at highest Q

### 3. Parameter Bounds
- Always set reasonable bounds
- Too wide: slow convergence
- Too narrow: may miss true minimum

### 4. Fitting Strategy
- Start by fitting only Rg and G (fix P and B)
- Then enable B fitting
- Finally enable P fitting
- Add levels incrementally

### 5. Check Fit Quality
- Look at residuals plot
- Reduced χ² should be near 1.0
- If >> 1: model doesn't fit data well
- If << 1: errors may be overestimated

### 6. Physical Constraints
- Rg > 0 (obviously)
- G > 0 (always)
- 0 < P < 6 (typical)
- B > 0 (usually)
- 0 < PACK < 0.64 (hard sphere maximum)

---

## Parameter Interpretation

### Radius of Gyration (Rg)
- Measure of particle size
- For sphere: R = Rg × √(5/3)
- Typical range: 1-10,000 Å

### Guinier Prefactor (G)
- Low-q scattering amplitude
- Proportional to: N × V² × Δρ²
  - N = number density
  - V = particle volume
  - Δρ = contrast

### Power Law Slope (P)
- **P = 4**: Porod (sharp interface)
- **P = 3**: Surface fractal
- **P < 3**: Mass fractal (e.g., 1.8-2.5)
- **P > 4**: Diffuse interface

### Porod Constant (B)
- Related to surface area
- When P = 4: S/V = πB/Q (where Q is invariant)

### Correlation Parameters
- **ETA**: Average spacing between particles
- **PACK**: Strength of correlation (0 = no correlation)

---

## Troubleshooting

### "Fit doesn't converge"
1. Check initial guesses (try `estimate_initial_parameters`)
2. Reduce number of fitted parameters
3. Adjust bounds
4. Increase `max_iterations`

### "Unphysical parameters"
1. Tighten parameter bounds
2. Enable `link_B` for automatic B calculation
3. Fix problematic parameters

### "High chi-squared"
1. Try adding another level
2. Enable correlation function
3. Check data quality

### "Negative intensities in model"
1. Check background level
2. Adjust G and B to positive values
3. Review parameter bounds

---

## Advanced Features

### Invariant Calculation
```python
inv = model.calculate_invariant(
    level_idx=0,
    q_min=0.001,
    q_max=1.0
)
print(f"Invariant: {inv['invariant']:.4e} cm⁻¹ Å⁻³")
```

### Size Distribution Moments
```python
from unified_utils import calculate_size_distribution_moments

size_params = calculate_size_distribution_moments(model, level_idx=0)
print(f"Average radius (sphere): {size_params['R_avg_sphere']:.2f} Å")
print(f"Volume: {size_params['volume_sphere']:.2e} Å³")
```

### Slit Smearing (USAXS)
```python
model.use_slit_smearing = True
model.slit_length = 0.01  # [1/Angstrom]
results = model.fit(q, intensity, error)
```

### Multiple Plots
```python
from unified_utils import plot_guinier_analysis, plot_porod_analysis

# Guinier plot (ln I vs Q²)
plot_guinier_analysis(model, level_idx=0)

# Porod plot (I×Q⁴ vs Q)
plot_porod_analysis(model, level_idx=0)
```

---

## Running the Demo

To see complete examples:

```bash
python unified_demo.py
```

This runs four demonstrations:
1. Single-level fit
2. Two-level hierarchical structure
3. Correlation function
4. Real data from NXcanSAS (if you provide a file path)

---

## Quick Reference Card

| Task | Code |
|------|------|
| Create model | `model = UnifiedFitModel(num_levels=1)` |
| Set parameter | `model.levels[0].Rg = 50.0` |
| Enable fitting | `model.levels[0].fit_Rg = True` |
| Set bounds | `model.levels[0].Rg_limits = (1.0, 500.0)` |
| Fit data | `results = model.fit(q, intensity, error)` |
| Print results | `print(model.get_parameter_summary())` |
| Plot fit | `plot_fit_results(model)` |
| Export results | `export_fit_results(model, 'out.txt')` |
| Calculate invariant | `inv = model.calculate_invariant(0)` |
| Enable mass fractal | `model.levels[0].mass_fractal = True` |
| Enable correlations | `model.levels[0].correlations = True` |
| Link levels | `model.levels[1].link_RGCO = True` |

---

## Need Help?

1. Read the full documentation: [README_unified.md](README_unified.md)
2. Run the demos: `python unified_demo.py`
3. Check the references for theory details
4. Examine the code comments in `unified.py`

Happy fitting!
