# NXcanSAS Unified Fit Results Format

## Overview

This document describes the HDF5/NXcanSAS format used by pyirena to store Unified Fit analysis results. The format follows the NeXus/NXcanSAS standard for small-angle scattering data and extends it with a dedicated group for Unified Fit model parameters and results.

## File Structure

```
file.h5 (or file_NX.h5)
├── @default = "entry"
├── @creator = "pyirena"
├── @NeXus_version = "4.3.0"
├── @HDF5_version
├── @h5py_version
├── @file_name
└── @file_time
└── entry/                          [NXentry]
    ├── @NX_class = "NXentry"
    ├── @canSAS_class = "SASentry"
    ├── definition = "NXsas"
    │
    ├── {sample_name}/              [NXsubentry] (optional, original data)
    │   ├── @NX_class = "NXsubentry"
    │   ├── @canSAS_class = "SASentry"
    │   ├── @default = "sasdata"
    │   ├── definition = "NXcanSAS"
    │   └── sasdata/                [NXdata]
    │       ├── @signal = "I"
    │       ├── @I_axes = "Q"
    │       ├── I                   [dataset] Intensity (1/cm)
    │       ├── Q                   [dataset] Q vector (1/Å)
    │       └── Idev               [dataset] Uncertainties (1/cm)
    │
    └── unified_fit_results/        [NXprocess]
        ├── @NX_class = "NXprocess"
        ├── @analysis_type = "Unified Fit"
        ├── @program = "pyirena"
        ├── @timestamp
        ├── @num_levels             Number of structural levels
        ├── @background             Background value (1/cm)
        ├── @chi_squared            Chi-squared from fit
        │
        ├── Q                       [dataset] Q vector (1/Å)
        ├── intensity_data          [dataset] Experimental I(Q) (1/cm)
        ├── intensity_error         [dataset] Experimental σ(Q) (1/cm)
        ├── intensity_model         [dataset] Unified model I(Q) (1/cm)
        ├── residuals               [dataset] Normalized residuals
        │
        ├── level_1/                [group]
        │   ├── @level_number = 1
        │   ├── @G                  Guinier prefactor
        │   ├── @Rg                 Radius of gyration (Å)
        │   ├── @B                  Porod constant
        │   ├── @P                  Porod exponent
        │   ├── @RgCutoff           Rg cutoff (Å)
        │   ├── @ETA                Correlation hole size (Å)
        │   ├── @PACK               Packing factor
        │   ├── @correlated         Correlations enabled (bool)
        │   ├── @Sv                 Surface to volume ratio (m²/cm³)
        │   └── @Invariant          Invariant (cm⁻⁴)
        │
        ├── level_2/                [group]
        │   └── ... (same structure as level_1)
        │
        └── level_N/                [group]
            └── ... (same structure as level_1)
```

## Unified Fit Parameters

Each level stores the following parameters as HDF5 attributes:

### Core Parameters
- **G**: Guinier prefactor - amplitude of Guinier region
- **Rg**: Radius of gyration (Å) - characteristic size at this level
- **B**: Porod constant - amplitude of power-law region
- **P**: Porod exponent - power-law slope (typically 1-5)

### Optional Parameters
- **RgCutoff**: Cutoff radius (Å) - limits structure factor at high-Q
- **ETA**: Correlation hole size (Å) - for locally correlated structures
- **PACK**: Packing factor - degree of correlation (0-16)
- **correlated**: Boolean - whether correlation corrections are applied

### Calculated Values
- **Sv**: Surface to volume ratio (m²/cm³) - calculated when P ≈ 4
- **Invariant**: Scattering invariant (cm⁻⁴)

## Data Arrays

All arrays are stored as HDF5 datasets with appropriate attributes:

- **Q**: Scattering vector magnitude (1/Å)
- **intensity_data**: Experimental intensity (1/cm or cm⁻¹)
- **intensity_error**: Experimental uncertainty (same units as intensity)
- **intensity_model**: Calculated Unified Fit model (same units as intensity)
- **residuals**: Normalized residuals = (data - model) / error

## Units

Following NXcanSAS conventions:
- **Q**: 1/Å (inverse angstroms)
- **Intensity**: 1/cm or cm⁻¹ (per centimeter)
- **Rg, ETA, RgCutoff**: Å (angstroms)
- **Sv**: m²/cm³
- **Invariant**: cm⁻⁴

## File Naming Convention

- **NXcanSAS input**: Results saved in same file
- **Other formats** (.dat, .txt, non-NXcanSAS .h5): New file created as `{original_name}_NX.h5`

Example:
- Input: `mydata.dat` → Output: `mydata_NX.h5`
- Input: `mydata.h5` (NXcanSAS) → Output: `mydata.h5` (appended)

## Python API

### Saving Results

```python
from pyirena.io.nxcansas_unified import save_unified_fit_results

save_unified_fit_results(
    filepath=Path("mydata_NX.h5"),
    q=q_array,
    intensity_data=intensity_exp,
    intensity_model=intensity_calc,
    residuals=residuals,
    levels=[
        {
            'G': 7.05e7, 'Rg': 1040, 'B': 3.47e-4, 'P': 4.0,
            'RgCutoff': 0, 'ETA': 2990, 'PACK': 1.84,
            'correlated': True, 'Sv': 37.7, 'Invariant': 2.1e10
        },
        # ... more levels
    ],
    background=1e-6,
    chi_squared=1.23,
    num_levels=2,
    error=error_array  # optional
)
```

### Loading Results

```python
from pyirena.io.nxcansas_unified import load_unified_fit_results

results = load_unified_fit_results(Path("mydata_NX.h5"))

# Access data
q = results['Q']
intensity_data = results['intensity_data']
intensity_model = results['intensity_model']
residuals = results['residuals']

# Access fit parameters
num_levels = results['num_levels']
background = results['background']
chi_squared = results['chi_squared']

# Access level parameters
for i, level in enumerate(results['levels']):
    print(f"Level {i+1}:")
    print(f"  G = {level['G']:.3e}")
    print(f"  Rg = {level['Rg']:.3f} Å")
    print(f"  B = {level['B']:.3e}")
    print(f"  P = {level['P']:.3f}")
    if level.get('correlated', False):
        print(f"  ETA = {level['ETA']:.1f} Å")
        print(f"  PACK = {level['PACK']:.2f}")
```

### Command-Line Reader

A command-line utility is provided for quick inspection:

```bash
python -m pyirena.io.nxcansas_unified mydata_NX.h5
```

Output:
```
Loaded Unified Fit results from mydata_NX.h5
  - 2 levels
  - χ² = 1.2345
  - Timestamp: 2026-02-17T14:30:00

Level 1:
  G = 7.050e+07
  Rg = 1040.000 Å
  B = 3.470e-04
  P = 4.000
  ETA = 2990.0 Å
  PACK = 1.84
  Sv = 37.7 m²/cm³
  Invariant = 2.10e+10 cm⁻⁴

Level 2:
  G = 0.000e+00
  Rg = 10000000000.000 Å
  B = 2.330e-08
  P = 4.410
  ...
```

## Compatibility

This format is compatible with:
- NeXus 4.3.0 standard
- NXcanSAS application definition
- pyirena (creator and reader)
- Any HDF5-compatible software (with custom interpretation of pyirena-specific groups)

## Version History

- **1.0** (2026-02-17): Initial implementation
  - NXcanSAS base structure
  - Unified Fit results storage
  - Level parameters as attributes
  - Python save/load API
