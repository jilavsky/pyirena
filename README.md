# pyIrena

**DO NOT use, not ready, not tested, used developement**
**Coded by Claude from SAXS_IgorCode Irena**

**Python tools for small-angle scattering data analysis**

[![Python Version](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

`pyIrena` is a comprehensive Python package for analyzing small-angle scattering (SAXS/SANS/USAXS) data. It provides powerful tools for fitting and interpreting hierarchical structures using the Unified Fit model (Beaucage method) and other advanced analysis techniques.

### Key Features

- **Unified Fit Model**: Full implementation of the Beaucage Unified fit model for hierarchical structures
- **Multi-Level Analysis**: Support for 1-5 structural levels with automatic parameter linking
- **Advanced Features**:
  - Mass fractal mode with automatic B calculation
  - Correlation functions (Born-Green approximation)
  - Slit smearing support for USAXS
  - Parameter constraints and bounds
- **Data I/O**: Native support for NXcanSAS HDF5 files
- **Analysis Tools**: Invariant calculation, Guinier/Porod analysis, size distribution moments
- **Visualization**: Comprehensive plotting utilities with matplotlib integration

## Installation

### From PyPI (when published)

```bash
pip install pyirena
```

### From source

```bash
git clone https://github.com/jilavsky/pyirena.git
cd pyirena
pip install -e .
```

### With optional dependencies

```bash
# Install with plotting support
pip install pyirena[plotting]

# Install with development tools
pip install pyirena[dev]

# Install everything
pip install pyirena[all]
```

### For Conda (when published)

```bash
conda install -c conda-forge pyirena
```

## Quick Start

### Basic Unified Fit Example

```python
import numpy as np
from pyirena.core.unified import UnifiedFitModel

# Create model with 1 structural level
model = UnifiedFitModel(num_levels=1)

# Set initial parameters
model.levels[0].Rg = 50.0    # Radius of gyration [Å]
model.levels[0].G = 1000.0   # Guinier prefactor [cm^-1]
model.levels[0].P = 4.0      # Power law slope
model.levels[0].B = 1e-3     # Porod constant
model.background = 0.01      # Flat background

# Fit to your data
results = model.fit(q_data, intensity_data, error_data)

# Display results
print(model.get_parameter_summary())
```

### Loading Data from NXcanSAS

```python
from pyirena.io.hdf5 import load_data_from_nxcansas

# Load data from HDF5 file
data = load_data_from_nxcansas('path/to/file.h5')

# Access arrays
q = data['Q']
intensity = data['Intensity']
error = data['Error']
```

### Multi-Level Hierarchical Fitting

```python
# Create model with 2 levels
model = UnifiedFitModel(num_levels=2)

# Level 1: Primary particles
model.levels[0].Rg = 20.0
model.levels[0].G = 100.0
model.levels[0].P = 4.0

# Level 2: Aggregates
model.levels[1].Rg = 200.0
model.levels[1].G = 10000.0
model.levels[1].P = 2.5  # Mass fractal
model.levels[1].link_RGCO = True  # Link cutoff to level 1

# Perform fit
results = model.fit(q, intensity, error)
```

## Documentation

Comprehensive documentation is available in the [docs](docs/) directory and includes:

- [Installation Guide](docs/installation.md)
- [User Guide](USAGE_GUIDE.md)
- [API Reference](docs/api.md)
- [Examples and Tutorials](pyirena/examples/)
- [Theory Background](docs/theory.md)

## Mathematical Model

The Unified fit model combines multiple structural levels using Guinier-Porod crossover functions:

For each structural level *i*:

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

## Package Structure

```
pyirena/
├── core/               # Core analysis modules
│   ├── __init__.py
│   └── unified.py     # Unified Fit model implementation
├── io/                # Data input/output
│   ├── __init__.py
│   └── hdf5.py        # HDF5/NXcanSAS support
├── plotting/          # Visualization tools
│   ├── __init__.py
│   └── unified_plots.py
├── examples/          # Example scripts and notebooks
│   ├── basic_fit.py
│   ├── multi_level.py
│   └── notebooks/
└── tests/            # Unit tests
    ├── test_unified.py
    └── test_io.py
```

## Examples

The package includes comprehensive examples in the [pyirena/examples](pyirena/examples/) directory:

- **basic_fit.py**: Simple single-level fitting
- **multi_level.py**: Hierarchical structure analysis
- **mass_fractal.py**: Mass fractal mode demonstration
- **correlation.py**: Correlation function fitting
- **notebooks/**: Jupyter notebooks with detailed walkthroughs

Run examples:

```bash
cd pyirena/examples
python basic_fit.py
```

## References

The Unified fit model is based on the work of Greg Beaucage:

- Beaucage, G. (1995). *J. Appl. Cryst.* **28**, 717-728
- Beaucage, G. (1996). *J. Appl. Cryst.* **29**, 134-146
- Beaucage, G. et al. (1997). *Macromolecules* **30**, 4158-4166

Online resources:
- http://www.eng.uc.edu/~gbeaucag/PDFPapers/Beaucage2.pdf
- http://www.eng.uc.edu/~gbeaucag/PDFPapers/Beaucage1.pdf

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. Areas for contribution:

- Additional analysis tools
- Proper slit smearing implementation (Lake integration)
- Additional correlation functions
- GUI interface
- Parallel fitting for batch processing
- MCMC uncertainty estimation

## Development

### Setting up development environment

```bash
git clone https://github.com/jilavsky/pyirena.git
cd pyirena
pip install -e ".[dev]"
```

### Running tests

```bash
pytest
```

### Code formatting

```bash
black pyirena/
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Greg Beaucage for developing the Unified fit model
- Igor Pro Irena package (original implementation by Jan Ilavsky)
- Small-angle scattering community at APS and beyond

## Citation

If you use pyIrena in your research, please cite:

```bibtex
@software{pyirena2024,
  author = {Ilavsky, Jan},
  title = {pyIrena: Python tools for small-angle scattering data analysis},
  year = {2024},
  url = {https://github.com/jilavsky/pyirena}
}
```

## Contact

- **Author**: Jan Ilavsky
- **Email**: ilavsky@aps.anl.gov
- **GitHub**: https://github.com/jilavsky/pyirena
- **Issues**: https://github.com/jilavsky/pyirena/issues

## Version

Current version: **0.1.0** (Alpha)

---

**Note**: This package is in active development. APIs may change between versions. For production use, please pin to a specific version.
