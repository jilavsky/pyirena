# pyIrena

> Coded by Claude from SAXS_IgorCode Irena. Planned, defined, debugged and validated by Jan Ilavsky.

Python tools for small-angle scattering (SAS) data analysis. A port of the Igor Pro
[Irena](https://usaxs.xray.aps.anl.gov/software/irena) package. Includes interactive
GUI tools for fitting, modeling, data merging, and visualization of SAXS/SANS/USAXS data.

**Current release: v0.4.6 (public beta)**

[![PyPI version](https://img.shields.io/pypi/v/pyirena.svg)](https://pypi.org/project/pyirena/)
[![Python Version](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

## Installation

**From PyPI (recommended):**

```bash
pip install pyirena[gui]
```

This installs pyirena with all GUI dependencies (PySide6, pyqtgraph, etc.). For the
core library only (no GUI), use `pip install pyirena`.

**From source (for development):**

```bash
git clone https://github.com/jilavsky/pyirena.git
cd pyirena
pip install -e ".[gui]"
```

**With conda:**

```bash
git clone https://github.com/jilavsky/pyirena.git
cd pyirena
conda env create -f environment.yml
conda activate pyirena
```

See [docs/installation.md](docs/installation.md) for full details, troubleshooting, and platform notes.

---

## Running the GUI

```bash
pyirena-gui
```

This launches the Data Selector, the main entry point for all analysis tools.
See [docs/gui_quickstart.md](docs/gui_quickstart.md) for a walkthrough.

Individual tools can also be launched directly:

| Command | Tool |
|---------|------|
| `pyirena-gui` | Data Selector (main entry point) |
| `pyirena-viewer` | HDF5 Viewer / Data Extractor |
| `pyirena-modeling` | Modeling tool (standalone) |
| `pyirena-datamerge` | Data Merge tool (standalone) |
| `pyirena-contrast` | Scattering Contrast Calculator |

---

## Analysis Tools

### Modeling
Parametric forward-modeling of SAS data. Combine up to 10 populations, each of which
can be a Size Distribution, Unified Fit Level, or Diffraction Peak. Fits the combined
model to experimental data using least-squares optimization.

- **5 distribution functions**: Gaussian, LogNormal, LSW, Schulz-Zimm, Ardell
- **9 form factors**: Sphere, Spheroid, Cylinder (aspect ratio / fixed length),
  Core-Shell Sphere and Spheroid (by core R / shell t / total R)
- **2 structure factors**: Born-Green interferences, Hard Sphere (Percus-Yevick)
- Monte Carlo uncertainty estimation
- [Modeling GUI guide](docs/modeling_gui.md)

### Unified Fit
Beaucage hierarchical model (1995, 1996) with 1-5 structural levels, each combining
Guinier and power-law contributions. Optional Born-Green correlation function.

- [Unified Fit GUI guide](docs/unified_fit_gui.md)
- [Unified Fit features & parameters](docs/unified_fit_features.md)

### Size Distribution
Indirect Fourier transform to recover particle size distributions from SAS data.
Four inversion methods: MaxEnt, Regularization, TNNLS, and Monte Carlo.

- [Size Distribution methods](docs/sizes_methods.md)

### Simple Fits
13 direct analytical models: Guinier, Guinier-Porod, Porod, Sphere, Spheroid,
Debye-Bueche, Treubner-Strey, Power Law, and more. Each with linearization plots
and Monte Carlo uncertainty estimation.

- [Simple Fits GUI guide](docs/simple_fits_gui.md)

### WAXS Peak Fit
Fit diffraction peaks (Gaussian, Lorentzian, Pseudo-Voigt, Log-Normal) on linear
I vs Q scale. Auto-detect peaks via Savitzky-Golay + `scipy.signal.find_peaks`.
Simultaneous background fitting (constant, linear, cubic, or 5th-order polynomial).

- [WAXS Peak Fit GUI guide](docs/waxs_peakfit_gui.md)

### Data Merge
Merge two SAS datasets (e.g., SAXS + WAXS) onto a common Q scale. Optimizes scale
factor, flat background, and optional Q-shift using Nelder-Mead.

- [Data Merge GUI guide](docs/data_merge_gui.md)

### HDF5 Viewer / Data Extractor
Browse NXcanSAS HDF5 files, inspect raw data and analysis results, extract and plot
datasets. Supports all pyIrena result types (Unified Fit, Sizes, Simple Fits, WAXS
Peaks, Modeling).

- [HDF5 Viewer guide](docs/hdf5_viewer_gui.md)

### Scattering Contrast Calculator
Look up X-ray and neutron scattering length densities for materials by chemical
formula. Computes contrast (Delta-rho-squared) between two materials.

- [Contrast Calculator guide](docs/scattering_contrast_gui.md)

### Data Selector
Central GUI panel for managing data files and launching analysis tools. Load HDF5
files, select datasets, tabulate results across files, and generate reports.

---

## Batch Scripting API

All analysis tools can be run headlessly from Python scripts or JSON configuration files:

```python
from pyirena.batch import fit_pyirena

results = fit_pyirena(
    data_file='sample.h5',
    config_file='pyirena_config.json',
)
```

Individual functions: `fit_unified`, `fit_sizes`, `fit_simple`, `fit_waxs`,
`fit_modeling`, `merge_data`.

See [docs/batch_api.md](docs/batch_api.md) for the full scripting guide.

---

## NXcanSAS I/O

All data and results are stored in HDF5 files using the
[NXcanSAS](https://www.nexusformat.org/NXcanSAS.html) format. Fit results are saved
alongside raw data, making files self-contained and shareable.

- [NXcanSAS format details](docs/NXcanSAS_UnifiedFit_Format.md)

---

## Documentation

| Topic | File |
|-------|------|
| Installation | [docs/installation.md](docs/installation.md) |
| Quick start (GUI) | [docs/gui_quickstart.md](docs/gui_quickstart.md) |
| Modeling GUI guide | [docs/modeling_gui.md](docs/modeling_gui.md) |
| Unified Fit GUI guide | [docs/unified_fit_gui.md](docs/unified_fit_gui.md) |
| Unified Fit features & parameters | [docs/unified_fit_features.md](docs/unified_fit_features.md) |
| Size Distribution methods | [docs/sizes_methods.md](docs/sizes_methods.md) |
| Simple Fits GUI guide | [docs/simple_fits_gui.md](docs/simple_fits_gui.md) |
| WAXS Peak Fit GUI guide | [docs/waxs_peakfit_gui.md](docs/waxs_peakfit_gui.md) |
| Data Merge GUI guide | [docs/data_merge_gui.md](docs/data_merge_gui.md) |
| HDF5 Viewer guide | [docs/hdf5_viewer_gui.md](docs/hdf5_viewer_gui.md) |
| Contrast Calculator guide | [docs/scattering_contrast_gui.md](docs/scattering_contrast_gui.md) |
| NXcanSAS file format | [docs/NXcanSAS_UnifiedFit_Format.md](docs/NXcanSAS_UnifiedFit_Format.md) |
| Batch fitting API | [docs/batch_api.md](docs/batch_api.md) |
| Usage guide (scripting) | [docs/usage_guide.md](docs/usage_guide.md) |
| Developer: adding form factors | [docs/developer_adding_form_factors.md](docs/developer_adding_form_factors.md) |
| Developer: adding structure factors | [docs/developer_adding_structure_factors.md](docs/developer_adding_structure_factors.md) |
| Testing | [docs/testing.md](docs/testing.md) |
| Distribution / packaging | [docs/distribution.md](docs/distribution.md) |
| Contributing | [CONTRIBUTING.md](CONTRIBUTING.md) |
| Changelog | [CHANGELOG.md](CHANGELOG.md) |

---

## License

MIT -- see [LICENSE](LICENSE).

## Contact

- **Author**: Jan Ilavsky -- ilavsky@aps.anl.gov
- **Affiliation**: X-ray Science Division, Argonne National Laboratory
- **Issues**: https://github.com/jilavsky/pyirena/issues
