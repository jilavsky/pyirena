# Installation Guide

## Requirements

- Python 3.10 or higher
- conda (Anaconda or Miniconda) **or** pip

## Recommended: Conda Environment (one command)

The repository includes `environment.yml` which creates a complete conda
environment with all dependencies — including the GUI stack — in one step:

```bash
# Clone the repository
git clone https://github.com/jilavsky/pyirena.git
cd pyirena

# Create the environment and install pyirena in editable mode
conda env create -f environment.yml

# Activate
conda activate pyirena
```

That's it. The environment file installs Python, all scientific packages,
PySide6, pyqtgraph, and pyirena itself (editable/development mode).

### Updating an existing environment

If the repository has been updated and new dependencies were added:

```bash
conda env update -f environment.yml --prune
```

### Removing the environment

```bash
conda env remove -n pyirena
```

---

## Alternative: pip install (no conda)

### From source (development mode)

```bash
git clone https://github.com/jilavsky/pyirena.git
cd pyirena

# GUI install (includes PySide6 + pyqtgraph)
pip install -e ".[gui]"
```

### Core only (no GUI)

```bash
pip install -e .
```

---

## Dependencies

### Core (always required)

| Package | Version | Purpose |
|---------|---------|---------|
| numpy | ≥ 1.20 | Numerical arrays |
| scipy | ≥ 1.7 | Optimization, fitting |
| h5py | ≥ 3.0 | HDF5 / NXcanSAS file I/O |

### GUI (required to run the interactive GUI)

| Package | Version | Purpose |
|---------|---------|---------|
| PySide6 | ≥ 6.4 | Qt6 Python bindings |
| pyqtgraph | ≥ 0.13 | Fast scientific plotting |
| matplotlib | ≥ 3.3 | Additional plot export |

### Development / testing (optional)

| Package | Purpose |
|---------|---------|
| pytest ≥ 7.0 | Test runner |
| pytest-cov | Coverage reports |

---

## Verifying the Installation

```bash
# Check the package is importable
python -c "import pyirena; print('pyirena', pyirena.__version__)"

# Launch the GUI
pyirena-gui
```

---

## Troubleshooting

### HDF5 / h5py build errors

`conda env create` installs pre-built binaries and avoids this entirely.
If using pip on Linux:

```bash
sudo apt-get install libhdf5-dev   # Debian/Ubuntu
pip install h5py
```

On macOS:
```bash
brew install hdf5
pip install h5py
```

### PySide6 / GUI not found

Make sure you installed with `[gui]` extras or that the conda environment
was created from `environment.yml`:

```bash
pip install -e ".[gui]"        # pip route
# or
conda env create -f environment.yml   # conda route
```

### ImportError after installation

Make sure the `pyirena` conda environment is active:

```bash
conda activate pyirena
python -c "import pyirena"
```

---

## Platform Notes

### Windows

Conda is strongly recommended — it provides pre-built PySide6 and h5py
binaries that avoid compiler requirements.

### macOS

- Intel and Apple Silicon both supported via conda-forge packages.
- If building from source with pip, Xcode Command Line Tools are required:
  `xcode-select --install`

### Linux

- Conda is the easiest path; all packages are available on conda-forge.
- With pip, install system HDF5 headers first (see Troubleshooting above).

---

## Next Steps

After installation, see:

- [Quick Start GUI guide](gui_quickstart.md) — open the GUI and load data
- [Unified Fit guide](unified_fit_gui.md)
- [WAXS Peak Fit guide](waxs_peakfit_gui.md)
- [Batch scripting API](batch_api.md)
