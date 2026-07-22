# Installation Guide

## Requirements

- Python 3.10 – 3.13 (Python 3.14 is not yet recommended — many compiled-binding wheels are still flaky)
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

That's it. The environment file installs Python and the scientific stack
via conda, then uses pip (inside conda) to install pyirena editable plus
the entire Qt6/GUI stack (PySide6, pyqtgraph, etc.) from PyPI. This
"conda for science, pip for Qt" split is intentional — see the warning
below.

### ⚠️ Do not mix conda-installed and pip-installed Qt

Never run `conda install pyside6` (or `pyqtgraph`) inside the pyirena env,
and do not add Qt packages to `environment.yml`'s conda dependencies. If
conda's `pyside6` and pip's `PySide6` end up coexisting, Windows will
fail at GUI launch with:

```
ImportError: DLL load failed while importing QtCore:
The specified procedure could not be found.
```

This is because `shiboken6` from one version pairs with Qt6 DLLs from the
other — the binaries don't match. Always let pip install the Qt stack
into a pyirena env. If you suspect a mix, see Troubleshooting below.

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

### Trying a pre-release (beta) version from PyPI

Beta releases (e.g. `1.1.0b1`) are published to PyPI ahead of a stable
release for early testing. `pip` ignores pre-releases by default, so pass
`--pre` explicitly:

```bash
pip install --pre "pyirena[gui]"
```

Or pin an exact beta version:

```bash
pip install "pyirena[gui]==1.1.0b1"
```

See [CHANGELOG.md](../CHANGELOG.md) for what changed, and report issues at
https://github.com/jilavsky/pyirena/issues.

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

### "DLL load failed while importing QtCore" (Windows)

Symptom — `pyirena-gui` prints:

```
Error: GUI dependencies not installed.
Details: Neither PySide6 nor PyQt6 found.
```

But `pip show PySide6` reports it is installed. Confirm with:

```bash
python -c "from PySide6 import QtCore"
```

If that gives `DLL load failed ... The specified procedure could not be
found`, you have a mixed conda + pip Qt install. Easiest fix is to
recreate the env so pip is the only Qt installer:

```bash
conda deactivate
conda env remove -n pyirena -y
conda env create -f environment.yml
conda activate pyirena
pyirena-gui
```

(The current `environment.yml` already routes the GUI stack through pip
to prevent this — older envs created before that change are the ones at
risk.)

### Python 3.14 wheels are flaky

Python 3.14 was released October 2025. Several scientific Python wheels
(including PySide6 on Windows) advertise cp314 support but in practice
fail with DLL or symbol-resolution errors. Stick with Python 3.10 – 3.13
until your stack of choice has caught up. The `environment.yml` caps
Python at `<3.14` for this reason.

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
