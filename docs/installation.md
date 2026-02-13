# Installation Guide

## Requirements

- Python 3.8 or higher
- pip or conda package manager

## Installation Methods

### 1. Install from PyPI (Recommended)

Once published, the easiest way to install pyIrena:

```bash
pip install pyirena
```

### 2. Install from Source

For the latest development version:

```bash
# Clone the repository
git clone https://github.com/ilavsky/pyirena.git
cd pyirena

# Install in development mode
pip install -e .
```

### 3. Install with Conda

Once published on conda-forge:

```bash
conda install -c conda-forge pyirena
```

## Optional Dependencies

### Plotting Support

For plotting capabilities with matplotlib:

```bash
pip install pyirena[plotting]
```

### Development Tools

For development and testing:

```bash
pip install pyirena[dev]
```

This includes:
- pytest for testing
- black for code formatting
- flake8 for linting
- mypy for type checking

### All Optional Dependencies

To install everything:

```bash
pip install pyirena[all]
```

## Verifying Installation

Test your installation:

```python
import pyirena
from pyirena.core.unified import UnifiedFitModel

print(f"pyIrena version: {pyirena.__version__}")

# Create a simple model
model = UnifiedFitModel(num_levels=1)
print("Installation successful!")
```

## Dependencies

### Core Dependencies (automatically installed)

- **numpy** (>=1.20.0): Numerical computing
- **scipy** (>=1.7.0): Optimization and scientific functions
- **h5py** (>=3.0.0): HDF5 file support

### Optional Dependencies

- **matplotlib** (>=3.3.0): Plotting and visualization
- **six** (>=1.15.0): Python 2/3 compatibility utilities

## Troubleshooting

### Common Issues

#### HDF5 Installation Issues

If you encounter problems installing h5py:

**On macOS:**
```bash
brew install hdf5
pip install h5py
```

**On Linux (Ubuntu/Debian):**
```bash
sudo apt-get install libhdf5-dev
pip install h5py
```

**On Windows:**
Use Anaconda/Miniconda which includes pre-built binaries:
```bash
conda install h5py
```

#### ImportError after installation

Make sure you're not in the source directory when importing:

```bash
cd ~
python -c "import pyirena; print(pyirena.__version__)"
```

#### Version conflicts

Create a fresh virtual environment:

```bash
python -m venv pyirena_env
source pyirena_env/bin/activate  # On Windows: pyirena_env\Scripts\activate
pip install pyirena
```

## Virtual Environments

### Using venv (Standard Library)

```bash
# Create environment
python -m venv myenv

# Activate (Linux/macOS)
source myenv/bin/activate

# Activate (Windows)
myenv\Scripts\activate

# Install pyirena
pip install pyirena
```

### Using conda

```bash
# Create environment
conda create -n pyirena_env python=3.10

# Activate
conda activate pyirena_env

# Install
pip install pyirena
# or
conda install -c conda-forge pyirena
```

## Updating

### Update from PyPI

```bash
pip install --upgrade pyirena
```

### Update from source

```bash
cd pyirena
git pull
pip install -e . --upgrade
```

## Uninstalling

```bash
pip uninstall pyirena
```

## Platform-Specific Notes

### Windows

- Recommend using Anaconda/Miniconda for easier dependency management
- Visual Studio Build Tools may be required for some dependencies

### macOS

- Xcode Command Line Tools recommended: `xcode-select --install`
- Use Homebrew for system dependencies: `brew install hdf5`

### Linux

- Install development headers for HDF5 and other libraries
- May need to install gcc/g++ compilers: `sudo apt-get install build-essential`

## Next Steps

After installation, check out:
- [Quick Start Guide](../README.md#quick-start)
- [Usage Guide](../USAGE_GUIDE.md)
- [Examples](../pyirena/examples/)
