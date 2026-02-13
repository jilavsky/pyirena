# pyIrena Package Structure

This document describes the organization of the pyIrena package repository.

## Directory Structure

```
pyirena/
├── .github/                    # GitHub-specific files
│   ├── workflows/
│   │   └── tests.yml          # CI/CD testing workflow
│   └── ISSUE_TEMPLATE/
│       ├── bug_report.md      # Bug report template
│       └── feature_request.md # Feature request template
│
├── conda/                      # Conda package recipe
│   └── meta.yaml              # Conda build configuration
│
├── docs/                       # Documentation
│   ├── installation.md        # Installation instructions
│   └── QUICK_START.md         # Quick start guide
│
├── pyirena/                    # Main package directory
│   ├── __init__.py            # Package initialization
│   ├── py.typed               # PEP 561 type marker
│   │
│   ├── core/                  # Core analysis modules
│   │   ├── __init__.py
│   │   └── unified.py         # Unified Fit model
│   │
│   ├── io/                    # Input/output utilities
│   │   ├── __init__.py
│   │   └── hdf5.py            # HDF5/NXcanSAS support
│   │
│   ├── plotting/              # Visualization tools
│   │   ├── __init__.py
│   │   └── unified_plots.py   # Plotting functions
│   │
│   ├── examples/              # Example scripts
│   │   ├── __init__.py
│   │   └── basic_demo.py      # Basic usage examples
│   │
│   └── tests/                 # Unit tests
│       ├── __init__.py
│       └── test_unified.py    # Tests for Unified model
│
├── .gitignore                 # Git ignore patterns
├── CHANGELOG.md               # Version history
├── CONTRIBUTING.md            # Contribution guidelines
├── LICENSE                    # MIT License
├── MANIFEST.in                # Package manifest
├── README.md                  # Main documentation
├── pyproject.toml             # Modern Python packaging config
├── requirements.txt           # Dependencies
└── setup.py                   # Setup script (backward compatibility)
```

## Old Files (Not in Package)

The following files remain in the root directory but are NOT part of the package:

- `unified.py` → Replaced by `pyirena/core/unified.py`
- `hdf5code.py` → Replaced by `pyirena/io/hdf5.py`
- `unified_utils.py` → Replaced by `pyirena/plotting/unified_plots.py`
- `unified_demo.py` → Replaced by `pyirena/examples/basic_demo.py`
- `README_unified.md` → Replaced by `README.md`
- `USAGE_GUIDE.md` → Kept for reference
- `FILE_SUMMARY.txt` → Archive
- `notes.txt` → Archive
- `*.ipf` files → Igor Pro source files (archive)

You can safely delete these old files or keep them as archives.

## Key Configuration Files

### pyproject.toml
Modern Python packaging configuration following PEP 517/518:
- Project metadata
- Dependencies
- Build system configuration
- Tool configurations (black, pytest, etc.)

### setup.py
Minimal setup script for backward compatibility. Actual configuration is in `pyproject.toml`.

### MANIFEST.in
Specifies which non-Python files to include in distributions.

### requirements.txt
Lists runtime dependencies. Useful for pip install -r requirements.txt.

## Installation Methods

### For Users (when published)

```bash
# From PyPI
pip install pyirena

# From GitHub
pip install git+https://github.com/ilavsky/pyirena.git

# From Conda
conda install -c conda-forge pyirena
```

### For Developers

```bash
# Clone repository
git clone https://github.com/ilavsky/pyirena.git
cd pyirena

# Install in editable mode
pip install -e ".[dev]"
```

## Distribution Workflow

### PyPI Distribution

1. Update version in `pyproject.toml`
2. Update `CHANGELOG.md`
3. Build distribution:
   ```bash
   python -m build
   ```
4. Upload to PyPI:
   ```bash
   python -m twine upload dist/*
   ```

### Conda Distribution

1. Update version in `conda/meta.yaml`
2. Create conda package:
   ```bash
   conda build conda/
   ```
3. Submit recipe to conda-forge feedstock

## Testing

Run tests with:

```bash
# All tests
pytest

# With coverage
pytest --cov=pyirena

# Specific test file
pytest pyirena/tests/test_unified.py
```

## Code Quality

```bash
# Format code
black pyirena/

# Check linting
flake8 pyirena/

# Type checking
mypy pyirena/
```

## Continuous Integration

GitHub Actions automatically:
- Runs tests on Python 3.8-3.12
- Tests on Linux, macOS, and Windows
- Generates coverage reports
- Triggered on push and pull requests

## Documentation

Documentation is organized as:
- `README.md`: Main package documentation
- `docs/installation.md`: Installation guide
- `docs/QUICK_START.md`: Quick start tutorial
- `USAGE_GUIDE.md`: Detailed usage guide (legacy)
- `CONTRIBUTING.md`: Contribution guidelines
- Code docstrings: API documentation

## Import Structure

Users can import from:

```python
# Main imports
from pyirena import UnifiedFitModel, UnifiedLevel

# Direct module imports
from pyirena.core.unified import UnifiedFitModel
from pyirena.io.hdf5 import readGenericNXcanSAS
from pyirena.plotting.unified_plots import plot_fit_results
```

## Version Management

Version is defined in `pyproject.toml` and accessed as:

```python
import pyirena
print(pyirena.__version__)  # "0.1.0"
```

## License

MIT License - See LICENSE file for details.

## Maintenance

Regular tasks:
1. Update dependencies in `pyproject.toml` and `requirements.txt`
2. Add new tests for new features
3. Update `CHANGELOG.md` for each release
4. Keep documentation synchronized with code
5. Review and merge pull requests
6. Tag releases in git

## Future Structure

Potential additions:
- `pyirena/gui/`: Graphical interface
- `pyirena/models/`: Additional scattering models
- `pyirena/utils/`: General utilities
- `docs/api/`: Auto-generated API docs
- `notebooks/`: Jupyter notebook tutorials
