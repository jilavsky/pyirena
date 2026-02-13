# Code Fragments Archive

This directory contains original/legacy code files that have been superseded by the organized package structure.

## Original Files (Pre-Restructuring)

These files were the original implementation before the package was restructured:

- **`unified.py`** - Original unified fit implementation
  - **Now**: `pyirena/core/unified.py`
  - Contains the UnifiedFitModel class and main fitting algorithms

- **`unified_utils.py`** - Original utility functions
  - **Now**: `pyirena/core/unified_utils.py`
  - Contains helper functions for unified fitting

- **`unified_demo.py`** - Original demonstration script
  - **Now**: Use `pyirena/examples/` for new examples
  - Shows how to use the unified fit model

- **`hdf5code.py`** - Original HDF5 I/O code
  - **Now**: `pyirena/io/hdf5.py`
  - Contains NXcanSAS file reading/writing functions

## Status

These files are **kept for reference only** and may be removed in future versions. All functionality has been migrated to the proper package structure under `pyirena/`.

## Current Active Files

The active, maintained code is now in:
- `pyirena/core/` - Core algorithms
- `pyirena/io/` - Input/Output handling
- `pyirena/plotting/` - Plotting utilities
- `pyirena/gui/` - GUI components

For new development, use the organized package structure, not these legacy files.
