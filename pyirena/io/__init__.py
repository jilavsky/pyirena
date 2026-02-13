"""
Data input/output utilities for pyIrena.

This module provides functions for reading and writing various data formats
commonly used in small-angle scattering, including NXcanSAS HDF5 files.

Functions:
    readGenericNXcanSAS: Read data from NXcanSAS HDF5 files
    load_data_from_nxcansas: Convenience wrapper for loading data
"""

try:
    from pyirena.io.hdf5 import readGenericNXcanSAS
    __all__ = ["readGenericNXcanSAS"]
except ImportError:
    # h5py might not be installed
    __all__ = []
