"""
Data input/output utilities for pyIrena.

This module provides functions for reading and writing various data formats
commonly used in small-angle scattering, including NXcanSAS HDF5 files.

Functions:
    readGenericNXcanSAS: Read data from NXcanSAS HDF5 files
    load_data_from_nxcansas: Convenience wrapper for loading data
    load_result: Load stored fit results from an NXcanSAS HDF5 file
"""

try:
    from pyirena.io.hdf5 import readGenericNXcanSAS
    from pyirena.io.results import load_result, SUPPORTED_ANALYSES
    __all__ = ["readGenericNXcanSAS", "load_result", "SUPPORTED_ANALYSES"]
except ImportError:
    # h5py might not be installed
    __all__ = []
