"""
Data input/output utilities for pyIrena.

This module provides functions for reading and writing various data formats
commonly used in small-angle scattering, including NXcanSAS HDF5 files.

Functions:
    readGenericNXcanSAS: Read data from NXcanSAS HDF5 files
    discover_scattering: Discover every curve in HDF5 or text sources
    load_scattering: Read one curve without modifying the source
    load_data_from_nxcansas: Convenience wrapper for loading data
    load_result: Load stored fit results from an NXcanSAS HDF5 file
"""

try:
    from pyirena.io.h5xp_writer import create_h5xp, write_iq_data
    from pyirena.io.hdf5 import readGenericNXcanSAS
    from pyirena.io.results import SUPPORTED_ANALYSES, load_result
    from pyirena.io.scattering import (
        ScatteringLocation,
        ScatteringRecord,
        discover_scattering,
        load_scattering,
    )
    __all__ = [
        "SUPPORTED_ANALYSES",
        "ScatteringLocation",
        "ScatteringRecord",
        "create_h5xp",
        "discover_scattering",
        "load_result",
        "load_scattering",
        "readGenericNXcanSAS",
        "write_iq_data",
    ]
except ImportError:
    # h5py might not be installed
    __all__ = []
