"""
pyIrena: Python tools for small-angle scattering data analysis

This package provides comprehensive tools for analyzing SAXS/SANS/USAXS data,
including the Unified Fit model (Beaucage method) for hierarchical structures.

Modules:
    core: Core analysis tools including the Unified Fit model
    io: Data input/output utilities for HDF5/NXcanSAS files
    plotting: Visualization and plotting functions

Example:
    >>> from pyirena.core.unified import UnifiedFitModel
    >>> model = UnifiedFitModel(num_levels=1)
    >>> model.levels[0].Rg = 50.0
    >>> results = model.fit(q_data, intensity_data)

References:
    Beaucage, G. (1995). J. Appl. Cryst. 28, 717-728
    Beaucage, G. (1996). J. Appl. Cryst. 29, 134-146
"""

__version__ = "0.1.1"
__author__ = "Jan Ilavsky"
__email__ = "ilavsky@aps.anl.gov"

from pyirena.core.unified import UnifiedFitModel, UnifiedLevel, load_data_from_nxcansas
from pyirena.core.sizes import SizesDistribution
from pyirena.batch import fit_unified, fit_sizes, fit_pyirena

__all__ = [
    "UnifiedFitModel",
    "UnifiedLevel",
    "load_data_from_nxcansas",
    "SizesDistribution",
    "fit_unified",
    "fit_sizes",
    "fit_pyirena",
]
