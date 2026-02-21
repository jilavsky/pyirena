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

Loading stored results example:
    >>> from pyirena import load_result
    >>> r = load_result("mydata.h5", "unified_fit")
    >>> if r["found"]:
    ...     print(f"chi² = {r['chi_squared']:.4f}")
    >>> r = load_result("mydata.h5", "size_distribution")
    >>> if r["found"]:
    ...     print(f"Vf = {r['volume_fraction']:.4g},  method = {r['method']}")

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
from pyirena.io.results import load_result, SUPPORTED_ANALYSES

try:
    from pyirena.plotting.plot_saxs import plot_saxs
except ImportError:
    pass  # matplotlib not installed

__all__ = [
    "UnifiedFitModel",
    "UnifiedLevel",
    "load_data_from_nxcansas",
    "SizesDistribution",
    "fit_unified",
    "fit_sizes",
    "fit_pyirena",
    "load_result",
    "SUPPORTED_ANALYSES",
    "plot_saxs",
]
