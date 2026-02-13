"""
Visualization and plotting utilities for pyIrena.

This module provides functions for plotting and visualizing:
- Fitted data with model overlays
- Guinier and Porod analyses
- Residual plots
- Multi-level decomposition

Note: Requires matplotlib to be installed (install with: pip install pyirena[plotting])
"""

try:
    from pyirena.plotting.unified_plots import (
        plot_fit_results,
        plot_guinier_analysis,
        plot_porod_analysis,
    )
    __all__ = [
        "plot_fit_results",
        "plot_guinier_analysis",
        "plot_porod_analysis",
    ]
except ImportError:
    # matplotlib might not be installed
    __all__ = []
