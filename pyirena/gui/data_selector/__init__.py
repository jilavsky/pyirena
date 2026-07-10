"""
Data Selector GUI for pyIrena.

This module provides a GUI panel for selecting data files and displaying
their content as graphs.
"""

# Package split of the former monolithic data_selector.py.
# All previous top-level names are re-exported so existing imports
# (pyirena.gui.data_selector.DataSelectorPanel, ._build_report, .main, ...)
# keep working unchanged.
from pyirena.gui.data_selector.config_dialogs import ConfigManagerDialog, DataSelectorConfigDialog
from pyirena.gui.data_selector.igor_import import _IgorImportDialog
from pyirena.gui.data_selector.panel import DataSelectorPanel, main
from pyirena.gui.data_selector.plot_utils import _gen_colors, _legend_indices, _LogDecadeAxis, _style_plot, _iq_error_bars, _add_jpeg_export, _rescaled_view
from pyirena.gui.data_selector.report import _quality_report_rows, _build_report
from pyirena.gui.data_selector.results_windows import GraphWindow, UnifiedFitResultsWindow, SizeDistResultsWindow, SimpleFitResultsWindow, WAXSPeakFitResultsWindow, TabulateResultsWindow
from pyirena.gui.data_selector.sorting import _sort_key_name, _sort_key_temperature, _sort_key_time, _sort_key_order, _sort_key_pressure, _SORT_KEYS
from pyirena.gui.data_selector.workers import BatchWorker

__all__ = [
    "_gen_colors",
    "_legend_indices",
    "_LogDecadeAxis",
    "_style_plot",
    "_iq_error_bars",
    "_sort_key_name",
    "_sort_key_temperature",
    "_sort_key_time",
    "_sort_key_order",
    "_sort_key_pressure",
    "_SORT_KEYS",
    "_add_jpeg_export",
    "_rescaled_view",
    "_quality_report_rows",
    "_build_report",
    "ConfigManagerDialog",
    "DataSelectorConfigDialog",
    "GraphWindow",
    "UnifiedFitResultsWindow",
    "SizeDistResultsWindow",
    "SimpleFitResultsWindow",
    "WAXSPeakFitResultsWindow",
    "TabulateResultsWindow",
    "BatchWorker",
    "_IgorImportDialog",
    "DataSelectorPanel",
    "main",
]
