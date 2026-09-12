"""pyirena.api — stable, AI-friendly facade over pyirena's HDF5 readers.

This package is the recommended entry point for non-GUI consumers (control
agents, automation scripts, LLM tools). Every public function:

- accepts ordinary str/path inputs,
- returns a JSON-serializable dict (numpy arrays decimated to a bounded
  length and NaN/inf replaced with None),
- never imports Qt / pyqtgraph,
- does not write to disk except for plot helpers (which write PNGs to a
  user-specified or temp directory).

Environment overrides
---------------------
PYIRENA_DATA_ROOT
    Restrict all file access to this subtree. Strongly recommended when
    exposing the api over a network/MCP boundary.
PYIRENA_MAX_ARRAY_POINTS
    Default decimation cap for returned arrays. Default 500.
PYIRENA_PLOT_CACHE
    Directory for plot PNGs. Default: ``<tempdir>/pyirena-mcp``.

Quick reference
---------------
Discovery
    list_files, summarize_folder, inspect_file

Reading
    read_reduced_data, read_metadata,
    read_simple_fit, read_unified_fit, read_size_distribution,
    read_modeling, read_saxs_morph, read_waxs_peakfit,
    read_fractals, read_merge_provenance, read_manipulation_provenance

Aggregation
    tabulate_parameter, summarize_sample

Plotting
    plot_iq, plot_parameter_trend

Calculators (stateless; no dataset, no file access)
    calc_contrast, calc_compound, calc_contrast_energy_scan,
    lookup_element, list_compound_library, load_compound

Data operations (the only group that WRITES data files)
    average_data, subtract_data, divide_data, scale_data, trim_data,
    rebin_data, merge_datasets, match_merge_files
"""
from __future__ import annotations

from pyirena.api.aggregate import summarize_sample, tabulate_parameter
from pyirena.api.calculators import (
    calc_compound,
    calc_contrast,
    calc_contrast_energy_scan,
    list_compound_library,
    load_compound,
    lookup_element,
)
from pyirena.api.data import read_metadata, read_reduced_data
from pyirena.api.data_ops import (
    average_data,
    divide_data,
    match_merge_files,
    merge_datasets,
    rebin_data,
    scale_data,
    subtract_data,
    trim_data,
)
from pyirena.api.discovery import inspect_file, list_files, summarize_folder
from pyirena.api.plotting import plot_iq, plot_parameter_trend
from pyirena.api.results import (
    read_fractals,
    read_manipulation_provenance,
    read_merge_provenance,
    read_modeling,
    read_saxs_morph,
    read_simple_fit,
    read_size_distribution,
    read_unified_fit,
    read_waxs_peakfit,
)

__all__ = [
    # discovery
    "list_files", "summarize_folder", "inspect_file",
    # data
    "read_reduced_data", "read_metadata",
    # per-tool results
    "read_simple_fit", "read_unified_fit", "read_size_distribution",
    "read_modeling", "read_saxs_morph", "read_waxs_peakfit",
    "read_fractals", "read_merge_provenance", "read_manipulation_provenance",
    # aggregation
    "tabulate_parameter", "summarize_sample",
    # plotting
    "plot_iq", "plot_parameter_trend",
    # calculators
    "calc_contrast", "calc_compound", "calc_contrast_energy_scan",
    "lookup_element", "list_compound_library", "load_compound",
    # data operations (write files)
    "average_data", "subtract_data", "divide_data", "scale_data",
    "trim_data", "rebin_data", "merge_datasets", "match_merge_files",
]
