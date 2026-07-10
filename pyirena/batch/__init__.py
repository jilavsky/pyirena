"""
pyirena.batch — headless (no-GUI) fitting API for scripting and automation.

Typical usage
-------------
Single file, one tool:
    from pyirena.batch import fit_unified
    result = fit_unified("data.h5", "pyirena_config.json")
    if result:
        print(result['parameters']['chi_squared'])

WAXS peak fitting (linear/linear, diffraction peaks):
    from pyirena.batch import fit_waxs
    result = fit_waxs("waxs_data.h5", "pyirena_config.json")
    if result and result['success']:
        print(result['n_peaks'], "peaks fitted")

Single file, all configured tools:
    from pyirena.batch import fit_pyirena
    results = fit_pyirena("data.h5", "pyirena_config.json")

Batch over many files:
    results = [fit_pyirena(f, cfg) for f in data_files]
    results = [r for r in results if r is not None]  # filter failures
"""

# Re-exports preserve the original `from pyirena.batch import X` API.
from pyirena.batch._common import _load_config, _load_data
from pyirena.batch.convert import igor_to_nexus, pxp_to_nexus
from pyirena.batch.manipulate import manipulate_data, average_data
from pyirena.batch.merge import merge_data
from pyirena.batch.modeling import fit_modeling
from pyirena.batch.pipeline import fit_pyirena
from pyirena.batch.saxs_morph import fit_saxs_morph
from pyirena.batch.simple import fit_simple, fit_simple_from_config
from pyirena.batch.sizes import _mc_uncertainty_sizes, fit_sizes
from pyirena.batch.unified import _state_to_model, _compute_invariant_sv, _mc_uncertainty_unified, _build_setup_state, _save_to_nexus, fit_unified
from pyirena.batch.waxs import fit_waxs_peaks_from_config, fit_waxs, fit_waxs_peaks

__all__ = [
    "fit_unified",
    "fit_sizes",
    "fit_pyirena",
    "fit_simple",
    "fit_simple_from_config",
    "fit_waxs_peaks_from_config",
    "fit_waxs",
    "fit_waxs_peaks",
    "merge_data",
    "fit_modeling",
    "fit_saxs_morph",
    "manipulate_data",
    "average_data",
    "igor_to_nexus",
    "pxp_to_nexus",
]
