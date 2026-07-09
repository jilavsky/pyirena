"""
pyirena.batch.pipeline — fit_pyirena — run every tool configured in a pyirena_config.json.

Split from the original monolithic pyirena/batch.py (no behavior change).
"""

from __future__ import annotations

import logging
import traceback
from pathlib import Path
from typing import TYPE_CHECKING, Dict, List, Optional, Union


from pyirena.logging_setup import ensure_console_output as _ensure_console

if TYPE_CHECKING:
    pass

log = logging.getLogger(__name__)

from pyirena.batch._common import _load_config
from pyirena.batch.modeling import fit_modeling
from pyirena.batch.simple import fit_simple_from_config
from pyirena.batch.sizes import fit_sizes
from pyirena.batch.unified import fit_unified
from pyirena.batch.waxs import fit_waxs


def fit_pyirena(
    data_file: Union[str, Path],
    config_file: Union[str, Path],
    save_to_nexus: bool = True,
    with_uncertainty: bool = False,
    n_mc_runs: int = 10,
    tools: Optional[List[str]] = None,
) -> Optional[Dict]:
    """Run all analysis tools that have a configuration group in the config file.

    This is the top-level batch entry point.  It reads the config file, detects
    which tool sections are present, and dispatches to the appropriate fitting
    function for each tool.  Currently supported tools:

    ``unified_fit``
        Runs :func:`fit_unified`.
    ``sizes``
        Runs :func:`fit_sizes`.
    ``simple_fits``
        Runs :func:`fit_simple_from_config`.
    ``waxs_peakfit``
        Runs :func:`fit_waxs` (alias for :func:`fit_waxs_peaks_from_config`).
    ``modeling``
        Runs :func:`fit_modeling`.

    Unknown sections in the config file are silently skipped.

    Parameters
    ----------
    data_file : str or Path
        Path to SAS data file.
    config_file : str or Path
        Path to a pyIrena JSON configuration file.
    save_to_nexus : bool, optional
        Passed to each individual tool's fitting function (default True).
    with_uncertainty : bool, optional
        Passed to each individual tool's fitting function (default False).
    n_mc_runs : int, optional
        Passed to each individual tool's fitting function (default 10).
    tools : list of str or None, optional
        If given, only the named tools are run (e.g. ``['unified_fit', 'sizes']``).
        Tools not in the list are skipped even if present in the config file.
        When None (default), every recognised tool section in the config is run —
        this is identical to the behaviour before this parameter was added, so
        existing callers are unaffected.

    Returns
    -------
    dict or None
        On success, a dictionary with:

        ``'input_file'``   Path
        ``'config_file'``  Path
        ``'tools_run'``    list of str — tool names that were executed.
        ``'results'``      dict mapping tool name → result dict (or None on failure).

        Returns None if the config file cannot be loaded or contains no known tools.

    Examples
    --------
    >>> # Run all tools present in the config (original behaviour)
    >>> results = fit_pyirena("sample.h5", "pyirena_config.json")

    >>> # Run only unified_fit, skip everything else in the config
    >>> results = fit_pyirena("sample.h5", "pyirena_config.json",
    ...                       tools=["unified_fit"])

    >>> if results:
    ...     uf = results['results'].get('unified_fit')
    ...     if uf and uf['success']:
    ...         print(uf['parameters']['chi_squared'])
    """
    _ensure_console()
    data_file = Path(data_file)
    config_file = Path(config_file)

    # Registry: config key → callable
    _TOOL_REGISTRY: Dict[str, callable] = {
        'unified_fit': lambda: fit_unified(
            data_file, config_file, save_to_nexus, with_uncertainty, n_mc_runs
        ),
        'sizes': lambda: fit_sizes(
            data_file, config_file, save_to_nexus, with_uncertainty, n_mc_runs
        ),
        'simple_fits': lambda: fit_simple_from_config(
            data_file, config_file, save_to_nexus, with_uncertainty, n_mc_runs
        ),
        'waxs_peakfit': lambda: fit_waxs(
            data_file, config_file, save_to_nexus, with_uncertainty, n_mc_runs
        ),
        'modeling': lambda: fit_modeling(
            data_file, config_file, save_to_nexus, with_uncertainty, n_mc_runs
        ),
    }

    # --- Load config to discover which tools are present ---
    try:
        config = _load_config(config_file)
        if config is None:
            return None
    except Exception:
        log.error(f"[pyirena.batch] Unexpected error reading config:\n{traceback.format_exc()}")
        return None

    tools_to_run = [key for key in _TOOL_REGISTRY if key in config]
    if tools is not None:
        tools_to_run = [t for t in tools_to_run if t in tools]

    if not tools_to_run:
        log.info(f"[pyirena.batch] Config file '{config_file}' contains no recognised "
              f"tool sections. Known tools: {list(_TOOL_REGISTRY)}")
        return None

    # --- Run each tool ---
    all_results: Dict[str, Optional[Dict]] = {}
    for tool in tools_to_run:
        log.info(f"[pyirena.batch] Running '{tool}' on '{data_file.name}' ...")
        try:
            all_results[tool] = _TOOL_REGISTRY[tool]()
        except Exception:
            log.error(f"[pyirena.batch] Unhandled error in '{tool}':\n{traceback.format_exc()}")
            all_results[tool] = None

    return {
        'input_file': data_file,
        'config_file': config_file,
        'tools_run': tools_to_run,
        'results': all_results,
    }
