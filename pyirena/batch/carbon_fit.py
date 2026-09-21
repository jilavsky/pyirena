"""
pyirena.batch.carbon_fit — headless Carbon model fitting.

Runs exactly the same :class:`~pyirena.core.carbon_fit.CarbonFitModel` the GUI
panel drives, from a JSON config section instead of a panel.  Because the whole
model — sections, peak list, links, fit flags and bounds — is one
``to_dict()``, the config section *is* the model and there is no second
translation layer to keep in step.

The ``carbon_fit`` section of a pyIrena config file is therefore whatever
``CarbonFitModel.to_dict()`` produced, optionally wrapped by the panel's own
state (the panel nests it under ``"model"``).  Both shapes load.

Example:
    >>> from pyirena.batch import fit_carbon                     # doctest: +SKIP
    >>> res = fit_carbon("carbon_powder.h5", "pyirena_config.json")
    >>> if res['success']:                                       # doctest: +SKIP
    ...     print(res['derived']['S_part_m2_g'], "m2/g")
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, Optional, Union

import numpy as np

from pyirena.batch._common import _load_config
from pyirena.logging_setup import ensure_console_output as _ensure_console

log = logging.getLogger(__name__)


def fit_carbon_from_config(
    data_file: Union[str, Path],
    config_file: Union[str, Path],
    save_to_nexus: bool = True,
    with_uncertainty: bool = False,
    n_mc_runs: int = 10,
) -> Optional[Dict]:
    """Fit the Carbon model using the ``carbon_fit`` section of a config file.

    Args:
        data_file: SAS data file — text (.dat/.txt) or NXcanSAS HDF5.
        config_file: pyIrena JSON config with a ``carbon_fit`` section, written
            by the panel's "Save params to JSON" button or by hand.
        save_to_nexus: Write results into ``entry/carbon_fit_results`` in the
            HDF5 file.
        with_uncertainty: Run the Monte-Carlo uncertainty pass.
        n_mc_runs: Number of Monte-Carlo passes when ``with_uncertainty``.

    Returns:
        Result dict, always with ``success`` and ``message``.  On success it
        also carries ``params``, ``errors``, ``derived``, ``chi_squared``,
        ``reduced_chi_squared``, ``n_points``, ``n_params`` and the arrays
        ``q``, ``I_data``, ``I_model``, ``I_porod``, ``I_mp``, ``I_waxs``.
    """
    _ensure_console()
    config_file = Path(config_file)
    config = _load_config(config_file)
    if config is None:
        return {'success': False, 'message': f"Cannot load config: {config_file}"}

    section = config.get('carbon_fit')
    if section is None:
        msg = f"No 'carbon_fit' section in '{config_file.name}'"
        log.info("[pyirena.batch.fit_carbon] %s", msg)
        return {'success': False, 'message': msg}

    return fit_carbon_model(
        data_file=data_file,
        config=section,
        save_to_nexus=save_to_nexus,
        with_uncertainty=with_uncertainty,
        n_mc_runs=n_mc_runs,
        setup_state=section,
    )


#: Short alias, matching ``fit_unified`` / ``fit_sizes`` / ``fit_waxs``.
fit_carbon = fit_carbon_from_config


def fit_carbon_model(
    data_file: Union[str, Path],
    config: Dict,
    save_to_nexus: bool = True,
    with_uncertainty: bool = False,
    n_mc_runs: int = 10,
    verbose: bool = True,
    setup_state: Optional[Dict] = None,
) -> Dict:
    """Fit one data file with a Carbon model config dict.

    Args:
        data_file: Text or NXcanSAS HDF5 data file.
        config: The model's ``to_dict()``, or a panel state dict containing it
            under ``"model"``.  Missing keys take their defaults, so a config
            written by an older pyIrena still runs.
        save_to_nexus: Save into ``entry/carbon_fit_results``.
        with_uncertainty: Run the Monte-Carlo uncertainty pass, overriding the
            config's own ``n_mc_runs``.
        n_mc_runs: Passes to use when ``with_uncertainty`` is set.
        verbose: Log progress.
        setup_state: State embedded as ``_pyirena_config`` for setup restore;
            defaults to the model's own dict.

    Returns:
        Result dict — see :func:`fit_carbon_from_config`.
    """
    _ensure_console()
    from pyirena.core.carbon_fit import CarbonFitModel

    data_file = Path(data_file)
    if verbose:
        log.info("[pyirena.batch.fit_carbon] %s", data_file.name)

    loaded = _read_data(data_file, verbose=verbose)
    if loaded is None:
        return {'success': False,
                'message': f"Could not load data from {data_file.name}"}
    q, I, err, h5_path = loaded

    model = CarbonFitModel.from_dict(config.get('model', config))
    if with_uncertainty:
        model.n_mc_runs = int(n_mc_runs)

    try:
        result = model.fit(q, I, err)
    except Exception as exc:
        log.error("[pyirena.batch.fit_carbon] fit failed: %s", exc)
        return {'success': False, 'message': f"Fit failed: {exc}"}

    if save_to_nexus and h5_path is not None:
        try:
            from pyirena.io.nxcansas_carbon_fit import save_carbon_fit_results
            save_carbon_fit_results(
                h5_path, result, model,
                setup_state=setup_state if setup_state is not None else model.to_dict())
            if verbose:
                log.info("[pyirena.batch.fit_carbon] saved to %s:entry/carbon_fit_results",
                         h5_path.name)
        except Exception as exc:
            log.error("[pyirena.batch.fit_carbon] could not save results: %s", exc)

    out = result.to_dict()
    out.update({
        'q': result.q, 'I_data': result.I_data, 'I_model': result.I_model,
        'I_porod': result.I_porod, 'I_mp': result.I_mp, 'I_waxs': result.I_waxs,
        'residuals': result.residuals,
        'model': model.to_dict(),
    })
    return out


def _read_data(data_file: Path, verbose: bool = True):
    """Load (q, I, error, hdf5_path) from a text or NXcanSAS file.

    Text files go through ``ensure_nxcansas_sibling`` first, exactly as every
    other batch module does, so results always have somewhere NXcanSAS-shaped
    to be written back to.  Returns ``None`` if the file cannot be read.
    """
    from pyirena.io.hdf5 import readGenericNXcanSAS

    try:
        if data_file.suffix.lower() in ('.txt', '.dat'):
            from pyirena.io.text_import import ensure_nxcansas_sibling
            from pyirena.state.state_manager import StateManager
            q_unit = StateManager().get('data_selector', 'q_unit', '1/A')
            h5_file = ensure_nxcansas_sibling(data_file, q_unit=q_unit)
        else:
            h5_file = data_file
        data = readGenericNXcanSAS(str(h5_file.parent), h5_file.name)
    except Exception as exc:
        if verbose:
            log.error("[pyirena.batch.fit_carbon] load error: %s", exc)
        return None
    if data is None:
        return None

    q = np.asarray(data['Q'], dtype=float)
    I = np.asarray(data['Intensity'], dtype=float)
    err = data.get('Error')
    err = np.asarray(err, dtype=float) if err is not None else None
    return q, I, err, Path(h5_file)
