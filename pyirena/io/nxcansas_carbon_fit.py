"""
HDF5 save / load for Carbon model results (NXcanSAS format extension).

Group path inside the NXcanSAS file::

    entry/carbon_fit_results   (NXprocess)

Structure
---------
::

    carbon_fit_results/
        (attrs)  NX_class, program, timestamp, schema_version, success,
                 message, saxs_mode, waxs_envelope, n_peaks, formula,
                 q_min, q_max, _pyirena_config (JSON setup envelope)
        chi_squared, reduced_chi_squared, n_points, n_params — float64 scalars
        Q                — 1-D float64, Å⁻¹
        intensity_data   — 1-D float64, cm⁻¹
        intensity_error  — 1-D float64, cm⁻¹ (optional)
        I_model          — 1-D float64, the total fit
        I_porod          — 1-D float64, grain Porod term + flat background
        I_mp             — 1-D float64, micropore term
        I_waxs           — 1-D float64, diffraction term
        residuals        — 1-D float64, (I_data − I_model)/σ
        params/          — one float64 scalar per fitted quantity
            (attrs on each) param_key, label, units, fit, limit_low, limit_high
        params_std/      — 1-σ uncertainty, same names
        derived/         — one float64 scalar per derived quantity
            (attrs on each) units

Dataset naming
--------------
The model's parameter keys are dotted (``background.S_macro``,
``peak.002.Q0``) because that is what makes them stable across a changing peak
list.  Inside HDF5 the dots become underscores — ``background_S_macro``,
``peak_002_Q0`` — so the names match the flat style the derived group and the
Data Explorer's trend plots already use.  The original dotted key is kept as
the ``param_key`` attribute, so the mapping is not something a reader has to
guess.

Every parameter lives in one flat ``params/`` group rather than in per-section
sub-groups.  A Carbon model fit has a handful of scalars per section and the
peak list is variable-length; a flat namespace keyed by a stable name means a
trend plot across files does not break when a peak is added, which a
positional ``peak_01/`` scheme would.
"""

from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path
from typing import Dict, Optional

import numpy as np

try:
    import h5py
except ImportError:                                    # pragma: no cover
    h5py = None                                        # type: ignore

log = logging.getLogger(__name__)

_GROUP = "entry/carbon_fit_results"
_PROGRAM = "pyIrena carbon_fit"
_TOOL = "carbon_fit"

#: Units for the derived quantities, by key prefix/name.  Written as an attr so
#: the Data Explorer, the report and any external reader all label the number
#: the same way without a second table.
_DERIVED_UNITS: Dict[str, str] = {
    'd002': 'angstrom', 'd100': 'angstrom',
    'rho_struc': 'g/cm^3', 'rho_sample': 'g/cm^3', 'porosity': '',
    'sld_struc': '10^10 cm^-2', 'sld_sample': '10^10 cm^-2',
    'contrast_porod': '10^20 cm^-4', 'contrast_micropore': '10^20 cm^-4',
    'S_macro': 'cm^2/cm^3', 'S_rough': 'cm^2/cm^3', 'S_part': 'cm^2/cm^3',
    'S_macro_m2_g': 'm^2/g', 'S_rough_m2_g': 'm^2/g', 'S_part_m2_g': 'm^2/g',
    'mp_phi': '', 'mp_I0': '1/cm', 'mp_radius': 'angstrom',
    'S_mp': 'cm^2/cm^3', 'S_mp_m2_g': 'm^2/g',
    'ts_xi': 'angstrom', 'ts_d': 'angstrom', 'ts_fa': '',
    'w_pore': 'angstrom', 'w_carbon': 'angstrom', 'ts_r_spheroid': 'angstrom',
    'delta_z2': 'angstrom^2',
    'crumple_D': '', 'crumple_sigma': 'angstrom', 'crumple_R': 'angstrom',
    'L_c': 'angstrom', 'L_a': 'angstrom', 'N_layers': '',
}


def _derived_units(key: str) -> str:
    """Unit string for a derived key, including the per-peak ones."""
    if key in _DERIVED_UNITS:
        return _DERIVED_UNITS[key]
    if key.startswith('peak_'):
        if key.endswith('_Q0') or key.endswith('_FWHM'):
            return '1/angstrom'
        if key.endswith('_d') or key.endswith('_L'):
            return 'angstrom'
        if key.endswith('_height'):
            return '1/cm'
    return ''


def _h5_name(param_key: str) -> str:
    """Dotted model key → flat HDF5 dataset name."""
    return param_key.replace('.', '_')


# ===========================================================================
# Save
# ===========================================================================

def save_carbon_fit_results(
    filepath: Path,
    result,
    model=None,
    setup_state: Optional[Dict] = None,
) -> None:
    """Write one Carbon model fit into an NXcanSAS HDF5 file.

    Args:
        filepath: HDF5 file to write into; created if missing, and any previous
            ``entry/carbon_fit_results`` group is replaced rather than merged,
            so a re-fit never leaves half of the old answer behind.
        result: A :class:`~pyirena.core.carbon_fit.CarbonFitResult`.
        model: The :class:`~pyirena.core.carbon_fit.CarbonFitModel` that
            produced it.  Used for the parameter metadata (labels, units,
            bounds, Fit? flags) and, when ``setup_state`` is not supplied, for
            the embedded setup config.
        setup_state: The panel's state dict.  Embedded as the
            ``_pyirena_config`` attribute so "Load Setup from File…" can put
            every control back — including the ones the fit did not touch.
            Defaults to ``model.to_dict()``, which is the headless equivalent.

    Raises:
        ImportError: if h5py is not installed.
    """
    if h5py is None:
        raise ImportError("h5py is required to save results.")

    filepath = Path(filepath)
    with h5py.File(filepath, "a") as f:
        if _GROUP in f:
            del f[_GROUP]
        grp = f.require_group(_GROUP)

        grp.attrs["NX_class"] = "NXprocess"
        grp.attrs["program"] = _PROGRAM
        grp.attrs["timestamp"] = str(result.timestamp
                                     or datetime.now().isoformat(timespec="seconds"))
        grp.attrs["success"] = bool(result.success)
        grp.attrs["message"] = str(result.message)
        if model is not None:
            grp.attrs["schema_version"] = int(model.SCHEMA_VERSION)
            grp.attrs["saxs_mode"] = str(model.saxs.mode)
            grp.attrs["waxs_envelope"] = str(model.waxs.envelope)
            grp.attrs["n_peaks"] = int(sum(1 for p in model.peaks if p.enabled))
            grp.attrs["formula"] = str(model.material.formula)
            grp.attrs["q_min"] = float(model.q_min)
            grp.attrs["q_max"] = float(model.q_max)

        for name, value in (("chi_squared", result.chi_squared),
                            ("reduced_chi_squared", result.reduced_chi_squared),
                            ("n_points", result.n_points),
                            ("n_params", result.n_params)):
            grp.create_dataset(name, data=float(value), dtype="float64")

        for name, arr, units in (
            ("Q", result.q, "1/angstrom"),
            ("intensity_data", result.I_data, "1/cm"),
            ("intensity_error", result.I_error, "1/cm"),
            ("I_model", result.I_model, "1/cm"),
            ("I_porod", result.I_porod, "1/cm"),
            ("I_mp", result.I_mp, "1/cm"),
            ("I_waxs", result.I_waxs, "1/cm"),
            ("residuals", result.residuals, "dimensionless"),
        ):
            if arr is None:
                continue
            ds = grp.create_dataset(name, data=np.asarray(arr, dtype=float),
                                    dtype="float64")
            ds.attrs["units"] = units

        # ── Parameters, with everything a reader needs to interpret them ──
        meta = {}
        if model is not None:
            meta = {r.key: r for r in model.parameter_refs(active_only=False)}

        p_grp = grp.create_group("params")
        s_grp = grp.create_group("params_std")
        for key, value in (result.params or {}).items():
            ds = p_grp.create_dataset(_h5_name(key), data=float(value),
                                      dtype="float64")
            ds.attrs["param_key"] = key
            ref = meta.get(key)
            if ref is not None:
                lo, hi = ref.limits
                ds.attrs["label"] = ref.label
                ds.attrs["units"] = ref.unit
                ds.attrs["fit"] = bool(ref.fit)
                ds.attrs["limit_low"] = float(lo)
                ds.attrs["limit_high"] = float(hi)
            std = (result.errors or {}).get(key)
            sd = s_grp.create_dataset(_h5_name(key),
                                      data=float(std if std is not None else np.nan),
                                      dtype="float64")
            sd.attrs["param_key"] = key

        d_grp = grp.create_group("derived")
        for key, value in (result.derived or {}).items():
            ds = d_grp.create_dataset(key, data=float(value), dtype="float64")
            ds.attrs["units"] = _derived_units(key)

        state = setup_state
        if state is None and model is not None:
            state = model.to_dict()
        if state is not None:
            try:
                from pyirena.io.setup_config import write_setup_config
                write_setup_config(grp, _TOOL, state)
            except Exception:
                log.warning("carbon_fit: could not embed setup config",
                            exc_info=True)


# ===========================================================================
# Load
# ===========================================================================

def load_carbon_fit_results(filepath: Path) -> Dict:
    """Read a Carbon model result group back into a plain dict.

    Args:
        filepath: HDF5/NXcanSAS file written by :func:`save_carbon_fit_results`.

    Returns:
        dict with the group's attributes, the arrays (``Q``, ``I_model``,
        ``I_porod``, ``I_mp``, ``I_waxs``, ``intensity_data``,
        ``intensity_error``, ``residuals``) as numpy arrays or ``None``, and
        ``params`` / ``params_std`` / ``derived`` as ``{key: float}`` maps.
        Parameter keys come back **dotted**, the way the model names them, so a
        caller can feed them straight back without knowing about the
        underscore spelling on disk.

    Raises:
        ImportError: if h5py is not installed.
        KeyError: if the file has no Carbon model results.
    """
    if h5py is None:
        raise ImportError("h5py is required to load results.")

    filepath = Path(filepath)
    with h5py.File(filepath, "r") as f:
        if _GROUP not in f:
            raise KeyError(f"No carbon_fit_results group in {filepath}")
        grp = f[_GROUP]
        attrs = dict(grp.attrs)

        def _arr(name):
            return np.array(grp[name], dtype=float) if name in grp else None

        def _num(name, default=np.nan):
            if name in grp and isinstance(grp[name], h5py.Dataset):
                try:
                    return float(grp[name][()])
                except Exception:
                    log.debug("carbon_fit: %s unreadable", name, exc_info=True)
            return float(attrs.get(name, default))

        def _scalars(sub, dotted: bool):
            out: Dict[str, float] = {}
            if sub not in grp:
                return out
            for name, ds in grp[sub].items():
                try:
                    value = float(ds[()])
                except Exception:
                    continue
                key = str(ds.attrs.get("param_key", name)) if dotted else name
                out[key] = value
            return out

        return {
            "success": bool(attrs.get("success", False)),
            "message": str(attrs.get("message", "")),
            "timestamp": str(attrs.get("timestamp", "")),
            "schema_version": int(attrs.get("schema_version", 1)),
            "saxs_mode": str(attrs.get("saxs_mode", "fractal")),
            "waxs_envelope": str(attrs.get("waxs_envelope", "none")),
            "n_peaks": int(attrs.get("n_peaks", 0)),
            "formula": str(attrs.get("formula", "C")),
            "q_min": float(attrs.get("q_min", np.nan)),
            "q_max": float(attrs.get("q_max", np.nan)),
            "chi_squared": _num("chi_squared"),
            "reduced_chi_squared": _num("reduced_chi_squared"),
            "n_points": int(_num("n_points", 0)),
            "n_params": int(_num("n_params", 0)),
            "Q": _arr("Q"),
            "intensity_data": _arr("intensity_data"),
            "intensity_error": _arr("intensity_error"),
            "I_model": _arr("I_model"),
            "I_porod": _arr("I_porod"),
            "I_mp": _arr("I_mp"),
            "I_waxs": _arr("I_waxs"),
            "residuals": _arr("residuals"),
            "params": _scalars("params", dotted=True),
            "params_std": _scalars("params_std", dotted=True),
            "derived": _scalars("derived", dotted=False),
        }


def load_carbon_fit_model(filepath: Path):
    """Rebuild the :class:`CarbonFitModel` stored alongside a saved fit.

    Reads the embedded ``_pyirena_config`` setup, which is the model's own
    ``to_dict()`` (or the panel state that wraps it), so a saved file can be
    re-opened, re-plotted and re-fitted without the user re-entering anything.

    Returns:
        A ``CarbonFitModel``, or ``None`` if the file carries no setup.
    """
    from pyirena.core.carbon_fit import CarbonFitModel
    from pyirena.io.setup_config import read_setup_config

    state = read_setup_config(filepath, _GROUP, _TOOL)
    if state is None:
        return None
    # The panel nests the model under 'model'; a headless save stores it flat.
    return CarbonFitModel.from_dict(state.get('model', state))
