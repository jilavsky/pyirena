"""
NXcanSAS I/O for Simple Fits results.

Saves and loads single-model fitting results (Guinier, Porod, Sphere, …) in
the same NXcanSAS HDF5 file that already contains experimental SAS data.
Results are stored under:

    entry/simple_fit_results   (NXprocess group)

Group attributes (fit metadata)
--------------------------------
model, success, chi_squared, reduced_chi_squared, dof,
q_min, q_max, use_complex_bg, timestamp, program

Sub-groups
----------
params/
    One scalar dataset per fitted parameter (e.g. I0, Rg, …).

params_std/
    Standard deviation for each parameter — from the covariance matrix
    immediately after fitting, or updated from Monte Carlo uncertainty runs.

derived/
    Model-specific computed quantities:
    - Guinier Sheet: Thickness  [Å]
    - Treubner-Strey: CorrLength [Å], RepeatDist [Å]
    - Unified Born Green: Rad [Å], G1, Rg2 [Å]
    - Invariant: Invariant [cm⁻⁴], VolumeFraction, PhiOneMinusPhi,
      QmaxUsed [Å⁻¹], and (optional) PorodTail [cm⁻⁴], PorodKp [cm⁻¹Å⁻⁴]

Datasets
--------
Q             — scattering vector used for I_model [Å⁻¹]
I_model       — fitted model intensity
residuals     — (I_data − I_model) / error  (stored if available)
intensity_data — measured intensity (stored if provided)
intensity_error — measurement uncertainty (stored if provided)

Extra arrays (calculation models — Invariant)
---------------------------------------------
Q_integral        — Q grid of the running invariant integral [Å⁻¹]
running_integral  — cumulative ∫q²·I_corr dq  [1/cm⁴]
I_corrected       — background-corrected intensity [1/cm]
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Optional

import h5py
import numpy as np


_GROUP = 'entry/simple_fit_results'
_PROGRAM = 'pyirena.core.simple_fits'


def save_simple_fit_results(
    filepath: Path,
    result: dict,
    model_obj=None,
    intensity_data: Optional[np.ndarray] = None,
    intensity_error: Optional[np.ndarray] = None,
    setup_state: Optional[dict] = None,
    fit_quality: Optional[dict] = None,
) -> None:
    """
    Save simple-model fitting results to an NXcanSAS HDF5 file.

    The file must already exist.  Any pre-existing ``simple_fit_results``
    group is deleted and recreated.

    Parameters
    ----------
    filepath : Path
        Path to the existing NXcanSAS HDF5 file.
    result : dict
        Return value of ``SimpleFitModel.fit()``.  Must contain:
        ``model``, ``success``, ``params``, ``params_std``, ``I_model``,
        ``q``, ``residuals``, ``chi2``, ``dof``, ``reduced_chi2``,
        ``derived``.
    model_obj : SimpleFitModel, optional
        Model object — used to persist ``use_complex_bg``, ``n_mc_runs``,
        and ``limits``.
    intensity_data : array, optional
        Measured I(Q) values to store alongside the fit.
    intensity_error : array, optional
        Measurement uncertainties [same units as intensity_data].
    setup_state : dict, optional
        Full GUI state dict (the ``simple_fits`` section from StateManager).
        Embedded as the ``_pyirena_config`` JSON attribute on
        ``simple_fit_results`` so the GUI can restore every control after an
        AI-driven run.
    """
    filepath = Path(filepath)
    timestamp = datetime.now().isoformat()

    mode = 'a' if filepath.exists() else 'w'
    with h5py.File(filepath, mode) as f:
        if _GROUP in f:
            del f[_GROUP]

        grp = f.require_group(_GROUP)
        grp.attrs['NX_class']         = 'NXprocess'
        grp.attrs['program']          = _PROGRAM
        grp.attrs['timestamp']        = timestamp
        grp.attrs['model']            = str(result.get('model', ''))
        grp.attrs['success']          = bool(result.get('success', False))
        if result.get('dof') is not None:
            grp.attrs['dof'] = int(result['dof'])
        # Fit quality stored as scalar datasets (browseable/collectable in HDF5 viewer)
        if result.get('chi2') is not None:
            grp.create_dataset('chi_squared', data=float(result['chi2']))
        if result.get('reduced_chi2') is not None:
            grp.create_dataset('reduced_chi_squared', data=float(result['reduced_chi2']))

        if model_obj is not None:
            grp.attrs['use_complex_bg'] = bool(model_obj.use_complex_bg)
            grp.attrs['n_mc_runs']      = int(model_obj.n_mc_runs)
            grp.attrs['invariant_porod_tail'] = bool(
                getattr(model_obj, 'invariant_porod_tail', False))

        if result.get('warning'):
            grp.attrs['warning'] = str(result['warning'])

        q = result.get('q')
        if q is not None:
            if q.size > 0:
                grp.attrs['q_min'] = float(q[q > 0].min()) if np.any(q > 0) else float(q.min())
                grp.attrs['q_max'] = float(q.max())

        # ── Arrays ────────────────────────────────────────────────────────────
        if q is not None:
            grp.create_dataset('Q', data=q.astype('f8'), compression='gzip')
            grp['Q'].attrs['units'] = '1/angstrom'

        I_model = result.get('I_model')
        if I_model is not None:
            grp.create_dataset('I_model', data=I_model.astype('f8'), compression='gzip')
            grp['I_model'].attrs['units'] = 'arb'

        # Slit-smearing provenance + ideal (pinhole) model curve.
        grp.attrs['slit_length'] = float(result.get('slit_length', 0.0) or 0.0)
        grp.attrs['data_is_slit_smeared'] = bool(result.get('data_is_slit_smeared', False))
        I_model_ideal = result.get('I_model_ideal')
        if I_model_ideal is not None:
            grp.create_dataset('I_model_ideal',
                               data=np.asarray(I_model_ideal, dtype='f8'),
                               compression='gzip')
            grp['I_model_ideal'].attrs['units'] = 'arb'

        residuals = result.get('residuals')
        if residuals is not None:
            grp.create_dataset('residuals', data=residuals.astype('f8'), compression='gzip')
            grp['residuals'].attrs['units'] = 'dimensionless'

        if intensity_data is not None:
            grp.create_dataset('intensity_data',
                               data=np.asarray(intensity_data, dtype='f8'),
                               compression='gzip')
            grp['intensity_data'].attrs['units'] = '1/cm'

        if intensity_error is not None:
            grp.create_dataset('intensity_error',
                               data=np.asarray(intensity_error, dtype='f8'),
                               compression='gzip')
            grp['intensity_error'].attrs['units'] = '1/cm'

        # Extra arrays from calculation models (e.g. Invariant running integral)
        _EXTRA_UNITS = {
            'Q_integral':       '1/angstrom',
            'running_integral': '1/cm^4',
            'I_corrected':      '1/cm',
        }
        for name, arr in (result.get('extra_arrays') or {}).items():
            if arr is None:
                continue
            ds = grp.create_dataset(name, data=np.asarray(arr, dtype='f8'),
                                    compression='gzip')
            if name in _EXTRA_UNITS:
                ds.attrs['units'] = _EXTRA_UNITS[name]

        # ── Parameters ────────────────────────────────────────────────────────
        params = result.get('params', {})
        if params:
            pgrp = grp.require_group('params')
            for name, val in params.items():
                if val is not None:
                    ds = pgrp.create_dataset(name, data=float(val))
                    # Store limits alongside if model_obj available
                    if model_obj is not None:
                        lo, hi = model_obj.limits.get(name, (None, None))
                        if lo is not None:
                            ds.attrs['limit_low'] = float(lo)
                        if hi is not None:
                            ds.attrs['limit_high'] = float(hi)

        params_std = result.get('params_std', {})
        if params_std:
            sgrp = grp.require_group('params_std')
            for name, val in params_std.items():
                if val is not None and np.isfinite(float(val)):
                    sgrp.create_dataset(name, data=float(val))

        # ── Derived quantities ─────────────────────────────────────────────────
        derived = result.get('derived', {})
        if derived:
            dgrp = grp.require_group('derived')
            for name, val in derived.items():
                if val is not None and np.isfinite(float(val)):
                    dgrp.create_dataset(name, data=float(val))

        # Robust fit-quality metrics under fit_quality/.  Computed from the
        # consistent fit-range (data, model, error) triple if not supplied.
        from pyirena.io.nxcansas_fit_quality import write_fit_quality
        # Calculation models (Invariant) have no least-squares fit: I_model is
        # just the background curve, so fit-quality metrics would be
        # meaningless.  They are recognisable by chi2 being None.
        _is_calculation = result.get('chi2') is None and result.get('success')
        if (fit_quality is None and not _is_calculation
                and q is not None and I_model is not None
                and intensity_data is not None
                and len(intensity_data) == len(I_model) == len(q)):
            from pyirena.core.fit_metrics import fit_quality_metrics
            n_params = len(result.get('params', {}) or {})
            fit_quality = fit_quality_metrics(
                q, np.asarray(intensity_data, float), I_model,
                np.asarray(intensity_error, float) if intensity_error is not None else None,
                n_params=max(1, n_params))
        write_fit_quality(grp, fit_quality)

        # Embed the full GUI setup so the panel can round-trip from this file.
        if setup_state is not None:
            from pyirena.io.setup_config import write_setup_config
            write_setup_config(grp, "simple_fits", setup_state)


def load_simple_fit_results(filepath: Path) -> dict:
    """
    Load simple-model fitting results from an NXcanSAS HDF5 file.

    Parameters
    ----------
    filepath : Path
        Path to the NXcanSAS HDF5 file.

    Returns
    -------
    dict with keys:
        ``model``, ``success``, ``chi_squared``, ``reduced_chi_squared``,
        ``dof``, ``q_min``, ``q_max``, ``use_complex_bg``, ``n_mc_runs``,
        ``timestamp``, ``program``,
        ``Q``, ``I_model``, ``residuals``,
        ``intensity_data`` (may be None), ``intensity_error`` (may be None),
        ``params`` (dict), ``params_std`` (dict), ``derived`` (dict).

    Raises
    ------
    KeyError
        If the file does not contain a ``simple_fit_results`` group.
    OSError
        If the file cannot be opened.
    """
    filepath = Path(filepath)
    with h5py.File(filepath, 'r') as f:
        if _GROUP not in f:
            raise KeyError(
                f"No simple_fit_results group found in {filepath.name}. "
                "Run the Simple Fits tool first."
            )
        grp = f[_GROUP]

        result: dict = {}

        # ── Scalar metadata ────────────────────────────────────────────────────
        for k in (
            'model', 'success', 'dof', 'q_min', 'q_max',
            'use_complex_bg', 'n_mc_runs', 'timestamp', 'program',
            'invariant_porod_tail', 'warning',
        ):
            result[k] = grp.attrs.get(k)
        # chi_squared / reduced_chi_squared: try dataset (new) then attr (old)
        for k in ('chi_squared', 'reduced_chi_squared'):
            if k in grp and isinstance(grp[k], h5py.Dataset):
                try:
                    result[k] = float(grp[k][()])
                except Exception:
                    result[k] = grp.attrs.get(k)
            else:
                result[k] = grp.attrs.get(k)

        # ── Arrays ────────────────────────────────────────────────────────────
        result['Q']       = grp['Q'][:]       if 'Q'       in grp else None
        result['I_model'] = grp['I_model'][:] if 'I_model' in grp else None
        result['residuals'] = grp['residuals'][:] if 'residuals' in grp else None
        result['intensity_data']  = grp['intensity_data'][:]  if 'intensity_data'  in grp else None
        result['intensity_error'] = grp['intensity_error'][:] if 'intensity_error' in grp else None

        # Extra arrays (Invariant running integral etc.)
        for name in ('Q_integral', 'running_integral', 'I_corrected'):
            result[name] = grp[name][:] if name in grp else None

        # Robust fit-quality metrics (None for files written before this feature)
        from pyirena.io.nxcansas_fit_quality import read_fit_quality
        result['fit_quality'] = read_fit_quality(grp)

        # ── Sub-groups ─────────────────────────────────────────────────────────
        result['params'] = {}
        if 'params' in grp:
            for name, ds in grp['params'].items():
                result['params'][name] = float(ds[()])

        result['params_std'] = {}
        if 'params_std' in grp:
            for name, ds in grp['params_std'].items():
                result['params_std'][name] = float(ds[()])

        result['derived'] = {}
        if 'derived' in grp:
            for name, ds in grp['derived'].items():
                result['derived'][name] = float(ds[()])

        # Slit-smearing provenance (absent in legacy files → pinhole defaults)
        result['slit_length'] = float(grp.attrs.get('slit_length', 0.0))
        result['data_is_slit_smeared'] = bool(grp.attrs.get('data_is_slit_smeared', False))

    return result


def print_simple_fit_results(result: dict) -> None:
    """Pretty-print simple fit results to the console."""
    print('=' * 60)
    print('  Simple Fits Results')
    print('=' * 60)
    print(f"  Model:            {result.get('model', 'N/A')}")
    print(f"  Success:          {result.get('success', 'N/A')}")
    chi2 = result.get('chi_squared')
    rchi2 = result.get('reduced_chi_squared')
    if chi2 is not None:
        print(f"  Chi-squared:      {chi2:.4g}")
    if rchi2 is not None:
        print(f"  Reduced chi²:     {rchi2:.4g}")
    params = result.get('params', {})
    std = result.get('params_std', {})
    if params:
        print('  Parameters:')
        for name, val in params.items():
            err = std.get(name)
            if err is not None and np.isfinite(err):
                print(f'    {name:20s} = {val:.6g}  ±  {err:.3g}')
            else:
                print(f'    {name:20s} = {val:.6g}')
    derived = result.get('derived', {})
    if derived:
        print('  Derived:')
        for name, val in derived.items():
            print(f'    {name:20s} = {val:.6g}')
    print(f"  Timestamp:        {result.get('timestamp', 'N/A')}")
    print('=' * 60)
