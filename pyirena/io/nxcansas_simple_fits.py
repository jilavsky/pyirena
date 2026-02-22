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

Datasets
--------
Q             — scattering vector used for I_model [Å⁻¹]
I_model       — fitted model intensity
residuals     — (I_data − I_model) / error  (stored if available)
intensity_data — measured intensity (stored if provided)
intensity_error — measurement uncertainty (stored if provided)
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

        if result.get('chi2') is not None:
            grp.attrs['chi_squared'] = float(result['chi2'])
        if result.get('reduced_chi2') is not None:
            grp.attrs['reduced_chi_squared'] = float(result['reduced_chi2'])
        if result.get('dof') is not None:
            grp.attrs['dof'] = int(result['dof'])

        if model_obj is not None:
            grp.attrs['use_complex_bg'] = bool(model_obj.use_complex_bg)
            grp.attrs['n_mc_runs']      = int(model_obj.n_mc_runs)

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
            'model', 'success', 'chi_squared', 'reduced_chi_squared', 'dof',
            'q_min', 'q_max', 'use_complex_bg', 'n_mc_runs',
            'timestamp', 'program',
        ):
            result[k] = grp.attrs.get(k)

        # ── Arrays ────────────────────────────────────────────────────────────
        result['Q']       = grp['Q'][:]       if 'Q'       in grp else None
        result['I_model'] = grp['I_model'][:] if 'I_model' in grp else None
        result['residuals'] = grp['residuals'][:] if 'residuals' in grp else None
        result['intensity_data']  = grp['intensity_data'][:]  if 'intensity_data'  in grp else None
        result['intensity_error'] = grp['intensity_error'][:] if 'intensity_error' in grp else None

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
