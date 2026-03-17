"""
NXcanSAS I/O for Modeling tool results.

Saves and loads ModelingResult objects to/from an HDF5 file group
``entry/modeling_results``.  The structure stores:

  - Global scalars (chi², background, q_min/max, timestamp, …)
  - Per-population sub-groups (pop_01/ … pop_10/) with all parameters,
    radius grids, volume/number distributions, and derived quantities.
  - Model arrays: model_q, model_I, plus per-population I(Q) curves.

Functions
---------
save_modeling_results(filepath, result)  → None
load_modeling_results(filepath)          → dict
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path
from typing import Optional

import h5py
import numpy as np

from pyirena.core.modeling import ModelingResult, ModelingConfig, SizeDistPopulation


# ──────────────────────────────────────────────────────────────────────────────
# Save
# ──────────────────────────────────────────────────────────────────────────────

def save_modeling_results(
    filepath: Path,
    result: ModelingResult,
    group_name: str = 'modeling_results',
) -> None:
    """Save Modeling fit results to an HDF5 file.

    If the file does not exist it is created with a minimal NXcanSAS shell.
    If a previous ``modeling_results`` group exists it is deleted and replaced.

    Args:
        filepath:   Path to HDF5 file (created if absent, opened 'a' if present).
        result:     ModelingResult from ModelingEngine.fit().
        group_name: HDF5 group name relative to 'entry/' (default 'modeling_results').
    """
    filepath = Path(filepath)
    timestamp = result.timestamp

    with h5py.File(filepath, 'a') as f:
        # ── Ensure top-level entry exists ─────────────────────────────────
        if 'entry' not in f:
            f.attrs.update({
                'default': 'entry',
                'file_name': filepath.name,
                'file_time': timestamp,
                'creator': 'pyirena',
                'NeXus_version': '4.3.0',
            })
            nxe = f.create_group('entry')
            nxe.attrs['NX_class'] = 'NXentry'
            nxe.attrs['canSAS_class'] = 'SASentry'
            nxe.create_dataset('definition', data='NXsas')

        full_path = f'entry/{group_name}'
        if full_path in f:
            del f[full_path]

        grp = f.create_group(full_path)
        grp.attrs['NX_class'] = 'NXprocess'
        grp.attrs['analysis_type'] = 'Modeling'
        grp.attrs['program'] = 'pyirena'
        grp.attrs['timestamp'] = timestamp

        # ── Global scalars ────────────────────────────────────────────────
        grp.create_dataset('chi_squared', data=float(result.chi_squared))
        grp.create_dataset('reduced_chi_squared', data=float(result.reduced_chi_squared))
        grp.create_dataset('dof', data=int(result.dof))
        grp.create_dataset('background', data=float(result.config.background))
        grp.create_dataset('q_min', data=float(result.config.q_min))
        grp.create_dataset('q_max', data=float(result.config.q_max))

        # ── Model arrays ──────────────────────────────────────────────────
        ds_q = grp.create_dataset('model_q', data=result.model_q)
        ds_q.attrs['units'] = '1/angstrom'
        ds_I = grp.create_dataset('model_I', data=result.model_I)
        ds_I.attrs['units'] = '1/cm'

        # ── Per-population groups ─────────────────────────────────────────
        cfg = result.config
        for k, pi in enumerate(result.pop_indices):
            pop: SizeDistPopulation = cfg.populations[pi]
            pname = f'pop_{pi+1:02d}'
            pg = grp.create_group(pname)
            pg.attrs['population_index'] = pi + 1
            pg.attrs['enabled'] = True
            pg.attrs['dist_type'] = pop.dist_type
            pg.attrs['form_factor'] = pop.form_factor
            pg.attrs['structure_factor'] = pop.structure_factor

            # Distribution parameters
            dg = pg.create_group('dist_params')
            for pn, pv in pop.dist_params.items():
                dg.create_dataset(pn, data=float(pv))
                dg[pn].attrs['fit'] = bool(pop.dist_params_fit.get(pn, False))
                lim = pop.dist_params_limits.get(pn, (0.0, 1e10))
                dg[pn].attrs['limit_lo'] = float(lim[0])
                dg[pn].attrs['limit_hi'] = float(lim[1])

            # Form-factor parameters
            fg = pg.create_group('ff_params')
            for pn, pv in pop.ff_params.items():
                fg.create_dataset(pn, data=float(pv))
                fg[pn].attrs['fit'] = bool(pop.ff_params_fit.get(pn, False))

            # Structure-factor parameters
            sg = pg.create_group('sf_params')
            for pn, pv in pop.sf_params.items():
                sg.create_dataset(pn, data=float(pv))
                sg[pn].attrs['fit'] = bool(pop.sf_params_fit.get(pn, False))

            # Scale / contrast
            pg.create_dataset('scale', data=float(pop.scale))
            pg['scale'].attrs['fit'] = bool(pop.fit_scale)
            pg.create_dataset('contrast', data=float(pop.contrast))
            pg['contrast'].attrs['fit'] = bool(pop.fit_contrast)
            pg.create_dataset('use_number_dist', data=bool(pop.use_number_dist))
            pg.create_dataset('n_bins', data=int(pop.n_bins))

            # MC uncertainties (if present)
            stds = result.params_std
            for label, stdval in stds.items():
                if label.startswith(f'pop{pi+1}_'):
                    pg.create_dataset(label + '_err', data=float(stdval))

            # Model I(Q) for this population
            ds_pI = pg.create_dataset('model_I', data=result.pop_model_I[k])
            ds_pI.attrs['units'] = '1/cm'

            # Radius grid and distributions
            pg.create_dataset('radius_grid', data=result.radius_grids[k])
            pg['radius_grid'].attrs['units'] = 'angstrom'
            pg.create_dataset('volume_dist', data=result.volume_dists[k])
            pg.create_dataset('number_dist', data=result.number_dists[k])

            # Derived quantities
            if k < len(result.derived):
                dv = result.derived[k]
                dervg = pg.create_group('derived')
                for dn, dval in dv.items():
                    dervg.create_dataset(dn, data=float(dval))

        # ── MC uncertainties for background ───────────────────────────────
        bg_err = result.params_std.get('background', 0.0)
        if bg_err > 0.0:
            grp.create_dataset('background_err', data=float(bg_err))

    print(f"Saved Modeling results to {filepath}")


# ──────────────────────────────────────────────────────────────────────────────
# Load
# ──────────────────────────────────────────────────────────────────────────────

def load_modeling_results(
    filepath: Path,
    group_name: str = 'modeling_results',
) -> Optional[dict]:
    """Load Modeling results from an HDF5 file.

    Returns a dict with the same structure as saved by save_modeling_results(),
    or None if no modeling results group is found.

    Returned dict keys (partial):
        chi_squared, reduced_chi_squared, dof, background, q_min, q_max,
        timestamp, model_q, model_I,
        populations: list of per-population dicts with sub-keys:
            dist_type, dist_params, form_factor, ff_params, structure_factor,
            sf_params, scale, contrast, use_number_dist, n_bins,
            model_I, radius_grid, volume_dist, number_dist, derived.
    """
    filepath = Path(filepath)
    if not filepath.exists():
        return None

    with h5py.File(filepath, 'r') as f:
        full_path = f'entry/{group_name}'
        if full_path not in f:
            return None

        grp = f[full_path]

        def _scalar(name, default=None):
            if name in grp:
                v = grp[name][()]
                return float(v) if np.isscalar(v) else v
            return default

        result = {
            'chi_squared':         _scalar('chi_squared'),
            'reduced_chi_squared': _scalar('reduced_chi_squared'),
            'dof':                 int(_scalar('dof', 0)),
            'background':          _scalar('background', 0.0),
            'q_min':               _scalar('q_min'),
            'q_max':               _scalar('q_max'),
            'timestamp':           grp.attrs.get('timestamp', ''),
            'model_q':             grp['model_q'][()] if 'model_q' in grp else None,
            'model_I':             grp['model_I'][()] if 'model_I' in grp else None,
            'populations':         [],
        }

        # ── Per-population groups ──────────────────────────────────────────
        for key in sorted(grp.keys()):
            if not key.startswith('pop_'):
                continue
            pg = grp[key]

            def _arr(subgrp, name):
                return subgrp[name][()] if name in subgrp else None

            def _load_param_group(subgrp):
                d = {}
                if subgrp is None:
                    return d
                for pn in subgrp:
                    d[pn] = float(subgrp[pn][()])
                return d

            pop_dict = {
                'population_index': int(pg.attrs.get('population_index', 0)),
                'enabled':          bool(pg.attrs.get('enabled', True)),
                'dist_type':        pg.attrs.get('dist_type', 'lognormal'),
                'form_factor':      pg.attrs.get('form_factor', 'sphere'),
                'structure_factor': pg.attrs.get('structure_factor', 'none'),
                'dist_params':      _load_param_group(pg.get('dist_params')),
                'ff_params':        _load_param_group(pg.get('ff_params')),
                'sf_params':        _load_param_group(pg.get('sf_params')),
                'scale':            float(pg['scale'][()]) if 'scale' in pg else 0.001,
                'contrast':         float(pg['contrast'][()]) if 'contrast' in pg else 1.0,
                'use_number_dist':  bool(pg['use_number_dist'][()]) if 'use_number_dist' in pg else False,
                'n_bins':           int(pg['n_bins'][()]) if 'n_bins' in pg else 200,
                'model_I':          _arr(pg, 'model_I'),
                'radius_grid':      _arr(pg, 'radius_grid'),
                'volume_dist':      _arr(pg, 'volume_dist'),
                'number_dist':      _arr(pg, 'number_dist'),
                'derived':          {},
            }

            dervg = pg.get('derived')
            if dervg is not None:
                for dn in dervg:
                    pop_dict['derived'][dn] = float(dervg[dn][()])

            result['populations'].append(pop_dict)

    return result
