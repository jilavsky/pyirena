"""
NXcanSAS I/O for SAXS Morph tool results.

Saves and loads SaxsMorphResult objects to/from an HDF5 file group
``entry/saxs_morph_results``.

Layout
------
entry/saxs_morph_results/                 (NXprocess)
  @NX_class            = "NXprocess"
  @analysis_type       = "SAXS Morph"
  @program             = "pyirena"
  @timestamp           (ISO date)
  chi_squared, reduced_chi_squared, dof    (scalars)
  q_min, q_max                             (scalars)
  volume_fraction, contrast, link_phi_contrast (scalars)
  power_law_B, power_law_P, background    (scalars)
  rng_seed, voxel_size, box_size_A, voxel_pitch_A, phi_actual  (scalars)
  data_q, data_I, data_dI, data_I_corr   (1-D arrays)
  model_q, model_I                       (1-D arrays)
  r_grid, gamma_r                        (1-D arrays)
  spectral_k, spectral_F                 (1-D arrays)
  voxelgram                              (3-D uint8, gzip-compressed,
                                          chunks=(N, N, 1) for cheap slice loads)
  <param>_err                            (scalars, only when MC was run)
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Optional

import h5py
import numpy as np

from pyirena.core.saxs_morph import (
    SaxsMorphConfig, SaxsMorphResult,
)


GROUP_NAME = 'saxs_morph_results'


# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------

def save_saxs_morph_results(
    filepath: Path,
    result: SaxsMorphResult,
    group_name: str = GROUP_NAME,
) -> None:
    """Save a SaxsMorphResult to ``entry/<group_name>``.

    The file is created with a minimal NXcanSAS shell if it does not exist;
    an existing group of the same name is replaced.
    """
    filepath = Path(filepath)
    timestamp = result.timestamp or datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    with h5py.File(filepath, 'a') as f:
        # ── Top-level entry ───────────────────────────────────────────────
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

        full = f'entry/{group_name}'
        if full in f:
            del f[full]

        grp = f.create_group(full)
        grp.attrs['NX_class'] = 'NXprocess'
        grp.attrs['analysis_type'] = 'SAXS Morph'
        grp.attrs['program'] = 'pyirena'
        grp.attrs['timestamp'] = timestamp

        cfg = result.config

        # ── Scalars ──────────────────────────────────────────────────────
        for name, val in [
            ('chi_squared', result.chi_squared),
            ('reduced_chi_squared', result.reduced_chi_squared),
            ('dof', result.dof),
            ('q_min', cfg.q_min),
            ('q_max', cfg.q_max),
            ('volume_fraction', cfg.volume_fraction),
            ('contrast', cfg.contrast),
            ('power_law_B', cfg.power_law_B),
            ('power_law_P', cfg.power_law_P),
            ('background', cfg.background),
            ('voxel_size', result.voxel_size),
            ('box_size_A', result.box_size_A),
            ('voxel_pitch_A', result.voxel_pitch_A),
            ('phi_actual', result.phi_actual),
            ('rng_seed', result.rng_seed_used),
            ('rg_A', getattr(result, 'rg_A', float('nan'))),
        ]:
            grp.create_dataset(name, data=_h5_scalar(val))

        grp.create_dataset('link_phi_contrast', data=bool(cfg.link_phi_contrast))

        # Fit-flag attrs alongside each fittable param
        for name in ('volume_fraction', 'contrast',
                     'power_law_B', 'power_law_P', 'background'):
            ds = grp[name]
            ds.attrs['fit'] = bool(getattr(cfg, f'fit_{name}'))
            lim = getattr(cfg, f'{name}_limits')
            ds.attrs['limit_lo'] = float(lim[0])
            ds.attrs['limit_hi'] = float(lim[1])

        # ── 1-D arrays ───────────────────────────────────────────────────
        _save_1d(grp, 'data_q', result.data_q, units='1/angstrom')
        _save_1d(grp, 'data_I', result.data_I, units='1/cm')
        _save_1d(grp, 'data_dI', result.data_dI, units='1/cm')
        _save_1d(grp, 'data_I_corr', result.data_I_corr, units='1/cm')
        _save_1d(grp, 'model_q', result.model_q, units='1/angstrom')
        _save_1d(grp, 'model_I', result.model_I, units='1/cm')
        _save_1d(grp, 'r_grid', result.r_grid, units='angstrom')
        _save_1d(grp, 'gamma_r', result.gamma_r, units='dimensionless')
        _save_1d(grp, 'spectral_k', result.spectral_k, units='1/angstrom')
        _save_1d(grp, 'spectral_F', result.spectral_F, units='angstrom**3')

        # ── Voxelgram (compressed) ───────────────────────────────────────
        N = result.voxelgram.shape[0]
        chunk_z = min(N, 1)
        ds_vox = grp.create_dataset(
            'voxelgram',
            data=np.ascontiguousarray(result.voxelgram, dtype=np.uint8),
            chunks=(N, N, chunk_z),
            compression='gzip', compression_opts=4,
        )
        ds_vox.attrs['pitch_A'] = float(result.voxel_pitch_A)
        ds_vox.attrs['box_size_A'] = float(result.box_size_A)
        ds_vox.attrs['phi_actual'] = float(result.phi_actual)
        ds_vox.attrs['description'] = (
            'Binary phase indicator: 0 = phase A (background), '
            '1 = phase B (scattering phase). Cube of side voxel_size.'
        )

        # ── MC uncertainties ─────────────────────────────────────────────
        for name, std in result.params_std.items():
            grp.create_dataset(f'{name}_err', data=float(std))

    print(f"Saved SAXS Morph results to {filepath}")


def _save_1d(grp, name, arr, units: str):
    if arr is None:
        return
    ds = grp.create_dataset(name, data=np.asarray(arr, dtype=float))
    ds.attrs['units'] = units


def _h5_scalar(v):
    """Coerce v into something h5py.create_dataset accepts as a scalar."""
    if v is None:
        return float('nan')
    if isinstance(v, bool):
        return bool(v)
    if isinstance(v, (int, np.integer)):
        return int(v)
    return float(v)


# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------

def load_saxs_morph_results(
    filepath: Path,
    group_name: str = GROUP_NAME,
) -> Optional[dict]:
    """Load SAXS Morph results from HDF5; returns None if absent.

    The returned dict has keys mirroring SaxsMorphResult fields plus the flat
    config scalars.  The 3-D voxelgram is loaded eagerly as a numpy.uint8
    array (16 MB at 256**3, 125 MB at 512**3 — callers that want lazy access
    should use h5py directly with the chunked dataset).
    """
    filepath = Path(filepath)
    if not filepath.exists():
        return None

    with h5py.File(filepath, 'r') as f:
        full = f'entry/{group_name}'
        if full not in f:
            return None
        grp = f[full]

        def _scal(name, default=None):
            if name not in grp:
                return default
            v = grp[name][()]
            if isinstance(v, (np.bool_, bool)):
                return bool(v)
            if hasattr(v, 'item'):
                return v.item()
            return v

        def _arr(name):
            return grp[name][()] if name in grp else None

        out = {
            # Administrative
            'timestamp':           grp.attrs.get('timestamp', ''),
            'program':             grp.attrs.get('program', ''),
            'analysis_type':       grp.attrs.get('analysis_type', ''),
            # Scalars
            'chi_squared':         _scal('chi_squared'),
            'reduced_chi_squared': _scal('reduced_chi_squared'),
            'dof':                 _scal('dof'),
            'q_min':               _scal('q_min'),
            'q_max':               _scal('q_max'),
            'volume_fraction':     _scal('volume_fraction'),
            'contrast':            _scal('contrast'),
            'power_law_B':         _scal('power_law_B'),
            'power_law_P':         _scal('power_law_P'),
            'background':          _scal('background'),
            'link_phi_contrast':   _scal('link_phi_contrast', False),
            'voxel_size':          _scal('voxel_size'),
            'box_size_A':          _scal('box_size_A'),
            'voxel_pitch_A':       _scal('voxel_pitch_A'),
            'phi_actual':          _scal('phi_actual'),
            'rng_seed':            _scal('rng_seed'),
            'rg_A':                _scal('rg_A', float('nan')),
            # Fit flags / limits (read attrs)
            'fit_flags':   {n: bool(grp[n].attrs.get('fit', False))
                            for n in ('volume_fraction', 'contrast',
                                      'power_law_B', 'power_law_P', 'background')
                            if n in grp},
            # 1-D arrays
            'data_q':       _arr('data_q'),
            'data_I':       _arr('data_I'),
            'data_dI':      _arr('data_dI'),
            'data_I_corr':  _arr('data_I_corr'),
            'model_q':      _arr('model_q'),
            'model_I':      _arr('model_I'),
            'r_grid':       _arr('r_grid'),
            'gamma_r':      _arr('gamma_r'),
            'spectral_k':   _arr('spectral_k'),
            'spectral_F':   _arr('spectral_F'),
            # Voxelgram (eager)
            'voxelgram':    _arr('voxelgram'),
            # MC uncertainties
            'params_std':   {n.replace('_err', ''): float(grp[n][()])
                             for n in grp
                             if n.endswith('_err') and n != 'reduced_chi_squared'},
        }
        return out


def result_from_loaded_dict(d: dict) -> SaxsMorphResult:
    """Reconstruct a SaxsMorphResult from a load_saxs_morph_results dict.

    Useful for the GUI to re-display an already-saved result without
    rerunning the engine.
    """
    cfg = SaxsMorphConfig(
        q_min=d.get('q_min'),
        q_max=d.get('q_max'),
        voxel_size_fit=int(d.get('voxel_size') or 128),
        voxel_size_render=int(d.get('voxel_size') or 256),
        box_size_A=float(d.get('box_size_A') or 1000.0),
        volume_fraction=float(d.get('volume_fraction') or 0.3),
        contrast=float(d.get('contrast') or 1.0),
        link_phi_contrast=bool(d.get('link_phi_contrast', False)),
        power_law_B=float(d.get('power_law_B') or 0.0),
        power_law_P=float(d.get('power_law_P') or 4.0),
        background=float(d.get('background') or 0.0),
    )
    flags = d.get('fit_flags') or {}
    for name, on in flags.items():
        setattr(cfg, f'fit_{name}', bool(on))

    return SaxsMorphResult(
        config=cfg,
        chi_squared=float(d.get('chi_squared') or 0.0),
        reduced_chi_squared=float(d.get('reduced_chi_squared') or 0.0),
        dof=int(d.get('dof') or 0),
        timestamp=str(d.get('timestamp') or ''),
        data_q=np.asarray(d.get('data_q')) if d.get('data_q') is not None else np.array([]),
        data_I=np.asarray(d.get('data_I')) if d.get('data_I') is not None else np.array([]),
        data_dI=np.asarray(d.get('data_dI')) if d.get('data_dI') is not None else np.array([]),
        data_I_corr=np.asarray(d.get('data_I_corr')) if d.get('data_I_corr') is not None else np.array([]),
        model_q=np.asarray(d.get('model_q')) if d.get('model_q') is not None else np.array([]),
        model_I=np.asarray(d.get('model_I')) if d.get('model_I') is not None else np.array([]),
        r_grid=np.asarray(d.get('r_grid')) if d.get('r_grid') is not None else np.array([]),
        gamma_r=np.asarray(d.get('gamma_r')) if d.get('gamma_r') is not None else np.array([]),
        spectral_k=np.asarray(d.get('spectral_k')) if d.get('spectral_k') is not None else np.array([]),
        spectral_F=np.asarray(d.get('spectral_F')) if d.get('spectral_F') is not None else np.array([]),
        voxelgram=np.asarray(d.get('voxelgram'), dtype=np.uint8) if d.get('voxelgram') is not None else np.zeros((0,0,0), dtype=np.uint8),
        voxel_size=int(d.get('voxel_size') or 0),
        box_size_A=float(d.get('box_size_A') or 0.0),
        voxel_pitch_A=float(d.get('voxel_pitch_A') or 0.0),
        phi_actual=float(d.get('phi_actual') or 0.0),
        rng_seed_used=int(d.get('rng_seed') or 0),
        rg_A=float(d.get('rg_A') if d.get('rg_A') is not None else float('nan')),
        params_std=dict(d.get('params_std') or {}),
    )
