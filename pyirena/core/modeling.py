"""
Modeling tool engine for pyIrena.

Computes small-angle scattering from a sum of up to 10 user-defined
"populations".  In Step 1 every population is a parametric size distribution
(Gauss, Log-Normal, LSW, Schulz-Zimm, or Ardell); future steps add Unified Fit
levels and diffraction peaks.

Each population:
  * Picks a bell-shaped distribution of particle radii (→ bin centres).
  * Multiplies a pre-built G-matrix by the volume-fraction distribution to get
    I_pop(Q).
  * Optionally applies a structure factor (None / Interferences / Hard-Sphere).

The total model is the sum over all enabled populations plus a flat background.

Fitting uses ``scipy.optimize.least_squares`` with bounds (or Nelder-Mead when
``no_limits=True``).  Monte-Carlo uncertainty is estimated by repeating the fit
on data perturbed by Gaussian noise ~ σ_I.

G-matrix caching
----------------
The G-matrix for a given (q-array, dist params, form-factor, contrast) is
expensive to rebuild but does *not* change during a "Graph Model" call.  We
cache it in ``ModelingEngine._g_cache`` keyed on a float-hash of
(q_bytes, dist_type, sorted(params.items()), form_factor, sorted(ff_params.items()), contrast).
The cache is intentionally *not* used during fitting because every optimizer
step changes the params.

References
----------
IR2L_Modeling*.ipf — Igor Pro source by Jan Ilavsky, Argonne National Laboratory.
IR2S_HardSphereStruct — Percus-Yevick hard-sphere structure factor (same file).
"""

from __future__ import annotations

import contextlib
import hashlib
import json
import logging
import os
import time
import warnings
from copy import deepcopy
from dataclasses import dataclass, field
from dataclasses import fields as dataclass_fields
from datetime import datetime
from typing import Optional

import numpy as np
from scipy.optimize import differential_evolution, least_squares, minimize
from scipy.special import erf as _scipy_erf

from pyirena.core.unified import sphere_amplitude as _sphere_amplitude

log = logging.getLogger(__name__)


class PopulationEvaluationError(RuntimeError):
    """A population's I(Q) could not be computed.

    Raised only on paths that must not silently continue — chiefly fitting,
    where dropping a population from the model changes chi-squared and the fit
    would otherwise report success for a model that is missing one of its
    components. Display paths keep warning and skipping so a half-typed
    parameter cannot take the panel down.
    """


# Weighted residual substituted for one failed trial evaluation during fitting.
# Large enough that the optimiser rejects the point outright, small enough that
# the sum of squares over a full Q range stays finite.
_FAILED_TRIAL_RESIDUAL = 1e8

# Thread-count environment variables read by the common BLAS / OpenMP backends.
# When running differential_evolution across worker processes we pin these to 1
# so N workers don't each spawn a full BLAS thread pool and oversubscribe the
# CPU. macOS Accelerate uses VECLIB_MAXIMUM_THREADS.
_BLAS_THREAD_ENV = (
    'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
    'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS',
)


@contextlib.contextmanager
def _limit_blas_threads(n: int = 1):
    """Temporarily pin BLAS/OpenMP thread counts (best-effort, restored after).

    Spawned worker processes (macOS/Windows 'spawn') re-import numpy and read
    these vars at import time, so setting them here before launching the DE
    worker pool keeps each worker single-threaded.
    """
    saved = {k: os.environ.get(k) for k in _BLAS_THREAD_ENV}
    try:
        for k in _BLAS_THREAD_ENV:
            os.environ[k] = str(n)
        yield
    finally:
        for k, v in saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v


class _DEObjective:
    """Picklable differential-evolution objective for parallel (workers>1) runs.

    A nested closure cannot be pickled to worker processes, and the GUI
    monkey-patches ``engine._chi2`` for cancellation (also unpicklable). This
    class carries only picklable state (a clean engine, the deep-copied config,
    the data arrays, and the log-sampling mask) and converts the internally
    log-scaled trial vector back to real units before evaluating χ². It must be
    module-level so ``pickle`` can reference it by qualified name.
    """

    def __init__(self, engine, keys, cfg, q, I, sigma, log_mask):
        self.engine = engine
        self.keys = keys
        self.cfg = cfg
        self.q = q
        self.I = I
        self.sigma = sigma
        self.log_mask = np.asarray(log_mask, dtype=bool)

    def __call__(self, t_x):
        x = np.array(t_x, dtype=float)
        if self.log_mask.any():
            x[self.log_mask] = 10.0 ** x[self.log_mask]
        return self.engine._chi2(x, self.keys, self.cfg, self.q, self.I, self.sigma)


# ── Monte-Carlo uncertainty: parallel machinery ───────────────────────────────
#
# One MC pass = perturb the data, refit, keep the fitted parameters. Passes are
# independent, so they distribute over worker processes with no coordination.
# A pool costs roughly a second to start (spawn re-imports numpy/scipy in every
# worker), so it only pays off once the remaining passes are slow enough; the
# first pass is timed in-process to decide. See ``_MC_PARALLEL_MIN_SECONDS``.

# Projected serial cost of the remaining passes, below which the pool is skipped.
_MC_PARALLEL_MIN_SECONDS = 5.0

# Process-local engine, created on first use inside each worker.
_MC_ENGINE: Optional["ModelingEngine"] = None


def _mc_engine() -> "ModelingEngine":
    """Return this process's reusable ``ModelingEngine``.

    One engine per worker process rather than one per pass, so the cached
    ``SlitSmearer`` (see :meth:`ModelingEngine._get_smearer`) is built once per
    process instead of rebuilt for every pass.

    The caller's live engine is deliberately never shipped across the process
    boundary: the Modeling GUI monkey-patches ``_chi2``/``_residuals`` with
    closures for "Cancel Fit", and closures cannot be pickled.
    """
    global _MC_ENGINE
    if _MC_ENGINE is None:
        _MC_ENGINE = ModelingEngine()
    return _MC_ENGINE


class _MCRun:
    """Picklable single Monte-Carlo pass, for parallel uncertainty estimation.

    Carries only picklable state (a config, the q grid, the uncertainties, and
    the expected parameter count) and is called with one pre-perturbed intensity
    array. Must be module-level so ``pickle`` can reference it by qualified name.

    Returns the flat list of fitted parameter values, or ``None`` when the pass
    failed — matching the serial path, which skips failed passes rather than
    aborting the whole estimate.
    """

    def __init__(self, cfg, q, dI, n_keys):
        self.cfg = cfg
        self.q = q
        self.dI = dI
        self.n_keys = n_keys

    def __call__(self, I_pert):
        engine = _mc_engine()
        # Every pass adds one entry to the q-keyed G cache (via the ideal-curve
        # evaluation in fit()). The engine outlives the pass, so drop them here
        # rather than growing without bound over a long run. This does not touch
        # the cached SlitSmearer, which is the whole point of reusing an engine.
        engine.clear_cache()
        # fit() writes best-fit values back into the config it is given, so each
        # pass needs its own copy to read the result out of.
        cfg_run = deepcopy(self.cfg)
        try:
            engine.fit(cfg_run, self.q, I_pert, self.dI)
            x_final, _, _, _ = engine._pack_params(deepcopy(cfg_run))
        except Exception:
            return None
        if len(x_final) != self.n_keys:
            return None
        return [float(v) for v in x_final]


def _mc_cancelled(cancel_cb) -> bool:
    """True when the caller has asked to stop (``cancel_cb`` may be None)."""
    return cancel_cb is not None and bool(cancel_cb())


def _resolve_mc_workers(requested) -> int:
    """Map the ``mc_workers`` convention onto a real process count.

    ``0`` (or None / negative) → auto: ``cpu_count - 2``, leaving a couple of
    cores for the GUI and the OS. ``1`` → serial. ``N > 1`` → N, capped at
    ``cpu_count``.
    """
    n_cpu = os.cpu_count() or 1
    if requested is None or int(requested) <= 0:
        return max(1, n_cpu - 2)
    return max(1, min(int(requested), n_cpu))


def _run_mc_parallel(run, batch, n_workers, record, cancel_cb):
    """Distribute MC passes over a process pool.

    ``record(x_final)`` is called once per completed pass, in completion order.

    Returns the entries of *batch* that were **not** completed — empty when the
    pool handled everything (including a clean cancel). Raises if the pool could
    not be created or the passes could not be submitted, which happens before
    anything is recorded so the caller can safely fall back to serial.
    """
    from concurrent.futures import ProcessPoolExecutor, as_completed
    from concurrent.futures.process import BrokenProcessPool
    from multiprocessing import get_context

    # 'spawn' everywhere: uniform semantics across platforms, and forking a
    # process that already has Qt loaded is not safe. Workers are created during
    # submit(), so the BLAS pin has to cover the submission loop — otherwise
    # every worker re-imports numpy with a full thread pool and oversubscribes.
    with _limit_blas_threads(1):
        executor = ProcessPoolExecutor(
            max_workers=n_workers, mp_context=get_context('spawn'),
        )
        try:
            futures = {executor.submit(run, arr): i for i, arr in enumerate(batch)}
        except Exception:
            executor.shutdown(wait=False)
            raise

    pending = set(range(len(batch)))
    try:
        for fut in as_completed(futures):
            try:
                x_final = fut.result()
            except BrokenProcessPool:
                raise
            except Exception:
                x_final = None          # failed pass — skipped, as in serial
            pending.discard(futures[fut])
            record(x_final)
            if _mc_cancelled(cancel_cb):
                # Passes already running cannot be interrupted; stop handing out
                # new ones and drop whatever is still queued.
                pending.clear()
                break
    except BrokenProcessPool:
        warnings.warn(
            "Monte-Carlo worker pool died; finishing the remaining passes serially.",
            RuntimeWarning,
        )
    finally:
        executor.shutdown(wait=False, cancel_futures=True)

    return [batch[i] for i in sorted(pending)]


from pyirena.core import distributions as D
from pyirena.core.form_factors import bin_widths, build_g_matrix

# ──────────────────────────────────────────────────────────────────────────────
# Serialisation
# ──────────────────────────────────────────────────────────────────────────────

def _jsonify(value):
    """Plain-Python form of a value, for JSON and HDF5.

    Tuples become lists (JSON has no tuple and would not give one back), numpy
    scalars become Python numbers (``json`` refuses them outright), and both
    conversions recurse into the nested dicts a population uses for its
    parameters.
    """
    if isinstance(value, tuple):
        return [_jsonify(v) for v in value]
    if isinstance(value, list):
        return [_jsonify(v) for v in value]
    if isinstance(value, dict):
        return {str(k): _jsonify(v) for k, v in value.items()}
    if isinstance(value, np.generic):
        return value.item()
    return value


def _restore(value, default):
    """Undo :func:`_jsonify` for one field, guided by the field's default.

    The default says what shape the field is meant to have: a tuple default
    means bounds and the two-element list from JSON has to become a tuple
    again, because the fitter unpacks them.  A dict default whose own values
    are tuples (``dist_params_limits`` and friends) gets the same treatment one
    level down.
    """
    if isinstance(default, tuple) and isinstance(value, (list, tuple)):
        return tuple(value)
    if isinstance(default, dict) and isinstance(value, dict):
        out = {}
        for key, item in value.items():
            proto = next(iter(default.values()), None)
            out[key] = (tuple(item)
                        if isinstance(proto, tuple) and isinstance(item, (list, tuple))
                        else item)
        return out
    return value


class _SerialisableDataclass:
    """``to_dict`` / ``from_dict`` driven by the dataclass fields themselves.

    Modeling has six population types, each with fifteen-odd fields; writing
    the two conversions out per type is how the io layer ended up with six
    branches in each direction that had to be kept in step by hand.  Walking
    ``dataclasses.fields()`` instead means a new field is serialised the moment
    it is declared, and an unknown key in an old file is ignored rather than
    raising.
    """

    def to_dict(self) -> dict:
        """This object's fields as a plain, JSON-serialisable dict."""
        return {f.name: _jsonify(getattr(self, f.name)) for f in dataclass_fields(self)}

    @classmethod
    def from_dict(cls, d: dict):
        """Rebuild from :meth:`to_dict`; every missing field takes its default."""
        obj = cls()
        for f in dataclass_fields(cls):
            if d and f.name in d and d[f.name] is not None:
                setattr(obj, f.name, _restore(d[f.name], getattr(obj, f.name)))
        return obj


# ──────────────────────────────────────────────────────────────────────────────
# Dataclasses
# ──────────────────────────────────────────────────────────────────────────────

@dataclass
class SizeDistPopulation(_SerialisableDataclass):
    """All parameters for one size-distribution population.

    Attributes
    ----------
    enabled:       Include this population in the model.
    dist_type:     One of 'gauss', 'lognormal', 'lsw', 'schulz_zimm', 'ardell'.
    dist_params:   Distribution parameter values (keys depend on dist_type).
    dist_params_fit:    {param_name: bool} — True = fit this parameter.
    dist_params_limits: {param_name: (min, max)} — fitting bounds.
    form_factor:   'sphere' or 'spheroid'.
    ff_params:     Form-factor extra parameters (e.g. {'aspect_ratio': 1.0}).
    ff_params_fit: {param_name: bool}.
    ff_params_limits: {param_name: (min, max)}.
    structure_factor: 'none', 'interferences', or 'hard_sphere'.
    sf_params:     Structure-factor parameters.
                   Interferences: {'eta': float, 'pack': float}
                   Hard-Sphere:   {'radius': float, 'volume_fraction': float}
    sf_params_fit, sf_params_limits: analogous to dist_params_*.
    contrast:      (Δρ)² in units of 10²⁰ cm⁻⁴.
    scale:         Volume-fraction proxy: scale ≈ Vf*(1-Vf).
                   The GUI shows Vf = 0.5*(1 - sqrt(1 - 4*scale)).
    use_number_dist: If True the user-supplied distribution is a number
                   distribution; otherwise it is a volume distribution.
    n_bins:        Number of radius bins for the radius grid (default 200).
    """
    pop_type: str = 'size_dist'
    enabled: bool = True
    dist_type: str = 'lognormal'
    dist_params: dict = field(default_factory=lambda: {
        'min_size': 10.0, 'mean_size': 100.0, 'sdeviation': 0.3
    })
    dist_params_fit: dict = field(default_factory=lambda: {
        'min_size': False, 'mean_size': True, 'sdeviation': True
    })
    dist_params_limits: dict = field(default_factory=lambda: {
        'min_size': (0.0, 1e6), 'mean_size': (1.0, 1e6), 'sdeviation': (0.01, 5.0)
    })

    form_factor: str = 'sphere'
    ff_params: dict = field(default_factory=dict)
    ff_params_fit: dict = field(default_factory=dict)
    ff_params_limits: dict = field(default_factory=dict)

    structure_factor: str = 'none'
    sf_params: dict = field(default_factory=lambda: {
        'eta': 100.0, 'pack': 0.1,          # interferences
        'radius': 50.0, 'volume_fraction': 0.1,  # hard sphere
    })
    sf_params_fit: dict = field(default_factory=lambda: {
        'eta': False, 'pack': False,
        'radius': False, 'volume_fraction': False,
    })
    sf_params_limits: dict = field(default_factory=lambda: {
        'eta': (1.0, 1e6), 'pack': (0.0, 16.0),
        'radius': (1.0, 1e6), 'volume_fraction': (0.0, 0.74),
    })

    contrast: float = 1.0
    fit_contrast: bool = False
    contrast_limits: tuple = (0.0, 1e10)

    scale: float = 0.001
    fit_scale: bool = True
    scale_limits: tuple = (1e-8, 1.0)

    use_number_dist: bool = False
    n_bins: int = 200
    label: str = ''


@dataclass
class UnifiedLevelPopulation(_SerialisableDataclass):
    """Single Beaucage Unified Fit level as a Modeling population.

    Computes: I(q) = G·exp(-q²Rg²/3) + B·Q*⁻ᴾ·exp(-q²RgCO²/3)
    where Q* = q / erf(K·q·Rg/√6)³  and  K = 1.0 if P>3 else 1.06.
    Optional Born-Green correlation factor: divide by 1 + PACK·f(q,ETA).
    """
    pop_type: str = 'unified_level'
    enabled: bool = True
    G: float = 1.0
    fit_G: bool = True
    G_limits: tuple = (1e-10, 1e10)
    Rg: float = 10.0
    fit_Rg: bool = True
    Rg_limits: tuple = (0.1, 1e6)
    B: float = 1e-4
    fit_B: bool = True
    B_limits: tuple = (1e-20, 1e10)
    P: float = 4.0
    fit_P: bool = False
    P_limits: tuple = (0.0, 6.0)
    RgCO: float = 0.0
    fit_RgCO: bool = False
    RgCO_limits: tuple = (0.0, 1e6)
    correlations: bool = False
    ETA: float = 10.0
    fit_ETA: bool = False
    ETA_limits: tuple = (0.1, 1e6)
    PACK: float = 0.0
    fit_PACK: bool = False
    PACK_limits: tuple = (0.0, 16.0)
    label: str = ''


@dataclass
class DiffractionPeakPopulation(_SerialisableDataclass):
    """Diffraction peak (Gaussian, Lorentzian, or pseudo-Voigt) as a Modeling population.

    Gaussian:    A · exp(-(q-q₀)² / (2σ²))
    Lorentzian:  A / (1 + ((q-q₀)/σ)²)
    pseudo-Voigt:  A · (η·Lorentzian + (1-η)·Gaussian)
    """
    pop_type: str = 'diffraction_peak'
    enabled: bool = True
    peak_type: str = 'gaussian'      # 'gaussian' | 'lorentzian' | 'voigt'
    position: float = 0.1
    fit_position: bool = True
    position_limits: tuple = (0.001, 10.0)
    amplitude: float = 1.0
    fit_amplitude: bool = True
    amplitude_limits: tuple = (0.0, 1e10)
    width: float = 0.01
    fit_width: bool = True
    width_limits: tuple = (1e-6, 10.0)
    eta_voigt: float = 0.5           # mixing: 0 = pure Gaussian, 1 = pure Lorentzian
    fit_eta_voigt: bool = False
    eta_voigt_limits: tuple = (0.0, 1.0)
    label: str = ''


@dataclass
class GuinierPorodPopulation(_SerialisableDataclass):
    """Guinier-Porod piecewise scattering model (Hammouda 2010, J. Appl. Cryst. 43, 716).

    Always 6 shape parameters: G, Rg1, s1, P, Rg2, s2.
    Default Rg2=1e10 collapses the low-Q regime (single-level behaviour).
    Optional: RgCO cutoff; Born-Green correlations (ETA, PACK).
    """
    pop_type: str = 'guinier_porod'
    enabled: bool = True
    G:   float = 1.0
    fit_G: bool = True
    G_limits: tuple = (1e-10, 1e10)
    Rg1: float = 10.0
    fit_Rg1: bool = True
    Rg1_limits: tuple = (0.1, 1e6)
    s1:  float = 0.0
    fit_s1: bool = False
    s1_limits: tuple = (0.0, 3.0)
    P:   float = 4.0
    fit_P: bool = False
    P_limits: tuple = (0.0, 6.0)
    Rg2: float = 1e10
    fit_Rg2: bool = False
    Rg2_limits: tuple = (0.1, 1e12)
    s2:  float = 0.0
    fit_s2: bool = False
    s2_limits: tuple = (0.0, 3.0)
    RgCO: float = 0.0
    fit_RgCO: bool = False
    RgCO_limits: tuple = (0.0, 1e6)
    correlations: bool = False
    ETA:  float = 10.0
    fit_ETA: bool = False
    ETA_limits: tuple = (0.1, 1e6)
    PACK: float = 0.0
    fit_PACK: bool = False
    PACK_limits: tuple = (0.0, 16.0)
    label: str = ''


@dataclass
class MassFractalPopulation(_SerialisableDataclass):
    """Mass fractal aggregate scattering (Teixeira 1988, J. Appl. Cryst. 21, 781).

    Sphere primary particles (Beta=1) with a fractal structure factor.
    Phi · Contrast · V · [S_fractal(Q) + (1-Eta)²] · |F(qR)|²
    """
    pop_type: str = 'mass_fractal'
    enabled: bool = True
    Phi:      float = 0.001
    fit_Phi: bool = True
    Phi_limits: tuple = (1e-8, 1.0)
    Radius:   float = 50.0
    fit_Radius: bool = True
    Radius_limits: tuple = (0.1, 1e6)
    Beta:     float = 1.0    # aspect ratio (1=sphere, <1=oblate, >1=prolate)
    fit_Beta: bool = False
    Beta_limits: tuple = (0.01, 100.0)
    Dv:       float = 2.5
    fit_Dv: bool = True
    Dv_limits: tuple = (1.0, 3.0)
    Ksi:      float = 500.0
    fit_Ksi: bool = True
    Ksi_limits: tuple = (1.0, 1e7)
    Eta:      float = 0.5
    fit_Eta: bool = False
    Eta_limits: tuple = (0.3, 0.8)
    Contrast: float = 1.0
    fit_Contrast: bool = False
    Contrast_limits: tuple = (0.0, 1e10)
    label: str = ''


@dataclass
class SurfaceFractalPopulation(_SerialisableDataclass):
    """Surface fractal scattering (Teixeira 1988, J. Appl. Cryst. 21, 781).

    π · Contrast · Ksi⁴ · Surface · Γ(5-Ds) · sin((3-Ds)·atan(Q·Ksi))
      / ((1+(Q·Ksi)²)^((5-Ds)/2) · Q·Ksi)
    Optional smooth Porod transition to Q⁻⁴ above Qc.
    """
    pop_type: str = 'surface_fractal'
    enabled: bool = True
    Surface:  float = 1e4
    fit_Surface: bool = True
    Surface_limits: tuple = (1.0, 1e12)
    Ds:       float = 2.5
    fit_Ds: bool = True
    Ds_limits: tuple = (2.0, 3.0)
    Ksi:      float = 500.0
    fit_Ksi: bool = True
    Ksi_limits: tuple = (1.0, 1e7)
    Contrast: float = 1.0
    fit_Contrast: bool = False
    Contrast_limits: tuple = (0.0, 1e10)
    use_porod_transition: bool = False
    Qc:       float = 0.1
    fit_Qc: bool = False
    Qc_limits: tuple = (0.001, 10.0)
    QcWidth:  float = 0.1
    fit_QcWidth: bool = False
    QcWidth_limits: tuple = (0.01, 1.0)
    label: str = ''


@dataclass
class ModelingConfig:
    """Global configuration for a Modeling fit."""
    populations: list = field(default_factory=list)   # list of SizeDistPopulation
    background: float = 0.0
    fit_background: bool = True
    background_limits: tuple = (0.0, 1e10)
    q_min: float = 0.001
    q_max: float = 1.0
    no_limits: bool = False
    n_mc_runs: int = 10
    # 'local'  → scipy least_squares (TRF) / Nelder-Mead in no-limits mode (default).
    # 'global' → differential_evolution to locate the basin, then a TRF local
    #            polish. Requires finite bounds, so it is forced to 'local'
    #            whenever no_limits is set. Best for monodisperse core-shell /
    #            core-shell-shell models whose χ² surface is highly multimodal.
    fit_method: str = 'local'
    # Worker processes for the global (DE) search. 1 = serial (default); N > 1
    # or -1 (all cores) evaluates the DE population in parallel. Only affects
    # the global method; the serial path is unchanged. Falls back to serial
    # automatically if multiprocessing fails on the host.
    de_workers: int = 1
    # Worker processes for the Monte-Carlo uncertainty passes (each pass is an
    # independent refit of noise-perturbed data, so they parallelise cleanly).
    # 0 = auto (cpu_count - 2, the default); 1 = serial; N > 1 = explicit count.
    # The first pass always runs in-process and is timed: when the remaining
    # work is too small to repay pool startup the rest also runs serially.
    # Falls back to serial automatically if multiprocessing fails on the host.
    mc_workers: int = 0
    # Slit smearing: when enabled the total model intensity is slit smeared
    # before comparison with (slit-smeared) data, so fitted parameters are
    # ideal-space.  Structure factors are part of the ideal model and are
    # applied before smearing (they live inside total_intensity).
    use_slit_smearing: bool = False
    slit_length: float = 0.0

    def to_dict(self) -> dict:
        """The whole fit setup as a plain dict — populations included.

        Settings only, in the same spirit as the other cores: the data, the
        fitted curve and the χ² of the last run are results and live in
        :class:`ModelingResult`.
        """
        return {
            'populations':       [pop.to_dict() for pop in self.populations],
            'background':        float(self.background),
            'fit_background':    bool(self.fit_background),
            'background_limits': [float(v) for v in self.background_limits],
            'q_min':             float(self.q_min),
            'q_max':             float(self.q_max),
            'no_limits':         bool(self.no_limits),
            'n_mc_runs':         int(self.n_mc_runs),
            'fit_method':        str(self.fit_method),
            'de_workers':        int(self.de_workers),
            'mc_workers':        int(self.mc_workers),
            'use_slit_smearing': bool(self.use_slit_smearing),
            'slit_length':       float(self.slit_length),
        }

    @classmethod
    def from_dict(cls, d: dict) -> 'ModelingConfig':
        """Rebuild a configuration from :meth:`to_dict`.

        Every field falls back to its default, so a file written before
        ``fit_method`` or the worker counts existed opens with the serial local
        fit it was actually run with.
        """
        d = d or {}
        cfg = cls()
        cfg.populations = [population_from_dict(pd) for pd in (d.get('populations') or [])]
        for name in ('background', 'q_min', 'q_max', 'slit_length'):
            if d.get(name) is not None:
                setattr(cfg, name, float(d[name]))
        for name in ('fit_background', 'no_limits', 'use_slit_smearing'):
            if d.get(name) is not None:
                setattr(cfg, name, bool(d[name]))
        for name in ('n_mc_runs', 'de_workers', 'mc_workers'):
            if d.get(name) is not None:
                setattr(cfg, name, int(d[name]))
        if d.get('fit_method'):
            cfg.fit_method = str(d['fit_method'])
        limits = d.get('background_limits')
        if limits is not None and len(limits) == 2:
            cfg.background_limits = (float(limits[0]), float(limits[1]))
        return cfg


#: Every population type, keyed by the ``pop_type`` written to file.  One
#: mapping instead of the ``if pop_type == …`` ladders that had to be repeated
#: wherever a population was rebuilt.
POPULATION_CLASSES = {
    'size_dist':        SizeDistPopulation,
    'unified_level':    UnifiedLevelPopulation,
    'diffraction_peak': DiffractionPeakPopulation,
    'guinier_porod':    GuinierPorodPopulation,
    'mass_fractal':     MassFractalPopulation,
    'surface_fractal':  SurfaceFractalPopulation,
}


def population_from_dict(d: dict):
    """Rebuild the right population class from a saved dict.

    An unrecognised ``pop_type`` — a file from a newer pyIrena — falls back to
    a size distribution rather than raising, so the rest of the fit still
    opens and the user sees something they can correct.
    """
    d = d or {}
    pop_type = d.get('pop_type', 'size_dist')
    cls = POPULATION_CLASSES.get(pop_type)
    if cls is None:
        log.warning("unknown population type %r; loading it as a size distribution",
                    pop_type)
        cls = SizeDistPopulation
    return cls.from_dict(d)


@dataclass
class ModelingResult:
    """Results from a Modeling fit."""
    config: ModelingConfig
    chi_squared: float
    reduced_chi_squared: float
    dof: int
    timestamp: str

    model_q: np.ndarray              # Q-array used for fitting
    model_I: np.ndarray              # Total model intensity (slit smeared when smearing on)

    # Per-population arrays (only enabled populations; order matches enabled pops)
    pop_indices: list                # list of int — original population indices
    pop_model_I: list                # list of np.ndarray
    radius_grids: list               # list of np.ndarray   [Å]
    volume_dists: list               # list of np.ndarray   [Å⁻¹ normalised]
    number_dists: list               # list of np.ndarray   [Å⁻¹ normalised]
    derived: list                    # list of dict per active population

    # MC uncertainty: {pop_i_paramname: std} and {'background': std}
    params_std: dict = field(default_factory=dict)

    # Ideal (pinhole) total on model_q; only set when slit smearing was used
    # (None otherwise).  Saved alongside model_I so downstream plots/readers can
    # offer a smeared/ideal toggle (parity with Unified/Sizes/Simple).
    model_I_ideal: Optional[np.ndarray] = None

    # Plain-language diagnostics from the fit (parameters pinned at limits,
    # non-convergence, unphysically broad distribution, huge reduced χ²).
    # Transient — shown by the GUI/batch layers, not persisted to HDF5.
    fit_warnings: list = field(default_factory=list)


def result_to_report_dict(result: "ModelingResult",
                          fit_quality: Optional[dict] = None) -> dict:
    """Flatten a :class:`ModelingResult` into the saved-results dict shape.

    ``pyirena.core.reporting`` (and everything else that reads modeling results)
    works on the dict ``load_modeling_results()`` returns, not on the live
    population objects.  This converter lets a panel report the fit it is
    holding without saving to HDF5 and loading it back.

    Populations are flattened with ``dataclasses.asdict``, so a new population
    type or parameter appears in the report automatically — there is no
    per-type list here to fall out of date.  Only enabled populations that took
    part in the fit are included, each carrying its ``derived`` quantities.

    Args:
        result: The fit result to describe.
        fit_quality: Optional robust fit-quality metrics dict.

    Returns:
        A dict with ``chi_squared``, ``reduced_chi_squared``, ``dof``,
        ``background``, ``q_min``, ``q_max``, ``slit_length``, ``timestamp``
        and a ``populations`` list.
    """
    from dataclasses import asdict, is_dataclass

    cfg = result.config
    populations = []
    for k, pop_index in enumerate(result.pop_indices):
        pop = cfg.populations[pop_index]
        flat = asdict(pop) if is_dataclass(pop) else dict(vars(pop))
        flat["population_index"] = pop_index + 1
        flat["enabled"] = True
        flat["derived"] = result.derived[k] if k < len(result.derived) else {}
        populations.append(flat)

    slit_length = float(getattr(cfg, "slit_length", 0.0) or 0.0)
    if not getattr(cfg, "use_slit_smearing", False):
        slit_length = 0.0

    return {
        "chi_squared":         float(result.chi_squared),
        "reduced_chi_squared": float(result.reduced_chi_squared),
        "dof":                 int(result.dof),
        "background":          float(cfg.background),
        "q_min":               float(cfg.q_min),
        "q_max":               float(cfg.q_max),
        "slit_length":         slit_length,
        "timestamp":           result.timestamp,
        "populations":         populations,
        "fit_quality":         fit_quality,
    }


# ──────────────────────────────────────────────────────────────────────────────
# Structure-factor helpers  (stand-alone, no class state needed)
# ──────────────────────────────────────────────────────────────────────────────

def _interferences_sf(q: np.ndarray, eta: float, pack: float) -> np.ndarray:
    """Born-Green interference structure factor: S(q) = 1/(1 + pack*f(q,eta))."""
    return 1.0 / (1.0 + pack * _sphere_amplitude(q, eta))


def _hard_sphere_sf(q: np.ndarray, radius: float,
                    volume_fraction: float) -> np.ndarray:
    """Percus-Yevick hard-sphere structure factor.

    Direct port of ``IR2S_HardSphereStruct`` from IR2_StructureFactors.ipf.

    Parameters
    ----------
    q:               Scattering vector [Å⁻¹]
    radius:          Hard-sphere radius [Å]
    volume_fraction: Volume fraction φ  (0 … ~0.74)
    """
    phi = float(volume_fraction)
    r = float(radius)

    if phi <= 0 or r <= 0:
        return np.ones_like(q)

    DENOM = (1.0 - phi) ** 4
    DNUM = (1.0 + 2.0 * phi) ** 2
    ALPHA = DNUM / DENOM
    BETA = -6.0 * phi * (1.0 + phi / 2.0) ** 2 / DENOM
    GAMMA = 0.5 * phi * DNUM / DENOM

    A = 2.0 * q * r
    # Protect against A=0
    A_safe = np.maximum(A, 1e-10)
    ASQ = A_safe ** 2
    ATH = ASQ * A_safe
    AFOR = ATH * A_safe

    RCA = np.cos(A_safe)
    RSA = np.sin(A_safe)

    CALP = ALPHA * (RSA / ASQ - RCA / A_safe)
    CBETA = BETA * (
        2.0 * RSA / ASQ - (ASQ - 2.0) * RCA / ATH - 2.0 / ATH
    )
    CGAM = GAMMA * (
        -RCA / A_safe
        + (4.0 / A_safe) * (
            (3.0 * ASQ - 6.0) * RCA / AFOR
            + (ASQ - 6.0) * RSA / ATH
            + 6.0 / AFOR
        )
    )

    PREFAC = -24.0 * phi / A_safe
    C = PREFAC * (CALP + CBETA + CGAM)
    S = 1.0 / (1.0 - C)

    # At q=0 (A→0) S → 1 numerically; set explicitly to be safe
    if np.ndim(q) == 0:
        return float(S)
    return np.where(A < 1e-10, 1.0, S)


def _unified_level_intensity(
    q: np.ndarray, pop: UnifiedLevelPopulation,
) -> np.ndarray:
    """Single Beaucage Unified level intensity.  Port of unified.py calculate_level_intensity."""
    K = 1.0 if pop.P > 3.0 else 1.06
    Rg = max(pop.Rg, 1e-10)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        erf_term = np.maximum(_scipy_erf(K * q * Rg / np.sqrt(6.0)), 1e-10)
        Q_star   = q / (erf_term ** 3)
        guinier  = pop.G * np.exp(-q ** 2 * Rg ** 2 / 3.0)
        cutoff   = np.exp(-pop.RgCO ** 2 * q ** 2 / 3.0) if pop.RgCO > 0 else 1.0
        powerlaw = pop.B / np.maximum(Q_star, 1e-100) ** pop.P * cutoff
    I = guinier + powerlaw
    if pop.correlations and pop.PACK > 0:
        I /= (1.0 + pop.PACK * _sphere_amplitude(q, pop.ETA))
    return I


def _diffraction_peak_intensity(
    q: np.ndarray, pop: DiffractionPeakPopulation,
) -> np.ndarray:
    """Gaussian, Lorentzian, or pseudo-Voigt diffraction peak."""
    dq = q - pop.position
    w  = max(pop.width, 1e-10)
    g  = np.exp(-dq ** 2 / (2.0 * w ** 2))
    l  = 1.0 / (1.0 + (dq / w) ** 2)
    if pop.peak_type == 'gaussian':
        return pop.amplitude * g
    if pop.peak_type == 'lorentzian':
        return pop.amplitude * l
    # pseudo-Voigt
    eta = float(np.clip(pop.eta_voigt, 0.0, 1.0))
    return pop.amplitude * (eta * l + (1.0 - eta) * g)


def _guinier_porod_intensity(
    q: np.ndarray, pop: 'GuinierPorodPopulation',
) -> np.ndarray:
    """Guinier-Porod piecewise model (Hammouda 2010).

    6 shape params: G, Rg1, s1, P, Rg2, s2.
    Rg2=1e10 default collapses the Q2 crossover (single-level behaviour).
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        Rg1 = max(pop.Rg1, 1e-10)
        # Guard against s1 == P (D→0/0) or P <= s1
        P_eff = pop.P if pop.P > pop.s1 + 1e-6 else pop.s1 + 1e-6

        Q1 = np.sqrt((P_eff - pop.s1) * (3.0 - pop.s1) / 2.0) / Rg1
        Q1 = max(Q1, 1e-30)
        D  = (pop.G * np.exp(-(P_eff - pop.s1) / 2.0)
              * ((3.0 - pop.s1) * (P_eff - pop.s1) / 2.0) ** ((P_eff - pop.s1) / 2.0)
              / Rg1 ** (P_eff - pop.s1))

        # Build middle + high-Q regime
        q_s1 = np.where(q < 1e-20, 1.0, q ** pop.s1)
        guinier1  = pop.G / np.maximum(q_s1, 1e-100) * np.exp(-q ** 2 * Rg1 ** 2 / (3.0 - pop.s1))
        powerlaw  = D / np.maximum(q, 1e-100) ** P_eff
        I = np.where(q < Q1, guinier1, powerlaw)

        # Second Guinier-Porod regime (active only when Q2 > 0 and Q2 < Q1)
        Rg2 = max(pop.Rg2, 1e-10)
        denom = 2.0 / (3.0 - pop.s2) * Rg2 ** 2 - 2.0 / (3.0 - pop.s1) * Rg1 ** 2
        if abs(denom) > 1e-30:
            Q2sq = (1.0 - pop.s2) / denom
            if Q2sq > 0:
                Q2 = np.sqrt(Q2sq)
                if 0 < Q2 < Q1:
                    G2 = (pop.G
                          * np.exp(-Q2 ** 2 * (Rg1 ** 2 / (3.0 - pop.s1)
                                               - Rg2 ** 2 / (3.0 - pop.s2)))
                          * Q2 ** (pop.s2 - pop.s1))
                    q_s2 = np.where(q < 1e-20, 1.0, q ** pop.s2)
                    guinier2 = G2 / np.maximum(q_s2, 1e-100) * np.exp(
                        -q ** 2 * Rg2 ** 2 / (3.0 - pop.s2))
                    I = np.where(q < Q2, guinier2, I)

    if pop.RgCO > 0:
        I = I * np.exp(-pop.RgCO ** 2 * q ** 2 / 3.0)
    if pop.correlations and pop.PACK > 0:
        I = I / (1.0 + pop.PACK * _sphere_amplitude(q, pop.ETA))
    return I


def _chi_s_mf(Beta: float) -> float:
    """CHiS factor for a spheroid (IR1V_CaculateChiS from IR1_Fractals.ipf)."""
    if Beta < 1.0 - 1e-6:
        return (1.0 / (2.0 * Beta)) * (
            1.0 + (Beta ** 2 / np.sqrt(1.0 - Beta ** 2))
            * np.log((1.0 + np.sqrt(1.0 - Beta ** 2)) / Beta)
        )
    if Beta > 1.0 + 1e-6:
        return (1.0 / (2.0 * Beta)) * (
            1.0 + (Beta ** 2 / np.sqrt(Beta ** 2 - 1.0))
            * np.arcsin(np.sqrt(Beta ** 2 - 1.0) / Beta)
        )
    return 1.0   # sphere


def _spheroid_ff2_avg(
    q: np.ndarray, R: float, Beta: float, n_pts: int = 50,
) -> np.ndarray:
    """Orientation-averaged spheroid |F(q)|² via Gauss-Legendre quadrature.

    Ports IR1V_CalculateFSquared from IR1_Fractals.ipf.
    Integration variable x = cos(theta), 0..1.
    P(q) = 9 * [∫₀¹ h(qR√(1+(β²-1)x²)) dx]²  where h(z) = (sin z - z cos z)/z³
    """
    from scipy.special import roots_legendre
    xi, wi = roots_legendre(n_pts)
    x = (xi + 1.0) / 2.0          # map [-1,1] → [0,1]
    w = wi / 2.0                   # Jacobian

    # effective radius factor u = sqrt(1 + (Beta²-1)·x²), shape (n_pts,)
    u = np.sqrt(np.maximum(1.0 + (Beta ** 2 - 1.0) * x ** 2, 0.0))

    # qRu: shape (n_q, n_pts)
    qRu = np.outer(q, R * u)

    # h(z) = (sin z - z cos z)/z³,  h(0) = 1/3
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        h = np.where(
            qRu < 1e-6,
            1.0 / 3.0,
            (np.sin(qRu) - qRu * np.cos(qRu)) / np.maximum(qRu ** 3, 1e-200),
        )

    integral = h @ w          # (n_q, n_pts) @ (n_pts,) → (n_q,)
    return 9.0 * integral ** 2


def _mass_fractal_intensity(
    q: np.ndarray, pop: 'MassFractalPopulation',
) -> np.ndarray:
    """Mass fractal aggregate scattering (Teixeira 1988).

    Primary particles are spheroids with aspect ratio Beta (1 = sphere).
    For Beta ≠ 1 uses orientation-averaged spheroid form factor via GL quadrature.
    """
    R    = max(pop.Radius, 1e-10)
    Dv   = float(np.clip(pop.Dv, 1.001, 2.999))
    Ksi  = max(pop.Ksi, 1e-10)
    Beta = max(getattr(pop, 'Beta', 1.0), 1e-3)

    # Spheroid volume (4/3)·π·R²·(Beta·R) = (4/3)·π·R³·Beta
    V = (4.0 / 3.0) * np.pi * R ** 3 * Beta

    # CHiS and RC generalise the sphere case (Beta=1 → ChiS=1, RC=2R)
    ChiS = _chi_s_mf(Beta)
    RC   = R * np.sqrt(2.0) / ChiS * np.sqrt(1.0 + ((2.0 + Beta ** 2) / 3.0) * ChiS ** 2)
    RC   = max(RC, 1e-10)

    Bracket = pop.Eta * RC ** 3 / (Beta * R ** 3) * (Ksi / RC) ** Dv

    qKsi = q * Ksi
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sf_num   = np.sin((Dv - 1.0) * np.arctan(qKsi))
        sf_denom = (Dv - 1.0) * qKsi * (1.0 + qKsi ** 2) ** ((Dv - 1.0) / 2.0)
        frac_sf  = np.where(qKsi < 1e-6, 1.0, sf_num / np.maximum(sf_denom, 1e-100))

    # Form factor: sphere path avoids integration overhead for Beta ≈ 1
    if abs(Beta - 1.0) < 0.01:
        qR  = q * R
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ff2 = np.where(
                qR < 1e-6,
                1.0,
                (3.0 * (np.sin(qR) - qR * np.cos(qR)) / np.maximum(qR ** 3, 1e-100)) ** 2,
            )
    else:
        ff2 = _spheroid_ff2_avg(q, R, Beta)

    return pop.Phi * pop.Contrast * 1e-4 * V * (Bracket * frac_sf + (1.0 - pop.Eta) ** 2) * ff2


def _surface_fractal_intensity(
    q: np.ndarray, pop: 'SurfaceFractalPopulation',
) -> np.ndarray:
    """Surface fractal scattering (Teixeira 1988).

    Uses gammaln for numerical stability near Ds→3.
    Optional smooth Porod transition above Qc (erf blending).
    """
    from scipy.special import erf as scipy_erf
    from scipy.special import gammaln

    Ksi = max(pop.Ksi, 1e-10)
    Ds  = float(np.clip(pop.Ds, 2.001, 2.999))

    qKsi = q * Ksi
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sf_num   = np.sin((3.0 - Ds) * np.arctan(qKsi))
        sf_denom = (1.0 + qKsi ** 2) ** ((5.0 - Ds) / 2.0) * qKsi
        gamma_val = np.exp(gammaln(5.0 - Ds))
        I = (np.pi * pop.Contrast * 1e20 * Ksi ** 4 * 1e-32
             * pop.Surface * gamma_val
             * np.where(qKsi < 1e-6, 0.0,
                        sf_num / np.maximum(sf_denom, 1e-100)))

    if pop.use_porod_transition and pop.Qc > 0:
        sigma = pop.Qc * pop.QcWidth / 2.3548
        if sigma > 1e-10:
            step1 = 0.5 * (1.0 + scipy_erf((pop.Qc - q) / (np.sqrt(2.0) * sigma)))
            # Continuity: A = I(Qc) * Qc^4
            i_qc  = np.interp(pop.Qc, q, I)
            A_por = i_qc * pop.Qc ** 4
            I = I * step1 + A_por * np.maximum(q, 1e-30) ** -4 * (1.0 - step1)

    return I


def apply_structure_factor(
    I_raw: np.ndarray,
    q: np.ndarray,
    pop: SizeDistPopulation,
) -> np.ndarray:
    """Multiply I_raw by S(Q) for the population's structure factor."""
    sf = pop.structure_factor.lower()
    if sf == 'none':
        return I_raw
    if sf == 'interferences':
        eta = float(pop.sf_params.get('eta', 100.0))
        pack = float(pop.sf_params.get('pack', 0.1))
        return I_raw * _interferences_sf(q, eta, pack)
    if sf == 'hard_sphere':
        rad = float(pop.sf_params.get('radius', 50.0))
        vf = float(pop.sf_params.get('volume_fraction', 0.1))
        return I_raw * _hard_sphere_sf(q, rad, vf)
    raise ValueError(f"Unknown structure factor: {pop.structure_factor!r}")


# ──────────────────────────────────────────────────────────────────────────────
# Fit diagnostics — plain-language warnings returned on ModelingResult
# ──────────────────────────────────────────────────────────────────────────────

def _format_param_key(key: tuple) -> str:
    """Human-readable name for a ``_pack_params`` key tuple."""
    if tuple(key) == ('background',):
        return 'background'
    return f"population {key[1] + 1} {key[-1]}"


def _pinned_parameter_warnings(
    x: np.ndarray, lo: np.ndarray, hi: np.ndarray, keys: list,
) -> list:
    """Warn for every fitted parameter that ended pinned at a fit limit.

    A parameter riding its limit means the optimizer wanted to go further —
    the reported value is dictated by the limit, not by the data. On badly
    misspecified models (see testData/bad1.h5, where ``scale`` pinned at its
    upper limit while the model sat ~600× below the low-Q data) this is the
    clearest available signal of what is actually wrong.

    The closeness test is relative to the limit's own magnitude (not the
    lo–hi span): many limits span decades (e.g. 0 … 1e10), where a
    span-relative test would flag every small value as "at the lower limit".
    """
    out = []
    for xv, lv, hv, key in zip(np.atleast_1d(x), lo, hi, keys):
        name = _format_param_key(key)
        if xv <= lv + 1e-3 * max(abs(lv), 1e-12):
            out.append(
                f"{name} = {xv:.4g} is pinned at its lower fit limit ({lv:g})")
        elif xv >= hv - 1e-3 * max(abs(hv), 1e-12):
            out.append(
                f"{name} = {xv:.4g} is pinned at its upper fit limit ({hv:g})")
    return out


def _model_sanity_warnings(cfg: 'ModelingConfig', rchi2: float) -> list:
    """Post-fit sanity checks on the fitted configuration.

    Flags a lognormal ``sdeviation`` above 3 (the radius grid then spans from
    Å to mm — the population is almost certainly being abused to mimic a
    power law) and a grossly large reduced χ² (the model cannot reach the
    data anywhere; more restarts will not fix that).
    """
    out = []
    for i, pop in enumerate(cfg.populations):
        if not getattr(pop, 'enabled', False):
            continue
        if getattr(pop, 'pop_type', 'size_dist') != 'size_dist':
            continue
        if getattr(pop, 'dist_type', '') != 'lognormal':
            continue
        s = float(pop.dist_params.get('sdeviation', 0.0))
        if s > 3.0:
            out.append(
                f"population {i + 1} lognormal sdeviation = {s:.3g} — the "
                f"distribution spans an unphysically wide size range; the "
                f"population is likely compensating for scattering the model "
                f"cannot otherwise describe (try adding a power-law/unified "
                f"level or another population)"
            )
    if rchi2 > 1e3:
        out.append(
            f"reduced χ² = {rchi2:.3g} — the model is far from the data; "
            f"check whether the populations can reach the measured intensity "
            f"over the whole fitted Q range"
        )
    return out


# ──────────────────────────────────────────────────────────────────────────────
# Derived-quantity helpers
# ──────────────────────────────────────────────────────────────────────────────

def _volume_sphere(r: np.ndarray) -> np.ndarray:
    return (4.0 / 3.0) * np.pi * r ** 3


def _volume_cylinder(r_grid: np.ndarray, pop) -> np.ndarray:
    """Cylinder volume per radius bin, based on the FF parameterisation.

    cylinder_ar:     V = π r² · 2L = 2π·AR·r³    (L = AR·r varies with r)
    cylinder_length: V = π r² · H                 (H = length is constant)
    """
    ff = getattr(pop, 'form_factor', 'sphere')
    if ff == 'cylinder_ar':
        AR = float(pop.ff_params.get('aspect_ratio', 1.0))
        return 2.0 * np.pi * AR * r_grid ** 3
    if ff == 'cylinder_length':
        length = float(pop.ff_params.get('length', 100.0))
        return np.pi * r_grid ** 2 * length
    return _volume_sphere(r_grid)


def _volume_coreshell(r_grid: np.ndarray, pop) -> np.ndarray:
    """Total particle volume per bin for core-shell FF modes.

    Volume is always computed from the *outer* (total) radius R_total so
    that ``pop.scale`` correctly encodes the total particle volume fraction.

    by_core / by_spheroid_core:  r_grid = R_core;  R_total = R_core + t_shell
    by_shell:                    r_grid = t_shell;  R_total = r_core_fixed + t_shell
    by_total / by_spheroid_total: r_grid = R_total
    """
    ff = getattr(pop, 'form_factor', 'sphere')
    if ff in ('cs_sphere_by_core', 'cs_spheroid_by_core'):
        t = float(pop.ff_params.get('t_shell', 0.0))
        r_total = r_grid + t
    elif ff == 'cs_sphere_by_shell':
        r_cf = float(pop.ff_params.get('r_core_fixed', 0.0))
        r_total = r_cf + r_grid          # r_grid = t_shell here
    else:  # by_total variants
        r_total = r_grid
    return _volume_sphere(r_total)


def _uf_invariant(pop: 'UnifiedLevelPopulation') -> float:
    """Compute the Porod invariant for a Unified Fit Level.

    Q_inv = ∫₀^∞ I(Q)·Q² dQ

    Numerically integrates from 0 to maxQ = 20π/Rg, then adds the Porod-tail
    analytical correction  -B·maxQ^(3-P)/(3-P)  when P > 3 and no RgCO cutoff.

    Returns NaN when P ≤ 3 (integral diverges) or the result is negative.
    Units: cm⁻⁴  (converted from cm⁻¹·Å⁻³ via factor 10²⁴).
    """
    import warnings as _w
    Rg = max(pop.Rg, 1e-10)
    P  = pop.P
    B  = pop.B
    maxQ = 20.0 * np.pi / Rg
    n = 2000
    q = np.linspace(0.0, maxQ, n)
    q[0] = q[1]  # avoid Q=0
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        I = _unified_level_intensity(q, pop)
    I[0] = I[1]
    invariant = float(np.trapezoid(I * q ** 2, q))
    # Porod tail correction  ∫_maxQ^∞ B·Q^(2-P) dQ = -B·maxQ^(3-P)/(3-P)  for P>3
    if pop.RgCO < 0.1 and P > 3.0:
        invariant += -B * maxQ ** (3.0 - P) / (3.0 - P)
    if invariant < 0.0:
        return float('nan')
    return invariant * 1e24   # cm⁻¹·Å⁻³ → cm⁻⁴


def compute_derived(radius_grid, vol_dist, num_dist, pop) -> dict:
    """Compute convenient derived quantities from a fitted population.

    For size-distribution populations returns: vol_mean_r, num_mean_r,
    volume_fraction, Rg, specific_surface.
    For Unified-level populations returns: Rg, G, B, P, invariant.
    For diffraction-peak populations returns: position, amplitude, width.
    """
    pop_type = getattr(pop, 'pop_type', 'size_dist')
    if pop_type == 'unified_level':
        inv = _uf_invariant(pop)
        return {'Rg': pop.Rg, 'G': pop.G, 'B': pop.B, 'P': pop.P,
                'invariant': inv if not np.isnan(inv) else None}
    if pop_type == 'diffraction_peak':
        return {
            'position': pop.position,
            'amplitude': pop.amplitude,
            'width': pop.width,
        }
    if pop_type == 'guinier_porod':
        return {'G': pop.G, 'Rg1': pop.Rg1, 's1': pop.s1, 'P': pop.P,
                'Rg2': pop.Rg2, 's2': pop.s2}
    if pop_type == 'mass_fractal':
        return {'Phi': pop.Phi, 'Radius': pop.Radius,
                'Beta': getattr(pop, 'Beta', 1.0),
                'Dv': pop.Dv, 'Ksi': pop.Ksi, 'Eta': pop.Eta}
    if pop_type == 'surface_fractal':
        return {'Surface': pop.Surface, 'Ds': pop.Ds, 'Ksi': pop.Ksi}
    # size_dist
    dr = bin_widths(radius_grid)
    vol_tot = np.sum(vol_dist * dr)
    num_tot = np.sum(num_dist * dr)

    vol_mean_r = (np.sum(radius_grid * vol_dist * dr) / vol_tot
                  if vol_tot > 0 else 0.0)
    num_mean_r = (np.sum(radius_grid * num_dist * dr) / num_tot
                  if num_tot > 0 else 0.0)

    # Volume fraction from scale ≈ Vf*(1-Vf)  → quadratic
    scale = pop.scale
    if 4.0 * scale > 1.0:
        vf = 0.5
    else:
        vf = 0.5 * (1.0 - np.sqrt(max(1.0 - 4.0 * scale, 0.0)))

    # Rg² (volume-weighted) — geometry depends on form factor
    ff = getattr(pop, 'form_factor', 'sphere')
    is_cs = ff.startswith('cs_')

    if ff == 'cylinder_ar':
        AR = float(pop.ff_params.get('aspect_ratio', 1.0))
        # Rg²(r) = r²/2 + L²/3  where L = AR·r
        rg2_integrand = radius_grid ** 2 * (0.5 + AR ** 2 / 3.0)
    elif ff == 'cylinder_length':
        length = float(pop.ff_params.get('length', 100.0))
        # Rg²(r) = r²/2 + (length/2)²/3 = r²/2 + length²/12
        rg2_integrand = radius_grid ** 2 / 2.0 + (length / 2.0) ** 2 / 3.0
    elif is_cs:
        # Core-shell: SLD-weighted Rg per bin, then volume-averaged.
        # Rg²(R_c, R_t) = [Δρ_c·V_c·3R_c²/5 + Δρ_s·(V_t·3R_t²/5 − V_c·3R_c²/5)]
        #                 / [Δρ_c·V_c + Δρ_s·V_t]
        # Fall back to geometric Rg² = 3/5·R_t² when denominator ≤ 0.
        sld_c = float(pop.ff_params.get('sld_core',    10.0))
        sld_s = float(pop.ff_params.get('sld_shell',    1.0))
        sld_0 = float(pop.ff_params.get('sld_solvent',  9.46))
        d_rho_c = sld_c - sld_s   # 10⁻⁶ Å⁻²
        d_rho_s = sld_s - sld_0
        # Compute R_core and R_total per bin
        t = float(pop.ff_params.get('t_shell', 0.0))
        r_cf = float(pop.ff_params.get('r_core_fixed', 0.0))
        if ff in ('cs_sphere_by_core', 'cs_spheroid_by_core'):
            r_c_arr = radius_grid
            r_t_arr = radius_grid + t
        elif ff == 'cs_sphere_by_shell':
            r_c_arr = np.full_like(radius_grid, r_cf)
            r_t_arr = r_cf + radius_grid
        else:  # by_total
            r_t_arr = radius_grid
            r_c_arr = np.maximum(radius_grid - t, 0.0)
        V_c = _volume_sphere(r_c_arr)
        V_t = _volume_sphere(r_t_arr)
        num_cs = d_rho_c * V_c * (3.0 / 5.0) * r_c_arr ** 2 + \
                 d_rho_s * (V_t * (3.0 / 5.0) * r_t_arr ** 2 -
                            V_c * (3.0 / 5.0) * r_c_arr ** 2)
        den_cs = d_rho_c * V_c + d_rho_s * V_t
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            rg2_integrand = np.where(den_cs > 0, num_cs / den_cs,
                                     (3.0 / 5.0) * r_t_arr ** 2)
    else:
        # Sphere: Rg²(r) = 3/5 · r²
        rg2_integrand = (3.0 / 5.0) * radius_grid ** 2

    if vol_tot > 0:
        Rg2 = np.sum(rg2_integrand * vol_dist * dr) / vol_tot
        Rg = np.sqrt(max(Rg2, 0.0))
    else:
        Rg = 0.0

    # Specific surface Sv = S/V (volume-averaged)  [Å⁻¹]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        safe_r = np.maximum(radius_grid, 1e-30)
        if ff == 'cylinder_ar':
            AR = float(pop.ff_params.get('aspect_ratio', 1.0))
            # Sv(r) = 2/L + 2/r = 2/(AR·r) + 2/r
            sv_integrand = 2.0 / (max(AR, 1e-30) * safe_r) + 2.0 / safe_r
        elif ff == 'cylinder_length':
            length = float(pop.ff_params.get('length', 100.0))
            # Sv(r) = 4/H + 2/r  where H = total height = length
            sv_integrand = 4.0 / max(length, 1e-30) + 2.0 / safe_r
        elif is_cs:
            # Sv based on outer surface only: Sv = 3/R_total (Porod law outer interface)
            safe_rt = np.maximum(r_t_arr, 1e-30)
            sv_integrand = 3.0 / safe_rt
        else:
            # Sphere: Sv(r) = 3/r
            sv_integrand = 3.0 / safe_r

        if vol_tot > 0:
            specific_surface = np.sum(sv_integrand * vol_dist * dr) / vol_tot
        else:
            specific_surface = 0.0

    # For core-shell: also report volume-weighted mean outer radius
    result = {
        'vol_mean_r':       vol_mean_r,
        'num_mean_r':       num_mean_r,
        'volume_fraction':  vf,
        'Rg':               Rg,
        'specific_surface': specific_surface,
    }
    if is_cs and vol_tot > 0:
        result['r_total_mean'] = float(np.sum(r_t_arr * vol_dist * dr) / vol_tot)
    return result


# ──────────────────────────────────────────────────────────────────────────────
# Main engine
# ──────────────────────────────────────────────────────────────────────────────

class ModelingEngine:
    """Compute and fit multi-population size-distribution models for SAS data."""

    def __init__(self) -> None:
        # "Graph Model" cache — keyed on full param hash, persists across calls
        self._g_cache: dict[str, np.ndarray] = {}
        # Fit-time smart cache — keyed on (r_grid, form_factor, ff_params, contrast)
        # only; reused when scale/background/SF params change but G inputs don't.
        # Cleared at the start of every fit() call.
        self._fit_g_cache: dict[int, np.ndarray] = {}
        # pop_index -> (key string, the q array G was built on)
        self._fit_g_keys:  dict[int, tuple] = {}
        # Cached slit smearer, rebuilt when (q, slit_length) change.
        self._smearer = None
        # Most recent radius grid and its bin widths, keyed on the distribution
        # inputs that define them. Building the grid inverts the CDF numerically
        # and is expensive; an optimiser step that leaves dist_params alone (most
        # of a finite-difference Jacobian) reuses it.
        self._rg_cache: dict = {}
        # Trial evaluations rejected during the current fit because a population
        # could not be computed, and the message from the most recent one.
        self._eval_failures: int = 0
        self._eval_failure_last: str = ''
        # Fit-time per-population intensity memo, keyed on every parameter the
        # population's I(Q) depends on. Lets a step that touched only the
        # background or another population skip this one outright.
        # Cleared at the start of every fit() call.
        self._pop_memo: dict = {}

    # ── Slit smearing ─────────────────────────────────────────────────────────

    def _get_smearer(self, q: np.ndarray, config: "ModelingConfig"):
        """Cached SlitSmearer for (q, slit_length); None when smearing is off."""
        if not getattr(config, 'use_slit_smearing', False) or config.slit_length <= 0:
            return None
        q = np.asarray(q, dtype=float)
        sm = self._smearer
        if (sm is not None and sm.slit_length == config.slit_length
                and sm.q.shape == q.shape and np.array_equal(sm.q, q)):
            return sm
        from pyirena.core.smearing import SlitSmearer
        self._smearer = SlitSmearer(q, config.slit_length)
        return self._smearer

    def total_intensity_maybe_smeared(
        self, config: "ModelingConfig", q: np.ndarray, use_cache: bool = True,
        strict: bool = False,
    ):
        """total_intensity, slit-smeared when config requests it.

        Evaluates all populations + background on the extended grid and applies
        the smearing operator to the total and to each per-population curve, so
        both the fit and the plotted overlays are consistent with the data.

        *strict* is forwarded to :meth:`total_intensity`; see there.
        """
        smearer = self._get_smearer(q, config)
        if smearer is None:
            return self.total_intensity(config, q, use_cache=use_cache, strict=strict)
        I_ext, pop_idx, pop_I_ext, pop_dist = self.total_intensity(
            config, smearer.q_ext, use_cache=False, strict=strict)
        I_total = smearer.smear(I_ext)
        pop_I = [smearer.smear(Ip) for Ip in pop_I_ext]
        return I_total, pop_idx, pop_I, pop_dist

    def total_intensity_ideal(self, config: "ModelingConfig", q: np.ndarray):
        """Ideal (pinhole) total on the data grid, or None when smearing is off.

        Cheap companion to :meth:`total_intensity_maybe_smeared` for saving the
        ``model_I_ideal`` curve: a single unsmeared evaluation on *q*.

        Uses the q-keyed full cache (``use_cache=True``), NOT the fit-time smart
        cache: that smart cache keys G on (radius grid, population) but not on
        *q*, so mixing it with the smeared path's ``q_ext`` evaluation in the
        same fit would reuse a wrong-length G.
        """
        if self._get_smearer(q, config) is None:
            return None
        I_ideal, _, _, _ = self.total_intensity(config, q, use_cache=True)
        return I_ideal

    # ── Radius grid ──────────────────────────────────────────────────────────

    def build_radius_grid(self, pop: SizeDistPopulation) -> np.ndarray:
        """Return the non-linear radius grid for *pop*."""
        return self._radius_grid_and_widths(pop)[0]

    def _radius_grid_and_widths(self, pop: SizeDistPopulation):
        """(radius grid, bin widths) for *pop*, memoised on their inputs.

        The grid depends only on the distribution type, its parameters and the
        bin count, so this holds a single entry and returns it whenever those
        are unchanged. Both arrays are treated as read-only by every caller.
        """
        key = (pop.dist_type, tuple(sorted(pop.dist_params.items())), pop.n_bins)
        cache = self._rg_cache
        if cache.get('key') == key:
            return cache['rg'], cache['dr']
        rg = D.generate_radius_grid(
            pop.dist_type, pop.dist_params,
            n_bins=pop.n_bins, tail_precision=0.01,
        )
        dr = bin_widths(rg)
        cache['key'], cache['rg'], cache['dr'] = key, rg, dr
        return rg, dr

    @staticmethod
    def _pop_memo_key(pop: SizeDistPopulation) -> tuple:
        """Every parameter a size-distribution population's I(Q) depends on.

        Used by the fit-time memo in :meth:`total_intensity`. Anything that can
        change the returned curve must appear here.
        """
        return (
            pop.dist_type, tuple(sorted(pop.dist_params.items())), pop.n_bins,
            pop.form_factor, tuple(sorted(pop.ff_params.items())),
            pop.contrast, pop.scale, pop.use_number_dist,
            getattr(pop, 'structure_factor', 'none'),
            tuple(sorted(getattr(pop, 'sf_params', {}).items())),
        )

    # ── Distribution ─────────────────────────────────────────────────────────

    def build_distributions(
        self, pop: SizeDistPopulation, radius_grid: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        """Return (vol_dist, num_dist) arrays normalised so that
        ∫ vol_dist·dr = pop.scale.

        If pop.use_number_dist is True, the analytical PDF is treated as a
        number distribution and converted to volume by multiplying by V(r).
        """
        dr = bin_widths(radius_grid)
        raw_pdf = D.pdf(pop.dist_type, radius_grid, pop.dist_params)

        # Guard against all-zero / NaN
        if not np.any(np.isfinite(raw_pdf)) or np.sum(raw_pdf * dr) <= 0:
            raw_pdf = np.ones_like(radius_grid)

        ff = getattr(pop, 'form_factor', 'sphere')
        if ff.startswith('cylinder'):
            V_r = _volume_cylinder(radius_grid, pop)
        elif ff.startswith('cs_'):
            V_r = _volume_coreshell(radius_grid, pop)
        else:
            V_r = _volume_sphere(radius_grid)      # particle volume per bin  [Å³]

        if pop.use_number_dist:
            num_raw = raw_pdf
            vol_raw = num_raw * V_r
        else:
            vol_raw = raw_pdf
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                num_raw = np.where(V_r > 0, vol_raw / V_r, 0.0)

        # Normalise so ∫ vol_dist dr = pop.scale
        vol_integral = np.sum(vol_raw * dr)
        if vol_integral > 0:
            norm = pop.scale / vol_integral
        else:
            norm = 1.0

        vol_dist = vol_raw * norm
        num_dist = num_raw * norm
        return vol_dist, num_dist

    # ── G-matrix (with cache for "Graph Model") ──────────────────────────────

    @staticmethod
    def _cache_key(pop: SizeDistPopulation, q: np.ndarray) -> str:
        data = {
            'q': q.tobytes().hex()[:32],        # first 32 hex chars as proxy
            'dist_type': pop.dist_type,
            'dist_params': sorted(pop.dist_params.items()),
            'form_factor': pop.form_factor,
            'ff_params': sorted(pop.ff_params.items()),
            'contrast': pop.contrast,
            'n_bins': pop.n_bins,
        }
        return hashlib.md5(json.dumps(data, sort_keys=True).encode()).hexdigest()

    def build_g_matrix(
        self, pop: SizeDistPopulation, q: np.ndarray,
        radius_grid: np.ndarray, use_cache: bool = True,
    ) -> np.ndarray:
        """Build (or retrieve cached) G matrix for *pop*.

        G[i,j] = FF(q_i, r_j) * contrast * 1e-4   [cm⁻¹ per unit vol-frac]
        """
        if use_cache:
            key = self._cache_key(pop, q)
            if key in self._g_cache:
                return self._g_cache[key]

        G = build_g_matrix(
            q=q,
            r_grid=radius_grid,
            shape=pop.form_factor,
            contrast=pop.contrast,
            **pop.ff_params,
        )

        if use_cache:
            self._g_cache[key] = G
        return G

    def clear_cache(self) -> None:
        """Discard every cached intermediate — G matrices, radius grid, memo."""
        self._g_cache.clear()
        self._rg_cache.clear()
        self._pop_memo.clear()

    @staticmethod
    def _fit_g_cache_key(r_grid: np.ndarray, pop: SizeDistPopulation) -> str:
        """Cache key for the fit-time smart G cache.

        Keyed only on what G actually depends on — the radius grid (as a hash
        of its bytes), form factor name, ff_params, and contrast.  Deliberately
        excludes scale, background, and sf_params so that optimizer steps that
        only touch those parameters can reuse the cached G matrix.

        G also depends on *q*, but hashing it on every evaluation would cost
        more than it saves; the caller pairs this key with the q array itself
        and compares that by identity.  See :meth:`_pop_intensity_and_grid`.
        """
        data = {
            'r': hashlib.md5(r_grid.tobytes()).hexdigest(),
            'ff': pop.form_factor,
            'ff_p': sorted(pop.ff_params.items()),
            'c': pop.contrast,
        }
        return hashlib.md5(json.dumps(data, sort_keys=True).encode()).hexdigest()

    # ── Population intensity ─────────────────────────────────────────────────

    def calculate_pop_intensity(
        self, pop: SizeDistPopulation, q: np.ndarray,
        use_cache: bool = True,
        pop_index: Optional[int] = None,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Calculate I(Q) for a single enabled population.

        Returns (I_pop, vol_dist, num_dist).

        Args:
            use_cache:  True  → "Graph Model" path: full hash cache persists
                                 across calls (params are fixed).
                        False → fitting path: use smart fit-time cache that
                                 reuses G when only scale/bg/SF params changed.
            pop_index:  Population index (0-based) used to key the fit-time
                        cache.  Must be supplied when use_cache=False to enable
                        smart reuse; ignored when use_cache=True.
        """
        return self._pop_intensity_and_grid(pop, q, use_cache, pop_index)[:3]

    def _pop_intensity_and_grid(
        self, pop: SizeDistPopulation, q: np.ndarray,
        use_cache: bool = True,
        pop_index: Optional[int] = None,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """As :meth:`calculate_pop_intensity`, plus the radius grid it used.

        ``total_intensity`` needs the grid too, and rebuilding it there cost as
        much as the intensity itself.
        """
        radius_grid, dr = self._radius_grid_and_widths(pop)
        vol_dist, num_dist = self.build_distributions(pop, radius_grid)

        if use_cache:
            # "Graph Model": use full-hash cache (params won't change mid-call)
            G = self.build_g_matrix(pop, q, radius_grid, use_cache=True)
        elif pop_index is not None:
            # Fitting: rebuild G only when its inputs (q / r_grid / FF /
            # contrast) actually changed; reuse it when only scale/bg/SF params
            # moved. q is compared by identity rather than hashed — a fit holds
            # its grids alive throughout, and with slit smearing on it evaluates
            # on two of them (the data grid and the smearer's extended grid).
            fit_key = self._fit_g_cache_key(radius_grid, pop)
            cached_key = self._fit_g_keys.get(pop_index)
            if cached_key is not None and cached_key[1] is q and cached_key[0] == fit_key:
                G = self._fit_g_cache[pop_index]
            else:
                G = build_g_matrix(
                    q=q, r_grid=radius_grid,
                    shape=pop.form_factor, contrast=pop.contrast,
                    **pop.ff_params,
                )
                self._fit_g_cache[pop_index] = G
                self._fit_g_keys[pop_index]  = (fit_key, q)
        else:
            # Fallback: always build fresh (no caching at all)
            G = build_g_matrix(
                q=q, r_grid=radius_grid,
                shape=pop.form_factor, contrast=pop.contrast,
                **pop.ff_params,
            )

        I_raw = G @ (vol_dist * dr)   # [cm⁻¹]
        I_pop = apply_structure_factor(I_raw, q, pop)
        return I_pop, vol_dist, num_dist, radius_grid

    # ── Total intensity ──────────────────────────────────────────────────────

    def total_intensity(
        self, config: ModelingConfig, q: np.ndarray,
        use_cache: bool = True,
        strict: bool = False,
    ) -> tuple[np.ndarray, list, list, list]:
        """Calculate summed model intensity over all enabled populations.

        Returns (I_total, pop_indices, pop_I_list, pop_dist_list)
        where pop_dist_list is list of (radius_grid, vol_dist, num_dist).

        Args:
            strict: what to do when a population cannot be evaluated.
                False (default, display paths) warns and leaves that population
                out of the sum, so a half-typed parameter cannot take a panel
                down. True (fitting) raises
                :class:`PopulationEvaluationError` instead — a silently dropped
                population changes chi-squared, and the fit would go on to
                report success for a model missing one of its components.
        """
        I_total = np.zeros(len(q), dtype=float)
        pop_indices = []
        pop_I_list = []
        pop_dist_list = []

        for i, pop in enumerate(config.populations):
            if not pop.enabled:
                continue
            try:
                pop_type = getattr(pop, 'pop_type', 'size_dist')
                if pop_type == 'unified_level':
                    I_pop = _unified_level_intensity(q, pop)
                    rg, vd, nd = None, None, None
                elif pop_type == 'diffraction_peak':
                    I_pop = _diffraction_peak_intensity(q, pop)
                    rg, vd, nd = None, None, None
                elif pop_type == 'guinier_porod':
                    I_pop = _guinier_porod_intensity(q, pop)
                    rg, vd, nd = None, None, None
                elif pop_type == 'mass_fractal':
                    I_pop = _mass_fractal_intensity(q, pop)
                    rg, vd, nd = None, None, None
                elif pop_type == 'surface_fractal':
                    I_pop = _surface_fractal_intensity(q, pop)
                    rg, vd, nd = None, None, None
                else:  # size_dist
                    # Fitting path only: reuse this population's curve outright
                    # when the optimiser step moved none of its parameters.
                    memo_key = None
                    hit = None
                    if not use_cache:
                        memo_key = self._pop_memo_key(pop)
                        cached = self._pop_memo.get(i)
                        # The entry keeps its own reference to q, so identity is
                        # a sound test — the smeared path evaluates on q_ext and
                        # must not collide with the data grid.
                        if (cached is not None
                                and cached[1] is q
                                and cached[0] == memo_key):
                            hit = cached
                    if hit is not None:
                        _, _, rg, vd, nd, I_pop = hit
                    else:
                        I_pop, vd, nd, rg = self._pop_intensity_and_grid(
                            pop, q,
                            use_cache=use_cache,
                            pop_index=(i if not use_cache else None),
                        )
                        if memo_key is not None:
                            self._pop_memo[i] = (memo_key, q, rg, vd, nd, I_pop)
            except Exception as e:
                if strict:
                    raise PopulationEvaluationError(
                        f"Population {i+1} "
                        f"({getattr(pop, 'pop_type', 'size_dist')}) failed: {e}"
                    ) from e
                warnings.warn(f"Population {i+1} failed: {e}", RuntimeWarning)
                continue

            I_total += I_pop
            pop_indices.append(i)
            pop_I_list.append(I_pop)
            pop_dist_list.append((rg, vd, nd))

        I_total += config.background
        return I_total, pop_indices, pop_I_list, pop_dist_list

    # ── Parameter packing / unpacking ────────────────────────────────────────

    def _pack_params(self, config: ModelingConfig) -> tuple[list, list, list, list]:
        """Flatten fittable parameters into a 1-D vector.

        Returns (x0, lo, hi, keys) where keys is a list of
        ('pop', pop_idx, 'dist'|'ff'|'sf'|'contrast'|'scale', param_name)
        or ('background',).
        """
        x0, lo, hi, keys = [], [], [], []

        for i, pop in enumerate(config.populations):
            if not pop.enabled:
                continue

            pop_type = getattr(pop, 'pop_type', 'size_dist')

            if pop_type == 'unified_level':
                for name in ['G', 'Rg', 'B', 'P', 'RgCO']:
                    if getattr(pop, f'fit_{name}', False):
                        lim = getattr(pop, f'{name}_limits')
                        x0.append(getattr(pop, name))
                        lo.append(lim[0]); hi.append(lim[1])
                        keys.append(('pop', i, 'uf', name))
                if pop.correlations:
                    for name in ['ETA', 'PACK']:
                        if getattr(pop, f'fit_{name}', False):
                            lim = getattr(pop, f'{name}_limits')
                            x0.append(getattr(pop, name))
                            lo.append(lim[0]); hi.append(lim[1])
                            keys.append(('pop', i, 'uf', name))
                continue

            if pop_type == 'diffraction_peak':
                for name in ['position', 'amplitude', 'width']:
                    if getattr(pop, f'fit_{name}', False):
                        lim = getattr(pop, f'{name}_limits')
                        x0.append(getattr(pop, name))
                        lo.append(lim[0]); hi.append(lim[1])
                        keys.append(('pop', i, 'peak', name))
                if pop.peak_type == 'voigt' and pop.fit_eta_voigt:
                    lim = pop.eta_voigt_limits
                    x0.append(pop.eta_voigt)
                    lo.append(lim[0]); hi.append(lim[1])
                    keys.append(('pop', i, 'peak', 'eta_voigt'))
                continue

            if pop_type == 'guinier_porod':
                for name in ['G', 'Rg1', 's1', 'P', 'Rg2', 's2', 'RgCO']:
                    if getattr(pop, f'fit_{name}', False):
                        lim = getattr(pop, f'{name}_limits')
                        x0.append(getattr(pop, name))
                        lo.append(lim[0]); hi.append(lim[1])
                        keys.append(('pop', i, 'gp', name))
                if pop.correlations:
                    for name in ['ETA', 'PACK']:
                        if getattr(pop, f'fit_{name}', False):
                            lim = getattr(pop, f'{name}_limits')
                            x0.append(getattr(pop, name))
                            lo.append(lim[0]); hi.append(lim[1])
                            keys.append(('pop', i, 'gp', name))
                continue

            if pop_type == 'mass_fractal':
                for name in ['Phi', 'Radius', 'Beta', 'Dv', 'Ksi', 'Eta', 'Contrast']:
                    if getattr(pop, f'fit_{name}', False):
                        lim = getattr(pop, f'{name}_limits')
                        x0.append(getattr(pop, name))
                        lo.append(lim[0]); hi.append(lim[1])
                        keys.append(('pop', i, 'mf', name))
                continue

            if pop_type == 'surface_fractal':
                for name in ['Surface', 'Ds', 'Ksi', 'Contrast']:
                    if getattr(pop, f'fit_{name}', False):
                        lim = getattr(pop, f'{name}_limits')
                        x0.append(getattr(pop, name))
                        lo.append(lim[0]); hi.append(lim[1])
                        keys.append(('pop', i, 'sf', name))
                if pop.use_porod_transition:
                    for name in ['Qc', 'QcWidth']:
                        if getattr(pop, f'fit_{name}', False):
                            lim = getattr(pop, f'{name}_limits')
                            x0.append(getattr(pop, name))
                            lo.append(lim[0]); hi.append(lim[1])
                            keys.append(('pop', i, 'sf', name))
                continue

            # size_dist — Distribution params
            for name in D.DIST_PARAM_NAMES.get(pop.dist_type, []):
                if pop.dist_params_fit.get(name, False):
                    val = float(pop.dist_params.get(name, 1.0))
                    lim = pop.dist_params_limits.get(name, (0.0, 1e10))
                    x0.append(val); lo.append(lim[0]); hi.append(lim[1])
                    keys.append(('pop', i, 'dist', name))

            # Form-factor params
            for name, val in pop.ff_params.items():
                if pop.ff_params_fit.get(name, False):
                    lim = pop.ff_params_limits.get(name, (1e-6, 1e6))
                    x0.append(float(val)); lo.append(lim[0]); hi.append(lim[1])
                    keys.append(('pop', i, 'ff', name))

            # Structure-factor params
            sf = pop.structure_factor.lower()
            sf_active_keys = {
                'interferences': ['eta', 'pack'],
                'hard_sphere':   ['radius', 'volume_fraction'],
            }.get(sf, [])
            for name in sf_active_keys:
                if pop.sf_params_fit.get(name, False):
                    val = float(pop.sf_params.get(name, 1.0))
                    lim = pop.sf_params_limits.get(name, (0.0, 1e10))
                    x0.append(val); lo.append(lim[0]); hi.append(lim[1])
                    keys.append(('pop', i, 'sf', name))

            # Scale
            if pop.fit_scale:
                lim = pop.scale_limits
                x0.append(pop.scale); lo.append(lim[0]); hi.append(lim[1])
                keys.append(('pop', i, 'scale', 'scale'))

            # Contrast
            if pop.fit_contrast:
                lim = pop.contrast_limits
                x0.append(pop.contrast); lo.append(lim[0]); hi.append(lim[1])
                keys.append(('pop', i, 'contrast', 'contrast'))

        # Background
        if config.fit_background:
            lim = config.background_limits
            x0.append(config.background); lo.append(lim[0]); hi.append(lim[1])
            keys.append(('background',))

        return x0, lo, hi, keys

    def _unpack_params(
        self, x: np.ndarray, keys: list, config: ModelingConfig
    ) -> None:
        """Write flat vector *x* back into *config* in-place."""
        for val, key in zip(x, keys):
            if key[0] == 'background':
                config.background = float(val)
            else:
                _, i, group, name = key
                pop = config.populations[i]
                if group in ('uf', 'gp', 'mf', 'sf', 'peak'):
                    setattr(pop, name, float(val))
                elif group == 'dist':
                    pop.dist_params[name] = float(val)
                elif group == 'ff':
                    pop.ff_params[name] = float(val)
                elif group == 'sf':
                    pop.sf_params[name] = float(val)
                elif group == 'scale':
                    pop.scale = float(val)
                elif group == 'contrast':
                    pop.contrast = float(val)

    # ── Fitting ──────────────────────────────────────────────────────────────

    def _residuals(
        self, x: np.ndarray, keys: list, config: ModelingConfig,
        q: np.ndarray, I: np.ndarray, sigma: np.ndarray,
    ) -> np.ndarray:
        """Weighted residuals for least_squares.

        Evaluated strictly: a population that cannot be computed must never be
        quietly dropped from a trial, because the resulting chi-squared would be
        that of a smaller model and the optimiser would be free to prefer it.
        Instead the trial is scored as maximally bad, so the optimiser walks
        away from that corner of parameter space, and the failure is counted for
        :meth:`fit` to report. A model that fails at the *starting* parameters
        is a different matter — ``fit`` checks that up front and raises.
        """
        self._unpack_params(x, keys, config)
        try:
            I_model, _, _, _ = self.total_intensity_maybe_smeared(
                config, q, use_cache=False, strict=True)
        except PopulationEvaluationError as exc:
            self._eval_failures += 1
            self._eval_failure_last = str(exc)
            return np.full(len(q), _FAILED_TRIAL_RESIDUAL, dtype=float)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            resid = (I - I_model) / np.maximum(sigma, 1e-30)
        return np.where(np.isfinite(resid), resid, 0.0)

    def _chi2(
        self, x: np.ndarray, keys: list, config: ModelingConfig,
        q: np.ndarray, I: np.ndarray, sigma: np.ndarray,
    ) -> float:
        """Scalar chi² for Nelder-Mead (no-limits mode) and global search."""
        resid = self._residuals(x, keys, config, q, I, sigma)
        return float(np.sum(resid ** 2))

    def _run_global_de(
        self, lo_arr: np.ndarray, hi_arr: np.ndarray, keys: list,
        cfg: ModelingConfig, q_fit: np.ndarray, I_fit: np.ndarray,
        dI_fit: np.ndarray, x0_arr: np.ndarray, seed=None,
        workers: int = 1, cancel_cb=None,
    ) -> np.ndarray:
        """Differential-evolution global search; returns the best parameter vector.

        This locates the correct basin of a multimodal χ² surface (e.g. the
        right Bessel lobe for a monodisperse core-shell radius); the caller
        then polishes the result with a local TRF least-squares step.

        Parameters that span many decades (lo > 0 and hi/lo > 100) are searched
        in log10 space so DE samples small and large values evenly. Without
        this, wide bounds such as scale ∈ [0, 1e10] or background ∈ [0, 1e10]
        make uniform DE sampling almost never probe physically relevant small
        values. The transform is internal to this method — it never touches
        _pack_params / _unpack_params or the stored parameter values.

        Cancellation:
          * Serial (workers == 1): the GUI worker wraps ``_chi2``, so a "Cancel
            Fit" raises out of each evaluation, as before.
          * Parallel (workers != 1): the wrapped ``_chi2`` cannot be pickled to
            worker processes, so a picklable objective on a clean engine is used
            instead and cancellation is polled via the DE ``callback`` (which
            runs in the main process once per generation). ``cancel_cb`` returns
            True to stop. Any multiprocessing failure falls back to a serial run.
        """
        lo_arr = np.asarray(lo_arr, dtype=float)
        hi_arr = np.asarray(hi_arr, dtype=float)

        # Decide which parameters to search in log10 space.
        with np.errstate(divide='ignore', invalid='ignore'):
            ratio = np.where(lo_arr > 0, hi_arr / np.maximum(lo_arr, 1e-300), 1.0)
        log_mask = (lo_arr > 0) & (ratio > 100.0)

        def _to_internal(x):
            t = np.array(x, dtype=float)
            t[log_mask] = np.log10(np.clip(t[log_mask], 1e-300, None))
            return t

        def _to_real(t):
            x = np.array(t, dtype=float)
            x[log_mask] = 10.0 ** x[log_mask]
            return x

        t_lo = _to_internal(lo_arr)
        t_hi = _to_internal(hi_arr)
        t_bounds = list(zip(t_lo, t_hi))

        # Seed DE with the current guess (clamped into the transformed box) so
        # the global search can never return something worse than the start.
        t_x0 = np.clip(_to_internal(x0_arr), t_lo, t_hi)

        common = dict(
            strategy='best1bin', popsize=15, maxiter=300, tol=1e-2,
            mutation=(0.5, 1.0), recombination=0.7, init='latinhypercube',
            polish=False, x0=t_x0, seed=seed,
        )

        def _serial_de():
            # Serial objective uses self._chi2 (possibly the GUI-wrapped one),
            # so per-evaluation cancellation keeps working. scipy appends the
            # tuple passed via `args` to each call: _de_obj(x, *args).
            def _de_obj(t_x, *args):
                return self._chi2(_to_real(t_x), *args)
            return differential_evolution(
                _de_obj, t_bounds,
                args=(keys, cfg, q_fit, I_fit, dI_fit),
                updating='immediate', workers=1, **common,
            )

        if workers == 1:
            de = _serial_de()
        else:
            # Parallel: build a picklable objective on a *clean* engine (the
            # live engine may carry a GUI cancellation monkey-patch that cannot
            # be pickled), and poll cancellation via the per-generation callback.
            objective = _DEObjective(
                ModelingEngine(), keys, deepcopy(cfg),
                q_fit, I_fit, dI_fit, log_mask,
            )

            # Legacy-style callback (no `intermediate_result` param) → scipy
            # calls it as (xk, convergence=…) and stops when it returns True.
            def _cb(*_a, **_k):
                return bool(cancel_cb()) if cancel_cb is not None else False

            try:
                # _DEObjective carries keys/cfg/data itself, so no `args` here.
                with _limit_blas_threads(1):
                    de = differential_evolution(
                        objective, t_bounds,
                        updating='deferred', workers=workers, callback=_cb,
                        **common,
                    )
            except Exception as exc:
                # Multiprocessing can fail on some hosts (spawn/pickling/frozen
                # apps). Never fail the fit for it — fall back to serial.
                warnings.warn(
                    f"Parallel global fit failed ({exc!r}); running serially.",
                    RuntimeWarning,
                )
                de = _serial_de()

        return np.clip(_to_real(de.x), lo_arr, hi_arr)

    def fit(
        self,
        config: ModelingConfig,
        q: np.ndarray,
        I: np.ndarray,
        dI: np.ndarray,
        cancel_cb=None,
    ) -> ModelingResult:
        """Fit the multi-population model to SAS data.

        The q-range is cropped to [config.q_min, config.q_max] before fitting.

        Args:
            config: ModelingConfig (modified in-place with best-fit params).
            q:      Scattering vector [Å⁻¹].
            I:      Measured intensity [cm⁻¹].
            dI:     Measurement uncertainty (σ) [cm⁻¹].
            cancel_cb: optional callable returning True to abort. Used by the
                parallel global search (workers > 1) to poll for "Cancel Fit"
                via the DE per-generation callback (the serial path cancels
                through the GUI's per-evaluation _chi2 wrapper instead).

        Returns:
            ModelingResult containing fitted parameters and derived quantities.

        Raises:
            PopulationEvaluationError: the model cannot be evaluated at the
                starting parameters. Continuing would fit a model with that
                population missing and report success for it.
        """
        # Fresh fit — discard anything cached from a previous fit() call
        self._fit_g_cache.clear()
        self._fit_g_keys.clear()
        self._rg_cache.clear()
        self._pop_memo.clear()
        self._eval_failures = 0
        self._eval_failure_last = ''

        # Crop q range
        mask = (q >= config.q_min) & (q <= config.q_max)
        q_fit = q[mask]
        I_fit = I[mask]
        dI_fit = np.maximum(dI[mask], 1e-30)

        if len(q_fit) == 0:
            raise ValueError("No data points within [q_min, q_max].")

        cfg = deepcopy(config)
        x0, lo, hi, keys = self._pack_params(cfg)

        # Fail fast on a model that cannot be evaluated at all. Without this a
        # broken population is dropped from every trial alike, so the optimiser
        # converges happily and the fit reports success for a model that is
        # quietly missing one of its components.
        self.total_intensity_maybe_smeared(cfg, q_fit, use_cache=False, strict=True)

        # Plain-language diagnostics collected during the fit and returned on
        # the result (GUI shows them in the status line, batch logs them).
        fit_warn: list = []
        converged = True

        if len(x0) == 0:
            # Nothing to fit — just evaluate
            I_model, pop_idx, pop_I, pop_dist = self.total_intensity_maybe_smeared(
                cfg, q_fit, use_cache=True, strict=True
            )
            resid = (I_fit - I_model) / dI_fit
            chi2 = float(np.sum(resid ** 2))
            dof = max(len(q_fit) - 1, 1)
        else:
            x0_arr = np.array(x0, dtype=float)

            if config.no_limits:
                result = minimize(
                    self._chi2, x0_arr,
                    args=(keys, cfg, q_fit, I_fit, dI_fit),
                    method='Nelder-Mead',
                    options={'maxiter': 500, 'xatol': 1e-4, 'fatol': 1e-4},
                )
                x_best = result.x
                converged = bool(result.success)
            else:
                lo_arr = np.array(lo, dtype=float)
                hi_arr = np.array(hi, dtype=float)

                if config.fit_method == 'global':
                    # Global search locates the basin; TRF then polishes it to
                    # the exact minimum (and yields the Jacobian for errors).
                    x_seed = self._run_global_de(
                        lo_arr, hi_arr, keys, cfg,
                        q_fit, I_fit, dI_fit, x0_arr,
                        workers=getattr(config, 'de_workers', 1),
                        cancel_cb=cancel_cb,
                    )
                else:
                    x_seed = x0_arr

                # x_scale='jac' auto-rescales each parameter by its Jacobian-
                # column norm every iteration. Modeling parameters span many
                # decades (unified-level G ~ 10³, B ~ 10⁻⁴, Rg ~ 10¹; scale
                # ~ 10⁻³; background ~ 10⁻²). Without it the TRF trust region
                # and the xtol/ftol tests act on the raw vector and terminate
                # far short of the minimum. Auto-scaling gives every parameter a
                # natural unit step (the scipy equivalent of Igor's
                # per-parameter fit-step on log-dependent parameters).
                #
                # Tight ftol/xtol/gtol (1e-12): with 'jac' scaling the
                # convergence tests run in SCALED space, where loose tolerances
                # trigger a spurious "converged" on the first small step while
                # still far from the minimum — the cause of the "press Fit 2–3
                # times" behaviour.
                #
                # Internal restart loop: re-seed the solver from its own result
                # until χ² stops improving, so scripts (which cannot re-press
                # Fit) get the fully-settled result on the first call. The
                # Modeling GUI worker wraps self._residuals to raise on each
                # evaluation for "Cancel Fit"; because every restart calls
                # self._residuals, cancellation stays responsive across the loop.
                #
                # The restart threshold is RELATIVE (1e-4·χ²). A model that
                # structurally cannot reach the data (e.g. a huge low-Q upturn
                # with `scale` pinned at its limit and a lognormal σ driven to
                # extremes) leaves TRF wandering a nearly flat χ² valley: every
                # round exhausts max_nfev, and with χ² ~ 10⁶ even numerical
                # wiggles of a few counts exceeded the old 1e-8·χ² threshold —
                # so all 5 restarts ran, ≈ 10 minutes of fitting with nothing
                # to show for it (see testData/bad1.h5). Genuine "press Fit
                # again" improvements are orders of magnitude above 1e-4·χ²,
                # so well-behaved fits are unaffected.
                # diff_step=1e-3: finite-difference Jacobian step of 0.1% of
                # each parameter's value — the scipy equivalent of Igor's
                # per-parameter epsilon (starting step as a fraction of the
                # value), which Irena needed for exactly this reason. scipy's
                # default step (√eps ≈ 1.5e-8 relative) yields a locally exact
                # Jacobian, and on stiff, strongly correlated surfaces (a
                # mean_size of 5×10⁴ fitted next to peak widths of 4×10⁻³)
                # TRF then advances in micro-steps: ~2400 evaluations per
                # round without meeting any tolerance. A 1e-3 secant slope
                # probes the surface at the scale of steps worth taking; the
                # same fits then converge in ~30–90 evaluations to an equal or
                # better χ² (verified on testData/bad1.h5 and Example2:
                # ~20–35× fewer evaluations).
                ls_common = dict(
                    args=(keys, cfg, q_fit, I_fit, dI_fit),
                    bounds=(lo_arr, hi_arr),
                    method='trf',
                    x_scale='jac',
                    diff_step=1e-3,
                    max_nfev=300,
                    ftol=1e-12, xtol=1e-12, gtol=1e-12,
                )
                prev_chi2 = np.inf
                max_restarts = 5
                result = least_squares(self._residuals, x_seed, **ls_common)
                for _ in range(max_restarts - 1):
                    chi2 = float(np.sum(result.fun ** 2))
                    if prev_chi2 - chi2 <= 1e-4 * max(chi2, 1.0):
                        break
                    prev_chi2 = chi2
                    result = least_squares(
                        self._residuals, result.x, **ls_common)
                x_best = result.x
                # status 0 = max_nfev exhausted without meeting any tolerance
                converged = (result.status != 0)
                fit_warn.extend(
                    _pinned_parameter_warnings(x_best, lo_arr, hi_arr, keys))

            self._unpack_params(x_best, keys, cfg)

            # Write best-fit params back to original config
            self._unpack_params(x_best, keys, config)

            I_model, pop_idx, pop_I, pop_dist = self.total_intensity_maybe_smeared(
                cfg, q_fit, use_cache=False, strict=True
            )
            resid = (I_fit - I_model) / dI_fit
            chi2 = float(np.sum(resid ** 2))
            dof = max(len(q_fit) - len(x0), 1)

            if self._eval_failures:
                # Counted in this process only: the parallel global search
                # evaluates on worker engines, which are still strict but keep
                # their own tallies.
                fit_warn.append(
                    f"{self._eval_failures} trial parameter set(s) were rejected "
                    f"because a population could not be evaluated there "
                    f"(last: {self._eval_failure_last}) — the reported fit does "
                    f"not use them, but check the parameter limits"
                )

        rchi2 = chi2 / dof

        if not converged:
            fit_warn.append(
                "fit stopped at the evaluation budget without converging — "
                "the parameters are the best found so far, not a settled minimum"
            )
        fit_warn.extend(_model_sanity_warnings(cfg, rchi2))

        # Ideal (pinhole) total for saving — None unless smearing was used.
        model_I_ideal = self.total_intensity_ideal(cfg, q_fit)

        # ── Derived quantities per population ─────────────────────────────
        derived = []
        radius_grids, vol_dists, num_dists = [], [], []
        for (rg, vd, nd), pi in zip(pop_dist, pop_idx):
            radius_grids.append(rg)
            vol_dists.append(vd)
            num_dists.append(nd)
            derived.append(compute_derived(rg, vd, nd, cfg.populations[pi]))

        return ModelingResult(
            config=deepcopy(cfg),
            chi_squared=chi2,
            reduced_chi_squared=rchi2,
            dof=dof,
            timestamp=datetime.now().isoformat(timespec='seconds'),
            model_q=q_fit,
            model_I=I_model,
            pop_indices=pop_idx,
            pop_model_I=pop_I,
            radius_grids=radius_grids,
            volume_dists=vol_dists,
            number_dists=num_dists,
            derived=derived,
            model_I_ideal=model_I_ideal,
            fit_warnings=fit_warn,
        )

    # ── Monte-Carlo uncertainty ──────────────────────────────────────────────

    def calculate_uncertainty_mc(
        self,
        config: ModelingConfig,
        q: np.ndarray,
        I: np.ndarray,
        dI: np.ndarray,
        n_runs: int = 10,
        progress_cb=None,
        cancel_cb=None,
        workers: Optional[int] = None,
        seed=None,
    ) -> dict[str, float]:
        """Estimate parameter uncertainty by repeated fitting on perturbed data.

        Data is perturbed as I_perturbed = I + dI * N(0,1) for each run.
        Returns a dict of std-dev for each fittable parameter.

        Each pass is an independent refit, so passes are spread over worker
        processes when that is worth doing (see *workers*). All the noise is
        drawn in *this* process before any pass starts, and workers receive
        finished arrays rather than a seed. That keeps parallel results
        bit-identical to serial ones for a given *seed*, and rules out the
        failure mode where forked workers inherit one RNG state and return N
        identical passes — which would silently drive every uncertainty to zero.

        Args:
            progress_cb: optional callable(n_done, n_runs) called as passes
                complete. Passes finish out of order when running in parallel,
                so this counts completions rather than reporting a pass index.
            cancel_cb: optional callable returning True to stop early. Passes
                already in flight are not interrupted, so a cancel takes effect
                within roughly one pass; uncertainties from the passes that did
                finish are still returned.
            workers: process-count override. None (default) uses
                ``config.mc_workers``: 0 = auto (cpu_count - 2), 1 = serial,
                N > 1 = explicit. Falls back to serial if the pool cannot start.
            seed: seed for the perturbation RNG. None (default) keeps passes
                nondeterministic; pass an int for reproducible uncertainties.
        """
        # Collect best-fit values per run
        _, _, _, keys = self._pack_params(deepcopy(config))
        if not keys:
            return {}

        n_runs = int(n_runs)
        if n_runs < 1:
            return {}

        # MC perturbs data around the already-found solution, so a fast local
        # refinement per pass is appropriate — running the global search n_runs
        # times would be needlessly slow. de_workers is pinned for the same
        # reason: with fit_method='local' the DE pool is already unreachable,
        # and pinning it stops a future edit from nesting pools inside this one.
        cfg_mc = deepcopy(config)
        cfg_mc.fit_method = 'local'
        cfg_mc.de_workers = 1

        rng = np.random.default_rng(seed)
        perturbed = [I + dI * rng.standard_normal(len(I)) for _ in range(n_runs)]

        run = _MCRun(cfg_mc, q, dI, len(keys))
        all_vals: list[list[float]] = [[] for _ in keys]
        n_done = 0

        def _record(x_final):
            nonlocal n_done
            n_done += 1
            if x_final is not None:
                for j, v in enumerate(x_final):
                    all_vals[j].append(v)
            if progress_cb is not None:
                progress_cb(n_done, n_runs)

        # Pass 1 runs here so its cost can be measured: it decides whether the
        # rest are worth a process pool at all. No work is wasted — it counts
        # towards n_runs like any other pass — but it does serialise one pass
        # ahead of the pool. That costs at most one pass of wall time, which is
        # a shrinking fraction as n_runs grows past the worker count, and it is
        # what lets cheap models skip the pool entirely.
        t0 = time.perf_counter()
        _record(run(perturbed[0]))
        first_pass_s = time.perf_counter() - t0

        remaining = perturbed[1:]
        n_workers = _resolve_mc_workers(
            workers if workers is not None else getattr(config, 'mc_workers', 0)
        )
        use_pool = (
            bool(remaining)
            and n_workers > 1
            and first_pass_s * len(remaining) >= _MC_PARALLEL_MIN_SECONDS
        )

        if remaining and not _mc_cancelled(cancel_cb):
            leftovers = remaining
            if use_pool:
                try:
                    leftovers = _run_mc_parallel(
                        run, remaining, n_workers, _record, cancel_cb)
                except Exception as exc:
                    # Multiprocessing can fail on some hosts (spawn/pickling/
                    # frozen apps). Never fail the estimate for it — the pool
                    # raises before recording anything, so a serial retry of the
                    # whole batch cannot double-count.
                    warnings.warn(
                        f"Parallel Monte-Carlo failed ({exc!r}); running serially.",
                        RuntimeWarning,
                    )
                    leftovers = remaining
            for I_pert in leftovers:
                if _mc_cancelled(cancel_cb):
                    break
                _record(run(I_pert))

        if len(all_vals[0]) < 2:
            return {}

        stds: dict[str, float] = {}
        for key, vals in zip(keys, all_vals):
            if len(vals) < 2:
                continue
            label = _key_label(key)
            stds[label] = float(np.std(vals, ddof=1))

        return stds


def _key_label(key: tuple) -> str:
    """Human-readable label for a parameter key tuple."""
    if key[0] == 'background':
        return 'background'
    _, i, group, name = key
    return f'pop{i+1}_{group}_{name}'
