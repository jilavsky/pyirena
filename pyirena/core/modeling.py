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

import warnings
import hashlib
import json
from copy import deepcopy
from dataclasses import dataclass, field
from datetime import datetime
from typing import Any, Optional

import numpy as np
from scipy.optimize import least_squares, minimize
from scipy.special import erf as _scipy_erf

from pyirena.core import distributions as D
from pyirena.core.form_factors import build_g_matrix, bin_widths


# ──────────────────────────────────────────────────────────────────────────────
# Dataclasses
# ──────────────────────────────────────────────────────────────────────────────

@dataclass
class SizeDistPopulation:
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
class UnifiedLevelPopulation:
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
class DiffractionPeakPopulation:
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


@dataclass
class ModelingResult:
    """Results from a Modeling fit."""
    config: ModelingConfig
    chi_squared: float
    reduced_chi_squared: float
    dof: int
    timestamp: str

    model_q: np.ndarray              # Q-array used for fitting
    model_I: np.ndarray              # Total model intensity

    # Per-population arrays (only enabled populations; order matches enabled pops)
    pop_indices: list                # list of int — original population indices
    pop_model_I: list                # list of np.ndarray
    radius_grids: list               # list of np.ndarray   [Å]
    volume_dists: list               # list of np.ndarray   [Å⁻¹ normalised]
    number_dists: list               # list of np.ndarray   [Å⁻¹ normalised]
    derived: list                    # list of dict per active population

    # MC uncertainty: {pop_i_paramname: std} and {'background': std}
    params_std: dict = field(default_factory=dict)


# ──────────────────────────────────────────────────────────────────────────────
# Structure-factor helpers  (stand-alone, no class state needed)
# ──────────────────────────────────────────────────────────────────────────────

def _sphere_amplitude(q: np.ndarray, eta: float) -> np.ndarray:
    """Born-Green sphere amplitude  f(q,η) = 3[sin(qη) - qη·cos(qη)] / (qη)³."""
    q_eta = q * eta
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return np.where(
            q_eta < 1e-10,
            1.0,
            3.0 * (np.sin(q_eta) - q_eta * np.cos(q_eta)) / (q_eta ** 3),
        )


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
# Derived-quantity helpers
# ──────────────────────────────────────────────────────────────────────────────

def _volume_sphere(r: np.ndarray) -> np.ndarray:
    return (4.0 / 3.0) * np.pi * r ** 3


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
    invariant = float(np.trapz(I * q ** 2, q))
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

    # Rg² (volume-weighted): ∫ (3/5)*r² * vol_dist dr / vol_tot
    if vol_tot > 0:
        Rg2 = np.sum((3.0 / 5.0) * radius_grid ** 2 * vol_dist * dr) / vol_tot
        Rg = np.sqrt(max(Rg2, 0.0))
    else:
        Rg = 0.0

    # Specific surface (sphere): ∫ (3/r) * vol_dist dr / vol_tot  [Å⁻¹]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        safe_r = np.maximum(radius_grid, 1e-30)
        if vol_tot > 0:
            specific_surface = np.sum((3.0 / safe_r) * vol_dist * dr) / vol_tot
        else:
            specific_surface = 0.0

    return {
        'vol_mean_r':       vol_mean_r,
        'num_mean_r':       num_mean_r,
        'volume_fraction':  vf,
        'Rg':               Rg,
        'specific_surface': specific_surface,
    }


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
        self._fit_g_keys:  dict[int, str] = {}

    # ── Radius grid ──────────────────────────────────────────────────────────

    def build_radius_grid(self, pop: SizeDistPopulation) -> np.ndarray:
        """Return the non-linear radius grid for *pop*."""
        return D.generate_radius_grid(
            pop.dist_type, pop.dist_params,
            n_bins=pop.n_bins, tail_precision=0.01,
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

        V_r = _volume_sphere(radius_grid)   # particle volume per bin  [Å³]

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
        """Discard all cached G matrices."""
        self._g_cache.clear()

    @staticmethod
    def _fit_g_cache_key(r_grid: np.ndarray, pop: SizeDistPopulation) -> str:
        """Cache key for the fit-time smart G cache.

        Keyed only on what G actually depends on — the radius grid (as a hash
        of its bytes), form factor name, ff_params, and contrast.  Deliberately
        excludes scale, background, and sf_params so that optimizer steps that
        only touch those parameters can reuse the cached G matrix.
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
        radius_grid = self.build_radius_grid(pop)
        vol_dist, num_dist = self.build_distributions(pop, radius_grid)

        if use_cache:
            # "Graph Model": use full-hash cache (params won't change mid-call)
            G = self.build_g_matrix(pop, q, radius_grid, use_cache=True)
        elif pop_index is not None:
            # Fitting: rebuild G only when its inputs (r_grid / FF / contrast)
            # actually changed; reuse it when only scale/bg/SF params moved.
            fit_key = self._fit_g_cache_key(radius_grid, pop)
            if self._fit_g_keys.get(pop_index) == fit_key:
                G = self._fit_g_cache[pop_index]
            else:
                G = build_g_matrix(
                    q=q, r_grid=radius_grid,
                    shape=pop.form_factor, contrast=pop.contrast,
                    **pop.ff_params,
                )
                self._fit_g_cache[pop_index] = G
                self._fit_g_keys[pop_index]  = fit_key
        else:
            # Fallback: always build fresh (no caching at all)
            G = build_g_matrix(
                q=q, r_grid=radius_grid,
                shape=pop.form_factor, contrast=pop.contrast,
                **pop.ff_params,
            )

        dr = bin_widths(radius_grid)
        I_raw = G @ (vol_dist * dr)   # [cm⁻¹]
        I_pop = apply_structure_factor(I_raw, q, pop)
        return I_pop, vol_dist, num_dist

    # ── Total intensity ──────────────────────────────────────────────────────

    def total_intensity(
        self, config: ModelingConfig, q: np.ndarray,
        use_cache: bool = True,
    ) -> tuple[np.ndarray, list, list, list]:
        """Calculate summed model intensity over all enabled populations.

        Returns (I_total, pop_indices, pop_I_list, pop_dist_list)
        where pop_dist_list is list of (radius_grid, vol_dist, num_dist).
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
                else:  # size_dist
                    I_pop, vd, nd = self.calculate_pop_intensity(
                        pop, q,
                        use_cache=use_cache,
                        pop_index=(i if not use_cache else None),
                    )
                    rg = self.build_radius_grid(pop)
            except Exception as e:
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
                if group == 'uf':
                    setattr(pop, name, float(val))
                elif group == 'peak':
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
        """Weighted residuals for least_squares."""
        self._unpack_params(x, keys, config)
        I_model, _, _, _ = self.total_intensity(config, q, use_cache=False)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            resid = (I - I_model) / np.maximum(sigma, 1e-30)
        return np.where(np.isfinite(resid), resid, 0.0)

    def _chi2(
        self, x: np.ndarray, keys: list, config: ModelingConfig,
        q: np.ndarray, I: np.ndarray, sigma: np.ndarray,
    ) -> float:
        """Scalar chi² for Nelder-Mead (no-limits mode)."""
        resid = self._residuals(x, keys, config, q, I, sigma)
        return float(np.sum(resid ** 2))

    def fit(
        self,
        config: ModelingConfig,
        q: np.ndarray,
        I: np.ndarray,
        dI: np.ndarray,
    ) -> ModelingResult:
        """Fit the multi-population model to SAS data.

        The q-range is cropped to [config.q_min, config.q_max] before fitting.

        Args:
            config: ModelingConfig (modified in-place with best-fit params).
            q:      Scattering vector [Å⁻¹].
            I:      Measured intensity [cm⁻¹].
            dI:     Measurement uncertainty (σ) [cm⁻¹].

        Returns:
            ModelingResult containing fitted parameters and derived quantities.
        """
        # Fresh fit — discard any G matrices cached from a previous fit() call
        self._fit_g_cache.clear()
        self._fit_g_keys.clear()

        # Crop q range
        mask = (q >= config.q_min) & (q <= config.q_max)
        q_fit = q[mask]
        I_fit = I[mask]
        dI_fit = np.maximum(dI[mask], 1e-30)

        if len(q_fit) == 0:
            raise ValueError("No data points within [q_min, q_max].")

        cfg = deepcopy(config)
        x0, lo, hi, keys = self._pack_params(cfg)

        if len(x0) == 0:
            # Nothing to fit — just evaluate
            I_model, pop_idx, pop_I, pop_dist = self.total_intensity(
                cfg, q_fit, use_cache=True
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
                    options={'maxiter': 100_000, 'xatol': 1e-6, 'fatol': 1e-6},
                )
                x_best = result.x
            else:
                lo_arr = np.array(lo, dtype=float)
                hi_arr = np.array(hi, dtype=float)
                result = least_squares(
                    self._residuals, x0_arr,
                    args=(keys, cfg, q_fit, I_fit, dI_fit),
                    bounds=(lo_arr, hi_arr),
                    method='trf',
                    max_nfev=50_000,
                    ftol=1e-8, xtol=1e-8, gtol=1e-8,
                )
                x_best = result.x

            self._unpack_params(x_best, keys, cfg)

            # Write best-fit params back to original config
            self._unpack_params(x_best, keys, config)

            I_model, pop_idx, pop_I, pop_dist = self.total_intensity(
                cfg, q_fit, use_cache=False
            )
            resid = (I_fit - I_model) / dI_fit
            chi2 = float(np.sum(resid ** 2))
            dof = max(len(q_fit) - len(x0), 1)

        rchi2 = chi2 / dof

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
    ) -> dict[str, float]:
        """Estimate parameter uncertainty by repeated fitting on perturbed data.

        Data is perturbed as I_perturbed = I + dI * N(0,1) for each run.
        Returns a dict of std-dev for each fittable parameter.

        Args:
            progress_cb: optional callable(run_index, n_runs) called before each run.
        """
        rng = np.random.default_rng()

        # Collect best-fit values per run
        _, _, _, keys = self._pack_params(deepcopy(config))
        if not keys:
            return {}

        all_vals: list[list[float]] = [[] for _ in keys]
        good_runs = 0

        for run_i in range(n_runs):
            if progress_cb is not None:
                progress_cb(run_i + 1, n_runs)
            I_pert = I + dI * rng.standard_normal(len(I))
            cfg_run = deepcopy(config)
            try:
                result = self.fit(cfg_run, q, I_pert, dI)
                x_final, _, _, keys_run = self._pack_params(deepcopy(cfg_run))
                if len(x_final) == len(keys):
                    for j, v in enumerate(x_final):
                        all_vals[j].append(v)
                    good_runs += 1
            except Exception:
                continue

        if good_runs < 2:
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
