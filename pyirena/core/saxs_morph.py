"""
SAXS Morph engine — 3D two-phase voxelgram model via Gaussian Random Fields.

This module is a Python port of Igor Pro's IR3_3DTwoPhaseSolid.ipf (J. Ilavsky,
Argonne National Laboratory), itself based on:

  * Berk, N. F. (1991). Phys. Rev. E 51, 4141.
  * Roberts, A. P. (1997).
  * Levitz, P. (2007). Modelling Simul. Mater. Sci. Eng. 15 S2.
  * Cherny, A. Yu., Anitas, E. M., Kuklin, A. I., Balasoiu, M., Osipov, V.
    A. (2013). J. Appl. Cryst. 46, 365 — DOI 10.1107/S0021889813003816.

Algorithm
---------
For a desired volume fraction phi of phase 1 and a target small-angle
scattering profile I(Q):

  1. Subtract user background (Power-law + Flat) from the data:
        I_corr(Q) = I(Q) - (B * Q**-P + flat)
  2. Compute the Debye autocorrelation function gamma(r) by sinc transform of
     Q**2 * I_corr(Q); normalise gamma(0) = 1.
  3. Compute the spectral density F(k) by sinc transform of gamma(r) * r**2.
     Clip to F >= 0 (negative bins are numerical noise).
  4. Compute the threshold parameter alfa = sqrt(2) * erfinv(1 - 2*phi).
  5. Generate a 3D Gaussian white-noise cube W, FFT, multiply by sqrt(F(k_3d))
     resampled onto the 3D radial k-grid, inverse FFT -> real scalar field.
  6. Threshold the scalar field at alfa * sigma -> binary uint8 voxelgram.
  7. Compute the model intensity I_model(Q) by FFT of the voxelgram and
     spherical averaging of |F(k)|**2.
  8. Add background back; compare to data; iterate parameters.

Engine class
------------
SaxsMorphEngine has the same surface as ModelingEngine: ``compute_voxelgram``
for one-shot evaluation, ``fit`` for least_squares / Nelder-Mead, and
``calculate_uncertainty_mc`` for Gaussian-perturbation Monte Carlo.

Memory note
-----------
A 256**3 binary cube is 16 MB; 512**3 is 125 MB; the intermediate complex128
FFT is 16x larger.  ``fit`` hard-clamps voxel size to <= 256 to keep iteration
memory within reach of a typical laptop; only the final post-convergence
voxelgram is rendered at the user-selected ``voxel_size_render`` (up to 512).
"""

from __future__ import annotations

import warnings
from copy import deepcopy
from dataclasses import dataclass, field
from datetime import datetime
from typing import Optional

import numpy as np
from scipy.optimize import least_squares, minimize
from scipy.special import erfinv as _erfinv
from scipy.integrate import simpson as _simpson


# Hard memory ceiling for the per-iteration voxelgram during fitting.
# 256**3 complex128 FFT ~= 256 MB transient; safe on 8 GB laptops.
MAX_FIT_VOXEL_SIZE = 256

# Allowed cube sizes (cube only in v1).
ALLOWED_VOXEL_SIZES = (64, 128, 256, 384, 512)


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------

@dataclass
class SaxsMorphConfig:
    """Configuration for one SAXS Morph run.

    Q range
    -------
    q_min, q_max          Cursor-driven Q range used for fitting; outside this
                          window data are ignored.

    Voxel grid
    ----------
    voxel_size_fit        Cube side used during the fit loop (clamped <= 256).
    voxel_size_render     Cube side used for the final post-convergence
                          voxelgram (up to 512).
    box_size_A            Physical cube side length [A].

    Two-phase parameters
    --------------------
    volume_fraction       phi of phase 1 (0..1). Always fittable in principle,
                          but ``fit_volume_fraction`` enables it.
    contrast              (Delta rho)**2 in 10**20 cm**-4. Used as a global
                          intensity scale once the structure factor S(Q) is
                          computed from the voxelgram.
    link_phi_contrast     If True, contrast is derived inside ``compute_voxelgram``
                          from the data invariant Q* = 2 pi**2 * phi (1-phi) * Delta rho**2,
                          so the fit only varies phi (contrast is auto-set).

    Background
    ----------
    power_law_B, power_law_P, background    Power-law amplitude / exponent
                          and flat background; subtracted from data before the
                          GRF inversion and re-added when computing model I(Q).

    Fit flags / limits
    ------------------
    *_fit                 Per-parameter boolean: include in least-squares.
    *_limits              Per-parameter (lo, hi) tuple.

    Misc
    ----
    no_limits             If True, switch to Nelder-Mead (no bounds).
    n_mc_runs             Number of Monte Carlo passes for uncertainty.
    rng_seed              Optional integer seed for reproducible voxelgrams;
                          None -> fresh per call.
    """
    # Q range
    q_min: Optional[float] = None
    q_max: Optional[float] = None

    # Voxel grid
    voxel_size_fit: int = 128
    voxel_size_render: int = 256
    box_size_A: float = 1000.0

    # Two-phase parameters
    volume_fraction: float = 0.3
    fit_volume_fraction: bool = True
    volume_fraction_limits: tuple = (0.05, 0.95)

    contrast: float = 1.0
    fit_contrast: bool = False
    contrast_limits: tuple = (0.0, 1e10)

    link_phi_contrast: bool = True

    # Background (Power-law + Flat)
    power_law_B: float = 0.0
    fit_power_law_B: bool = False
    power_law_B_limits: tuple = (0.0, 1e10)

    power_law_P: float = 4.0
    fit_power_law_P: bool = False
    power_law_P_limits: tuple = (0.0, 6.0)

    background: float = 0.0
    fit_background: bool = False
    background_limits: tuple = (0.0, 1e10)

    # Misc
    no_limits: bool = False
    n_mc_runs: int = 10
    rng_seed: Optional[int] = None


@dataclass
class SaxsMorphResult:
    """Outputs of one SaxsMorphEngine.compute_voxelgram or .fit call."""
    config: SaxsMorphConfig
    chi_squared: float
    reduced_chi_squared: float
    dof: int
    timestamp: str

    # Q grid + intensities (over the fit window only)
    data_q: np.ndarray
    data_I: np.ndarray
    data_dI: np.ndarray
    data_I_corr: np.ndarray            # data minus user background
    model_q: np.ndarray
    model_I: np.ndarray                # full model = structure*contrast + bg

    # Intermediate functions (1D)
    r_grid: np.ndarray
    gamma_r: np.ndarray
    spectral_k: np.ndarray
    spectral_F: np.ndarray

    # The voxelgram itself
    voxelgram: np.ndarray              # uint8, shape (N, N, N)
    voxel_size: int                    # N
    box_size_A: float
    voxel_pitch_A: float               # box_size_A / N
    phi_actual: float                  # realised volume fraction of voxelgram

    rng_seed_used: int

    # MC uncertainties (param_name -> std). Empty if not run.
    params_std: dict = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Pure-math free functions (no class state)
# ---------------------------------------------------------------------------

def subtract_background(
    q: np.ndarray, I: np.ndarray,
    power_law_B: float, power_law_P: float, background: float,
) -> np.ndarray:
    """Return I - (B * q**-P + flat). q must be > 0."""
    q_safe = np.maximum(q, 1e-30)
    return I - (power_law_B * q_safe ** -power_law_P + background)


def add_background(
    q: np.ndarray, I_struct: np.ndarray,
    power_law_B: float, power_law_P: float, background: float,
) -> np.ndarray:
    """Inverse of subtract_background."""
    q_safe = np.maximum(q, 1e-30)
    return I_struct + (power_law_B * q_safe ** -power_law_P + background)


def alfa_threshold(phi: float) -> float:
    """Threshold level for thresholded Gaussian field given target volume fraction.

    alfa = sqrt(2) * erfinv(1 - 2*phi)

    phi is clamped to [1e-3, 1-1e-3] to avoid +/-inf from erfinv.
    """
    phi_safe = float(np.clip(phi, 1e-3, 1.0 - 1e-3))
    return float(np.sqrt(2.0) * _erfinv(1.0 - 2.0 * phi_safe))


def debye_autocorr(
    q: np.ndarray, I_corr: np.ndarray,
    r_grid: Optional[np.ndarray] = None,
    n_r: int = 256,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute the Debye phase autocorrelation gamma(r) from corrected I(Q).

    gamma(r) = ( 1 / (2*pi**2) ) * integral_0^inf  I_corr(Q) * Q**2 * sinc(Q*r) dQ

    where sinc(x) = sin(x)/x.  The result is normalised so gamma(0) = 1.

    Parameters
    ----------
    q       : 1-D Q array (sorted ascending, > 0).
    I_corr  : Background-subtracted intensity at each q.
    r_grid  : Optional 1-D r grid [A]. Default: linear from 0 to ~pi/q_min,
              n_r points.
    n_r     : Default r-grid resolution if r_grid is None.

    Returns
    -------
    (r, gamma) — both length n_r.
    """
    q = np.asarray(q, dtype=float)
    I_corr = np.asarray(I_corr, dtype=float)
    order = np.argsort(q)
    q_s = q[order]
    I_s = I_corr[order]

    if r_grid is None:
        r_max = float(np.pi / max(q_s[0], 1e-12))
        r_grid = np.linspace(0.0, r_max, n_r)

    qsq_I = (q_s ** 2) * I_s

    gamma = np.empty_like(r_grid, dtype=float)
    for i, r in enumerate(r_grid):
        if r <= 0.0:
            integrand = qsq_I  # sinc(0) = 1
        else:
            qr = q_s * r
            integrand = qsq_I * np.sin(qr) / qr
        gamma[i] = _simpson(integrand, x=q_s)

    if abs(gamma[0]) > 0:
        gamma /= gamma[0]
    return r_grid, gamma


def spectral_function(
    r: np.ndarray, gamma: np.ndarray,
    k_grid: Optional[np.ndarray] = None,
    n_k: int = 256,
) -> tuple[np.ndarray, np.ndarray]:
    """Spectral density F(k) of gamma(r): sinc transform of gamma * r**2.

    F(k) = ( 1 / (2*pi**2) ) * integral_0^inf  gamma(r) * r**2 * sinc(k*r) dr

    Negative bins are clipped to zero (numerical noise).

    Returns
    -------
    (k, F) — both length n_k.
    """
    r = np.asarray(r, dtype=float)
    gamma = np.asarray(gamma, dtype=float)

    if k_grid is None:
        r_max = float(r[-1] if r[-1] > 0 else 1.0)
        k_max = float(np.pi / max(r[1] - r[0], 1e-12)) if len(r) > 1 else 10.0
        k_grid = np.linspace(0.0, k_max, n_k)

    g_rsq = gamma * r ** 2

    F = np.empty_like(k_grid, dtype=float)
    for i, k in enumerate(k_grid):
        if k <= 0.0:
            integrand = g_rsq  # sinc(0) = 1
        else:
            kr = k * r
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                sinc = np.where(kr > 1e-12, np.sin(kr) / np.where(kr > 1e-12, kr, 1.0), 1.0)
            integrand = g_rsq * sinc
        F[i] = _simpson(integrand, x=r)

    F = np.maximum(F, 0.0)
    return k_grid, F


def generate_voxelgram(
    spectral_k: np.ndarray, spectral_F: np.ndarray,
    voxel_size: int, box_size_A: float,
    alfa: float,
    rng_seed: Optional[int] = None,
) -> tuple[np.ndarray, float, int]:
    """Generate the binary uint8 voxelgram via the GRF FFT pipeline.

    Steps:
      1. Build 3D radial k-grid (magnitudes) from FFT freqs of an N**3 cube
         with spacing pitch = box_size_A / N.
      2. Resample sqrt(spectral_F) onto each |k| via 1-D linear interpolation.
      3. Multiply FFT of white-noise cube by the resampled spectrum.
      4. Inverse FFT, take real part, threshold at alfa * sigma.

    Returns
    -------
    voxelgram   uint8 cube, shape (N, N, N), values 0 or 1.
    sigma       Standard deviation of the underlying scalar field.
    seed_used   The integer seed actually used (so it can be saved for
                reproducibility).
    """
    N = int(voxel_size)
    pitch = float(box_size_A) / N

    # k-vectors in physical units (rad/A): k = 2*pi * fftfreq(N, d=pitch).
    k_axis = 2.0 * np.pi * np.fft.fftfreq(N, d=pitch)
    kx, ky, kz = np.meshgrid(k_axis, k_axis, k_axis, indexing='ij')
    k_mag = np.sqrt(kx * kx + ky * ky + kz * kz)
    del kx, ky, kz

    # Resample sqrt(F) onto |k|. For |k| outside the supplied table we use 0.
    sqrt_F = np.sqrt(np.maximum(spectral_F, 0.0))
    spectrum_3d = np.interp(
        k_mag.ravel(), spectral_k, sqrt_F, left=sqrt_F[0], right=0.0
    ).reshape(k_mag.shape)
    del k_mag

    # White-noise cube.
    if rng_seed is None:
        seed_used = int(np.random.SeedSequence().entropy & 0xFFFFFFFF)
    else:
        seed_used = int(rng_seed)
    rng = np.random.default_rng(seed_used)
    noise = rng.standard_normal(size=(N, N, N)).astype(np.float32)

    # Filter in k-space.
    noise_k = np.fft.fftn(noise)
    del noise
    field_k = noise_k * spectrum_3d
    del noise_k, spectrum_3d
    field = np.real(np.fft.ifftn(field_k))
    del field_k

    sigma = float(field.std())
    if sigma <= 0.0:
        sigma = 1.0  # degenerate (e.g. all-zero spectrum) — return all-zeros voxelgram

    voxelgram = (field > alfa * sigma).astype(np.uint8)
    return voxelgram, sigma, seed_used


def voxelgram_to_iq(
    voxelgram: np.ndarray, voxel_pitch_A: float,
    q_target: np.ndarray,
) -> np.ndarray:
    """Compute spherically-averaged structure factor of a binary voxelgram.

    Returns I_struct(Q) such that the full SAS intensity is
    contrast * I_struct(Q) where contrast is the (Delta rho)**2 of phase
    against background.

    Steps:
      1. FFT of (voxelgram - <voxelgram>) -> F_k.  Subtracting the mean kills
         the DC spike that would otherwise dominate the lowest q bin.
      2. |F_k|**2 -> 3D power spectrum.
      3. Spherically average into 1-D bins of |k|.
      4. Convert to per-volume intensity: divide by N**3.
      5. Resample onto user q_target by log-log linear interpolation.

    Returns I_struct of shape q_target.shape.
    """
    N = voxelgram.shape[0]
    pitch = float(voxel_pitch_A)

    centered = voxelgram.astype(np.float32) - float(voxelgram.mean())
    F_k = np.fft.fftn(centered)
    pwr = (F_k.real * F_k.real + F_k.imag * F_k.imag).astype(np.float64)
    del F_k, centered

    # Build radial bins in k-space (rad/A).
    k_axis = 2.0 * np.pi * np.fft.fftfreq(N, d=pitch)
    kx, ky, kz = np.meshgrid(k_axis, k_axis, k_axis, indexing='ij')
    k_mag = np.sqrt(kx * kx + ky * ky + kz * kz)
    del kx, ky, kz

    # Bin edges: linear from 0 to k_max, ~N/2 bins.
    k_max = float(k_axis.max())
    n_bins = max(N // 2, 32)
    edges = np.linspace(0.0, k_max, n_bins + 1)
    bin_idx = np.digitize(k_mag.ravel(), edges) - 1
    valid = (bin_idx >= 0) & (bin_idx < n_bins)
    bin_idx = bin_idx[valid]
    pwr_flat = pwr.ravel()[valid]
    k_flat = k_mag.ravel()[valid]
    del k_mag, pwr

    counts = np.bincount(bin_idx, minlength=n_bins).astype(np.float64)
    sums = np.bincount(bin_idx, weights=pwr_flat, minlength=n_bins)
    k_sums = np.bincount(bin_idx, weights=k_flat, minlength=n_bins)

    nonzero = counts > 0
    I_radial = np.zeros(n_bins, dtype=np.float64)
    k_radial = np.zeros(n_bins, dtype=np.float64)
    I_radial[nonzero] = sums[nonzero] / counts[nonzero] / (N ** 3)
    k_radial[nonzero] = k_sums[nonzero] / counts[nonzero]

    # Drop the DC bin (k=0).
    if k_radial[0] == 0.0:
        nonzero[0] = False

    k_keep = k_radial[nonzero]
    I_keep = I_radial[nonzero]

    # Multiply by box volume so I_struct(Q) has units of 1/cm when contrast
    # is in 10**20 cm**-4 and lengths are in A.  (See e.g. Levitz 2007 and
    # IR3T_AvgIntensityFromBox in the Igor source: factor V = (N*pitch)**3
    # is needed because the FFT defines the structure factor without
    # normalisation by sample volume.)
    V_box = (N * pitch) ** 3
    I_keep = I_keep * V_box

    # Log-log interpolation onto q_target.  Outside the table return small
    # but positive values to keep ratios well-behaved.
    if len(k_keep) < 2:
        return np.zeros_like(q_target)

    log_k = np.log(np.maximum(k_keep, 1e-30))
    log_I = np.log(np.maximum(I_keep, 1e-30))
    log_qt = np.log(np.maximum(q_target, 1e-30))
    log_It = np.interp(log_qt, log_k, log_I, left=log_I[0], right=log_I[-1])
    return np.exp(log_It)


def derive_contrast_from_invariant(
    q: np.ndarray, I_corr: np.ndarray, phi: float,
) -> float:
    """Compute (Delta rho)**2 from the Porod invariant.

    Q* = integral_0^inf I_corr(Q) * Q**2 dQ = 2 * pi**2 * phi * (1 - phi) * (Delta rho)**2

    so (Delta rho)**2 = Q* / ( 2 * pi**2 * phi * (1 - phi) )
    """
    phi = float(np.clip(phi, 1e-3, 1.0 - 1e-3))
    order = np.argsort(q)
    q_s = q[order]
    I_s = I_corr[order]
    Qstar = float(_simpson((q_s ** 2) * I_s, x=q_s))
    denom = 2.0 * np.pi ** 2 * phi * (1.0 - phi)
    if denom <= 0:
        return 0.0
    return Qstar / denom


# ---------------------------------------------------------------------------
# Engine
# ---------------------------------------------------------------------------

class SaxsMorphEngine:
    """SAXS Morph engine: compute_voxelgram + fit + MC uncertainty."""

    def __init__(self):
        # Cache of (k, F) keyed by hash of (q_corr_bytes, B, P, flat).
        # The data are constant during a fit; only voxel size changes per call.
        self._spectral_cache: dict = {}
        # Optional cancellation hook installed by GUI workers; raises an
        # exception when set to abort the current fit between iterations.
        self._cancel_check = None

    # ----- one-shot evaluation ---------------------------------------------

    def compute_voxelgram(
        self,
        config: SaxsMorphConfig,
        q: np.ndarray,
        I: np.ndarray,
        dI: Optional[np.ndarray] = None,
        voxel_size_override: Optional[int] = None,
    ) -> SaxsMorphResult:
        """One full evaluation: data -> voxelgram -> model I(Q) -> chi**2.

        The voxel size used is ``voxel_size_override`` if given, else
        ``config.voxel_size_fit`` (callers that want the high-resolution
        rendering after a fit should pass voxel_size_override = voxel_size_render).
        """
        q = np.asarray(q, dtype=float)
        I = np.asarray(I, dtype=float)
        if dI is None:
            dI = np.maximum(I * 0.05, 1e-30)
        else:
            dI = np.asarray(dI, dtype=float)

        # Crop to fit window
        q_min = config.q_min if config.q_min is not None else float(q.min())
        q_max = config.q_max if config.q_max is not None else float(q.max())
        mask = (q >= q_min) & (q <= q_max)
        if not np.any(mask):
            raise ValueError("No data points within [q_min, q_max].")
        q_fit = q[mask]
        I_fit = I[mask]
        dI_fit = np.maximum(dI[mask], 1e-30)

        # Subtract background
        I_corr = subtract_background(
            q_fit, I_fit,
            config.power_law_B, config.power_law_P, config.background,
        )

        # Optionally derive contrast from invariant
        contrast = config.contrast
        if config.link_phi_contrast:
            contrast = derive_contrast_from_invariant(
                q_fit, np.maximum(I_corr, 0.0), config.volume_fraction,
            )

        # Cached spectral function
        r_grid, gamma_r, spectral_k, spectral_F = self._spectral(q_fit, I_corr)

        # Voxelgram
        N = int(voxel_size_override) if voxel_size_override else config.voxel_size_fit
        if N < 8:
            raise ValueError(f"voxel_size must be >= 8 (got {N})")
        alfa = alfa_threshold(config.volume_fraction)
        voxelgram, _sigma, seed_used = generate_voxelgram(
            spectral_k, spectral_F,
            voxel_size=N, box_size_A=config.box_size_A,
            alfa=alfa, rng_seed=config.rng_seed,
        )
        pitch = config.box_size_A / N
        phi_actual = float(voxelgram.mean())

        # Model intensity
        I_struct = voxelgram_to_iq(voxelgram, pitch, q_fit)
        I_model = add_background(
            q_fit, contrast * I_struct,
            config.power_law_B, config.power_law_P, config.background,
        )

        # Chi**2
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            resid = (I_fit - I_model) / dI_fit
        resid = np.where(np.isfinite(resid), resid, 0.0)
        chi2 = float(np.sum(resid * resid))
        dof = max(len(q_fit) - 1, 1)

        # Reflect derived contrast back so callers/UI can read it
        cfg_out = deepcopy(config)
        cfg_out.contrast = float(contrast)
        cfg_out.q_min = float(q_min)
        cfg_out.q_max = float(q_max)

        return SaxsMorphResult(
            config=cfg_out,
            chi_squared=chi2,
            reduced_chi_squared=chi2 / dof,
            dof=dof,
            timestamp=datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            data_q=q_fit,
            data_I=I_fit,
            data_dI=dI_fit,
            data_I_corr=I_corr,
            model_q=q_fit,
            model_I=I_model,
            r_grid=r_grid,
            gamma_r=gamma_r,
            spectral_k=spectral_k,
            spectral_F=spectral_F,
            voxelgram=voxelgram,
            voxel_size=N,
            box_size_A=float(config.box_size_A),
            voxel_pitch_A=float(pitch),
            phi_actual=phi_actual,
            rng_seed_used=seed_used,
        )

    # ----- fitting ---------------------------------------------------------

    def fit(
        self,
        config: SaxsMorphConfig,
        q: np.ndarray,
        I: np.ndarray,
        dI: Optional[np.ndarray] = None,
    ) -> SaxsMorphResult:
        """Refine fittable parameters by least_squares (or Nelder-Mead).

        Parameters that may be fit (via fit_* flags):
          volume_fraction, contrast, power_law_B, power_law_P, background.

        The fit loop is hard-clamped to ``MAX_FIT_VOXEL_SIZE`` to bound memory.
        After convergence, one final ``compute_voxelgram`` is run at
        ``config.voxel_size_render`` and that high-resolution voxelgram is
        returned in the result.
        """
        cfg = deepcopy(config)

        # Clamp fit-time voxel size
        cfg.voxel_size_fit = int(min(cfg.voxel_size_fit, MAX_FIT_VOXEL_SIZE))

        x0, lo, hi, keys = self._pack_params(cfg)

        if len(x0) == 0:
            # Nothing to fit — just evaluate at render resolution
            return self.compute_voxelgram(
                cfg, q, I, dI, voxel_size_override=cfg.voxel_size_render,
            )

        x0_arr = np.array(x0, dtype=float)

        def _residuals(x):
            if self._cancel_check is not None:
                self._cancel_check()
            self._unpack_params(x, keys, cfg)
            res = self.compute_voxelgram(
                cfg, q, I, dI, voxel_size_override=cfg.voxel_size_fit,
            )
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                r = (res.data_I - res.model_I) / res.data_dI
            return np.where(np.isfinite(r), r, 0.0)

        def _chi2(x):
            r = _residuals(x)
            return float(np.sum(r * r))

        # Stash hooks so a GUI worker can monkey-patch _residuals/_chi2 on
        # this engine instance for cancellation (mirror ModelingEngine pattern).
        self._residuals = _residuals
        self._chi2 = _chi2

        if cfg.no_limits:
            res_opt = minimize(
                self._chi2, x0_arr,
                method='Nelder-Mead',
                options={'maxiter': 200, 'xatol': 1e-3, 'fatol': 1e-3},
            )
            x_best = res_opt.x
        else:
            lo_arr = np.array(lo, dtype=float)
            hi_arr = np.array(hi, dtype=float)
            # least_squares wants strict inequalities lo < x0 < hi
            x0_clipped = np.clip(x0_arr, lo_arr + 1e-12, hi_arr - 1e-12)
            res_opt = least_squares(
                self._residuals, x0_clipped,
                bounds=(lo_arr, hi_arr),
                method='trf',
                max_nfev=100,
                xtol=1e-4, ftol=1e-4,
            )
            x_best = res_opt.x

        self._unpack_params(x_best, keys, cfg)

        # Final high-res evaluation
        return self.compute_voxelgram(
            cfg, q, I, dI, voxel_size_override=cfg.voxel_size_render,
        )

    # ----- Monte-Carlo uncertainty -----------------------------------------

    def calculate_uncertainty_mc(
        self,
        config: SaxsMorphConfig,
        q: np.ndarray, I: np.ndarray, dI: np.ndarray,
        n_runs: int = 10,
        progress_cb=None,
    ) -> dict:
        """Estimate parameter standard deviations by Monte Carlo perturbation."""
        rng = np.random.default_rng(config.rng_seed)
        collected: dict = {name: [] for name in
                           ('volume_fraction', 'contrast',
                            'power_law_B', 'power_law_P', 'background')}

        for k in range(n_runs):
            if progress_cb:
                progress_cb(k, n_runs)
            I_perturbed = I + rng.standard_normal(I.shape) * np.maximum(dI, 1e-30)
            cfg = deepcopy(config)
            try:
                res = self.fit(cfg, q, I_perturbed, dI)
            except Exception:
                continue
            for name in collected:
                collected[name].append(getattr(res.config, name))

        if progress_cb:
            progress_cb(n_runs, n_runs)

        out = {}
        for name, vals in collected.items():
            if len(vals) > 1:
                out[name] = float(np.std(vals, ddof=1))
        return out

    # ----- internals -------------------------------------------------------

    def _spectral(self, q_fit, I_corr):
        """Compute and cache (r, gamma, k, F)."""
        key = (q_fit.tobytes(), I_corr.tobytes())
        cached = self._spectral_cache.get(key)
        if cached is not None:
            return cached
        r, gamma = debye_autocorr(q_fit, I_corr)
        k, F = spectral_function(r, gamma)
        # Bound cache size
        if len(self._spectral_cache) > 8:
            self._spectral_cache.clear()
        self._spectral_cache[key] = (r, gamma, k, F)
        return r, gamma, k, F

    def _pack_params(self, cfg: SaxsMorphConfig):
        """Flatten fit_* booleans into x0 / lo / hi / keys lists."""
        names = ['volume_fraction', 'contrast',
                 'power_law_B', 'power_law_P', 'background']
        # Skip contrast if linked to phi via invariant (it's derived, not fit).
        x0, lo, hi, keys = [], [], [], []
        for name in names:
            if name == 'contrast' and cfg.link_phi_contrast:
                continue
            if not getattr(cfg, f'fit_{name}', False):
                continue
            val = float(getattr(cfg, name))
            limits = getattr(cfg, f'{name}_limits')
            x0.append(val)
            lo.append(float(limits[0]))
            hi.append(float(limits[1]))
            keys.append(name)
        return x0, lo, hi, keys

    def _unpack_params(self, x, keys, cfg: SaxsMorphConfig):
        for val, name in zip(x, keys):
            setattr(cfg, name, float(val))
