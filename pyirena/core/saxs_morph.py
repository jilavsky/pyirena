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

    Workflow
    --------
    The SAXS Morph method does NOT iteratively fit the model parameters
    (phi, contrast).  Given the data invariant and one of (phi, Delta rho^2),
    the other is uniquely determined; given those plus a spectral function
    derived from the data, the voxelgram is generated *deterministically*
    (modulo the RNG seed).  The only true fits are the two background
    pre-fits: a Power-law fit at low Q and a flat-background fit at high Q.

    Q ranges
    --------
    q_min, q_max                        Range used for the GRF modelling
                                         step (read from the main I(Q)
                                         cursors at "Calculate 3D" time).
    power_law_q_min, power_law_q_max    Range for the Power-law pre-fit
                                         (typically a low-Q window).
    background_q_min, background_q_max  Range for the flat-background
                                         pre-fit (typically a high-Q window).

    Input mode
    ----------
    input_mode  One of:
        'phi'      User supplies phi; contrast is derived from the invariant.
        'contrast' User supplies contrast; phi is derived from the invariant.
        'both'     User supplies both; no derivation (invariant ignored).

    The fit_* flags and *_limits fields are kept for backward compatibility
    with the (deprecated) Engine.fit() method but are no longer used by the
    GUI workflow.
    """
    # Modelling Q range (driven by main cursors)
    q_min: Optional[float] = None
    q_max: Optional[float] = None

    # Background pre-fit Q ranges
    power_law_q_min: Optional[float] = None
    power_law_q_max: Optional[float] = None
    background_q_min: Optional[float] = None
    background_q_max: Optional[float] = None

    # Voxel grid
    voxel_size_fit: int = 128
    voxel_size_render: int = 256
    box_size_A: float = 1000.0

    # Two-phase parameters
    input_mode: str = 'phi'           # 'phi' | 'contrast' | 'both'
    volume_fraction: float = 0.3
    contrast: float = 1.0

    # Background (Power-law + Flat) — values populated by pre-fits
    power_law_B: float = 0.0
    power_law_P: float = 4.0
    background: float = 0.0

    # Post-threshold smoothing of the binary voxelgram (in voxels).
    # Mimics Igor Pro's ``ImageFilter/N=5 gauss3d`` step that knocks
    # down the per-voxel numerical noise in the GRF realisation.
    # 0 = no smoothing (output stays uint8 binary).
    # > 0 = scipy.ndimage.gaussian_filter with this sigma (output is
    #       float32 in [0, 1]).  Default 1.0 matches Igor's default.
    smooth_sigma: float = 1.0

    # ----- Backward-compatibility (deprecated, used only by Engine.fit) ----
    fit_volume_fraction: bool = False
    volume_fraction_limits: tuple = (0.05, 0.95)
    fit_contrast: bool = False
    contrast_limits: tuple = (0.0, 1e10)
    link_phi_contrast: bool = True   # legacy alias: True == input_mode 'phi'
    fit_power_law_B: bool = False
    power_law_B_limits: tuple = (0.0, 1e10)
    fit_power_law_P: bool = False
    power_law_P_limits: tuple = (0.0, 6.0)
    fit_background: bool = False
    background_limits: tuple = (0.0, 1e10)
    no_limits: bool = False
    n_mc_runs: int = 10
    # -----------------------------------------------------------------------

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

    # Radius of gyration derived from the autocorrelation function gamma(r).
    # nan if it cannot be computed (e.g. degenerate autocorrelation).
    rg_A: float = float('nan')

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


def berk_lut(alfa: float, n: int = 2000) -> tuple[np.ndarray, np.ndarray]:
    """Build a forward lookup table for the Berk thresholding transformation.

    For a Gaussian field with autocorrelation g(r), thresholded at
    level alfa to give binary indicator eta(r) with mean phi, the
    indicator's autocorrelation is

        gamma_eta(r) = T(g(r), alfa)
        T(g, alfa) = (1/(2 pi)) * integral_0^g exp(-alfa^2/(1+t)) / sqrt(1-t^2) dt

    See Berk (1991), Phys. Rev. E 51, 4141.

    The integrand is integrably singular at g = +/- 1, where 1/sqrt(1-t^2)
    diverges.  We substitute t = sin(theta) to remove the singularity:
    sqrt(1-t^2) = cos(theta), dt = cos(theta) d(theta), so

        T(g, alfa) = (1/(2 pi)) * integral_0^arcsin(g) exp(-alfa^2/(1+sin theta)) d theta

    The integrand on theta is bounded in (0, 1) and trapezoidal
    integration converges quickly.

    Returns
    -------
    (g_grid, T_grid) : two 1-D arrays of length ``n`` covering
        g in (-1, 1).  T_grid is monotonically increasing in g (the
        integrand is positive), so use ``np.interp(target_T, T_grid, g_grid)``
        to invert T.
    """
    from scipy.integrate import cumulative_trapezoid
    # theta ranges over (-pi/2, pi/2) excluding the singular endpoints
    theta = np.linspace(-np.pi / 2 * 0.999, np.pi / 2 * 0.999, n)
    g_grid = np.sin(theta)
    integrand = np.exp(-alfa ** 2 / (1.0 + g_grid))   # smooth, in (0, 1]
    cum = cumulative_trapezoid(integrand, theta, initial=0.0)
    # cum[i] = integral from theta[0] to theta[i].  We want integral from 0,
    # so subtract the value at theta=0.
    cum_at_zero = float(np.interp(0.0, theta, cum))
    T_grid = (cum - cum_at_zero) / (2.0 * np.pi)
    return g_grid, T_grid


def berk_invert(target_T, alfa: float, lut: tuple = None) -> np.ndarray:
    """Solve T(g, alfa) = target_T for g, vectorised over target_T.

    Uses the forward LUT from :func:`berk_lut` and ``np.interp`` to
    invert.  Values of ``target_T`` outside the LUT range are clamped
    to ``g = +/- (1 - epsilon)``.

    Parameters
    ----------
    target_T : scalar or array of desired indicator autocorrelation values.
               Should be in [-T_max, T_max] where T_max = T(1, alfa) ~= phi(1-phi).
    alfa     : threshold parameter (sqrt(2) * erfinv(1 - 2 phi)).
    lut      : optional pre-built (g_grid, T_grid) tuple from berk_lut.
               Saves the LUT-build cost when inverting many points.

    Returns
    -------
    g : array of the same shape as target_T.
    """
    if lut is None:
        g_grid, T_grid = berk_lut(alfa)
    else:
        g_grid, T_grid = lut
    return np.interp(np.asarray(target_T, dtype=float), T_grid, g_grid)


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


def compute_rg_from_autocorr(r_grid: np.ndarray, gamma_r: np.ndarray) -> float:
    """Estimate the radius of gyration from the Debye autocorrelation function.

    Uses the standard PDDF formula:
        Rg² = ∫r⁴ γ(r) dr / (2 ∫r² γ(r) dr)

    where p(r) ∝ r² γ(r) is the pair distance distribution function.
    Integration is truncated at the first downward zero-crossing of γ(r) to
    avoid negative contributions from interference oscillations at large r.

    Returns nan if the denominator is non-positive or if there are too few
    r-grid points after truncation.
    """
    r = np.asarray(r_grid, dtype=float)
    g = np.asarray(gamma_r, dtype=float)
    # Truncate at first downward zero crossing (γ goes from 1 toward 0)
    crossings = np.where((g[:-1] > 0) & (g[1:] <= 0))[0]
    n = int(crossings[0]) + 1 if len(crossings) > 0 else len(r)
    r = r[:n]
    g = np.maximum(g[:n], 0.0)
    if len(r) < 3:
        return float('nan')
    num = _simpson(r ** 4 * g, x=r)
    den = _simpson(r ** 2 * g, x=r)
    if den <= 0.0:
        return float('nan')
    return float(np.sqrt(0.5 * num / den))


def generate_voxelgram(
    spectral_k: np.ndarray, spectral_F: np.ndarray,
    voxel_size: int, box_size_A: float,
    alfa: float,
    rng_seed: Optional[int] = None,
    smooth_sigma: float = 0.0,
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

    if smooth_sigma > 0:
        # Mimic Igor's ``ImageFilter/N=5 gauss3d`` post-processing step.
        # Output is float32 in [0, 1] (no longer strictly binary, but
        # the 3D viewer still renders an isosurface at 0.5 cleanly and
        # the spectrum has its high-Q noise damped by the Gaussian's
        # frequency response).
        from scipy.ndimage import gaussian_filter
        voxelgram = gaussian_filter(
            voxelgram.astype(np.float32), sigma=float(smooth_sigma),
        )

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

    # Convert from <|F_DFT|^2>/N**3 (dimensionless) to a structure-factor
    # pre-factor with the right units so that I [cm**-1] = contrast
    # [10**20 cm**-4] * I_struct.
    #
    # Continuous Fourier transform vs DFT:
    #   F_continuous(q) ~= pitch**3 * F_DFT(k)
    #   |F_continuous|**2 ~= pitch**6 * |F_DFT|**2
    #   S(q) [length**3] = |F_continuous|**2 / V_sample
    #                     = pitch**3 / N**3 * |F_DFT|**2
    #
    # In cm units: S(q) [cm**3] = pitch_A**3 * 10**(-24) / N**3 * <|F_DFT|^2>
    #   I [cm**-1] = (Delta rho)**2 [cm**-4] * S [cm**3]
    #              = (contrast_e20cm4 * 10**20) * pitch_A**3 * 10**(-24)/N**3 * <|F|^2>
    #              = contrast_e20cm4 * pitch_A**3 * 10**(-4)/N**3 * <|F|^2>
    #              = contrast_e20cm4 * pitch_A**3 * 10**(-4) * (<|F|^2>/N**3)
    #
    # I_keep currently equals <|F|^2>/N**3, so multiply by pitch_A**3 * 1e-4
    # (NOT by the (N*pitch)**3 box volume — that double-counted N**3 and
    #  was missing the 10**(-24) Angstrom-to-cm conversion).
    I_keep = I_keep * (pitch ** 3) * 1e-4

    # Log-log interpolation onto q_target.  Outside the table return small
    # but positive values to keep ratios well-behaved.
    if len(k_keep) < 2:
        return np.zeros_like(q_target)

    log_k = np.log(np.maximum(k_keep, 1e-30))
    log_I = np.log(np.maximum(I_keep, 1e-30))
    log_qt = np.log(np.maximum(q_target, 1e-30))
    log_It = np.interp(log_qt, log_k, log_I, left=log_I[0], right=log_I[-1])
    return np.exp(log_It)


def compute_invariant_extrapolated(
    q: np.ndarray, I_corr: np.ndarray, n_avg: int = 5,
) -> float:
    """Porod invariant Q* = integral_0^inf I(Q) * Q**2 dQ with extrapolations
    beyond the data range.

    The numerical integral covers ``q[0]..q[-1]``.  Beyond the data ends we
    add closed-form contributions matching the standard SAS extrapolations
    used by the Igor reference (``IR3T_ExtendDataCalcParams``):

      Low-Q  (Q < q_min)
        Approximate I(Q) ~ I0 (constant) where
        I0 = mean of the first ``n_avg`` data points.
        Contribution: I0 * q_min**3 / 3.

      High-Q (Q > q_max)
        Approximate I(Q) ~ K * Q**-4 (Porod tail) where
        K = mean of (I * Q**4) over the last ``n_avg`` data points.
        Contribution: integral_qmax^inf K * Q**-2 dQ = K / q_max.

    Returns Q*_numerical in units of [Å**-3 * cm**-1] when ``q`` is in
    Å**-1 and ``I_corr`` is in cm**-1.

    Without extrapolation (e.g. n_avg=0) the function reduces to the
    bare Simpson integral over the data range.
    """
    order = np.argsort(q)
    q_s = q[order]
    I_s = I_corr[order]

    # Numerical part
    Q_data = float(_simpson((q_s ** 2) * I_s, x=q_s))

    if n_avg <= 0 or len(q_s) < 2:
        return Q_data

    n = min(int(n_avg), len(q_s))

    # Low-Q extrapolation: constant I = mean of first n points
    q_min = float(q_s[0])
    I0 = float(np.mean(I_s[:n]))
    Q_low = I0 * q_min ** 3 / 3.0 if q_min > 0 else 0.0

    # High-Q extrapolation: Porod tail I = K * Q^-4
    q_max = float(q_s[-1])
    K = float(np.mean(I_s[-n:] * q_s[-n:] ** 4))
    Q_high = K / q_max if q_max > 0 else 0.0

    return Q_data + Q_low + Q_high


def derive_contrast_from_invariant(
    q: np.ndarray, I_corr: np.ndarray, phi: float,
) -> float:
    """Compute (Delta rho)**2 from the Porod invariant.

    Q* = integral_0^inf I(Q) * Q**2 dQ = 2 * pi**2 * phi * (1 - phi) * (Delta rho)**2

    Unit conversion
    ---------------
    With I in cm**-1 and Q in A**-1 (the pyirena standard), the *numerical*
    integral has units of A**-3 * cm**-1.  But the right-hand side of the
    invariant identity is in cm**-4 (or 10**20 cm**-4 in the SAS-conventional
    "1e20" reporting unit).  Since 1 A**-1 = 10**8 cm**-1, we have

        Q*_numerical [A**-3 * cm**-1] * 10**24 = Q*_proper [cm**-4]

    and dividing by 10**20 to express the result in 10**20 cm**-4 leaves a
    net factor of 10**4 multiplying the numerical integral.  This is the
    (10**4) constant in the formula below; it is the difference between
    the value Igor reports and the unit-naive integration we used in
    pyirena <= 0.5.1.

    Returns
    -------
    Delta rho**2 in 10**20 cm**-4 (the standard SAS reporting unit, also
    used by the Igor reference implementation).
    """
    phi = float(np.clip(phi, 1e-3, 1.0 - 1e-3))
    Qstar_num = compute_invariant_extrapolated(q, I_corr)   # A**-3 * cm**-1
    Qstar_e20cm4 = Qstar_num * 1e4                          # 10**20 cm**-4
    denom = 2.0 * np.pi ** 2 * phi * (1.0 - phi)
    if denom <= 0:
        return 0.0
    return Qstar_e20cm4 / denom


def derive_phi_from_invariant(
    q: np.ndarray, I_corr: np.ndarray, contrast: float,
) -> float:
    """Inverse of :func:`derive_contrast_from_invariant`.

    Given Q* (computed from data) and contrast Delta rho^2, solve

        Q* = 2 * pi^2 * phi * (1 - phi) * Delta rho^2

    for phi.  This is a quadratic in phi:

        phi^2 - phi + K = 0,  where K = Q* / (2 * pi^2 * Delta rho^2)

    with roots phi = 0.5 * (1 +/- sqrt(1 - 4K)).  We return the *smaller*
    root (phi < 0.5) which corresponds to the minority phase.

    If 1 - 4K < 0 the data is inconsistent with the supplied contrast (the
    invariant is too large to be explained by any two-phase split at this
    contrast) — we return 0.5 in that case as a sentinel.
    """
    if contrast <= 0:
        return 0.5
    Qstar_num = compute_invariant_extrapolated(q, I_corr)   # A**-3 * cm**-1
    Qstar_e20cm4 = Qstar_num * 1e4                          # 10**20 cm**-4
    # contrast is in 10**20 cm**-4 to match the inverse function above
    K = Qstar_e20cm4 / (2.0 * np.pi ** 2 * contrast)
    disc = 1.0 - 4.0 * K
    if disc < 0:
        return 0.5
    return float(0.5 * (1.0 - np.sqrt(disc)))


def fit_power_law_bg(
    q: np.ndarray, I: np.ndarray,
    q_min: float, q_max: float,
) -> tuple[float, float]:
    """Linear least-squares fit of log10(I) = log10(B) - P * log10(Q) over [q_min, q_max].

    Returns (B, P) where I_pl(Q) = B * Q**(-P).  Suitable for a low-Q
    window where the data is dominated by a power-law tail (e.g. Porod's
    Q^-4 in a sharp-interface system).

    Falls back to (0.0, 4.0) if the window contains fewer than 2 valid
    (positive Q, positive I) points.
    """
    q = np.asarray(q, dtype=float)
    I = np.asarray(I, dtype=float)
    mask = (q >= q_min) & (q <= q_max) & (q > 0) & (I > 0)
    if mask.sum() < 2:
        return 0.0, 4.0
    lq = np.log10(q[mask])
    lI = np.log10(I[mask])
    # lI = log10(B) - P * lq    →    polyfit returns [slope, intercept]
    slope, intercept = np.polyfit(lq, lI, 1)
    P = float(-slope)
    B = float(10.0 ** intercept)
    return B, P


def fit_flat_bg(
    q: np.ndarray, I: np.ndarray,
    q_min: float, q_max: float,
    power_law_B: float = 0.0, power_law_P: float = 4.0,
) -> float:
    """Estimate flat background as the median of (I - power_law) over [q_min, q_max].

    Use this in a high-Q window where the structural scattering has
    decayed and the residual is dominated by detector / electronic noise.

    Median is preferred over mean because high-Q data points often have a
    long-tailed noise distribution.
    """
    q = np.asarray(q, dtype=float)
    I = np.asarray(I, dtype=float)
    mask = (q >= q_min) & (q <= q_max)
    if not np.any(mask):
        return 0.0
    q_w = q[mask]
    I_w = I[mask]
    pl = power_law_B * np.maximum(q_w, 1e-30) ** -power_law_P
    return float(np.median(I_w - pl))


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

        # Resolve input mode: legacy `link_phi_contrast=True` maps to 'phi'.
        mode = (config.input_mode or '').lower()
        if mode not in ('phi', 'contrast', 'both'):
            mode = 'phi' if config.link_phi_contrast else 'both'

        phi = float(config.volume_fraction)
        contrast = float(config.contrast)
        I_corr_pos = np.maximum(I_corr, 0.0)

        if mode == 'phi':
            # User input: phi.  Derive contrast from invariant.
            contrast = derive_contrast_from_invariant(q_fit, I_corr_pos, phi)
        elif mode == 'contrast':
            # User input: contrast.  Derive phi from invariant.
            phi = derive_phi_from_invariant(q_fit, I_corr_pos, contrast)
        # else 'both': use both as supplied.

        # gamma_r_norm is the indicator (binary phase) autocorrelation
        # normalised to gamma(0) = 1.  Its absolute scale should be
        # phi*(1-phi), so we rescale before the Berk inversion below.
        r_grid, gamma_r_norm, _spectral_k_unused, _spectral_F_unused = (
            self._spectral(q_fit, I_corr)
        )

        # Voxelgram size
        N = int(voxel_size_override) if voxel_size_override else config.voxel_size_fit
        if N < 8:
            raise ValueError(f"voxel_size must be >= 8 (got {N})")
        alfa = alfa_threshold(phi)

        # ── Berk inversion ──────────────────────────────────────────────
        # The Gaussian field generated by FFT-filtering white noise with
        # sqrt(F(k)) has autocorrelation = inverse-FT of F(k).  After
        # thresholding at alfa, the BINARY indicator's autocorrelation is
        # T(g(r), alfa), NOT g(r) itself.  To match the data's measured
        # indicator autocorrelation gamma(r), we must:
        #   1. Rescale gamma_norm to physical: gamma_phys = gamma_norm * phi(1-phi)
        #   2. Invert the Berk relation: find g(r) such that T(g(r), alfa) = gamma_phys(r)
        #   3. Use FT of g(r) -- not gamma -- as the spectral function for the field.
        # Without this step the model I(Q) is systematically too low
        # (by a factor T'(0,alfa) at small gamma, which can be 5-10x).
        gamma_phys = gamma_r_norm * (phi * (1.0 - phi))
        lut = berk_lut(alfa)
        g_r = berk_invert(gamma_phys, alfa, lut=lut)

        # Spectral function of the underlying GAUSSIAN field (FT of g(r)).
        # Compute on a k-grid that covers the 3D FFT range used below
        # so we don't lose high-k content during interpolation.
        pitch = config.box_size_A / N
        k_max_3d = float(np.sqrt(3.0) * np.pi / pitch) * 1.05  # +5% safety
        spectral_k, spectral_F = spectral_function(
            r_grid, g_r,
            k_grid=np.linspace(0.0, k_max_3d, max(512, N)),
        )
        # ALWAYS generate binary voxelgram for I(Q) computation.
        # Smoothing is applied SEPARATELY for display only (see below).
        # Using a smoothed field for I(Q) would reduce the variance
        # (smoothing damps high-frequency content), causing the model
        # intensity to be systematically lower than the data.
        voxelgram, _sigma, seed_used = generate_voxelgram(
            spectral_k, spectral_F,
            voxel_size=N, box_size_A=config.box_size_A,
            alfa=alfa, rng_seed=config.rng_seed,
            smooth_sigma=0.0,   # binary for physics — display smoothing happens in GUI
        )
        # pitch already computed above for the spectral k-grid sizing
        phi_actual = float(voxelgram.mean())
        rg_A = compute_rg_from_autocorr(r_grid, gamma_r_norm)

        # Model intensity (from binary voxelgram — correct physics)
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

        # Reflect resolved phi/contrast back so callers/UI can read them
        cfg_out = deepcopy(config)
        cfg_out.volume_fraction = float(phi)
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
            gamma_r=gamma_r_norm,
            spectral_k=spectral_k,
            spectral_F=spectral_F,
            voxelgram=voxelgram,
            voxel_size=N,
            box_size_A=float(config.box_size_A),
            voxel_pitch_A=float(pitch),
            phi_actual=phi_actual,
            rng_seed_used=seed_used,
            rg_A=rg_A,
        )

    # ----- fitting ---------------------------------------------------------

    def fit(
        self,
        config: SaxsMorphConfig,
        q: np.ndarray,
        I: np.ndarray,
        dI: Optional[np.ndarray] = None,
    ) -> SaxsMorphResult:
        """**DEPRECATED**: kept only for backward compatibility / regression tests.

        The standard SAXS Morph workflow does NOT iteratively fit the model
        parameters; volume fraction and contrast are linked through the data
        invariant and the voxelgram is computed deterministically from them
        plus the spectral function.  Use
        :meth:`compute_voxelgram` instead, after pre-fitting the background
        with :func:`fit_power_law_bg` and :func:`fit_flat_bg`.

        This method still varies any parameters whose ``fit_*`` flag is True
        (volume_fraction, contrast, power_law_B, power_law_P, background)
        via least_squares/Nelder-Mead and was the basis of the original
        implementation; it is preserved only to avoid breaking existing
        scripts and tests.
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
