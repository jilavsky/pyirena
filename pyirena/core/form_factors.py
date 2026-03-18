"""
Form factors for SAS size distribution analysis.

All functions return the form factor contribution per unit volume fraction,
in units suitable for direct use in the G-matrix:

    sphere_ff(q, r) → V(r) × F_norm²(Q, r)   [Å³]

The G-matrix is then:

    G[i, j] = sphere_ff(q[i], r[j]) × contrast × 1e-4   [cm^-1]

where contrast = (Δρ)² in units of 10^20 cm^-4, and the factor 1e-4 converts
the volume from Å³ to cm³ (1 Å³ = 1e-24 cm³; 1e-24 × 1e20 = 1e-4).

The G-matrix columns represent the scattered intensity per unit volume fraction
from monodisperse particles of radius r[j]:

    I(Q) [cm^-1] = Σ_j G[i,j] × φ_j

where φ_j is the volume fraction of particles in bin j (dimensionless).
"""

from __future__ import annotations

import numpy as np
from numpy.polynomial.legendre import leggauss


# ──────────────────────────────────────────────────────────────────────────────
# Sphere form factor
# ──────────────────────────────────────────────────────────────────────────────

def _sphere_amplitude(qr: np.ndarray) -> np.ndarray:
    """
    Normalised sphere amplitude F(Qr).

    F(Qr) = 3 [sin(Qr) - Qr·cos(Qr)] / (Qr)³

    Small-Qr branch uses a Taylor expansion to avoid catastrophic cancellation.
    Large-Qr branch uses pure NumPy trig (sin/cos), which is AVX-vectorised and
    much faster than scipy.special.spherical_jn for large arrays.
    The identity used is:  3·j₁(x)/x = 3·(sin x − x·cos x) / x³.
    """
    out = np.empty_like(qr, dtype=float)
    small = qr < 1e-3
    # Taylor: F ≈ 1 - (Qr)²/10 + (Qr)⁴/280 + …
    qr2 = qr[small] ** 2
    out[small] = 1.0 - qr2 / 10.0 + qr2 ** 2 / 280.0
    x = qr[~small]
    # 3·(sin x − x·cos x) / x³  —  algebraically identical to 3·j₁(x)/x
    out[~small] = 3.0 * (np.sin(x) - x * np.cos(x)) / (x ** 3)
    return out


def sphere_ff(q: np.ndarray, r: float) -> np.ndarray:
    """
    Sphere form factor per unit volume fraction.

    Returns V(r) × F_norm²(Q, r)   [Å³]

    where V(r) = (4/3)π r³  and  F_norm = 3[sin(Qr) - Qr·cos(Qr)]/(Qr)³.

    At Q=0:  sphere_ff → V(r).
    Multiply by (Δρ)² × 1e-4 to get the scattering intensity per unit
    volume fraction in cm^-1.

    Args:
        q:  1-D array of Q values [Å^-1]
        r:  Sphere radius [Å]

    Returns:
        1-D array of length len(q)  [Å³]
    """
    q = np.asarray(q, dtype=float)
    V = (4.0 / 3.0) * np.pi * r ** 3
    F = _sphere_amplitude(q * r)
    return V * F ** 2


# ──────────────────────────────────────────────────────────────────────────────
# Spheroid / Ellipsoid form factor
# ──────────────────────────────────────────────────────────────────────────────

# Gauss-Legendre nodes and weights cached at module level.
# 50 points match the quadrature order used in Igor IR1R_CalcSpheroidFormFactor.
_GL_NODES, _GL_WEIGHTS = leggauss(50)


def spheroid_ff(q: np.ndarray, r: float, aspect_ratio: float = 1.0) -> np.ndarray:
    """
    Spheroid (oblate or prolate) form factor per unit volume fraction.

    Returns V_spheroid × ∫₀¹ F_norm²(Q, r_eff(θ)) d(cosθ)   [Å³]

    The spheroid semi-axes are (r, r, r·AR).  The amplitude at polar angle θ
    relative to Q is V_spheroid × F_norm(Q, r_eff(θ)) where:

        r_eff(θ) = r · √[1 + (AR²-1)·cos²(θ)]

    Integration over cos(θ) ∈ [0, 1] uses 50-point Gauss-Legendre quadrature,
    matching Igor IR1R_CalcSpheroidFormFactor.

    At AR=1 the result equals sphere_ff(q, r).

    Args:
        q:            1-D array of Q values [Å^-1]
        r:            Equatorial radius [Å]
        aspect_ratio: AR = polar/equatorial axis ratio (default 1 → sphere)

    Returns:
        1-D array of length len(q)  [Å³]
    """
    q = np.asarray(q, dtype=float)
    AR = float(aspect_ratio)

    V = (4.0 / 3.0) * np.pi * r ** 3 * AR

    # Map Gauss-Legendre [-1, 1] nodes to cos(θ) ∈ [0, 1]
    cos_t = 0.5 * (_GL_NODES + 1.0)   # shape (50,)
    weights = 0.5 * _GL_WEIGHTS        # Jacobian factor ½ included

    # r_eff for each (q, cos_theta) pair → shape (len(q), 50)
    r_eff = r * np.sqrt(1.0 + (AR ** 2 - 1.0) * cos_t[np.newaxis, :] ** 2)
    qr_eff = q[:, np.newaxis] * r_eff  # (len(q), 50)

    F_eff = _sphere_amplitude(qr_eff)  # (len(q), 50)  — F_norm for each θ

    # Integrand: V × F_norm²(Q, r_eff(θ))  [Å³]
    integrand = F_eff ** 2 * V         # (len(q), 50)
    result = integrand @ weights        # (len(q),)

    return result


# ──────────────────────────────────────────────────────────────────────────────
# G-matrix builder — fully vectorized (no Python loop over r bins)
# ──────────────────────────────────────────────────────────────────────────────

def _build_g_sphere(
    q: np.ndarray,
    r_grid: np.ndarray,
    contrast: float,
) -> np.ndarray:
    """Vectorized G matrix for sphere form factor.

    Computes the entire (M × N) matrix with a single NumPy broadcast, avoiding
    the Python loop over radius bins that the old per-column approach required.

    qr[i,j] = q[i] * r[j]   → _sphere_amplitude works on any-shaped array.
    G[i,j]  = V(r[j]) * F(qr[i,j])² * contrast * 1e-4
    """
    qr = q[:, np.newaxis] * r_grid[np.newaxis, :]   # (M, N)
    V  = (4.0 / 3.0) * np.pi * r_grid ** 3          # (N,)
    F  = _sphere_amplitude(qr)                        # (M, N)
    return V[np.newaxis, :] * F ** 2 * (contrast * 1e-4)   # (M, N)


def _build_g_spheroid(
    q: np.ndarray,
    r_grid: np.ndarray,
    contrast: float,
    aspect_ratio: float = 1.0,
) -> np.ndarray:
    """G matrix for spheroid form factor, per-r-bin loop over Gauss-Legendre.

    A fully-fused (M × N × K) tensor would be memory-bandwidth-bound at
    ~16 MB.  Instead we loop over N radius bins so each (M × K) = ~80 KB
    slice fits in L2 cache.  The per-bin NumPy overhead is small because
    _sphere_amplitude uses AVX-vectorised sin/cos (no scipy overhead).

    G[i,j] = V(r[j]) * ∫₀¹ F²(Q_i, r_eff(r_j, θ)) d(cosθ)
    """
    AR = float(aspect_ratio)
    cos_t   = 0.5 * (_GL_NODES + 1.0)   # (K,)
    weights = 0.5 * _GL_WEIGHTS          # (K,)

    M, N = len(q), len(r_grid)
    G = np.empty((M, N), dtype=float)

    for j, r in enumerate(r_grid):
        V = (4.0 / 3.0) * np.pi * r ** 3 * AR
        # r_eff(θ) for this radius bin — shape (K,)
        r_eff  = r * np.sqrt(1.0 + (AR ** 2 - 1.0) * cos_t ** 2)
        # qr_eff[i, k] = q[i] * r_eff[k] — shape (M, K)
        qr_eff = q[:, np.newaxis] * r_eff[np.newaxis, :]
        F_eff  = _sphere_amplitude(qr_eff)   # (M, K)
        G[:, j] = V * (F_eff ** 2 @ weights)

    return G * (contrast * 1e-4)


#: Vectorized G-matrix builders, keyed by shape name.
_G_BUILDERS: dict[str, callable] = {
    'sphere':   _build_g_sphere,
    'spheroid': _build_g_spheroid,
}


def build_g_matrix(
    q: np.ndarray,
    r_grid: np.ndarray,
    shape: str = 'sphere',
    contrast: float = 1.0,
    **shape_params,
) -> np.ndarray:
    """
    Build the instrument/shape-independent G matrix for size distribution fitting.

    G[i, j] = FF²(q[i], r[j]) × contrast × 1e-4   [cm^-1 per unit vol-frac]

    Fully vectorised — no Python loop over radius bins.

    Args:
        q:           1-D array of Q values [Å^-1], shape (M,)
        r_grid:      1-D array of radius bin centres [Å], shape (N,)
        shape:       Particle shape: 'sphere' or 'spheroid'
        contrast:    (Δρ)² in units of 10^20 cm^-4
        **shape_params: Extra parameters for the form factor.
                        For spheroid: aspect_ratio=<float>

    Returns:
        G matrix of shape (M, N)

    Raises:
        ValueError: if shape is not recognised.
    """
    q      = np.asarray(q,      dtype=float)
    r_grid = np.asarray(r_grid, dtype=float)

    if shape not in _G_BUILDERS:
        raise ValueError(
            f"Unknown shape '{shape}'. Supported: {sorted(_G_BUILDERS)}"
        )

    return _G_BUILDERS[shape](q, r_grid, contrast, **shape_params)


# ──────────────────────────────────────────────────────────────────────────────
# Radius grid helpers
# ──────────────────────────────────────────────────────────────────────────────

def make_r_grid(
    r_min: float,
    r_max: float,
    n_bins: int,
    log_spacing: bool = False,
) -> np.ndarray:
    """
    Create a radius bin-centre grid.

    Args:
        r_min:       Minimum radius [Å]
        r_max:       Maximum radius [Å]
        n_bins:      Number of bins
        log_spacing: If True, bins are equally spaced in log(r)

    Returns:
        1-D array of bin centres [Å], shape (n_bins,)
    """
    if log_spacing:
        return np.logspace(np.log10(r_min), np.log10(r_max), n_bins)
    return np.linspace(r_min, r_max, n_bins)


def bin_widths(r_grid: np.ndarray) -> np.ndarray:
    """
    Compute trapezoidal bin widths at each radius, matching the Igor Pro
    function IR1R_BinWidthInRadia.

    Rules:
        First point:    Δr[0]    = r[1]   - r[0]
        Last point:     Δr[-1]   = r[-1]  - r[-2]
        Interior:       Δr[i]    = (r[i+1] - r[i-1]) / 2

    Args:
        r_grid: 1-D array of radius bin centres [Å]

    Returns:
        1-D array of bin widths [Å], same shape as r_grid
    """
    r = np.asarray(r_grid, dtype=float)
    n = len(r)
    dw = np.empty(n, dtype=float)

    if n == 1:
        dw[0] = 1.0
        return dw

    dw[0] = r[1] - r[0]
    dw[-1] = r[-1] - r[-2]
    if n > 2:
        dw[1:-1] = (r[2:] - r[:-2]) / 2.0

    return dw
