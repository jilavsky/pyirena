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

    Handles the limit Qr → 0 via a Taylor expansion to avoid division by zero.
    """
    out = np.empty_like(qr, dtype=float)
    small = qr < 1e-3   # Taylor expansion avoids catastrophic cancellation
    # Taylor: F ≈ 1 - (Qr)²/10 + ...
    qr2 = qr[small] ** 2
    out[small] = 1.0 - qr2 / 10.0 + qr2 ** 2 / 280.0
    x = qr[~small]
    out[~small] = 3.0 * (np.sin(x) - x * np.cos(x)) / x ** 3
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
# G-matrix builder
# ──────────────────────────────────────────────────────────────────────────────

#: Registry of supported shapes and their form factor functions.
#: Each entry maps shape name → callable(q, r, **shape_params) → array [Å^6].
_FF_REGISTRY: dict[str, callable] = {
    'sphere':   lambda q, r, **kw: sphere_ff(q, r),
    'spheroid': lambda q, r, aspect_ratio=1.0, **kw: spheroid_ff(q, r, aspect_ratio),
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

    G[i, j] = FF²(q[i], r[j]) × contrast × 1e20

    Units:  [Å^6] × [10^20 cm^-4] × [1e20] = cm^-1 (per unit volume fraction)

    The contrast should be provided in the same units as (Δρ)², i.e. in
    10^20 cm^-4.  Multiplying by 1e20 converts to absolute cm^-1 units.

    Args:
        q:           1-D array of Q values [Å^-1], shape (M,)
        r_grid:      1-D array of radius bin centres [Å], shape (N,)
        shape:       Particle shape: 'sphere' or 'spheroid'
        contrast:    (Δρ)² in units of 10^20 cm^-4
        **shape_params: Extra parameters forwarded to the form factor function.
                        For spheroid: aspect_ratio=<float>

    Returns:
        G matrix of shape (M, N)

    Raises:
        ValueError: if shape is not in the registry.
    """
    q = np.asarray(q, dtype=float)
    r_grid = np.asarray(r_grid, dtype=float)

    if shape not in _FF_REGISTRY:
        raise ValueError(
            f"Unknown shape '{shape}'. Supported: {sorted(_FF_REGISTRY)}"
        )

    ff_func = _FF_REGISTRY[shape]
    M, N = len(q), len(r_grid)
    G = np.empty((M, N), dtype=float)

    for j, r in enumerate(r_grid):
        G[:, j] = ff_func(q, r, **shape_params)

    # Unit conversion:
    #   ff_func returns [Å³] = V(r) × F_norm²
    #   contrast in [10^20 cm^-4]
    #   1 Å³ = 1e-24 cm³,  so contrast × 1e20 × 1e-24 = contrast × 1e-4
    G *= contrast * 1e-4   # → [cm^-1] per unit volume fraction

    return G


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
