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
from scipy.special import j1 as _scipy_j1


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


# ──────────────────────────────────────────────────────────────────────────────
# Cylinder form factor
# ──────────────────────────────────────────────────────────────────────────────

def _j1_over_x(x: np.ndarray) -> np.ndarray:
    """J₁(x)/x with Taylor expansion for x < 1e-3 to avoid cancellation.

    As x → 0:  J₁(x)/x → 1/2 − x²/16 + x⁴/384 − …
    Multiply by 2 to get the normalised cylinder radial amplitude, which → 1.
    """
    out = np.empty_like(x, dtype=float)
    small = x < 1e-3
    x2 = x[small] ** 2
    out[small] = 0.5 - x2 / 16.0 + x2 ** 2 / 384.0
    out[~small] = _scipy_j1(x[~small]) / x[~small]
    return out


def _build_g_cylinder_ar(
    q: np.ndarray,
    r_grid: np.ndarray,
    contrast: float,
    aspect_ratio: float = 1.0,
) -> np.ndarray:
    """G matrix for a finite cylinder with aspect-ratio parameterisation.

    The cylinder has equatorial radius r (from the distribution) and
    half-length L = AR * r, so total height H = 2·AR·r.

    Volume per bin:  V(r) = π r² · 2L = 2π · AR · r³

    The orientationally-averaged form factor squared is integrated via
    50-point Gauss-Legendre quadrature over cos(α) ∈ [0, 1]:

        <F²(q)> = ∫₀¹ [2J₁(qr√(1−u²)) / (qr√(1−u²))]² · sinc²(qLu/π) du

    where  sinc(x) = sin(πx)/(πx)  (NumPy convention).

    At AR → 0:  thin disk limit.
    At AR ≫ 1:  long rod limit.
    At AR = 1:  cylinder with H = 2R.
    """
    AR = float(aspect_ratio)
    cos_t   = 0.5 * (_GL_NODES + 1.0)   # cos(α) nodes ∈ [0, 1],  shape (K,)
    weights = 0.5 * _GL_WEIGHTS          # GL weights with Jacobian,  shape (K,)
    sin_t   = np.sqrt(np.maximum(1.0 - cos_t ** 2, 0.0))   # sin(α)

    M, N = len(q), len(r_grid)
    G = np.empty((M, N), dtype=float)

    for j, r in enumerate(r_grid):
        L = AR * r                            # half-length [Å]
        V = 2.0 * np.pi * AR * r ** 3        # cylinder volume [Å³]

        # Radial (perpendicular) amplitude: 2 J₁(qr sinα) / (qr sinα)
        qr_perp = q[:, np.newaxis] * (r * sin_t)        # (M, K)
        F_perp  = 2.0 * _j1_over_x(qr_perp)            # (M, K); → 1 at qr_perp=0

        # Axial (parallel) amplitude: sin(qL cosα) / (qL cosα) = sinc(qLu/π)
        qL_par = q[:, np.newaxis] * (L * cos_t)         # (M, K)
        F_par  = np.sinc(qL_par / np.pi)                # (M, K); → 1 at qL_par=0

        G[:, j] = V * ((F_perp * F_par) ** 2 @ weights)

    return G * (contrast * 1e-4)


def _build_g_cylinder_length(
    q: np.ndarray,
    r_grid: np.ndarray,
    contrast: float,
    length: float = 100.0,
) -> np.ndarray:
    """G matrix for a finite cylinder with fixed absolute length.

    The cylinder has equatorial radius r (from the distribution) and a
    fixed total height H = length [Å] (half-length L = length/2), so the
    volume is V(r) = π r² · length (scales as r², not r³).

    The form factor is identical to cylinder_ar except L is constant:

        <F²(q)> = ∫₀¹ [2J₁(qr√(1−u²)) / (qr√(1−u²))]² · sinc²(qLu/π) du

    Useful for disk-like scatterers (small length) or rod-like (large length)
    when the physical thickness/length is known independently of the radius.
    """
    L = float(length) / 2.0              # half-length [Å]
    cos_t   = 0.5 * (_GL_NODES + 1.0)
    weights = 0.5 * _GL_WEIGHTS
    sin_t   = np.sqrt(np.maximum(1.0 - cos_t ** 2, 0.0))

    M, N = len(q), len(r_grid)
    G = np.empty((M, N), dtype=float)

    # Axial factor is constant across all r bins (L is fixed)
    qL_par = q[:, np.newaxis] * (L * cos_t)             # (M, K)
    F_par  = np.sinc(qL_par / np.pi)                    # (M, K)

    for j, r in enumerate(r_grid):
        V = np.pi * r ** 2 * float(length)              # cylinder volume [Å³]

        qr_perp = q[:, np.newaxis] * (r * sin_t)        # (M, K)
        F_perp  = 2.0 * _j1_over_x(qr_perp)            # (M, K)

        G[:, j] = V * ((F_perp * F_par) ** 2 @ weights)

    return G * (contrast * 1e-4)


# ──────────────────────────────────────────────────────────────────────────────
# Core-Shell form factors
# ──────────────────────────────────────────────────────────────────────────────
#
# SLD convention: sld_core, sld_shell, sld_solvent are in units of 10⁻⁶ Å⁻².
# The contrast embedded in the amplitude is:
#   Δρ_c = (sld_core − sld_shell)   × 1e-6  [Å⁻²]
#   Δρ_s = (sld_shell − sld_solvent) × 1e-6  [Å⁻²]
#
# The scattering amplitude for a core-shell sphere (real, centrosymmetric):
#   F_cs = Δρ_c · V_core · f_sph(q·R_core)
#         + Δρ_s · V_total · f_sph(q·R_total)
#
# G[i,j] = F_cs²(q_i, R_c[j], R_t[j]) × 1e-4   [cm⁻¹ per unit vol-frac]
# The factor 1e-4 converts Å⁶ × Å⁻⁴ × 1e-12 × 1e24 = cm⁻¹.
# (Å³ × Å⁻² = Å, squared → Å²; 1e-6² = 1e-12; 1 Å² = 1e-20 cm²;
#  1 cm⁻¹ = 1 cm⁻¹. Net: [Å³ · 10⁻⁶Å⁻²]² = 10⁻¹² Å⁶ Å⁻⁴ = 10⁻¹² Å²
#  → × 1e-20+24 = ×1e4 → need ×1e-4 for Å³→cm³, giving ×1e-16 total
#  Combined pre-factor: 1e-4 [Å³→cm³] × 1e-12 [SLD²] = 1e-16 cm⁻¹)
#
# The `contrast` argument is accepted for API consistency but ignored.
#
# Polydispersity modes — what the r_grid axis represents:
#   by_core:  r_grid = R_core;  R_total = R_core + t_shell
#   by_shell: r_grid = t_shell; R_total = r_core_fixed + t_shell
#   by_total: r_grid = R_total; R_core  = R_total − t_shell


def _coreshell_f(
    q: np.ndarray,
    r_c: float,
    r_t: float,
    d_rho_c: float,
    d_rho_s: float,
) -> np.ndarray:
    """Core-shell sphere amplitude for a single (R_core, R_total) pair.

    Returns F_cs(q) [Å³ · Å⁻² · 10⁻⁶] — a 1-D array over q values.

    Args:
        q:       Q values [Å⁻¹], shape (M,)
        r_c:     Core radius [Å]
        r_t:     Total (outer) radius [Å]; must be ≥ r_c
        d_rho_c: ρ_core − ρ_shell  in 10⁻⁶ Å⁻²
        d_rho_s: ρ_shell − ρ_solvent in 10⁻⁶ Å⁻²
    """
    V_c = (4.0 / 3.0) * np.pi * r_c ** 3
    V_t = (4.0 / 3.0) * np.pi * r_t ** 3
    F_c = _sphere_amplitude(q * r_c)   # (M,)
    F_t = _sphere_amplitude(q * r_t)   # (M,)
    return d_rho_c * V_c * F_c + d_rho_s * V_t * F_t


def _cs_g_from_pairs(
    q: np.ndarray,
    r_c_arr: np.ndarray,
    r_t_arr: np.ndarray,
    sld_core: float,
    sld_shell: float,
    sld_solvent: float,
) -> np.ndarray:
    """Build G matrix from pre-computed (R_core, R_total) arrays.

    G[i,j] = F_cs²(q_i, r_c[j], r_t[j]) / V_total(r_t[j]) × 1e-4  [cm⁻¹]

    Unit derivation (SLDs in 10⁻⁶ Å⁻², volumes in Å³):
      F_cs has implicit units 10⁻⁶ Å⁻² × Å³ = 10⁻⁶ Å
      Physical amplitude F_A [cm] = F_cs × 10⁻¹⁴
      Intensity per Vf = |F_A|² / V_total = F_cs² × 10⁻²⁸ / (V_t × 10⁻²⁴)
                       = F_cs² / V_t × 10⁻⁴  [cm⁻¹]
    """
    d_rho_c = float(sld_core)   - float(sld_shell)    # 10⁻⁶ Å⁻²
    d_rho_s = float(sld_shell)  - float(sld_solvent)  # 10⁻⁶ Å⁻²

    M, N = len(q), len(r_c_arr)
    G = np.empty((M, N), dtype=float)
    for j in range(N):
        r_t = float(r_t_arr[j])
        V_t = (4.0 / 3.0) * np.pi * r_t ** 3          # total sphere volume [Å³]
        F = _coreshell_f(q, float(r_c_arr[j]), r_t, d_rho_c, d_rho_s)
        G[:, j] = F ** 2 / V_t
    return G * 1e-4   # convert to cm⁻¹ per unit vol-frac


def _build_g_cs_sphere_by_core(
    q: np.ndarray, r_grid: np.ndarray, contrast: float,
    sld_core: float = 10.0, sld_shell: float = 1.0,
    sld_solvent: float = 9.46, t_shell: float = 20.0,
) -> np.ndarray:
    """Core-Shell Sphere G matrix — distribution over core radius.

    r_grid = R_core;  R_total = R_core + t_shell (constant)
    """
    r_t = r_grid + float(t_shell)
    return _cs_g_from_pairs(q, r_grid, r_t, sld_core, sld_shell, sld_solvent)


def _build_g_cs_sphere_by_shell(
    q: np.ndarray, r_grid: np.ndarray, contrast: float,
    sld_core: float = 10.0, sld_shell: float = 1.0,
    sld_solvent: float = 9.46, r_core_fixed: float = 50.0,
) -> np.ndarray:
    """Core-Shell Sphere G matrix — distribution over shell thickness.

    r_grid = t_shell;  R_core = r_core_fixed (constant);  R_total = R_core + t_shell
    """
    r_c = np.full_like(r_grid, float(r_core_fixed))
    r_t = r_c + r_grid
    return _cs_g_from_pairs(q, r_c, r_t, sld_core, sld_shell, sld_solvent)


def _build_g_cs_sphere_by_total(
    q: np.ndarray, r_grid: np.ndarray, contrast: float,
    sld_core: float = 10.0, sld_shell: float = 1.0,
    sld_solvent: float = 9.46, t_shell: float = 20.0,
) -> np.ndarray:
    """Core-Shell Sphere G matrix — distribution over total outer radius.

    r_grid = R_total;  R_core = R_total − t_shell (constant shell thickness)
    Bins where R_total ≤ t_shell are treated as zero-shell (R_core = 0).
    """
    t = float(t_shell)
    r_c = np.maximum(r_grid - t, 0.0)
    return _cs_g_from_pairs(q, r_c, r_grid, sld_core, sld_shell, sld_solvent)


def _cs_spheroid_g_from_pairs(
    q: np.ndarray,
    r_c_arr: np.ndarray,
    r_t_arr: np.ndarray,
    sld_core: float,
    sld_shell: float,
    sld_solvent: float,
    aspect_ratio: float,
) -> np.ndarray:
    """Core-Shell Spheroid G matrix from (R_core, R_total) arrays.

    Orientational average via 50-point Gauss-Legendre quadrature over
    cos(α) ∈ [0, 1].  The spheroid semi-axes are (R, R, R·AR) for
    both the core and total particle, sharing the same AR.

    The core-shell spheroid amplitude at angle α:
      F(q, α) = Δρ_c · V_core · f_sph(q · r_eff(R_core, α))
               + Δρ_s · V_total· f_sph(q · r_eff(R_total, α))
    where r_eff(R, α) = R · √[1 + (AR²−1)·cos²(α)]

    The exact orientation average of the squared amplitude is:
      ⟨F²⟩ = ∫₀¹ F²(q, α) d(cosα)  ≈  Σₖ wₖ · F²(q, αₖ)
    """
    AR = float(aspect_ratio)
    cos_t   = 0.5 * (_GL_NODES + 1.0)   # (K,)
    weights = 0.5 * _GL_WEIGHTS          # (K,)

    d_rho_c = float(sld_core)  - float(sld_shell)
    d_rho_s = float(sld_shell) - float(sld_solvent)

    M, N = len(q), len(r_c_arr)
    G = np.empty((M, N), dtype=float)

    for j in range(N):
        R_c = float(r_c_arr[j])
        R_t = float(r_t_arr[j])
        V_c = (4.0 / 3.0) * np.pi * R_c ** 3 * AR
        V_t = (4.0 / 3.0) * np.pi * R_t ** 3 * AR

        # r_eff at each GL angle for core and total radii — shape (K,)
        stretch = np.sqrt(1.0 + (AR ** 2 - 1.0) * cos_t ** 2)
        r_eff_c = R_c * stretch   # (K,)
        r_eff_t = R_t * stretch   # (K,)

        # Sphere amplitudes — shape (M, K)
        F_c = _sphere_amplitude(q[:, np.newaxis] * r_eff_c)
        F_t = _sphere_amplitude(q[:, np.newaxis] * r_eff_t)

        # Core-shell amplitude at each angle — (M, K)
        F_cs = d_rho_c * V_c * F_c + d_rho_s * V_t * F_t

        G[:, j] = (F_cs ** 2 @ weights) / V_t

    return G * 1e-4


def _build_g_cs_spheroid_by_core(
    q: np.ndarray, r_grid: np.ndarray, contrast: float,
    sld_core: float = 10.0, sld_shell: float = 1.0,
    sld_solvent: float = 9.46, t_shell: float = 20.0,
    aspect_ratio: float = 1.0,
) -> np.ndarray:
    """Core-Shell Spheroid G matrix — distribution over core radius.

    r_grid = R_core;  R_total = R_core + t_shell.
    Both core and shell share the same aspect ratio AR.
    """
    r_t = r_grid + float(t_shell)
    return _cs_spheroid_g_from_pairs(
        q, r_grid, r_t, sld_core, sld_shell, sld_solvent, aspect_ratio
    )


def _build_g_cs_spheroid_by_total(
    q: np.ndarray, r_grid: np.ndarray, contrast: float,
    sld_core: float = 10.0, sld_shell: float = 1.0,
    sld_solvent: float = 9.46, t_shell: float = 20.0,
    aspect_ratio: float = 1.0,
) -> np.ndarray:
    """Core-Shell Spheroid G matrix — distribution over total outer radius.

    r_grid = R_total;  R_core = R_total − t_shell.
    """
    t = float(t_shell)
    r_c = np.maximum(r_grid - t, 0.0)
    return _cs_spheroid_g_from_pairs(
        q, r_c, r_grid, sld_core, sld_shell, sld_solvent, aspect_ratio
    )


#: Vectorized G-matrix builders, keyed by shape name.
_G_BUILDERS: dict[str, callable] = {
    'sphere':               _build_g_sphere,
    'spheroid':             _build_g_spheroid,
    'cylinder_ar':          _build_g_cylinder_ar,
    'cylinder_length':      _build_g_cylinder_length,
    'cs_sphere_by_core':    _build_g_cs_sphere_by_core,
    'cs_sphere_by_shell':   _build_g_cs_sphere_by_shell,
    'cs_sphere_by_total':   _build_g_cs_sphere_by_total,
    'cs_spheroid_by_core':  _build_g_cs_spheroid_by_core,
    'cs_spheroid_by_total': _build_g_cs_spheroid_by_total,
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
        shape:       Particle shape: 'sphere', 'spheroid', 'cylinder_ar',
                     'cylinder_length', or any 'cs_*' core-shell variant.
        contrast:    (Δρ)² in units of 10^20 cm^-4.  Ignored for core-shell
                     shapes (SLDs encode the contrast internally).
        **shape_params: Extra parameters for the form factor.
                        For spheroid:           aspect_ratio=<float>
                        For cylinder_ar:        aspect_ratio=<float>  (L = AR·r)
                        For cylinder_length:    length=<float> [Å]
                        For cs_sphere_by_core:  sld_core, sld_shell, sld_solvent
                                                [10⁻⁶ Å⁻²], t_shell [Å]
                        For cs_sphere_by_shell: …, r_core_fixed [Å]
                        For cs_sphere_by_total: …, t_shell [Å]
                        For cs_spheroid_*:      above + aspect_ratio

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
