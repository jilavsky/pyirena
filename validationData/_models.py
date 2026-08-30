"""Independent analytic scattering models used to synthesise validation data.

**This module deliberately does not import pyIrena.**  Every expression below is
written directly from the published literature (references in each docstring),
so that fitting the resulting data with pyIrena tests pyIrena's mathematics
rather than merely testing that pyIrena can invert its own forward model.

Conventions (identical to pyIrena and Igor Irena)
------------------------------------------------
* Q             : scattering vector, 1/Angstrom
* I             : absolute intensity, 1/cm
* radii, Rg     : Angstrom
* contrast      : (delta-rho)^2 in units of 1e20 cm^-4
* volume        : Angstrom^3, converted to cm^3 by the 1e-4 factor that appears
                  with contrast (1 A^3 = 1e-24 cm^3; 1e-24 * 1e20 = 1e-4)

So for a dilute assembly of particles with volume-fraction density f_V(r):

    I(Q) = 1e-4 * contrast * INTEGRAL V(r) * F^2(Q, r) * f_V(r) dr      [1/cm]
"""

from __future__ import annotations

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import erf as _erf
from scipy.special import gammaln as _gammaln

__all__ = [
    "sphere_ff2", "spheroid_ff2", "sphere_volume", "spheroid_volume",
    "lognormal_pdf", "gauss_pdf", "schulz_zimm_pdf",
    "size_distribution_intensity", "hard_sphere_sf",
    "unified_level", "guinier", "guinier_rod", "guinier_sheet",
    "porod", "power_law", "sphere_model", "debye_polymer_chain",
    "debye_bueche", "guinier_porod", "mass_fractal", "surface_fractal",
    "gauss_peak", "lorentz_peak", "pseudo_voigt_peak",
    "slit_smear",
]


# ---------------------------------------------------------------------------
# Form factors
# ---------------------------------------------------------------------------

def sphere_ff2(q, r):
    """Normalised sphere form factor squared, F^2(x) with x = Q*r.

    F(x) = 3 (sin x - x cos x) / x^3,  F(0) = 1.
    Rayleigh (1911); see e.g. Guinier & Fournet (1955).
    """
    x = np.asarray(q, dtype=float) * float(r)
    out = np.empty_like(x)
    small = x < 1e-3
    x2 = x[small] ** 2
    out[small] = 1.0 - x2 / 10.0 + x2 * x2 / 280.0      # Taylor of F(x)
    xb = x[~small]
    out[~small] = 3.0 * (np.sin(xb) - xb * np.cos(xb)) / xb ** 3
    return out ** 2


def spheroid_ff2(q, r, beta, n_quad=200):
    """Orientationally averaged spheroid F^2, semi-axes (r, r, beta*r).

    P(Q) = INTEGRAL_0^1 F^2(Q * r * sqrt(1 + (beta^2-1) u^2)) du,  u = cos(theta).
    Feigin & Svergun (1987), eq. for ellipsoids of revolution.
    Evaluated by Gauss-Legendre quadrature on [0, 1].
    """
    xi, wi = leggauss(int(n_quad))
    u = (xi + 1.0) / 2.0
    w = wi / 2.0
    q = np.asarray(q, dtype=float)
    reff = float(r) * np.sqrt(1.0 + (float(beta) ** 2 - 1.0) * u ** 2)
    x = q[:, None] * reff[None, :]
    f = np.where(x < 1e-3,
                 1.0 - x ** 2 / 10.0 + x ** 4 / 280.0,
                 3.0 * (np.sin(x) - x * np.cos(x)) / np.where(x == 0, 1.0, x) ** 3)
    return (f ** 2) @ w


def sphere_volume(r):
    """V = (4/3) pi r^3  [A^3]."""
    return (4.0 / 3.0) * np.pi * np.asarray(r, dtype=float) ** 3


def spheroid_volume(r, beta):
    """V = (4/3) pi r^3 beta for semi-axes (r, r, beta*r)  [A^3]."""
    return (4.0 / 3.0) * np.pi * np.asarray(r, dtype=float) ** 3 * float(beta)


# ---------------------------------------------------------------------------
# Size distributions (probability densities in r, 1/Angstrom)
# ---------------------------------------------------------------------------

def lognormal_pdf(r, median, sigma, shift=0.0):
    """Two-parameter log-normal shifted by `shift`.

    f(r) = 1 / ((r-s) sigma sqrt(2 pi)) * exp(-[ln((r-s)/median)]^2 / (2 sigma^2))

    `median` is the median of (r - shift); `sigma` is the log-space standard
    deviation.  With shift = 0 this is the ordinary log-normal, and it is
    exactly the distribution Irena / pyIrena call "LogNormal" with
    (min_size = shift, mean_size = median, sdeviation = sigma).
    """
    r = np.asarray(r, dtype=float)
    x = r - float(shift)
    out = np.zeros_like(r)
    m = x > 0
    s = float(sigma)
    out[m] = (1.0 / (x[m] * s * np.sqrt(2.0 * np.pi))
              * np.exp(-np.log(x[m] / float(median)) ** 2 / (2.0 * s ** 2)))
    return out


def gauss_pdf(r, mean, sd):
    """Normal density, mean and standard deviation in Angstrom."""
    r = np.asarray(r, dtype=float)
    return np.exp(-(r - float(mean)) ** 2 / (2.0 * float(sd) ** 2)) / (
        float(sd) * np.sqrt(2.0 * np.pi))


def schulz_zimm_pdf(r, mean, width):
    """Schulz-Zimm (gamma) density with mean `mean` and std-dev `width`."""
    r = np.asarray(r, dtype=float)
    z = (float(mean) / float(width)) ** 2 - 1.0
    a = float(mean) / (z + 1.0)
    out = np.zeros_like(r)
    m = r > 0
    out[m] = np.exp((z + 1.0) * np.log(r[m] / a) - r[m] / a - np.log(r[m])
                    - _gammaln(z + 1.0))
    return out


def size_distribution_intensity(q, r_grid, fv, contrast, shape="sphere",
                                aspect_ratio=1.0):
    """I(Q) from a polydisperse assembly.

        I(Q) = 1e-4 * contrast * INTEGRAL V(r) F^2(Q,r) f_V(r) dr

    Args:
        q:        Q values [1/A]
        r_grid:   integration grid in radius [A] (dense; trapezoid rule)
        fv:       volume-fraction density f_V(r) on r_grid [1/A]
        contrast: (delta-rho)^2 in 1e20 cm^-4
        shape:    'sphere' or 'spheroid'
        aspect_ratio: spheroid beta
    Returns:
        I(Q) [1/cm]
    """
    q = np.asarray(q, dtype=float)
    r_grid = np.asarray(r_grid, dtype=float)
    fv = np.asarray(fv, dtype=float)

    if shape == "sphere":
        V = sphere_volume(r_grid)
        x = q[:, None] * r_grid[None, :]
        f = np.where(x < 1e-3,
                     1.0 - x ** 2 / 10.0 + x ** 4 / 280.0,
                     3.0 * (np.sin(x) - x * np.cos(x)) / np.where(x == 0, 1.0, x) ** 3)
        ff2 = f ** 2
    elif shape == "spheroid":
        V = spheroid_volume(r_grid, aspect_ratio)
        ff2 = np.stack([spheroid_ff2(q, rr, aspect_ratio) for rr in r_grid], axis=1)
    else:
        raise ValueError(f"unknown shape {shape!r}")

    integrand = ff2 * (V * fv)[None, :]
    return 1e-4 * float(contrast) * np.trapezoid(integrand, r_grid, axis=1)


def hard_sphere_sf(q, radius, phi):
    """Percus-Yevick hard-sphere structure factor.

    Analytic PY solution as given by Ashcroft & Lehner, Phys. Rev. 145 (1966) 83:

        S(Q) = 1 / (1 + 24 phi G(A) / A),   A = 2 Q R
        G(A) = alpha (sinA - A cosA)/A^2
             + beta  (2A sinA + (2-A^2) cosA - 2)/A^3
             + gamma {-A^4 cosA + 4[(3A^2-6) cosA + (A^3-6A) sinA + 6]}/A^5
        alpha = (1+2phi)^2/(1-phi)^4
        beta  = -6 phi (1+phi/2)^2/(1-phi)^4
        gamma = phi alpha / 2
    """
    q = np.asarray(q, dtype=float)
    phi = float(phi)
    if phi <= 0 or radius <= 0:
        return np.ones_like(q)
    d = (1.0 - phi) ** 4
    alpha = (1.0 + 2.0 * phi) ** 2 / d
    beta = -6.0 * phi * (1.0 + phi / 2.0) ** 2 / d
    gamma = phi * alpha / 2.0

    A = np.maximum(2.0 * q * float(radius), 1e-10)
    sA, cA = np.sin(A), np.cos(A)
    G = (alpha * (sA - A * cA) / A ** 2
         + beta * (2.0 * A * sA + (2.0 - A ** 2) * cA - 2.0) / A ** 3
         + gamma * (-A ** 4 * cA
                    + 4.0 * ((3.0 * A ** 2 - 6.0) * cA
                             + (A ** 3 - 6.0 * A) * sA + 6.0)) / A ** 5)
    return 1.0 / (1.0 + 24.0 * phi * G / A)


# ---------------------------------------------------------------------------
# Beaucage Unified model
# ---------------------------------------------------------------------------

def unified_level(q, G, Rg, B, P, RgCO=0.0, K=None):
    """One Beaucage Unified level.

    Beaucage, J. Appl. Cryst. 28 (1995) 717; 29 (1996) 134.

        I(Q) = G exp(-Q^2 Rg^2/3) + B / Q*^P  * exp(-Q^2 RgCO^2/3)
        Q*   = Q / [erf(K Q Rg / sqrt(6))]^3
        K    = 1.00 for P > 3, 1.06 otherwise (Irena convention)

    RgCO is the low-Q cut-off radius that suppresses this level's power law
    below the size of the next-smaller structural level.
    """
    q = np.asarray(q, dtype=float)
    if K is None:
        K = 1.0 if P > 3.0 else 1.06
    e = np.maximum(_erf(K * q * Rg / np.sqrt(6.0)), 1e-10)
    q_star = q / e ** 3
    guinier_term = G * np.exp(-q ** 2 * Rg ** 2 / 3.0)
    cutoff = np.exp(-RgCO ** 2 * q ** 2 / 3.0) if RgCO > 0 else 1.0
    return guinier_term + B / np.maximum(q_star, 1e-100) ** P * cutoff


def unified_B_from_G(G, Rg, P):
    """Hammouda's "linked B": the B that makes Guinier and power law join smoothly.

        B = G exp(-P/2) (3P/2)^(P/2) / Rg^P
    """
    return G * np.exp(-P / 2.0) * (1.5 * P) ** (P / 2.0) / Rg ** P


# ---------------------------------------------------------------------------
# Simple-fit models
# ---------------------------------------------------------------------------

def guinier(q, I0, Rg):
    """I = I0 exp(-Q^2 Rg^2 / 3).  Guinier (1939)."""
    q = np.asarray(q, dtype=float)
    return I0 * np.exp(-q ** 2 * Rg ** 2 / 3.0)


def guinier_rod(q, I0, Rc):
    """I = I0 exp(-Q^2 Rc^2 / 2) / Q — rod cross-section Guinier law."""
    q = np.asarray(q, dtype=float)
    return I0 * np.exp(-q ** 2 * Rc ** 2 / 2.0) / q


def guinier_sheet(q, I0, Rg):
    """I = I0 exp(-Q^2 Rg^2) / Q^2 — lamella thickness Guinier law."""
    q = np.asarray(q, dtype=float)
    return I0 * np.exp(-q ** 2 * Rg ** 2) / q ** 2


def porod(q, Kp, background=0.0):
    """I = Kp Q^-4 + bg.  Porod (1951)."""
    return Kp * np.asarray(q, dtype=float) ** -4 + background


def power_law(q, prefactor, exponent, background=0.0):
    """I = prefactor Q^-exponent + bg."""
    return prefactor * np.asarray(q, dtype=float) ** -float(exponent) + background


def sphere_model(q, scale, R):
    """I = scale |F(QR)|^2 — monodisperse sphere."""
    return float(scale) * sphere_ff2(q, R)


def debye_polymer_chain(q, scale, Rg):
    """Debye (1947) Gaussian-coil form factor: 2(exp(-x)-1+x)/x^2, x = Q^2 Rg^2."""
    x = np.asarray(q, dtype=float) ** 2 * float(Rg) ** 2
    return float(scale) * np.where(x < 1e-6, 1.0 - x / 3.0,
                                   2.0 * (np.exp(-x) - 1.0 + x) / np.where(x == 0, 1.0, x) ** 2)


def debye_bueche(q, prefactor, eta, corr_length):
    """Debye-Bueche (1949): I = prefactor eta^2 xi^3 / (1 + Q^2 xi^2)^2."""
    q = np.asarray(q, dtype=float)
    return (float(prefactor) * float(eta) ** 2 * float(corr_length) ** 3
            / (1.0 + q ** 2 * float(corr_length) ** 2) ** 2)


def guinier_porod(q, G, Rg1, s1, P):
    """Single-level Guinier-Porod model, Hammouda, J. Appl. Cryst. 43 (2010) 716.

        Q1 = sqrt((P - s1)(3 - s1)/2) / Rg1
        I(Q < Q1) = G Q^-s1 exp(-Q^2 Rg1^2 / (3 - s1))
        I(Q >= Q1) = D Q^-P
        D = G exp(-(P-s1)/2) [ (3-s1)(P-s1)/2 ]^((P-s1)/2) / Rg1^(P-s1)
    """
    q = np.asarray(q, dtype=float)
    Q1 = np.sqrt((P - s1) * (3.0 - s1) / 2.0) / Rg1
    D = (G * np.exp(-(P - s1) / 2.0)
         * ((3.0 - s1) * (P - s1) / 2.0) ** ((P - s1) / 2.0) / Rg1 ** (P - s1))
    lo = G / q ** s1 * np.exp(-q ** 2 * Rg1 ** 2 / (3.0 - s1))
    hi = D / q ** P
    return np.where(q < Q1, lo, hi)


# ---------------------------------------------------------------------------
# Fractal models (Teixeira, J. Appl. Cryst. 21 (1988) 781)
# ---------------------------------------------------------------------------

def mass_fractal(q, phi, R, Dv, ksi, eta, contrast):
    """Mass-fractal aggregate of spherical primary particles.

        I(Q) = phi contrast 1e-4 V [ bracket S_f(Q) + (1-eta)^2 ] F^2(QR)
        S_f(Q) = sin[(Dv-1) atan(Q ksi)] /
                 [ (Dv-1) Q ksi (1 + Q^2 ksi^2)^((Dv-1)/2) ]
        bracket = eta (RC/R)^3 (ksi/RC)^Dv,   RC = R sqrt(2) sqrt(1 + ChiS^2)
    with ChiS = 1 for spheres, so RC = 2R.
    """
    q = np.asarray(q, dtype=float)
    V = sphere_volume(R)
    RC = R * np.sqrt(2.0) * np.sqrt(1.0 + (2.0 + 1.0) / 3.0)     # = 2R for spheres
    bracket = eta * RC ** 3 / R ** 3 * (ksi / RC) ** Dv
    x = q * ksi
    sf = np.where(x < 1e-6, 1.0,
                  np.sin((Dv - 1.0) * np.arctan(x))
                  / ((Dv - 1.0) * np.where(x == 0, 1.0, x)
                     * (1.0 + x ** 2) ** ((Dv - 1.0) / 2.0)))
    return phi * contrast * 1e-4 * V * (bracket * sf + (1.0 - eta) ** 2) * sphere_ff2(q, R)


def surface_fractal(q, surface, Ds, ksi, contrast):
    """Surface-fractal scattering.

        I(Q) = pi contrast*1e20 ksi^4 *1e-32* surface * Gamma(5-Ds)
               * sin[(3-Ds) atan(Q ksi)]
               / [ (1 + Q^2 ksi^2)^((5-Ds)/2) Q ksi ]
    (unit bookkeeping as in Irena: contrast in 1e20 cm^-4, surface in cm^-1.)
    """
    q = np.asarray(q, dtype=float)
    x = q * ksi
    return (np.pi * contrast * 1e20 * ksi ** 4 * 1e-32 * surface
            * np.exp(_gammaln(5.0 - Ds))
            * np.where(x < 1e-6, 0.0,
                       np.sin((3.0 - Ds) * np.arctan(x))
                       / ((1.0 + x ** 2) ** ((5.0 - Ds) / 2.0)
                          * np.where(x == 0, 1.0, x))))


# ---------------------------------------------------------------------------
# Diffraction peaks
# ---------------------------------------------------------------------------

def gauss_peak(q, A, Q0, FWHM):
    """Gaussian of height A at Q0 with full width at half maximum FWHM."""
    sigma = FWHM / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    return A * np.exp(-0.5 * ((np.asarray(q, dtype=float) - Q0) / sigma) ** 2)


def lorentz_peak(q, A, Q0, FWHM):
    """Lorentzian of height A at Q0 with full width at half maximum FWHM."""
    g = FWHM / 2.0
    return A * g ** 2 / ((np.asarray(q, dtype=float) - Q0) ** 2 + g ** 2)


def pseudo_voigt_peak(q, A, Q0, FWHM, eta):
    """eta * Lorentz + (1 - eta) * Gauss, both of height A and width FWHM."""
    return eta * lorentz_peak(q, A, Q0, FWHM) + (1.0 - eta) * gauss_peak(q, A, Q0, FWHM)


# ---------------------------------------------------------------------------
# Slit smearing
# ---------------------------------------------------------------------------

def slit_smear(model_fn, q, slit_length, n_l=4001):
    """Infinite-slit-length (Lake) smearing of an analytic model.

        I_smeared(Q) = (1/SL) INTEGRAL_0^SL I_ideal(sqrt(Q^2 + l^2)) dl

    Because `model_fn` is analytic, the integrand is evaluated exactly at every
    quadrature node -- no interpolation or extrapolation of a tabulated curve is
    involved, which makes this an independent reference for pyIrena's smearing.
    """
    q = np.asarray(q, dtype=float)
    ll = np.linspace(0.0, float(slit_length), int(n_l))
    qq = np.sqrt(q[:, None] ** 2 + ll[None, :] ** 2)
    vals = model_fn(qq.ravel()).reshape(qq.shape)
    return np.trapezoid(vals, ll, axis=1) / float(slit_length)
