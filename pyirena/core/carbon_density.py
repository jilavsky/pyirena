"""
Carbon model — density / SLD / contrast convenience layer.

The Carbon model tool's value beyond the raw scattering formulas is this
layer: the user types a chemical formula once, the WAXS tab's fitted peak
positions give the interlayer and in-plane spacings, and everything else —
structural density, sample density, both scattering length densities, both
contrasts — follows without anyone re-typing a number between tabs.

Source
------
D. Saurel, J. Segalini, M. Jauregui, A. Pendashteh, B. Daffos, P. Simon,
M. Casas-Cabanas, *Energy Storage Materials* **21** (2019) 162–173,
supplementary information eq. (S1)–(S2), as corrected by the corrigendum,
*Energy Storage Materials* **28** (2020) 418.

Units (pyIrena convention throughout — see AGENTS.md §5)
--------------------------------------------------------
lengths Å · density g/cm³ · SLD 10¹⁰ cm⁻² · contrast (Δρ)² 10²⁰ cm⁻⁴.

No SLD physics is implemented here.  Everything goes through
:mod:`pyirena.core.scattering_contrast`, which already does formula + density
→ X-ray SLD against the NIST/Chantler tables.
"""

from __future__ import annotations

import logging
from functools import lru_cache
from typing import Optional

log = logging.getLogger(__name__)

# ── Graphite reference values ───────────────────────────────────────────────
#
# These are the denominators of the rho_struc ratio, so what matters is that
# they are the *same kind of quantity* as the numbers the user's WAXS fit
# produces — which are d-spacings, d = 2π/Q_c.
#
# d002: half the c axis of hexagonal graphite (c = 6.708 Å).  Franklin 1951.
# d100: the (100) d-spacing, a·√3/2 with a = 2.4612 Å  →  2.1315 Å.
#
#   The source paper prints "0.246 nm" for d100_graphite.  That is the a
#   lattice parameter, not the (100) d-spacing — the graphite (100) reflection
#   sits at Q ≈ 2.948 Å⁻¹, i.e. d = 2.131 Å.  Since rho_struc uses the ratio
#   d100/d100_graphite and the numerator is read off a fitted peak as 2π/Q_c,
#   the denominator has to be the d-spacing too, or the ratio is wrong by
#   √3/2 ≈ 0.866 (a 33 % error in the density once squared).  Hence 2.1315.
#   Users comparing against the paper's table should keep this in mind.
RHO_GRAPHITE = 2.26       # g/cm³
D002_GRAPHITE = 3.354     # Å   (c/2, c = 6.708 Å)
D100_GRAPHITE = 2.1315    # Å   (a·√3/2, a = 2.4612 Å)


def rho_struc_from_spacings(
    d002: float,
    d100: float,
    rho_graphite: float = RHO_GRAPHITE,
    d002_graphite: float = D002_GRAPHITE,
    d100_graphite: float = D100_GRAPHITE,
) -> float:
    """Structural (pore-free) density of a turbostratic carbon, g/cm³.

    ::

        ρ_struc = ρ_graphite · (d002_graphite / d002) · (d100_graphite / d100)²

    Derivation: the mass per hexagonal unit cell is fixed (4 C atoms), and the
    cell volume is V = (√3/2)·a²·c with a ∝ d100 and c = 2·d002.  Density is
    mass/volume, so ρ ∝ 1/(d100² · d002) — both spacing ratios sit in the
    numerator with the *graphite* value on top.  A carbon whose layers are
    further apart than graphite's (d002 > 3.354 Å, always true for disordered
    carbons) is correspondingly less dense.

    .. note::
       The planning transcription of the corrigendum had the in-plane ratio
       inverted, ``(d100/d100_graphite)²``, which would make a swollen lattice
       *denser*.  The form above is the one the unit-cell geometry gives and is
       what this code uses; ``pyirena/tests/test_carbon_fit.py`` pins both the
       graphite self-consistency case (ρ_struc = ρ_graphite when both spacings
       match) and the monotonic direction.

    Args:
        d002: Interlayer (002) spacing [Å], from the fitted (002) peak, 2π/Q_c.
        d100: In-plane (100) spacing [Å], from the fitted (100) peak, 2π/Q_c.
        rho_graphite: Reference graphite density [g/cm³].
        d002_graphite: Reference graphite (002) spacing [Å].
        d100_graphite: Reference graphite (100) spacing [Å].

    Returns:
        Structural density [g/cm³], or ``nan`` if either spacing is non-positive.
    """
    if not (d002 > 0 and d100 > 0 and d002_graphite > 0 and d100_graphite > 0):
        return float('nan')
    return float(rho_graphite * (d002_graphite / d002) * (d100_graphite / d100) ** 2)


def rho_sample(rho_struc: float, porosity: float) -> float:
    """Apparent (skeletal-plus-pores) density, g/cm³.

    ``ρ_sample = (1 − φ) · ρ_struc`` — SI eq. (S2).  ``φ`` is the micropore
    volume fraction from the SAXS region of the fit.
    """
    if not rho_struc > 0:
        return float('nan')
    return float(max(1.0 - float(porosity), 0.0) * rho_struc)


@lru_cache(maxsize=64)
def xray_sld_per_gram(formula: str) -> float:
    """X-ray SLD per unit density, [10¹⁰ cm⁻²] / (g/cm³), for one formula.

    The free-electron X-ray SLD is exactly proportional to density, so the
    whole formula-parsing and periodic-table lookup only has to happen once per
    formula — after which any density is a multiply.  That matters: the Carbon
    model recomputes its contrasts on every fit iteration, because they depend
    on the fitted peak positions and porosity.

    Returns ``nan`` rather than raising when the formula will not parse, since
    this is called live from the GUI while the user is still typing.
    """
    try:
        from pyirena.core.scattering_contrast import compute_compound
        # Any positive reference density works; the ratio is density-independent.
        return float(compute_compound(formula, 1.0).xray_sld_per_gram)
    except Exception:
        log.debug("xray_sld_per_gram(%r) failed", formula, exc_info=True)
        return float('nan')


def xray_sld(formula: str, density: float) -> float:
    """X-ray SLD [10¹⁰ cm⁻²] for a formula at a density, via NIST/Chantler.

    Thin wrapper over :func:`pyirena.core.scattering_contrast.compute_compound`
    (through the cached :func:`xray_sld_per_gram`) so the Carbon model never
    grows its own copy of the SLD physics.  Returns ``nan`` rather than raising
    when the formula will not parse or the density is non-positive.
    """
    if not (density and density > 0):
        return float('nan')
    per_gram = xray_sld_per_gram(formula)
    return float(per_gram * float(density))


def contrast_from_slds(sld_a: float, sld_b: float = 0.0) -> float:
    """(Δρ)² in 10²⁰ cm⁻⁴ between two SLDs given in 10¹⁰ cm⁻².

    The second phase defaults to vacuum (SLD = 0), which is what both of the
    Carbon model's contrasts need: grain-vs-air for the Porod term, and
    matrix-vs-empty-pore for the micropore term.
    """
    try:
        d = float(sld_a) - float(sld_b)
    except (TypeError, ValueError):
        return float('nan')
    return float(d * d)


def specific_surface_area_m2_g(S_per_volume: float, density: float) -> float:
    """Convert a per-volume surface area to the per-gram number BET reports.

    Args:
        S_per_volume: Surface area per unit volume [cm²/cm³ ≡ cm⁻¹], which is
            what a Porod term fitted against absolute intensity in cm⁻¹ gives.
        density: Density of the material the area belongs to [g/cm³].

    Returns:
        Specific surface area [m²/g].  ``1e-4`` converts cm²/g → m²/g.
    """
    if not (density and density > 0):
        return float('nan')
    return float(S_per_volume / density * 1e-4)


def scherrer_size(fwhm_q: float, shape_factor: float = 0.9) -> float:
    """Crystallite size [Å] from a peak's Gaussian FWHM in Q [Å⁻¹].

    Scherrer's relation written directly in Q rather than 2θ::

        L = 2π·K / ΔQ

    with ``K`` the shape factor (0.9 for the usual FWHM convention).  For the
    (002) reflection of a carbon this is L_c, the stack height.
    """
    if not (fwhm_q and fwhm_q > 0):
        return float('nan')
    return float(2.0 * 3.141592653589793 * shape_factor / fwhm_q)


def d_spacing(q_c: Optional[float]) -> float:
    """Bragg spacing d = 2π/Q_c [Å] from a peak centre in Å⁻¹."""
    if not (q_c and q_c > 0):
        return float('nan')
    return float(2.0 * 3.141592653589793 / q_c)
