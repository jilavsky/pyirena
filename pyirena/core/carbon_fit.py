"""
Carbon model — full-range SAXS+WAXS fitting of disordered carbonaceous materials.

One model, fitted in one pass across the whole measured range (USAXS → SAXS →
WAXS, often five decades in Q), as the sum of three physically distinct
contributions::

    I(Q) = I_Porod(Q) + I_mp(Q) + I_waxs(Q) + background

* ``I_Porod`` — scattering from the outer surface of the powder grains, plus an
  optional nanometre-scale surface-roughness term.  Dominates the lowest Q.
* ``I_mp``   — the micropore structure, either dilute pores with an optional
  fractal-aggregate structure factor, or the Teubner-Strey two-phase model.
  Dominates the middle of the range.
* ``I_waxs`` — turbostratic stacking: one or more Voigt diffraction peaks with
  a Debye-Waller intensity factor, a powder-orientation factor, and optionally
  a crumpled-layer envelope.  Dominates the high-Q end.

The three cannot be fitted separately: they share the contrast, which is
computed from the structural density, which is computed from the *fitted* WAXS
peak positions, and their tails overlap in the middle of the range.  Fitting
them together is the point of the tool.

Source
------
D. Saurel, J. Segalini, M. Jauregui, A. Pendashteh, B. Daffos, P. Simon,
M. Casas-Cabanas, "A SAXS outlook on disordered carbonaceous materials for
electrochemical energy storage", *Energy Storage Materials* **21** (2019)
162–173, ``doi:10.1016/j.ensm.2019.05.007``, with the corrigendum
*Energy Storage Materials* **28** (2020) 418, ``doi:10.1016/j.ensm.2020.03.013``.
Equation numbers in the docstrings below are the paper's own; "A3.N" refers to
Annex 3 of the supplementary information.

Units
-----
Å and Å⁻¹ throughout, as everywhere else in pyIrena — the paper's nm values
need ×10 (lengths) or ×0.1 (Q) before they can be compared with anything here.
Intensity cm⁻¹, SLD 10¹⁰ cm⁻², contrast 10²⁰ cm⁻⁴, surface area cm²/cm³ (the
m²/g number BET reports is a derived quantity), ⟨δz²⟩ Å², density g/cm³.

Reuse
-----
Nothing here reimplements maths pyIrena already has.  The Teubner-Strey
formula comes from :mod:`pyirena.core.simple_fits`, the Voigt peak from
:mod:`pyirena.core.waxs_peakfit` (which gained a true Voigt for this tool),
SLD/contrast from :mod:`pyirena.core.scattering_contrast` via
:mod:`pyirena.core.carbon_density`, and fit quality from
:mod:`pyirena.core.fit_metrics`.  The Teixeira structure factor and the
Beaucage globule/discoid form factors *are* written out here, because the
paper's parameterisation differs from Irena's existing one (Irena's
``_mass_fractal_intensity`` uses an η-based normalisation with no Γ(D−1)) and
bending one onto the other would change numbers rather than share code.
"""

from __future__ import annotations

import logging
import warnings
from dataclasses import dataclass, field
from datetime import datetime
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.optimize import least_squares
from scipy.special import erf, gammaln

from pyirena.core import carbon_density as _cd
from pyirena.core.modeling import _SerialisableDataclass
from pyirena.core.simple_fits import _teubner_strey
from pyirena.core.waxs_peakfit import voigt_fwhm, voigt_peak

log = logging.getLogger(__name__)


class CarbonFitAborted(RuntimeError):
    """Raised out of a fit's progress callback to stop it early.

    The GUI's Stop button works by raising this from the callback it hands to
    :meth:`CarbonFitModel.fit`; the fit catches it, puts the starting
    parameters back, and returns an unsuccessful result rather than leaving the
    model holding whatever the optimiser happened to be trying.
    """


__all__ = [
    'CARBON_SAXS_MODES', 'CARBON_WAXS_ENVELOPES', 'CarbonFitAborted',
    'CarbonBackground', 'CarbonSaxsRegion', 'CarbonWaxsPeak',
    'CarbonWaxsRegion', 'CarbonMaterial', 'CarbonFitModel', 'CarbonFitResult',
    'porod_roughness_factor', 'teixeira_structure_factor',
    'globule_form_factor', 'discoid_form_factor', 'default_carbon_peaks',
]


# ─────────────────────────────────────────────────────────────────────────────
# Unit bookkeeping
# ─────────────────────────────────────────────────────────────────────────────
#
# Porod's law in cgs is I(Q) = 2π·(Δρ)²·S/Q⁴ with every quantity in cm.  This
# package keeps Q in Å⁻¹, contrast in 10²⁰ cm⁻⁴ and S in cm²/cm³, so
#
#   I [cm⁻¹] = 2π · Δρ²[10²⁰cm⁻⁴]·10²⁰ · S[cm⁻¹] · (Q[Å⁻¹]·10⁸)⁻⁴
#            = 2π · Δρ² · S · 10⁻¹² / Q⁴
#
# and the same 10⁻¹² appears wherever a surface area is read back out of a
# fitted Porod prefactor.  ``_VOL_UNIT`` is the matching factor for a volume:
# Δρ²[10²⁰cm⁻⁴] × V[Å³] → cm⁻¹ costs 10²⁰ × 10⁻²⁴ = 10⁻⁴, which is the same
# convention the Modeling G-matrix uses.
_POROD_UNIT = 1.0e-12
_VOL_UNIT = 1.0e-4

#: Selectable micropore (SAXS-region) models.  A registry rather than an
#: if-chain so adding a model is one entry here plus one evaluator branch —
#: more are expected, and which one fits is sample-dependent.
CARBON_SAXS_MODES: Dict[str, str] = {
    'fractal': 'Pores + fractal aggregation',
    'teubner_strey': 'Teubner-Strey (two-phase)',
}

#: Selectable WAXS-region envelopes multiplying the summed diffraction peaks.
CARBON_WAXS_ENVELOPES: Dict[str, str] = {
    'none': 'Flat layers (peaks only)',
    'crumpled': 'Crumpled / curved layers',
}


# ─────────────────────────────────────────────────────────────────────────────
# Building blocks — pure functions, each independently testable
# ─────────────────────────────────────────────────────────────────────────────

def porod_roughness_factor(q: np.ndarray, R_rough: float) -> np.ndarray:
    """Surface-roughness envelope of the grain Porod term — eq. (3).

    ::

        f_rough(Q, R) = (2/9)·R⁴ / [ 1 + (1/5)(QR)² + (2/9)(QR)⁴ ]

    Limits that make this the right function: ``f_rough(0, R) = (2/9)R⁴``
    (finite, so the roughness adds nothing at Q → 0) and, for QR ≫ 1, the
    quartic term dominates the denominator and ``f_rough → Q⁻⁴``, so the full
    bracket of eq. (3) becomes ``(S_macro + S_rough)·Q⁻⁴`` — the paper's own
    eq. (4).  The roughness is therefore invisible below Q ≈ 1/R and simply
    adds its surface area to the Porod law above it.

    Args:
        q:       Scattering vector [Å⁻¹].
        R_rough: Roughness correlation length [Å].

    Returns:
        Envelope in Å⁴, same shape as ``q``.
    """
    R = max(float(R_rough), 1e-10)
    x2 = (np.asarray(q, dtype=float) * R) ** 2
    denom = 1.0 + x2 / 5.0 + (2.0 / 9.0) * x2 * x2
    return (2.0 / 9.0) * R ** 4 / denom


def teixeira_structure_factor(
    q: np.ndarray, D: float, sigma: float, R: float,
) -> np.ndarray:
    """Fractal-aggregate structure factor — eq. (5), Annex 3 eq. (A3.17).

    Teixeira (1988) *J. Appl. Cryst.* **21**, 781, in the form the paper
    writes it::

        S(Q) = 1 + (D·Σ^D / R^D)·Γ(D−1)·sin[(D−1)·arctan(QΣ)]
                   / [ QΣ·(1 + (QΣ)²)^((D−1)/2) ]

    ``D`` is the mass-fractal dimension of the aggregate, ``Σ`` the cutoff
    length above which the fractal ordering is lost, and ``R`` the radius of
    the primary unit being aggregated (the pore radius in the SAXS region, the
    layer transition radius in the crumpled WAXS envelope).  ``S → 1`` at high
    Q, and ``S(0) = 1 + D·Γ(D)·(Σ/R)^D``.

    .. note::
       This is *not* the same normalisation as Irena's
       ``modeling._mass_fractal_intensity``, which folds in a packing factor η
       and drops Γ(D−1).  Both are "Teixeira"; only this one reproduces the
       numbers in the source paper, so the two are deliberately separate.

    Args:
        q:     Scattering vector [Å⁻¹].
        D:     Mass fractal dimension, clipped to (1, 3).
        sigma: Fractal cutoff length Σ [Å].
        R:     Primary-unit radius [Å].

    Returns:
        S(Q), same shape as ``q``.
    """
    q = np.asarray(q, dtype=float)
    D = float(np.clip(D, 1.001, 2.999))
    S = max(float(sigma), 1e-10)
    R = max(float(R), 1e-10)

    x = q * S
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        # Γ(D−1) through gammaln: D−1 can approach 0 where Γ blows up.
        pref = D * np.exp(gammaln(D - 1.0)) * (S / R) ** D
        num = np.sin((D - 1.0) * np.arctan(x))
        den = x * (1.0 + x * x) ** ((D - 1.0) / 2.0)
        # x → 0: sin((D−1)arctan x)/x → (D−1), so the term → D·Γ(D)·(Σ/R)^D.
        ratio = np.where(x < 1e-8, D - 1.0, num / np.maximum(den, 1e-300))
    out = 1.0 + pref * ratio
    return np.where(np.isfinite(out), out, 1.0)


def globule_form_factor(q: np.ndarray, r: float, k: float = 1.0) -> np.ndarray:
    """"Algebraic globule" pore form factor — eq. (6), Annex 3 eq. (A3.8/A3.9).

    ::

        P(Q) = exp(−(Qr)²/5) + [erf(Qr/√10)]¹² · (9/2)·k / (Qr)⁴

    This is Beaucage's unified level for a compact 3-D object written out for a
    sphere of radius ``r``: the Guinier term uses Rg² = 3r²/5, and the
    erf-blended Porod term has exponent P = 4 with prefactor B = (9/2)k/r⁴.
    ``k`` is the paper's shape allowance — ``k = 1`` is a monodisperse sphere,
    ``k > 1`` a globule of less well-defined shape.  ``P(0) = 1``, so the
    overall scale sits entirely in ``I₀``.

    Args:
        q: Scattering vector [Å⁻¹].
        r: Pore radius [Å].
        k: Globule shape factor (1 = sphere).

    Returns:
        P(Q), same shape as ``q``, normalised to 1 at Q = 0.
    """
    q = np.asarray(q, dtype=float)
    r = max(float(r), 1e-10)
    x = q * r
    guinier = np.exp(-(x ** 2) / 5.0)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        porod = np.where(
            x < 1e-6,
            0.0,
            erf(x / np.sqrt(10.0)) ** 12 * 4.5 * float(k) / np.maximum(x ** 4, 1e-300),
        )
    out = guinier + porod
    return np.where(np.isfinite(out), out, 0.0)


def discoid_form_factor(q: np.ndarray, R: float) -> np.ndarray:
    """"Unified discoid" layer form factor — eq. (17), Annex 3 eq. (A3.14/A3.47).

    ::

        P(Q) = exp(−(QR)²/6) + [erf(1.06·QR/√12)]⁶ · 2/(QR)²

    Beaucage's unified level for a thin disc of radius ``R``: Guinier with
    Rg² = R²/2, then an erf-blended Q⁻² power law — the signature of a locally
    two-dimensional (sheet-like) object.  Used only by the crumpled-layer WAXS
    envelope, where ``R`` is the radius over which a graphene layer stays flat
    before it bends.

    Args:
        q: Scattering vector [Å⁻¹].
        R: Layer transition radius [Å].

    Returns:
        P(Q), same shape as ``q``, normalised to 1 at Q = 0.
    """
    q = np.asarray(q, dtype=float)
    R = max(float(R), 1e-10)
    x = q * R
    guinier = np.exp(-(x ** 2) / 6.0)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        power = np.where(
            x < 1e-6,
            0.0,
            erf(1.06 * x / np.sqrt(12.0)) ** 6 * 2.0 / np.maximum(x ** 2, 1e-300),
        )
    out = guinier + power
    return np.where(np.isfinite(out), out, 0.0)


def _coherence_length(fwhm_g: float, fwhm_l: float,
                      shape_factor: float = 0.9) -> float:
    """Scherrer size [Å] from a Voigt peak's Gaussian component.

    Only the Gaussian part is finite-size broadening; the Lorentzian part is
    layer curvature, so the size must not be read off the observed width.

    Returns ``nan`` when the Gaussian component is less than a tenth of the
    observed width.  The two Voigt components are strongly correlated for a
    weak or poorly resolved peak, and the optimiser will happily collapse the
    Gaussian to zero while the Lorentzian absorbs the whole profile — the
    total width stays right, but ``2π·K/FWHM_G`` then reports a coherence
    length of tens of thousands of Å for a carbon whose layers are nanometres
    across.  A blank is the honest answer: the size is simply not determined
    by that peak.

    Args:
        fwhm_g: Gaussian component FWHM [Å⁻¹].
        fwhm_l: Lorentzian component FWHM [Å⁻¹].
        shape_factor: Scherrer K — 0.9 for a stack height L_c, 1.84 (Warren)
            for the two-dimensional hk band that gives the layer extent L_a.
    """
    total = voigt_fwhm(fwhm_g, fwhm_l)
    if not (total > 0) or fwhm_g <= 0 or fwhm_g < 0.1 * total:
        return float('nan')
    return _cd.scherrer_size(fwhm_g, shape_factor)


# ─────────────────────────────────────────────────────────────────────────────
# Model sections
# ─────────────────────────────────────────────────────────────────────────────
#
# Every fittable quantity follows Modeling's triple convention —
# ``X`` / ``fit_X`` / ``X_limits`` — so one generic packer can walk the whole
# model, and a new parameter is serialised the moment it is declared
# (_SerialisableDataclass walks dataclasses.fields()).

@dataclass
class CarbonBackground(_SerialisableDataclass):
    """Complex background: grain Porod scattering + optional surface roughness.

    ``I_Porod(Q) = 2π(Δρ)²·[ S_macro·Q⁻ⁿ + S_rough·f_rough(Q, R_rough) ]``
    (eq. 3), plus a flat instrumental background.

    ``S_macro`` is a specific surface area in cm²/cm³ only while the exponent
    ``n`` is exactly 4; away from 4 it is a Porod prefactor and the derived
    surface areas are reported as ``nan``.  The exponent is exposed because
    real grain surfaces are rarely perfectly sharp, but it is fixed at 4 by
    default — freeing it and ``S_rough`` together is usually over-fitting.
    """
    enabled: bool = True

    S_macro: float = 1.0e4                 # cm²/cm³
    fit_S_macro: bool = True
    S_macro_limits: tuple = (0.0, 1.0e14)

    porod_exponent: float = 4.0
    fit_porod_exponent: bool = False
    porod_exponent_limits: tuple = (2.0, 4.5)

    use_roughness: bool = False
    S_rough: float = 1.0e4                 # cm²/cm³
    fit_S_rough: bool = True
    S_rough_limits: tuple = (0.0, 1.0e14)
    R_rough: float = 50.0                  # Å
    fit_R_rough: bool = True
    R_rough_limits: tuple = (1.0, 1.0e5)

    flat_background: float = 0.0           # cm⁻¹
    fit_flat_background: bool = True
    flat_background_limits: tuple = (0.0, 1.0e6)


@dataclass
class CarbonSaxsRegion(_SerialisableDataclass):
    """Micropore scattering, ``I_mp`` — one of two mutually exclusive models.

    ``mode='fractal'`` (eq. 5–6): dilute pores of radius ``r`` described by the
    algebraic-globule form factor, optionally aggregated with a Teixeira mass
    fractal of dimension ``D`` and cutoff ``Σ``.  The scale is physical:
    ``I₀ = φ·(Δρ)²·(4/3)πr³``.

    ``mode='teubner_strey'`` (eq. 8): the two-phase semi-empirical form
    ``I₀/(1 + C1·Q² + C2·Q⁴)``, evaluated by ``simple_fits._teubner_strey``.
    Correlation length, repeat distance, amphiphilicity and the pore/wall
    widths all come out as derived quantities.

    The two are alternative descriptions of the same pore population and are
    never summed — switching modes switches which parameters are live.
    """
    enabled: bool = True
    mode: str = 'fractal'

    # ── fractal + globule branch ──
    phi: float = 0.10                      # pore volume fraction
    fit_phi: bool = True
    phi_limits: tuple = (1.0e-6, 0.95)

    pore_radius: float = 5.0               # r [Å]
    fit_pore_radius: bool = True
    pore_radius_limits: tuple = (1.0, 1.0e4)

    globule_k: float = 1.0
    fit_globule_k: bool = False
    globule_k_limits: tuple = (0.1, 20.0)

    use_fractal: bool = False              # False → dilute pores, S(Q) ≡ 1
    fractal_D: float = 2.5
    fit_fractal_D: bool = True
    fractal_D_limits: tuple = (1.01, 2.99)
    fractal_sigma: float = 100.0           # Σ [Å]
    fit_fractal_sigma: bool = True
    fractal_sigma_limits: tuple = (5.0, 1.0e5)

    # ── Teubner-Strey branch ──
    ts_I0: float = 1.0                     # cm⁻¹
    fit_ts_I0: bool = True
    ts_I0_limits: tuple = (1.0e-30, 1.0e10)
    ts_C1: float = -30.0                   # Å²
    fit_ts_C1: bool = True
    ts_C1_limits: tuple = (-1.0e8, 1.0e8)
    ts_C2: float = 5000.0                  # Å⁴
    fit_ts_C2: bool = True
    ts_C2_limits: tuple = (1.0e-30, 1.0e14)


@dataclass
class CarbonWaxsPeak(_SerialisableDataclass):
    """One turbostratic diffraction peak — eq. (2).

    A true Voigt (Lorentzian ⊗ Gaussian, not a linear pseudo-Voigt mix),
    evaluated by :func:`pyirena.core.waxs_peakfit.voigt_peak`.  The two widths
    have different physical origins and are therefore separate parameters:
    ``FWHM_G`` approximates finite-crystallite-size broadening (Annex 3 eq.
    A3.32 — the true shape is a squared sinc, of which a Gaussian is an
    excellent fit over the upper two thirds), and ``FWHM_L`` comes from
    distortions of the second kind, i.e. layer bending and curvature
    (Annex 3 eq. A3.38–A3.41, Vonk).

    ``K`` is the peak amplitude *before* the 1/Q² orientation factor and the
    Debye-Waller factor are applied; the observed height is reported as a
    derived quantity.  ``label`` is the Miller index the user recognises and is
    what the Material section looks up when it needs d002 or d100.
    """
    enabled: bool = True
    label: str = '002'

    Q0: float = 1.873                      # Å⁻¹ (d = 3.354 Å)
    fit_Q0: bool = True
    Q0_limits: tuple = (0.05, 25.0)

    K: float = 1.0
    fit_K: bool = True
    K_limits: tuple = (0.0, 1.0e12)

    FWHM_G: float = 0.30                   # Å⁻¹, crystallite size broadening
    fit_FWHM_G: bool = True
    FWHM_G_limits: tuple = (1.0e-4, 10.0)

    FWHM_L: float = 0.15                   # Å⁻¹, curvature broadening
    fit_FWHM_L: bool = True
    FWHM_L_limits: tuple = (0.0, 10.0)


@dataclass
class CarbonWaxsRegion(_SerialisableDataclass):
    """Everything the WAXS peaks share — eq. (2) and eq. (16)–(17).

    ``I_waxs(Q) = envelope(Q) · exp(−Q²⟨δz²⟩/3) · Σ_peaks V_i(Q) / Q²``

    * ``delta_z2`` — ⟨δz²⟩, the mean-square fluctuation of the interlayer
      spacing (distortions of the *first* kind).  It attenuates intensity
      without touching peak width or position, so it is a property of the
      material and is shared by every peak: that is precisely what makes it
      identifiable, because it sets the (004)/(002) intensity ratio.  Fitted
      per-peak it would be degenerate with that peak's ``K``.
    * ``use_orientation_factor`` — the 1/Q² powder average that ties a locally
      one-dimensional stacking correlation to the measured 3-D isotropic
      intensity (Annex 3 eq. A3.25).  It is easy to drop silently and a fit
      without it looks fine while reporting a wrong ``K``, so it is an explicit,
      testable switch rather than something buried in the peak shape.
    * ``envelope='crumpled'`` — eq. (16)–(17), for layers that are bent rather
      than forming flat nanocrystallites: a Teixeira fractal structure factor
      times a unified discoid form factor.  This is the same crumpling geometry
      the SAXS region may be seeing, hence the three ``link_*`` switches, which
      make the WAXS parameter follow the SAXS one instead of being fitted
      independently.  ``link_R_to_pore`` equates the layer transition radius
      with the pore radius — a strong physical assumption, off by default.
    """
    enabled: bool = True

    delta_z2: float = 0.0                  # ⟨δz²⟩ [Å²]
    fit_delta_z2: bool = False
    delta_z2_limits: tuple = (0.0, 10.0)

    use_orientation_factor: bool = True

    envelope: str = 'none'
    R_layer: float = 20.0                  # Å
    fit_R_layer: bool = True
    R_layer_limits: tuple = (1.0, 1.0e4)
    fractal_D: float = 2.5
    fit_fractal_D: bool = True
    fractal_D_limits: tuple = (1.01, 2.99)
    fractal_sigma: float = 100.0           # Å
    fit_fractal_sigma: bool = True
    fractal_sigma_limits: tuple = (5.0, 1.0e5)

    link_D_to_saxs: bool = False
    link_sigma_to_saxs: bool = False
    link_R_to_pore: bool = False


@dataclass
class CarbonMaterial(_SerialisableDataclass):
    """Composition and density → the two contrasts the other sections need.

    This is the wiring that makes the tool more than three formulas.  The WAXS
    tab's fitted (002) and (100) positions give the spacings, the spacings give
    the structural density (corrigendum eq. S1), the SAXS region's porosity
    gives the sample density (eq. S2), and
    :mod:`pyirena.core.scattering_contrast` turns each density plus the
    chemical formula into an X-ray SLD.  Nothing is re-typed between tabs.

    Two contrasts come out, and they are different:

    * ``contrast_porod`` — grain against vacuum, so it uses ``ρ_sample``
      (the grain *including* its pores).  Used by the background section.
    * ``contrast_micropore`` — empty pore against the carbon matrix, so it uses
      ``ρ_struc`` (the pore-free skeleton).  Used by the SAXS region.

    Each of the three stages can be overridden: ``rho_struc_mode='manual'`` when
    there is no usable (100) peak, ``porosity_mode='manual'`` for the
    Teubner-Strey branch (whose φ is derived, not fitted, so taking it live
    would make the contrast depend on its own output), and
    ``contrast_mode='manual'`` to bypass the chain entirely.
    """
    formula: str = 'C'

    rho_struc_mode: str = 'from_peaks'     # 'from_peaks' | 'manual'
    rho_struc_manual: float = 2.0          # g/cm³
    d002_label: str = '002'
    d100_label: str = '100'
    rho_graphite: float = _cd.RHO_GRAPHITE
    d002_graphite: float = _cd.D002_GRAPHITE
    d100_graphite: float = _cd.D100_GRAPHITE

    porosity_mode: str = 'auto'            # 'auto' | 'manual'
    porosity_manual: float = 0.0

    contrast_mode: str = 'auto'            # 'auto' | 'manual'
    contrast_porod_manual: float = 287.0       # 10²⁰ cm⁻⁴
    contrast_micropore_manual: float = 287.0   # 10²⁰ cm⁻⁴


def default_carbon_peaks() -> List[CarbonWaxsPeak]:
    """The three reflections a disordered-carbon pattern usually shows.

    (002) at d = 3.354 Å, (100) at d = 2.1315 Å and (004) at d = 1.677 Å —
    graphite positions, which is where a turbostratic carbon's peaks start
    before the fit moves them outwards.  (004) is disabled by default because
    it is often too weak to see.
    """
    return [
        CarbonWaxsPeak(label='002', Q0=2.0 * np.pi / _cd.D002_GRAPHITE,
                       K=1.0, FWHM_G=0.30, FWHM_L=0.15),
        CarbonWaxsPeak(label='100', Q0=2.0 * np.pi / _cd.D100_GRAPHITE,
                       K=0.3, FWHM_G=0.25, FWHM_L=0.10),
        CarbonWaxsPeak(label='004', Q0=4.0 * np.pi / _cd.D002_GRAPHITE,
                       K=0.05, FWHM_G=0.40, FWHM_L=0.20, enabled=False),
    ]


# ─────────────────────────────────────────────────────────────────────────────
# Parameter references — one generic walker over the whole model
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class _ParamRef:
    """One fittable scalar: where it lives, what it is called, what bounds it has.

    Holding a reference to the owning dataclass rather than a copy of the value
    means the fitter writes straight back into the model, so the model is
    always the single source of truth and the GUI needs no separate sync step.
    """
    owner: object
    attr: str
    key: str            # dotted, stable, used in JSON/HDF5: 'background.S_macro'
    label: str          # what the user reads
    unit: str = ''

    @property
    def value(self) -> float:
        return float(getattr(self.owner, self.attr))

    @value.setter
    def value(self, v: float) -> None:
        setattr(self.owner, self.attr, float(v))

    @property
    def fit(self) -> bool:
        return bool(getattr(self.owner, f'fit_{self.attr}', False))

    @property
    def limits(self) -> Tuple[float, float]:
        lo, hi = getattr(self.owner, f'{self.attr}_limits', (-np.inf, np.inf))
        return float(lo), float(hi)


@dataclass
class CarbonFitResult:
    """Outcome of one Carbon model fit — settings-free, results only.

    Arrays are kept here (not in the model) because the model is the *setup*:
    it is what round-trips to JSON and to the HDF5 ``_pyirena_config``, and it
    must stay free of data.
    """
    success: bool = False
    message: str = ''
    params: Dict[str, float] = field(default_factory=dict)
    errors: Dict[str, float] = field(default_factory=dict)
    derived: Dict[str, float] = field(default_factory=dict)
    chi_squared: float = float('nan')
    reduced_chi_squared: float = float('nan')
    n_points: int = 0
    n_params: int = 0
    n_iterations: int = 0
    timestamp: str = ''
    quality: Dict[str, object] = field(default_factory=dict)

    q: Optional[np.ndarray] = None
    I_data: Optional[np.ndarray] = None
    I_error: Optional[np.ndarray] = None
    I_model: Optional[np.ndarray] = None
    I_porod: Optional[np.ndarray] = None
    I_mp: Optional[np.ndarray] = None
    I_waxs: Optional[np.ndarray] = None
    residuals: Optional[np.ndarray] = None

    def to_dict(self) -> dict:
        """JSON-serialisable summary — scalars only, no arrays."""
        return {
            'success': bool(self.success),
            'message': str(self.message),
            'params': {k: float(v) for k, v in self.params.items()},
            'errors': {k: float(v) for k, v in self.errors.items()},
            'derived': {k: float(v) for k, v in self.derived.items()},
            'chi_squared': float(self.chi_squared),
            'reduced_chi_squared': float(self.reduced_chi_squared),
            'n_points': int(self.n_points),
            'n_params': int(self.n_params),
            'n_iterations': int(self.n_iterations),
            'timestamp': str(self.timestamp),
        }


# ─────────────────────────────────────────────────────────────────────────────
# The model
# ─────────────────────────────────────────────────────────────────────────────

class CarbonFitModel:
    """The whole Carbon model: four sections, N diffraction peaks, one fit.

    Composition rather than a flat parameter bag, because the three scattering
    contributions are genuinely separate physics that happen to share a
    contrast and a Q range.  ``to_dict()``/``from_dict()`` compose the sections'
    own serialisation, so a new field anywhere is persisted, scriptable and
    batch-runnable the moment it is declared.

    Example:
        >>> m = CarbonFitModel()
        >>> m.saxs.mode = 'teubner_strey'
        >>> q = np.logspace(-3, 0.5, 400)
        >>> I = m.evaluate(q)
        >>> comps = m.evaluate_components(q)
        >>> sorted(comps)
        ['background', 'mp', 'porod', 'total', 'waxs']
    """

    #: Bumped when a stored field changes meaning.  ``from_dict`` migrates.
    SCHEMA_VERSION = 1

    def __init__(self) -> None:
        self.background = CarbonBackground()
        self.saxs = CarbonSaxsRegion()
        self.waxs = CarbonWaxsRegion()
        self.peaks: List[CarbonWaxsPeak] = default_carbon_peaks()
        self.material = CarbonMaterial()

        # Fit range and numerics (not physics — kept flat on the model).
        self.q_min: float = 0.0
        self.q_max: float = 0.0            # 0 = "no limit", use all data
        self.weighting: str = 'auto'       # 'auto' | 'sigma' | 'relative' | 'log'
        self.max_iterations: int = 400
        self.n_mc_runs: int = 0            # 0 = no Monte-Carlo uncertainty pass

    # ── Peak bookkeeping ────────────────────────────────────────────────────

    def peak_by_label(self, label: str) -> Optional[CarbonWaxsPeak]:
        """First enabled peak whose label matches, case-insensitively.

        Igor names and Miller indices are matched loosely on purpose: users
        type ``002``, ``(002)`` and ``d002`` interchangeably.
        """
        if not label:
            return None
        want = str(label).strip().lower().lstrip('d').strip('()')
        for pk in self.peaks:
            if not pk.enabled:
                continue
            got = str(pk.label).strip().lower().lstrip('d').strip('()')
            if got == want:
                return pk
        return None

    def add_peak(self, label: str = '', Q0: float = 1.9) -> CarbonWaxsPeak:
        """Append a peak and return it."""
        pk = CarbonWaxsPeak(label=label or f'peak{len(self.peaks) + 1}', Q0=Q0)
        self.peaks.append(pk)
        return pk

    def remove_peak(self, index: int) -> None:
        """Drop peak ``index`` if it exists."""
        if 0 <= index < len(self.peaks):
            del self.peaks[index]

    # ── Linked parameters ───────────────────────────────────────────────────

    def effective_waxs_geometry(self) -> Tuple[float, float, float]:
        """(D, Σ, R) actually used by the crumpled-layer envelope.

        Resolves the three ``link_*`` switches: a linked parameter takes the
        SAXS region's value instead of its own, and is excluded from the fit
        vector so the optimiser never sees a duplicate degree of freedom.
        """
        w, s = self.waxs, self.saxs
        D = s.fractal_D if w.link_D_to_saxs else w.fractal_D
        sigma = s.fractal_sigma if w.link_sigma_to_saxs else w.fractal_sigma
        R = s.pore_radius if w.link_R_to_pore else w.R_layer
        return float(D), float(sigma), float(R)

    # ── Contrast chain ──────────────────────────────────────────────────────

    def resolve_material(self) -> Dict[str, float]:
        """Run the spacing → density → SLD → contrast chain once.

        Returns every intermediate as well as the two contrasts, because the
        Material tab displays all of them and the report lists all of them.
        Called on every fit iteration (the peak positions move), which is why
        the SLD lookup is cached per formula.

        Returns:
            dict with keys ``d002``, ``d100``, ``rho_struc``, ``porosity``,
            ``rho_sample``, ``sld_struc``, ``sld_sample``, ``contrast_porod``,
            ``contrast_micropore``.  Any link in the chain that cannot be
            evaluated is ``nan``; the contrasts then fall back to the manual
            values so a fit is still possible.
        """
        mat = self.material
        pk002 = self.peak_by_label(mat.d002_label)
        pk100 = self.peak_by_label(mat.d100_label)
        d002 = _cd.d_spacing(pk002.Q0) if pk002 is not None else float('nan')
        d100 = _cd.d_spacing(pk100.Q0) if pk100 is not None else float('nan')

        if mat.rho_struc_mode == 'manual':
            rho_struc = float(mat.rho_struc_manual)
        else:
            rho_struc = _cd.rho_struc_from_spacings(
                d002, d100, mat.rho_graphite, mat.d002_graphite, mat.d100_graphite)
            if not np.isfinite(rho_struc):
                # No usable (002)/(100) pair yet — fall back rather than fail,
                # so the panel still shows a curve before the peaks are set up.
                rho_struc = float(mat.rho_struc_manual)

        if mat.porosity_mode == 'manual':
            porosity = float(mat.porosity_manual)
        elif self.saxs.enabled and self.saxs.mode == 'fractal':
            porosity = float(self.saxs.phi)
        else:
            # Teubner-Strey's φ is *derived* from I₀, which itself depends on
            # the contrast — taking it live here would close a feedback loop
            # around the fit.  The user's estimate breaks it.
            porosity = float(mat.porosity_manual)
        porosity = float(np.clip(porosity, 0.0, 0.999))

        rho_sample = _cd.rho_sample(rho_struc, porosity)
        sld_struc = _cd.xray_sld(mat.formula, rho_struc)
        sld_sample = _cd.xray_sld(mat.formula, rho_sample)

        if mat.contrast_mode == 'manual':
            c_porod = float(mat.contrast_porod_manual)
            c_mp = float(mat.contrast_micropore_manual)
        else:
            c_porod = _cd.contrast_from_slds(sld_sample, 0.0)
            c_mp = _cd.contrast_from_slds(sld_struc, 0.0)
            if not np.isfinite(c_porod):
                c_porod = float(mat.contrast_porod_manual)
            if not np.isfinite(c_mp):
                c_mp = float(mat.contrast_micropore_manual)

        return {
            'd002': float(d002), 'd100': float(d100),
            'rho_struc': float(rho_struc), 'porosity': float(porosity),
            'rho_sample': float(rho_sample),
            'sld_struc': float(sld_struc), 'sld_sample': float(sld_sample),
            'contrast_porod': float(c_porod), 'contrast_micropore': float(c_mp),
        }

    # ── Forward model ───────────────────────────────────────────────────────

    def evaluate_components(self, q: np.ndarray,
                            material: Optional[Dict[str, float]] = None,
                            ) -> Dict[str, np.ndarray]:
        """The three contributions and their sum, on grid ``q``.

        Args:
            q: Scattering vector [Å⁻¹].
            material: Pre-resolved output of :meth:`resolve_material`, passed in
                during a fit so the chain runs once per iteration rather than
                three times.

        Returns:
            dict with ``porod``, ``mp``, ``waxs``, ``background`` and ``total``.
            A disabled section contributes an array of zeros rather than being
            absent, so callers can plot unconditionally.
        """
        q = np.asarray(q, dtype=float)
        mat = material if material is not None else self.resolve_material()
        zeros = np.zeros_like(q)

        I_porod = self._eval_porod(q, mat['contrast_porod'])
        I_mp = self._eval_micropore(q, mat['contrast_micropore'])
        I_waxs = self._eval_waxs(q)
        bg = (np.full_like(q, float(self.background.flat_background))
              if self.background.enabled else zeros)

        total = I_porod + I_mp + I_waxs + bg
        return {'porod': I_porod, 'mp': I_mp, 'waxs': I_waxs,
                'background': bg, 'total': total}

    def evaluate(self, q: np.ndarray,
                 material: Optional[Dict[str, float]] = None) -> np.ndarray:
        """Total model intensity [cm⁻¹] on grid ``q`` — eq. (1)."""
        return self.evaluate_components(q, material)['total']

    def _eval_porod(self, q: np.ndarray, contrast: float) -> np.ndarray:
        """``2π(Δρ)²[S_macro·Q⁻ⁿ + S_rough·f_rough(Q,R_rough)]`` — eq. (3).

        The flat background is *not* included here; it is added once in
        :meth:`evaluate_components` so it stays a single instrumental term
        rather than part of the grain morphology.
        """
        bg = self.background
        if not bg.enabled:
            return np.zeros_like(q)
        qsafe = np.maximum(q, 1e-30)
        bracket = bg.S_macro * qsafe ** (-float(bg.porod_exponent))
        if bg.use_roughness:
            bracket = bracket + bg.S_rough * porod_roughness_factor(q, bg.R_rough)
        out = 2.0 * np.pi * float(contrast) * _POROD_UNIT * bracket
        return np.where(np.isfinite(out), out, 0.0)

    def _eval_micropore(self, q: np.ndarray, contrast: float) -> np.ndarray:
        """``I_mp`` for whichever SAXS-region model is selected."""
        s = self.saxs
        if not s.enabled:
            return np.zeros_like(q)

        if s.mode == 'teubner_strey':
            # Exactly the Simple Fits model, A ≡ 1 (eq. 8).
            out = _teubner_strey(q, float(s.ts_I0), 1.0, float(s.ts_C1), float(s.ts_C2))
            return np.where(np.isfinite(out), out, 0.0)

        if s.mode != 'fractal':
            raise ValueError(
                f"Unknown SAXS-region mode {s.mode!r}; expected one of "
                f"{sorted(CARBON_SAXS_MODES)}")

        # I₀ = φ·(Δρ)²·V_pore — eq. (6), with V in Å³ and the 1e-4 that turns
        # 10²⁰cm⁻⁴ × Å³ into cm⁻¹.
        r = max(float(s.pore_radius), 1e-10)
        V = (4.0 / 3.0) * np.pi * r ** 3
        I0 = float(s.phi) * float(contrast) * V * _VOL_UNIT
        out = I0 * globule_form_factor(q, r, s.globule_k)
        if s.use_fractal:
            out = out * teixeira_structure_factor(q, s.fractal_D, s.fractal_sigma, r)
        return np.where(np.isfinite(out), out, 0.0)

    def _eval_waxs(self, q: np.ndarray) -> np.ndarray:
        """``I_waxs`` — summed Voigt peaks, orientation-averaged and damped.

        ``envelope · exp(−Q²⟨δz²⟩/3) · Σᵢ Vᵢ(Q) / Q²`` (eq. 2, and eq. 16–17
        when the envelope is the crumpled-layer one).  The three factors are
        kept separate and named, rather than folded into the peak shape,
        because each is independently testable and the paper writes them that
        way.
        """
        w = self.waxs
        if not w.enabled:
            return np.zeros_like(q)

        peaks = np.zeros_like(q)
        for pk in self.peaks:
            if not pk.enabled:
                continue
            peaks = peaks + voigt_peak(q, float(pk.K), float(pk.Q0),
                                       float(pk.FWHM_G), float(pk.FWHM_L))
        if not np.any(peaks):
            return np.zeros_like(q)

        # Distortions of the first kind: intensity only, never width or position.
        if w.delta_z2 > 0:
            peaks = peaks * np.exp(-(q ** 2) * float(w.delta_z2) / 3.0)

        # Powder average of a locally 1-D stacking correlation (Annex 3 A3.25).
        if w.use_orientation_factor:
            peaks = peaks / np.maximum(q, 1e-30) ** 2

        if w.envelope == 'crumpled':
            D, sigma, R = self.effective_waxs_geometry()
            peaks = peaks * teixeira_structure_factor(q, D, sigma, R) \
                * discoid_form_factor(q, R)
        elif w.envelope != 'none':
            raise ValueError(
                f"Unknown WAXS envelope {w.envelope!r}; expected one of "
                f"{sorted(CARBON_WAXS_ENVELOPES)}")

        return np.where(np.isfinite(peaks), peaks, 0.0)

    # ── Parameter table ─────────────────────────────────────────────────────

    def parameter_refs(self, active_only: bool = True) -> List[_ParamRef]:
        """Every scalar parameter of the model, in a stable order.

        Args:
            active_only: Drop parameters that the current model *shape* makes
                meaningless — a disabled section, the branch of the SAXS-region
                selector that is not chosen, roughness when it is switched off,
                envelope geometry when the envelope is flat, and any parameter
                the user has linked to another one.  Passing ``False`` returns
                the full table, which is what the GUI needs in order to draw
                greyed-out rows.

        Returns:
            list of :class:`_ParamRef`.  The ``key`` is the stable dotted name
            used in JSON, HDF5 and the report — ``background.S_macro``,
            ``saxs.pore_radius``, ``peak.002.Q0`` — and peak keys use the peak
            *label*, not its index, so inserting a peak does not rename the
            others' results.
        """
        refs: List[_ParamRef] = []

        def add(owner, attr, label, unit=''):
            prefix = {id(self.background): 'background', id(self.saxs): 'saxs',
                      id(self.waxs): 'waxs'}.get(id(owner))
            refs.append(_ParamRef(owner, attr, f'{prefix}.{attr}', label, unit))

        bg = self.background
        if bg.enabled or not active_only:
            add(bg, 'S_macro', 'Grain surface area S_macro', 'cm²/cm³')
            add(bg, 'porod_exponent', 'Porod exponent n')
            if bg.use_roughness or not active_only:
                add(bg, 'S_rough', 'Roughness surface area S_rough', 'cm²/cm³')
                add(bg, 'R_rough', 'Roughness length R_rough', 'Å')
            add(bg, 'flat_background', 'Flat background', 'cm⁻¹')

        s = self.saxs
        if s.enabled or not active_only:
            if s.mode == 'fractal' or not active_only:
                add(s, 'phi', 'Pore volume fraction φ')
                add(s, 'pore_radius', 'Pore radius r', 'Å')
                add(s, 'globule_k', 'Globule shape factor k')
                if s.use_fractal or not active_only:
                    add(s, 'fractal_D', 'Fractal dimension D')
                    add(s, 'fractal_sigma', 'Fractal cutoff Σ', 'Å')
            if s.mode == 'teubner_strey' or not active_only:
                add(s, 'ts_I0', 'Teubner-Strey I₀', 'cm⁻¹')
                add(s, 'ts_C1', 'Teubner-Strey C₁', 'Å²')
                add(s, 'ts_C2', 'Teubner-Strey C₂', 'Å⁴')

        w = self.waxs
        if w.enabled or not active_only:
            add(w, 'delta_z2', 'Stacking disorder ⟨δz²⟩', 'Å²')
            if w.envelope == 'crumpled' or not active_only:
                if not (w.link_R_to_pore and active_only):
                    add(w, 'R_layer', 'Layer transition radius R', 'Å')
                if not (w.link_D_to_saxs and active_only):
                    add(w, 'fractal_D', 'Crumpling fractal dimension D')
                if not (w.link_sigma_to_saxs and active_only):
                    add(w, 'fractal_sigma', 'Crumpling cutoff Σ', 'Å')

            for pk in self.peaks:
                if not (pk.enabled or not active_only):
                    continue
                for attr, label, unit in (
                    ('Q0', 'centre Q₀', 'Å⁻¹'), ('K', 'amplitude K', ''),
                    ('FWHM_G', 'Gaussian FWHM', 'Å⁻¹'),
                    ('FWHM_L', 'Lorentzian FWHM', 'Å⁻¹'),
                ):
                    refs.append(_ParamRef(
                        pk, attr, f'peak.{pk.label}.{attr}',
                        f'({pk.label}) {label}', unit))
        return refs

    def fitted_refs(self) -> List[_ParamRef]:
        """Active parameters whose Fit? box is ticked — the fit vector."""
        return [r for r in self.parameter_refs(active_only=True) if r.fit]

    def parameter_values(self) -> Dict[str, float]:
        """``{key: value}`` for every active parameter."""
        return {r.key: r.value for r in self.parameter_refs(active_only=True)}

    # ── Derived quantities ──────────────────────────────────────────────────

    def compute_derived(self) -> Dict[str, float]:
        """The materials-science layer: ~20 numbers the fit parameters imply.

        These are what the user actually wants — a BET-comparable specific
        surface area, a pore width, a stack height — rather than the
        coefficients that produced them.  Everything is a plain float in the
        units listed in the module docstring, so the same dict goes to the
        panel, the report, the HDF5 ``derived/`` group and the api layer
        unchanged.

        Non-physical or not-applicable entries are ``nan`` rather than missing,
        so the set of keys depends only on the model's *mode*, not on whether a
        particular fit happened to succeed.
        """
        out: Dict[str, float] = {}
        mat = self.resolve_material()
        out.update(mat)

        nan = float('nan')
        rho_sample = mat['rho_sample']

        # ── Background: grain surface areas ──
        bg = self.background
        if bg.enabled and abs(bg.porod_exponent - 4.0) < 1e-9:
            s_macro = float(bg.S_macro)
            s_rough = float(bg.S_rough) if bg.use_roughness else 0.0
            out['S_macro'] = s_macro
            out['S_rough'] = s_rough
            out['S_part'] = s_macro + s_rough
            out['S_macro_m2_g'] = _cd.specific_surface_area_m2_g(s_macro, rho_sample)
            out['S_rough_m2_g'] = _cd.specific_surface_area_m2_g(s_rough, rho_sample)
            out['S_part_m2_g'] = _cd.specific_surface_area_m2_g(
                s_macro + s_rough, rho_sample)
        else:
            # A free Porod exponent makes the prefactor stop being an area.
            for k in ('S_macro', 'S_rough', 'S_part',
                      'S_macro_m2_g', 'S_rough_m2_g', 'S_part_m2_g'):
                out[k] = nan

        # ── SAXS region: porosity, pore size, micropore surface area ──
        s = self.saxs
        for k in ('mp_phi', 'mp_I0', 'mp_radius', 'S_mp', 'S_mp_m2_g',
                  'ts_xi', 'ts_d', 'ts_fa', 'w_pore', 'w_carbon',
                  'ts_r_spheroid'):
            out[k] = nan

        if s.enabled and s.mode == 'fractal':
            r = max(float(s.pore_radius), 1e-30)
            phi = float(s.phi)
            k = float(s.globule_k)
            out['mp_phi'] = phi
            out['mp_radius'] = r
            out['mp_I0'] = (phi * mat['contrast_micropore']
                            * (4.0 / 3.0) * np.pi * r ** 3 * _VOL_UNIT)
            # Porod limit of I_mp is I₀·(9/2)k/(Qr)⁴, which against
            # 2π(Δρ)²S/Q⁴ gives S = 3φk/r — the sphere result for k = 1.
            # ×1e8 converts Å⁻¹ to cm⁻¹.
            out['S_mp'] = 3.0 * phi * k / r * 1.0e8
            out['S_mp_m2_g'] = _cd.specific_surface_area_m2_g(out['S_mp'], rho_sample)

        elif s.enabled and s.mode == 'teubner_strey':
            C1, C2, I0 = float(s.ts_C1), float(s.ts_C2), float(s.ts_I0)
            if C2 > 0:
                half = 0.5 / np.sqrt(C2)
                xi_arg = half + C1 / (4.0 * C2)
                d_arg = half - C1 / (4.0 * C2)
                xi = 1.0 / np.sqrt(xi_arg) if xi_arg > 0 else nan
                d = 2.0 * np.pi / np.sqrt(d_arg) if d_arg > 0 else nan
                fa = C1 / (2.0 * np.sqrt(C2))
                out['ts_xi'] = float(xi)
                out['ts_d'] = float(d)
                out['ts_fa'] = float(fa)
                # φ(1−φ) from I₀ = 8π·φ(1−φ)(Δρ)²·ξ³/(1+(2πξ/d)²)² (eq. 8).
                c_mp = mat['contrast_micropore']
                if np.isfinite(xi) and np.isfinite(d) and c_mp > 0 and xi > 0:
                    denom = 8.0 * np.pi * c_mp * xi ** 3 * _VOL_UNIT
                    x = I0 * (1.0 + (2.0 * np.pi * xi / d) ** 2) ** 2 / denom
                    # Take the low-porosity root, as the Modeling tool does for
                    # its own Vf ≈ scale inversion.
                    phi = 0.5 * (1.0 - np.sqrt(1.0 - 4.0 * x)) if x <= 0.25 else nan
                    out['mp_phi'] = float(phi)
                    if np.isfinite(phi) and 0.0 < phi < 1.0:
                        # Corrigendum (2020): general for any f_a, unlike the
                        # original main-text r = sqrt(5·C1).
                        out['w_pore'] = float(xi / (1.0 - phi))
                        out['w_carbon'] = float(xi / phi)
                # Porod limit I₀/(C₂Q⁴) against 2π(Δρ)²S·1e-12/Q⁴ (eq. 15).
                if C2 > 0 and c_mp > 0:
                    out['S_mp'] = float(I0 / (C2 * 2.0 * np.pi * c_mp * _POROD_UNIT))
                    out['S_mp_m2_g'] = _cd.specific_surface_area_m2_g(
                        out['S_mp'], rho_sample)
                # The pre-corrigendum spheroid estimate, only where it is valid.
                if C1 < 0 and abs(fa - 0.4) < 0.15:
                    out['ts_r_spheroid'] = float(np.sqrt(5.0 * abs(C1)))

        # ── WAXS region ──
        w = self.waxs
        out['delta_z2'] = float(w.delta_z2) if w.enabled else nan
        if w.enabled and w.envelope == 'crumpled':
            D, sigma, R = self.effective_waxs_geometry()
            out['crumple_D'], out['crumple_sigma'], out['crumple_R'] = D, sigma, R
        else:
            out['crumple_D'] = out['crumple_sigma'] = out['crumple_R'] = nan

        for pk in self.peaks:
            if not pk.enabled:
                continue
            tag = f'peak_{pk.label}'
            d = _cd.d_spacing(pk.Q0)
            fw = voigt_fwhm(pk.FWHM_G, pk.FWHM_L)
            out[f'{tag}_Q0'] = float(pk.Q0)
            out[f'{tag}_d'] = float(d)
            out[f'{tag}_FWHM'] = float(fw)
            out[f'{tag}_L'] = _coherence_length(pk.FWHM_G, pk.FWHM_L)
            out[f'{tag}_height'] = float(
                pk.K * np.exp(-(pk.Q0 ** 2) * w.delta_z2 / 3.0)
                / (pk.Q0 ** 2 if w.use_orientation_factor else 1.0))

        pk002 = self.peak_by_label(self.material.d002_label)
        if pk002 is not None:
            Lc = _coherence_length(pk002.FWHM_G, pk002.FWHM_L)
            d002 = _cd.d_spacing(pk002.Q0)
            out['L_c'] = float(Lc)
            out['N_layers'] = float(Lc / d002) if (d002 and d002 > 0) else nan
        else:
            out['L_c'] = out['N_layers'] = nan
        pk100 = self.peak_by_label(self.material.d100_label)
        out['L_a'] = (_coherence_length(pk100.FWHM_G, pk100.FWHM_L,
                                        shape_factor=1.84)
                       if pk100 is not None else nan)

        return {k: float(v) for k, v in out.items()}

    # ── Serialisation ───────────────────────────────────────────────────────

    def to_dict(self) -> dict:
        """Settings — not results — as a plain dict.

        This shape is the contract the io, batch, api and GUI layers all build
        on: ``dict → JSON → HDF5`` unchanged.  No data arrays, no χ², nothing
        that belongs on :class:`CarbonFitResult`.
        """
        return {
            'schema_version': self.SCHEMA_VERSION,
            'background': self.background.to_dict(),
            'saxs': self.saxs.to_dict(),
            'waxs': self.waxs.to_dict(),
            'peaks': [pk.to_dict() for pk in self.peaks],
            'material': self.material.to_dict(),
            'q_min': float(self.q_min),
            'q_max': float(self.q_max),
            'weighting': str(self.weighting),
            'max_iterations': int(self.max_iterations),
            'n_mc_runs': int(self.n_mc_runs),
        }

    @classmethod
    def from_dict(cls, d: Optional[dict]) -> 'CarbonFitModel':
        """Rebuild from :meth:`to_dict`; every missing field takes its default.

        A file written before a field existed has to open with that field at its
        default — that is the backwards-compatibility guarantee, not politeness.
        An absent ``peaks`` list restores the default (002)/(100)/(004) set,
        while an explicitly empty one stays empty, because "no peaks" is a
        legitimate saved state (a pure SAXS fit).
        """
        obj = cls()
        if not d:
            return obj
        obj.background = CarbonBackground.from_dict(d.get('background') or {})
        obj.saxs = CarbonSaxsRegion.from_dict(d.get('saxs') or {})
        obj.waxs = CarbonWaxsRegion.from_dict(d.get('waxs') or {})
        obj.material = CarbonMaterial.from_dict(d.get('material') or {})
        if 'peaks' in d and d['peaks'] is not None:
            obj.peaks = [CarbonWaxsPeak.from_dict(p or {}) for p in d['peaks']]
        for key, cast in (('q_min', float), ('q_max', float),
                          ('weighting', str), ('max_iterations', int),
                          ('n_mc_runs', int)):
            if d.get(key) is not None:
                try:
                    setattr(obj, key, cast(d[key]))
                except (TypeError, ValueError):
                    log.debug("carbon_fit: ignoring bad %s=%r", key, d[key])
        if obj.saxs.mode not in CARBON_SAXS_MODES:
            log.warning("carbon_fit: unknown SAXS mode %r, using 'fractal'",
                        obj.saxs.mode)
            obj.saxs.mode = 'fractal'
        if obj.waxs.envelope not in CARBON_WAXS_ENVELOPES:
            log.warning("carbon_fit: unknown WAXS envelope %r, using 'none'",
                        obj.waxs.envelope)
            obj.waxs.envelope = 'none'
        return obj

    def copy(self) -> 'CarbonFitModel':
        """Deep copy, via the serialisation round-trip."""
        return CarbonFitModel.from_dict(self.to_dict())

    # ── Fitting ─────────────────────────────────────────────────────────────

    def q_mask(self, q: np.ndarray) -> np.ndarray:
        """Boolean mask for the model's fit range, ``0`` meaning "no limit"."""
        q = np.asarray(q, dtype=float)
        mask = np.isfinite(q) & (q > 0)
        if self.q_min > 0:
            mask &= q >= self.q_min
        if self.q_max > 0:
            mask &= q <= self.q_max
        return mask

    def _sigma(self, I: np.ndarray, error: Optional[np.ndarray]) -> np.ndarray:
        """Residual weights for the selected weighting mode.

        ``'auto'`` uses the measured uncertainties when they are present and
        usable, and relative weighting otherwise.  Relative weighting is the
        right default for this tool even when errors exist but are optimistic:
        the fit spans five decades in Q and ten or more in intensity, and
        absolute weighting would let the low-Q Porod region set every
        parameter while the diffraction peaks contributed nothing.
        """
        I = np.asarray(I, dtype=float)
        mode = self.weighting
        if mode == 'auto':
            ok = (error is not None and np.all(np.isfinite(error))
                  and np.all(np.asarray(error) > 0))
            mode = 'sigma' if ok else 'relative'
        if mode == 'sigma' and error is not None:
            sig = np.abs(np.asarray(error, dtype=float))
            sig = np.where(np.isfinite(sig) & (sig > 0), sig, np.abs(I))
        else:
            sig = np.abs(I)
        floor = np.nanmax(np.abs(I)) * 1e-12 if np.any(np.isfinite(I)) else 1e-30
        return np.maximum(sig, max(floor, 1e-30))

    def fit(
        self,
        q: np.ndarray,
        I: np.ndarray,
        error: Optional[np.ndarray] = None,
        progress: Optional[object] = None,
    ) -> CarbonFitResult:
        """Least-squares fit of the whole model to the whole Q range.

        All three contributions are refined together.  That is not a
        convenience: the contrast that scales the background and the micropore
        term is computed from the WAXS peak positions and the porosity, so the
        regions are coupled through the physics and fitting them in sequence
        would converge to a different — wrong — answer.

        Args:
            q:        Scattering vector [Å⁻¹].
            I:        Measured intensity [cm⁻¹].
            error:    1-σ uncertainties, optional.
            progress: Optional callable ``progress(iteration, chi2)`` for the
                GUI.  Raising :class:`CarbonFitAborted` from it stops the fit
                and restores the starting parameters — that is how the Stop
                button works.

        Returns:
            :class:`CarbonFitResult`.  The model itself is left holding the
            best-fit parameters, so a second call refines from there.

        Raises:
            ValueError: if no parameter is marked for fitting, or the fit range
                leaves fewer points than parameters.
        """
        q = np.asarray(q, dtype=float)
        I = np.asarray(I, dtype=float)
        mask = self.q_mask(q) & np.isfinite(I)
        qf, If = q[mask], I[mask]
        ef = np.asarray(error, dtype=float)[mask] if error is not None else None

        refs = self.fitted_refs()
        if not refs:
            raise ValueError(
                "No parameters are marked for fitting — tick at least one Fit? box.")
        if qf.size <= len(refs):
            raise ValueError(
                f"Fit range holds {qf.size} points but the model has "
                f"{len(refs)} free parameters.")

        sigma = self._sigma(If, ef)
        x_start = np.array([r.value for r in refs], dtype=float)
        lo = np.array([r.limits[0] for r in refs], dtype=float)
        hi = np.array([r.limits[1] for r in refs], dtype=float)
        # A start value sitting exactly on a bound stops TRF dead, so the
        # optimiser is handed a nudged copy — but an abort restores the values
        # the user actually typed, not the nudged ones.
        x0 = np.clip(x_start, lo + 1e-12 * np.maximum(np.abs(lo), 1.0),
                     hi - 1e-12 * np.maximum(np.abs(hi), 1.0))

        state = {'n': 0}

        def residuals(x: np.ndarray) -> np.ndarray:
            for ref, v in zip(refs, x):
                ref.value = v
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                model = self.evaluate(qf)
                res = (If - model) / sigma
            res = np.where(np.isfinite(res), res, 1.0e6)
            state['n'] += 1
            if progress is not None:
                progress(state['n'], float(np.sum(res ** 2)))
            return res

        try:
            sol = least_squares(
                residuals, x0, bounds=(lo, hi), method='trf',
                max_nfev=int(self.max_iterations) * max(len(refs), 1),
                xtol=1e-10, ftol=1e-10, gtol=1e-10)
        except Exception as exc:                # CarbonFitAborted, or a genuine
            for ref, v in zip(refs, x_start):   # numerical failure inside scipy
                ref.value = v
            return self._result(qf, If, ef, success=False,
                                message=f"Fit aborted: {exc}")

        for ref, v in zip(refs, sol.x):
            ref.value = v

        errors = self._parameter_errors(sol, refs, qf.size)
        if self.n_mc_runs and self.n_mc_runs > 1:
            mc = self._monte_carlo(qf, If, sigma, refs, sol.x)
            errors.update(mc)

        res = self._result(qf, If, ef, success=bool(sol.success),
                           message=str(sol.message), errors=errors,
                           n_iterations=int(state['n']))
        return res

    def _parameter_errors(self, sol, refs: List[_ParamRef],
                          n_points: int) -> Dict[str, float]:
        """1-σ uncertainties from the Gauss-Newton covariance ``(JᵀJ)⁻¹·s²``.

        Rank-deficient Jacobians are common here — a linked or badly-determined
        parameter makes a column vanish — so a singular decomposition yields
        ``nan`` for that parameter rather than an exception or a fabricated
        number.
        """
        out = {r.key: float('nan') for r in refs}
        try:
            J = np.asarray(sol.jac, dtype=float)
            dof = max(n_points - len(refs), 1)
            s2 = float(np.sum(np.asarray(sol.fun) ** 2)) / dof
            _, sv, VT = np.linalg.svd(J, full_matrices=False)
            keep = sv > max(sv[0], 1e-300) * 1e-12 if sv.size else np.array([], bool)
            if not np.any(keep):
                return out
            inv = (VT[keep].T / sv[keep] ** 2) @ VT[keep]
            var = np.diag(inv) * s2
            for r, v in zip(refs, var):
                out[r.key] = float(np.sqrt(v)) if v >= 0 else float('nan')
        except Exception:
            log.debug("carbon_fit: covariance estimate failed", exc_info=True)
        return out

    def _monte_carlo(self, q, I, sigma, refs, x_best) -> Dict[str, float]:
        """Uncertainties by refitting noise-perturbed copies of the data.

        More honest than the covariance estimate when parameters are
        correlated, which in this model they always are (contrast couples the
        regions).  Serial: the fit is fast enough that process startup would
        dominate, unlike Modeling's much heavier per-iteration cost.
        """
        rng = np.random.default_rng(12345)
        samples: List[np.ndarray] = []
        saved = [r.value for r in refs]
        lo = np.array([r.limits[0] for r in refs], dtype=float)
        hi = np.array([r.limits[1] for r in refs], dtype=float)
        for _ in range(int(self.n_mc_runs)):
            I_pert = I + rng.normal(0.0, sigma)

            def residuals(x, _Ip=I_pert):
                for ref, v in zip(refs, x):
                    ref.value = v
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    res = (_Ip - self.evaluate(q)) / sigma
                return np.where(np.isfinite(res), res, 1.0e6)

            try:
                sol = least_squares(residuals, x_best, bounds=(lo, hi),
                                    method='trf', max_nfev=100 * len(refs))
                samples.append(np.asarray(sol.x, dtype=float))
            except Exception:
                log.debug("carbon_fit: one Monte-Carlo pass failed", exc_info=True)
        for ref, v in zip(refs, saved):
            ref.value = v
        if len(samples) < 2:
            return {}
        arr = np.vstack(samples)
        return {r.key: float(np.std(arr[:, i], ddof=1))
                for i, r in enumerate(refs)}

    def _result(self, q, I, error, *, success: bool, message: str,
                errors: Optional[Dict[str, float]] = None,
                n_iterations: int = 0) -> CarbonFitResult:
        """Assemble a :class:`CarbonFitResult` from the model's current state."""
        mat = self.resolve_material()
        comps = self.evaluate_components(q, mat)
        model = comps['total']
        sigma = self._sigma(I, error)
        resid = (I - model) / sigma
        n_par = len(self.fitted_refs())
        chi2 = float(np.sum(resid ** 2))

        quality: Dict[str, object] = {}
        try:
            from pyirena.core.fit_metrics import fit_quality_metrics
            quality = fit_quality_metrics(q, I, model, error, n_par)
        except Exception:
            log.debug("carbon_fit: quality metrics unavailable", exc_info=True)

        return CarbonFitResult(
            success=success, message=message,
            params=self.parameter_values(),
            errors=errors or {},
            derived=self.compute_derived(),
            chi_squared=chi2,
            reduced_chi_squared=chi2 / max(q.size - n_par, 1),
            n_points=int(q.size), n_params=n_par, n_iterations=n_iterations,
            timestamp=datetime.now().isoformat(timespec='seconds'),
            quality=quality,
            q=q, I_data=I, I_error=error, I_model=model,
            I_porod=comps['porod'] + comps['background'],
            I_mp=comps['mp'], I_waxs=comps['waxs'], residuals=resid,
        )
