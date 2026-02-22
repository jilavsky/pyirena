"""
Simple Fits model library for SAS data.

Provides single-model fits ported from the Irena Igor Pro package
(IR3_SimpleFits.ipf and IR3_SystemSpecificModels.ipf).

Models
------
From SimpleFits:
  * Guinier          — I0·exp(−Q²Rg²/3)
  * Guinier Rod      — I0·exp(−Q²Rc²/2)/Q
  * Guinier Sheet    — I0·exp(−Q²Rg²)/Q²
  * Porod            — Kp·Q⁻⁴ + Background
  * Power Law        — Prefactor·Q⁻Exponent + Background
  * Sphere           — Scale·|FF(QR)|² + Background
  * Spheroid         — orientationally-averaged sphere FF with aspect ratio

From SystemSpecificModels:
  * Debye-Bueche         — Prefactor·Eta²·ξ³/(1+Q²ξ²)²
  * Treubner-Strey        — Prefactor/(A+C1·Q²+C2·Q⁴)
  * Benedetti-Ciccariello — coated-interface Porod law
  * Hermans               — lamellar structure (complex exponential form)
  * Hybrid Hermans        — Hermans + two Unified-level contributions
  * Unified Born Green    — two-level Unified with structure factor

Usage example
-------------
>>> from pyirena.core.simple_fits import SimpleFitModel
>>> import numpy as np
>>> model = SimpleFitModel()
>>> model.model = 'Guinier'
>>> model.params['I0'] = 10.0
>>> model.params['Rg'] = 80.0
>>> result = model.fit(q, intensity, error)
>>> print(result['chi2'], result['params'])
"""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import curve_fit
from scipy.special import erf

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Gauss-Legendre quadrature nodes on [0, 1] — used for spheroid orientation
# averaging (precomputed at module load, shared across all calls).
# ---------------------------------------------------------------------------
_GL_N = 50
_GL_X_11, _GL_W_11 = leggauss(_GL_N)        # nodes on [-1, 1]
_GL_X = (_GL_X_11 + 1.0) / 2.0              # remapped to [0, 1]
_GL_W = _GL_W_11 / 2.0                      # scaled weights

# ---------------------------------------------------------------------------
# UBG angular quadrature (matches Igor 20-point Riemann sum)
# ---------------------------------------------------------------------------
_UBG_N_ANGLES = 20
_UBG_ANGLES = np.arange(_UBG_N_ANGLES) * (np.pi / (2 * _UBG_N_ANGLES))
_UBG_COS_A = np.cos(_UBG_ANGLES)            # shape (20,)

# ---------------------------------------------------------------------------
# Complex background parameter spec
# ---------------------------------------------------------------------------
# (name, default, lower_bound, upper_bound)
_BG_PARAMS: list[tuple] = [
    ('BG_A',    0.0,  0.0,  None),
    ('BG_n',    4.0,  0.1,  10.0),
    ('BG_flat', 0.0,  None, None),
]


# ===========================================================================
# Private model formula functions
# ===========================================================================

def _guinier(q: np.ndarray, I0: float, Rg: float) -> np.ndarray:
    """I(Q) = I0·exp(−Q²·Rg²/3)"""
    return I0 * np.exp(-q**2 * Rg**2 / 3.0)


def _guinier_rod(q: np.ndarray, I0: float, Rc: float) -> np.ndarray:
    """I(Q) = I0·exp(−Q²·Rc²/2) / Q"""
    with np.errstate(divide='ignore', invalid='ignore'):
        val = I0 * np.exp(-q**2 * Rc**2 / 2.0) / q
    return np.where(q > 0, val, 0.0)


def _guinier_sheet(q: np.ndarray, I0: float, Rg: float) -> np.ndarray:
    """I(Q) = I0·exp(−Q²·Rg²) / Q²"""
    with np.errstate(divide='ignore', invalid='ignore'):
        val = I0 * np.exp(-q**2 * Rg**2) / q**2
    return np.where(q > 0, val, 0.0)


def _porod(q: np.ndarray, Kp: float, Background: float) -> np.ndarray:
    """I(Q) = Kp·Q⁻⁴ + Background"""
    with np.errstate(divide='ignore', invalid='ignore'):
        val = Kp / q**4
    return np.where(q > 0, val, 0.0) + Background


def _power_law(q: np.ndarray, Prefactor: float, Exponent: float,
               Background: float) -> np.ndarray:
    """I(Q) = Prefactor·Q⁻Exponent + Background"""
    with np.errstate(divide='ignore', invalid='ignore'):
        val = Prefactor * q**(-Exponent)
    return np.where(q > 0, val, 0.0) + Background


def _sphere_ff_sq(u: np.ndarray) -> np.ndarray:
    """Sphere form factor squared: [3(sin u − u cos u)/u³]²"""
    with np.errstate(divide='ignore', invalid='ignore'):
        ff = np.where(
            np.abs(u) < 1e-6,
            1.0,
            3.0 * (np.sin(u) - u * np.cos(u)) / u**3,
        )
    return ff**2


def _sphere(q: np.ndarray, Scale: float, R: float,
            Background: float) -> np.ndarray:
    """I(Q) = Scale·|FF(QR)|² + Background"""
    return Scale * _sphere_ff_sq(q * R) + Background


def _spheroid(q: np.ndarray, Scale: float, R: float, Beta: float,
              Background: float) -> np.ndarray:
    """Orientationally-averaged spheroid form factor.

    Uses Gauss-Legendre quadrature over cos(θ) from 0 to 1.
    The semi-axes are (R, R, Beta·R); Beta=1 gives a sphere.
    """
    # x = cos(θ), r(x) = R·sqrt(1 + (Beta²−1)·x²)
    x = _GL_X                                     # shape (N,)
    w = _GL_W
    r = R * np.sqrt(1.0 + (Beta**2 - 1.0) * x**2)  # (N,)

    # ff_sq[i, j] = sphere FF² at q[i] * r[j]
    u = q[:, None] * r[None, :]                   # (n_q, N)
    ff_sq = _sphere_ff_sq(u)                       # (n_q, N)
    avg = np.einsum('ij,j->i', ff_sq, w)          # weighted sum over j

    return Scale * avg + Background


def _debye_bueche(q: np.ndarray, Prefactor: float, Eta: float,
                  CorrLength: float) -> np.ndarray:
    """I(Q) = Prefactor·Eta²·ξ³/(1+Q²ξ²)²

    The Prefactor absorbs the wavelength-dependent constant from the
    original Debye-Bueche formula: 32π³·λ⁻⁴ (Igor-compatible form).
    """
    return Prefactor * Eta**2 * CorrLength**3 / (1.0 + q**2 * CorrLength**2)**2


def _treubner_strey(q: np.ndarray, Prefactor: float, A: float,
                    C1: float, C2: float) -> np.ndarray:
    """I(Q) = Prefactor/(A + C1·Q² + C2·Q⁴)"""
    denom = A + C1 * q**2 + C2 * q**4
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(denom != 0.0, Prefactor / denom, 0.0)


def _benedetti_ciccariello(q: np.ndarray, SolidSLD: float, VoidSLD: float,
                            LayerSLD: float, Sp: float, t: float) -> np.ndarray:
    """Benedetti-Ciccariello model for coated interfaces.

    Parameters
    ----------
    SolidSLD, VoidSLD, LayerSLD : float
        Scattering length densities [10¹⁰ cm⁻²].
    Sp : float
        Porod's specific surface area [cm⁻¹].
    t : float
        Coating layer thickness [Å].
    """
    # Contrast (convert SLD from 10^10 cm^-2 → cm^-2)
    n12 = (SolidSLD - VoidSLD) * 1e10
    denom = SolidSLD - VoidSLD
    if abs(denom) < 1e-30:
        NuValue = 0.0
    else:
        NuValue = (SolidSLD - 2.0 * LayerSLD + VoidSLD) / denom
    Alpha = (1.0 + NuValue**2) / 2.0
    Rnu = (1.0 - NuValue**2) / (1.0 + NuValue**2)
    # Q in Å⁻¹; Q^4 in Å⁻⁴; factor 1e32 converts to cm⁻⁴ for unit consistency
    with np.errstate(divide='ignore', invalid='ignore'):
        prefact = 2.0 * np.pi * n12**2 * Alpha * Sp / (q**4 * 1e32)
    return np.where(q > 0, prefact * (1.0 + Rnu * np.cos(q * t)), 0.0)


def _hermans_core(q: np.ndarray, B: float, d1: float, sigma1: float,
                  d2: float, sigma2: float) -> np.ndarray:
    """Hermans lamellar model (complex-plane form).

    I(Q) = B/(2π)⁴·s⁴ · Re[(1−H1)(1−H2)/(1−H1·H2)]
    where s = Q/(2π), H1 and H2 are complex exponentials.

    Reference: https://doi.org/10.1016/j.polymer.2021.124281
    """
    s = q / (2.0 * np.pi)
    H1 = np.exp(2j * np.pi * d1 * s - 2.0 * np.pi**2 * sigma1**2 * s**2)
    H2 = np.exp(2j * np.pi * d2 * s - 2.0 * np.pi**2 * sigma2**2 * s**2)
    Bs = B / (2.0 * np.pi)**4
    denom = 1.0 - H1 * H2
    # Guard against near-zero denominator
    denom = np.where(np.abs(denom) < 1e-30, 1e-30 + 0j, denom)
    with np.errstate(divide='ignore', invalid='ignore'):
        result = (Bs / s**4) * np.real((1.0 - H1) * (1.0 - H2) / denom)
    return np.maximum(result, 0.0)


def _hermans(q: np.ndarray, B: float, AmorphThick: float, AmorphSigma: float,
             LamThick: float, LamSigma: float) -> np.ndarray:
    """Hermans model.  See _hermans_core for formula."""
    return _hermans_core(q, B, AmorphThick, AmorphSigma, LamThick, LamSigma)


def _qstar(q: np.ndarray, Rg: float) -> np.ndarray:
    """Unified Q* = Q / erf(Q·|Rg|/√6)³ (Beaucage cut-off function)."""
    arg = q * abs(Rg) / np.sqrt(6.0)
    erf_val = erf(arg)
    with np.errstate(divide='ignore', invalid='ignore'):
        qs = np.where(erf_val > 1e-30, q / erf_val**3, q)
    return qs


def _hybrid_hermans(q: np.ndarray,
                    B: float, AmorphThick: float, AmorphSigma: float,
                    LamThick: float, LamSigma: float,
                    G2: float, Rg2: float,
                    G3: float, Rg3: float, B3: float, P3: float) -> np.ndarray:
    """Hybrid Hermans model.

    Combines the Hermans lamellar model with two additional Unified levels:
      I = I_Hermans + G2·exp(−Q²Rg2²/3) + G3·exp(−Q²Rg3²/3) + B3·Qstar2⁻P3
    where Qstar2 uses Rg2 as the cut-off length.
    """
    I_herm = _hermans_core(q, B, AmorphThick, AmorphSigma, LamThick, LamSigma)
    qs2 = _qstar(q, Rg2)
    level2 = G2 * np.exp(-q**2 * Rg2**2 / 3.0)
    level3 = (G3 * np.exp(-q**2 * Rg3**2 / 3.0)
              + B3 * np.exp(-q**2 * Rg2**2 / 3.0) * qs2**(-abs(P3)))
    return I_herm + level2 + level3


def _ubg(q: np.ndarray, Rg1: float, B1: float, Pack: float,
         CorrDist: float, StackIrreg: float, kI: float) -> np.ndarray:
    """Unified Born Green model (Formula 8).

    Two-level Unified low-Q fit with a Born-Green structure factor.
    Parameters:  Rg1, B1, Pack (packing factor), CorrDist (correlation
    distance ξ), StackIrreg (stacking irregularity δ), kI (damping).

    Reference: https://doi.org/10.1016/j.polymer.2021.124281, Formula 8
    """
    if Rg1 <= 0 or B1 <= 0 or CorrDist <= 0:
        return np.zeros_like(q)

    # Derived structural parameters
    Rad = np.sqrt(abs(Pack * CorrDist * Rg1) / 4.0)
    denom12 = 2.0 * Rg1 + Rad
    G1 = B1 * Rg1**4 * Rad / denom12 if denom12 > 0 else 0.0
    B2 = 2.0 * G1 / Rg1**2
    Rg2 = np.sqrt(Rad**2 / 2.0 + Rg1**2 / 3.0)
    G2 = G1 * (Rad / Rg1)**2

    qs1 = _qstar(q, Rg1)
    qs2 = _qstar(q, Rg2)

    with np.errstate(divide='ignore', invalid='ignore'):
        base = (G1 * np.exp(-q**2 * Rg1**2 / 3.0) + B1 * qs1**(-4.0)
                + G2 * np.exp(-q**2 * Rg2**2 / 3.0)
                + np.exp(-q**2 * Rg1**2 / 3.0) * B2 * qs2**(-2.0))

    if Pack == 0.0:
        return np.maximum(base, 0.0)

    # Structure factor (vectorised over q and angles, matching Igor 20-pt sum)
    # Apply stacking irregularity: Q' = Q·exp(δ·(Q − 2π/ξ)/Q)
    two_pi_over_xi = 2.0 * np.pi / CorrDist
    with np.errstate(divide='ignore', invalid='ignore'):
        q_corr = np.where(
            q > 1e-30,
            q * np.exp(StackIrreg * (q - two_pi_over_xi) / q),
            q,
        )
    corrnum = q_corr * CorrDist                   # shape (n_q,)

    # (n_q, 20) sinc matrix
    args = corrnum[:, None] * _UBG_COS_A[None, :]  # (n_q, 20)
    sinc = np.where(np.abs(args) < 1e-10, 1.0, np.sin(args) / args)

    # Average with π/2 factor matching Igor: sumtheta += sinc * π/2; sumtheta/N
    avg = np.sum(sinc * (np.pi / 2.0), axis=1) / _UBG_N_ANGLES   # (n_q,)
    S = 1.0 + Pack * avg * np.exp(-corrnum**2 * kI)

    with np.errstate(divide='ignore', invalid='ignore'):
        result = np.where(S != 0.0, base / S, base)

    return np.maximum(result, 0.0)


# ===========================================================================
# MODEL_REGISTRY
# ===========================================================================
# Each entry:
#   'params'        : list of (name, default, lower_bound, upper_bound)
#                     lower/upper_bound: None → −inf / +inf
#   'formula'       : callable (q, *params) → np.ndarray
#   'linearization' : str key or None
#                     'guinier' | 'guinier_rod' | 'guinier_sheet' | 'porod'
#   'complex_bg'    : True if a complex background (A·Q^-n + flat) can be
#                     added; False if the formula already includes Background

MODEL_REGISTRY: dict[str, dict] = {
    'Guinier': {
        'params': [
            ('I0', 1.0,   1e-30, None),
            ('Rg', 50.0,  0.1,   10_000.0),
        ],
        'formula': _guinier,
        'linearization': 'guinier',
        'complex_bg': True,
    },
    'Guinier Rod': {
        'params': [
            ('I0', 1.0,  1e-30, None),
            ('Rc', 10.0, 0.1,   10_000.0),
        ],
        'formula': _guinier_rod,
        'linearization': 'guinier_rod',
        'complex_bg': True,
    },
    'Guinier Sheet': {
        'params': [
            ('I0', 1.0,  1e-30, None),
            ('Rg', 10.0, 0.1,   10_000.0),
        ],
        'formula': _guinier_sheet,
        'linearization': 'guinier_sheet',
        'complex_bg': True,
    },
    'Porod': {
        'params': [
            ('Kp',         1.0, 1e-30, None),
            ('Background', 0.0, None,  None),
        ],
        'formula': _porod,
        'linearization': 'porod',
        'complex_bg': False,
    },
    'Power Law': {
        'params': [
            ('Prefactor',  1.0, 1e-30, None),
            ('Exponent',   4.0, 1.0,   6.0),
            ('Background', 0.0, None,  None),
        ],
        'formula': _power_law,
        'linearization': None,
        'complex_bg': False,
    },
    'Sphere': {
        'params': [
            ('Scale',      1.0,  1e-30, None),
            ('R',          50.0, 0.1,   10_000.0),
            ('Background', 0.0,  None,  None),
        ],
        'formula': _sphere,
        'linearization': None,
        'complex_bg': True,
    },
    'Spheroid': {
        'params': [
            ('Scale',      1.0,  1e-30,  None),
            ('R',          50.0, 0.1,    10_000.0),
            ('Beta',       1.0,  0.001,  1000.0),
            ('Background', 0.0,  None,   None),
        ],
        'formula': _spheroid,
        'linearization': None,
        'complex_bg': True,
    },
    'Debye-Bueche': {
        'params': [
            ('Prefactor',  1.0,   1e-30, None),
            ('Eta',        0.01,  1e-30, None),
            ('CorrLength', 100.0, 0.1,   100_000.0),
        ],
        'formula': _debye_bueche,
        'linearization': None,
        'complex_bg': True,
    },
    'Treubner-Strey': {
        'params': [
            ('Prefactor', 1.0,    1e-30, None),
            ('A',         0.1,    1e-30, None),
            ('C1',        -30.0,  None,  None),
            ('C2',        5000.0, 1e-30, None),
        ],
        'formula': _treubner_strey,
        'linearization': None,
        'complex_bg': True,
    },
    'Benedetti-Ciccariello': {
        'params': [
            ('SolidSLD', 1.0,  None,  None),
            ('VoidSLD',  0.0,  None,  None),
            ('LayerSLD', 0.5,  None,  None),
            ('Sp',       1e5,  1e-30, None),
            ('t',        10.0, 0.1,   10_000.0),
        ],
        'formula': _benedetti_ciccariello,
        'linearization': None,
        'complex_bg': True,
    },
    'Hermans': {
        'params': [
            ('B',           1e-4, 1e-30, None),
            ('AmorphThick', 47.0, 0.1,   10_000.0),
            ('AmorphSigma', 23.0, 0.0,   10_000.0),
            ('LamThick',   146.0, 0.1,   10_000.0),
            ('LamSigma',    58.0, 0.0,   10_000.0),
        ],
        'formula': _hermans,
        'linearization': None,
        'complex_bg': True,
    },
    'Hybrid Hermans': {
        'params': [
            ('B',           1e-4,  1e-30, None),
            ('AmorphThick', 47.0,  0.1,   10_000.0),
            ('AmorphSigma', 23.0,  0.0,   10_000.0),
            ('LamThick',   146.0,  0.1,   10_000.0),
            ('LamSigma',    65.0,  0.0,   10_000.0),
            ('G2',           1.0,  1e-30, None),
            ('Rg2',         50.0,  0.1,   10_000.0),
            ('G3',           0.1,  1e-30, None),
            ('Rg3',        200.0,  0.1,   100_000.0),
            ('B3',          1e-4,  1e-30, None),
            ('P3',           2.0,  0.1,   6.0),
        ],
        'formula': _hybrid_hermans,
        'linearization': None,
        'complex_bg': True,
    },
    'Unified Born Green': {
        'params': [
            ('Rg1',        50.0,  0.1,   10_000.0),
            ('B1',         1e-3,  1e-30, None),
            ('Pack',      100.0,  0.0,   None),
            ('CorrDist',  200.0,  0.1,   100_000.0),
            ('StackIrreg', 0.17,  None,  None),
            ('kI',         0.02,  0.0,   None),
        ],
        'formula': _ubg,
        'linearization': None,
        'complex_bg': True,
    },
}

# Ordered list of all model names (for combo boxes etc.)
MODEL_NAMES: list[str] = list(MODEL_REGISTRY.keys())


# ===========================================================================
# SimpleFitModel
# ===========================================================================

class SimpleFitModel:
    """
    Single-model fitting for SAS data.

    Attributes
    ----------
    model : str
        Active model name (must be a key in ``MODEL_REGISTRY``).
    params : dict[str, float]
        Current parameter values.  Populated from registry defaults when
        the model is changed via ``set_model()``.
    limits : dict[str, tuple]
        Fitting bounds ``(lower, upper)`` for each parameter.
        ``None`` entries map to ±∞ for ``scipy.optimize.curve_fit``.
    use_complex_bg : bool
        When True (and the model supports it), add a complex background
        ``BG_A·Q^(−BG_n) + BG_flat`` as extra fit parameters.
    n_mc_runs : int
        Number of Monte Carlo runs used by the batch ``fit_simple()``
        uncertainty estimation.
    """

    def __init__(self) -> None:
        self.model: str = 'Guinier'
        self.params: dict[str, float] = {}
        self.limits: dict[str, tuple] = {}
        self.use_complex_bg: bool = False
        self.n_mc_runs: int = 50
        self._reset_to_defaults()

    # ── Model switching ───────────────────────────────────────────────────────

    def set_model(self, model_name: str) -> None:
        """Switch to a different model and reset params to registry defaults."""
        if model_name not in MODEL_REGISTRY:
            raise ValueError(
                f"Unknown model '{model_name}'.  "
                f"Available: {MODEL_NAMES}"
            )
        self.model = model_name
        self._reset_to_defaults()

    def _reset_to_defaults(self) -> None:
        """Populate params/limits from registry defaults for the current model."""
        entry = MODEL_REGISTRY[self.model]
        self.params = {name: default for name, default, lo, hi in entry['params']}
        self.limits = {name: (lo, hi) for name, default, lo, hi in entry['params']}
        if self.use_complex_bg and entry['complex_bg']:
            for name, default, lo, hi in _BG_PARAMS:
                self.params.setdefault(name, default)
                self.limits.setdefault(name, (lo, hi))

    # ── Forward model evaluation ───────────────────────────────────────────────

    def compute(self, q: np.ndarray) -> np.ndarray:
        """Evaluate the model at *q* with current parameter values (no fitting).

        Returns the model intensity array (same shape as *q*).
        Negative or non-finite values are replaced with NaN.
        """
        q = np.asarray(q, dtype=float)
        func = self._build_fit_func()
        specs = self._active_param_specs()
        vals = [s[1] for s in specs]
        with np.errstate(divide='ignore', invalid='ignore'):
            result = func(q, *vals)
        return np.where(np.isfinite(result) & (result > 0), result, np.nan)

    # ── Helpers ───────────────────────────────────────────────────────────────

    def _active_param_specs(self) -> list[tuple]:
        """Return ordered list of (name, value, lo, hi) for the active model."""
        entry = MODEL_REGISTRY[self.model]
        specs = []
        for name, default, lo_def, hi_def in entry['params']:
            val = self.params.get(name, default)
            lo, hi = self.limits.get(name, (lo_def, hi_def))
            specs.append((name, val, lo, hi))
        if self.use_complex_bg and entry['complex_bg']:
            for name, default, lo_def, hi_def in _BG_PARAMS:
                val = self.params.get(name, default)
                lo, hi = self.limits.get(name, (lo_def, hi_def))
                specs.append((name, val, lo, hi))
        return specs

    def _build_fit_func(self):
        """Return a callable ``f(q, *params)`` for ``scipy.optimize.curve_fit``."""
        entry = MODEL_REGISTRY[self.model]
        formula = entry['formula']
        n_model = len(entry['params'])
        use_bg = self.use_complex_bg and entry['complex_bg']

        if use_bg:
            def func(q, *all_params):
                model_params = all_params[:n_model]
                BG_A, BG_n, BG_flat = all_params[n_model:]
                with np.errstate(divide='ignore', invalid='ignore'):
                    bg = np.where(q > 0, BG_A * q**(-BG_n), 0.0) + BG_flat
                return formula(q, *model_params) + bg
        else:
            def func(q, *all_params):
                return formula(q, *all_params)

        return func

    # ── Fitting ───────────────────────────────────────────────────────────────

    def fit(
        self,
        q: np.ndarray,
        intensity: np.ndarray,
        error: Optional[np.ndarray] = None,
        fixed_params: Optional[dict] = None,
        no_limits: bool = False,
    ) -> dict:
        """
        Fit the active model to experimental SAS data.

        Parameters
        ----------
        q : array
            Scattering vector [Å⁻¹].
        intensity : array
            Measured intensity.
        error : array or None
            Measurement uncertainty.  If None, set to ``intensity × 0.05``.
        fixed_params : dict or None
            Parameters that should be held fixed during fitting.
            Keys are parameter names; values are the fixed values to use.
            Parameters not in this dict are fitted normally.
        no_limits : bool
            When True, ignore all lower/upper bounds (fit unconstrained).

        Returns
        -------
        dict with keys:
            success, model, params, params_std, I_model, q, residuals,
            chi2, dof, derived, error (error message if success is False).
        """
        q = np.asarray(q, dtype=float)
        I = np.asarray(intensity, dtype=float)
        if error is None:
            dI = np.maximum(I * 0.05, 1e-30)
        else:
            dI = np.asarray(error, dtype=float)
            dI = np.maximum(dI, 1e-30)

        # Filter to finite, positive q points
        mask = np.isfinite(q) & np.isfinite(I) & np.isfinite(dI) & (q > 0)
        qf, If, dIf = q[mask], I[mask], dI[mask]

        if len(qf) < 2:
            return self._failure('Not enough valid data points (need ≥ 2).')

        fixed = dict(fixed_params) if fixed_params else {}
        specs = self._active_param_specs()
        all_names = [s[0] for s in specs]

        # Split parameters into free (to be fitted) and fixed
        free_specs = [(n, v, lo, hi) for (n, v, lo, hi) in specs if n not in fixed]
        free_names = [s[0] for s in free_specs]

        base_func = self._build_fit_func()

        # ── All parameters fixed: just evaluate the model ─────────────────────
        if not free_names:
            all_vals = [fixed.get(n, v) for (n, v, lo, hi) in specs]
            try:
                I_model_fit = base_func(qf, *all_vals)
                residuals_fit = (If - I_model_fit) / dIf
                chi2 = float(np.sum(residuals_fit**2))
                dof = max(len(qf), 1)
                I_model = base_func(q, *all_vals)
                residuals_full = (I - I_model) / dI
            except Exception as exc:
                return self._failure(str(exc))
            all_fitted = {n: fixed.get(n, v) for (n, v, lo, hi) in specs}
            self.params.update(all_fitted)
            derived = self._compute_derived(all_fitted)
            return {
                'success': True,
                'model': self.model,
                'params': all_fitted,
                'params_std': {n: 0.0 for n in all_names},
                'I_model': I_model,
                'q': q,
                'residuals': residuals_full,
                'chi2': chi2,
                'dof': dof,
                'reduced_chi2': chi2 / dof,
                'derived': derived,
                'error': None,
            }

        # ── Build wrapper that only accepts free parameters ───────────────────
        def func(q_arr, *free_vals):
            free_map = dict(zip(free_names, free_vals))
            all_vals = [
                fixed[n] if n in fixed else free_map[n]
                for n in all_names
            ]
            return base_func(q_arr, *all_vals)

        p0 = [s[1] for s in free_specs]
        if no_limits:
            bounds_lo = [-np.inf] * len(free_names)
            bounds_hi = [np.inf] * len(free_names)
        else:
            bounds_lo = [-np.inf if s[2] is None else s[2] for s in free_specs]
            bounds_hi = [np.inf if s[3] is None else s[3] for s in free_specs]

        # Clamp p0 to bounds
        p0 = [max(lo, min(hi, v)) if not (np.isinf(lo) and np.isinf(hi)) else v
              for v, lo, hi in zip(p0, bounds_lo, bounds_hi)]

        try:
            popt, pcov = curve_fit(
                func, qf, If,
                p0=p0,
                sigma=dIf,
                absolute_sigma=True,
                bounds=(bounds_lo, bounds_hi),
                maxfev=100_000,
            )
        except Exception as exc:
            return self._failure(str(exc))

        # Combine free fitted values + fixed values into full result dict
        free_fitted = dict(zip(free_names, popt))
        fitted_params = {}
        for (n, v, lo, hi) in specs:
            fitted_params[n] = fixed[n] if n in fixed else free_fitted[n]
        fitted_params = {k: float(v) for k, v in fitted_params.items()}
        self.params.update(fitted_params)

        # Parameter uncertainties (0 for fixed params, from pcov for free)
        param_std: dict[str, float] = {n: 0.0 for n in fixed}
        try:
            diag = np.diag(pcov)
            for name, d in zip(free_names, diag):
                param_std[name] = float(np.sqrt(max(d, 0.0)))
        except Exception:
            for name in free_names:
                param_std[name] = float('nan')

        # Model intensity and residuals (using free_vals order)
        free_vals_opt = [free_fitted[n] for n in free_names]
        I_model = func(q, *free_vals_opt)
        residuals_full = (I - I_model) / dI
        I_model_fit = func(qf, *free_vals_opt)
        residuals_fit = (If - I_model_fit) / dIf
        chi2 = float(np.sum(residuals_fit**2))
        dof = max(len(qf) - len(free_names), 1)

        derived = self._compute_derived(fitted_params)

        return {
            'success': True,
            'model': self.model,
            'params': fitted_params,
            'params_std': param_std,
            'I_model': I_model,
            'q': q,
            'residuals': residuals_full,
            'chi2': chi2,
            'dof': dof,
            'reduced_chi2': chi2 / dof,
            'derived': derived,
            'error': None,
        }

    def _failure(self, message: str) -> dict:
        return {
            'success': False,
            'model': self.model,
            'params': dict(self.params),
            'params_std': {},
            'I_model': None,
            'q': None,
            'residuals': None,
            'chi2': None,
            'dof': None,
            'reduced_chi2': None,
            'derived': {},
            'error': message,
        }

    def _compute_derived(self, fitted: dict) -> dict:
        """Compute model-specific derived quantities from fitted parameters."""
        derived: dict[str, float] = {}
        if self.model == 'Guinier Sheet':
            Rg = fitted.get('Rg', 0.0)
            derived['Thickness'] = Rg * 2.0 * np.sqrt(3.0)
        elif self.model == 'Treubner-Strey':
            A = fitted.get('A', 1e-4)
            C1 = fitted.get('C1', 0.0)
            C2 = fitted.get('C2', 1e-8)
            try:
                xi_arg = 0.5 * np.sqrt(abs(A / C2)) + C1 / (4.0 * C2)
                derived['CorrLength'] = 1.0 / np.sqrt(abs(xi_arg)) if xi_arg != 0 else float('nan')
                d_arg = 0.5 * np.sqrt(abs(A / C2)) - C1 / (4.0 * C2)
                derived['RepeatDist'] = (2.0 * np.pi / np.sqrt(abs(d_arg))
                                         if d_arg > 0 else float('nan'))
            except Exception:
                pass
        elif self.model == 'Unified Born Green':
            Rg1 = fitted.get('Rg1', 50.0)
            B1 = fitted.get('B1', 1e-3)
            Pack = fitted.get('Pack', 0.5)
            CorrDist = fitted.get('CorrDist', 100.0)
            if Rg1 > 0 and CorrDist > 0:
                Rad = np.sqrt(abs(Pack * CorrDist * Rg1) / 4.0)
                derived['Rad'] = float(Rad)
                denom12 = 2.0 * Rg1 + Rad
                if denom12 > 0:
                    derived['G1'] = float(B1 * Rg1**4 * Rad / denom12)
                Rg2 = np.sqrt(Rad**2 / 2.0 + Rg1**2 / 3.0)
                derived['Rg2'] = float(Rg2)
        return derived

    # ── Linearization helpers ─────────────────────────────────────────────────

    def linearize(
        self,
        q: np.ndarray,
        intensity: np.ndarray,
        error: Optional[np.ndarray] = None,
        q_min: Optional[float] = None,
        q_max: Optional[float] = None,
    ) -> Optional[dict]:
        """
        Return linearized (X, Y, dY) arrays and a model line for the current
        model.  Returns None if the model has no linearization.

        Data points are background-corrected: when ``use_complex_bg`` is True
        the fitted complex background (BG_A·Q⁻BG_n + BG_flat) is subtracted
        before applying the log/Porod transform so that the resulting data
        points are linear in the transformed space.

        The model line (``Y_fit``) is computed **analytically** from the
        currently fitted model parameters, not by linear regression.  This
        ensures the line matches the fit shown in the main I(Q) plot even
        when background is present.

        slope/intercept are derived analytically from the model params and
        displayed in the panel title; R² is computed by comparing the
        (background-subtracted) **in-range** data points to the analytic model
        line.  Out-of-range points (outside [q_min, q_max]) are included in the
        returned X/Y arrays for display (plotted grey in the panel) but are
        excluded from the R² calculation so that curvature outside the Guinier
        regime does not corrupt the goodness-of-fit metric.

        Returns
        -------
        dict with keys:
            x_label, y_label, X, Y, dY, X_fit, Y_fit,
            slope, intercept, r_squared
        """
        lin = MODEL_REGISTRY[self.model]['linearization']
        if lin is None:
            return None

        q = np.asarray(q, dtype=float)
        I = np.asarray(intensity, dtype=float)
        if error is None:
            dI = np.maximum(I * 0.05, 1e-30)
        else:
            dI = np.asarray(error, dtype=float)
            dI = np.maximum(dI, 1e-30)

        mask = np.isfinite(q) & np.isfinite(I) & (q > 0) & (I > 0)
        q, I, dI = q[mask], I[mask], dI[mask]

        # ── Subtract complex background (Guinier-family only) ─────────────────
        # Porod and Power Law have an explicit Background param in their formula;
        # they are not combined with the complex background widget.
        I_corr = I.copy()
        if lin != 'porod' and self.use_complex_bg and MODEL_REGISTRY[self.model].get('complex_bg', False):
            BG_A    = self.params.get('BG_A',    0.0)
            BG_n    = self.params.get('BG_n',    4.0)
            BG_flat = self.params.get('BG_flat', 0.0)
            with np.errstate(divide='ignore', invalid='ignore'):
                bg = np.where(q > 0, BG_A * q**(-BG_n), 0.0) + BG_flat
            I_corr = I - bg

        # Only keep points where background-corrected intensity is positive
        pos = I_corr > 0
        q, I_corr, dI = q[pos], I_corr[pos], dI[pos]

        if len(q) < 2:
            return None

        # ── Transform data ────────────────────────────────────────────────────
        if lin == 'guinier':
            # I_model = I0·exp(−Q²·Rg²/3)  →  ln(I_model) = ln(I0) − Q²·Rg²/3
            X = q**2
            Y = np.log(I_corr)
            dY = dI / I_corr           # d(ln I) = dI/I
            x_label, y_label = 'Q²  [Å⁻²]', 'ln(I)'
            I0 = self.params.get('I0', 1.0)
            Rg = self.params.get('Rg', 50.0)
            slope     = -Rg**2 / 3.0
            intercept = np.log(max(I0, 1e-300))
        elif lin == 'guinier_rod':
            # Q·I_model = I0·exp(−Q²·Rc²/2)  →  ln(Q·I) = ln(I0) − Q²·Rc²/2
            X = q**2
            Y = np.log(q * I_corr)
            dY = dI / I_corr
            x_label, y_label = 'Q²  [Å⁻²]', 'ln(Q·I)'
            I0 = self.params.get('I0', 1.0)
            Rc = self.params.get('Rc', 10.0)
            slope     = -Rc**2 / 2.0
            intercept = np.log(max(I0, 1e-300))
        elif lin == 'guinier_sheet':
            # Q²·I_model = I0·exp(−Q²·Rg²)  →  ln(Q²·I) = ln(I0) − Q²·Rg²
            X = q**2
            Y = np.log(q**2 * I_corr)
            dY = dI / I_corr
            x_label, y_label = 'Q²  [Å⁻²]', 'ln(Q²·I)'
            I0 = self.params.get('I0', 1.0)
            Rg = self.params.get('Rg', 10.0)
            slope     = -Rg**2
            intercept = np.log(max(I0, 1e-300))
        elif lin == 'porod':
            # I = Kp·Q⁻⁴ + Background  →  I·Q⁴ = Kp + Background·Q⁴
            X = q**4
            Y = I_corr * q**4          # I_corr == I for Porod (no bg subtraction)
            dY = dI * q**4
            x_label, y_label = 'Q⁴  [Å⁻⁴]', 'I·Q⁴'
            Kp         = self.params.get('Kp',         1.0)
            Background = self.params.get('Background', 0.0)
            slope     = Background      # coefficient of Q⁴
            intercept = Kp             # constant term
        else:
            return None

        # ── Analytic model line ───────────────────────────────────────────────
        Y_fit = slope * X + intercept

        # ── R² of in-range data vs analytic model ─────────────────────────────
        # R² is computed ONLY for points inside [q_min, q_max] (the fitting
        # range).  Including out-of-range points would give a misleadingly bad
        # (even negative) R² because the Guinier/Porod approximation naturally
        # breaks down outside that region.
        in_range = np.ones(len(q), dtype=bool)
        if q_min is not None:
            in_range &= (q >= q_min)
        if q_max is not None:
            in_range &= (q <= q_max)

        Y_r  = Y[in_range]
        Yf_r = Y_fit[in_range]
        dY_r = dY[in_range]

        if len(Y_r) < 2:
            r_sq = float('nan')
        else:
            w = 1.0 / np.maximum(dY_r, 1e-30)**2
            sw   = np.sum(w)
            Y_mean_w = np.sum(w * Y_r) / sw
            ss_res = np.sum(w * (Y_r - Yf_r)**2)
            ss_tot = np.sum(w * (Y_r - Y_mean_w)**2)
            r_sq = 1.0 - ss_res / ss_tot if ss_tot > 0 else float('nan')

        return {
            'x_label':   x_label,
            'y_label':   y_label,
            'q':         q,         # original Q values (for Q-range highlighting)
            'X':         X,
            'Y':         Y,
            'dY':        dY,
            'X_fit':     X,
            'Y_fit':     Y_fit,
            'slope':     float(slope),
            'intercept': float(intercept),
            'r_squared': float(r_sq),
        }

    # ── Serialisation ─────────────────────────────────────────────────────────

    def to_dict(self) -> dict:
        """Serialise to a plain dict (for state persistence and JSON export)."""
        return {
            'model': self.model,
            'params': dict(self.params),
            'limits': {k: list(v) for k, v in self.limits.items()},
            'use_complex_bg': self.use_complex_bg,
            'n_mc_runs': self.n_mc_runs,
        }

    @classmethod
    def from_dict(cls, d: dict) -> 'SimpleFitModel':
        """Deserialise from a plain dict."""
        obj = cls.__new__(cls)
        obj.model = d.get('model', 'Guinier')
        if obj.model not in MODEL_REGISTRY:
            obj.model = 'Guinier'
        obj.params = dict(d.get('params', {}))
        raw_lims = d.get('limits', {})
        obj.limits = {k: tuple(v) for k, v in raw_lims.items()}
        obj.use_complex_bg = bool(d.get('use_complex_bg', False))
        obj.n_mc_runs = int(d.get('n_mc_runs', 50))
        # Fill in any missing params/limits from registry defaults
        entry = MODEL_REGISTRY[obj.model]
        for name, default, lo, hi in entry['params']:
            obj.params.setdefault(name, default)
            obj.limits.setdefault(name, (lo, hi))
        if obj.use_complex_bg and entry['complex_bg']:
            for name, default, lo, hi in _BG_PARAMS:
                obj.params.setdefault(name, default)
                obj.limits.setdefault(name, (lo, hi))
        return obj
