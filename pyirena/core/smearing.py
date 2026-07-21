"""
Slit-smearing engine for pyIrena.

Single home for all slit-smearing math — no tool implements its own.  The core
operation is the APS/USAXS (Lake) infinite-slit-length smearing of an *ideal*
(pinhole) model:

    I_smeared(q) = (1 / SL) * ∫_0^SL  I_ideal( sqrt(q² + l²) )  dl

where ``SL`` is the slit (half-)length in 1/Å — the ``dQl`` scalar in NXcanSAS
slit-smeared data.  pyIrena smears the *model*, never the data (the data are
what the instrument measured); the smeared model is compared to slit-smeared
data during fitting and display.

Design contract
---------------
* **SL <= 0 is a strict no-op.**  Every entry point returns the input unchanged
  (identity matrix / identity array) so tools can call unconditionally without
  ``if use_smearing:`` branches.
* :func:`build_smearing_matrix` returns a *linear operator* ``W`` such that
  ``I_smeared(q) = W @ I_ideal(q_ext)``.  Because the smearing is linear in I,
  Size Distribution can smear its whole G matrix once (``G_sm = W @ G(q_ext)``)
  and the inversion machinery is untouched; fit loops cost one matvec per
  iteration.  Interpolation inside ``W`` is linear in *linear* intensity space
  (not log-log) precisely so ``W`` stays a fixed matrix independent of I.
* :func:`smear_curve` smears a *tabulated* curve (validation tests and the
  Fractals Monte-Carlo intensity).  It interpolates in log-log space and
  extrapolates explicitly, so it is a one-shot (non-linear-in-I) helper, not a
  matrix.

Igor reference: ``IR1B_SmearData`` in ``IR1_Desmearing.ipf`` and the model-side
extension in ``IR1A_UnifiedCalculateIntensity`` (``IR1_UnifiedFitFncts.ipf``).
The pyIrena improvement over Igor: because we control the model, we evaluate it
directly on an extended q grid (:func:`build_extended_q`) instead of
extrapolating a tabulated curve, so the "artifact at qmax" problem largely
disappears.

Pure NumPy/SciPy — no h5py, no Qt — unit-testable in isolation.
"""

from __future__ import annotations

from typing import Callable, Optional, Tuple

import numpy as np
from scipy import sparse

__all__ = [
    "build_extended_q",
    "build_smearing_matrix",
    "smear_model",
    "smear_curve",
    "SlitSmearer",
]

# Default number of trapezoid points across the slit length [0, SL] used to
# discretise the l-integral.  200 is ample for the accuracy validated against
# the DSM->SMR test-file pair (median |rel. dev.| ~0.1% mid-range); confirmed
# adequate by Jan (2026-07-20).
DEFAULT_N_L = 200

# Default geometric refinement of each interior q segment when building the
# extended grid used for W.  The l-integral samples I at sqrt(q² + l²), which
# for low q hops across several decades of the (log-spaced) data grid; refining
# keeps the linear-in-linear interpolation inside W accurate.  Model evaluation
# on the refined grid is analytic and cheap, so refinement is nearly free.
DEFAULT_REFINE = 8

# Irena target: extend the model grid to at least this multiple of SL so the
# smearing integral is well sampled even when the data barely span the slit.
DEFAULT_MIN_SPAN_FACTOR = 3.0


def build_extended_q(
    q: np.ndarray,
    slit_length: float,
    min_span_factor: float = DEFAULT_MIN_SPAN_FACTOR,
    refine: int = DEFAULT_REFINE,
) -> np.ndarray:
    """Extend (and optionally refine) a q grid for model-side slit smearing.

    The smearing integral evaluates the model at ``sqrt(q² + l²)`` for
    ``l`` in ``[0, SL]``, so the model must be known up to
    ``sqrt(qmax² + SL²)``.  Following Irena, we also require the grid to span at
    least ``min_span_factor * SL``.  Points beyond ``qmax`` are appended with a
    constant spacing of twice the last data step (Igor-compatible); the model is
    analytic there so a coarse tail is fine.

    The original ``q`` values are always a subset of the returned grid, so the
    corresponding rows of :func:`build_smearing_matrix` map cleanly back onto
    ``q``.

    Parameters
    ----------
    q : array-like
        Monotonically increasing data q grid [1/Å].
    slit_length : float
        Slit (half-)length ``SL`` [1/Å].  ``<= 0`` returns ``q`` unchanged.
    min_span_factor : float
        Extend at least to ``min_span_factor * SL`` past 0 (Irena uses 3).
    refine : int
        Insert ``refine - 1`` geometrically spaced points inside every interior
        segment of ``q`` (``1`` = no interior refinement).

    Returns
    -------
    np.ndarray
        The extended (and refined) q grid, strictly increasing, containing the
        original ``q`` values.
    """
    q = np.asarray(q, dtype=float)
    if slit_length <= 0 or q.size < 2:
        return q.copy()

    # --- interior geometric refinement (keeps original nodes) ----------------
    if refine and refine > 1:
        segments = []
        for i in range(q.size - 1):
            sub = np.geomspace(q[i], q[i + 1], refine + 1)[:-1]
            segments.append(sub)
        segments.append(q[-1:])
        base = np.concatenate(segments)
    else:
        base = q.copy()

    # --- extension beyond qmax ----------------------------------------------
    qmax = float(q[-1])
    target = max(np.sqrt(qmax ** 2 + slit_length ** 2), min_span_factor * slit_length)
    if target <= base[-1]:
        return base

    # Igor extends with 2*Δq_last for its tabulated curve; because we evaluate
    # an analytic model on this grid we can afford the (already-refined) last
    # step itself, which keeps the linear-in-linear interpolation of steep
    # (e.g. Guinier) tails accurate right up to sqrt(qmax²+SL²).
    step = base[-1] - base[-2]
    if step <= 0:
        step = qmax * 0.01
    n_ext = int(np.ceil((target - base[-1]) / step))
    ext = base[-1] + step * np.arange(1, n_ext + 1)
    return np.concatenate([base, ext])


def build_smearing_matrix(
    q: np.ndarray,
    slit_length: float,
    q_ext: Optional[np.ndarray] = None,
    n_l: int = DEFAULT_N_L,
    min_span_factor: float = DEFAULT_MIN_SPAN_FACTOR,
    refine: int = DEFAULT_REFINE,
) -> Tuple[sparse.csr_matrix, np.ndarray]:
    """Build the sparse slit-smearing operator ``W``.

    ``I_smeared(q) = W @ I_ideal(q_ext)`` where the rows of ``W`` are the
    trapezoid weights of ``(1/SL) ∫_0^SL I(sqrt(q² + l²)) dl`` after linear
    interpolation of the eval points onto ``q_ext``.

    Parameters
    ----------
    q : array-like
        Data q grid [1/Å] (the output grid; one row of ``W`` per point).
    slit_length : float
        Slit length ``SL`` [1/Å].  ``<= 0`` yields the identity operator with
        ``q_ext == q``.
    q_ext : array-like, optional
        Extended grid the model is evaluated on.  Built via
        :func:`build_extended_q` when omitted.
    n_l : int
        Number of trapezoid points across ``[0, SL]``.
    min_span_factor, refine :
        Forwarded to :func:`build_extended_q` when ``q_ext`` is built here.

    Returns
    -------
    (W, q_ext) : (scipy.sparse.csr_matrix, np.ndarray)
        ``W`` has shape ``(len(q), len(q_ext))``.
    """
    q = np.asarray(q, dtype=float)
    n = q.size

    if slit_length <= 0:
        # Strict no-op: identity onto q itself.
        return sparse.identity(n, format="csr", dtype=float), q.copy()

    if q_ext is None:
        q_ext = build_extended_q(
            q, slit_length, min_span_factor=min_span_factor, refine=refine
        )
    else:
        q_ext = np.asarray(q_ext, dtype=float)

    # Trapezoid weights over l in [0, SL], normalised by 1/SL.
    l = np.linspace(0.0, slit_length, n_l)
    dl = l[1] - l[0]
    w_l = np.full(n_l, dl, dtype=float)
    w_l[0] *= 0.5
    w_l[-1] *= 0.5
    w_l /= slit_length  # the (1/SL) prefactor

    # Eval points r[i, j] = sqrt(q_i² + l_j²).  Shape (n, n_l).
    r = np.sqrt(q[:, None] ** 2 + l[None, :] ** 2)

    # Linear-interpolation brackets of every r onto q_ext.
    idx = np.searchsorted(q_ext, r)
    idx = np.clip(idx, 1, q_ext.size - 1)
    x0 = q_ext[idx - 1]
    x1 = q_ext[idx]
    denom = x1 - x0
    frac = np.where(denom > 0, (r - x0) / denom, 0.0)
    frac = np.clip(frac, 0.0, 1.0)

    # Accumulate weights: each (i, j) contributes w_l[j] split between the two
    # bracketing columns of q_ext.
    row = np.repeat(np.arange(n), n_l)
    w_l_full = np.broadcast_to(w_l, (n, n_l)).ravel()
    frac_flat = frac.ravel()
    lo_idx = (idx - 1).ravel()
    hi_idx = idx.ravel()

    rows = np.concatenate([row, row])
    cols = np.concatenate([lo_idx, hi_idx])
    vals = np.concatenate([w_l_full * (1.0 - frac_flat), w_l_full * frac_flat])

    W = sparse.coo_matrix(
        (vals, (rows, cols)), shape=(n, q_ext.size), dtype=float
    ).tocsr()
    W.sum_duplicates()
    return W, q_ext


def smear_model(
    model_fn: Callable[[np.ndarray], np.ndarray],
    q: np.ndarray,
    slit_length: float,
    min_span_factor: float = DEFAULT_MIN_SPAN_FACTOR,
    refine: int = DEFAULT_REFINE,
    n_l: int = DEFAULT_N_L,
) -> np.ndarray:
    """Convenience: slit-smear an analytic model evaluated on the fly.

    Builds the extended grid and ``W``, evaluates ``model_fn(q_ext)``, and
    returns ``W @ model_fn(q_ext)`` on the original ``q`` grid.  For repeated
    calls with the same ``(q, SL)`` (e.g. a fit loop) hold a :class:`SlitSmearer`
    instead so ``W`` is built once.

    ``slit_length <= 0`` returns ``model_fn(q)`` unchanged.
    """
    q = np.asarray(q, dtype=float)
    if slit_length <= 0:
        return np.asarray(model_fn(q), dtype=float)
    W, q_ext = build_smearing_matrix(
        q, slit_length, n_l=n_l, min_span_factor=min_span_factor, refine=refine
    )
    return W @ np.asarray(model_fn(q_ext), dtype=float)


def smear_curve(
    q: np.ndarray,
    intensity: np.ndarray,
    slit_length: float,
    extrapolation: str = "power_law",
    n_l: int = DEFAULT_N_L,
) -> np.ndarray:
    """Slit-smear a *tabulated* curve (no model available).

    Used for validation (smear the desmeared data and compare to the file's SMR
    entry) and for the Fractals Monte-Carlo intensity, which is tabulated on its
    own grid.  Interpolation is in log-log space (intensities span decades) and
    values beyond ``qmax`` are supplied by ``extrapolation``:

    * ``"power_law"`` (default) — fit a power law to the last decade and extend.
    * ``"flat"``       — hold the last value.
    * ``"linear_zero"``— linear ramp to zero at ``20 * qmax`` (Igor-compatible).
    * ``"none"``       — raise if any eval point exceeds ``qmax``.

    ``slit_length <= 0`` returns ``intensity`` unchanged.

    .. warning::
       Log-log interpolation requires positive intensities, so any non-positive
       (or non-finite) points are dropped before interpolation.  For curves that
       cross zero — e.g. an over-subtracted background near the noise floor —
       this biases the smeared result slightly upward there.  Prefer
       :func:`smear_model` / :class:`SlitSmearer` when an analytic model is
       available (they smear in linear space with no such bias).
    """
    q = np.asarray(q, dtype=float)
    intensity = np.asarray(intensity, dtype=float)
    if slit_length <= 0:
        return intensity.copy()

    qmax = float(q[-1])
    l = np.linspace(0.0, slit_length, n_l)
    dl = l[1] - l[0]
    w_l = np.full(n_l, dl, dtype=float)
    w_l[0] *= 0.5
    w_l[-1] *= 0.5

    r = np.sqrt(q[:, None] ** 2 + l[None, :] ** 2)  # (n, n_l)
    r_max = float(r.max())

    # Interior: log-log linear interpolation on positive samples.
    pos = intensity > 0
    if np.count_nonzero(pos) < 2:
        raise ValueError("smear_curve needs at least two positive intensity points.")
    lq = np.log(q[pos])
    lI = np.log(intensity[pos])

    def eval_curve(x):
        x = np.asarray(x, dtype=float)
        out = np.empty_like(x)
        interior = x <= qmax
        # log-log interpolation within the tabulated range
        out[interior] = np.exp(np.interp(np.log(x[interior]), lq, lI))
        beyond = ~interior
        if np.any(beyond):
            out[beyond] = _extrapolate(x[beyond], q, intensity, qmax, extrapolation)
        return out

    if extrapolation == "none" and r_max > qmax * (1.0 + 1e-9):
        raise ValueError(
            f"smear_curve: evaluation reaches q={r_max:.4g} > qmax={qmax:.4g} "
            "but extrapolation='none'."
        )

    I_r = eval_curve(r)  # (n, n_l)
    return (I_r * w_l[None, :]).sum(axis=1) / slit_length


def _extrapolate(x, q, intensity, qmax, method):
    """Values of the tabulated curve beyond ``qmax`` per ``method``."""
    if method == "flat":
        return np.full_like(x, intensity[-1])
    if method == "linear_zero":
        # Linear ramp from I(qmax) at qmax to 0 at 20*qmax (Igor convention).
        q_zero = 20.0 * qmax
        frac = np.clip((q_zero - x) / (q_zero - qmax), 0.0, 1.0)
        return intensity[-1] * frac
    if method == "power_law":
        # Fit I = B q^-P over the last decade in log-log, then extend.
        last = q >= qmax / 10.0
        if np.count_nonzero(last & (intensity > 0)) >= 2:
            sel = last & (intensity > 0)
            slope, intercept = np.polyfit(np.log(q[sel]), np.log(intensity[sel]), 1)
            return np.exp(intercept + slope * np.log(x))
        return np.full_like(x, max(intensity[-1], 0.0))
    raise ValueError(f"Unknown extrapolation method: {method!r}")


class SlitSmearer:
    """Cached slit-smearing operator for a fixed ``(q, slit_length)``.

    The object a fitting loop holds: builds ``W`` and ``q_ext`` once, then
    :meth:`smear` / :meth:`smear_model` are cheap.  When ``slit_length <= 0``
    the smearer is a strict no-op (``is_noop`` is True), so tools construct one
    unconditionally and never branch on whether smearing is active.

    Examples
    --------
    >>> smearer = SlitSmearer(q, slit_length)
    >>> I_sm = smearer.smear_model(model.calculate_intensity)   # analytic model
    >>> I_sm = smearer.smear(I_on_q_ext)                        # pre-evaluated
    """

    def __init__(
        self,
        q: np.ndarray,
        slit_length: float,
        min_span_factor: float = DEFAULT_MIN_SPAN_FACTOR,
        refine: int = DEFAULT_REFINE,
        n_l: int = DEFAULT_N_L,
    ):
        self.q = np.asarray(q, dtype=float)
        self.slit_length = float(slit_length)
        self._W, self._q_ext = build_smearing_matrix(
            self.q,
            self.slit_length,
            n_l=n_l,
            min_span_factor=min_span_factor,
            refine=refine,
        )

    @property
    def is_noop(self) -> bool:
        """True when there is no slit length, so smearing is the identity."""
        return self.slit_length <= 0

    @property
    def q_ext(self) -> np.ndarray:
        """Extended grid the model must be evaluated on before :meth:`smear`."""
        return self._q_ext

    @property
    def matrix(self) -> sparse.csr_matrix:
        """The sparse smearing operator ``W`` (``I_sm(q) = W @ I(q_ext)``)."""
        return self._W

    def smear(self, intensity_on_q_ext: np.ndarray) -> np.ndarray:
        """Smear a model already evaluated on :attr:`q_ext`."""
        return self._W @ np.asarray(intensity_on_q_ext, dtype=float)

    def smear_model(
        self, model_fn: Callable[[np.ndarray], np.ndarray]
    ) -> np.ndarray:
        """Evaluate ``model_fn`` on :attr:`q_ext` and return the smeared result."""
        return self._W @ np.asarray(model_fn(self._q_ext), dtype=float)

    def smear_columns(self, matrix_on_q_ext: np.ndarray) -> np.ndarray:
        """Smear every column of a matrix evaluated on :attr:`q_ext`.

        Used by Size Distribution: ``G_smeared = smearer.smear_columns(G_ext)``
        where ``G_ext`` has shape ``(len(q_ext), n_r)``.
        """
        return self._W @ np.asarray(matrix_on_q_ext, dtype=float)
