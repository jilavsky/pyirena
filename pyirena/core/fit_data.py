"""
Which points a fit sees, and with what σ — one rule for every path.

Decided for issue #38 (``planning/fit-runner/README.md``, D2/D3), after an
audit found every path doing something different:

* **I ≤ 0 points are kept.** Background-subtracted data have real negative
  values at high Q; dropping only the negative half of the noise biases a
  fitted background upward. (Text-file import keeps its own cleaning, and
  core math that needs I > 0 — e.g. the Size Distribution inversion — masks
  internally as before.)
* **σ ≤ 0 or non-finite σ is repaired** to ``|I| × error_fraction`` — the
  Data Selector's "Generated uncertainty fraction", the same rule text import
  uses. It used to be clamped to 1e-30 on some paths, giving that one point
  ~1e60 weight so it dominated χ².
* A point whose σ cannot be repaired (I = 0 and no valid σ) is dropped and
  counted — the alternative is an infinite weight.
* Non-finite Q or I, and Q ≤ 0, are dropped (nothing can be fitted there).
* ``err=None`` (no uncertainties at all) stays None: each tool keeps its own
  weighting for that case, and Monte-Carlo uses the generated fraction.

numpy only.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np


@dataclass
class FitArrays:
    """The prepared arrays plus what was done to get them."""

    q: np.ndarray
    intensity: np.ndarray
    error: Optional[np.ndarray]
    #: Boolean mask over the *input* arrays selecting the kept points.
    mask: np.ndarray
    n_sigma_repaired: int = 0
    n_dropped: int = 0

    def notes(self) -> list:
        out = []
        if self.n_sigma_repaired:
            out.append(f"{self.n_sigma_repaired} point(s) with uncertainty <= 0 "
                       "were given σ = I × the generated-uncertainty fraction.")
        if self.n_dropped:
            out.append(f"{self.n_dropped} point(s) dropped (non-finite Q/I, Q <= 0, "
                       "or I = 0 with no valid uncertainty).")
        return out


def prepare_fit_arrays(q, intensity, error=None, error_fraction: float = 0.05,
                       base_mask=None) -> FitArrays:
    """Apply the shared point/σ rule. See the module docstring.

    Args:
        q, intensity, error: Input arrays (error may be None).
        error_fraction: Fraction of |I| used for σ ≤ 0 points.
        base_mask: Optional extra selection (e.g. the fit Q range), ANDed in.
    """
    q = np.asarray(q, dtype=float)
    I = np.asarray(intensity, dtype=float)
    mask = np.isfinite(q) & (q > 0) & np.isfinite(I)
    if base_mask is not None:
        mask &= np.asarray(base_mask, dtype=bool)
    n_repaired = 0
    err_out = None
    if error is not None:
        e = np.array(error, dtype=float)          # copy: never edit the caller's
        bad = ~np.isfinite(e) | (e <= 0)
        repair = bad & mask
        e[repair] = np.abs(I[repair]) * float(error_fraction)
        still_bad = mask & (~np.isfinite(e) | (e <= 0))
        n_repaired = int(np.sum(repair & ~still_bad))
        mask &= ~still_bad
        err_out = e[mask]
    n_in = int(np.sum(base_mask)) if base_mask is not None else len(q)
    return FitArrays(q=q[mask], intensity=I[mask], error=err_out, mask=mask,
                     n_sigma_repaired=n_repaired,
                     n_dropped=max(0, n_in - int(np.sum(mask))))


def q_range_mask(q, q_min=None, q_max=None) -> np.ndarray:
    """Inclusive Q-range selection; either bound may be None (open)."""
    q = np.asarray(q, dtype=float)
    m = np.ones(q.shape, dtype=bool)
    if q_min is not None:
        m &= q >= float(q_min)
    if q_max is not None:
        m &= q <= float(q_max)
    return m
