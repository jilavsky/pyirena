"""
Similarity analysis for SAXS datasets — outlier / radiation-damage detection.

Implements the CorMap test (Franke et al., Nature Methods 12, 419-422, 2015),
based on the Schilling (1990) "Longest Run of Heads" probability formula.

The module is intentionally free of intra-package imports so it can be used
independently and extended with additional methods without coupling to the
averaging engine.

Extension pattern
-----------------
Add a new method function with signature::

    def _run_my_method(datasets, filenames, reference, p_min) -> list[SimilarityResult]

then register it in :data:`SIMILARITY_METHODS`::

    SIMILARITY_METHODS['my_method'] = _run_my_method

The GUI method dropdown is populated from ``SIMILARITY_METHODS.keys()``.
"""

from __future__ import annotations

import math
from functools import lru_cache
from typing import NamedTuple, Optional

import numpy as np


# ---------------------------------------------------------------------------
#  Public result type
# ---------------------------------------------------------------------------

class SimilarityResult(NamedTuple):
    idx: int        # position in the input datasets list
    filename: str   # display name; '' when not provided
    p_value: float  # 0–1; low p → curves differ → likely radiation damaged
    longest_run: int  # C: observed longest run of same-sign residuals
    n_points: int     # N: number of Q points compared
    accepted: bool    # True when p_value >= p_min


# ---------------------------------------------------------------------------
#  CorMap internals
# ---------------------------------------------------------------------------

def _longest_run(arr: np.ndarray) -> int:
    """Return the longest run of consecutive same-sign finite values."""
    finite = arr[np.isfinite(arr)]
    signs = np.sign(finite).astype(int)
    if len(signs) == 0:
        return 0
    changes = np.where(np.diff(signs) != 0)[0] + 1
    runs = np.diff(np.concatenate(([0], changes, [len(signs)])))
    return int(runs.max()) if len(runs) else 0


@lru_cache(maxsize=None)
def _A(n: int, c: int) -> int:
    """Count binary sequences of length n with no run of same symbol >= c.

    Schilling (1990), Recurrence 2.  Memoised via lru_cache.
    """
    if n <= 0:
        return 1
    if n < c:
        return 1 << n  # 2^n
    return sum(_A(n - 1 - j, c) for j in range(c))


def _cormap_p_value(n: int, c: int) -> float:
    """P(longest run of same-sign residuals >= c | n i.i.d. ±1 trials).

    Low p-value means the observed run is unlikely to arise by chance, i.e.
    the two curves are statistically different.
    """
    if n <= 0 or c <= 1:
        return 1.0
    # B(n, c) = number of sequences with no run of + OR - of length >= c
    # delta = sequences where the longest run IS >= c
    delta = (1 << n) - 2 * _A(n - 1, c - 1)
    if delta <= 0:
        return 0.0
    # Convert to float in log space to handle very large integers safely
    try:
        return min(math.exp2(math.log2(float(delta)) - n), 1.0)
    except (ValueError, OverflowError):
        return 0.0


def _interp_loglog(q_ref: np.ndarray, q_src: np.ndarray, I_src: np.ndarray) -> np.ndarray:
    """Interpolate I_src(q_src) onto q_ref in log-log space."""
    return np.exp(
        np.interp(
            np.log(q_ref),
            np.log(q_src),
            np.log(np.maximum(I_src, 1e-300)),
            left=np.nan,
            right=np.nan,
        )
    )


def _compare_pair(
    q1: np.ndarray, I1: np.ndarray,
    q2: np.ndarray, I2: np.ndarray,
) -> tuple[int, int, float]:
    """Compare two SAXS curves with CorMap on their common Q range.

    Returns (N, C, p_value).
    """
    q_lo = max(float(q1[0]), float(q2[0]))
    q_hi = min(float(q1[-1]), float(q2[-1]))
    if q_lo >= q_hi:
        return 0, 0, 1.0

    mask = (q1 >= q_lo) & (q1 <= q_hi)
    q_ref = q1[mask]
    I1_ref = I1[mask]
    if len(q_ref) < 2:
        return 0, 0, 1.0

    I2_ref = _interp_loglog(q_ref, q2, I2)
    diff = I2_ref - I1_ref

    C = _longest_run(diff)
    N = int(np.sum(np.isfinite(diff)))
    return N, C, _cormap_p_value(N, C)


def _build_majority_reference(
    datasets: list,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute a median reference curve from all datasets on the first Q grid.

    Median is used instead of mean so that a single large outlier does not
    bias the reference and cause all other datasets to be rejected (the
    "swamping" problem with outlier detection).
    """
    q_ref = datasets[0][0]
    I_stack = []
    for i, (q, I, *_) in enumerate(datasets):
        if i == 0:
            I_stack.append(I.copy())
        else:
            I_interp = _interp_loglog(q_ref, q, I)
            I_stack.append(I_interp)
    I_median = np.nanmedian(np.vstack(I_stack), axis=0)
    return q_ref, I_median


# ---------------------------------------------------------------------------
#  CorMap method driver
# ---------------------------------------------------------------------------

def _run_cormap(
    datasets: list,
    filenames: list[str],
    reference: str,
    p_min: float,
) -> list[SimilarityResult]:
    """Run CorMap similarity check on *datasets*.

    Parameters
    ----------
    datasets:
        List of ``(q, I, dI, dQ)`` tuples (same format as
        ``DataManipulation.average()``).
    filenames:
        Display names, parallel to *datasets*.
    reference:
        ``'first'`` — compare every frame against frame 0 (frame 0 always
        accepted with p=1.0); ``'majority'`` — compare every frame against the
        mean of all frames (frame 0 is also tested).
    p_min:
        Rejection threshold; frame rejected when ``p_value < p_min``.
    """
    results: list[SimilarityResult] = []

    if reference == 'majority':
        q_ref, I_ref = _build_majority_reference(datasets)
        for i, (q, I, *_) in enumerate(datasets):
            N, C, p = _compare_pair(q_ref, I_ref, q, I)
            results.append(SimilarityResult(
                idx=i, filename=filenames[i],
                p_value=p, longest_run=C, n_points=N,
                accepted=(p >= p_min),
            ))
    else:  # 'first'
        # Frame 0 is the reference — always accepted
        q0, I0 = datasets[0][0], datasets[0][1]
        results.append(SimilarityResult(
            idx=0, filename=filenames[0],
            p_value=1.0, longest_run=0, n_points=0,
            accepted=True,
        ))
        for i, (q, I, *_) in enumerate(datasets[1:], start=1):
            N, C, p = _compare_pair(q0, I0, q, I)
            results.append(SimilarityResult(
                idx=i, filename=filenames[i],
                p_value=p, longest_run=C, n_points=N,
                accepted=(p >= p_min),
            ))

    return results


# ---------------------------------------------------------------------------
#  Extensibility registry  (add new methods here)
# ---------------------------------------------------------------------------

SIMILARITY_METHODS: dict[str, callable] = {
    'cormap': _run_cormap,
    # 'my_method': _run_my_method,
}

#: Human-readable labels for the GUI dropdown, keyed by method name.
SIMILARITY_METHOD_LABELS: dict[str, str] = {
    'cormap': 'CorMap (Franke 2015)',
}


# ---------------------------------------------------------------------------
#  Public entry point
# ---------------------------------------------------------------------------

def check_similarity(
    datasets: list,
    filenames: Optional[list[str]] = None,
    method: str = 'cormap',
    reference: str = 'first',
    p_min: float = 0.01,
) -> list[SimilarityResult]:
    """Assess pairwise similarity of SAXS datasets and flag outliers.

    Parameters
    ----------
    datasets:
        List of ``(q, I, dI, dQ)`` tuples — the same format accepted by
        ``DataManipulation.average()``.
    filenames:
        Optional display names, one per dataset.  Defaults to ``''``.
    method:
        Key into :data:`SIMILARITY_METHODS`.  Currently ``'cormap'``.
    reference:
        ``'first'`` (compare each frame vs. frame 0; frame 0 always accepted)
        or ``'majority'`` (compare each frame vs. the mean of all frames).
    p_min:
        P-value threshold.  Frames with ``p_value < p_min`` are marked as
        rejected (``accepted=False``).  Typical range: 0.001–0.05.

    Returns
    -------
    list[SimilarityResult]
        One entry per dataset, in input order.
    """
    if filenames is None:
        filenames = [''] * len(datasets)
    if len(filenames) != len(datasets):
        raise ValueError("filenames must have the same length as datasets")
    if method not in SIMILARITY_METHODS:
        raise ValueError(f"Unknown similarity method {method!r}; "
                         f"available: {list(SIMILARITY_METHODS)}")
    if reference not in ('first', 'majority'):
        raise ValueError(f"reference must be 'first' or 'majority', got {reference!r}")
    if not datasets:
        return []
    if len(datasets) == 1:
        q, I, *_ = datasets[0]
        return [SimilarityResult(
            idx=0, filename=filenames[0],
            p_value=1.0, longest_run=0, n_points=len(q),
            accepted=True,
        )]

    return SIMILARITY_METHODS[method](datasets, filenames, reference, p_min)
