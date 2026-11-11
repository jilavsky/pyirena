"""Shared JSON helpers hoisted out of the control modules (issue #38)."""
import math

import numpy as np
import pytest

from pyirena.api._json import decimate_for_json, json_float


@pytest.mark.parametrize("value, expected", [
    (1, 1.0), (np.float32(2.5), 2.5), ("3", 3.0),
    (None, None), ("x", None), (float("nan"), None), (np.inf, None),
])
def test_json_float(value, expected):
    assert json_float(value) == expected


@pytest.mark.parametrize("n", [10, 500, 501, 999, 1000, 1001, 5000])
def test_decimate_never_exceeds_cap(n):
    out = decimate_for_json(np.arange(n, dtype=float), max_points=500)
    assert len(out) <= 500
    assert out[0] == 0.0
    if n <= 500:
        assert len(out) == n


def test_decimate_keeps_alignment_and_nulls_nonfinite():
    q = np.linspace(1, 2, 999)
    I = q.copy()
    I[3] = np.nan
    a, b = decimate_for_json(q, 300), decimate_for_json(I, 300)
    assert len(a) == len(b)
    assert all(x == y for x, y in zip(a, b, strict=True) if y is not None)
    assert decimate_for_json([1.0, math.inf]) == [1.0, None]
