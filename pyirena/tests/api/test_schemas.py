"""Tests for pyirena.api.schemas — array sanitization, JSON-safe output."""
from __future__ import annotations

import math

import numpy as np

from pyirena.api.schemas import (
    FileEntry, ReducedData, array_to_list, asdict_clean,
)


def test_array_to_list_decimates():
    arr = np.arange(1000.0)
    out = array_to_list(arr, max_points=50)
    assert len(out) == 50
    # Endpoints preserved
    assert out[0] == 0.0
    assert out[-1] == 999.0


def test_array_to_list_handles_nan_inf():
    arr = np.array([1.0, np.nan, np.inf, -np.inf, 2.0])
    out = array_to_list(arr, max_points=10, include_full=True)
    assert out == [1.0, None, None, None, 2.0]


def test_array_to_list_none_passthrough():
    assert array_to_list(None) is None


def test_array_to_list_empty():
    out = array_to_list(np.array([]))
    assert out == []


def test_asdict_clean_replaces_nan():
    fe = FileEntry(path="/x", name="x.h5", scan_number=None,
                   mtime=None, size_bytes=10)
    d = asdict_clean(fe)
    assert d["path"] == "/x"
    assert d["scan_number"] is None
    assert d["analyses"] == []


def test_reduced_data_to_dict_is_json_safe():
    rd = ReducedData(path="/x", found=True, n_points=3,
                     q_min=0.001, q_max=1.0,
                     Q=[0.001, 0.1, 1.0],
                     I=[1.0, math.nan, 0.5])
    d = rd.to_dict()
    # NaN must be None after to_dict()
    assert d["I"][1] is None
    # JSON serialization round-trips cleanly
    import json
    s = json.dumps(d)
    assert json.loads(s)["found"] is True
