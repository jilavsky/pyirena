"""Strict-JSON coercion for everything ``pyirena.api`` returns.

Layering invariant 3 says the api layer returns JSON-serialisable dicts only.
In practice that has been enforced by each function being careful, which is
enough for MCP (the stdio server tolerates NaN, and a human reads the result)
but not for a socket protocol. ``json.dumps`` emits bare ``NaN`` and
``Infinity`` by default — tokens that are **not** valid JSON. Python's own
``json.loads`` accepts them, so a pyIrena-to-pyIrena round trip looks fine and
a strict parser on the other end (Go, Rust, a JS ``JSON.parse``, most
orchestrators) rejects the whole reply. A single un-converged uncertainty is
enough to do it.

So the ZMQ service serialises with ``allow_nan=False`` and everything that
reaches the wire goes through :func:`to_strict_json` first:

* numpy scalars and arrays become Python numbers and lists,
* NaN and ±Inf become ``null`` — the honest JSON spelling of "no value",
* Paths, bytes, sets and tuples become strings and lists,
* dict keys become strings.

``null`` for NaN loses the distinction between "not a number" and "not
computed", which is the right trade here: both mean "do not use this value",
and the alternative is a reply the caller cannot parse at all.
"""
from __future__ import annotations

import json
import math
from pathlib import Path, PurePath
from typing import Any

import numpy as np


def to_strict_json(obj: Any) -> Any:
    """Return *obj* rebuilt from types that ``json.dumps(allow_nan=False)`` accepts.

    Non-finite floats become None. Unknown objects fall back to ``str(obj)``
    rather than raising: a reply that names an unexpected object is more use
    than no reply at all.
    """
    # Order matters: bool before int (bool is an int), and the numpy scalar
    # check before the plain-number check.
    if obj is None or isinstance(obj, (str, bool)):
        return obj

    if isinstance(obj, (np.bool_,)):
        return bool(obj)

    if isinstance(obj, (int, np.integer)):
        return int(obj)

    if isinstance(obj, (float, np.floating)):
        value = float(obj)
        return value if math.isfinite(value) else None

    if isinstance(obj, (np.complexfloating, complex)):
        return {"real": to_strict_json(obj.real), "imag": to_strict_json(obj.imag)}

    if isinstance(obj, np.ndarray):
        return [to_strict_json(v) for v in obj.tolist()]

    if isinstance(obj, dict):
        return {str(k): to_strict_json(v) for k, v in obj.items()}

    if isinstance(obj, (list, tuple, set, frozenset)):
        return [to_strict_json(v) for v in obj]

    if isinstance(obj, (Path, PurePath)):
        return str(obj)

    if isinstance(obj, (bytes, bytearray)):
        return obj.decode("utf-8", errors="replace")

    if hasattr(obj, "item") and callable(obj.item):     # any stray numpy scalar
        try:
            return to_strict_json(obj.item())
        except Exception:
            pass

    return str(obj)


def dumps_strict(obj: Any, **kwargs: Any) -> str:
    """Serialise *obj* as JSON that a strict parser will accept.

    Raises ``ValueError`` only if :func:`to_strict_json` missed something,
    which is a bug in this module rather than in the caller.
    """
    kwargs.setdefault("allow_nan", False)
    return json.dumps(to_strict_json(obj), **kwargs)


def is_strict_json(obj: Any) -> bool:
    """True if *obj* can already be serialised without coercion. For tests."""
    try:
        json.dumps(obj, allow_nan=False)
    except (TypeError, ValueError):
        return False
    return True


def json_float(value: Any):
    """Float, or None for missing / non-numeric / non-finite — JSON has no NaN.

    The one copy of a helper every control module used to define for itself
    (``_f`` / ``_finite``; issue #38).
    """
    if value is None:
        return None
    try:
        v = float(value)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def decimate_for_json(arr: Any, max_points: int = 500) -> list:
    """Thin an array to **at most** *max_points* values for JSON transport.

    Every ``ceil(n / max_points)``-th point is kept (the same stride
    ``export_results`` uses), and NaN/inf become None. Arrays of equal length
    decimate identically, so a Q array and its intensities stay aligned.
    The control modules' private copies used ``n // max_points``, which
    returned up to twice the cap.
    """
    a = np.asarray(arr, dtype=float)
    n = len(a)
    if max_points and max_points > 0 and n > max_points:
        a = a[::int(math.ceil(n / max_points))]
    return [v if math.isfinite(v) else None for v in a.tolist()]
