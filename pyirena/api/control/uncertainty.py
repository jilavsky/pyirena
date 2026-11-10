"""Parameter uncertainties on the agent surface (issues #28 and #38).

Monte-Carlo uncertainties used to exist only in the GUI and in
``pyirena.batch``; an agent driving the same fit over MCP or ZMQ got none, and
``get_parameter_uncertainties`` was a placeholder. The ``run_*_fit`` control
functions for Unified Fit, Size Distribution, Simple Fits and Modeling now
take ``n_mc_runs`` (0 = off, the default) and call the same core
implementation the GUI and batch use, so the three surfaces return the same
numbers for the same setup and seed. Carbon already had ``n_mc_runs``; WAXS
Peak Fit reports ``curve_fit`` covariance uncertainties, which every fit
computes.

**Off by default, deliberately.** Monte-Carlo is the one setting that reliably
outruns a network client's timeout (a ZMQ request has a 55 s budget), so it
must be asked for. Every result says how many refits succeeded.

What is stored on the session, and returned by
:func:`get_parameter_uncertainties` and in the ``export_results`` envelope, is
one JSON-safe block::

    {available, reason, method, n_runs, n_ok, seed,
     parameters: [{name, std}, ...],
     scalars: {...}            # Size Distribution: rg_std, volume_fraction_std, …
     distribution_std: [...]}  # Size Distribution, with include_arrays only
"""
from __future__ import annotations

import math
from typing import Optional

from pyirena.api.control.errors import make_error, no_fit, no_session
from pyirena.api.control.session import get_session

#: Methods reported in the ``method`` field.
MONTE_CARLO = "monte_carlo"
COVARIANCE = "covariance"


def _f(value):
    """Float, or None for missing/non-finite — JSON has no NaN."""
    if value is None:
        return None
    try:
        v = float(value)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def check_n_mc_runs(n_mc_runs) -> tuple[int, Optional[dict]]:
    """Validate the ``n_mc_runs`` argument: 0 (off) or at least 2."""
    try:
        n = int(n_mc_runs or 0)
    except (TypeError, ValueError):
        return 0, make_error(
            f"n_mc_runs must be an integer, got {n_mc_runs!r}.",
            suggestion="Pass 0 for no uncertainties, or e.g. 20.",
            code="BAD_ARGUMENTS",
        )
    if n < 0 or n == 1:
        return 0, make_error(
            f"n_mc_runs must be 0 (off) or at least 2, got {n}.",
            suggestion="A standard deviation needs at least two refits; 10-50 is usual.",
            code="BAD_ARGUMENTS",
        )
    return n, None


def summarize_mc(mc: dict, parameters: list[dict], scalars: Optional[dict] = None,
                 distribution_std=None) -> dict:
    """Turn a core ``calculate_uncertainty_mc`` result into the stored block."""
    block = {
        "available": bool(mc.get("available")),
        "reason": mc.get("reason"),
        "method": MONTE_CARLO,
        "n_runs": int(mc.get("n_runs", 0)),
        "n_ok": int(mc.get("n_ok", 0)),
        "seed": mc.get("seed"),
        "parameters": parameters if mc.get("available") else [],
    }
    if scalars is not None and mc.get("available"):
        block["scalars"] = {k: _f(v) for k, v in scalars.items()}
    if distribution_std is not None and mc.get("available"):
        block["_distribution_std"] = [_f(v) for v in distribution_std]
    return block


def mc_reply(block: dict) -> dict:
    """The compact form a ``run_*_fit`` reply carries."""
    out = {k: block[k] for k in ("available", "reason", "method", "n_runs", "n_ok", "seed")}
    if "scalars" in block:
        out["scalars"] = block["scalars"]
    return out


def _public(block: dict, include_arrays: bool) -> dict:
    out = {k: v for k, v in block.items() if not k.startswith("_")}
    if include_arrays and "_distribution_std" in block:
        out["distribution_std"] = block["_distribution_std"]
    return out


def _not_computed(tool: str) -> dict:
    return {
        "available": False,
        "method": None,
        "reason": (
            "No uncertainties were computed for the last fit. Re-run "
            f"{_RUNNER.get(tool, 'the fit')} with n_mc_runs > 0 (e.g. 20) for "
            "Monte-Carlo uncertainties."
        ),
        "parameters": [],
    }


_RUNNER = {
    "unified_fit": "run_fit",
    "sizes": "run_sizes_fit",
    "simple_fits": "run_simple_fit",
    "modeling": "run_modeling_fit",
    "carbon_fit": "run_carbon_fit",
}


def _covariance_simple(s) -> dict:
    stds = (s.last_fit_result or {}).get("params_std") or {}
    return {
        "available": any(_f(v) is not None for v in stds.values()),
        "reason": None if stds else "the fit returned no covariance",
        "method": COVARIANCE,
        "parameters": [{"name": k, "std": _f(v)} for k, v in stds.items()],
    }


def _covariance_waxs(s) -> dict:
    r = s.last_fit_result or {}
    rows = [{"name": f"background {k}", "std": _f(v)}
            for k, v in (r.get("bg_params_std") or {}).items()]
    for i, peak_std in enumerate(r.get("peaks_std") or []):
        for k, v in (peak_std or {}).items():
            rows.append({"name": f"peak {i + 1} {k}", "std": _f(v)})
    return {
        "available": any(row["std"] is not None for row in rows),
        "reason": None if rows else "the fit returned no covariance",
        "method": COVARIANCE,
        "parameters": rows,
    }


def _carbon(s) -> dict:
    r = s.last_fit_result
    errors = getattr(r, "errors", None) or {}
    n_mc = int(getattr(s.model, "n_mc_runs", 0) or 0)
    return {
        "available": any(_f(v) is not None for v in errors.values()),
        "reason": None if errors else "the fit returned no uncertainties",
        "method": MONTE_CARLO if n_mc > 1 else COVARIANCE,
        "n_runs": n_mc or None,
        "parameters": [{"name": k, "std": _f(v)} for k, v in errors.items()],
    }


def get_parameter_uncertainties(session_id: str, include_arrays: bool = False) -> dict:
    """Parameter uncertainties from the last fit, for any fitting tool.

    Monte-Carlo uncertainties are returned when the last ``run_*_fit`` was
    called with ``n_mc_runs > 0``. Without them, Simple Fits and WAXS Peak Fit
    fall back to the covariance σ every fit computes, and the other tools
    report ``available: false`` with how to get them.

    Args:
        session_id: The session.
        include_arrays: Size Distribution only — also return the per-bin
            ``distribution_std``.

    Returns:
        ``{available, reason, method, parameters: [{name, std}], …}`` with
        ``n_runs``/``n_ok``/``seed`` for Monte-Carlo results.
    """
    s = get_session(session_id)
    if s is None:
        return no_session(session_id)
    if s.last_fit_result is None:
        return no_fit(session_id)

    tool = s.model_name
    block = getattr(s, "last_uncertainty", None)
    if block is not None:
        out = _public(block, include_arrays)
    elif tool == "simple_fits":
        out = _covariance_simple(s)
    elif tool == "waxs_peakfit":
        out = _covariance_waxs(s)
    elif tool == "carbon_fit":
        out = _carbon(s)
    else:
        out = _not_computed(tool)
    out = dict(out)
    out.update({"ok": True, "session_id": session_id, "tool": tool})
    return out


def export_block(s) -> Optional[dict]:
    """The uncertainty block for the ``export_results`` envelope (no arrays)."""
    if s.last_fit_result is None:
        return None
    out = get_parameter_uncertainties(s.session_id)
    if "error" in out:
        return None
    return {k: v for k, v in out.items() if k not in ("ok", "session_id")}
