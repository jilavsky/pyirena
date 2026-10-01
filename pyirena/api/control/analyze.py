"""``analyze`` — fit a curve with a saved configuration, in one call.

The coarse counterpart to the session tools. Those let an agent explore a
fit step by step; this is for the case where the answer is already known —
a scientist set the fit up once in the GUI, exported the parameters, and now
every new measurement should be fitted the same way.

    analyze(data={"q": [...], "intensity": [...], "error": [...]},
            config=json.load(open("modeling.json")))

That is one request instead of six, which matters over a synchronous
transport: each round trip costs latency, and a session left open by a client
that died costs memory until it is evicted. ``analyze`` opens its own session
and closes it in a ``finally``, so it cannot leak one.

**A config is more than model state.** The fitted Q range, the slit settings
and (for Simple Fits) which parameters are held fixed live alongside the
model in a config file, not inside ``model.to_dict()``. They are applied
here too — a replay that silently used the full Q range instead of the one
the scientist chose would return a plausible wrong number rather than an
error. ``pyirena.core.tool_config`` owns that translation and is the same
code ``pyirena.batch`` uses, so a fit replayed here and the same fit run from
the command line build an identical model.

Uncertainties are off by default and deliberately so: Monte-Carlo
uncertainties are the one setting that reliably outruns a network client's
timeout (a 10-population Modeling config at its exported ``n_mc_runs=50`` is
several minutes), so they must be asked for, and the number of runs is
reported back in case it was clamped.
"""
from __future__ import annotations

from typing import Any, Dict, Optional

from pyirena.api.control.errors import make_error
from pyirena.core.tool_config import (
    TOOL_SECTIONS,
    ConfigError,
    apply_prefits,
    build_setup,
)

#: Which control-API fit function drives each tool, and which config keys
#: (if any) feed its run-time arguments.
_RUNNERS = {
    "unified_fit":  ("run_fit", {"no_limits": "no_limits"}),
    "sizes":        ("run_sizes_fit", {}),
    "simple_fits":  ("run_simple_fit", {"no_limits": "no_limits"}),
    "modeling":     ("run_modeling_fit", {"fit_method": "fit_method",
                                          "no_limits": "no_limits"}),
    "waxs_peakfit": ("run_waxs_fit", {"weight_mode": "weight_mode",
                                      "no_limits": "no_limits"}),
    "carbon_fit":   ("run_carbon_fit", {"weighting": "weighting",
                                        "no_limits": "no_limits"}),
}


def analyze(
    data: Dict[str, Any],
    config: Dict[str, Any],
    tool: Optional[str] = None,
    include_arrays: bool = False,
    max_points: Optional[int] = 2000,
) -> dict:
    """Fit one curve with a saved pyIrena configuration and return the results.

    Parameters
    ----------
    data : dict
        ``{q, intensity, error?, dq?, label?, is_slit_smeared?, slit_length?}``
        — the same arguments :func:`open_dataset_from_data` takes. Q in Å⁻¹.
    config : dict
        Any configuration pyIrena writes, in any of its wrappers:

        * the GUI's *Export Parameters* sidecar,
          ``{"_pyirena_config": …, "<tool>": {…}}``;
        * the ``_pyirena_config`` attribute read out of a result ``.h5``,
          ``{"_pyirena_config": …, "state": {…}}``;
        * **an :func:`export_results` reply handed straight back** — the
          usual way to replay the fit you just ran on the next measurement,
          and the one that also carries the fitted Q range and slit settings;
        * ``{"tool": …, "model": {…}}``, or one tool's bare section with
          *tool* given.
    tool : str, optional
        Which tool to run. Inferred from the config when it holds exactly one
        tool section, which is what the GUI writes.
    include_arrays : bool
        Also return the model curve, residuals and per-component curves.
    max_points : int or None
        Cap on each returned array; None for no cap.

    Returns
    -------
    dict
        The :func:`export_results` payload, plus an ``analyze`` block naming
        the tool, the Q range actually fitted and anything that was adjusted.
        Errors are the usual ``{error, code, suggestion}`` dicts — including
        ``BAD_CONFIG`` when the configuration cannot be turned into a model.
    """
    from pyirena.api.control import close_session, open_dataset_from_data
    from pyirena.api.control.export import export_results
    from pyirena.api.control.session import get_session

    if not isinstance(data, dict):
        return make_error(
            f"'data' must be an object with q and intensity, got {type(data).__name__}.",
            suggestion='Pass {"q": [...], "intensity": [...], "error": [...]}.',
            code="BAD_ARGUMENTS",
        )
    if not isinstance(config, dict):
        return make_error(
            f"'config' must be an object, got {type(config).__name__}.",
            suggestion="Pass the JSON the GUI's Export Parameters wrote.",
            code="BAD_CONFIG",
        )

    # Resolve the config first: a bad config should cost nothing and should
    # not leave a session behind.
    try:
        setup = build_setup(
            config, tool,
            data_is_slit_smeared=bool(data.get("is_slit_smeared", False)),
            data_slit_length=float(data.get("slit_length", 0.0) or 0.0),
        )
    except ConfigError as exc:
        return make_error(
            str(exc),
            suggestion=f"Expected a config for one of: {', '.join(TOOL_SECTIONS)}.",
            code="BAD_CONFIG",
        )

    opened = open_dataset_from_data(
        q=data.get("q"),
        intensity=data.get("intensity"),
        error=data.get("error"),
        dq=data.get("dq"),
        label=data.get("label", ""),
        is_slit_smeared=bool(data.get("is_slit_smeared", False)),
        slit_length=float(data.get("slit_length", 0.0) or 0.0),
        error_fraction=float(data.get("error_fraction", 0.05)),
    )
    if "error" in opened:
        return opened

    session_id = opened["session_id"]
    notes: list = list(setup.warnings)
    try:
        s = get_session(session_id)
        s.model = setup.model
        s.model_name = setup.tool

        # The fit range is part of the configuration, not a detail: a replay
        # over the full curve is a different fit, and a silent one.
        applied_q = _apply_q_range(s, setup, notes)

        # Simple Fits carries its held-fixed parameters separately.
        if setup.tool == "simple_fits" and setup.fixed_params:
            s._simple_fixed = set(setup.fixed_params)
            for name, value in setup.fixed_params.items():
                if name in s.model.params:
                    s.model.params[name] = value

        # Pre-fit steps the config asks for (Sizes background windows, the
        # Simple Fits background the Invariant integrates on top of, WAXS
        # peak re-centring). Which data each sees is not the same — see
        # apply_prefits — so hand it both the full curve and the fitted range.
        from pyirena.api.control.session import fit_mask
        mask = fit_mask(s)
        notes.extend(apply_prefits(setup, setup.section,
                                   s.q, s.intensity,
                                   s.q[mask], s.intensity[mask]))

        fit_result = _run_fit(setup, session_id)
        if isinstance(fit_result, dict) and "error" in fit_result:
            return fit_result

        report = export_results(session_id, include_arrays=include_arrays,
                                max_points=max_points)
        if isinstance(report, dict) and "error" in report:
            return report

        report["analyze"] = {
            "tool": setup.tool,
            "config_applied": True,
            "fit_q_min": applied_q[0],
            "fit_q_max": applied_q[1],
            "fixed_parameters": sorted(setup.fixed_params),
            "n_points_input": opened["summary"]["n_points"],
            "cleaning": opened["summary"]["cleaning"],
            "notes": notes,
        }
        return report
    finally:
        # One call, one session, no leak — whatever happened above.
        close_session(session_id)


def _apply_q_range(session, setup, notes: list):
    """Apply the config's Q range, reporting rather than silently ignoring it."""
    import numpy as np

    q_min, q_max = setup.fit_q_min, setup.fit_q_max
    if q_min is None and q_max is None:
        return None, None

    data_min, data_max = float(np.min(session.q)), float(np.max(session.q))
    lo = data_min if q_min is None else max(q_min, data_min)
    hi = data_max if q_max is None else min(q_max, data_max)

    if lo >= hi or not np.any((session.q >= lo) & (session.q <= hi)):
        # A one-sided range is legitimate, so either bound may be absent —
        # say so rather than formatting None into the note that explains it.
        def _bound(value):
            return "open" if value is None else f"{value:.4g}"

        notes.append(
            f"The config's Q range [{_bound(q_min)}, {_bound(q_max)}] does not "
            f"overlap this curve [{data_min:.4g}, {data_max:.4g}]; fitting the "
            "full range instead."
        )
        return None, None

    # Only say "clipped" when it means something. A config's Q range is
    # usually the data's own limits written back out, so it lands a few float
    # ulps outside them; reporting that as a clip trains the reader to ignore
    # the note that matters.
    def _materially_inside(bound, edge):
        return bound is not None and abs(bound - edge) > 1e-9 * max(abs(edge), 1e-30)

    if _materially_inside(q_min, lo) or _materially_inside(q_max, hi):
        notes.append(
            f"The config's Q range was clipped to this curve: "
            f"[{lo:.4g}, {hi:.4g}]."
        )
    session.fit_q_min, session.fit_q_max = lo, hi
    return lo, hi


def _run_fit(setup, session_id: str):
    """Call the tool's own run_* function, passing what the config asks for.

    **Unified Fit keeps ``walk_limits`` on here, on purpose.** It looks wrong
    for a tool whose job is to replay a setup exactly — the bounds travel
    correctly and are then widened when the fit reaches them — and it is the
    behaviour the tool exists for. A saved config is replayed across a series
    of measurements in which the structure is *growing*, so a limit that was
    right for the first scan is passed by the tenth; a fit that stopped dead
    at it would report the bound rather than the size. The limits cannot
    simply be dropped instead, because several fitting methods require finite
    bounds. Walking them is what lets one configuration follow a sample
    through its whole run.

    Pass ``no_limits`` in the config for the other case — a bound you want
    ignored for this fit rather than followed.
    """
    from pyirena.api import control as ctrl

    name, arg_map = _RUNNERS[setup.tool]
    runner = getattr(ctrl, name)
    section = setup.section if isinstance(setup.section, dict) else {}
    kwargs = {}
    for arg, config_key in arg_map.items():
        value = section.get(config_key)
        if value is None:
            # Carbon keeps its weighting inside the model block, where model
            # construction has already restored it; the runner's default
            # would overwrite that with 'auto' and change the answer. Only
            # the model's own value is a fallback — never an invented one.
            value = getattr(setup.model, arg, None) if arg == "weighting" else None
        if value is not None:
            kwargs[arg] = value
    try:
        return runner(session_id, **kwargs)
    except TypeError:
        # A config naming an argument this version does not take should not
        # be fatal — run the fit with the defaults and say so.
        return runner(session_id)
