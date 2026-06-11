"""pyirena.api.control — bidirectional fitting control surface.

Version
-------
``__version__`` mirrors ``pyirena.__version__`` and is intended to be stamped
into audit JSON produced by pyirena-ai so that replays can verify the control
API version that was used.

    >>> import pyirena.api.control as ctrl
    >>> ctrl.__version__
    '0.8.2'


Extends the read-only pyirena.api facade with tools that allow an LLM agent
(or any automation script) to configure models, run fits, and evaluate results.

Quick start
-----------
>>> import pyirena.api.control as ctrl
>>> r = ctrl.open_dataset("/data/scan_007.h5")
>>> sid = r["session_id"]
>>> ctrl.select_model(sid, model_name="unified_fit", nlevels=1)
>>> ctrl.fix_all_except(sid, ["Rg_1", "G_1", "background"])
>>> ctrl.run_fit(sid)
>>> print(ctrl.get_chi_squared(sid))
>>> ctrl.get_fit_image(sid)   # returns {"image_base64": "...", ...}
>>> ctrl.save_fit(sid)

All functions return plain dicts.  Errors are {"error": ..., "code": ...,
"suggestion": ...} dicts rather than exceptions.

JSON schemas for LLM tool-use
------------------------------
>>> from pyirena.api.control.schemas import TOOL_SCHEMAS
>>> # pass to client.messages.create(tools=TOOL_SCHEMAS, ...)
"""
from __future__ import annotations

from pyirena import __version__  # re-exported for audit-trail stamping

from pyirena.api.control.unified_fit import (
    # A — session lifecycle
    open_dataset,
    list_open_sessions,
    close_session,
    get_session_summary,
    # B — model
    list_available_models,
    select_model,
    get_model_parameters,
    get_model_description,
    # C — parameter control
    set_parameter_value,
    set_parameter_bounds,
    fix_parameter,
    free_parameter,
    fix_all_except,
    reset_parameters_to_defaults,
    # C' — level management
    add_unified_level,
    remove_unified_level,
    # C''' — per-level boolean options (correlations, mass_fractal, …)
    get_level_options,
    set_level_option,
    # C'' — Q range
    get_data_q_range,
    get_fit_q_range,
    set_fit_q_range,
    reset_fit_q_range,
    # D — fit execution
    run_fit,
    # E — quality
    get_chi_squared,
    get_residuals,
    get_fit_image,
    get_residuals_image,
    get_parameter_uncertainties,
    # F — persistence
    save_fit,
    export_fit_report,
)

__all__ = [
    "__version__",
    # A
    "open_dataset", "list_open_sessions", "close_session", "get_session_summary",
    # B
    "list_available_models", "select_model", "get_model_parameters", "get_model_description",
    # C
    "set_parameter_value", "set_parameter_bounds",
    "fix_parameter", "free_parameter", "fix_all_except",
    "reset_parameters_to_defaults",
    # C'
    "add_unified_level", "remove_unified_level",
    # C'''
    "get_level_options", "set_level_option",
    # C''
    "get_data_q_range", "get_fit_q_range", "set_fit_q_range", "reset_fit_q_range",
    # D
    "run_fit",
    # E
    "get_chi_squared", "get_residuals",
    "get_fit_image", "get_residuals_image",
    "get_parameter_uncertainties",
    # F
    "save_fit", "export_fit_report",
]
