"""JSON schemas for pyirena.api.control tools in Anthropic tool-use format.

Each entry is a dict suitable for passing directly as a tool in the
Anthropic messages API (anthropic.types.Tool / the tools= list).

Usage
-----
from pyirena.api.control.schemas import TOOL_SCHEMAS
# pass to client.messages.create(tools=TOOL_SCHEMAS, ...)
"""
from __future__ import annotations

TOOL_SCHEMAS: list[dict] = [

    # -----------------------------------------------------------------------
    # Category A — Session lifecycle
    # -----------------------------------------------------------------------
    {
        "name": "open_dataset",
        "description": (
            "Load a NXcanSAS HDF5 file and create a fitting session. "
            "Returns a session_id you must pass to all other tools, plus a "
            "summary of the data (Q range, number of points)."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "file_path": {
                    "type": "string",
                    "description": "Absolute or relative path to the NXcanSAS HDF5 file.",
                },
            },
            "required": ["file_path"],
        },
    },
    {
        "name": "list_open_sessions",
        "description": "List all currently open fitting sessions.",
        "input_schema": {"type": "object", "properties": {}, "required": []},
    },
    {
        "name": "close_session",
        "description": "Close and discard a fitting session, freeing memory.",
        "input_schema": {
            "type": "object",
            "properties": {
                "session_id": {"type": "string", "description": "Session ID from open_dataset."},
            },
            "required": ["session_id"],
        },
    },
    {
        "name": "get_session_summary",
        "description": "Return a summary of the current session state: file, model, Q range, fit status.",
        "input_schema": {
            "type": "object",
            "properties": {
                "session_id": {"type": "string"},
            },
            "required": ["session_id"],
        },
    },

    # -----------------------------------------------------------------------
    # Category B — Model selection
    # -----------------------------------------------------------------------
    {
        "name": "list_available_models",
        "description": "List the fitting models available in this version.",
        "input_schema": {"type": "object", "properties": {}, "required": []},
    },
    {
        "name": "select_model",
        "description": (
            "Select and initialise a fitting model. For 'unified_fit', use "
            "nlevels to set the number of structural levels (1–5). "
            "Replaces any previously selected model and clears prior fit results. "
            "Returns the initial parameter table."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "session_id": {"type": "string"},
                "model_name": {
                    "type": "string",
                    "enum": ["unified_fit"],
                    "description": "Model to use.",
                    "default": "unified_fit",
                },
                "nlevels": {
                    "type": "integer",
                    "minimum": 1,
                    "maximum": 5,
                    "description": "Number of Unified Fit levels. Start with 1.",
                    "default": 1,
                },
            },
            "required": ["session_id"],
        },
    },
    {
        "name": "get_model_parameters",
        "description": (
            "Return the full parameter table for the current model. "
            "Each parameter has: name, value, fixed, lo (lower bound), "
            "hi (upper bound), units, description."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "session_id": {"type": "string"},
            },
            "required": ["session_id"],
        },
    },
    {
        "name": "get_model_description",
        "description": (
            "Return a text description of the Unified Fit model and its parameters, "
            "including physical meaning and fitting tips. "
            "Useful before starting a new fit."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "session_id": {"type": "string"},
            },
            "required": ["session_id"],
        },
    },

    # -----------------------------------------------------------------------
    # Category C — Parameter control
    # -----------------------------------------------------------------------
    {
        "name": "set_parameter_value",
        "description": (
            "Set the starting/current value of a named parameter. "
            "Use get_model_parameters() to see valid parameter names "
            "(e.g. 'Rg_1', 'G_2', 'background')."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "session_id":  {"type": "string"},
                "param_name":  {"type": "string", "description": "Parameter name, e.g. 'Rg_1'."},
                "value":       {"type": "number"},
            },
            "required": ["session_id", "param_name", "value"],
        },
    },
    {
        "name": "set_parameter_bounds",
        "description": "Set lower and upper bounds for a parameter during fitting.",
        "input_schema": {
            "type": "object",
            "properties": {
                "session_id": {"type": "string"},
                "param_name": {"type": "string"},
                "lo":         {"type": "number", "description": "Lower bound."},
                "hi":         {"type": "number", "description": "Upper bound."},
            },
            "required": ["session_id", "param_name", "lo", "hi"],
        },
    },
    {
        "name": "fix_parameter",
        "description": "Hold a parameter fixed at its current value during fitting.",
        "input_schema": {
            "type": "object",
            "properties": {
                "session_id": {"type": "string"},
                "param_name": {"type": "string"},
            },
            "required": ["session_id", "param_name"],
        },
    },
    {
        "name": "free_parameter",
        "description": "Allow a parameter to vary during fitting.",
        "input_schema": {
            "type": "object",
            "properties": {
                "session_id": {"type": "string"},
                "param_name": {"type": "string"},
            },
            "required": ["session_id", "param_name"],
        },
    },
    {
        "name": "fix_all_except",
        "description": (
            "Fix every parameter except those in free_list. "
            "Useful for staged fitting: e.g. fix all, then free only 'Rg_1' and 'G_1'."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "session_id": {"type": "string"},
                "free_list": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "Parameter names to leave free.",
                },
            },
            "required": ["session_id", "free_list"],
        },
    },
    {
        "name": "reset_parameters_to_defaults",
        "description": (
            "Reset all parameters to their factory defaults. "
            "Equivalent to calling select_model() again with the same number of levels."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "session_id": {"type": "string"},
            },
            "required": ["session_id"],
        },
    },

    # -----------------------------------------------------------------------
    # Category C' — Unified Fit level management
    # -----------------------------------------------------------------------
    {
        "name": "add_unified_level",
        "description": (
            "Add a new structural level to the Unified Fit model. "
            "position=-1 appends; position=1 inserts before the first level. "
            "Returns the updated parameter table."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "session_id": {"type": "string"},
                "position":   {
                    "type": "integer",
                    "description": "-1 = append (default); 1–n = insert before that level.",
                    "default": -1,
                },
            },
            "required": ["session_id"],
        },
    },
    {
        "name": "remove_unified_level",
        "description": (
            "Remove a structural level (1-based) from the Unified Fit model. "
            "Levels are renumbered from 1 after removal."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "session_id": {"type": "string"},
                "level":      {
                    "type": "integer",
                    "minimum": 1,
                    "description": "1-based level number to remove.",
                },
            },
            "required": ["session_id", "level"],
        },
    },

    # -----------------------------------------------------------------------
    # Category C''' — Per-level boolean options (correlations, mass_fractal, …)
    # -----------------------------------------------------------------------
    {
        "name": "get_level_options",
        "description": (
            "Return per-level boolean flag state: correlations, mass_fractal, "
            "link_B, link_RGCO.  These flags switch entire features of the "
            "intensity formula on or off and are SEPARATE from numeric "
            "parameters (use set_parameter_value for those, "
            "set_level_option for these).  Omit `level` to get all levels."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "session_id": {"type": "string"},
                "level": {
                    "type": ["integer", "null"],
                    "minimum": 1,
                    "description": "1-based level number; omit for all levels.",
                },
            },
            "required": ["session_id"],
        },
    },
    {
        "name": "set_level_option",
        "description": (
            "Toggle a per-level boolean flag.  Critical: setting numeric "
            "parameters ETA and PACK has NO EFFECT unless the level's "
            "`correlations` option is True — use this tool to enable it.  "
            "Other options: `mass_fractal` (auto-compute B), `link_B` (estimate "
            "B from Porod invariant), `link_RGCO` (link RgCO to previous "
            "level's Rg)."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "session_id": {"type": "string"},
                "level": {
                    "type": "integer",
                    "minimum": 1,
                    "description": "1-based level number.",
                },
                "option": {
                    "type": "string",
                    "enum": ["correlations", "mass_fractal", "link_B", "link_RGCO"],
                    "description": "Which boolean flag to toggle.",
                },
                "enabled": {
                    "type": "boolean",
                    "description": "True to turn on, False to turn off.",
                },
            },
            "required": ["session_id", "level", "option", "enabled"],
        },
    },

    # -----------------------------------------------------------------------
    # Category C'' — Q range
    # -----------------------------------------------------------------------
    {
        "name": "get_data_q_range",
        "description": "Return the Q range of the loaded dataset.",
        "input_schema": {
            "type": "object",
            "properties": {"session_id": {"type": "string"}},
            "required": ["session_id"],
        },
    },
    {
        "name": "get_fit_q_range",
        "description": "Return the Q range currently used for fitting (may be narrower than the data).",
        "input_schema": {
            "type": "object",
            "properties": {"session_id": {"type": "string"}},
            "required": ["session_id"],
        },
    },
    {
        "name": "set_fit_q_range",
        "description": (
            "Restrict the Q range used for fitting. "
            "Pass q_min and/or q_max; omit either to leave that end unchanged. "
            "Useful to exclude low-Q beam stop artefacts or high-Q noise."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "session_id": {"type": "string"},
                "q_min": {
                    "type": ["number", "null"],
                    "description": "Minimum Q to include in the fit.",
                },
                "q_max": {
                    "type": ["number", "null"],
                    "description": "Maximum Q to include in the fit.",
                },
            },
            "required": ["session_id"],
        },
    },
    {
        "name": "reset_fit_q_range",
        "description": "Restore the fit Q range to the full data range.",
        "input_schema": {
            "type": "object",
            "properties": {"session_id": {"type": "string"}},
            "required": ["session_id"],
        },
    },

    # -----------------------------------------------------------------------
    # Category C'''' — Local estimators (one-term fits over a Q sub-range)
    # -----------------------------------------------------------------------
    {
        "name": "fit_local_guinier",
        "description": (
            "Fit a single Guinier term I(q) = G·exp(-q²·Rg²/3) on the Q "
            "range [q_min, q_max] and return G and Rg.  Equivalent to the "
            "GUI's 'Fit Rg/G btwn cursors' button.  Useful for estimating "
            "good starting values for a level before run_fit.  Does NOT "
            "modify the model — call set_parameter_value() afterwards if "
            "you want to apply the result."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "session_id": {"type": "string"},
                "q_min":      {"type": "number", "description": "Lower Q [Å⁻¹]."},
                "q_max":      {"type": "number", "description": "Upper Q [Å⁻¹]."},
            },
            "required": ["session_id", "q_min", "q_max"],
        },
    },
    {
        "name": "fit_local_power_law",
        "description": (
            "Fit a single power-law term I(q) = B·q⁻ᴾ on the Q range "
            "[q_min, q_max] and return P and B.  Equivalent to the GUI's "
            "'Fit P/B btwn cursors' button.  Pick a window over the linear "
            "power-law region of the log-log plot (beyond the Guinier knee). "
            "Does NOT modify the model — call set_parameter_value() "
            "afterwards if you want to apply the result."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "session_id": {"type": "string"},
                "q_min":      {"type": "number", "description": "Lower Q [Å⁻¹]."},
                "q_max":      {"type": "number", "description": "Upper Q [Å⁻¹]."},
            },
            "required": ["session_id", "q_min", "q_max"],
        },
    },

    # -----------------------------------------------------------------------
    # Category D — Fit execution
    # -----------------------------------------------------------------------
    {
        "name": "run_fit",
        "description": (
            "Run the fitting algorithm. Uses the current parameter values as "
            "starting point and the current Q range restriction. "
            "Re-running uses the previous fitted values as the new starting point "
            "(re-calling after partial convergence is a valid strategy). "
            "Returns chi_squared, reduced_chi_squared, and the updated parameter values."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "session_id": {"type": "string"},
                "max_iter": {
                    "type": ["integer", "null"],
                    "description": "Maximum number of iterations. Default: 1000.",
                },
                "tolerance": {
                    "type": ["number", "null"],
                    "description": "Convergence tolerance. Default: scipy default.",
                },
                "random_seed": {
                    "type": ["integer", "null"],
                    "description": (
                        "Seed numpy RNG before the fit. Pass the same value to "
                        "reproduce a fit from an audit trail. Stored in the result "
                        "and in session.last_fit_result for stamping into audit JSON."
                    ),
                },
            },
            "required": ["session_id"],
        },
    },

    # -----------------------------------------------------------------------
    # Category E — Quality assessment
    # -----------------------------------------------------------------------
    {
        "name": "get_chi_squared",
        "description": "Return χ² and reduced χ² from the last fit.",
        "input_schema": {
            "type": "object",
            "properties": {"session_id": {"type": "string"}},
            "required": ["session_id"],
        },
    },
    {
        "name": "get_residuals",
        "description": (
            "Return normalised residuals and summary statistics from the last fit. "
            "rms close to 1 means the fit matches the error bars well."
        ),
        "input_schema": {
            "type": "object",
            "properties": {"session_id": {"type": "string"}},
            "required": ["session_id"],
        },
    },
    {
        "name": "get_fit_image",
        "description": (
            "Capture the current fit as a base64-encoded PNG. "
            "Works before a fit (shows data + model at current parameter values) "
            "and after a fit (adds a residuals subplot below). "
            "Examine the image to assess fit quality visually."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "session_id": {"type": "string"},
                "width":  {"type": "integer", "default": 1024, "description": "Image width in pixels."},
                "height": {"type": "integer", "default": 768,  "description": "Image height in pixels."},
            },
            "required": ["session_id"],
        },
    },
    {
        "name": "get_residuals_image",
        "description": (
            "Capture the fit + residuals plot as a base64-encoded PNG. "
            "Requires a completed fit."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "session_id": {"type": "string"},
                "width":  {"type": "integer", "default": 1024},
                "height": {"type": "integer", "default": 768},
            },
            "required": ["session_id"],
        },
    },
    {
        "name": "get_parameter_uncertainties",
        "description": (
            "Return parameter uncertainties from the last fit. "
            "Note: not available in Phase 1 — returns a placeholder."
        ),
        "input_schema": {
            "type": "object",
            "properties": {"session_id": {"type": "string"}},
            "required": ["session_id"],
        },
    },

    # -----------------------------------------------------------------------
    # Category F — Persistence
    # -----------------------------------------------------------------------
    {
        "name": "save_fit",
        "description": (
            "Save the fitted result back to NXcanSAS HDF5. "
            "Defaults to overwriting the original file. "
            "Pass output_path to save to a different location."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "session_id":  {"type": "string"},
                "output_path": {
                    "type": ["string", "null"],
                    "description": "Output file path. Defaults to the input file.",
                },
            },
            "required": ["session_id"],
        },
    },
    {
        "name": "export_fit_report",
        "description": (
            "Export a human-readable fit report. "
            "format='json' for machine-readable; format='markdown' for human review."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "session_id": {"type": "string"},
                "format": {
                    "type": "string",
                    "enum": ["json", "markdown"],
                    "default": "json",
                },
            },
            "required": ["session_id"],
        },
    },
]

# Convenience: look up a schema by name
TOOL_SCHEMA_BY_NAME: dict[str, dict] = {t["name"]: t for t in TOOL_SCHEMAS}
