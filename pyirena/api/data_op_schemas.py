"""JSON schemas for pyirena.api.data_ops in Anthropic tool-use format.

Mirrors pyirena/api/calculator_schemas.py. Kept as its own registry so the
control-surface schema count (TOOL_SCHEMA_BY_NAME) stays untouched;
pyirena/mcp/dispatch.py merges all registries to build the dispatcher.

These are the only pyirena tools that WRITE data files, so every
description says what gets written and where.

Usage
-----
from pyirena.api.data_op_schemas import DATA_OP_TOOL_SCHEMAS
# pass to client.messages.create(tools=DATA_OP_TOOL_SCHEMAS, ...)
"""
from __future__ import annotations

from pyirena.api.data_ops import REBIN_MODES, SIMILARITY_REFERENCES

_OUTPUT_FOLDER_PROPERTY = {
    "type": ["string", "null"],
    "description": (
        "Where to write the result. Defaults to a sibling of the source "
        "folder with '_manip' appended (e.g. /data/run42 -> "
        "/data/run42_manip). Pass an explicit folder when the default is "
        "refused as outside PYIRENA_DATA_ROOT."
    ),
}

_MERGE_OUTPUT_FOLDER_PROPERTY = {
    "type": ["string", "null"],
    "description": (
        "Where to write the merged file. Defaults to a sibling of file1's "
        "folder with '_merged' appended (e.g. /data/usaxs -> "
        "/data/usaxs_merged)."
    ),
}

_OVERWRITE_NOTE = (
    "Repeating the same operation on the same input overwrites the previous "
    "result; different operations use different filename suffixes and do not "
    "collide."
)


DATA_OP_TOOL_SCHEMAS: list[dict] = [

    # -----------------------------------------------------------------------
    # Many-dataset
    # -----------------------------------------------------------------------
    {
        "name": "average_data",
        "description": (
            "Average two or more SAS datasets and write the result as "
            "<first file stem>_avg.h5. All frames are interpolated onto the "
            "first file's Q grid; the uncertainty is the larger of the "
            "propagated error and the point-to-point standard deviation. "
            "Set similarity_check=true to screen for radiation damage — a "
            "cormap test drops frames whose shape has changed, and the "
            "discarded ones are listed in 'rejected'. Datasets must all have "
            "the same slit-smearing status. " + _OVERWRITE_NOTE
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "files": {
                    "type": "array",
                    "items": {"type": "string"},
                    "minItems": 2,
                    "description": "Paths of the datasets to average (at least 2).",
                },
                "output_folder": _OUTPUT_FOLDER_PROPERTY,
                "similarity_check": {
                    "type": "boolean",
                    "description": (
                        "Screen frames with a cormap similarity test and drop "
                        "outliers before averaging. Use for a time series "
                        "where radiation damage is possible."
                    ),
                    "default": False,
                },
                "similarity_p_min": {
                    "type": "number",
                    "description": (
                        "P-value threshold; frames below it are discarded. "
                        "Typical range 0.001-0.05."
                    ),
                    "default": 0.01,
                },
                "similarity_method": {
                    "type": "string",
                    "description": "Similarity algorithm. Currently 'cormap'.",
                    "default": "cormap",
                },
                "similarity_reference": {
                    "type": "string",
                    "enum": list(SIMILARITY_REFERENCES),
                    "description": (
                        "'first' compares every frame with frame 0 (always "
                        "kept); 'majority' compares with the median of all "
                        "frames."
                    ),
                    "default": "first",
                },
                "similarity_normalize_scale": {
                    "type": "boolean",
                    "description": (
                        "Rescale each frame to the reference before comparing "
                        "so flux drift is not mistaken for a shape change."
                    ),
                    "default": True,
                },
            },
            "required": ["files"],
        },
    },

    # -----------------------------------------------------------------------
    # Two-dataset
    # -----------------------------------------------------------------------
    {
        "name": "subtract_data",
        "description": (
            "Subtract a buffer/background dataset from a sample and write "
            "<sample stem>_sub.h5. Computes I_sample - buffer_scale * "
            "I_buffer, interpolating the buffer onto the sample's Q grid. "
            "Check 'n_dropped_nonpositive' in the result: points that go to "
            "zero or negative are stripped on save, which usually means "
            "over-subtraction — lower buffer_scale. The two datasets must "
            "have matching slit smearing. " + _OVERWRITE_NOTE
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "sample_file": {"type": "string", "description": "The sample dataset."},
                "buffer_file": {
                    "type": "string",
                    "description": (
                        "The buffer, solvent or background dataset to subtract."
                    ),
                },
                "buffer_scale": {
                    "type": "number",
                    "description": "Multiplier applied to the buffer before subtracting.",
                    "default": 1.0,
                },
                "auto_scale": {
                    "type": "boolean",
                    "description": (
                        "Fit buffer_scale from the intensity integral ratio "
                        "over a Q window. Requires BOTH auto_q_min and "
                        "auto_q_max."
                    ),
                    "default": False,
                },
                "auto_q_min": {
                    "type": ["number", "null"],
                    "description": (
                        "Lower bound (1/A) of the window where sample and "
                        "buffer should match. Required when auto_scale is on."
                    ),
                },
                "auto_q_max": {
                    "type": ["number", "null"],
                    "description": (
                        "Upper bound (1/A) of the auto-scale window. Required "
                        "when auto_scale is on."
                    ),
                },
                "output_folder": _OUTPUT_FOLDER_PROPERTY,
            },
            "required": ["sample_file", "buffer_file"],
        },
    },

    {
        "name": "divide_data",
        "description": (
            "Divide one dataset by another and write <numerator stem>_div.h5. "
            "Computes I_num / (denominator_scale * I_den - "
            "denominator_background), interpolating the denominator onto the "
            "numerator's Q grid. Points where the denominator is zero come "
            "back non-finite and are reported in 'n_nonfinite', then stripped "
            "on save. Slit smearing must match. " + _OVERWRITE_NOTE
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "numerator_file": {"type": "string", "description": "Dataset to divide."},
                "denominator_file": {"type": "string", "description": "Dataset to divide by."},
                "denominator_scale": {
                    "type": "number",
                    "description": "Multiplier applied to the denominator before dividing.",
                    "default": 1.0,
                },
                "denominator_background": {
                    "type": "number",
                    "description": (
                        "Flat background subtracted from the scaled denominator, "
                        "in the file's intensity units (1/cm when absolute)."
                    ),
                    "default": 0.0,
                },
                "output_folder": _OUTPUT_FOLDER_PROPERTY,
            },
            "required": ["numerator_file", "denominator_file"],
        },
    },

    # -----------------------------------------------------------------------
    # Single-dataset
    # -----------------------------------------------------------------------
    {
        "name": "scale_data",
        "description": (
            "Scale a dataset's intensity and/or subtract a flat background, "
            "writing <stem>_scaled.h5. Computes scale_I * I - background — "
            "the background is subtracted AFTER scaling. Use to put a dataset "
            "on absolute scale, or to remove a known constant background. "
            + _OVERWRITE_NOTE
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "file": {"type": "string", "description": "Input dataset."},
                "scale_I": {
                    "type": "number",
                    "description": "Multiplicative intensity factor.",
                    "default": 1.0,
                },
                "background": {
                    "type": "number",
                    "description": (
                        "Flat background subtracted after scaling, in the "
                        "file's intensity units (1/cm when absolute)."
                    ),
                    "default": 0.0,
                },
                "scale_uncertainty": {
                    "type": ["number", "null"],
                    "description": (
                        "Factor applied to the uncertainties. Defaults to "
                        "scale_I."
                    ),
                },
                "output_folder": _OUTPUT_FOLDER_PROPERTY,
            },
            "required": ["file"],
        },
    },

    {
        "name": "trim_data",
        "description": (
            "Keep only the points inside a Q window, writing <stem>_trimmed.h5. "
            "Use to drop a noisy high-Q tail or a beamstop-contaminated low-Q "
            "region before fitting. Both bounds are inclusive. "
            + _OVERWRITE_NOTE
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "file": {"type": "string", "description": "Input dataset."},
                "q_min": {
                    "type": "number",
                    "description": "Lower Q bound in 1/A, inclusive.",
                    "default": 0.0,
                },
                "q_max": {
                    "type": ["number", "null"],
                    "description": (
                        "Upper Q bound in 1/A, inclusive. Null means no upper "
                        "limit."
                    ),
                },
                "output_folder": _OUTPUT_FOLDER_PROPERTY,
            },
            "required": ["file"],
        },
    },

    {
        "name": "rebin_data",
        "description": (
            "Resample a dataset onto a new Q grid, writing <stem>_rebinned.h5. "
            "Use to thin a very dense curve before fitting, or to put two "
            "datasets on a common grid. Interpolation is linear in log(I) vs "
            "log(Q) and does NOT extrapolate, so the result can have fewer "
            "points than requested. " + _OVERWRITE_NOTE
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "file": {"type": "string", "description": "Input dataset."},
                "mode": {
                    "type": "string",
                    "enum": list(REBIN_MODES),
                    "description": (
                        "'log' for geometric spacing (the usual choice for "
                        "SAS), 'linear' for even spacing, 'reference' to reuse "
                        "another file's Q grid."
                    ),
                    "default": "log",
                },
                "n_points": {
                    "type": "integer",
                    "minimum": 2,
                    "description": "Target number of points for 'log' and 'linear'.",
                    "default": 200,
                },
                "q_min": {
                    "type": ["number", "null"],
                    "description": (
                        "Grid lower bound in 1/A; defaults to the data's own "
                        "minimum. Must be > 0 for mode='log'."
                    ),
                },
                "q_max": {
                    "type": ["number", "null"],
                    "description": "Grid upper bound in 1/A; defaults to the data's maximum.",
                },
                "reference_file": {
                    "type": ["string", "null"],
                    "description": (
                        "Dataset whose Q grid to reuse. Required for "
                        "mode='reference'."
                    ),
                },
                "output_folder": _OUTPUT_FOLDER_PROPERTY,
            },
            "required": ["file"],
        },
    },

    # -----------------------------------------------------------------------
    # Merge
    # -----------------------------------------------------------------------
    {
        "name": "merge_datasets",
        "description": (
            "Merge two datasets that overlap in Q into one curve, writing "
            "<file1 stem>_merged.h5. file1 must be the LOWER-Q dataset "
            "(typically USAXS) and is the absolute-intensity reference; file2 "
            "is the higher-Q dataset (typically SAXS) brought onto it. A "
            "scale factor is fitted in the overlap region. Note a background "
            "is ALWAYS fitted and subtracted from file1 and cannot be forced "
            "to zero — check the returned 'background' is small relative to "
            "your intensities. Mixing slit-smeared USAXS with pinhole SAXS is "
            "normal and allowed. " + _OVERWRITE_NOTE
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "file1": {
                    "type": "string",
                    "description": (
                        "Lower-Q dataset (usually USAXS). Keeps its intensity "
                        "scale; the fitted background is subtracted from it."
                    ),
                },
                "file2": {
                    "type": "string",
                    "description": (
                        "Higher-Q dataset (usually SAXS), rescaled onto file1."
                    ),
                },
                "q_overlap_min": {
                    "type": ["number", "null"],
                    "description": (
                        "Lower bound (1/A) of the overlap region used for the "
                        "fit. Omit both bounds to auto-detect."
                    ),
                },
                "q_overlap_max": {
                    "type": ["number", "null"],
                    "description": "Upper bound (1/A) of the overlap region.",
                },
                "fit_scale": {
                    "type": "boolean",
                    "description": (
                        "Fit the scale factor. When false, fixed_scale_value "
                        "is used instead."
                    ),
                    "default": True,
                },
                "scale_dataset": {
                    "type": "integer",
                    "enum": [1, 2],
                    "description": (
                        "Which dataset is rescaled. 2 (default) keeps file1 "
                        "absolute and is far faster — it has a closed-form "
                        "solution, while 1 falls back to iterative fitting."
                    ),
                    "default": 2,
                },
                "fixed_scale_value": {
                    "type": "number",
                    "description": "Scale used when fit_scale is false.",
                    "default": 1.0,
                },
                "fit_qshift": {
                    "type": "boolean",
                    "description": (
                        "Also fit a Q offset between the datasets. Usually "
                        "unnecessary; leave off unless the overlap clearly "
                        "misaligns in Q."
                    ),
                    "default": False,
                },
                "fixed_qshift_value": {
                    "type": "number",
                    "description": "Q shift in 1/A used when fit_qshift is false.",
                    "default": 0.0,
                },
                "qshift_dataset": {
                    "type": "integer",
                    "enum": [0, 1, 2],
                    "description": "Which dataset the Q shift applies to; 0 means none.",
                    "default": 0,
                },
                "split_at_left_cursor": {
                    "type": "boolean",
                    "description": (
                        "True: hard split at q_overlap_min so no duplicate Q "
                        "values remain. False (default): keep both datasets' "
                        "points across the overlap."
                    ),
                    "default": False,
                },
                "output_folder": _MERGE_OUTPUT_FOLDER_PROPERTY,
            },
            "required": ["file1", "file2"],
        },
    },

    {
        "name": "match_merge_files",
        "description": (
            "Pair up files from two folders for batch merging. Reads only — "
            "writes nothing. Matches on the text before the first underscore "
            "plus the last integer in the filename, so sampleA_usaxs_007.h5 "
            "pairs with sampleA_saxs_007.h5. Call this before looping "
            "merge_datasets over a run, and check 'unmatched_1' / "
            "'unmatched_2' for files with no partner."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "folder1": {
                    "type": "string",
                    "description": "Folder of lower-Q datasets (usually USAXS).",
                },
                "folder2": {
                    "type": "string",
                    "description": "Folder of higher-Q datasets (usually SAXS).",
                },
            },
            "required": ["folder1", "folder2"],
        },
    },
]

# Convenience: look up a schema by name
DATA_OP_SCHEMA_BY_NAME: dict[str, dict] = {
    t["name"]: t for t in DATA_OP_TOOL_SCHEMAS
}
