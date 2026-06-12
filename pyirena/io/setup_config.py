"""Embed the full pyIrena GUI setup for one tool inside a NXcanSAS results group.

Background
----------
NXcanSAS files written by pyirena store *fit results* (Q, I_data, I_model,
fitted values, …) but historically have not stored the *setup* — the controls,
fit flags, bounds, link flags, cursor positions, etc. that the user (or the
pyirena-ai agent) chose before running the fit. Without that, a downstream
user who opens the result has no way to "pick up where the agent left off"
in the GUI.

This module provides a tiny, uniform way to round-trip the full setup
without inventing one NXcanSAS sub-group per control: the entire state dict
(the same dict the panel's ``get_current_state()`` / ``apply_state()``
methods already use) is JSON-serialised and stored as a single string
attribute (``_pyirena_config``) on the tool's results group. The envelope
mirrors the existing ``_pyirena_config`` header that ``ModelingPanel.export_json``
writes to sidecar JSON files, so the two channels (sidecar file, embedded
attribute) share one shape.

Envelope shape on disk
----------------------
``_pyirena_config`` is a JSON object of the form::

    {
      "_pyirena_config": {
        "tool":           "unified_fit",
        "schema_version": 1,
        "saved_by":       "pyirena <__version__>",
        "saved_at":       "<ISO-8601 timestamp>"
      },
      "state": { ... the tool state dict as understood by apply_state() ... }
    }

Schema versioning is delegated to ``StateManager._migrate_state()`` — there is
no second migration table here. The header's ``schema_version`` is copied
from ``state["schema_version"]`` if present (most tools have one); panels
that don't version their state (``unified_fit``) omit the field.
"""
from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Optional

import h5py


PYIRENA_CONFIG_ATTR = "_pyirena_config"


class SetupConfigError(Exception):
    """Base class for setup-config read errors."""


class SetupConfigToolMismatch(SetupConfigError):
    """Raised when the on-disk attribute belongs to a different tool."""

    def __init__(self, expected: str, actual: str, file_path: str):
        super().__init__(
            f"Setup in '{file_path}' belongs to tool '{actual}', "
            f"not '{expected}'."
        )
        self.expected = expected
        self.actual = actual
        self.file_path = file_path


def _build_envelope(tool: str, state: Dict[str, Any]) -> Dict[str, Any]:
    try:
        from pyirena import __version__ as pyirena_version
    except Exception:
        pyirena_version = "unknown"

    header: Dict[str, Any] = {
        "tool":     tool,
        "saved_by": f"pyirena {pyirena_version}",
        "saved_at": datetime.now().isoformat(),
    }
    if isinstance(state, dict) and "schema_version" in state:
        header["schema_version"] = state["schema_version"]

    return {
        PYIRENA_CONFIG_ATTR: header,
        "state":             state,
    }


def write_setup_config(group: h5py.Group, tool: str, state: Dict[str, Any]) -> None:
    """Serialise ``state`` and store it as the ``_pyirena_config`` attribute.

    Parameters
    ----------
    group
        An open ``h5py.Group`` for the tool's results (e.g. the group at
        ``/entry/unified_fit_results``).  The writer overwrites any existing
        attribute with the same name.
    tool
        Short tool identifier matching the state-manager section name
        (e.g. ``"unified_fit"``, ``"sizes"``, ``"modeling"``,
        ``"simple_fits"``, ``"waxs_peakfit"``).
    state
        The full state dict (same shape the panel's ``apply_state`` expects).
        Values must be JSON-serialisable (no numpy arrays / objects).
    """
    envelope = _build_envelope(tool, state)
    try:
        payload = json.dumps(envelope, default=_json_default)
    except (TypeError, ValueError) as exc:
        raise SetupConfigError(
            f"Could not JSON-encode setup state for tool '{tool}': {exc}"
        ) from exc
    group.attrs[PYIRENA_CONFIG_ATTR] = payload


def read_setup_config(file_path: Path | str, group_path: str,
                      expected_tool: str) -> Optional[Dict[str, Any]]:
    """Read the embedded setup state from an existing NXcanSAS file.

    Parameters
    ----------
    file_path
        Path to a HDF5/NXcanSAS file.
    group_path
        Internal path to the tool's results group
        (e.g. ``"entry/unified_fit_results"``).
    expected_tool
        The tool the caller expects.  Raises ``SetupConfigToolMismatch`` if
        the attribute exists but identifies a different tool — guards against
        loading a Sizes setup into the Unified Fit panel.

    Returns
    -------
    The ``state`` dict (after stripping the envelope header), or ``None`` if
    the file or group exists but no ``_pyirena_config`` attribute is present.
    """
    fp = Path(file_path)
    if not fp.exists():
        raise FileNotFoundError(fp)

    with h5py.File(fp, "r") as f:
        if group_path not in f:
            return None
        group = f[group_path]
        raw = group.attrs.get(PYIRENA_CONFIG_ATTR)
        if raw is None:
            return None
        # h5py returns bytes for older files; normalise to str
        if isinstance(raw, bytes):
            raw = raw.decode("utf-8")

    try:
        envelope = json.loads(raw)
    except (TypeError, ValueError) as exc:
        raise SetupConfigError(
            f"'{PYIRENA_CONFIG_ATTR}' on {fp}:{group_path} is not valid JSON: {exc}"
        ) from exc

    header = envelope.get(PYIRENA_CONFIG_ATTR, {})
    actual_tool = header.get("tool")
    if actual_tool and actual_tool != expected_tool:
        raise SetupConfigToolMismatch(expected_tool, actual_tool, str(fp))

    state = envelope.get("state")
    if not isinstance(state, dict):
        raise SetupConfigError(
            f"'{PYIRENA_CONFIG_ATTR}' on {fp}:{group_path} has no 'state' object."
        )
    return state


def _json_default(obj: Any) -> Any:
    """Best-effort coercion for non-JSON-native values (numpy scalars, paths, …)."""
    import numpy as np

    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, Path):
        return str(obj)
    raise TypeError(f"Object of type {type(obj).__name__} is not JSON serialisable")
