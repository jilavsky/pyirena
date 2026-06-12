"""Shared "Load Setup from File…" dialog flow used by every fit panel.

Each panel exposes a small wrapper that supplies:
  - the tool name and NXcanSAS group path,
  - a default folder for the file picker,
  - a callback that applies a restored state dict to the GUI,
  - a status-reporting hook (one for success, one for warnings).

This module centralises the file-dialog, the error-handling, and the
tool-mismatch guard so the five panels don't repeat 30 lines of identical
boilerplate each.  Keeping it in one place also makes the UX uniform.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Callable, Dict, Optional, Tuple

try:
    from PySide6.QtWidgets import QFileDialog, QMessageBox, QWidget
except ImportError:
    try:
        from PyQt6.QtWidgets import QFileDialog, QMessageBox, QWidget
    except ImportError:
        from PyQt5.QtWidgets import QFileDialog, QMessageBox, QWidget  # type: ignore

from pyirena.io.setup_config import (
    SetupConfigError, SetupConfigToolMismatch, read_setup_config,
)


# Map tool name → results-group path inside the NXcanSAS file.
TOOL_GROUP_PATH = {
    "unified_fit":  "entry/unified_fit_results",
    "sizes":        "entry/sizes_results",
    "modeling":     "entry/modeling_results",
    "simple_fits":  "entry/simple_fit_results",
    "waxs_peakfit": "entry/waxs_peakfit_results",
}

# Human-readable label per tool — used in dialog titles and status messages.
TOOL_LABEL = {
    "unified_fit":  "Unified Fit",
    "sizes":        "Sizes",
    "modeling":     "Modeling",
    "simple_fits":  "Simple Fits",
    "waxs_peakfit": "WAXS Peak Fit",
}


def prompt_and_load_setup(
    parent:        QWidget,
    tool:          str,
    default_folder: str,
    apply_state:   Callable[[Dict[str, Any]], None],
    on_status:     Optional[Callable[[str], None]] = None,
    suggested_path: Optional[str] = None,
) -> Optional[Path]:
    """Run the "Load Setup from File…" interaction.

    Parameters
    ----------
    parent
        Qt parent for dialogs.
    tool
        One of the keys in ``TOOL_GROUP_PATH`` (e.g. ``"unified_fit"``).
    default_folder
        Folder the file dialog opens in.  Use the currently loaded data
        file's folder so the AI-produced .h5 is one click away.
    apply_state
        Callback invoked with the restored state dict.  Typically wraps
        the panel's existing ``apply_state`` / ``_apply_state`` / ``load_state``
        method.  Any exception is caught and surfaced to the user.
    on_status
        Optional sink for short success/info messages (panel status bar).
    suggested_path
        Pre-select this file in the dialog.  Default = no pre-selection.

    Returns
    -------
    Path
        The path the user picked, even when no setup was found / on mismatch.
    None
        User cancelled the dialog.
    """
    if tool not in TOOL_GROUP_PATH:
        raise ValueError(f"Unknown tool '{tool}' — must be one of {list(TOOL_GROUP_PATH)}")

    label = TOOL_LABEL.get(tool, tool)
    title = f"Load {label} Setup from File"
    start = suggested_path or default_folder
    path_str, _ = QFileDialog.getOpenFileName(
        parent, title, start,
        "NXcanSAS HDF5 (*.h5 *.hdf5 *.nxs);;All Files (*)",
    )
    if not path_str:
        return None

    path = Path(path_str)
    try:
        state = read_setup_config(path, TOOL_GROUP_PATH[tool], expected_tool=tool)
    except FileNotFoundError:
        QMessageBox.warning(parent, "File not found", f"Could not find:\n{path}")
        return path
    except SetupConfigToolMismatch as exc:
        other_label = TOOL_LABEL.get(exc.actual, exc.actual)
        QMessageBox.warning(
            parent, "Wrong tool",
            f"This file contains a {other_label} setup, not {label}.\n\n"
            f"Open the {other_label} panel to restore it from this file.",
        )
        return path
    except SetupConfigError as exc:
        QMessageBox.warning(parent, "Setup unreadable", str(exc))
        return path
    except Exception as exc:
        QMessageBox.critical(parent, "Load failed",
                             f"Could not read setup from:\n{path}\n\n{exc}")
        return path

    if state is None:
        QMessageBox.information(
            parent, "No setup stored",
            f"No {label} setup is stored in:\n{path.name}\n\n"
            "The file was likely produced by an older pyirena version, by a "
            "different tool, or not by pyirena at all.",
        )
        return path

    try:
        apply_state(state)
    except Exception as exc:
        QMessageBox.critical(
            parent, "Could not apply setup",
            f"The setup was read but could not be applied:\n{exc}",
        )
        return path

    if on_status is not None:
        try:
            on_status(f"{label} setup restored from {path.name}")
        except Exception:
            pass

    return path
