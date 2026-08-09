"""
Source-level guard for the "every plot exports the same way" contract.

Before consolidation there were eight parallel implementations of "add export
to a plot", so each new export feature had to be written eight times and some
panels silently missed out.  These checks read source text only — no Qt, no
display — so they run in every CI job and fail the build if a ninth one
appears.

If this test fails on a plot you added, the fix is one call::

    from pyirena.gui.plot_export import attach_plot_export

    attach_plot_export(my_plot, self, 'my_tool_iq', window=self.graphics_layout)

or, if the plot is built with ``make_sas_plot``, simply pass ``parent_widget``.
"""

from __future__ import annotations

import re
from pathlib import Path

GUI_DIR = Path(__file__).resolve().parents[1] / "gui"

# The one module allowed to drive pyqtgraph's exporters directly.
_EXPORTER_OWNER = "plot_export.py"

# Modules with a deliberate, specialised exporter of their own.  The HDF5
# viewer's graph window writes NXcanSAS metadata into its PNG/HDF5/ITX output
# that the generic curve exporters do not carry; it still uses the shared
# clipboard copy and folder memory.
_ALLOWED_OWN_SAVE = {"export.py"}


def _gui_sources() -> list[Path]:
    return sorted(p for p in GUI_DIR.rglob("*.py"))


def test_only_plot_export_touches_pyqtgraph_exporters():
    offenders = []
    for path in _gui_sources():
        if path.name == _EXPORTER_OWNER:
            continue
        text = path.read_text(encoding="utf-8")
        if re.search(r"\b(ImageExporter|SVGExporter)\b", text):
            offenders.append(str(path.relative_to(GUI_DIR.parent.parent)))
    assert not offenders, (
        "These modules use a pyqtgraph exporter directly instead of "
        "pyirena.gui.plot_export:\n  " + "\n  ".join(offenders)
    )


def test_no_panel_rolls_its_own_jpeg_menu_entry():
    """The menu wording is shared, so it cannot drift between panels again."""
    offenders = []
    pattern = re.compile(r"addAction\(\s*[\"'][^\"']*(JPEG|PNG)[^\"']*[\"']")
    for path in _gui_sources():
        if path.name in _ALLOWED_OWN_SAVE | {_EXPORTER_OWNER}:
            continue
        if pattern.search(path.read_text(encoding="utf-8")):
            offenders.append(str(path.relative_to(GUI_DIR.parent.parent)))
    assert not offenders, (
        "These modules add their own image-export menu entry; call "
        "attach_plot_export() instead:\n  " + "\n  ".join(offenders)
    )


def test_export_dialogs_do_not_default_to_the_home_directory():
    """``Path.home()`` in a save dialog is the bug A8 was about.

    Export dialogs must start from the shared remembered folder
    (``plot_export.export_folder``) so the saved graph lands next to the data.
    """
    offenders = []
    pattern = re.compile(
        r"getSaveFileName\((?:[^)]|\n)*?Path\.home\(\)", re.MULTILINE
    )
    for path in _gui_sources():
        if pattern.search(path.read_text(encoding="utf-8")):
            offenders.append(str(path.relative_to(GUI_DIR.parent.parent)))
    assert not offenders, (
        "Save dialogs defaulting to the home directory instead of the "
        "remembered export folder:\n  " + "\n  ".join(offenders)
    )


def test_plot_export_exposes_the_public_api():
    src = (GUI_DIR / "plot_export.py").read_text(encoding="utf-8")
    for name in (
        "attach_plot_export",
        "collect_plot_curves",
        "copy_plot_to_clipboard",
        "curves_to_csv_text",
        "export_folder",
        "save_curves_as_csv",
        "save_itx_from_plot",
        "save_plot_image",
        "save_widget_image",
    ):
        assert f'"{name}"' in src, f"{name} missing from plot_export.__all__"
