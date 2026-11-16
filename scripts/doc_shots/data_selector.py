"""Screenshots for docs/tools/data-selector.md.

Dataset: the same folder of three stored Unified Fits as the Data Explorer
shots (copies of validationData/unified/, fitted here). The Data Selector
shows the full folder path, so the folder is made at a fixed, neutral place
(``/tmp/pyirena-docs/unified``; the system temp folder on Windows) instead of
the script's randomly named scratch folder, which keeps the image the same
from run to run. It is rebuilt on every run.
"""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

from ._common import Context, draw_callouts, shot, widget_rect
from .data_explorer import demo_folder

SIZE = (1000, 720)  # the right column fits without scrolling
RESULT_SIZE = (800, 600)


def _base() -> Path:
    root = Path("/tmp") if os.name == "posix" and Path("/tmp").is_dir() else Path(
        tempfile.gettempdir())
    return root / "pyirena-docs"


def _open(ctx: Context, select_all: bool = True):
    from pyirena.state.state_manager import StateManager

    # No GitHub update check (network) and no "new version" banner in the shot.
    sm = StateManager()
    sm.set("data_selector", "check_for_updates", False)
    sm.save()

    folder = demo_folder(ctx, base=_base())
    from pyirena.gui.data_selector.panel import DataSelectorPanel

    panel = DataSelectorPanel()
    ctx.show(panel, SIZE)
    panel._load_folder(str(folder))
    if select_all:
        panel.file_list.selectAll()
    panel.unified_fit_result_checkbox.setChecked(True)
    ctx.settle()
    return panel


@shot("data-selector-01-window")
def window(ctx: Context):
    panel = _open(ctx)
    r = lambda w: widget_rect(w, panel)  # noqa: E731
    marks = [
        (1, r(panel.folder_button), "tl", (0, 0)),
        (2, r(panel.file_type_combo), "r", (-14, 0)),
        (3, r(panel.sort_combo), "r", (-14, 0)),
        (4, r(panel.filter_input), "r", (-14, 0)),
        (5, r(panel.file_list), "tr", (-16, 16)),
        (6, r(panel.locate_logs_button), "r", (16, 0)),
        (7, r(panel.data_checkbox), "l", (-14, 0)),
        (8, r(panel.plot_button), "l", (-12, 0)),
        (9, r(panel.hdf5_viewer_button), "l", (-12, 0)),
        (10, r(panel.unified_fit_button), "l", (-12, 0)),
        (11, r(panel.data_merge_button), "l", (-12, 0)),
        (12, r(panel.status_label), "l", (230, 0)),
    ]
    return draw_callouts(panel.grab(), marks)


@shot("data-selector-02-graph")
def graph(ctx: Context):
    panel = _open(ctx)
    panel.plot_selected_files()
    ctx.settle()
    win = panel.unified_fit_results_window or panel.graph_window
    ctx.show(win, RESULT_SIZE)
    return win.grab()


@shot("data-selector-03-tabulate")
def tabulate(ctx: Context):
    panel = _open(ctx)
    panel.tabulate_results()
    ctx.settle()
    win = panel.tabulate_results_window
    ctx.show(win, (900, 420))
    return win.grab()
