"""Screenshots for docs/tools/data-manipulation.md.

Datasets: validationData/manipulation/manipulate_input.h5 is
manipulate_reference.h5 multiplied by 2.5 with 0.3 cm⁻¹ added. So
Scale I = 0.4 with Background = 0.12 (I_out = scale·I − background) gives
back the reference exactly, and subtracting the reference as buffer with
Buffer scale = 2.5 leaves a flat 0.3 cm⁻¹. Both files are copied into one
scratch folder named "manipulation" (the panel shows the folder name, so it
must not be a random temp name); outputs never reach validationData/.
"""

from __future__ import annotations

import shutil
import tempfile
from pathlib import Path

from ._common import Context, draw_callouts, shot, widget_rect
from ._package import VALIDATION

SIZE = (1350, 760)
INPUT = "manipulate_input.h5"
REFERENCE = "manipulate_reference.h5"
TAB_SCALE, TAB_SUBTRACT = 0, 4


def _open(ctx: Context, tab: int):
    from pyirena.gui.data_manipulation_panel import DataManipulationPanel

    folder = Path(tempfile.mkdtemp(dir=ctx.scratch)) / "manipulation"
    folder.mkdir()
    for name in (INPUT, REFERENCE):
        shutil.copy2(VALIDATION / "manipulation" / name, folder)

    panel = DataManipulationPanel()
    ctx.show(panel, SIZE)
    panel.set_folder(str(folder))
    panel._tabs.setCurrentIndex(tab)
    ctx.settle()
    return panel


def _select(panel, name: str):
    from pyirena.gui._qt import Qt

    item = panel._fb.file_list.findItems(name, Qt.MatchFlag.MatchExactly)[0]
    panel._fb.file_list.setCurrentItem(item)
    return item


def _scaled(ctx: Context):
    panel = _open(ctx, TAB_SCALE)
    item = _select(panel, INPUT)
    panel._on_file_double_clicked(item)
    panel._scale_I_edit.setText("0.4")
    panel._scale_bg_edit.setText("0.12")
    panel._preview_scale()
    ctx.settle()
    return panel


@shot("data-manipulation-01-window")
def window(ctx: Context):
    panel = _scaled(ctx)
    r = lambda w: widget_rect(w, panel)  # noqa: E731
    fb = panel._fb
    marks = [
        (1, r(fb.folder_btn), "r", (-14, 0)),
        (2, r(fb.file_list), "c", (0, 40)),
        (3, r(panel._tabs.tabBar()), "r", (16, 0)),
        (4, r(panel._scale_unc_edit), "r", (16, 0)),
        (5, r(panel._graph), "tr", (-60, 40)),
        (6, r(panel._out_folder_label), "r", (-14, 0)),
        (7, r(panel._save_btn), "r", (-14, 0)),
        (8, r(panel._status), "l", (420, 0)),
    ]
    return draw_callouts(panel.grab(), marks)


@shot("data-manipulation-02-scale-result")
def scale_result(ctx: Context):
    panel = _scaled(ctx)
    return panel._graph.grab()


@shot("data-manipulation-03-subtract")
def subtract(ctx: Context):
    panel = _open(ctx, TAB_SUBTRACT)
    panel._set_buffer(REFERENCE)
    item = _select(panel, INPUT)
    panel._on_file_double_clicked(item)
    panel._sub_scale_edit.setText("2.5")
    panel._preview_subtract()
    ctx.settle()
    return panel.grab()
