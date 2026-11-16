"""Screenshots for docs/tools/data-merge.md.

Datasets: validationData/merge/merge_usaxs.h5 (DS1, low Q, correct scale) and
merge_saxs.h5 (DS2, high Q, the same curve deliberately scaled by 1.15). With
"Scale: DS2" the merge should recover 1/1.15 = 0.8696 and a background near 0.
Each file is copied into its own scratch folder, which the panel then uses as
the DS1 / DS2 folder, so "Create Default Folder" and "Save Merged Data" write
next to the copies, never into validationData/.
"""

from __future__ import annotations

import os
import shutil
import tempfile
from pathlib import Path

from ._common import VALIDATION, Context, draw_callouts, shot, widget_rect

SIZE = (1500, 820)
DS1 = "merge/merge_usaxs.h5"
DS2 = "merge/merge_saxs.h5"
OVERLAP = (0.01, 0.04)  # Å⁻¹; both branches have data here


def _load(panel, number: int, path: str) -> None:
    selector = panel._ds1 if number == 1 else panel._ds2
    selector.set_folder(os.path.dirname(path))
    items = selector.file_list.findItems(os.path.basename(path),
                                         _match_exactly())
    selector.file_list.setCurrentItem(items[0])
    (panel._on_ds1_double_clicked if number == 1 else panel._on_ds2_double_clicked)(items[0])


def _match_exactly():
    from pyirena.gui._qt import Qt

    return Qt.MatchFlag.MatchExactly


def _copy(ctx: Context, relpath: str, folder: str) -> str:
    """Copy into <scratch>/<random>/<folder>/ — the panel shows the folder name,
    so it must not be a random temp name (the image would change every run)."""
    dst = Path(tempfile.mkdtemp(dir=ctx.scratch)) / folder
    dst.mkdir()
    return str(shutil.copy2(VALIDATION / relpath, dst))


def _open(ctx: Context):
    from pyirena.gui.data_merge_panel import DataMergePanel

    panel = DataMergePanel()
    ctx.show(panel, SIZE)
    _load(panel, 1, _copy(ctx, DS1, "USAXS"))
    _load(panel, 2, _copy(ctx, DS2, "SAXS"))
    panel._graph.init_cursors(*OVERLAP)
    panel._update_cursor_display()
    ctx.settle()
    return panel


def _merged(ctx: Context):
    panel = _open(ctx)
    panel._scale_ds_combo.setCurrentText("DS2")
    panel._optimize_merge()
    ctx.settle()
    return panel


@shot("data-merge-01-window")
def window(ctx: Context):
    panel = _merged(ctx)
    r = lambda w: widget_rect(w, panel)  # noqa: E731
    marks = [
        (1, r(panel._ds1.file_list), "c", (0, 0)),
        (2, r(panel._ds2.file_list), "c", (0, 0)),
        (3, r(panel._mode_combo), "r", (-10, 0)),
        (4, r(panel._scale_ds_combo), "r", (-10, 0)),
        (5, r(panel._optimize_btn), "l", (-6, 0)),
        (6, r(panel._q_min_edit), "r", (-12, 0)),
        (7, r(panel._graph), "tr", (-60, 50)),
        (8, r(panel._out_folder_label), "r", (-14, 0)),
        (9, r(panel._save_btn), "r", (-14, 0)),
        (10, r(panel._match_chk), "r", (14, 0)),
        (11, r(panel._status), "l", (480, 0)),
    ]
    return draw_callouts(panel.grab(), marks)


@shot("data-merge-02-loaded")
def loaded(ctx: Context):
    panel = _open(ctx)
    return panel._graph.grab()
