"""Screenshots for docs/tools/data-explorer.md.

Dataset: a scratch folder named ``unified`` holding copies of the three files
in validationData/unified/, each with a Unified Fit stored in it (fitted here
through pyirena.api.control from rough start values, as in the Unified Fit
page). Only the folder's base name appears in the window.
"""

from __future__ import annotations

import shutil
from pathlib import Path

from ._common import VALIDATION, Context, draw_callouts, shot, widget_rect

SIZE = (1150, 700)  # the window's own default size
GRAPH_SIZE = (760, 560)

# file -> (number of levels, start values, level options)
_FITS = {
    "unified_one_level.h5": (1, {"G_1": 80, "Rg_1": 150, "B_1": 1e-7}, {}),
    "unified_one_level_slitsmeared.h5": (1, {"G_1": 80, "Rg_1": 150, "B_1": 1e-7}, {}),
    "unified_two_level.h5": (
        2,
        {"G_1": 10, "Rg_1": 100, "B_1": 1e-7, "G_2": 3000, "Rg_2": 1000, "B_2": 1e-6, "P_2": 3},
        {2: "link_RGCO"},
    ),
}


def demo_folder(ctx: Context) -> Path:
    """Fresh ``<scratch>/explorer/unified`` with the three files, each fitted and saved."""
    from pyirena.api.control import (
        close_session,
        free_parameter,
        open_dataset,
        run_fit,
        save_fit,
        select_model,
        set_level_option,
        set_parameter_value,
    )

    folder = ctx.scratch / "explorer" / "unified"
    if folder.exists():
        shutil.rmtree(folder)
    folder.mkdir(parents=True)
    for name, (nlevels, start, options) in _FITS.items():
        dst = folder / name
        shutil.copy2(VALIDATION / "unified" / name, dst)
        sid = open_dataset(str(dst))["session_id"]
        select_model(sid, nlevels=nlevels)
        for level, option in options.items():
            set_level_option(sid, level, option, True)
        for pname, value in start.items():
            set_parameter_value(sid, pname, value)
        if nlevels == 2:
            free_parameter(sid, "P_2")
        run_fit(sid)
        save_fit(sid)          # into the scratch copy
        close_session(sid)
    return folder


def _open(ctx: Context):
    from pyirena.gui.hdf5viewer.main_window import HDF5ViewerWindow

    folder = demo_folder(ctx)
    win = HDF5ViewerWindow(initial_folder=str(folder))
    ctx.show(win, SIZE)
    tree = win._file_tree._tree
    tree.selectAll()
    ctx.settle()
    return win


def _expand_results(ctx: Context, win, filename: str = "unified_two_level.h5") -> None:
    """Show *filename* in the centre panel with entry/unified_fit_results expanded."""
    win._hdf5_browser.load_file(str(Path(win._file_tree.root_folder) / filename))
    ctx.settle()
    tree = win._hdf5_browser._tree

    def find(parent_items, name):
        for item in parent_items:
            if item.text(0) == name:
                return item
        return None

    entry = find([tree.topLevelItem(i) for i in range(tree.topLevelItemCount())], "entry")
    if entry is not None:
        entry.setExpanded(True)       # children load lazily on expand
        ctx.settle(5)
        res = find([entry.child(i) for i in range(entry.childCount())], "unified_fit_results")
        if res is not None:
            res.setExpanded(True)
            ctx.settle(5)
            tree.setCurrentItem(res)
    ctx.settle()


@shot("data-explorer-01-window")
def window(ctx: Context):
    win = _open(ctx)
    _expand_results(ctx, win)
    pc = win._plot_controls
    pc._cb_nxcansas.setChecked(True)
    pc._cb_unified.setChecked(True)
    ctx.settle()
    r = lambda w: widget_rect(w, win)  # noqa: E731
    ft = win._file_tree
    marks = [
        (1, r(ft._folder_btn), "tl", (0, 2)),
        (2, r(ft._sort_combo), "r", (-14, 0)),
        (3, r(ft._filter_edit), "r", (-14, 0)),
        (4, r(ft._tree), "tr", (-16, 16)),
        (5, r(win._hdf5_browser._tree), "tr", (-16, 16)),
        (6, r(win._hdf5_browser._value_edit), "r", (-14, 0)),
        (7, r(pc._tabs.tabBar()), "r", (14, 0)),
        (8, r(pc._cb_nxcansas), "l", (-6, 0)),
        (9, r(pc._new_graph_btn), "l", (-2, 0)),
        (10, r(win._status_bar), "l", (130, 0)),
    ]
    return draw_callouts(win.grab(), marks)


@shot("data-explorer-02-graph")
def graph(ctx: Context):
    win = _open(ctx)
    pc = win._plot_controls
    pc._cb_nxcansas.setChecked(True)
    pc._cb_unified.setChecked(True)
    pc._rb_first.setChecked(False)
    _all_files_radio(pc).setChecked(True)
    ctx.settle()
    pc._on_new_graph()
    ctx.settle()
    gw = win._graph_windows[-1]
    ctx.show(gw, GRAPH_SIZE)
    return gw.grab()


def _all_files_radio(pc):
    from pyirena.gui._qt import QRadioButton

    return next(b for b in pc.findChildren(QRadioButton) if b.text().startswith("All selected"))


@shot("data-explorer-03-collect")
def collect(ctx: Context):
    win = _open(ctx)
    pc = win._plot_controls
    pc._tabs.setCurrentIndex(1)
    pc._collect_type.setCurrentIndex(0)          # Unified Fit
    ctx.settle()
    idx = pc._collect_item.findText("Rg")
    pc._collect_item.setCurrentIndex(max(idx, 0))
    pc._collect_index.setValue(1)
    ctx.settle()
    pc._on_collect()
    ctx.settle()
    cw = win._collect_windows[-1]
    ctx.show(cw, GRAPH_SIZE)
    return cw.grab()
