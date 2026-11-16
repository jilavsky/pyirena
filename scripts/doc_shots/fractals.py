"""Screenshots for docs/tools/fractals.md.

Reference file: a scratch copy of validationData/modeling/modeling_mass_fractal.h5
(Teixeira mass fractal of R = 40 Å spheres, Dv = 2.4, ξ = 800 Å) with a
two-level Unified fit stored in it first, so the panel shows the data and the
Unified curve. The levels are deliberately NOT RgCutoff-linked, so the panel
detects no fractal pair and shows no targets: a two-level Unified fit of this
synthetic set is not a valid fractal reference (free B: c = 0.67 < 1; with the
mass-fractal option c = 1.00 and level 1's B collapses to its bound; 16-11-2026),
and the targets picture (fractals-05-targets) is left for a real dataset.
Aggregates are grown with a fixed random seed so the shots are reproducible.

The PyVista/VTK 3D view crashes the offscreen Qt platform (segfault), so the
script switches it off before the window is built: the 3D pane then shows
its "not available" placeholder, and the 3D picture is a manual shot.
"""

from __future__ import annotations

from pathlib import Path

from ._common import Context, draw_callouts, scene_item_rect, shot, widget_rect

SIZE = (1400, 1450)
REFERENCE = "modeling/modeling_mass_fractal.h5"
SEED = 12345


def _reference_with_unified(ctx: Context) -> str:
    """Copy the mass-fractal dataset and store an unlinked two-level Unified fit in it."""
    from pyirena.api.control import (
        free_parameter,
        open_dataset,
        run_fit,
        save_fit,
        select_model,
        set_parameter_value,
    )

    path = ctx.dataset(REFERENCE)
    sid = open_dataset(path)["session_id"]
    select_model(sid, nlevels=2)   # no link_RGCO: see the module docstring
    for name, value in {"G_1": 2, "Rg_1": 30, "B_1": 1e-5, "P_1": 4,
                        "G_2": 500, "Rg_2": 1500, "B_2": 1e-5, "P_2": 2.4}.items():
        set_parameter_value(sid, name, value)
    free_parameter(sid, "P_2")
    run_fit(sid)
    save_fit(sid, output_path=path)
    return path


def _open(ctx: Context, with_reference: bool = True):
    from pyirena.gui import saxs_morph_3d
    from pyirena.gui.fractals_panel import FractalsGraphWindow

    saxs_morph_3d.HAS_PYVISTA = False  # VTK cannot render offscreen (see module doc)

    win = FractalsGraphWindow()
    ctx.show(win, SIZE)
    if with_reference:
        path = _reference_with_unified(ctx)
        win.panel._load_unified_from_nexus(Path(path))
        # Show the file name, not the scratch path (display only).
        win.panel.nexus_path_edit.setText(Path(path).name)
        ctx.settle()
    return win


def _mode_tabs(panel):
    """The Grow One / Grow Many / Optimizer QTabWidget (no attribute on the panel)."""
    from pyirena.gui._qt import QTabWidget

    w = panel.btn_grow
    while not isinstance(w, QTabWidget):
        w = w.parentWidget()
    return w


def _grow_one(ctx: Context, win, z: int = 250) -> None:
    p = win.panel
    p.z_spin.setValue(z)
    p.seed_spin.setValue(SEED)
    n = p.agg_list.count()
    p._on_grow_one()
    ctx.wait_until(lambda: p.agg_list.count() > n, timeout=300)
    # Let the "done" job entry disappear (removed 1.5 s after completion).
    ctx.wait_until(lambda: p.jobs_list.count() == 0, timeout=10)


@shot("fractals-01-window")
def window(ctx: Context):
    win = _open(ctx)
    _grow_one(ctx, win)
    p = win.panel
    r = lambda w: widget_rect(w, win)  # noqa: E731
    marks = [
        (1, r(p.nexus_path_edit), "r", (-14, 0)),
        (2, r(p.z_spin.parentWidget()), "tr", (-14, 14)),
        (3, r(p.btn_grow.parentWidget().parentWidget()), "tr", (-14, 14)),
        (4, r(p.jobs_list), "tr", (-14, 14)),
        (5, r(p.agg_list), "tr", (-14, 14)),
        (6, r(p.targets_summary_label), "r", (-14, 0)),
        (7, r(p.btn_mc), "r", (-14, 0)),
        (8, r(p.btn_save), "r", (-14, 0)),
        (9, scene_item_rect(p.iq_plot.scene().views()[0], p.iq_plot, win), "tr", (-16, 16)),
        (10, r(p.slice_viewer), "tr", (-16, 16)),
        (11, r(p.voxel3d_viewer), "tr", (-16, 16)),
    ]
    return draw_callouts(win.grab(), marks)


@shot("fractals-02-montecarlo-iq")
def montecarlo_iq(ctx: Context):
    win = _open(ctx)
    _grow_one(ctx, win)
    p = win.panel
    mc = p._mc_worker
    configure = mc.configure
    # The panel passes no seed (random primary-sphere polydispersity); pin it
    # so the shot is reproducible.
    mc.configure = lambda agg, q, **kw: configure(agg, q, **{**kw, "seed": SEED})
    p._on_compute_mc()
    ctx.wait_until(lambda: not mc.isRunning() and p._active.i_montecarlo is not None,
                   timeout=300)
    return p.iq_plot.scene().views()[0].grab()


@shot("fractals-03-optimizer")
def optimizer(ctx: Context):
    win = _open(ctx)
    p = win.panel
    p.z_spin.setValue(100)          # small aggregates keep the run short
    p.seed_spin.setValue(SEED)
    _mode_tabs(p).setCurrentIndex(2)  # Optimizer
    p.opt_dmin_spin.setValue(1.3)
    p.opt_c_spin.setValue(1.4)
    p.opt_iter_spin.setValue(2)
    ctx.settle()
    p._on_optimize()
    ctx.wait_until(lambda: any("BEST" in (a.label or "") for a in p._aggregates), timeout=600)
    ctx.wait_until(lambda: p.jobs_list.count() == 0, timeout=10)
    # The trial labels are long: widen the controls pane and show only it.
    from pyirena.gui._qt import QSplitter

    splitter = p.findChildren(QSplitter)[0]
    splitter.setSizes([700, 700])
    ctx.settle()
    return splitter.widget(0).grab()
