"""Screenshots for docs/tools/size-distribution.md.

Datasets (validationData/README.md):

* sizes/sizes_sphere_lognormal.h5 — spheres, log-normal volume distribution,
  R_median = 100 Å, sigma_log = 0.25, volume fraction 0.01, contrast 100,
  no background.
* sizes/sizes_sphere_flat_background.h5 — R_median = 120 Å, sigma_log = 0.22,
  volume fraction 0.006, contrast 100, flat background 0.02 cm⁻¹ dominating
  above Q ≈ 0.15 Å⁻¹.

The Monte-Carlo uncertainty button is not shot: the GUI runs it unseeded, so
the image would change on every run.
"""

from __future__ import annotations

from ._common import Context, draw_callouts, scene_item_rect, shot, widget_rect

SIZE = (1200, 960)  # the panel's own default size
LOGNORMAL = "sizes/sizes_sphere_lognormal.h5"
FLAT_BG = "sizes/sizes_sphere_flat_background.h5"


def _open(ctx: Context, dataset: str):
    from pyirena.gui.sizes_panel import SizesFitPanel

    panel = SizesFitPanel()
    ctx.show(panel, SIZE)
    panel.data_loader.load_path(ctx.dataset(dataset))
    ctx.settle()
    return panel


def _setup(ctx: Context, panel, method: str = "regularization") -> None:
    """Contrast of the datasets (100), grid 10–1000 Å, the chosen method."""
    panel.contrast_edit.setText("100")
    panel.contrast_edit.editingFinished.emit()
    panel.rmin_edit.setText("10")
    panel.rmax_edit.setText("1000")
    panel.method_combo.setCurrentText(method)
    ctx.settle()


def _tabs(panel):
    from pyirena.gui._qt import QTabWidget

    # The Sizes / Background tab widget of the control panel.
    for tw in panel.findChildren(QTabWidget):
        if tw.count() == 2 and tw.tabText(1) == "Background":
            return tw
    raise LookupError("Sizes/Background tabs not found")


def _flat_background_fitted(ctx: Context):
    """Flat-background dataset: background fitted at Q 0.2–0.35, sizes at 0.0008–0.15."""
    panel = _open(ctx, FLAT_BG)
    _setup(ctx, panel)
    panel.bg_q_min_edit.setText("0.2")
    panel.bg_q_max_edit.setText("0.35")
    panel.fit_background_action()
    ctx.settle()
    return panel


@shot("size-distribution-01-window")
def window(ctx: Context):
    panel = _open(ctx, LOGNORMAL)
    _setup(ctx, panel)
    panel.run_fit()
    ctx.wait_jobs(panel)
    gw = panel.graph_window
    r = lambda w: widget_rect(w, panel)  # noqa: E731
    gl = gw.graphics_layout
    return draw_callouts(panel.grab(), [
        (1, r(panel.data_loader._edit), "r", (-14, 0)),
        (2, r(panel.slit_length_edit), "r", (16, 0)),
        (3, r(_tabs(panel).tabBar()), "l", (165, 0)),
        (4, r(panel.contrast_edit), "r", (100, 0)),
        (5, r(panel.qmax_display), "r", (40, 0)),
        (6, r(panel.rmax_edit), "r", (100, 0)),
        (7, r(panel.error_scale_edit), "r", (100, 0)),
        (8, r(panel.method_combo), "r", (16, 0)),
        (9, r(panel.fit_button), "t", (0, -2)),
        (10, r(panel.result_vf), "r", (-14, 0)),
        (11, r(panel.results_to_graph_button), "l", (-2, 0)),
        (12, r(panel.uncertainty_button), "tr", (-6, 2)),
        (13, scene_item_rect(gl, gw.main_plot, panel), "tr", (-16, 16)),
        (14, scene_item_rect(gl, gw.residuals_plot, panel), "tr", (-16, 16)),
        (15, scene_item_rect(gl, gw.distribution_plot, panel), "tr", (-16, 30)),
        (16, r(gw.status_message), "l", (2, 0)),
    ])


@shot("size-distribution-02-background-tab")
def background_tab(ctx: Context):
    panel = _flat_background_fitted(ctx)
    _tabs(panel).setCurrentIndex(1)
    ctx.settle()
    return panel.grab()


@shot("size-distribution-03-fit-result")
def fit_result(ctx: Context):
    panel = _flat_background_fitted(ctx)
    panel.graph_window.set_cursor_range(0.0008, 0.15)
    ctx.settle()
    panel.run_fit()
    ctx.wait_jobs(panel)
    panel._results_to_graph()
    ctx.settle()
    return panel.graph_window.grab()
