"""Screenshots for docs/tools/simple-fits.md.

Datasets (validationData/README.md):

* simple_fits/simple_guinier.h5 — Guinier, I0 = 250 cm⁻¹, Rg = 45 Å.
* simple_fits/simple_sphere.h5 — monodisperse sphere, Scale = 8 cm⁻¹,
  R = 120 Å, flat background 0.01 cm⁻¹.
* sizes/sizes_sphere_flat_background.h5 — for the Invariant: volume fraction
  0.006, contrast 100, flat background 0.02 cm⁻¹ above Q ≈ 0.15 Å⁻¹.
"""

from __future__ import annotations

from ._common import Context, draw_callouts, scene_item_rect, shot, widget_rect
from ._package import wait_jobs

SIZE = (1200, 900)  # the coded 1200×700 clips the parameter grid; taller shows it all


def _open(ctx: Context, dataset: str, model: str):
    from pyirena.gui.simple_fits_panel import SimpleFitsPanel

    panel = SimpleFitsPanel()
    ctx.show(panel, SIZE)
    panel.data_loader.load_path(ctx.dataset(dataset))
    ctx.settle()
    panel.model_combo.setCurrentText(model)
    ctx.settle()
    return panel


def _set(panel, name: str, value: float, fit: bool | None = None) -> None:
    panel._param_value_edits[name].setText(f"{value:g}")
    if fit is not None:
        panel._param_fit_checks[name].setChecked(fit)


def _fit(ctx: Context, panel) -> None:
    panel._run_fit()
    wait_jobs(ctx, panel)


def _guinier_fitted(ctx: Context):
    panel = _open(ctx, "simple_fits/simple_guinier.h5", "Guinier")
    _set(panel, "I0", 100.0)
    _set(panel, "Rg", 30.0)
    ctx.settle()
    _fit(ctx, panel)
    return panel


@shot("simple-fits-01-window")
def window(ctx: Context):
    panel = _guinier_fitted(ctx)
    gw = panel.graph_window
    gl = gw.graphics_layout
    r = lambda w: widget_rect(w, panel)  # noqa: E731
    return draw_callouts(panel.grab(), [
        (1, r(panel.data_loader._edit), "r", (-14, 0)),
        (2, r(panel.slit_length_edit), "r", (16, 0)),
        (3, r(panel.model_combo), "r", (16, 0)),
        (4, r(panel.q_range_fields), "c", (0, -22)),
        (5, r(panel.complex_bg_check), "r", (14, 0)),
        (6, r(panel.params_box), "tr", (-16, 8)),
        (7, r(panel.fit_btn), "r", (-16, 0)),
        (8, r(panel.chi2_label), "r", (16, 0)),
        (9, r(panel.results_to_graph_btn), "l", (-2, 0)),
        (10, r(panel.unc_btn), "tr", (-6, 2)),
        (11, scene_item_rect(gl, gw.main_plot, panel), "tr", (-16, 16)),
        (12, scene_item_rect(gl, gw.residuals_plot, panel), "tr", (-16, 16)),
        (13, scene_item_rect(gl, gw.lin_plot, panel), "tr", (-16, 16)),
        (14, r(panel.status_label), "l", (720, 0)),
    ])


@shot("simple-fits-02-sphere")
def sphere(ctx: Context):
    panel = _open(ctx, "simple_fits/simple_sphere.h5", "Sphere")
    panel.complex_bg_check.setChecked(True)
    ctx.settle()
    _set(panel, "Scale", 5.0)
    _set(panel, "R", 100.0)
    _set(panel, "BG_B", 0.0, fit=False)
    _set(panel, "BG_P", 4.0, fit=False)
    _set(panel, "BG_flat", 0.01, fit=True)
    ctx.settle()
    _fit(ctx, panel)
    return panel.graph_window.grab()


@shot("simple-fits-03-invariant")
def invariant(ctx: Context):
    panel = _open(ctx, "sizes/sizes_sphere_flat_background.h5", "Invariant")
    panel.complex_bg_check.setChecked(True)
    ctx.settle()
    _set(panel, "Contrast", 100.0)
    _set(panel, "BG_B", 0.0, fit=False)
    ctx.settle()
    # Flat background from the high-Q plateau, then integrate the whole curve.
    panel.graph_window.set_cursor_range(0.2, 0.35)
    ctx.settle()
    panel._prefit_bg_flat()
    ctx.settle()
    panel.graph_window.set_cursor_range(0.0008, 0.35)
    ctx.settle()
    _fit(ctx, panel)
    return panel.grab()
