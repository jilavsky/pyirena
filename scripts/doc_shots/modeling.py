"""Screenshots for docs/tools/modeling.md.

Dataset: validationData/modeling/modeling_sizedist_plus_unified.h5 — population 1
log-normal spheres (median 30 Å, sdeviation 0.3, scale 0.005, contrast 100);
population 2 a Unified level (G = 20000, Rg = 4000 Å, P = 3.5, B = 1.56e-8);
background 0.01 cm⁻¹ (validationData/README.md). Starting values below are
deliberately rough. The fit is the Standard (local) one, which is deterministic.
"""

from __future__ import annotations

from ._common import Context, draw_callouts, scene_item_rect, shot, widget_rect

SIZE = (1500, 1000)  # wider than the 1300×900 minimum, so the controls fit unclipped
DATASET = "modeling/modeling_sizedist_plus_unified.h5"


def _open(ctx: Context, dataset: str = DATASET):
    from pyirena.gui.modeling_panel import ModelingPanel

    panel = ModelingPanel()
    ctx.show(panel, SIZE)
    panel.data_loader.load_path(ctx.dataset(dataset))
    ctx.settle()
    return panel


def _two_populations(ctx: Context, panel) -> None:
    """P1 log-normal spheres, P2 Unified level, rough start values."""
    from pyirena.core.modeling import SizeDistPopulation, UnifiedLevelPopulation

    spheres = SizeDistPopulation(
        dist_type="lognormal",
        dist_params={"min_size": 0.0, "mean_size": 40.0, "sdeviation": 0.3},
        dist_params_fit={"min_size": False, "mean_size": True, "sdeviation": True},
        scale=0.002, contrast=100.0, form_factor="sphere")
    level = UnifiedLevelPopulation(G=1e4, Rg=3000.0, B=1e-8, P=3.5)
    panel._pop_widgets[0].from_population(spheres)
    panel._pop_widgets[1].from_population(level)
    panel.bg_edit.setText("0.01")
    panel.bg_fit_cb.setChecked(True)
    ctx.settle()


def _fitted(ctx: Context):
    panel = _open(ctx)
    _two_populations(ctx, panel)
    panel.run_fit()
    ctx.settle(5)
    ctx.wait_until(lambda: panel._fit_worker is None or not panel._fit_worker.isRunning(), 300)
    return panel


@shot("modeling-01-window")
def window(ctx: Context):
    panel = _fitted(ctx)
    g = panel.graph
    r = lambda w: widget_rect(w, panel)  # noqa: E731
    return draw_callouts(panel.grab(), [
        (1, r(panel.data_loader._edit), "r", (-14, 0)),
        (2, r(panel.slit_length_edit), "r", (16, 0)),
        (3, r(panel.q_range_fields), "c", (0, -22)),
        (4, r(panel.no_limits_cb), "r", (14, 0)),
        (5, r(panel.pop_tabs.tabBar()), "r", (16, 0)),
        (6, r(panel.bg_edit), "t", (0, -12)),
        (7, r(panel.fit_method_combo), "t", (0, -12)),
        (8, r(panel.btn_fit), "r", (-16, 0)),
        (9, r(panel.chi2_lbl), "t", (0, -14)),
        (10, r(panel.btn_revert), "l", (-2, 0)),
        (11, r(panel.btn_mc), "tr", (-6, 2)),
        (12, scene_item_rect(g.gl, g.iq_plot, panel), "tr", (-16, 16)),
        (13, scene_item_rect(g.gl, g.resid_plot, panel), "tr", (-16, 16)),
        (14, scene_item_rect(g.gl, g.dist_plot, panel), "tr", (-16, 16)),
        (15, r(g.status_lbl), "r", (-16, 0)),
    ])


@shot("modeling-02-population-tab")
def population_tab(ctx: Context):
    panel = _open(ctx)
    _two_populations(ctx, panel)
    panel.pop_tabs.setCurrentIndex(0)
    panel.graph_model()
    ctx.settle()
    return panel.grab()


@shot("modeling-03-fit-result")
def fit_result(ctx: Context):
    panel = _fitted(ctx)
    panel.show_individual_cb.setChecked(True)
    panel.pop_tabs.setCurrentIndex(1)
    ctx.settle()
    panel._results_to_graph()
    ctx.settle()
    return panel.graph.grab()
