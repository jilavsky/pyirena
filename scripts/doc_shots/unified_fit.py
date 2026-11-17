"""Screenshots for docs/tools/unified-fit.md.

Dataset: validationData/unified/unified_two_level.h5 — a two-level Unified
curve with known parameters (level 1: G = 8, Rg = 120 Å, P = 4; level 2:
G = 4000, Rg = 1200 Å, P = 3.2, RgCO linked to Rg of level 1; background
0.02 cm⁻¹). Starting values below are deliberately rough, as a user's would be.
"""

from __future__ import annotations

from ._common import Context, draw_callouts, scene_item_rect, shot, widget_rect
from ._package import wait_jobs

SIZE = (1200, 960)  # the panel's own default size
DATASET = "unified/unified_two_level.h5"

# Rough starting values, roughly what the cursor fits give.
START = {
    1: {"G": 10, "Rg": 100, "B": 1e-7, "P": 4},
    2: {"G": 3000, "Rg": 1000, "B": 1e-6, "P": 3},
}


def _open(ctx: Context, dataset: str = DATASET):
    from pyirena.gui.unified_fit import UnifiedFitPanel

    panel = UnifiedFitPanel()
    ctx.show(panel, SIZE)
    panel.data_loader.load_path(ctx.dataset(dataset))
    ctx.settle()
    return panel


def _two_level_start(ctx: Context, panel) -> None:
    """Two levels, rough start values, every parameter that should be fitted ticked."""
    panel.num_levels_spin.setValue(2)
    ctx.settle()
    for level, values in START.items():
        lw = panel.level_widgets[level - 1]
        lw.set_parameters({**values, "fit_G": True, "fit_Rg": True, "fit_B": True,
                           "fit_P": level == 2})
        lw.fix_limits()
    panel.level_widgets[1].link_rgco_check.setChecked(True)
    panel.fit_background_check.setChecked(True)
    ctx.settle()


def _fitted(ctx: Context):
    panel = _open(ctx)
    _two_level_start(ctx, panel)
    panel.run_fit()
    wait_jobs(ctx, panel)
    return panel


@shot("unified-fit-01-window")
def window(ctx: Context):
    panel = _fitted(ctx)
    gw = panel.graph_window
    r = lambda w: widget_rect(w, panel)  # noqa: E731
    out = (16, 0)  # just right of the widget's right edge
    marks = [
        (1, r(panel.data_loader._edit), "r", (-14, 0)),
        (2, r(panel.slit_length_edit), "r", out),
        (3, r(panel.q_range_fields), "r", (-4, -22)),
        (4, r(panel.num_levels_spin), "tr", (0, -2)),
        (5, r(panel.level_tabs.tabBar()), "r", out),
        (6, r(panel.fit_background_check), "r", (14, 0)),
        (7, r(panel.fit_button), "r", (-16, 0)),
        (8, r(panel.results_graphs_button), "l", (-2, 0)),
        (9, r(panel.analyze_results_button), "tr", (-6, 2)),
        (10, r(gw.tab_widget.tabBar()), "r", out),
        (11, scene_item_rect(gw.graphics_layout, gw.main_plot, panel), "tr", (-16, 16)),
        (12, scene_item_rect(gw.graphics_layout, gw.residual_plot, panel), "tr", (-16, 16)),
        (13, r(gw.status_message), "l", (2, 0)),
    ]
    return draw_callouts(panel.grab(), marks)


@shot("unified-fit-02-porod-tab")
def porod_tab(ctx: Context):
    panel = _open(ctx)
    panel.graph_window.tab_widget.setCurrentIndex(1)
    ctx.settle()
    return panel.graph_window.grab()


@shot("unified-fit-03-guinier-fit")
def guinier_fit(ctx: Context):
    panel = _open(ctx)
    _two_level_start(ctx, panel)
    # Level 2 (Rg ≈ 1200 Å): its Guinier knee is the clean one, near
    # Q ≈ 0.001 Å⁻¹; level 1 sits on level 2's power-law tail.
    panel.graph_window.set_cursor_range(0.0002, 0.0012)
    panel.level_tabs.setCurrentIndex(1)
    ctx.settle()
    panel.fit_local_guinier(2)
    wait_jobs(ctx, panel)
    return panel.grab()


@shot("unified-fit-04-fit-result")
def fit_result(ctx: Context):
    panel = _fitted(ctx)
    panel.display_results_on_graph()
    ctx.settle()
    return panel.graph_window.grab()


@shot("unified-fit-05-slit-smeared")
def slit_smeared(ctx: Context):
    # One Unified level, infinite-slit smeared (dQl = 0.03 Å⁻¹ in the file).
    panel = _open(ctx, "unified/unified_one_level_slitsmeared.h5")
    lw = panel.level_widgets[0]
    lw.set_parameters({"G": 80, "Rg": 150, "B": 1e-7, "P": 4,
                       "fit_G": True, "fit_Rg": True, "fit_B": True})
    lw.fix_limits()
    panel.fit_background_check.setChecked(True)
    ctx.settle()
    panel.run_fit()
    wait_jobs(ctx, panel)
    return panel.grab()
