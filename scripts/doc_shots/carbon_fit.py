"""Screenshots for docs/tools/carbon-fit.md.

There is no carbon dataset in validationData/ or testData/, so the shots
synthesise one, exactly as pyirena/tests/test_carbon_fit.py
(test_fit_recovers_parameters_from_noisy_data) does: the model's defaults
with S_macro = 3.0e4 cm²/cm³, φ = 0.20 and pore radius 7.0 Å, evaluated on
900 log-spaced points 0.001–4.47 Å⁻¹, 2 % Gaussian noise with seed 7. The
curve is written as a three-column text file into the scratch folder and
loaded like any user file. The panel then starts from its own defaults.
"""

from __future__ import annotations

from pathlib import Path

from ._common import Context, draw_callouts, scene_item_rect, shot, widget_rect

SIZE = (1300, 1000)
TRUTH = {"S_macro": 3.0e4, "phi": 0.20, "pore_radius": 7.0}


def synthetic_carbon(folder: Path) -> Path:
    """Write the synthetic carbon curve (Q, I, dI) to *folder*; return the path."""
    import numpy as np

    from pyirena.core.carbon_fit import CarbonFitModel

    truth = CarbonFitModel()
    truth.background.S_macro = TRUTH["S_macro"]
    truth.saxs.phi = TRUTH["phi"]
    truth.saxs.pore_radius = TRUTH["pore_radius"]
    q = np.logspace(-3, 0.65, 900)
    rng = np.random.default_rng(7)
    intensity = truth.evaluate(q) * (1.0 + rng.normal(0.0, 0.02, q.shape))
    path = folder / "synthetic_carbon.dat"
    np.savetxt(path, np.column_stack([q, intensity, 0.02 * np.abs(intensity)]),
               header="Q[1/A] I[1/cm] dI[1/cm]  synthetic carbon, seed 7")
    return path


def _open(ctx: Context):
    import tempfile

    from pyirena.gui.carbon_fit_panel import CarbonFitPanel

    panel = CarbonFitPanel()
    ctx.show(panel, SIZE)
    path = synthetic_carbon(Path(tempfile.mkdtemp(dir=ctx.scratch)))
    panel.data_loader.load_path(str(path))
    ctx.settle()
    return panel


def _fitted(ctx: Context):
    panel = _open(ctx)
    panel._on_fit()
    ctx.wait_jobs(panel)
    return panel


@shot("carbon-fit-01-window")
def window(ctx: Context):
    panel = _fitted(ctx)
    g = panel.graph_window
    r = lambda w: widget_rect(w, panel)  # noqa: E731
    marks = [
        (1, r(panel.data_loader), "r", (-80, 0)),
        (2, r(panel.tabs.tabBar()), "r", (14, 0)),
        (3, r(panel.q_range_fields), "tr", (-14, 0)),
        (4, r(panel.show_components_check), "r", (10, 0)),
        (5, r(panel.fit_btn), "r", (-14, 0)),
        (6, r(panel.mc_btn), "tr", (-4, 0)),
        (7, r(panel.save_btn), "r", (-14, 0)),
        (8, scene_item_rect(g.graphics_layout, g.main_plot, panel), "tr", (-16, 16)),
        (9, scene_item_rect(g.graphics_layout, g.residuals_plot, panel), "tr", (-16, 16)),
        (10, r(panel.status_label), "l", (2, 0)),
    ]
    return draw_callouts(panel.grab(), marks)


@shot("carbon-fit-02-waxs-zoom")
def waxs_zoom(ctx: Context):
    panel = _fitted(ctx)
    panel.tabs.setCurrentIndex(2)        # WAXS region
    panel.zoom_check.setChecked(True)
    ctx.settle()
    return panel.grab()


@shot("carbon-fit-03-results")
def results(ctx: Context):
    panel = _fitted(ctx)
    panel.tabs.setCurrentIndex(4)        # Results
    ctx.settle()
    return panel.grab()
