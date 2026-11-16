"""Screenshots for docs/tools/waxs-peak-fit.md.

Dataset: validationData/waxs/waxs_three_gauss.h5 — three Gaussian peaks
(Q0 = 1.85, 2.55, 3.05 Å⁻¹; FWHM 0.09, 0.12, 0.15 Å⁻¹) on a linear
background 8 + 1.5·Q cm⁻¹. The Diffraction Lines shot overlays
testData/Si.cif (Dans_Diffraction must be installed); silicon is not in the
synthetic data, the shot only shows what the overlay looks like.
"""

from __future__ import annotations

from ._common import TEST_DATA, Context, draw_callouts, scene_item_rect, shot, widget_rect

SIZE = (1200, 1020)  # the panel's own default size
DATASET = "waxs/waxs_three_gauss.h5"


def _open(ctx: Context, dataset: str = DATASET):
    from pyirena.gui.waxs_peakfit_panel import WAXSPeakFitPanel

    panel = WAXSPeakFitPanel()
    ctx.show(panel, SIZE)
    panel.data_loader.load_path(ctx.dataset(dataset))
    ctx.settle()
    return panel


def _remove_false_detections(ctx: Context, panel, fraction: float = 0.2) -> None:
    """Press "Remove" on detected peaks lower than *fraction* of the tallest one."""
    rows = list(panel._peak_rows)
    heights = [r.get_params()["A"]["value"] for r in rows]
    for row, a in zip(rows, heights, strict=True):
        if a < fraction * max(heights):
            panel._remove_peak_row(row)
    ctx.settle()


def _fitted(ctx: Context):
    panel = _open(ctx)
    panel._bg_combo.setCurrentText("Linear")   # the generating background
    ctx.settle()
    panel._find_peaks()
    ctx.settle()
    _remove_false_detections(ctx, panel)
    panel._run_fit()
    ctx.wait_jobs(panel)
    return panel


@shot("waxs-peak-fit-01-window")
def window(ctx: Context):
    panel = _fitted(ctx)
    g = panel._graph
    r = lambda w: widget_rect(w, panel)  # noqa: E731
    marks = [
        (1, r(panel.data_loader), "r", (-80, 0)),
        (2, r(panel._no_limits_chk), "r", (12, 0)),
        (3, r(panel._bg_combo), "r", (14, 0)),
        (4, r(panel._find_btn), "r", (-14, 0)),
        (5, r(panel._peaks_scroll), "tr", (-26, 14)),
        (6, r(panel._fit_btn), "t", (0, -14)),
        (7, r(panel._weight_combo), "r", (14, 0)),
        (8, r(panel._ps_crosscorr_cb), "r", (14, 0)),
        (9, r(panel._status), "tr", (-14, 10)),
        (10, scene_item_rect(g.gl, g.main_plot, panel), "tr", (-16, 16)),
        (11, scene_item_rect(g.gl, g.resid_plot, panel), "tr", (-16, 16)),
    ]
    return draw_callouts(panel.grab(), marks)


@shot("waxs-peak-fit-02-find-peaks")
def find_peaks(ctx: Context):
    panel = _open(ctx)
    panel._bg_combo.setCurrentText("Linear")
    ctx.settle()
    panel._find_peaks()
    ctx.settle()
    return panel.grab()


@shot("waxs-peak-fit-03-fit-result")
def fit_result(ctx: Context):
    panel = _fitted(ctx)
    panel._results_to_graphs()
    ctx.settle()
    return panel._graph.grab()


@shot("waxs-peak-fit-04-diffraction-lines")
def diffraction_lines(ctx: Context):
    panel = _fitted(ctx)
    dl = panel._diffraction_panel
    tabs = dl.parentWidget()
    while tabs is not None and not hasattr(tabs, "setCurrentWidget"):
        tabs = tabs.parentWidget()
    tabs.setCurrentWidget(dl)
    dl._add_entry(ctx.dataset("Si.cif", root=TEST_DATA), show_hkl=True)
    ctx.settle()
    return panel.grab()
