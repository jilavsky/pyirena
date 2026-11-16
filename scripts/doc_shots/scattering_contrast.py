"""Screenshots for docs/tools/scattering-contrast.md.

No dataset: the calculator works on compound definitions. The example is
magnetite (Fe3O4, 5.17 g/cm³) particles in water (H2O, 1.0 g/cm³); the energy
scan brackets the Fe K edge (7.11 keV), where the anomalous contrast changes.
"""

from __future__ import annotations

from ._common import Context, draw_callouts, shot, widget_rect

SIZE = (1120, 1140)  # default is 1120×740; taller so no control is scrolled away
SCAN_SIZE = (900, 760)


def _open(ctx: Context):
    # The compound library lives at a fixed path in the user's home folder
    # (pyirena/io/contrast_io.py), not under PYIRENA_STATE_DIR: point it at the
    # scratch folder so the shot never shows the user's own saved compounds.
    import pyirena.io.contrast_io as cio

    cio.DEFAULT_LIBRARY_PATH = ctx.scratch / "contrast_compounds.h5"
    from pyirena.gui.contrast_panel import ContrastPanel

    panel = ContrastPanel()
    ctx.show(panel, SIZE)
    for n, (name, formula, density) in {
            1: ("magnetite", "Fe3O4", 5.17), 2: ("water", "H2O", 1.0)}.items():
        refs = panel._refs(n)
        refs["name"].setText(name)
        refs["formula"].setText(formula)
        refs["mode"].setCurrentIndex(0)
        refs["density"].setValue(density)
        panel._update_isotope_table(n)
    ctx.settle()
    return panel


def _button(panel, text: str):
    from pyirena.gui._qt import QPushButton

    return next(b for b in panel.findChildren(QPushButton) if b.text() == text)


@shot("scattering-contrast-01-window")
def window(ctx: Context):
    panel = _open(ctx)
    panel._on_calculate()
    ctx.settle()
    r = lambda w: widget_rect(w, panel)  # noqa: E731
    c1, c2 = panel._refs(1), panel._refs(2)
    marks = [
        (1, r(c1["formula"]), "r", (-14, 0)),
        (2, r(c1["iso_table"]), "r", (-14, 0)),
        (3, r(c1["lib_combo"]), "r", (14, 0)),
        (4, r(c2["formula"]), "r", (-14, 0)),
        (5, r(panel._energy_spin), "r", (-14, 0)),
        (6, r(panel._e_start), "r", (-14, 0)),
        (7, r(_button(panel, "Calculate")), "l", (16, 0)),
        (8, r(panel._tbl), "bl", (16, -16)),
        (9, r(panel._status_lbl), "l", (-2, 0)),
        (10, r(_button(panel, "Export Results CSV")), "l", (-2, 0)),
        (11, r(_button(panel, "Show Graphs")), "r", (-14, 0)),
    ]
    return draw_callouts(panel.grab(), marks)


@shot("scattering-contrast-02-energy-scan")
def energy_scan(ctx: Context):
    panel = _open(ctx)
    panel._e_start.setValue(6.5)
    panel._e_end.setValue(8.0)
    panel._n_pts.setValue(300)
    panel._on_energy_scan()
    ctx.settle()
    win = panel._graph_win
    ctx.show(win, SCAN_SIZE)
    return win.grab()
