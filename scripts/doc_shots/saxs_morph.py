"""Screenshots for docs/tools/saxs-morph.md.

Dataset: validationData/simple_fits/simple_debye_bueche.h5 — a Debye-Bueche
curve (random two-phase structure, correlation length 80 Å, no background),
the kind of morphology the Gaussian-random-field method describes. With
φ = 0.3 the invariant Q* = Prefactor·η²·π/4 = 1.96e-3 cm⁻¹ Å⁻³ gives
Δρ² = Q*/(2π² φ(1−φ) 1e-4) = 4.74 × 10²⁰ cm⁻⁴, which checks the
"derive Δρ²" mode. The file has no background, so both pre-fits are skipped
(B = 0 and background = 0 are the panel defaults).

The 3D PyVista view crashes under QT_QPA_PLATFORM=offscreen, so PyVista is
switched off here (the panel then shows its "not installed" hint in that
pane). The automated window shot covers that pane with a neutral box; the 3D
view itself is a manual shot (saxs-morph-04-3d-view).
"""

from __future__ import annotations

from ._common import Context, draw_callouts, scene_item_rect, shot, widget_rect

SIZE = (1500, 1150)  # minimum is 1300 x 900; taller so the result block fits
DATASET = "simple_fits/simple_debye_bueche.h5"


def _open(ctx: Context):
    import pyirena.gui.saxs_morph_3d as morph3d

    morph3d.HAS_PYVISTA = False  # VTK segfaults offscreen
    from pyirena.gui.saxs_morph_panel import SaxsMorphPanel

    panel = SaxsMorphPanel()
    ctx.show(panel, SIZE)
    panel.data_loader.load_path(ctx.dataset(DATASET))
    ctx.settle()
    return panel


def _bg_tabs(panel):
    """The Power-law Bckg / Flat Bckg tab widget (the panel keeps no reference)."""
    from pyirena.gui._qt import QTabWidget

    w = panel.bg_widget
    while not isinstance(w, QTabWidget):
        w = w.parentWidget()
    return w


def _calculated(ctx: Context):
    panel = _open(ctx)
    panel.box_size_edit.setText("4000")
    panel.seed_edit.setText("42")
    panel.phi_edit.setText("0.3")  # Input φ → derive Δρ² (the default mode)
    # Modelling window: Q_min ≥ 3 × 2π/L (L = 4000 Å), or the panel asks first.
    panel.graph.set_q_range(0.005, 0.2)
    ctx.settle()
    panel._on_calculate()
    ctx.wait_until(lambda: not panel._calculating, timeout=300)
    return panel


def _cover_3d(image, panel):
    """Grey box over the 3D pane, which cannot render offscreen."""
    from PySide6.QtGui import QPainter

    from pyirena.gui._qt import QColor, Qt

    rect = widget_rect(panel.graph.voxel3d_viewer, panel)
    p = QPainter(image)
    p.fillRect(rect, QColor("#eef0f2"))
    p.setPen(QColor("#555555"))
    p.drawText(rect, int(Qt.AlignmentFlag.AlignCenter),
               "3D isosurface view\n(see the 3D viewer screenshot below)")
    p.end()
    return image


@shot("saxs-morph-01-window")
def window(ctx: Context):
    panel = _calculated(ctx)
    g = panel.graph
    r = lambda w: widget_rect(w, panel)  # noqa: E731
    marks = [
        (1, r(panel.data_loader._edit), "r", (-14, 0)),
        (2, r(panel.qmax_lbl), "r", (40, 0)),
        (3, r(panel.render_size_combo.parentWidget()), "tl", (100, 6)),
        (4, r(panel.mode_combo.parentWidget()), "tl", (172, 0)),
        (5, r(_bg_tabs(panel).tabBar()), "r", (16, 0)),
        (6, r(panel.btn_calc), "l", (16, 0)),
        (7, r(panel.result_lbl), "tr", (-12, 12)),
        (8, r(panel.btn_save), "r", (-14, 0)),
        (9, scene_item_rect(g.gl, g.iq_plot, panel), "tr", (-16, 16)),
        (10, r(g.slice_viewer), "tr", (-16, 16)),
        (11, r(g.voxel3d_viewer), "tr", (-16, 16)),
        (12, r(g.status_lbl), "r", (-16, 0)),
    ]
    return draw_callouts(_cover_3d(panel.grab(), panel), marks)


@shot("saxs-morph-02-calculated")
def calculated(ctx: Context):
    panel = _calculated(ctx)
    return panel.graph.gl.grab()


@shot("saxs-morph-03-slice")
def slice_view(ctx: Context):
    panel = _calculated(ctx)
    return panel.graph.slice_viewer.parentWidget().grab()
