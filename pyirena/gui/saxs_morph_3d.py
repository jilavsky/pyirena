"""
3D and 2D voxelgram viewers for the SAXS Morph tool.

PyVista is imported lazily — if VTK is missing, ``HAS_PYVISTA`` is False and
``Voxel3DViewer.__init__`` shows a single QLabel with an install hint instead
of crashing. The Slice2DViewer needs only pyqtgraph and works regardless.

The popout helper reparents a viewer widget into a top-level QDialog the user
can resize / maximise; closing returns it to its slot in the panel.

Public API
----------
HAS_PYVISTA : bool                    True when pyvistaqt could be imported.
PYVISTA_INSTALL_HINT : str            Helpful install hint string.
Voxel3DViewer(QWidget)                3D isosurface viewer.
Slice2DViewer(QWidget)                2D slice viewer with axis combo + slider.
make_popout_button(widget) -> QPushButton
"""

from __future__ import annotations

from typing import Optional

import numpy as np

try:
    from PySide6.QtWidgets import (
        QWidget, QMainWindow, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
        QSlider, QSplitter, QDialog, QComboBox, QColorDialog, QFileDialog,
        QMenu, QSizePolicy, QApplication,
    )
    from PySide6.QtCore import Qt, Signal
    from PySide6.QtGui import QAction, QColor
except ImportError:
    try:
        from PyQt6.QtWidgets import (
            QWidget, QMainWindow, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
            QSlider, QSplitter, QDialog, QComboBox, QColorDialog, QFileDialog,
            QMenu, QSizePolicy, QApplication,
        )
        from PyQt6.QtCore import Qt, pyqtSignal as Signal
        from PyQt6.QtGui import QAction, QColor
    except ImportError:
        from PyQt5.QtWidgets import (
            QWidget, QMainWindow, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
            QSlider, QSplitter, QDialog, QComboBox, QColorDialog, QFileDialog,
            QMenu, QSizePolicy, QApplication,
        )
        from PyQt5.QtCore import Qt, pyqtSignal as Signal
        from PyQt5.QtGui import QAction, QColor

import pyqtgraph as pg


# ---------------------------------------------------------------------------
# Lazy PyVista import
# ---------------------------------------------------------------------------

PYVISTA_INSTALL_HINT = (
    'PyVista is not installed.\n\n'
    'Install with:\n'
    '    pip install pyirena[gui3d]\n'
    'or\n'
    '    pip install pyvista pyvistaqt vtk\n\n'
    '2D slice and I(Q) views still work without it.'
)

try:
    import pyvista as pv
    from pyvistaqt import QtInteractor
    HAS_PYVISTA = True
    _IMPORT_ERROR: Optional[Exception] = None
except Exception as _e:
    HAS_PYVISTA = False
    _IMPORT_ERROR = _e
    pv = None  # type: ignore
    QtInteractor = QWidget  # type: ignore


# ---------------------------------------------------------------------------
# Silence VTK shader / driver warnings.
# ---------------------------------------------------------------------------
# Some shader-based post-processing passes (SSAO, cast shadows, occasionally
# EDL) emit verbose vtkOpenGLPolyDataMapper / vtkShaderProgram errors when
# their fragment shader fails to compile on the current GPU.  This is a
# well-known issue on macOS where Apple has deprecated OpenGL — the modern
# VTK shader templates miss substitutions on Apple's translation layer and
# the driver rejects the program.  The viewer recovers (it falls back to the
# previous valid shader) but the console fills with red noise per render.
#
# We silence VTK's global warning display once on import.  The cost is
# losing visibility into other VTK warnings that don't bubble up through
# Python exceptions — acceptable for an end-user tool.  Real exceptions
# coming through the Python layer are still surfaced normally.
if HAS_PYVISTA:
    try:
        import vtk
        vtk.vtkObject.GlobalWarningDisplayOff()
    except Exception:
        pass


# ---------------------------------------------------------------------------
# 3D viewer
# ---------------------------------------------------------------------------

class Voxel3DViewer(QWidget):
    """Renders a binary uint8 voxel cube as an isosurface (flying_edges).

    For binary 0/1 data, isosurface is dramatically faster than full volume
    rendering on integrated GPUs (1-3 s mesh build at 256**3, then a static
    triangle mesh that rotates at 60 fps).  Right-click context menu gives
    the user control over view, color, outline, and screenshot.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumSize(300, 300)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

        self._actor = None
        self._bounds_actor = None
        self._show_outline = True
        # Default isosurface color: emerald green (#2ecc71 ≈ RGB 46/204/113).
        # Matches the convention used by molecular-visualisation tools
        # (PyMOL / VMD default highlight) and reads well against the
        # white background under any of our lighting modes.
        self._iso_color = (0.18, 0.80, 0.44)
        self._voxelgram = None
        self._pitch_A = 1.0
        self._mesh = None

        # ── Visual-tuning state ──────────────────────────────────────────
        # All toggles default to OFF / "default" so the initial render
        # matches the current behaviour; users opt in via the right-click
        # menu.  State persists for the lifetime of the viewer instance.
        self._lighting_mode = 'default'      # 'default' (headlight) | 'lightkit'
        self._edl_on = False                 # eye-dome lighting
        self._ssao_on = False                # screen-space ambient occlusion
        self._shadows_on = False             # cast shadows from positional lights
        self._parallel_projection = False    # False=perspective, True=orthographic

        if HAS_PYVISTA:
            self._build_ui()
        else:
            self._build_disabled()

    # ----- enabled UI ------------------------------------------------------

    def _build_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        self.plotter = QtInteractor(self)
        self.plotter.set_background('white')
        self.plotter.show_axes()
        layout.addWidget(self.plotter.interactor)

        self.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.customContextMenuRequested.connect(self._show_context_menu)

    def _show_context_menu(self, pos):
        menu = QMenu(self)

        # ── Reset view ───────────────────────────────────────────────────
        a_reset = QAction('Reset view', self)
        a_reset.triggered.connect(self.reset_view)
        menu.addAction(a_reset)

        menu.addSeparator()

        # ── View submenu ─────────────────────────────────────────────────
        view_menu = menu.addMenu('View')
        for label, callback in (
            ('XY plane (top)',   self._view_xy),
            ('XZ plane (front)', self._view_xz),
            ('YZ plane (side)',  self._view_yz),
            ('Isometric',        self._view_isometric),
        ):
            a = QAction(label, self)
            a.triggered.connect(callback)
            view_menu.addAction(a)
        view_menu.addSeparator()
        a_persp = QAction('Perspective projection', self)
        a_persp.setCheckable(True)
        a_persp.setChecked(not self._parallel_projection)
        a_persp.triggered.connect(lambda: self._set_projection(parallel=False))
        view_menu.addAction(a_persp)
        a_ortho = QAction('Orthographic projection', self)
        a_ortho.setCheckable(True)
        a_ortho.setChecked(self._parallel_projection)
        a_ortho.triggered.connect(lambda: self._set_projection(parallel=True))
        view_menu.addAction(a_ortho)

        # ── Lighting submenu ─────────────────────────────────────────────
        light_menu = menu.addMenu('Lighting')
        a_default = QAction('Default headlight', self)
        a_default.setCheckable(True)
        a_default.setChecked(self._lighting_mode == 'default')
        a_default.triggered.connect(lambda: self._set_lighting('default'))
        light_menu.addAction(a_default)
        a_kit = QAction('3-point light kit', self)
        a_kit.setCheckable(True)
        a_kit.setChecked(self._lighting_mode == 'lightkit')
        a_kit.triggered.connect(lambda: self._set_lighting('lightkit'))
        light_menu.addAction(a_kit)
        light_menu.addSeparator()
        a_edl = QAction('Eye-dome lighting (silhouette edges)', self)
        a_edl.setCheckable(True)
        a_edl.setChecked(self._edl_on)
        a_edl.triggered.connect(self._toggle_edl)
        light_menu.addAction(a_edl)
        a_ssao = QAction('SSAO (depth shading in crevices)', self)
        a_ssao.setCheckable(True)
        a_ssao.setChecked(self._ssao_on)
        a_ssao.triggered.connect(self._toggle_ssao)
        light_menu.addAction(a_ssao)
        a_shadows = QAction('Cast shadows', self)
        a_shadows.setCheckable(True)
        a_shadows.setChecked(self._shadows_on)
        a_shadows.triggered.connect(self._toggle_shadows)
        a_shadows.setToolTip(
            'Shadows are most visible with the 3-point light kit.\n'
            'With the default headlight (camera-position light) shadows\n'
            'are usually invisible because the light follows the camera.'
        )
        light_menu.addAction(a_shadows)

        menu.addSeparator()

        # ── Existing actions ─────────────────────────────────────────────
        a_color = QAction('Pick isosurface color…', self)
        a_color.triggered.connect(self._pick_color)
        menu.addAction(a_color)

        a_outline = QAction(
            'Hide bounding box' if self._show_outline else 'Show bounding box',
            self,
        )
        a_outline.triggered.connect(self._toggle_outline)
        menu.addAction(a_outline)

        menu.addSeparator()
        a_shot = QAction('Save screenshot…', self)
        a_shot.triggered.connect(self._save_screenshot)
        menu.addAction(a_shot)

        menu.exec(self.mapToGlobal(pos))

    # ----- View handlers ---------------------------------------------------

    def _view_xy(self):
        if HAS_PYVISTA and self.plotter is not None:
            try: self.plotter.view_xy()
            except Exception: pass

    def _view_xz(self):
        if HAS_PYVISTA and self.plotter is not None:
            try: self.plotter.view_xz()
            except Exception: pass

    def _view_yz(self):
        if HAS_PYVISTA and self.plotter is not None:
            try: self.plotter.view_yz()
            except Exception: pass

    def _view_isometric(self):
        if HAS_PYVISTA and self.plotter is not None:
            try: self.plotter.view_isometric()
            except Exception: pass

    def _set_projection(self, parallel: bool):
        if not (HAS_PYVISTA and self.plotter is not None):
            return
        self._parallel_projection = bool(parallel)
        try:
            if parallel:
                self.plotter.enable_parallel_projection()
            else:
                self.plotter.disable_parallel_projection()
            self.plotter.render()
        except Exception as exc:
            print(f'[Voxel3DViewer] projection toggle failed: {exc}')

    # ----- Lighting handlers -----------------------------------------------

    def _set_lighting(self, mode: str):
        """Switch between default headlight and 3-point lightkit.

        Removes all existing lights first to keep the lighting setup
        deterministic — calling ``enable_lightkit`` repeatedly without a
        clean slate would stack additional lights on top.
        """
        if not (HAS_PYVISTA and self.plotter is not None):
            return
        self._lighting_mode = mode
        try:
            self.plotter.remove_all_lights()
        except Exception:
            pass
        try:
            if mode == 'lightkit':
                # PyVista's enable_lightkit sets up a key/fill/back/headlight
                # rig that gives surfaces a much stronger sense of form than
                # the default single camera-mounted headlight.
                self.plotter.enable_lightkit()
            else:
                # Single camera-mounted light (PyVista default)
                light = pv.Light(light_type='headlight')
                self.plotter.add_light(light)
            self.plotter.render()
        except Exception as exc:
            print(f'[Voxel3DViewer] lighting "{mode}" failed: {exc}')

    def _toggle_edl(self):
        """Toggle Eye-dome lighting (post-process edge enhancement).

        EDL darkens silhouette edges so concave / convex features pop
        visually.  Massive readability boost for fractal aggregates and
        any porous structure.  Composes on top of any lighting mode.
        """
        if not (HAS_PYVISTA and self.plotter is not None):
            return
        self._edl_on = not self._edl_on
        try:
            if self._edl_on:
                self.plotter.enable_eye_dome_lighting()
            else:
                self.plotter.disable_eye_dome_lighting()
            self.plotter.render()
        except Exception as exc:
            print(f'[Voxel3DViewer] EDL toggle failed: {exc}')
            self._edl_on = not self._edl_on   # revert state on failure

    def _toggle_ssao(self):
        """Toggle Screen-Space Ambient Occlusion.

        SSAO darkens crevices and concavities by sampling neighbour
        depths.  Excellent for sponge-like saxsMorph volumes — gives a
        true 3-D-cavity feel that flat shading can't.
        """
        if not (HAS_PYVISTA and self.plotter is not None):
            return
        self._ssao_on = not self._ssao_on
        try:
            if self._ssao_on:
                self.plotter.enable_ssao()
            else:
                self.plotter.disable_ssao()
            self.plotter.render()
        except Exception as exc:
            print(f'[Voxel3DViewer] SSAO toggle failed: {exc}')
            self._ssao_on = not self._ssao_on

    def _toggle_shadows(self):
        """Toggle cast shadows from positional lights.

        Most visible with the 3-point light kit; with the default
        headlight (light moves with camera) shadows are usually
        invisible because every visible surface gets directly lit.
        """
        if not (HAS_PYVISTA and self.plotter is not None):
            return
        self._shadows_on = not self._shadows_on
        try:
            if self._shadows_on:
                self.plotter.enable_shadows()
            else:
                self.plotter.disable_shadows()
            self.plotter.render()
        except Exception as exc:
            print(f'[Voxel3DViewer] shadows toggle failed: {exc}')
            self._shadows_on = not self._shadows_on

    # ----- disabled UI -----------------------------------------------------

    def _build_disabled(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        msg = QLabel(PYVISTA_INSTALL_HINT)
        msg.setAlignment(Qt.AlignmentFlag.AlignCenter)
        msg.setWordWrap(True)
        msg.setStyleSheet(
            'background:#fff8e1; color:#856404; font-size:10pt;'
            'padding:20px; border:1px dashed #d4a017; border-radius:6px;'
        )
        layout.addWidget(msg)

    # ----- public API ------------------------------------------------------

    def set_voxelgram(self, voxelgram: np.ndarray, pitch_A: float):
        """Build an isosurface mesh from a binary uint8 voxel cube."""
        if not HAS_PYVISTA:
            return
        self._voxelgram = voxelgram
        self._pitch_A = float(pitch_A)

        self._rebuild_mesh()

    def _rebuild_mesh(self):
        if not HAS_PYVISTA or self._voxelgram is None:
            return

        N = self._voxelgram.shape[0]
        # pv.ImageData with point-data scalars — dimensions = N points per axis.
        grid = pv.ImageData(
            dimensions=(N, N, N),
            spacing=(self._pitch_A, self._pitch_A, self._pitch_A),
            origin=(0.0, 0.0, 0.0),
        )
        grid.point_data['phase'] = (
            self._voxelgram.astype(np.float32).ravel(order='F')
        )

        # Flying-edges isosurface is the fast path for binary data.
        try:
            mesh = grid.contour([0.5], scalars='phase', method='flying_edges')
        except Exception:
            mesh = grid.contour([0.5], scalars='phase')
        self._mesh = mesh

        # Replace any existing actors
        if self._actor is not None:
            try:
                self.plotter.remove_actor(self._actor)
            except Exception:
                pass
        if self._bounds_actor is not None:
            try:
                self.plotter.remove_actor(self._bounds_actor)
            except Exception:
                pass

        self._actor = self.plotter.add_mesh(
            mesh, color=self._iso_color, smooth_shading=True,
            specular=0.3, specular_power=15,
        )
        # Bounding box with Å tick labels (5 per side).
        # show_bounds() draws the cube ruler; if it fails (old VTK), fall back
        # to a plain outline.
        self._bounds_actor = None
        if self._show_outline:
            box_label = f'{N * self._pitch_A:.0f} A'  # noqa: F841 (used in log)
            # NOTE: VTK's font does not include the Å (U+00C5) glyph; using
            # it produces missing-glyph rectangles or "X()".  Use ASCII "[A]".
            try:
                self._bounds_actor = self.plotter.show_bounds(
                    mesh=grid,
                    xtitle='X [A]', ytitle='Y [A]', ztitle='Z [A]',
                    n_xlabels=5, n_ylabels=5, n_zlabels=5,
                    fmt='%.0f', font_size=9, bold=False,
                    ticks='outside', grid=False, all_edges=True,
                    color='black',
                )
            except Exception:
                # Older VTK / pyvistaqt that does not support show_bounds well
                self._bounds_actor = self.plotter.add_mesh(
                    grid.outline(), color='black', line_width=1,
                )
        self.plotter.reset_camera()

    def reset_view(self):
        if HAS_PYVISTA:
            self.plotter.reset_camera()

    def _pick_color(self):
        c0 = QColor(int(self._iso_color[0] * 255),
                    int(self._iso_color[1] * 255),
                    int(self._iso_color[2] * 255))
        c = QColorDialog.getColor(c0, self, 'Pick isosurface color')
        if c.isValid():
            self._iso_color = (c.redF(), c.greenF(), c.blueF())
            self._rebuild_mesh()

    def _toggle_outline(self):
        self._show_outline = not self._show_outline
        self._rebuild_mesh()

    def _save_screenshot(self):
        path, _ = QFileDialog.getSaveFileName(
            self, 'Save 3D screenshot', 'voxelgram_3d.png',
            'PNG (*.png);;JPEG (*.jpg *.jpeg);;All Files (*)',
        )
        if path and HAS_PYVISTA:
            try:
                self.plotter.screenshot(path)
            except Exception as e:
                from PySide6.QtWidgets import QMessageBox
                QMessageBox.warning(self, 'Screenshot failed', str(e))

    def screenshot(self, path: str):
        if HAS_PYVISTA:
            self.plotter.screenshot(path)

    def clear(self):
        if not HAS_PYVISTA:
            return
        if self._actor is not None:
            try:
                self.plotter.remove_actor(self._actor)
            except Exception:
                pass
            self._actor = None
        if self._bounds_actor is not None:
            try:
                self.plotter.remove_actor(self._bounds_actor)
            except Exception:
                pass
            self._bounds_actor = None
        self._voxelgram = None

    def shutdown(self):
        """Tear down the PyVista plotter and the underlying VTK render window
        before Qt destroys the widget hierarchy.

        On macOS, the VTK CocoaRenderWindow holds a reference to a layer-backed
        NSView; if Qt destroys the QNSView while VTK still owns it, a delayed
        signal (e.g. propagateBackingProperties) will invoke vtkCocoaRenderWindow::Render
        on a freed view and segfault the whole process.

        This method:
          1. removes all actors,
          2. calls Finalize() on the render window (releases the GL context),
          3. closes the plotter and nulls our reference (so any late signal
             that does reach Python finds None and bails out).
        Safe to call multiple times.
        """
        if not HAS_PYVISTA:
            return
        plotter = getattr(self, 'plotter', None)
        if plotter is None:
            return
        # Drop our own references first
        self._actor = None
        self._bounds_actor = None
        self._mesh = None
        self._voxelgram = None
        # Tell VTK to release the render window resources
        try:
            ren_win = plotter.render_window
            if ren_win is not None:
                ren_win.Finalize()
        except Exception:
            pass
        # Close the QtInteractor / Plotter
        try:
            plotter.close()
        except Exception:
            pass
        # Drop our reference so late signals see None
        self.plotter = None

    def closeEvent(self, evt):
        self.shutdown()
        super().closeEvent(evt)


# ---------------------------------------------------------------------------
# 2D slice viewer
# ---------------------------------------------------------------------------

class Slice2DViewer(QWidget):
    """One 2D slice of a 3D voxelgram, with axis combo and a position slider."""

    AXIS_LABELS = [
        ('XY plane (Z slice)', 'z'),
        ('XZ plane (Y slice)', 'y'),
        ('YZ plane (X slice)', 'x'),
    ]
    # Physical axis labels for each slice plane: (horizontal_label, vertical_label)
    # ASCII units (no Å glyph) so they render reliably in every Qt font.
    _PLANE_AXIS_LABELS = {
        'z': ('x [A]', 'y [A]'),
        'y': ('x [A]', 'z [A]'),
        'x': ('y [A]', 'z [A]'),
    }

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumSize(300, 300)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

        self._voxelgram: Optional[np.ndarray] = None
        self._pitch_A = 1.0
        self._axis = 'z'
        self._building = True

        self._build_ui()
        self._building = False

    def _build_ui(self):
        lay = QVBoxLayout(self)
        lay.setContentsMargins(2, 2, 2, 2)
        lay.setSpacing(3)

        # Axis selector
        top = QHBoxLayout()
        top.addWidget(QLabel('Slice plane:'))
        self.axis_combo = QComboBox()
        for label, key in self.AXIS_LABELS:
            self.axis_combo.addItem(label, key)
        self.axis_combo.currentIndexChanged.connect(self._on_axis_changed)
        top.addWidget(self.axis_combo)
        top.addStretch()
        lay.addLayout(top)

        # Image view: white = minority phase (void), dark grey = majority
        # phase (solid).  Sharp two-level rendering keeps structural edges
        # crisp without anti-aliasing blur.  Smoothing is NOT applied here —
        # it is reserved for the 3D isosurface viewer.
        plot_item = pg.PlotItem()
        self.image_view = pg.ImageView(view=plot_item)
        self.image_view.ui.histogram.hide()
        self.image_view.ui.roiBtn.hide()
        self.image_view.ui.menuBtn.hide()
        self.image_view.view.setAspectLocked(True)
        self.image_view.view.invertY(False)  # y increases upward (physical)
        # Black axes against white background — match the 3D viewer style.
        for ax_name in ('left', 'bottom', 'top', 'right'):
            ax = plot_item.getAxis(ax_name)
            if ax is not None:
                ax.setPen(pg.mkPen('k', width=1))
                ax.setTextPen(pg.mkPen('k'))
        # Physical axis labels (updated when the slice plane changes).
        # ASCII units so they always render.
        self.image_view.view.setLabel('bottom', 'x [A]', size='9pt', color='k')
        self.image_view.view.setLabel('left', 'y [A]', size='9pt', color='k')
        # White-background canvas (was default dark grey of pyqtgraph)
        plot_item.getViewBox().setBackgroundColor('w')
        # 0 → white (void / minority);  1 → dark grey (solid / majority).
        bw_lut = np.array(
            [[255, 255, 255, 255], [80, 80, 80, 255]], dtype=np.uint8,
        )
        self.image_view.imageItem.setLookupTable(bw_lut)
        self.image_view.imageItem.setLevels([0, 1])
        lay.addWidget(self.image_view)

        # Position slider
        bot = QHBoxLayout()
        self.slider = QSlider(Qt.Orientation.Horizontal)
        self.slider.setRange(0, 0)
        self.slider.valueChanged.connect(self._on_slice_changed)
        bot.addWidget(self.slider)
        self.slice_label = QLabel('—')
        self.slice_label.setMinimumWidth(160)
        self.slice_label.setStyleSheet('font-family: monospace; font-size: 9pt;')
        bot.addWidget(self.slice_label)
        lay.addLayout(bot)

    # ----- public API ------------------------------------------------------

    def set_voxelgram(self, voxelgram: np.ndarray, pitch_A: float):
        self._voxelgram = voxelgram
        self._pitch_A = float(pitch_A)
        N = voxelgram.shape[0]
        self._building = True
        self.slider.setRange(0, N - 1)
        self.slider.setValue(N // 2)
        self._building = False
        self._update_axis_labels()
        self._refresh()

    def clear(self):
        self._voxelgram = None
        self.image_view.clear()
        self.slice_label.setText('—')

    def current_slice_index(self) -> int:
        return int(self.slider.value())

    # ----- internals -------------------------------------------------------

    def _on_axis_changed(self, idx):
        if self._building:
            return
        self._axis = self.axis_combo.itemData(idx)
        self._update_axis_labels()
        self._refresh()

    def _update_axis_labels(self):
        h_label, v_label = self._PLANE_AXIS_LABELS.get(
            self._axis, ('Position [A]', 'Position [A]')
        )
        self.image_view.view.setLabel('bottom', h_label, size='9pt', color='k')
        self.image_view.view.setLabel('left', v_label, size='9pt', color='k')

    def _on_slice_changed(self, _val):
        if self._building:
            return
        self._refresh()

    def _refresh(self):
        if self._voxelgram is None:
            return
        N = self._voxelgram.shape[0]
        idx = int(np.clip(self.slider.value(), 0, N - 1))

        if self._axis == 'z':
            slc = self._voxelgram[:, :, idx]
            coord_label = f'z = {idx * self._pitch_A:.1f} A'
        elif self._axis == 'y':
            slc = self._voxelgram[:, idx, :]
            coord_label = f'y = {idx * self._pitch_A:.1f} A'
        else:
            slc = self._voxelgram[idx, :, :]
            coord_label = f'x = {idx * self._pitch_A:.1f} A'

        # Display with physical Å coordinates on the axes.
        # pos=(0,0): lower-left corner; scale=(pitch, pitch): Å per voxel.
        self.image_view.setImage(
            slc.astype(np.float32),
            autoLevels=False, levels=[0, 1],
            pos=(0.0, 0.0), scale=(self._pitch_A, self._pitch_A),
        )
        self.slice_label.setText(f'Slice {idx + 1} / {N} ({coord_label})')


# ---------------------------------------------------------------------------
# Popout: reparent a widget into a QDialog and back
# ---------------------------------------------------------------------------

class _PopoutDialog(QDialog):
    """Dialog hosting a single widget; on close, re-parents it back."""

    def __init__(self, widget: QWidget, restore_to_layout, restore_index: int,
                 title: str = 'Popout', parent=None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.resize(900, 700)
        self._widget = widget
        self._layout = restore_to_layout
        self._index = restore_index

        lay = QVBoxLayout(self)
        lay.setContentsMargins(2, 2, 2, 2)
        lay.addWidget(widget)

    def closeEvent(self, evt):
        # Return the widget to its original layout slot.  If the original
        # layout has been destroyed (e.g. parent panel closed first), the
        # widget is a 3D viewer that needs explicit shutdown to avoid the
        # macOS VTK CocoaRenderWindow crash.
        restored = False
        try:
            self._layout.insertWidget(self._index, self._widget)
            restored = True
        except Exception:
            try:
                self._layout.addWidget(self._widget)
                restored = True
            except Exception:
                pass
        if not restored:
            # Original layout is gone — shut down the widget explicitly.
            shutdown = getattr(self._widget, 'shutdown', None)
            if callable(shutdown):
                try:
                    shutdown()
                except Exception:
                    pass
        super().closeEvent(evt)


def make_popout_button(widget: QWidget, title: str = 'Popout') -> QPushButton:
    """Return a small button that opens a _PopoutDialog around ``widget``.

    Usage:
        viewer = Voxel3DViewer()
        btn = make_popout_button(viewer, '3D viewer')
        layout.addWidget(viewer)
        layout.addWidget(btn)
    """
    btn = QPushButton('Pop out ⤢')
    btn.setMaximumHeight(22)
    btn.setStyleSheet(
        'background:#ecf0f1;color:#2c3e50;'
        'border:1px solid #bdc3c7;border-radius:3px;'
        'padding:2px 8px;font-size:9pt;'
    )

    state = {'dlg': None}

    def _toggle():
        if state['dlg'] is not None:
            return
        parent_layout = widget.parentWidget().layout() if widget.parentWidget() else None
        if parent_layout is None:
            return
        index = parent_layout.indexOf(widget)
        # Remove from current layout (so QDialog can take ownership)
        parent_layout.removeWidget(widget)
        dlg = _PopoutDialog(widget, parent_layout, index, title=title)
        state['dlg'] = dlg
        btn.setEnabled(False)

        def _on_finished():
            state['dlg'] = None
            btn.setEnabled(True)

        dlg.finished.connect(_on_finished)
        dlg.show()

    btn.clicked.connect(_toggle)
    return btn


# ---------------------------------------------------------------------------
# VoxelViewerWindow — standalone 2D + 3D viewer for stored voxelgrams
# ---------------------------------------------------------------------------

class VoxelViewerWindow(QMainWindow):
    """Standalone window that displays one or more stored voxelgrams.

    Used by the Data Selector "Create Graph" button when the user has
    selected the **Fractals** or **3D saxsMorph** checkbox: instead of
    opening the full SAXS-Morph / Fractals tool (which would re-run the
    engine), this lightweight viewer just loads the saved 3D structure
    from the NeXus file and displays it in the same `Slice2DViewer` +
    `Voxel3DViewer` pair the parent tools use.

    No fitting, no Calculate button, no engine — purely visualisation
    of saved data.

    Parameters
    ----------
    items : list of dict
        Each dict has keys:
            ``voxelgram`` : (N, N, N) uint8 / float32 array
            ``pitch_A``   : float (voxel pitch in Å)
            ``label``     : str (shown in the dropdown / window title)
    title : str
        Window title (e.g. "Fractals — stored aggregates").
    """

    def __init__(self, items: list, title: str = 'Voxel viewer', parent=None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.setMinimumSize(1000, 650)
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, True)

        self._items = list(items)

        central = QWidget()
        v = QVBoxLayout(central)
        v.setContentsMargins(6, 6, 6, 6)
        v.setSpacing(6)

        # Dropdown to switch between items (only if more than one)
        if len(self._items) > 1:
            top = QHBoxLayout()
            top.addWidget(QLabel('Show:'))
            self._combo = QComboBox()
            for item in self._items:
                self._combo.addItem(str(item.get('label', '?')))
            self._combo.currentIndexChanged.connect(self._on_select)
            top.addWidget(self._combo, stretch=1)
            v.addLayout(top)
        else:
            self._combo = None
            if self._items:
                lbl = QLabel(f"<b>{self._items[0].get('label', '')}</b>")
                lbl.setStyleSheet('color:#2c3e50;font-size:10pt;padding:2px;')
                v.addWidget(lbl)

        # 2D + 3D side by side
        body = QWidget()
        h = QHBoxLayout(body)
        h.setContentsMargins(0, 0, 0, 0)
        h.setSpacing(6)

        slice_box = QWidget()
        sc = QVBoxLayout(slice_box)
        sc.setContentsMargins(0, 0, 0, 0)
        self.slice_viewer = Slice2DViewer()
        sc.addWidget(self.slice_viewer)
        sc.addWidget(make_popout_button(self.slice_viewer, '2D slice viewer'))
        h.addWidget(slice_box, 1)

        voxel_box = QWidget()
        vc = QVBoxLayout(voxel_box)
        vc.setContentsMargins(0, 0, 0, 0)
        self.voxel3d_viewer = Voxel3DViewer()
        vc.addWidget(self.voxel3d_viewer)
        vc.addWidget(make_popout_button(self.voxel3d_viewer, '3D viewer'))
        h.addWidget(voxel_box, 1)

        v.addWidget(body, stretch=1)
        self.setCentralWidget(central)

        if self._items:
            self._show_item(0)

    # ── Item switching ───────────────────────────────────────────────────

    def _on_select(self, idx: int):
        if 0 <= idx < len(self._items):
            self._show_item(idx)

    def _show_item(self, idx: int):
        item = self._items[idx]
        vox = item.get('voxelgram')
        pitch = float(item.get('pitch_A', 1.0))
        if vox is None or pitch <= 0:
            return
        try:
            self.slice_viewer.set_voxelgram(vox, pitch)
            self.voxel3d_viewer.set_voxelgram(vox, pitch)
        except Exception as exc:
            print(f"[VoxelViewerWindow] failed to display item {idx}: {exc}")

    # ── Lifecycle ────────────────────────────────────────────────────────

    def closeEvent(self, evt):
        # VTK / Cocoa cleanup — same pattern as SaxsMorphPanel.
        try:
            self.voxel3d_viewer.shutdown()
        except Exception:
            pass
        super().closeEvent(evt)
