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
        QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QSlider,
        QDialog, QComboBox, QColorDialog, QFileDialog, QMenu, QSizePolicy,
        QApplication,
    )
    from PySide6.QtCore import Qt, Signal
    from PySide6.QtGui import QAction, QColor
except ImportError:
    try:
        from PyQt6.QtWidgets import (
            QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QSlider,
            QDialog, QComboBox, QColorDialog, QFileDialog, QMenu, QSizePolicy,
            QApplication,
        )
        from PyQt6.QtCore import Qt, pyqtSignal as Signal
        from PyQt6.QtGui import QAction, QColor
    except ImportError:
        from PyQt5.QtWidgets import (
            QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QSlider,
            QDialog, QComboBox, QColorDialog, QFileDialog, QMenu, QSizePolicy,
            QApplication,
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
        self._outline_actor = None
        self._show_outline = True
        self._iso_color = (0.10, 0.50, 0.90)
        self._voxelgram = None
        self._pitch_A = 1.0
        self._mesh = None

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
        a_reset = QAction('Reset view', self)
        a_reset.triggered.connect(self.reset_view)
        menu.addAction(a_reset)

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
        if self._outline_actor is not None:
            try:
                self.plotter.remove_actor(self._outline_actor)
            except Exception:
                pass

        self._actor = self.plotter.add_mesh(
            mesh, color=self._iso_color, smooth_shading=True,
            specular=0.3, specular_power=15,
        )
        if self._show_outline:
            self._outline_actor = self.plotter.add_mesh(
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
        if self._outline_actor is not None:
            try:
                self.plotter.remove_actor(self._outline_actor)
            except Exception:
                pass
            self._outline_actor = None
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
        self._outline_actor = None
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

        # Image view: strictly binary black/white.  The 2D slice shows
        # the PHYSICAL microstructure — thresholded 0 (void) and 1 (solid)
        # as sharp black and white regions.  Smoothing is deliberately
        # NOT applied here; it goes only to the 3D isosurface viewer.
        self.image_view = pg.ImageView(view=pg.PlotItem())
        self.image_view.ui.histogram.hide()
        self.image_view.ui.roiBtn.hide()
        self.image_view.ui.menuBtn.hide()
        self.image_view.view.setAspectLocked(True)
        self.image_view.view.invertY(True)
        bw_lut = np.array([[0, 0, 0, 255], [255, 255, 255, 255]], dtype=np.uint8)
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
        self._refresh()

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
            coord_label = f'z = {idx * self._pitch_A:.2f} Å'
        elif self._axis == 'y':
            slc = self._voxelgram[:, idx, :]
            coord_label = f'y = {idx * self._pitch_A:.2f} Å'
        else:
            slc = self._voxelgram[idx, :, :]
            coord_label = f'x = {idx * self._pitch_A:.2f} Å'

        self.image_view.setImage(slc.astype(np.float32), autoLevels=False, levels=[0, 1])
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
