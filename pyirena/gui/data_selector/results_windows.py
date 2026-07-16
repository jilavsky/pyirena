"""
pyirena.gui.data_selector.results_windows — stored-fit viewer windows and the tabulate window.

Split from the original monolithic data_selector.py (no behavior change).
"""

import logging
import os
from pathlib import Path
from typing import List

import numpy as np
import pyqtgraph as pg

log = logging.getLogger(__name__)

from pyirena.gui.data_selector._qt import (
    QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QFileDialog, QAbstractItemView, QTableWidget, QTableWidgetItem, Qt,
)
from pyirena.gui.sas_plot import add_slope_line_menu
from pyirena.io.hdf5 import readGenericNXcanSAS
from pyirena.io.text_import import ensure_nxcansas_sibling
from pyirena.io.nxcansas_unified import load_unified_fit_results

from pyirena.gui.data_selector.plot_utils import _LogDecadeAxis, _add_jpeg_export, _gen_colors, _iq_error_bars, _legend_indices, _rescaled_view, _style_plot


class GraphWindow(QWidget):
    """
    Separate window for displaying raw SAS data.
    Uses pyqtgraph for fast, interactive rendering even with large datasets.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("pyIrena - Data Viewer")
        self.setGeometry(100, 100, 850, 620)
        self._itx_technique = 'Data'   # ITX export → root:Data: (multi-file viewer)

        self.gl = pg.GraphicsLayoutWidget()
        self.gl.setBackground('w')

        self.plot = self.gl.addPlot(
            row=0, col=0,
            axisItems={
                'left':   _LogDecadeAxis(orientation='left'),
                'bottom': _LogDecadeAxis(orientation='bottom'),
            }
        )
        self.plot.setLogMode(True, True)
        self.plot.setLabel('left', 'Intensity (cm⁻¹)')
        self.plot.setLabel('bottom', 'Q (Å⁻¹)')
        self.plot.setTitle('Small-Angle Scattering Data', size='13pt')
        self.plot.showGrid(x=True, y=True, alpha=0.3)
        self.plot.addLegend(offset=(-10, 10), labelTextSize='18pt', labelTextColor='k')
        _style_plot(self.plot)
        _add_jpeg_export(self, self.plot)
        add_slope_line_menu(self.plot)

        # ── Per-session state for right-click toggles ──────────────────────
        self._plot_cache: list = []          # list of dicts, one per file
        self._data_items: list = []          # scatter / line PlotDataItems
        self._error_bar_items: list = []     # error bar PlotDataItems
        self._show_errorbars: bool = True
        self._line_mode: bool = False

        vb = self.plot.getViewBox()
        vb.menu.addSeparator()
        self._errbar_action = vb.menu.addAction("Hide Error Bars")
        self._errbar_action.triggered.connect(self._toggle_errorbars)
        self._linemode_action = vb.menu.addAction("Switch to Line Mode")
        self._linemode_action.triggered.connect(self._toggle_linemode)

        layout = QVBoxLayout()
        layout.addWidget(self.gl)
        self.setLayout(layout)

    def plot_data(
        self,
        file_paths: List[str],
        error_fraction: float = 0.05,
        max_legend_items: int = 12,
    ):
        """
        Plot data from the selected files.

        Args:
            file_paths:       List of file paths to plot.
            error_fraction:   Fraction used to synthesise uncertainty for 2-column
                              text files.
            max_legend_items: Maximum number of entries shown in the legend.
                              When there are more files, only first, last, and
                              evenly-spaced intermediate files are labelled.
        """
        self._plot_cache = []
        legend_idx = _legend_indices(len(file_paths), max_legend_items)
        colors = _gen_colors(len(file_paths))

        for idx, file_path in enumerate(file_paths):
            color = colors[idx]
            try:
                path, filename = os.path.split(file_path)
                _, ext = os.path.splitext(filename)

                if ext.lower() in ('.txt', '.dat'):
                    # Convert-once: creates/reuses a cleaned NXcanSAS sibling
                    h5_path = ensure_nxcansas_sibling(
                        Path(file_path),
                        error_fraction=error_fraction,
                    )
                    data = readGenericNXcanSAS(str(h5_path.parent), h5_path.name)
                else:
                    data = readGenericNXcanSAS(path, filename)

                if data is None:
                    continue

                q   = np.asarray(data['Q'],        dtype=float)
                I   = np.asarray(data['Intensity'], dtype=float)
                err = data.get('Error')
                if err is not None:
                    err = np.asarray(err, dtype=float)

                name = os.path.basename(file_path) if idx in legend_idx else None
                self._plot_cache.append({
                    'q': q, 'I': I, 'err': err, 'color': color, 'name': name,
                })

            except Exception as e:
                log.warning("Error loading %s: %s", file_path, e)

        self._redraw_items()

        # Set x-axis to the actual Q range of all loaded data so data from
        # instruments with very different Q ranges is always in view.
        # Use raw Q values when the user has switched to linear mode via the
        # right-click Transform menu; otherwise use log10 (ViewBox coords are
        # log10 when log mode is active).
        all_q = []
        for entry in self._plot_cache:
            q = entry['q']
            mask = np.isfinite(q) & (q > 0)
            if mask.any():
                all_q.append(q[mask])
        if all_q:
            q_all = np.concatenate(all_q)
            if len(q_all) >= 2:
                q_min, q_max = float(q_all.min()), float(q_all.max())
                vb = self.plot.getViewBox()
                try:
                    x_log = vb.state['logMode'][0]
                except (KeyError, IndexError, TypeError):
                    x_log = True  # safe default: log mode
                self.plot.setXRange(
                    np.log10(q_min) if x_log else q_min,
                    np.log10(q_max) if x_log else q_max,
                    padding=0.05,
                )

        self.show()

    def _redraw_items(self):
        """Rebuild scatter/line and error-bar items from the cached data.

        Called on initial load and whenever a toggle (error bars / line mode)
        changes the display style.  Clears only the managed data items so that
        the axes, grid, legend widget, and JPEG action are preserved.
        """
        for item in self._data_items + self._error_bar_items:
            self.plot.removeItem(item)
        self._data_items = []
        self._error_bar_items = []
        if self.plot.legend is not None:
            self.plot.legend.clear()

        for entry in self._plot_cache:
            q, I, err, color, name = (
                entry['q'], entry['I'], entry['err'],
                entry['color'], entry['name'],
            )
            mask = np.isfinite(q) & np.isfinite(I) & (q > 0) & (I > 0)
            q_, I_ = q[mask], I[mask]
            if len(q_) == 0:
                continue

            if self._line_mode:
                item = self.plot.plot(
                    q_, I_, pen=pg.mkPen(color, width=1.5), name=name,
                )
            else:
                item = self.plot.plot(
                    q_, I_,
                    pen=None, symbol='o', symbolSize=4,
                    symbolPen=pg.mkPen(color, width=1),
                    symbolBrush=pg.mkBrush(color),
                    name=name,
                )
            self._data_items.append(item)

            if name is not None and self.plot.legend is not None and self.plot.legend.items:
                self.plot.legend.items[-1][1].setAttr('color', color.name())

            if err is not None and self._show_errorbars:
                xb, yb = _iq_error_bars(q_, I_, err[mask] if err.shape == q.shape else err)
                if len(xb):
                    eb = self.plot.plot(
                        xb, yb,
                        pen=pg.mkPen(color, width=1),
                        connect='finite',
                    )
                    self._error_bar_items.append(eb)

    def _toggle_errorbars(self):
        """Right-click handler: show or hide error bars."""
        self._show_errorbars = not self._show_errorbars
        self._errbar_action.setText(
            "Show Error Bars" if not self._show_errorbars else "Hide Error Bars"
        )
        self._redraw_items()

    def _toggle_linemode(self):
        """Right-click handler: switch between point and line display."""
        self._line_mode = not self._line_mode
        self._linemode_action.setText(
            "Switch to Point Mode" if self._line_mode else "Switch to Line Mode"
        )
        self._redraw_items()


class UnifiedFitResultsWindow(QWidget):
    """
    Separate window for displaying Unified Fit results stored in HDF5 files.

    Two pyqtgraph panels: I(Q) data + model fit (top), normalised residuals
    (bottom).  X-axes are linked so zoom/pan stays synchronised.
    Files that contain no unified fit group are silently skipped.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("pyIrena - Unified Fit Results")
        self.setGeometry(130, 130, 900, 700)
        self._itx_technique = 'UnifiedFit'   # ITX export → root:UnifiedFit: (batch viewer)

        self.gl = pg.GraphicsLayoutWidget()
        self.gl.setBackground('w')

        # ── Top panel: data + model ────────────────────────────────────────
        self.ax_main = self.gl.addPlot(
            row=0, col=0,
            axisItems={
                'left':   _LogDecadeAxis(orientation='left'),
                'bottom': _LogDecadeAxis(orientation='bottom'),
            }
        )
        self.ax_main.setLogMode(True, True)
        self.ax_main.setLabel('left', 'Intensity (cm⁻¹)')
        self.ax_main.setTitle('Unified Fit Results', size='13pt')
        self.ax_main.showGrid(x=True, y=True, alpha=0.3)
        self.ax_main.addLegend(offset=(-10, 10), labelTextSize='18pt', labelTextColor='k')
        _style_plot(self.ax_main)
        add_slope_line_menu(self.ax_main)

        # ── Bottom panel: residuals ────────────────────────────────────────
        self.ax_resid = self.gl.addPlot(
            row=1, col=0,
            axisItems={'bottom': _LogDecadeAxis(orientation='bottom')},
        )
        self.ax_resid.setLogMode(True, False)
        self.ax_resid.setLabel('bottom', 'Q (Å⁻¹)')
        self.ax_resid.setLabel('left', "Residuals r' (rescaled)")
        self.ax_resid.showGrid(x=True, y=True, alpha=0.3)
        self.ax_resid.addLegend(offset=(-10, 10), labelTextSize='18pt', labelTextColor='k')
        _style_plot(self.ax_resid)

        # Synchronise x-axes
        self.ax_resid.setXLink(self.ax_main)

        # Allocate more vertical space to the main panel
        self.gl.ci.layout.setRowStretchFactor(0, 3)
        self.gl.ci.layout.setRowStretchFactor(1, 1)
        _add_jpeg_export(self, self.ax_main, self.ax_resid)

        layout = QVBoxLayout()
        layout.addWidget(self.gl)
        self.setLayout(layout)

    def plot_results(self, file_paths: List[str], max_legend_items: int = 12):
        """
        Load and plot Unified Fit results from the given file paths.

        Only HDF5 files are considered.  Files without a unified fit group
        are skipped silently.

        Args:
            file_paths:       List of file paths to load.
            max_legend_items: Maximum number of files shown in the legend.
        """
        self.ax_main.clear()
        self.ax_resid.clear()
        if self.ax_main.legend is not None:
            self.ax_main.legend.clear()
        if self.ax_resid.legend is not None:
            self.ax_resid.legend.clear()

        # Restore the zero-line after clear()
        self.ax_resid.addLine(y=0, pen=pg.mkPen('k', width=1))

        found_any = False
        colors = _gen_colors(len(file_paths))
        legend_idx = _legend_indices(len(file_paths), max_legend_items)

        for idx, file_path in enumerate(file_paths):
            _, ext = os.path.splitext(file_path)
            if ext.lower() not in ('.h5', '.hdf5', '.hdf'):
                continue

            try:
                results = load_unified_fit_results(Path(file_path))
            except Exception:
                continue

            color = colors[idx]
            label = os.path.basename(file_path)
            in_legend = idx in legend_idx
            Q         = results['Q']
            I_data    = results['intensity_data']
            I_model   = results['intensity_model']
            residuals = results['residuals']
            I_error   = results.get('intensity_error')
            chi2      = results.get('chi_squared', float('nan'))

            # ── data points ────────────────────────────────────────────────
            data_name = f'{label}  data' if in_legend else None
            self.ax_main.plot(
                Q, I_data,
                pen=None, symbol='o', symbolSize=4,
                symbolPen=pg.mkPen(color, width=1),
                symbolBrush=pg.mkBrush(color),
                name=data_name,
            )
            if data_name is not None and self.ax_main.legend is not None and self.ax_main.legend.items:
                self.ax_main.legend.items[-1][1].setAttr('color', color.name())

            # ── error bars ─────────────────────────────────────────────────
            if I_error is not None:
                xb, yb = _iq_error_bars(Q, I_data, I_error)
                if len(xb):
                    self.ax_main.plot(xb, yb,
                                      pen=pg.mkPen(color, width=1),
                                      connect='finite')

            # ── model line ─────────────────────────────────────────────────
            fit_color = pg.mkColor(color)
            fit_color.setAlpha(210)
            fit_name = f'{label}  fit  χ²={chi2:.3f}' if in_legend else None
            self.ax_main.plot(
                Q, I_model,
                pen=pg.mkPen(fit_color, width=2.5),
                name=fit_name,
            )
            if fit_name is not None and self.ax_main.legend is not None and self.ax_main.legend.items:
                self.ax_main.legend.items[-1][1].setAttr('color', fit_color.name())

            # ── residuals (rescaled, robust — uniform with live panels) ──────
            resid_name = label if in_legend else None
            self.ax_resid.plot(
                Q, _rescaled_view(residuals),
                pen=None, symbol='o', symbolSize=3,
                symbolPen=pg.mkPen(color, width=1),
                symbolBrush=pg.mkBrush(color),
                name=resid_name,
            )
            if resid_name is not None and self.ax_resid.legend is not None and self.ax_resid.legend.items:
                self.ax_resid.legend.items[-1][1].setAttr('color', color.name())

            found_any = True

        if not found_any:
            self.ax_main.setTitle(
                'No Unified Fit results found in selected files',
                size='12pt', color='#7f8c8d',
            )

        self.show()


class SizeDistResultsWindow(QWidget):
    """
    Separate window for displaying Size Distribution fit results stored in HDF5 files.

    Three pyqtgraph panels (all x-axes linked):
      top    — I(Q) data + model fit (log-log)
      middle — normalised residuals vs Q (log x)
      bottom — P(r) size distribution vs r (log x)
    Files that contain no sizes_results group are silently skipped.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("pyIrena - Size Distribution Results")
        self.setGeometry(150, 150, 900, 900)
        self._itx_technique = 'Sizes'   # ITX export → root:Sizes: (batch viewer)

        self.gl = pg.GraphicsLayoutWidget()
        self.gl.setBackground('w')

        # ── Top: I(Q) data + model ─────────────────────────────────────────
        self.ax_main = self.gl.addPlot(
            row=0, col=0,
            axisItems={
                'left':   _LogDecadeAxis(orientation='left'),
                'bottom': _LogDecadeAxis(orientation='bottom'),
            }
        )
        self.ax_main.setLogMode(True, True)
        self.ax_main.setLabel('left', 'Intensity (cm⁻¹)')
        self.ax_main.setTitle('Size Distribution Fit Results', size='13pt')
        self.ax_main.showGrid(x=True, y=True, alpha=0.3)
        self.ax_main.addLegend(offset=(-10, 10), labelTextSize='18pt', labelTextColor='k')
        _style_plot(self.ax_main)
        add_slope_line_menu(self.ax_main)

        # ── Middle: residuals ──────────────────────────────────────────────
        self.ax_resid = self.gl.addPlot(
            row=1, col=0,
            axisItems={'bottom': _LogDecadeAxis(orientation='bottom')},
        )
        self.ax_resid.setLogMode(True, False)
        self.ax_resid.setLabel('bottom', 'Q (Å⁻¹)')
        self.ax_resid.setLabel('left', "Residuals r' (rescaled)")
        self.ax_resid.showGrid(x=True, y=True, alpha=0.3)
        _style_plot(self.ax_resid)
        self.ax_resid.setXLink(self.ax_main)

        # ── Bottom: P(r) distribution ──────────────────────────────────────
        self.ax_dist = self.gl.addPlot(
            row=2, col=0,
            axisItems={'bottom': _LogDecadeAxis(orientation='bottom')},
        )
        self.ax_dist.setLogMode(True, False)
        self.ax_dist.setLabel('bottom', 'r (Å)')
        self.ax_dist.setLabel('left', 'P(r)  (vol. frac. Å⁻¹)')
        self.ax_dist.showGrid(x=True, y=True, alpha=0.3)
        self.ax_dist.addLegend(offset=(-10, 10), labelTextSize='18pt', labelTextColor='k')
        _style_plot(self.ax_dist)

        # Allocate space: main 3×, residuals 1×, distribution 2×
        self.gl.ci.layout.setRowStretchFactor(0, 3)
        self.gl.ci.layout.setRowStretchFactor(1, 1)
        self.gl.ci.layout.setRowStretchFactor(2, 2)
        _add_jpeg_export(self, self.ax_main, self.ax_resid, self.ax_dist)

        layout = QVBoxLayout()
        layout.addWidget(self.gl)
        self.setLayout(layout)

    def plot_results(self, file_paths: List[str], max_legend_items: int = 12):
        """
        Load and plot Size Distribution results from the given file paths.

        Only HDF5 files are considered.  Files without a sizes_results group
        are skipped silently.

        Args:
            file_paths:       List of file paths to load.
            max_legend_items: Maximum number of files shown in the legend.
        """
        from pyirena.io.nxcansas_sizes import load_sizes_results

        self.ax_main.clear()
        self.ax_resid.clear()
        self.ax_dist.clear()
        for ax in (self.ax_main, self.ax_resid, self.ax_dist):
            if ax.legend is not None:
                ax.legend.clear()

        self.ax_resid.addLine(y=0, pen=pg.mkPen('k', width=1))

        found_any = False
        colors = _gen_colors(len(file_paths))
        legend_idx = _legend_indices(len(file_paths), max_legend_items)

        for idx, file_path in enumerate(file_paths):
            _, ext = os.path.splitext(file_path)
            if ext.lower() not in ('.h5', '.hdf5', '.hdf'):
                continue

            try:
                results = load_sizes_results(Path(file_path))
            except Exception:
                continue

            color = colors[idx]
            label = os.path.basename(file_path)
            in_legend = idx in legend_idx
            Q         = results.get('Q')
            I_data    = results.get('intensity_data')
            I_model   = results.get('intensity_model')
            I_error   = results.get('intensity_error')
            residuals = results.get('residuals')
            r_grid    = results.get('r_grid')
            dist      = results.get('distribution')
            dist_std  = results.get('distribution_std')
            chi2      = results.get('chi_squared', float('nan'))
            vf        = results.get('volume_fraction', float('nan'))

            if Q is None or I_data is None or I_model is None:
                continue

            # ── data points ────────────────────────────────────────────────
            data_name = f'{label}  data' if in_legend else None
            self.ax_main.plot(
                Q, I_data,
                pen=None, symbol='o', symbolSize=4,
                symbolPen=pg.mkPen(color, width=1),
                symbolBrush=pg.mkBrush(color),
                name=data_name,
            )
            if data_name is not None and self.ax_main.legend is not None and self.ax_main.legend.items:
                self.ax_main.legend.items[-1][1].setAttr('color', color.name())

            # ── error bars ─────────────────────────────────────────────────
            if I_error is not None:
                xb, yb = _iq_error_bars(Q, I_data, I_error)
                if len(xb):
                    self.ax_main.plot(xb, yb,
                                      pen=pg.mkPen(color, width=1),
                                      connect='finite')

            # ── model line ─────────────────────────────────────────────────
            fit_color = pg.mkColor(color)
            fit_color.setAlpha(210)
            vf_str = f'{vf:.3g}' if (vf == vf) else 'N/A'
            fit_name = f'{label}  fit  χ²={chi2:.3f}' if in_legend else None
            self.ax_main.plot(
                Q, I_model,
                pen=pg.mkPen(fit_color, width=2.5),
                name=fit_name,
            )
            if fit_name is not None and self.ax_main.legend is not None and self.ax_main.legend.items:
                self.ax_main.legend.items[-1][1].setAttr('color', fit_color.name())

            # ── residuals (rescaled, robust — uniform with live panels) ──────
            if residuals is not None:
                self.ax_resid.plot(
                    Q, _rescaled_view(residuals),
                    pen=None, symbol='o', symbolSize=3,
                    symbolPen=pg.mkPen(color, width=1),
                    symbolBrush=pg.mkBrush(color),
                )

            # ── size distribution ───────────────────────────────────────────
            if r_grid is not None and dist is not None:
                # ±1σ band as dashed upper/lower bounds
                if dist_std is not None:
                    band_color = pg.mkColor(color)
                    band_color.setAlpha(120)
                    band_pen = pg.mkPen(band_color, width=1,
                                        style=Qt.PenStyle.DashLine)
                    upper = np.maximum(dist + dist_std, 0.0)
                    lower = np.maximum(dist - dist_std, 0.0)
                    self.ax_dist.plot(r_grid, upper, pen=band_pen)
                    self.ax_dist.plot(r_grid, lower, pen=band_pen)

                dist_name = f'{label}  Vf={vf_str}' if in_legend else None
                self.ax_dist.plot(
                    r_grid, dist,
                    pen=pg.mkPen(color, width=2),
                    name=dist_name,
                )
                if dist_name is not None and self.ax_dist.legend is not None and self.ax_dist.legend.items:
                    self.ax_dist.legend.items[-1][1].setAttr('color', color.name())

            found_any = True

        if not found_any:
            self.ax_main.setTitle(
                'No Size Distribution results found in selected files',
                size='12pt', color='#7f8c8d',
            )

        self.show()


class SimpleFitResultsWindow(QWidget):
    """
    Separate window for displaying Simple Fits results stored in HDF5 files.

    Two pyqtgraph panels (x-axes linked):
      top    — I(Q) data + model fit (log-log)
      bottom — normalised residuals vs Q (log x)
    Files that contain no simple_fit_results group are silently skipped.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("pyIrena - Simple Fits Results")
        self.setGeometry(140, 140, 900, 700)
        self._itx_technique = 'SimpleFits'   # ITX export → root:SimpleFits: (batch viewer)

        self.gl = pg.GraphicsLayoutWidget()
        self.gl.setBackground('w')

        # ── Top: I(Q) data + model ─────────────────────────────────────────
        self.ax_main = self.gl.addPlot(
            row=0, col=0,
            axisItems={
                'left':   _LogDecadeAxis(orientation='left'),
                'bottom': _LogDecadeAxis(orientation='bottom'),
            }
        )
        self.ax_main.setLogMode(True, True)
        self.ax_main.setLabel('left', 'Intensity (cm⁻¹)')
        self.ax_main.setTitle('Simple Fits Results', size='13pt')
        self.ax_main.showGrid(x=True, y=True, alpha=0.3)
        self.ax_main.addLegend(offset=(-10, 10), labelTextSize='18pt', labelTextColor='k')
        _style_plot(self.ax_main)
        add_slope_line_menu(self.ax_main)

        # ── Bottom: residuals ──────────────────────────────────────────────
        self.ax_resid = self.gl.addPlot(
            row=1, col=0,
            axisItems={'bottom': _LogDecadeAxis(orientation='bottom')},
        )
        self.ax_resid.setLogMode(True, False)
        self.ax_resid.setLabel('bottom', 'Q (Å⁻¹)')
        self.ax_resid.setLabel('left', "Residuals r' (rescaled)")
        self.ax_resid.showGrid(x=True, y=True, alpha=0.3)
        _style_plot(self.ax_resid)
        self.ax_resid.setXLink(self.ax_main)

        self.gl.ci.layout.setRowStretchFactor(0, 3)
        self.gl.ci.layout.setRowStretchFactor(1, 1)
        _add_jpeg_export(self, self.ax_main, self.ax_resid)

        layout = QVBoxLayout()
        layout.addWidget(self.gl)
        self.setLayout(layout)

    def plot_results(self, file_paths: List[str], max_legend_items: int = 12):
        """
        Load and plot Simple Fits results from the given file paths.

        Only HDF5 files are considered.  Files without a simple_fit_results
        group are skipped silently.
        """
        from pyirena.io.nxcansas_simple_fits import load_simple_fit_results

        self.ax_main.clear()
        self.ax_resid.clear()
        for ax in (self.ax_main, self.ax_resid):
            if ax.legend is not None:
                ax.legend.clear()

        self.ax_resid.addLine(y=0, pen=pg.mkPen('k', width=1))

        found_any = False
        colors = _gen_colors(len(file_paths))
        legend_idx = _legend_indices(len(file_paths), max_legend_items)

        for idx, file_path in enumerate(file_paths):
            _, ext = os.path.splitext(file_path)
            if ext.lower() not in ('.h5', '.hdf5', '.hdf'):
                continue

            try:
                results = load_simple_fit_results(Path(file_path))
            except Exception:
                continue

            color     = colors[idx]
            label     = os.path.basename(file_path)
            in_legend = idx in legend_idx

            Q         = results.get('Q')
            I_data    = results.get('intensity_data')
            I_model   = results.get('I_model')
            I_error   = results.get('intensity_error')
            residuals = results.get('residuals')
            chi2      = results.get('chi_squared', float('nan'))
            model     = results.get('model', '')

            if Q is None or I_model is None:
                continue

            # ── data points ────────────────────────────────────────────────
            if I_data is not None:
                data_name = f'{label}  data' if in_legend else None
                self.ax_main.plot(
                    Q, I_data,
                    pen=None, symbol='o', symbolSize=4,
                    symbolPen=pg.mkPen(color, width=1),
                    symbolBrush=pg.mkBrush(color),
                    name=data_name,
                )
                if data_name is not None and self.ax_main.legend is not None and self.ax_main.legend.items:
                    self.ax_main.legend.items[-1][1].setAttr('color', color.name())

                # ── error bars ─────────────────────────────────────────────
                if I_error is not None:
                    xb, yb = _iq_error_bars(Q, I_data, I_error)
                    if len(xb):
                        self.ax_main.plot(xb, yb,
                                          pen=pg.mkPen(color, width=1),
                                          connect='finite')

            # ── model line — darker shade for clear contrast with data symbols
            fit_color = color.darker(280)   # ~36% brightness of the data colour
            derived = results.get('derived', {}) or {}
            if model == 'Invariant':
                # Calculation model: no fitted curve.  Show the background-
                # corrected intensity (same units as the data) and put the
                # invariant / volume fraction into the legend instead of χ².
                I_line = results.get('I_corrected')
                if I_line is None:
                    I_line = I_model
                phi = derived.get('VolumeFraction', float('nan'))
                inv = derived.get('Invariant', float('nan'))
                fit_name = (f'{label}  Invariant Q*={inv:.3g} cm⁻⁴  φ={phi:.3g}'
                            if in_legend else None)
                _mask = np.isfinite(I_line) & (np.asarray(I_line) > 0)
                self.ax_main.plot(
                    np.asarray(Q)[_mask], np.asarray(I_line)[_mask],
                    pen=pg.mkPen(fit_color, width=3.0),
                    name=fit_name,
                )
            else:
                chi2_str = f'{chi2:.3f}' if (chi2 == chi2) else 'N/A'
                fit_name = f'{label}  {model}  χ²={chi2_str}' if in_legend else None
                self.ax_main.plot(
                    Q, I_model,
                    pen=pg.mkPen(fit_color, width=3.0),
                    name=fit_name,
                )
            if fit_name is not None and self.ax_main.legend is not None and self.ax_main.legend.items:
                self.ax_main.legend.items[-1][1].setAttr('color', fit_color.name())

            # ── residuals (rescaled, robust — uniform with live panels) ──────
            if residuals is not None:
                self.ax_resid.plot(
                    Q, _rescaled_view(residuals),
                    pen=None, symbol='o', symbolSize=3,
                    symbolPen=pg.mkPen(color, width=1),
                    symbolBrush=pg.mkBrush(color),
                )

            found_any = True

        if not found_any:
            self.ax_main.setTitle(
                'No Simple Fits results found in selected files',
                size='12pt', color='#7f8c8d',
            )

        self.show()


class WAXSPeakFitResultsWindow(QWidget):
    """
    Separate window for displaying WAXS Peak Fit results stored in HDF5 files.

    Two pyqtgraph panels (x-axes linked), both LINEAR/LINEAR:
      top    — I(Q) data + model fit
      bottom — normalised residuals vs Q
    Files that contain no waxs_peakfit_results group are silently skipped.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("pyIrena - WAXS Peak Fit Results")
        self.setGeometry(145, 145, 900, 700)
        self._itx_technique = 'WAXSPeakFit'   # ITX export → root:WAXSPeakFit: (batch viewer)

        self.gl = pg.GraphicsLayoutWidget()
        self.gl.setBackground('w')

        # ── Top: data + model ──────────────────────────────────────────────
        self.ax_main = self.gl.addPlot(row=0, col=0)
        self.ax_main.setLogMode(False, False)     # linear / linear
        self.ax_main.setLabel('left', 'Intensity')
        self.ax_main.setLabel('bottom', 'Q  (Å⁻¹)')
        self.ax_main.setTitle('WAXS Peak Fit Results', size='13pt')
        self.ax_main.showGrid(x=True, y=True, alpha=0.3)
        self.ax_main.addLegend(offset=(-10, 10), labelTextSize='10pt', labelTextColor='k')
        _style_plot(self.ax_main)

        # ── Bottom: residuals ──────────────────────────────────────────────
        self.ax_resid = self.gl.addPlot(row=1, col=0)
        self.ax_resid.setLogMode(False, False)
        self.ax_resid.setLabel('bottom', 'Q  (Å⁻¹)')
        self.ax_resid.setLabel('left', "Residuals r' (rescaled)")
        self.ax_resid.showGrid(x=True, y=True, alpha=0.3)
        self.ax_resid.setXLink(self.ax_main)
        _style_plot(self.ax_resid)

        self.gl.ci.layout.setRowStretchFactor(0, 3)
        self.gl.ci.layout.setRowStretchFactor(1, 1)
        _add_jpeg_export(self, self.ax_main, self.ax_resid)

        # ── Display toggle (points ↔ line) ─────────────────────────────────
        self._waxs_data_as_lines = False
        self._waxs_data_items: list = []
        vb = self.ax_main.getViewBox()
        vb.menu.addSeparator()
        self._waxs_data_mode_act = vb.menu.addAction("Show data as line")
        self._waxs_data_mode_act.triggered.connect(self._toggle_waxs_data_mode)

        layout = QVBoxLayout()
        layout.addWidget(self.gl)
        self.setLayout(layout)

    def _toggle_waxs_data_mode(self):
        self._waxs_data_as_lines = not self._waxs_data_as_lines
        for item, color in self._waxs_data_items:
            if self._waxs_data_as_lines:
                item.setData(pen=pg.mkPen(color, width=1.5), symbol=None)
            else:
                item.setData(pen=None, symbol='o', symbolSize=3,
                             symbolPen=pg.mkPen(color, width=1),
                             symbolBrush=pg.mkBrush(color))
        self._waxs_data_mode_act.setText(
            "Show data as points" if self._waxs_data_as_lines else "Show data as line"
        )

    def plot_results(self, file_paths: List[str], max_legend_items: int = 12):
        """Load and plot WAXS peak-fit results from the given file paths."""
        from pyirena.io.nxcansas_waxs_peakfit import load_waxs_peakfit_results

        self.ax_main.clear()
        self.ax_resid.clear()
        self._waxs_data_items.clear()
        self._waxs_data_as_lines = False
        self._waxs_data_mode_act.setText("Show data as line")
        for ax in (self.ax_main, self.ax_resid):
            if ax.legend is not None:
                ax.legend.clear()

        self.ax_resid.addLine(y=0, pen=pg.mkPen('k', width=1))

        found_any = False
        colors = _gen_colors(len(file_paths))
        legend_idx = _legend_indices(len(file_paths), max_legend_items)

        for idx, file_path in enumerate(file_paths):
            _, ext = os.path.splitext(file_path)
            if ext.lower() not in ('.h5', '.hdf5', '.hdf'):
                continue

            try:
                results = load_waxs_peakfit_results(Path(file_path))
            except Exception:
                continue

            color     = colors[idx]
            label     = os.path.basename(file_path)
            in_legend = idx in legend_idx
            Q         = results.get('Q')
            I_data    = results.get('intensity_data')
            I_fit     = results.get('I_fit')
            residuals = results.get('residuals')

            if Q is None or I_fit is None:
                continue

            name = label if in_legend else None

            # Data scatter (tracked for points ↔ line toggle)
            if I_data is not None:
                mask = np.isfinite(Q) & np.isfinite(I_data)
                di = self.ax_main.plot(
                    Q[mask], I_data[mask],
                    pen=None, symbol='o', symbolSize=3,
                    symbolPen=pg.mkPen(color, width=1),
                    symbolBrush=pg.mkBrush(color),
                    name=name,
                )
                self._waxs_data_items.append((di, color))

            # Model line — darker than the data colour so they are distinguishable
            mask_f = np.isfinite(Q) & np.isfinite(I_fit)
            model_col = pg.mkColor(color).darker(175)
            self.ax_main.plot(
                Q[mask_f], I_fit[mask_f],
                pen=pg.mkPen(model_col, width=2),
                name=(f"{label} fit") if in_legend else None,
            )

            # Residuals (rescaled, robust — uniform with live panels)
            if residuals is not None:
                residuals_r = _rescaled_view(residuals)
                mask_r = np.isfinite(Q) & np.isfinite(residuals_r)
                self.ax_resid.plot(
                    Q[mask_r], residuals_r[mask_r],
                    pen=None, symbol='o', symbolSize=3,
                    symbolPen=pg.mkPen(color, width=1),
                    symbolBrush=pg.mkBrush(color),
                )

            found_any = True

        if not found_any:
            self.ax_main.setTitle(
                'No WAXS Peak Fit results found in selected files',
                size='12pt', color='#7f8c8d',
            )

        self.show()


class TabulateResultsWindow(QWidget):
    """
    Separate window that shows fit results for selected files in a spreadsheet-like
    table and lets the user save them as a CSV file.

    One row per file.  Columns include Unified Fit scalars + per-level parameters
    and/or Size Distribution scalars, depending on which results were found.
    CSV is plain ASCII with comma separators — readable by Excel, Igor Pro,
    Origin, SigmaPlot, and any scripting language.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("pyIrena - Tabulated Results")
        self.setGeometry(160, 160, 1300, 500)

        self._headers: list = []
        self._rows: list    = []

        self.table = QTableWidget()
        self.table.setAlternatingRowColors(True)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.horizontalHeader().setStretchLastSection(False)

        self.save_btn = QPushButton("Save as CSV…")
        self.save_btn.setMinimumHeight(32)
        self.save_btn.setStyleSheet("""
            QPushButton {
                background-color: #16a085; color: white;
                font-size: 12px; font-weight: bold;
                border-radius: 4px; padding: 5px 12px;
            }
            QPushButton:hover { background-color: #138d75; }
            QPushButton:disabled { background-color: #bdc3c7; }
        """)
        self.save_btn.setToolTip("Save the current results table to a CSV file for use in spreadsheets.")
        self.save_btn.clicked.connect(self._save_csv)
        self.save_btn.setEnabled(False)

        btn_row = QHBoxLayout()
        btn_row.addStretch()
        btn_row.addWidget(self.save_btn)

        layout = QVBoxLayout()
        layout.addWidget(self.table)
        layout.addLayout(btn_row)
        self.setLayout(layout)

    def set_data(self, headers: list, rows: list, default_save_path: str = ''):
        """Populate the table and remember data for CSV export."""
        self._headers = headers
        self._rows    = rows
        self._default_save_path = default_save_path

        self.table.setColumnCount(len(headers))
        self.table.setRowCount(len(rows))
        self.table.setHorizontalHeaderLabels(headers)

        for r, row in enumerate(rows):
            for c, val in enumerate(row):
                text = '' if val is None else str(val)
                item = QTableWidgetItem(text)
                item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                self.table.setItem(r, c, item)

        self.table.resizeColumnsToContents()
        self.save_btn.setEnabled(bool(rows))
        self.show()

    def _save_csv(self):
        import csv as _csv
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Save tabulated results as CSV",
            self._default_save_path,
            "CSV files (*.csv);;All files (*)",
        )
        if not path:
            return
        if not path.lower().endswith('.csv'):
            path += '.csv'
        with open(path, 'w', newline='', encoding='utf-8') as fh:
            writer = _csv.writer(fh)
            writer.writerow(self._headers)
            writer.writerows(
                ['' if v is None else v for v in row]
                for row in self._rows
            )
