"""
pyirena.gui.data_selector.plot_utils — shared pyqtgraph helpers (palette, axes, styling, JPEG export).

Split from the original monolithic data_selector.py (no behavior change).
"""

from pathlib import Path

import numpy as np
import pyqtgraph as pg

from pyirena.gui.data_selector._qt import (
    QFileDialog, QAction,
)
from pyirena.gui.sas_plot import save_itx_from_plot



# ── Shared colour palette for multi-file graphs ────────────────────────────
def _gen_colors(n: int) -> list:
    """
    Return a list of *n* QColor objects evenly distributed from blue (index 0)
    to red (index n-1) along the HSV hue wheel, covering the full visible
    spectrum (blue → cyan → green → yellow → red).

    With a single file a neutral blue is returned.  With N files the colours
    are guaranteed to be unique — the palette never repeats regardless of how
    many datasets are plotted simultaneously, making it easy to track where
    features appear across an ordered series.
    """
    if n <= 0:
        return []
    if n == 1:
        return [pg.hsvColor(0.60, sat=0.85, val=0.82)]   # single medium blue
    return [
        pg.hsvColor(0.667 * (1.0 - i / (n - 1)), sat=0.88, val=0.82)
        for i in range(n)
    ]


def _legend_indices(n: int, max_items: int) -> set:
    """Return the set of file indices that should receive a legend entry.

    Always includes the first and last index.  Middle indices are thinned so
    that the total does not exceed *max_items*.
    """
    if n <= max_items:
        return set(range(n))
    if max_items <= 2:
        return {0, n - 1}
    step = max(1, (n - 2) // (max_items - 2))
    middle = list(range(1, n - 1, step))[:(max_items - 2)]
    return {0} | set(middle) | {n - 1}


class _LogDecadeAxis(pg.AxisItem):
    """
    Decade-only tick labels for log-mode plots.
    Minor ticks at ×2 (log10≈0.301) and ×5 (log10≈0.699) per decade.
    Identical to LogDecadeAxis in sizes_panel.py.
    """
    def tickValues(self, minVal, maxVal, size):
        if not self.logMode:
            return super().tickValues(minVal, maxVal, size)
        import math
        lo = math.floor(minVal - 0.001)
        hi = math.ceil(maxVal + 0.001)
        major = [float(v) for v in range(lo, hi + 1)
                 if minVal - 0.01 <= float(v) <= maxVal + 0.01]
        minor = []
        for decade in range(lo, hi + 1):
            for offset in (0.301, 0.699):
                v = float(decade) + offset
                if minVal - 0.01 <= v <= maxVal + 0.01:
                    minor.append(v)
        return [(1.0, major), (0.301, minor)]

    def tickStrings(self, values, scale, spacing):
        if not self.logMode:
            return super().tickStrings(values, scale, spacing)
        strings = []
        for v in values:
            if abs(v - round(v)) < 0.05:
                pwr = int(round(v))
                val = 10.0 ** pwr
                strings.append(f'{val:g}' if 0.001 <= abs(val) <= 999 else f'1e{pwr}')
            else:
                strings.append('')
        return strings


def _style_plot(plot_item):
    """White background, black axes, SI-prefix scaling disabled."""
    for side in ('left', 'bottom'):
        ax = plot_item.getAxis(side)
        ax.setPen(pg.mkPen('k'))
        ax.setTextPen(pg.mkPen('k'))
        ax.enableAutoSIPrefix(False)


def _iq_error_bars(q, I, err, cap_frac=0.02):
    """
    Build NaN-separated (x, y) arrays for I(Q) error bars on a log-log plot.
    Caps are ±cap_frac in multiplicative (log) space.
    """
    x_lines, y_lines = [], []
    for i in range(len(q)):
        if q[i] <= 0 or I[i] <= 0 or not (err[i] > 0):
            continue
        y_top = I[i] + err[i]
        y_bot = max(I[i] - err[i], I[i] * 0.001)
        x_lines.extend([q[i], q[i], np.nan])
        y_lines.extend([y_bot, y_top, np.nan])
        x_lines.extend([q[i] / (1 + cap_frac), q[i] * (1 + cap_frac), np.nan])
        y_lines.extend([y_top, y_top, np.nan])
        x_lines.extend([q[i] / (1 + cap_frac), q[i] * (1 + cap_frac), np.nan])
        y_lines.extend([y_bot, y_bot, np.nan])
    return np.array(x_lines, dtype=float), np.array(y_lines, dtype=float)


# ── JPEG export helper ─────────────────────────────────────────────────────
def _add_jpeg_export(window, *plot_items):
    """
    Add 'Save as JPEG' and 'Save as Igor Pro ITX' actions to the ViewBox
    context menu of every PlotItem passed in.  JPEG captures the whole
    *window* widget; ITX exports named curves from the individual plot_item.
    """
    def _save():
        path, _ = QFileDialog.getSaveFileName(
            window, "Save as JPEG",
            str(Path.home()),
            "JPEG images (*.jpg *.jpeg);;All files (*)",
        )
        if not path:
            return
        if not path.lower().endswith(('.jpg', '.jpeg')):
            path += '.jpg'
        window.grab().save(path, 'JPEG', 95)

    for plot_item in plot_items:
        act = QAction("Save as JPEG…", window)
        act.triggered.connect(_save)
        plot_item.getViewBox().menu.addAction(act)

        act_itx = QAction("Save as Igor Pro ITX…", window)
        act_itx.triggered.connect(
            lambda checked=False, p=plot_item, w=window:
                save_itx_from_plot(p, w)
        )
        plot_item.getViewBox().menu.addAction(act_itx)


def _rescaled_view(residuals):
    """Rescale stored residuals by their robust (MAD-based) scale for display.

    Uniform with the live fit panels: the viewer shows r' = r/s so scatter is
    compared to the data's own robust noise floor. Returns the input unchanged
    if it is None or has no finite points.
    """
    if residuals is None:
        return residuals
    from pyirena.core.fit_metrics import rescale_residuals
    r_prime, _ = rescale_residuals(np.asarray(residuals, dtype=float))
    return r_prime
