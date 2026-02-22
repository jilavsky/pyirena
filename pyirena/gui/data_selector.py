"""
Data Selector GUI for pyIrena.

This module provides a GUI panel for selecting data files and displaying
their content as graphs.
"""

import os
import sys
import subprocess
from pathlib import Path
from typing import List, Optional

try:
    from PySide6.QtWidgets import (
        QApplication, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QPushButton,
        QListWidget, QLabel, QLineEdit, QFileDialog, QComboBox,
        QAbstractItemView, QMessageBox, QMenuBar, QMenu,
        QDialog, QFormLayout, QDialogButtonBox, QGroupBox, QCheckBox, QColorDialog,
        QTableWidget, QTableWidgetItem,
    )
    from PySide6.QtCore import Qt, QDir, QThread, Signal
    from PySide6.QtGui import QAction, QDoubleValidator
except ImportError:
    try:
        from PyQt6.QtWidgets import (
            QApplication, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QPushButton,
            QListWidget, QLabel, QLineEdit, QFileDialog, QComboBox,
            QAbstractItemView, QMessageBox, QMenuBar, QMenu,
            QDialog, QFormLayout, QDialogButtonBox, QGroupBox, QCheckBox, QColorDialog,
            QTableWidget, QTableWidgetItem,
        )
        from PyQt6.QtCore import Qt, QDir, QThread, pyqtSignal as Signal
        from PyQt6.QtGui import QAction, QDoubleValidator
    except ImportError:
        raise ImportError(
            "Neither PySide6 nor PyQt6 found. Install with: pip install PySide6"
        )

import re
import numpy as np
import pyqtgraph as pg

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


def _iq_error_bars(q, I, err, cap_frac=0.05):
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


# ── Filename sort-key extractors ───────────────────────────────────────────
def _sort_key_name(name: str) -> str:
    return name.lower()

def _sort_key_temperature(name: str) -> float:
    m = re.search(r'_(-?\d+(?:\.\d+)?)C(?=_|\.|$)', name, re.IGNORECASE)
    return float(m.group(1)) if m else float('inf')

def _sort_key_time(name: str) -> float:
    m = re.search(r'_(\d+(?:\.\d+)?)min(?=_|\.|$)', name, re.IGNORECASE)
    return float(m.group(1)) if m else float('inf')

def _sort_key_order(name: str) -> float:
    # last _NNN before the file extension
    m = re.search(r'_(\d+)(?:\.[^.]+)?$', name)
    return float(m.group(1)) if m else float('inf')

def _sort_key_pressure(name: str) -> float:
    m = re.search(r'_(\d+(?:\.\d+)?)PSI(?=_|\.|$)', name, re.IGNORECASE)
    return float(m.group(1)) if m else float('inf')

_SORT_KEYS = [
    _sort_key_name,        # 0 Filename A→Z
    _sort_key_name,        # 1 Filename Z→A
    _sort_key_temperature, # 2 Temperature ↑
    _sort_key_temperature, # 3 Temperature ↓
    _sort_key_time,        # 4 Time ↑
    _sort_key_time,        # 5 Time ↓
    _sort_key_order,       # 6 Order number ↑
    _sort_key_order,       # 7 Order number ↓
    _sort_key_pressure,    # 8 Pressure ↑
    _sort_key_pressure,    # 9 Pressure ↓
]


# ── JPEG export helper ─────────────────────────────────────────────────────
def _add_jpeg_export(window, *plot_items):
    """
    Add a 'Save as JPEG…' action to the ViewBox context menu of every
    PlotItem passed in.  Captures the whole *window* widget via grab().
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


from pyirena.io.hdf5 import readGenericNXcanSAS, readTextFile
from pyirena.io.nxcansas_unified import load_unified_fit_results
from pyirena.gui.unified_fit import UnifiedFitPanel
from pyirena.gui.sizes_panel import SizesFitPanel
from pyirena.state import StateManager
from pyirena.batch import fit_unified, fit_sizes, fit_simple_from_config


def _build_report(file_path: str,
                  data_info: Optional[dict] = None,
                  fit_results: Optional[dict] = None,
                  sizes_results: Optional[dict] = None,
                  simple_fit_results: Optional[dict] = None) -> str:
    """
    Build a Markdown report string.

    Args:
        file_path:          Absolute path to the source file.
        data_info:          Dict with keys 'Q', 'I', 'I_error' (optional array).
                            Pass None to omit the data section.
        fit_results:        Dict from load_unified_fit_results().
                            Pass None to omit the unified fit section.
        sizes_results:      Dict from load_sizes_results().
                            Pass None to omit the size distribution section.
        simple_fit_results: Dict from load_simple_fit_results().
                            Pass None to omit the simple fits section.

    Returns:
        Multi-line Markdown string ready to be written to a .md file.
    """
    from datetime import datetime as _dt

    filename = os.path.basename(file_path)
    now = _dt.now().strftime('%Y-%m-%d %H:%M:%S')

    L = []

    # ── Header ───────────────────────────────────────────────────────────────
    L += [
        "# pyIrena Report",
        "",
        "| | |",
        "|---|---|",
        f"| **File** | `{filename}` |",
        f"| **Report generated** | {now} |",
    ]
    if fit_results is not None:
        L += [
            f"| **Unified Fit timestamp** | {fit_results.get('timestamp', 'unknown')} |",
        ]
    if sizes_results is not None:
        L += [
            f"| **Size Dist. timestamp** | {sizes_results.get('timestamp', 'unknown')} |",
        ]
    if simple_fit_results is not None:
        L += [
            f"| **Simple Fits timestamp** | {simple_fit_results.get('timestamp', 'unknown')} |",
        ]
    L.append("")

    # ── Data summary ─────────────────────────────────────────────────────────
    if data_info is not None:
        Q = data_info['Q']
        I = data_info['I']
        I_error = data_info.get('I_error')
        L += [
            "## Data Summary",
            "",
            "| Parameter | Value |",
            "|-----------|-------|",
            f"| Q range | {Q.min():.4g} – {Q.max():.4g} Å⁻¹ |",
            f"| Intensity range | {I.min():.4e} – {I.max():.4e} cm⁻¹ |",
            f"| Data points | {len(Q)} |",
        ]
        if I_error is not None:
            L.append(
                f"| Uncertainty range | {I_error.min():.4e} – {I_error.max():.4e} cm⁻¹ |"
            )
        L.append("")

    # ── Fit quality ──────────────────────────────────────────────────────────
    if fit_results is not None:
        chi2      = fit_results['chi_squared']
        bg        = fit_results['background']
        n_levels  = fit_results['num_levels']
        Q_fit     = fit_results['Q']
        residuals = fit_results['residuals']
        levels    = fit_results.get('levels', [])

        L += [
            "## Fit Quality",
            "",
            "| Parameter | Value |",
            "|-----------|-------|",
            f"| Chi-squared (χ²) | {chi2:.4f} |",
            f"| Number of levels | {n_levels} |",
            f"| Background | {bg:.4e} cm⁻¹ |",
            f"| Q range (fit) | {Q_fit.min():.4g} – {Q_fit.max():.4g} Å⁻¹ |",
            f"| Data points (fit) | {len(Q_fit)} |",
            f"| Residuals mean | {np.mean(residuals):.4f} |",
            f"| Residuals std dev | {np.std(residuals):.4f} |",
            f"| Max \\|residual\\| | {np.max(np.abs(residuals)):.4f} |",
            "",
        ]

        # ── Level parameters ──────────────────────────────────────────────
        for i, level in enumerate(levels):
            lnum = i + 1
            L.append(f"## Level {lnum} Parameters")
            L.append("")

            has_mc = any(f'{p}_err' in level for p in ('G', 'Rg', 'B', 'P', 'ETA', 'PACK'))
            if has_mc:
                L += ["| Parameter | Value | Uncertainty (1σ) |",
                      "|-----------|-------|------------------|"]
            else:
                L += ["| Parameter | Value |",
                      "|-----------|-------|"]

            def _row(label, key, unit='', fmt='.4e'):
                val = level.get(key)
                if val is None:
                    return
                val_str = f"{val:{fmt}}{unit}"
                if has_mc:
                    err = level.get(f'{key}_err', 0.0)
                    err_str = f"± {err:{fmt}}{unit}" if err > 0 else "—"
                    L.append(f"| {label} | {val_str} | {err_str} |")
                else:
                    L.append(f"| {label} | {val_str} |")

            _row('G',  'G')
            _row('Rg', 'Rg', ' Å')
            _row('B',  'B')
            _row('P',  'P',  fmt='.4f')

            rgcut = level.get('RgCutoff', 0.0)
            if isinstance(rgcut, float) and rgcut > 0.01:
                _row('RgCutoff', 'RgCutoff', ' Å')

            if level.get('correlated', False):
                _row('ETA',  'ETA',  ' Å', fmt='.2f')
                _row('PACK', 'PACK', '',   fmt='.4f')

            _row('Sv',        'Sv',        ' m²/cm³')
            _row('Invariant', 'Invariant', ' cm⁻⁴')

            L.append("")

    # ── Size Distribution results ─────────────────────────────────────────────
    if sizes_results is not None:
        # All scalar metadata is at the top level of the dict returned by
        # load_sizes_results() — not nested under a 'params' key.
        r_grid    = sizes_results.get('r_grid')
        dist      = sizes_results.get('distribution')
        residuals = sizes_results.get('residuals')

        def _sv(key, default=float('nan')):
            """Return scalar from sizes_results, substituting default for None."""
            v = sizes_results.get(key)
            return v if v is not None else default

        chi2        = _sv('chi_squared')
        vf          = _sv('volume_fraction')
        rg          = _sv('rg')
        n_iter      = _sv('n_iterations', 'N/A')
        method      = _sv('method', 'unknown')
        shape       = _sv('shape', 'unknown')
        contrast    = _sv('contrast')
        aspect_ratio = _sv('aspect_ratio')
        n_bins      = _sv('n_bins', 0)
        r_min       = _sv('r_min')
        r_max       = _sv('r_max')
        log_spacing = _sv('log_spacing', False)
        background  = _sv('background')
        error_scale = _sv('error_scale')
        power_law_B = _sv('power_law_B')
        power_law_P = _sv('power_law_P')
        q_fit_min   = sizes_results.get('cursor_q_min')
        q_fit_max   = sizes_results.get('cursor_q_max')

        peak_r = float('nan')
        if dist is not None and r_grid is not None and len(dist) > 0:
            peak_r = float(r_grid[int(np.argmax(dist))])

        def _fmt(v, spec='.4g', suffix=''):
            """Format a value, returning 'N/A' for None/nan."""
            if v is None:
                return 'N/A'
            try:
                if np.isnan(float(v)):
                    return 'N/A'
            except (TypeError, ValueError):
                pass
            try:
                return f"{v:{spec}}{suffix}"
            except (TypeError, ValueError):
                return str(v)

        L += [
            "## Size Distribution",
            "",
            "**Fit results:**",
            "",
            "| Parameter | Value |",
            "|-----------|-------|",
            f"| Method | {method} |",
            f"| Particle shape | {shape} |",
            f"| Chi-squared (χ²) | {_fmt(chi2)} |",
            f"| Volume fraction | {_fmt(vf, '.4e')} |",
            f"| Rg | {_fmt(rg)} Å |",
            f"| Peak r | {_fmt(peak_r)} Å |",
            f"| Iterations | {n_iter} |",
        ]
        if residuals is not None:
            L += [
                f"| Residuals mean | {np.mean(residuals):.4f} |",
                f"| Residuals std dev | {np.std(residuals):.4f} |",
            ]
        L.append("")

        L += [
            "**Grid / model setup:**",
            "",
            "| Parameter | Value |",
            "|-----------|-------|",
            f"| r min | {_fmt(r_min)} Å |",
            f"| r max | {_fmt(r_max)} Å |",
            f"| Bins | {n_bins} |",
            f"| Log spacing | {log_spacing} |",
            f"| Contrast (Δρ)² | {_fmt(contrast)} ×10²⁰ cm⁻⁴ |",
        ]
        try:
            if aspect_ratio is not None and not np.isnan(float(aspect_ratio)) and float(aspect_ratio) != 1.0:
                L.append(f"| Aspect ratio | {_fmt(aspect_ratio, '.4f')} |")
        except (TypeError, ValueError):
            pass
        L += [
            f"| Background | {_fmt(background, '.4e')} cm⁻¹ |",
            f"| Error scale | {_fmt(error_scale)} |",
        ]
        try:
            if power_law_B is not None and not np.isnan(float(power_law_B)) and float(power_law_B) != 0.0:
                L.append(f"| Power law B | {_fmt(power_law_B, '.4e')} |")
                L.append(f"| Power law P | {_fmt(power_law_P, '.4f')} |")
        except (TypeError, ValueError):
            pass
        if q_fit_min is not None and q_fit_max is not None:
            L.append(f"| Q range (fit) | {_fmt(q_fit_min)} – {_fmt(q_fit_max)} Å⁻¹ |")
        L.append("")

        # Method-specific parameters
        if str(method).lower() == 'maxent':
            sky      = _sv('maxent_sky_background')
            stab     = _sv('maxent_stability')
            max_iter = _sv('maxent_max_iter', 'N/A')
            L += [
                "**MaxEnt parameters:**",
                "",
                "| Parameter | Value |",
                "|-----------|-------|",
                f"| Sky background | {_fmt(sky)} |",
                f"| Stability | {_fmt(stab)} |",
                f"| Max iterations | {max_iter} |",
                "",
            ]
        elif str(method).lower() == 'regularization':
            evalue    = _sv('regularization_evalue')
            min_ratio = _sv('regularization_min_ratio')
            L += [
                "**Regularization parameters:**",
                "",
                "| Parameter | Value |",
                "|-----------|-------|",
                f"| Eigenvalue weight | {_fmt(evalue)} |",
                f"| Min ratio | {_fmt(min_ratio)} |",
                "",
            ]
        elif str(method).lower() == 'tnnls':
            approach = _sv('tnnls_approach_param')
            max_iter = _sv('tnnls_max_iter', 'N/A')
            L += [
                "**TNNLS parameters:**",
                "",
                "| Parameter | Value |",
                "|-----------|-------|",
                f"| Approach parameter | {_fmt(approach)} |",
                f"| Max iterations | {max_iter} |",
                "",
            ]
        elif str(method).lower() in ('montecarlo', 'mcsas'):
            n_rep    = _sv('montecarlo_n_repetitions') or _sv('mcsas_n_repetitions', 'N/A')
            conv     = _sv('montecarlo_convergence') or _sv('mcsas_convergence')
            max_iter = _sv('montecarlo_max_iter') or _sv('mcsas_max_iter', 'N/A')
            L += [
                "**Monte Carlo parameters:**",
                "",
                "| Parameter | Value |",
                "|-----------|-------|",
                f"| Repetitions | {n_rep} |",
                f"| Convergence | {_fmt(conv)} |",
                f"| Max iterations | {max_iter} |",
                "",
            ]

    # ── Simple Fits results ───────────────────────────────────────────────────
    if simple_fit_results is not None:
        sf_model    = simple_fit_results.get('model', 'unknown')
        sf_chi2     = simple_fit_results.get('chi_squared')
        sf_rchi2    = simple_fit_results.get('reduced_chi_squared')
        sf_dof      = simple_fit_results.get('dof')
        sf_q_min    = simple_fit_results.get('q_min')
        sf_q_max    = simple_fit_results.get('q_max')
        sf_complex  = simple_fit_results.get('use_complex_bg', False)
        sf_params   = simple_fit_results.get('params', {})
        sf_std      = simple_fit_results.get('params_std', {})
        sf_derived  = simple_fit_results.get('derived', {})

        def _sf_fmt(v, spec='.4g', suffix=''):
            if v is None:
                return 'N/A'
            try:
                if np.isnan(float(v)):
                    return 'N/A'
            except (TypeError, ValueError):
                pass
            try:
                return f"{v:{spec}}{suffix}"
            except (TypeError, ValueError):
                return str(v)

        L += [
            "## Simple Fits",
            "",
            "| Parameter | Value |",
            "|-----------|-------|",
            f"| Model | {sf_model} |",
            f"| Chi-squared (χ²) | {_sf_fmt(sf_chi2)} |",
            f"| Reduced chi² | {_sf_fmt(sf_rchi2)} |",
            f"| DOF | {sf_dof if sf_dof is not None else 'N/A'} |",
            f"| Q range (fit) | {_sf_fmt(sf_q_min)} – {_sf_fmt(sf_q_max)} Å⁻¹ |",
            f"| Complex background | {sf_complex} |",
            "",
        ]

        if sf_params:
            has_std = bool(sf_std)
            if has_std:
                L += ["**Parameters:**", "",
                      "| Parameter | Value | Uncertainty (1σ) |",
                      "|-----------|-------|------------------|"]
            else:
                L += ["**Parameters:**", "",
                      "| Parameter | Value |",
                      "|-----------|-------|"]
            for name, val in sf_params.items():
                val_str = _sf_fmt(val, '.6g')
                if has_std:
                    err = sf_std.get(name)
                    err_str = f"± {_sf_fmt(err, '.3g')}" if err is not None else "—"
                    L.append(f"| {name} | {val_str} | {err_str} |")
                else:
                    L.append(f"| {name} | {val_str} |")
            L.append("")

        if sf_derived:
            L += ["**Derived quantities:**", "",
                  "| Quantity | Value |",
                  "|----------|-------|"]
            for name, val in sf_derived.items():
                L.append(f"| {name} | {_sf_fmt(val, '.6g')} |")
            L.append("")

    L += ["---", "*Generated by pyIrena*", ""]
    return "\n".join(L)


class DataSelectorConfigDialog(QDialog):
    """
    Extensible configuration dialog for the Data Selector.

    Settings are defined as a list of field specifications (FIELD_SPECS).
    Adding a new configurable parameter only requires adding one entry to that list.

    Supported field types
    ---------------------
    'float'  — QLineEdit with QDoubleValidator
    'int'    — QLineEdit with integer validation
    'str'    — plain QLineEdit
    'bool'   — QCheckBox
    'color'  — QPushButton that opens QColorDialog (stores hex color string)
    """

    # ---------------------------------------------------------------------------
    # Field specifications — add new settings here
    # ---------------------------------------------------------------------------
    FIELD_SPECS = [
        {
            'group':    'Text File Options',
            'key':      'error_fraction',
            'label':    'Generated uncertainty fraction',
            'tooltip':  (
                'When a text file has only two columns (Q and I) and no uncertainty\n'
                'column, the uncertainty is generated as:  σ = I × this_value\n'
                'Default: 0.05  (5 % of intensity)'
            ),
            'type':     'float',
            'default':  0.05,
            'min':      0.0,
            'max':      100.0,
            'decimals': 4,
        },
        {
            'group':   'Graph Options',
            'key':     'max_legend_items',
            'label':   'Maximum items in legend',
            'tooltip': (
                'When more datasets than this limit are plotted, the legend shows\n'
                'only the first, last, and evenly spaced items in between.\n'
                'Default: 12'
            ),
            'type':    'int',
            'default': 12,
            'min':     2,
            'max':     200,
            'decimals': 0,
        },
        # -----------------------------------------------------------------------
        # Future settings — just append a dict here, no other code changes needed
        # -----------------------------------------------------------------------
        # {
        #     'group':   'Text File Options',
        #     'key':     'q_units_scale',
        #     'label':   'Q unit scale factor',
        #     'tooltip': 'Multiply Q by this factor on load (1.0 = no change)',
        #     'type':    'float',
        #     'default': 1.0,
        #     'min':     0.0,
        #     'max':     1000.0,
        #     'decimals': 6,
        # },
        # {
        #     'group':   'Display',
        #     'key':     'plot_color',
        #     'label':   'Default plot color',
        #     'tooltip': 'Color used for single-file plots',
        #     'type':    'color',
        #     'default': '#3498db',
        # },
    ]

    def __init__(self, current_values: dict, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Data Selector — Configure")
        self.setMinimumWidth(420)
        self._widgets = {}   # key -> (field_type, widget)
        self._init_ui(current_values)

    def _init_ui(self, current_values: dict):
        outer = QVBoxLayout()
        outer.setSpacing(12)

        # Group fields by their 'group' key
        groups = {}
        for spec in self.FIELD_SPECS:
            g = spec.get('group', 'General')
            groups.setdefault(g, []).append(spec)

        for group_name, specs in groups.items():
            box = QGroupBox(group_name)
            form = QFormLayout()
            form.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapLongRows)

            for spec in specs:
                key       = spec['key']
                label_txt = spec['label']
                ftype     = spec['type']
                value     = current_values.get(key, spec.get('default', ''))
                tooltip   = spec.get('tooltip', '')

                if ftype in ('float', 'int'):
                    widget = QLineEdit(str(value))
                    validator = QDoubleValidator(
                        float(spec.get('min', -1e300)),
                        float(spec.get('max',  1e300)),
                        int(spec.get('decimals', 6)),
                    )
                    widget.setValidator(validator)
                    widget.setMaximumWidth(120)

                elif ftype == 'bool':
                    widget = QCheckBox()
                    widget.setChecked(bool(value))

                elif ftype == 'color':
                    widget = QPushButton()
                    widget._color = str(value)
                    widget.setStyleSheet(f"background-color: {value};")
                    widget.setFixedSize(60, 24)
                    widget.clicked.connect(
                        lambda checked, btn=widget: self._pick_color(btn)
                    )

                else:   # 'str'
                    widget = QLineEdit(str(value))

                if tooltip:
                    widget.setToolTip(tooltip)

                lbl = QLabel(label_txt)
                if tooltip:
                    lbl.setToolTip(tooltip)

                form.addRow(lbl, widget)
                self._widgets[key] = (ftype, widget)

            box.setLayout(form)
            outer.addWidget(box)

        # OK / Cancel buttons
        btn_box = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        btn_box.accepted.connect(self.accept)
        btn_box.rejected.connect(self.reject)
        outer.addWidget(btn_box)

        self.setLayout(outer)

    def _pick_color(self, btn):
        color = QColorDialog.getColor(parent=self)
        if color.isValid():
            btn._color = color.name()
            btn.setStyleSheet(f"background-color: {color.name()};")

    def get_values(self) -> dict:
        """Return validated values from all widgets keyed by field key."""
        result = {}
        for spec in self.FIELD_SPECS:
            key   = spec['key']
            ftype = spec['type']
            _, widget = self._widgets[key]

            if ftype in ('float', 'int'):
                try:
                    result[key] = float(widget.text()) if ftype == 'float' else int(widget.text())
                except ValueError:
                    result[key] = spec.get('default', 0)
            elif ftype == 'bool':
                result[key] = widget.isChecked()
            elif ftype == 'color':
                result[key] = widget._color
            else:
                result[key] = widget.text()
        return result


class GraphWindow(QWidget):
    """
    Separate window for displaying raw SAS data.
    Uses pyqtgraph for fast, interactive rendering even with large datasets.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("pyIrena - Data Viewer")
        self.setGeometry(100, 100, 850, 620)

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
        self.plot.addLegend(offset=(-10, 10), labelTextSize='18pt')
        _style_plot(self.plot)
        _add_jpeg_export(self, self.plot)

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
        self.plot.clear()
        if self.plot.legend is not None:
            self.plot.legend.clear()

        legend_idx = _legend_indices(len(file_paths), max_legend_items)
        colors = _gen_colors(len(file_paths))
        for idx, file_path in enumerate(file_paths):
            color = colors[idx]
            try:
                path, filename = os.path.split(file_path)
                _, ext = os.path.splitext(filename)

                if ext.lower() in ('.txt', '.dat'):
                    data = readTextFile(path, filename, error_fraction=error_fraction)
                else:
                    data = readGenericNXcanSAS(path, filename)

                if data is None:
                    continue

                q   = np.asarray(data['Q'],         dtype=float)
                I   = np.asarray(data['Intensity'],  dtype=float)
                err = data.get('Error')

                label = os.path.basename(file_path)
                name  = label if idx in legend_idx else None

                self.plot.plot(
                    q, I,
                    pen=None, symbol='o', symbolSize=4,
                    symbolPen=pg.mkPen(color, width=1),
                    symbolBrush=pg.mkBrush(color),
                    name=name,
                )
                if name is not None and self.plot.legend is not None and self.plot.legend.items:
                    self.plot.legend.items[-1][1].setAttr('color', color.name())

                if err is not None:
                    err = np.asarray(err, dtype=float)
                    xb, yb = _iq_error_bars(q, I, err)
                    if len(xb):
                        self.plot.plot(xb, yb,
                                       pen=pg.mkPen(color, width=1),
                                       connect='finite')

            except Exception as e:
                print(f"Error loading {file_path}: {e}")

        self.show()


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
        self.ax_main.addLegend(offset=(-10, 10), labelTextSize='18pt')
        _style_plot(self.ax_main)

        # ── Bottom panel: residuals ────────────────────────────────────────
        self.ax_resid = self.gl.addPlot(
            row=1, col=0,
            axisItems={'bottom': _LogDecadeAxis(orientation='bottom')},
        )
        self.ax_resid.setLogMode(True, False)
        self.ax_resid.setLabel('bottom', 'Q (Å⁻¹)')
        self.ax_resid.setLabel('left', 'Residuals (norm.)')
        self.ax_resid.showGrid(x=True, y=True, alpha=0.3)
        self.ax_resid.addLegend(offset=(-10, 10), labelTextSize='18pt')
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

            # ── residuals ──────────────────────────────────────────────────
            resid_name = label if in_legend else None
            self.ax_resid.plot(
                Q, residuals,
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
        self.ax_main.addLegend(offset=(-10, 10), labelTextSize='18pt')
        _style_plot(self.ax_main)

        # ── Middle: residuals ──────────────────────────────────────────────
        self.ax_resid = self.gl.addPlot(
            row=1, col=0,
            axisItems={'bottom': _LogDecadeAxis(orientation='bottom')},
        )
        self.ax_resid.setLogMode(True, False)
        self.ax_resid.setLabel('bottom', 'Q (Å⁻¹)')
        self.ax_resid.setLabel('left', 'Residuals (norm.)')
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
        self.ax_dist.addLegend(offset=(-10, 10), labelTextSize='18pt')
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

            # ── residuals ──────────────────────────────────────────────────
            if residuals is not None:
                self.ax_resid.plot(
                    Q, residuals,
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
        self.ax_main.addLegend(offset=(-10, 10), labelTextSize='18pt')
        _style_plot(self.ax_main)

        # ── Bottom: residuals ──────────────────────────────────────────────
        self.ax_resid = self.gl.addPlot(
            row=1, col=0,
            axisItems={'bottom': _LogDecadeAxis(orientation='bottom')},
        )
        self.ax_resid.setLogMode(True, False)
        self.ax_resid.setLabel('bottom', 'Q (Å⁻¹)')
        self.ax_resid.setLabel('left', 'Residuals (norm.)')
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
            chi2_str = f'{chi2:.3f}' if (chi2 == chi2) else 'N/A'
            fit_name = f'{label}  {model}  χ²={chi2_str}' if in_legend else None
            self.ax_main.plot(
                Q, I_model,
                pen=pg.mkPen(fit_color, width=3.0),
                name=fit_name,
            )
            if fit_name is not None and self.ax_main.legend is not None and self.ax_main.legend.items:
                self.ax_main.legend.items[-1][1].setAttr('color', fit_color.name())

            # ── residuals ──────────────────────────────────────────────────
            if residuals is not None:
                self.ax_resid.plot(
                    Q, residuals,
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


class BatchWorker(QThread):
    """
    Background thread for running batch fitting (Unified Fit or Size Distribution)
    on a list of files with a shared pyirena_config.json.
    """
    progress = Signal(str)              # emitted before each file: "Working: N/M — name"
    finished = Signal(int, int, list)   # n_ok, n_fail, per-file messages

    def __init__(self, tool: str, file_paths: list, config_file: str, parent=None):
        super().__init__(parent)
        self.tool = tool                # 'unified', 'sizes', or 'simple_fits'
        self.file_paths = file_paths
        self.config_file = config_file

    def run(self):
        n_ok, n_fail = 0, 0
        messages = []
        if self.tool == 'unified':
            fit_fn = fit_unified
        elif self.tool == 'simple_fits':
            fit_fn = fit_simple_from_config
        else:
            fit_fn = fit_sizes
        total = len(self.file_paths)

        for i, fp in enumerate(self.file_paths):
            fname = os.path.basename(fp)
            self.progress.emit(f"Working: {i + 1}/{total} — {fname}")
            try:
                result = fit_fn(fp, self.config_file, save_to_nexus=True)
                if result.get('success', False):
                    n_ok += 1
                    messages.append(f"✓ {fname}")
                else:
                    n_fail += 1
                    messages.append(f"✗ {fname}: {result.get('message', 'fit failed')}")
            except Exception as exc:
                n_fail += 1
                messages.append(f"✗ {fname}: {exc}")

        self.finished.emit(n_ok, n_fail, messages)


class DataSelectorPanel(QWidget):
    """
    Main data selector panel for pyIrena.

    Provides file browsing, filtering, and data visualization capabilities.
    """

    # File extensions for different data types
    HDF5_EXTENSIONS = ['.hdf', '.h5', '.hdf5']
    TEXT_EXTENSIONS = ['.txt', '.dat']

    def __init__(self):
        super().__init__()
        self.current_folder = None
        self.last_folder = None  # Remember last selected folder
        self.graph_window = None
        self.unified_fit_results_window = None  # Graph of stored fit results
        self.size_dist_results_window = None   # Graph of stored size dist results
        self.tabulate_results_window = None    # Tabulated results window
        self.unified_fit_window = None  # Unified fit panel
        self.sizes_fit_window = None   # Size distribution panel
        self.simple_fits_window = None         # Simple Fits panel
        self.simple_fits_results_window = None # Graph of stored simple fit results
        self._batch_worker = None      # Batch fitting thread

        # Initialize state manager
        self.state_manager = StateManager()

        # Load last used folder from state
        self.load_last_folder()

        self.init_ui()

        # Restore saved sort selection and re-sort the already-populated list
        saved_sort = int(self.state_manager.get('data_selector', 'sort_index', 0) or 0)
        self.sort_combo.blockSignals(True)
        self.sort_combo.setCurrentIndex(saved_sort)
        self.sort_combo.blockSignals(False)
        self.sort_file_list()   # apply restored sort order to the initial file list

    def init_ui(self):
        """Initialize the user interface."""
        self.setWindowTitle("pyIrena - Data Selector")

        # Main layout
        main_layout = QVBoxLayout()
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)

        # Menu bar
        menu_bar = self.create_menu_bar()
        main_layout.addWidget(menu_bar)

        # Content layout (with margins)
        content_layout = QVBoxLayout()
        content_layout.setContentsMargins(20, 20, 20, 20)
        content_layout.setSpacing(15)

        # Title
        title_label = QLabel("pyIrena")
        title_label.setStyleSheet("""
            QLabel {
                font-size: 24px;
                font-weight: bold;
                color: #2c3e50;
                padding: 10px;
            }
        """)
        title_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        content_layout.addWidget(title_label)

        # Folder selection section
        folder_layout = QHBoxLayout()
        self.folder_button = QPushButton("Select Data Folder")
        self.folder_button.setMinimumHeight(40)
        self.folder_button.clicked.connect(self.select_folder)
        folder_layout.addWidget(self.folder_button)

        self.refresh_button = QPushButton("Refresh")
        self.refresh_button.setMinimumHeight(40)
        self.refresh_button.setMaximumWidth(100)
        self.refresh_button.clicked.connect(self.refresh_file_list)
        self.refresh_button.setEnabled(False)
        folder_layout.addWidget(self.refresh_button)

        self.folder_label = QLabel("No folder selected")
        self.folder_label.setStyleSheet("color: #7f8c8d; font-style: italic;")
        folder_layout.addWidget(self.folder_label)
        folder_layout.addStretch()

        content_layout.addLayout(folder_layout)

        # File type + Sort selection (combined row)
        type_layout = QHBoxLayout()
        type_layout.addWidget(QLabel("File Type:"))

        self.file_type_combo = QComboBox()
        self.file_type_combo.addItem("HDF5 Files (.hdf, .h5, .hdf5)", "hdf5")
        self.file_type_combo.addItem("Text Files (.txt, .dat)", "text")
        self.file_type_combo.addItem("All Supported Files", "all")
        self.file_type_combo.currentIndexChanged.connect(self.refresh_file_list)
        type_layout.addWidget(self.file_type_combo)

        type_layout.addSpacing(12)
        type_layout.addWidget(QLabel("Sort:"))
        self.sort_combo = QComboBox()
        self.sort_combo.addItems([
            "Filename  A→Z",
            "Filename  Z→A",
            "Temperature  ↑",
            "Temperature  ↓",
            "Time  ↑",
            "Time  ↓",
            "Order number  ↑",
            "Order number  ↓",
            "Pressure  ↑",
            "Pressure  ↓",
        ])
        self.sort_combo.setToolTip(
            "Sort order for the file list.\n"
            "Patterns recognised in filenames:\n"
            "  Temperature : _25C\n"
            "  Time        : _50min\n"
            "  Order number: _354  (last underscore-number before extension)\n"
            "  Pressure    : _35PSI"
        )
        self.sort_combo.currentIndexChanged.connect(self._on_sort_changed)
        type_layout.addWidget(self.sort_combo)
        type_layout.addStretch()

        content_layout.addLayout(type_layout)

        # Content area (listbox + graph button)
        file_area_layout = QHBoxLayout()

        # Left side: file list section
        left_layout = QVBoxLayout()

        # File filter input
        filter_layout = QHBoxLayout()
        filter_layout.addWidget(QLabel("Filter:"))
        self.filter_input = QLineEdit()
        self.filter_input.setPlaceholderText("Enter text to filter files...")
        self.filter_input.textChanged.connect(self.filter_files)
        filter_layout.addWidget(self.filter_input)
        left_layout.addLayout(filter_layout)

        # File list
        self.file_list = QListWidget()
        self.file_list.setMinimumWidth(400)  # At least 30 characters wide
        self.file_list.setMinimumHeight(400)  # Show at least 15 items
        self.file_list.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.file_list.itemDoubleClicked.connect(self.plot_selected_files)
        self.file_list.itemSelectionChanged.connect(self.update_plot_button_state)
        left_layout.addWidget(self.file_list)

        # Configure button — small, sits below the file list
        configure_row = QHBoxLayout()
        self.configure_button = QPushButton("Configure...")
        self.configure_button.setMaximumWidth(110)
        self.configure_button.setMinimumHeight(24)
        self.configure_button.setToolTip("Configure data loading options")
        self.configure_button.clicked.connect(self.open_configure_dialog)
        configure_row.addWidget(self.configure_button)
        configure_row.addStretch()
        left_layout.addLayout(configure_row)

        file_area_layout.addLayout(left_layout, stretch=2)

        # Right side: action buttons — starts at same level as Filter row
        right_layout = QVBoxLayout()
        right_layout.setSpacing(6)

        # ── Batch script status display ────────────────────────────────────
        self.batch_status_label = QLabel("")
        self.batch_status_label.setWordWrap(True)
        self.batch_status_label.setMinimumHeight(28)
        self.batch_status_label.setMaximumHeight(56)
        self.batch_status_label.setVisible(False)
        self.batch_status_label.setStyleSheet(
            "QLabel { padding: 5px 8px; border-radius: 4px; font-size: 12px; }"
        )
        right_layout.addWidget(self.batch_status_label)

        # ── Graph content checkboxes ───────────────────────────────────────
        graph_content_label = QLabel("Show in graph/reports:")
        graph_content_label.setStyleSheet("font-weight: bold; color: #2c3e50;")
        right_layout.addWidget(graph_content_label)

        cb_row = QHBoxLayout()
        cb_row.setSpacing(10)
        self.data_checkbox = QCheckBox("Data")
        self.data_checkbox.setChecked(True)
        self.data_checkbox.setToolTip("Plot experimental data for selected files")
        cb_row.addWidget(self.data_checkbox)

        self.unified_fit_result_checkbox = QCheckBox("Unified Fit")
        self.unified_fit_result_checkbox.setChecked(False)
        self.unified_fit_result_checkbox.setToolTip(
            "Plot stored Unified Fit results (data + model + residuals).\n"
            "Only HDF5 files are checked; files without fit results are skipped."
        )
        cb_row.addWidget(self.unified_fit_result_checkbox)

        self.size_dist_checkbox = QCheckBox("Size Dist.")
        self.size_dist_checkbox.setChecked(False)
        self.size_dist_checkbox.setToolTip(
            "Open Size Distribution panel with selected data (Create Graph),\n"
            "or include stored sizes results in report (Create Report).\n"
            "Only HDF5 files with stored sizes results are used for reports."
        )
        cb_row.addWidget(self.size_dist_checkbox)

        self.simple_fits_checkbox = QCheckBox("Simple Fits")
        self.simple_fits_checkbox.setChecked(False)
        self.simple_fits_checkbox.setToolTip(
            "Plot stored Simple Fits results (data + model + residuals).\n"
            "Only HDF5 files with stored simple fit results are used."
        )
        cb_row.addWidget(self.simple_fits_checkbox)
        cb_row.addStretch()
        right_layout.addLayout(cb_row)

        right_layout.addSpacing(4)

        # ── Create Graph / Create Report side by side ──────────────────────
        _btn_style_blue = """
            QPushButton {
                background-color: #3498db; color: white;
                font-size: 12px; font-weight: bold;
                border-radius: 4px; padding: 5px 8px;
            }
            QPushButton:hover { background-color: #2980b9; }
            QPushButton:disabled { background-color: #bdc3c7; }
        """
        _btn_style_purple = """
            QPushButton {
                background-color: #8e44ad; color: white;
                font-size: 12px; font-weight: bold;
                border-radius: 4px; padding: 5px 8px;
            }
            QPushButton:hover { background-color: #7d3c98; }
            QPushButton:disabled { background-color: #bdc3c7; }
        """

        self.plot_button = QPushButton("Create Graph")
        self.plot_button.setMinimumHeight(34)
        self.plot_button.setStyleSheet(_btn_style_blue)
        self.plot_button.clicked.connect(self.plot_selected_files)
        self.plot_button.setEnabled(False)

        self.report_button = QPushButton("Create Report")
        self.report_button.setMinimumHeight(34)
        self.report_button.setStyleSheet(_btn_style_purple)
        self.report_button.setToolTip(
            "Generate a Markdown report (.md) summarising the Unified Fit\n"
            "results for each selected HDF5 file.  Files without stored\n"
            "fit results are skipped."
        )
        self.report_button.clicked.connect(self.create_report)
        self.report_button.setEnabled(False)

        _btn_style_teal = """
            QPushButton {
                background-color: #16a085; color: white;
                font-size: 12px; font-weight: bold;
                border-radius: 4px; padding: 5px 8px;
            }
            QPushButton:hover { background-color: #138d75; }
            QPushButton:disabled { background-color: #bdc3c7; }
        """
        self.tabulate_button = QPushButton("Tabulate Results")
        self.tabulate_button.setMinimumHeight(34)
        self.tabulate_button.setStyleSheet(_btn_style_teal)
        self.tabulate_button.setToolTip(
            "Build a table of fit results for selected files and display it.\n"
            "Results included depend on the checked checkboxes.\n"
            "You can save the table as a CSV file (Excel, Igor Pro, Origin, etc.)."
        )
        self.tabulate_button.clicked.connect(self.tabulate_results)
        self.tabulate_button.setEnabled(False)

        right_layout.addSpacing(10)

        # ── Unified Fit: GUI button + Script batch button side by side ─────
        _uf_gui_style = """
            QPushButton {
                background-color: #27ae60; color: white;
                font-size: 12px; font-weight: bold;
                border-radius: 4px; padding: 6px 8px;
            }
            QPushButton:hover { background-color: #229954; }
            QPushButton:disabled { background-color: #bdc3c7; }
        """
        _uf_script_style = """
            QPushButton {
                background-color: #1e8449; color: white;
                font-size: 12px; font-weight: bold;
                border-radius: 4px; padding: 6px 8px;
            }
            QPushButton:hover { background-color: #196f3d; }
            QPushButton:disabled { background-color: #bdc3c7; }
        """
        self.unified_fit_button = QPushButton("Unified Fit (GUI)")
        self.unified_fit_button.setMinimumHeight(38)
        self.unified_fit_button.setStyleSheet(_uf_gui_style)
        self.unified_fit_button.setToolTip(
            "Open Unified Fit panel for the first selected file."
        )
        self.unified_fit_button.clicked.connect(self.launch_unified_fit)
        self.unified_fit_button.setEnabled(False)

        self.unified_script_button = QPushButton("Unified Fit (script)")
        self.unified_script_button.setMinimumHeight(38)
        self.unified_script_button.setStyleSheet(_uf_script_style)
        self.unified_script_button.setToolTip(
            "Batch-fit all selected files with Unified Fit using pyirena_config.json.\n"
            "Results are saved into each file's NXcanSAS record."
        )
        self.unified_script_button.clicked.connect(self.run_unified_script)
        self.unified_script_button.setEnabled(False)

        # ── Size Distribution: GUI button + Script batch button ────────────
        _sz_gui_style = """
            QPushButton {
                background-color: #2980b9; color: white;
                font-size: 12px; font-weight: bold;
                border-radius: 4px; padding: 6px 8px;
            }
            QPushButton:hover { background-color: #2471a3; }
            QPushButton:disabled { background-color: #bdc3c7; }
        """
        _sz_script_style = """
            QPushButton {
                background-color: #1f618d; color: white;
                font-size: 12px; font-weight: bold;
                border-radius: 4px; padding: 6px 8px;
            }
            QPushButton:hover { background-color: #1a5276; }
            QPushButton:disabled { background-color: #bdc3c7; }
        """
        self.sizes_fit_button = QPushButton("Size Distribution (GUI)")
        self.sizes_fit_button.setMinimumHeight(38)
        self.sizes_fit_button.setStyleSheet(_sz_gui_style)
        self.sizes_fit_button.setToolTip(
            "Open Size Distribution panel for the first selected file."
        )
        self.sizes_fit_button.clicked.connect(self.launch_sizes_fit)
        self.sizes_fit_button.setEnabled(False)

        self.sizes_script_button = QPushButton("Size Distribution (script)")
        self.sizes_script_button.setMinimumHeight(38)
        self.sizes_script_button.setStyleSheet(_sz_script_style)
        self.sizes_script_button.setToolTip(
            "Batch-fit all selected files with Size Distribution using pyirena_config.json.\n"
            "Results are saved into each file's NXcanSAS record."
        )
        self.sizes_script_button.clicked.connect(self.run_sizes_script)
        self.sizes_script_button.setEnabled(False)

        # ── Single grid so all buttons share equal column widths ───────────
        # Row 0: Create Graph | Create Report
        # Row 1: Tabulate Results (spans both columns)
        # Row 2: spacer
        # Row 3: Unified Fit (GUI) | Unified Fit (script)
        # Row 4: Size Distribution (GUI) | Size Distribution (script)
        btn_grid = QGridLayout()
        btn_grid.setHorizontalSpacing(4)
        btn_grid.setVerticalSpacing(4)
        btn_grid.setColumnStretch(0, 1)
        btn_grid.setColumnStretch(1, 1)
        btn_grid.addWidget(self.plot_button,            0, 0)
        btn_grid.addWidget(self.report_button,          0, 1)
        btn_grid.addWidget(self.tabulate_button,        1, 0, 1, 2)  # full-width
        btn_grid.setRowMinimumHeight(2, 10)                          # visual separator
        btn_grid.addWidget(self.unified_fit_button,     3, 0)
        btn_grid.addWidget(self.unified_script_button,  3, 1)
        btn_grid.addWidget(self.sizes_fit_button,       4, 0)
        btn_grid.addWidget(self.sizes_script_button,    4, 1)

        # ── Simple Fits: GUI button + Script batch button ─────────────────
        _sf_gui_style = """
            QPushButton {
                background-color: #27ae60; color: white;
                font-size: 12px; font-weight: bold;
                border-radius: 4px; padding: 4px;
                border: none;
            }
            QPushButton:hover { background-color: #1e8449; }
            QPushButton:disabled { background-color: #bdc3c7; }
        """
        _sf_script_style = """
            QPushButton {
                background-color: #1e8449; color: white;
                font-size: 12px; font-weight: bold;
                border-radius: 4px; padding: 4px;
                border: none;
            }
            QPushButton:hover { background-color: #196f3d; }
            QPushButton:disabled { background-color: #bdc3c7; }
        """
        self.simple_fits_button = QPushButton("Simple Fits (GUI)")
        self.simple_fits_button.setMinimumHeight(38)
        self.simple_fits_button.setStyleSheet(_sf_gui_style)
        self.simple_fits_button.setToolTip(
            "Open Simple Fits panel for the first selected file."
        )
        self.simple_fits_button.clicked.connect(self.launch_simple_fits)
        self.simple_fits_button.setEnabled(False)

        self.simple_fits_script_button = QPushButton("Simple Fits (script)")
        self.simple_fits_script_button.setMinimumHeight(38)
        self.simple_fits_script_button.setStyleSheet(_sf_script_style)
        self.simple_fits_script_button.setToolTip(
            "Batch-fit all selected files with Simple Fits using pyirena_config.json.\n"
            "Results are saved into each file's NXcanSAS record."
        )
        self.simple_fits_script_button.clicked.connect(self.run_simple_fits_script)
        self.simple_fits_script_button.setEnabled(False)

        btn_grid.addWidget(self.simple_fits_button,        5, 0)
        btn_grid.addWidget(self.simple_fits_script_button, 5, 1)

        right_layout.addLayout(btn_grid)

        right_layout.addStretch()
        file_area_layout.addLayout(right_layout, stretch=1)

        content_layout.addLayout(file_area_layout)

        # Status bar
        self.status_label = QLabel("Ready - Select a folder to begin")
        self.status_label.setStyleSheet("""
            QLabel {
                color: #7f8c8d;
                padding: 5px;
                border-top: 1px solid #bdc3c7;
            }
        """)
        content_layout.addWidget(self.status_label)

        # Add content layout to main layout
        main_layout.addLayout(content_layout)

        self.setLayout(main_layout)

        # Set minimum window size (at least twice the listbox width)
        self.setMinimumSize(900, 600)

        # If we have a restored folder, display it and list files
        if self.current_folder:
            self.folder_label.setText(self.current_folder)
            self.folder_label.setStyleSheet("color: #2c3e50;")
            self.refresh_button.setEnabled(True)
            self.refresh_file_list()
            self.status_label.setText(f"Restored folder: {os.path.basename(self.current_folder)}")

    def create_menu_bar(self) -> QMenuBar:
        """Create the menu bar with Models menu."""
        menu_bar = QMenuBar()

        # Models menu
        models_menu = QMenu("&Models", self)

        # Unified Fit action
        unified_fit_action = QAction("&Unified Fit", self)
        unified_fit_action.setStatusTip("Open Unified Fit model panel")
        unified_fit_action.triggered.connect(self.launch_unified_fit)
        models_menu.addAction(unified_fit_action)

        # Size Distribution action
        sizes_fit_action = QAction("&Size Distribution", self)
        sizes_fit_action.setStatusTip("Open Size Distribution fitting panel")
        sizes_fit_action.triggered.connect(self.launch_sizes_fit)
        models_menu.addAction(sizes_fit_action)

        # Simple Fits action
        simple_fits_action = QAction("Simple &Fits", self)
        simple_fits_action.setStatusTip("Open Simple Fits panel")
        simple_fits_action.triggered.connect(self.launch_simple_fits)
        models_menu.addAction(simple_fits_action)

        # Add separator and future models placeholder
        models_menu.addSeparator()
        placeholder_action = QAction("More models coming soon...", self)
        placeholder_action.setEnabled(False)
        models_menu.addAction(placeholder_action)

        menu_bar.addMenu(models_menu)

        # Help menu
        help_menu = QMenu("&Help", self)

        about_action = QAction("&About pyIrena", self)
        about_action.setStatusTip("About pyIrena")
        about_action.triggered.connect(self.show_about)
        help_menu.addAction(about_action)

        menu_bar.addMenu(help_menu)

        return menu_bar

    def select_folder(self):
        """Open folder selection dialog."""
        # Use last folder if available, otherwise home directory
        start_dir = self.last_folder if self.last_folder else QDir.homePath()

        folder = QFileDialog.getExistingDirectory(
            self,
            "Select Data Folder",
            start_dir,
            QFileDialog.Option.ShowDirsOnly
        )

        if folder:
            self.current_folder = folder
            self.last_folder = folder  # Remember for next time
            self.save_last_folder(folder)  # Save to state
            self.folder_label.setText(folder)
            self.folder_label.setStyleSheet("color: #2c3e50;")
            self.refresh_button.setEnabled(True)
            self.refresh_file_list()
            self.status_label.setText(f"Folder: {os.path.basename(folder)}")

    def get_file_extensions(self) -> List[str]:
        """Get current file extensions based on selection."""
        file_type = self.file_type_combo.currentData()

        if file_type == "hdf5":
            return self.HDF5_EXTENSIONS
        elif file_type == "text":
            return self.TEXT_EXTENSIONS
        else:  # all
            return self.HDF5_EXTENSIONS + self.TEXT_EXTENSIONS

    def refresh_file_list(self):
        """Refresh the file list based on current folder and file type."""
        self.file_list.clear()

        if not self.current_folder:
            return

        extensions = self.get_file_extensions()

        try:
            files = []
            for file in os.listdir(self.current_folder):
                file_path = os.path.join(self.current_folder, file)
                if os.path.isfile(file_path):
                    _, ext = os.path.splitext(file)
                    if ext.lower() in extensions:
                        files.append(file)

            self.file_list.addItems(files)
            self.status_label.setText(f"Found {len(files)} files")
            self.sort_file_list()   # apply current sort-combo order
            self.update_plot_button_state()

        except Exception as e:
            self.status_label.setText(f"Error reading folder: {e}")

    def filter_files(self, filter_text: str):
        """Filter the file list based on search text."""
        if not filter_text:
            # Show all items
            for i in range(self.file_list.count()):
                self.file_list.item(i).setHidden(False)
            return

        # Use grep-like filtering (case-insensitive substring match)
        filter_text = filter_text.lower()
        visible_count = 0

        for i in range(self.file_list.count()):
            item = self.file_list.item(i)
            matches = filter_text in item.text().lower()
            item.setHidden(not matches)
            if matches:
                visible_count += 1

        self.status_label.setText(f"Showing {visible_count} of {self.file_list.count()} files")

    def _on_sort_changed(self, index: int):
        """Save sort selection to state, then apply it to the file list."""
        self.state_manager.set('data_selector', 'sort_index', index)
        self.state_manager.save()
        self.sort_file_list()

    def sort_file_list(self):
        """
        Re-order all items in the file list according to the sort combo selection.
        The current text filter is re-applied afterwards so hidden items stay hidden.
        """
        n = self.file_list.count()
        if n == 0:
            return

        items = [self.file_list.item(i).text() for i in range(n)]

        idx     = self.sort_combo.currentIndex()
        reverse = bool(idx % 2)                 # odd indices → descending
        key_fn  = _SORT_KEYS[min(idx, len(_SORT_KEYS) - 1)]

        items.sort(key=key_fn, reverse=reverse)

        self.file_list.clear()
        self.file_list.addItems(items)

        # Re-apply whatever text filter is active
        self.filter_files(self.filter_input.text())

    def update_plot_button_state(self):
        """Enable or disable buttons based on file selection."""
        has_selection = len(self.file_list.selectedItems()) > 0
        self.plot_button.setEnabled(has_selection)
        self.report_button.setEnabled(has_selection)
        self.tabulate_button.setEnabled(has_selection)
        self.unified_fit_button.setEnabled(has_selection)
        self.unified_script_button.setEnabled(has_selection)
        self.sizes_fit_button.setEnabled(has_selection)
        self.sizes_script_button.setEnabled(has_selection)
        self.simple_fits_button.setEnabled(has_selection)
        self.simple_fits_script_button.setEnabled(has_selection)

    def plot_selected_files(self):
        """Plot the selected files according to the Data / Unified Fit checkboxes."""
        selected_items = self.file_list.selectedItems()

        if not selected_items:
            QMessageBox.warning(
                self,
                "No Selection",
                "Please select one or more files to plot."
            )
            return

        file_paths = [
            os.path.join(self.current_folder, item.text())
            for item in selected_items
        ]

        error_fraction    = self.state_manager.get('data_selector', 'error_fraction', 0.05)
        max_legend_items  = int(self.state_manager.get('data_selector', 'max_legend_items', 12))
        plotted = []

        # ── Experimental data ──────────────────────────────────────────────
        if self.data_checkbox.isChecked():
            if self.graph_window is None:
                self.graph_window = GraphWindow()
            try:
                self.graph_window.plot_data(
                    file_paths,
                    error_fraction=error_fraction,
                    max_legend_items=max_legend_items,
                )
                plotted.append("data")
            except Exception as e:
                QMessageBox.critical(
                    self, "Plot Error", f"Error creating data plot:\n{str(e)}"
                )

        # ── Unified Fit results ────────────────────────────────────────────
        if self.unified_fit_result_checkbox.isChecked():
            if self.unified_fit_results_window is None:
                self.unified_fit_results_window = UnifiedFitResultsWindow()
            try:
                self.unified_fit_results_window.plot_results(
                    file_paths,
                    max_legend_items=max_legend_items,
                )
                plotted.append("Unified Fit results")
            except Exception as e:
                QMessageBox.critical(
                    self, "Plot Error", f"Error creating Unified Fit plot:\n{str(e)}"
                )

        # ── Size Distribution results ───────────────────────────────────────
        if self.size_dist_checkbox.isChecked():
            if self.size_dist_results_window is None:
                self.size_dist_results_window = SizeDistResultsWindow()
            try:
                self.size_dist_results_window.plot_results(
                    file_paths,
                    max_legend_items=max_legend_items,
                )
                plotted.append("Size Distribution results")
            except Exception as e:
                QMessageBox.critical(
                    self, "Plot Error", f"Error creating Size Distribution plot:\n{str(e)}"
                )

        # ── Simple Fits results ────────────────────────────────────────────
        if self.simple_fits_checkbox.isChecked():
            if self.simple_fits_results_window is None:
                self.simple_fits_results_window = SimpleFitResultsWindow()
            try:
                self.simple_fits_results_window.plot_results(
                    file_paths,
                    max_legend_items=max_legend_items,
                )
                plotted.append("Simple Fits results")
            except Exception as e:
                QMessageBox.critical(
                    self, "Plot Error", f"Error creating Simple Fits plot:\n{str(e)}"
                )

        if plotted:
            self.status_label.setText(
                f"Plotted {len(file_paths)} file(s): {', '.join(plotted)}"
            )
        else:
            self.status_label.setText(
                "Nothing to plot — check 'Data', 'Unified Fit', 'Size Dist.' or 'Simple Fits' checkbox"
            )

    def create_report(self):
        """
        Generate a Markdown report for each selected file.

        Content depends on which checkboxes are ticked:
          - 'Data'        → basic data summary (Q range, I range, N points)
          - 'Unified Fit' → fit quality metrics and level parameters (HDF5 only)
          - 'Size Dist.'  → size distribution results (HDF5 only)
        Multiple checkboxes can be active simultaneously.
        """
        selected_items = self.file_list.selectedItems()
        if not selected_items:
            QMessageBox.warning(self, "No Selection", "Please select files to report.")
            return

        show_data         = self.data_checkbox.isChecked()
        show_fit          = self.unified_fit_result_checkbox.isChecked()
        show_sizes        = self.size_dist_checkbox.isChecked()
        show_simple_fits  = self.simple_fits_checkbox.isChecked()

        if not show_data and not show_fit and not show_sizes and not show_simple_fits:
            self.status_label.setText(
                "Nothing to report — check 'Data', 'Unified Fit', 'Size Dist.' or 'Simple Fits' checkbox"
            )
            return

        file_paths = [
            os.path.join(self.current_folder, item.text())
            for item in selected_items
        ]

        error_fraction = self.state_manager.get('data_selector', 'error_fraction', 0.05)
        saved, skipped = [], []

        for file_path in file_paths:
            _, ext = os.path.splitext(file_path)
            is_hdf = ext.lower() in ['.h5', '.hdf5', '.hdf']

            data_info          = None
            fit_results        = None
            sizes_results      = None
            simple_fit_results = None

            # ── Load raw data ──────────────────────────────────────────────
            if show_data:
                try:
                    dir_path, filename = os.path.split(file_path)
                    if ext.lower() in ['.txt', '.dat']:
                        raw = readTextFile(dir_path, filename, error_fraction=error_fraction)
                    else:
                        raw = readGenericNXcanSAS(dir_path, filename)

                    if raw is not None:
                        data_info = {
                            'Q':       raw['Q'],
                            'I':       raw['Intensity'],
                            'I_error': raw.get('Error'),
                        }
                except Exception:
                    pass   # data load failure is non-fatal; section simply omitted

            # ── Load Unified Fit results (HDF5 only) ───────────────────────
            if show_fit and is_hdf:
                try:
                    fit_results = load_unified_fit_results(Path(file_path))
                except Exception:
                    pass   # no fit group or unreadable — section omitted silently

            # ── Load Size Distribution results (HDF5 only) ─────────────────
            if show_sizes and is_hdf:
                try:
                    from pyirena.io.nxcansas_sizes import load_sizes_results
                    sizes_results = load_sizes_results(Path(file_path))
                except Exception:
                    pass   # no sizes group or unreadable — section omitted silently

            # ── Load Simple Fits results (HDF5 only) ───────────────────────
            if show_simple_fits and is_hdf:
                try:
                    from pyirena.io.nxcansas_simple_fits import load_simple_fit_results
                    simple_fit_results = load_simple_fit_results(Path(file_path))
                except Exception:
                    pass   # no simple_fit_results group — section omitted silently

            # Nothing to write for this file?
            if data_info is None and fit_results is None and sizes_results is None and simple_fit_results is None:
                skipped.append(os.path.basename(file_path))
                continue

            md = _build_report(
                file_path,
                data_info=data_info,
                fit_results=fit_results,
                sizes_results=sizes_results,
                simple_fit_results=simple_fit_results,
            )
            out_path = Path(file_path).parent / (Path(file_path).stem + '_report.md')
            out_path.write_text(md, encoding='utf-8')
            saved.append(out_path.name)

            # Open in system default application for .md files
            try:
                if sys.platform == 'darwin':
                    subprocess.run(['open', str(out_path)], check=False)
                elif sys.platform == 'win32':
                    os.startfile(str(out_path))
                else:
                    subprocess.run(['xdg-open', str(out_path)], check=False)
            except Exception:
                pass   # Opening is best-effort; don't fail the whole save

        if saved:
            msg = f"Report(s) saved: {', '.join(saved)}"
            if skipped:
                msg += f"  (skipped: {', '.join(skipped)})"
            self.status_label.setText(msg)
        else:
            self.status_label.setText(
                "No reportable content found — no reports generated"
            )

    # ─────────────────────────────────────────────────────────────────────────
    def tabulate_results(self):
        """
        Collect fit results from selected HDF5 files and display them in a
        spreadsheet-like table.  The active checkboxes control which result
        sets are included:
          • 'Unified Fit' checkbox → Unified Fit scalars + per-level parameters
          • 'Size Dist.'  checkbox → Size Distribution scalars

        The table can be saved as a CSV file from within the result window.
        """
        selected_items = self.file_list.selectedItems()
        if not selected_items:
            QMessageBox.warning(self, "No Selection", "Please select files to tabulate.")
            return

        show_unified     = self.unified_fit_result_checkbox.isChecked()
        show_sizes       = self.size_dist_checkbox.isChecked()
        show_simple_fits = self.simple_fits_checkbox.isChecked()

        if not show_unified and not show_sizes and not show_simple_fits:
            self.status_label.setText(
                "Nothing to tabulate — check 'Unified Fit', 'Size Dist.' or 'Simple Fits' checkbox"
            )
            return

        file_paths = [
            os.path.join(self.current_folder, item.text())
            for item in selected_items
        ]

        from pyirena.io.nxcansas_unified import load_unified_fit_results
        from pyirena.io.nxcansas_sizes import load_sizes_results

        # ── First pass: load all data, determine max unified-fit levels ──────
        loaded = []   # list of (fname, uf_dict_or_None, sz_dict_or_None, sf_dict_or_None)
        max_levels = 0
        # Collect all SF param names seen across files (order preserved via dict)
        sf_param_names_seen: dict = {}
        sf_derived_names_seen: dict = {}

        for fp in file_paths:
            _, ext = os.path.splitext(fp)
            is_hdf = ext.lower() in ('.h5', '.hdf5', '.hdf')
            fname  = os.path.basename(fp)

            uf = None
            sz = None
            sf = None

            if is_hdf:
                if show_unified:
                    try:
                        uf = load_unified_fit_results(Path(fp))
                        max_levels = max(max_levels, len(uf.get('levels', [])))
                    except Exception:
                        pass
                if show_sizes:
                    try:
                        sz = load_sizes_results(Path(fp))
                    except Exception:
                        pass
                if show_simple_fits:
                    try:
                        from pyirena.io.nxcansas_simple_fits import load_simple_fit_results
                        sf = load_simple_fit_results(Path(fp))
                        for pname in sf.get('params', {}):
                            sf_param_names_seen[pname] = None
                        for dname in sf.get('derived', {}):
                            sf_derived_names_seen[dname] = None
                    except Exception:
                        pass

            loaded.append((fname, uf, sz, sf))

        # ── Build column headers ──────────────────────────────────────────────
        headers = ['filename']

        if show_unified:
            headers += ['UF_chi2', 'UF_background', 'UF_background_err', 'UF_num_levels']
            _uf_level_params = ['G', 'G_err', 'Rg', 'Rg_err',
                                 'B', 'B_err', 'P', 'P_err',
                                 'RgCutoff', 'ETA', 'ETA_err',
                                 'PACK', 'PACK_err', 'correlated',
                                 'Sv', 'Invariant']
            for lvl in range(1, max_levels + 1):
                for p in _uf_level_params:
                    headers.append(f'L{lvl}_{p}')

        if show_sizes:
            headers += [
                'SD_method', 'SD_shape', 'SD_chi2', 'SD_volume_fraction',
                'SD_rg', 'SD_n_iterations', 'SD_q_power',
                'SD_contrast', 'SD_aspect_ratio',
                'SD_background', 'SD_error_scale',
                'SD_power_law_B', 'SD_power_law_P',
                'SD_r_min', 'SD_r_max', 'SD_n_bins', 'SD_log_spacing',
                'SD_cursor_q_min', 'SD_cursor_q_max',
            ]

        _sf_param_names   = list(sf_param_names_seen.keys())
        _sf_derived_names = list(sf_derived_names_seen.keys())
        if show_simple_fits:
            headers += ['SF_model', 'SF_chi2', 'SF_reduced_chi2', 'SF_dof',
                        'SF_q_min', 'SF_q_max', 'SF_use_complex_bg']
            for pname in _sf_param_names:
                headers.append(f'SF_{pname}')
                headers.append(f'SF_{pname}_std')
            for dname in _sf_derived_names:
                headers.append(f'SF_derived_{dname}')

        # ── Build rows ────────────────────────────────────────────────────────
        def _fmt(v):
            """Format a value for CSV/table: round floats to 6 sig-figs."""
            if v is None:
                return None
            if isinstance(v, float):
                import math
                if math.isnan(v) or math.isinf(v):
                    return None
                return float(f'{v:.6g}')
            return v

        rows = []
        for fname, uf, sz, sf in loaded:
            row = [fname]

            if show_unified:
                if uf is not None:
                    row += [
                        _fmt(uf.get('chi_squared')),
                        _fmt(uf.get('background')),
                        _fmt(uf.get('background_err')),
                        uf.get('num_levels'),
                    ]
                    levels = uf.get('levels', [])
                    for lvl_idx in range(max_levels):
                        lv = levels[lvl_idx] if lvl_idx < len(levels) else {}
                        row += [
                            _fmt(lv.get('G')),
                            _fmt(lv.get('G_err')),
                            _fmt(lv.get('Rg')),
                            _fmt(lv.get('Rg_err')),
                            _fmt(lv.get('B')),
                            _fmt(lv.get('B_err')),
                            _fmt(lv.get('P')),
                            _fmt(lv.get('P_err')),
                            _fmt(lv.get('RgCutoff')),
                            _fmt(lv.get('ETA')),
                            _fmt(lv.get('ETA_err')),
                            _fmt(lv.get('PACK')),
                            _fmt(lv.get('PACK_err')),
                            lv.get('correlated'),
                            _fmt(lv.get('Sv')),
                            _fmt(lv.get('Invariant')),
                        ]
                else:
                    row += [None] * (4 + max_levels * 16)

            if show_sizes:
                if sz is not None:
                    row += [
                        sz.get('method'),
                        sz.get('shape'),
                        _fmt(sz.get('chi_squared')),
                        _fmt(sz.get('volume_fraction')),
                        _fmt(sz.get('rg')),
                        sz.get('n_iterations'),
                        _fmt(sz.get('q_power')),
                        _fmt(sz.get('contrast')),
                        _fmt(sz.get('aspect_ratio')),
                        _fmt(sz.get('background')),
                        _fmt(sz.get('error_scale')),
                        _fmt(sz.get('power_law_B')),
                        _fmt(sz.get('power_law_P')),
                        _fmt(sz.get('r_min')),
                        _fmt(sz.get('r_max')),
                        sz.get('n_bins'),
                        sz.get('log_spacing'),
                        _fmt(sz.get('cursor_q_min')),
                        _fmt(sz.get('cursor_q_max')),
                    ]
                else:
                    row += [None] * 19

            if show_simple_fits:
                if sf is not None:
                    sf_params_d  = sf.get('params', {})
                    sf_std_d     = sf.get('params_std', {})
                    sf_derived_d = sf.get('derived', {})
                    row += [
                        sf.get('model'),
                        _fmt(sf.get('chi_squared')),
                        _fmt(sf.get('reduced_chi_squared')),
                        sf.get('dof'),
                        _fmt(sf.get('q_min')),
                        _fmt(sf.get('q_max')),
                        sf.get('use_complex_bg'),
                    ]
                    for pname in _sf_param_names:
                        row.append(_fmt(sf_params_d.get(pname)))
                        row.append(_fmt(sf_std_d.get(pname)))
                    for dname in _sf_derived_names:
                        row.append(_fmt(sf_derived_d.get(dname)))
                else:
                    n_sf_cols = 7 + 2 * len(_sf_param_names) + len(_sf_derived_names)
                    row += [None] * n_sf_cols

            rows.append(row)

        # ── Display ───────────────────────────────────────────────────────────
        default_path = os.path.join(
            self.current_folder or '',
            'pyIrena_TableOfResults.csv',
        )

        if self.tabulate_results_window is None:
            self.tabulate_results_window = TabulateResultsWindow()

        self.tabulate_results_window.set_data(headers, rows, default_path)
        self.status_label.setText(
            f"Tabulated {len(rows)} file(s) — use 'Save as CSV' in the results window."
        )

    def launch_unified_fit(self):
        """Launch the Unified Fit model panel with selected data."""
        selected_items = self.file_list.selectedItems()

        if not selected_items:
            QMessageBox.warning(
                self,
                "No Selection",
                "Please select one or more files to analyze with Unified Fit."
            )
            return

        # Get full file paths
        file_paths = []
        for item in selected_items:
            file_path = os.path.join(self.current_folder, item.text())
            file_paths.append(file_path)

        # Load first selected file for fitting
        # (Multiple files can be loaded and fitted separately)
        file_path = file_paths[0]
        path, filename = os.path.split(file_path)
        _, ext = os.path.splitext(filename)

        error_fraction = self.state_manager.get('data_selector', 'error_fraction', 0.05)
        try:
            # Load data based on file extension
            if ext.lower() in ['.txt', '.dat']:
                data = readTextFile(path, filename, error_fraction=error_fraction)
                is_nxcansas = False
            else:
                data = readGenericNXcanSAS(path, filename)
                is_nxcansas = True  # HDF5 files loaded with NXcanSAS reader

            if data is None:
                QMessageBox.critical(
                    self,
                    "Load Error",
                    f"Could not load data from {filename}"
                )
                return

            # Create or show unified fit window
            if self.unified_fit_window is None:
                self.unified_fit_window = UnifiedFitPanel()

            # Set the data with filepath and format information
            self.unified_fit_window.set_data(
                data['Q'],
                data['Intensity'],
                data.get('Error'),
                filename,
                filepath=file_path,  # Pass full path to file
                is_nxcansas=is_nxcansas  # Pass format information
            )

            # Show the window
            self.unified_fit_window.show()
            self.unified_fit_window.raise_()
            self.unified_fit_window.activateWindow()

            self.status_label.setText(f"Opened Unified Fit for {filename}")

        except Exception as e:
            QMessageBox.critical(
                self,
                "Error",
                f"Error loading data for Unified Fit:\n{str(e)}"
            )
            self.status_label.setText(f"Error: {str(e)}")

    def launch_sizes_fit(self):
        """Launch the Size Distribution fitting panel with selected data."""
        selected_items = self.file_list.selectedItems()

        if not selected_items:
            QMessageBox.warning(
                self,
                "No Selection",
                "Please select one or more files to analyze with Size Distribution."
            )
            return

        file_path = os.path.join(self.current_folder, selected_items[0].text())
        path, filename = os.path.split(file_path)
        _, ext = os.path.splitext(filename)

        error_fraction = self.state_manager.get('data_selector', 'error_fraction', 0.05)
        try:
            if ext.lower() in ['.txt', '.dat']:
                data = readTextFile(path, filename, error_fraction=error_fraction)
                is_nxcansas = False
            else:
                data = readGenericNXcanSAS(path, filename)
                is_nxcansas = True

            if data is None:
                QMessageBox.critical(
                    self,
                    "Load Error",
                    f"Could not load data from {filename}"
                )
                return

            if self.sizes_fit_window is None:
                self.sizes_fit_window = SizesFitPanel()

            self.sizes_fit_window.set_data(
                data['Q'],
                data['Intensity'],
                data.get('Error'),
                filename,
                filepath=file_path,
                is_nxcansas=is_nxcansas,
            )

            self.sizes_fit_window.show()
            self.sizes_fit_window.raise_()
            self.sizes_fit_window.activateWindow()

            self.status_label.setText(f"Opened Size Distribution for {filename}")

        except Exception as e:
            QMessageBox.critical(
                self,
                "Error",
                f"Error loading data for Size Distribution:\n{str(e)}"
            )
            self.status_label.setText(f"Error: {str(e)}")

    # ── Batch script fitting ───────────────────────────────────────────────

    def run_unified_script(self):
        """Batch-fit all selected files with Unified Fit."""
        self._run_batch_fit('unified')

    def run_sizes_script(self):
        """Batch-fit all selected files with Size Distribution."""
        self._run_batch_fit('sizes')

    def launch_simple_fits(self):
        """Launch the Simple Fits panel with the first selected file."""
        from pyirena.gui.simple_fits_panel import SimpleFitsPanel
        selected_items = self.file_list.selectedItems()

        if not selected_items:
            QMessageBox.warning(
                self,
                "No Selection",
                "Please select a file to open in Simple Fits.",
            )
            return

        file_path = os.path.join(self.current_folder, selected_items[0].text())
        path, filename = os.path.split(file_path)
        _, ext = os.path.splitext(filename)

        error_fraction = self.state_manager.get('data_selector', 'error_fraction', 0.05)
        try:
            if ext.lower() in ['.txt', '.dat']:
                data = readTextFile(path, filename, error_fraction=error_fraction)
                is_nxcansas = False
            else:
                data = readGenericNXcanSAS(path, filename)
                is_nxcansas = True

            if data is None:
                QMessageBox.critical(
                    self, "Error",
                    f"Could not read data from file: {filename}",
                )
                return

            if self.simple_fits_window is None:
                self.simple_fits_window = SimpleFitsPanel()

            self.simple_fits_window.set_data(
                data['Q'],
                data['Intensity'],
                data.get('Error'),
                filename,
                filepath=file_path,
                is_nxcansas=is_nxcansas,
            )

            self.simple_fits_window.show()
            self.simple_fits_window.raise_()
            self.simple_fits_window.activateWindow()

            self.status_label.setText(f"Opened Simple Fits for {filename}")

        except Exception as e:
            QMessageBox.critical(
                self, "Error",
                f"Error loading data for Simple Fits:\n{str(e)}",
            )
            self.status_label.setText(f"Error: {str(e)}")

    def run_simple_fits_script(self):
        """Batch-fit all selected files with Simple Fits."""
        self._run_batch_fit('simple_fits')

    def _find_config_file(self) -> Optional[str]:
        """
        Return the path to pyirena_config.json.
        Looks first in the current folder; if not found, prompts the user.
        """
        if self.current_folder:
            candidate = os.path.join(self.current_folder, 'pyirena_config.json')
            if os.path.isfile(candidate):
                return candidate
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select pyirena_config.json",
            self.current_folder or str(Path.home()),
            "pyIrena Config (*.json);;All Files (*)",
        )
        return path or None

    def _run_batch_fit(self, tool: str):
        """Start a BatchWorker thread for the given tool ('unified' or 'sizes')."""
        if self._batch_worker is not None and self._batch_worker.isRunning():
            QMessageBox.information(
                self, "Busy", "A batch fit is already in progress — please wait."
            )
            return

        selected_items = self.file_list.selectedItems()
        if not selected_items:
            QMessageBox.warning(
                self, "No Selection", "Please select one or more files to fit."
            )
            return

        config_file = self._find_config_file()
        if not config_file:
            QMessageBox.warning(
                self,
                "No Config File",
                "No pyirena_config.json found.\n"
                "Export parameters from a GUI fit panel first, then try again.",
            )
            return

        file_paths = [
            os.path.join(self.current_folder, item.text())
            for item in selected_items
        ]
        _tool_display = {'unified': "Unified Fit", 'sizes': "Size Distribution",
                         'simple_fits': "Simple Fits"}
        tool_name = _tool_display.get(tool, tool)
        self._set_batch_status(
            f"⏳  Starting {tool_name} batch on {len(file_paths)} file(s)…", 'working'
        )

        self._batch_worker = BatchWorker(tool, file_paths, config_file, parent=self)
        self._batch_worker.progress.connect(self._on_batch_progress)
        self._batch_worker.finished.connect(self._on_batch_finished)
        self._batch_worker.start()

    def _set_batch_status(self, text: str, state: str):
        """Update the batch status label colour and text."""
        colour_map = {
            'working': ('#f39c12', 'white'),
            'done':    ('#27ae60', 'white'),
            'partial': ('#e67e22', 'white'),
            'error':   ('#e74c3c', 'white'),
        }
        bg, fg = colour_map.get(state, ('#7f8c8d', 'white'))
        self.batch_status_label.setStyleSheet(
            f"QLabel {{ padding: 5px 8px; border-radius: 4px; font-size: 12px; "
            f"background-color: {bg}; color: {fg}; }}"
        )
        self.batch_status_label.setText(text)
        self.batch_status_label.setVisible(True)

    def _on_batch_progress(self, msg: str):
        self._set_batch_status(f"⏳  {msg}", 'working')

    def _on_batch_finished(self, n_ok: int, n_fail: int, messages: list):
        total = n_ok + n_fail
        if n_fail == 0:
            summary = f"✓  Done: all {total} fit(s) succeeded"
            state = 'done'
        elif n_ok == 0:
            summary = f"✗  Done: all {total} fit(s) failed"
            state = 'error'
        else:
            summary = f"⚠  Done: {n_ok} succeeded, {n_fail} failed (out of {total})"
            state = 'partial'
        self._set_batch_status(summary, state)
        self.status_label.setText(f"Batch complete — {summary.lstrip('✓✗⚠ ')}")

        # Show detail list if there were any failures
        if n_fail > 0:
            detail = "\n".join(messages)
            QMessageBox.information(
                self,
                "Batch Fit — Details",
                f"{summary}\n\n{detail}",
            )

    def open_configure_dialog(self):
        """Open the extensible configuration dialog for data loading options."""
        current = {
            spec['key']: self.state_manager.get(
                'data_selector', spec['key'], spec.get('default')
            )
            for spec in DataSelectorConfigDialog.FIELD_SPECS
        }
        dialog = DataSelectorConfigDialog(current, parent=self)
        if dialog.exec() == QDialog.DialogCode.Accepted:
            values = dialog.get_values()
            for key, value in values.items():
                self.state_manager.set('data_selector', key, value)
            self.state_manager.save()

    def load_last_folder(self):
        """Load the last used folder from state and set it if it exists."""
        last_folder = self.state_manager.get('data_selector', 'last_folder')

        if last_folder and os.path.isdir(last_folder):
            # Folder exists, use it
            self.last_folder = last_folder
            self.current_folder = last_folder
            print(f"Restored last folder: {last_folder}")
        else:
            # Folder doesn't exist or wasn't saved, use home directory
            self.last_folder = str(Path.home())
            self.current_folder = None
            if last_folder:
                print(f"Last folder no longer exists: {last_folder}")
                print(f"Starting in home directory: {self.last_folder}")

    def save_last_folder(self, folder: str):
        """Save the current folder to state."""
        self.state_manager.set('data_selector', 'last_folder', folder)
        self.state_manager.save()

    def show_about(self):
        """Show about dialog."""
        QMessageBox.about(
            self,
            "About pyIrena",
            """<h3>pyIrena</h3>
            <p><b>Python tools for small-angle scattering data analysis</b></p>
            <p>Version: 0.1.0</p>
            <p>pyIrena provides tools for analyzing SAXS/SANS/USAXS data,
            including the Unified Fit model for hierarchical structures.</p>
            <p>Based on Irena SAS package for Igor Pro by Jan Ilavsky</p>
            <p><a href='https://github.com/jilavsky/SAXS_IgorCode'>
            Original Irena Package</a></p>
            """
        )


def main():
    """Main entry point for the data selector GUI."""
    app = QApplication(sys.argv)

    # Set application style
    app.setStyle('Fusion')

    # Create and show the main window
    window = DataSelectorPanel()
    window.show()

    sys.exit(app.exec())


if __name__ == "__main__":
    main()
