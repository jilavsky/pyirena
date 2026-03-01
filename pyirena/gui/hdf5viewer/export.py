"""
Export utilities for GraphWindow.

Provides: JPEG, PNG, CSV, HDF5, ITX (Igor Pro text format), matplotlib figure.
"""

from __future__ import annotations

import os
import re
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

try:
    from PySide6.QtWidgets import QFileDialog, QMessageBox
except ImportError:
    from PyQt6.QtWidgets import QFileDialog, QMessageBox  # type: ignore[no-redef]

if TYPE_CHECKING:
    from .graph_window import GraphWindow


def _default_save_path(title: str, ext: str) -> str:
    """Return CWD / sanitised_title + ext as the default save path."""
    safe = re.sub(r"[^\w\s-]", "", title).strip().replace(" ", "_")
    return str(Path.cwd() / ((safe or "graph") + ext))


# ── JPEG / PNG ─────────────────────────────────────────────────────────────

def save_jpeg(gw: "GraphWindow", filepath: str | None = None) -> bool:
    if filepath is None:
        filepath, _ = QFileDialog.getSaveFileName(
            gw, "Save as JPEG", _default_save_path(gw.get_title(), ".jpg"),
            "JPEG images (*.jpg *.jpeg);;All files (*)",
        )
    if not filepath:
        return False
    if not filepath.lower().endswith((".jpg", ".jpeg")):
        filepath += ".jpg"
    gw.grab().save(filepath, "JPEG", 95)
    return True


def save_png(gw: "GraphWindow", filepath: str | None = None) -> bool:
    if filepath is None:
        filepath, _ = QFileDialog.getSaveFileName(
            gw, "Save as PNG", _default_save_path(gw.get_title(), ".png"),
            "PNG images (*.png);;All files (*)",
        )
    if not filepath:
        return False
    if not filepath.lower().endswith(".png"):
        filepath += ".png"
    gw.grab().save(filepath, "PNG")
    return True


# ── CSV ────────────────────────────────────────────────────────────────────

def save_csv(gw: "GraphWindow", filepath: str | None = None) -> bool:
    curves = gw.get_curves()
    if not curves:
        QMessageBox.warning(gw, "No data", "No curves to export.")
        return False

    if filepath is None:
        filepath, _ = QFileDialog.getSaveFileName(
            gw, "Save as CSV", _default_save_path(gw.get_title(), ".csv"),
            "CSV files (*.csv);;All files (*)",
        )
    if not filepath:
        return False
    if not filepath.lower().endswith(".csv"):
        filepath += ".csv"

    # Build header and columns
    lines = []
    header_parts = []
    columns = []

    for curve in curves:
        lbl = curve["label"].replace(",", ";")
        x   = np.asarray(curve["x"], float)
        y   = np.asarray(curve["y"], float)
        header_parts.append(f"X_{lbl}")
        header_parts.append(f"Y_{lbl}")
        columns.append(x)
        columns.append(y)
        if curve.get("yerr") is not None:
            header_parts.append(f"Yerr_{lbl}")
            columns.append(np.asarray(curve["yerr"], float))
        if curve.get("xerr") is not None:
            header_parts.append(f"Xerr_{lbl}")
            columns.append(np.asarray(curve["xerr"], float))

    lines.append(",".join(header_parts))
    n_rows = max(len(c) for c in columns)
    for i in range(n_rows):
        row = []
        for col in columns:
            if i < len(col):
                row.append(f"{col[i]:.8g}")
            else:
                row.append("")
        lines.append(",".join(row))

    try:
        with open(filepath, "w", encoding="utf-8") as f:
            f.write("\n".join(lines))
        return True
    except Exception as exc:
        QMessageBox.critical(gw, "Save failed", str(exc))
        return False


# ── HDF5 ───────────────────────────────────────────────────────────────────

def save_hdf5(gw: "GraphWindow", filepath: str | None = None) -> bool:
    import h5py

    curves = gw.get_curves()
    if not curves:
        QMessageBox.warning(gw, "No data", "No curves to export.")
        return False

    if filepath is None:
        filepath, _ = QFileDialog.getSaveFileName(
            gw, "Save as HDF5", _default_save_path(gw.get_title(), ".h5"),
            "HDF5 files (*.h5 *.hdf5);;All files (*)",
        )
    if not filepath:
        return False

    try:
        with h5py.File(filepath, "w") as f:
            f.attrs["NX_class"] = "NXroot"
            entry = f.require_group("entry")
            entry.attrs["NX_class"] = "NXentry"
            for i, curve in enumerate(curves):
                name = f"curve_{i+1:02d}"
                grp = entry.require_group(name)
                grp.attrs["NX_class"] = "NXdata"
                grp.attrs["signal"] = "y"
                grp.attrs["axes"] = "x"
                grp.attrs["label"] = curve["label"]
                grp.create_dataset("x", data=np.asarray(curve["x"], float), compression="gzip")
                grp.create_dataset("y", data=np.asarray(curve["y"], float), compression="gzip")
                if curve.get("yerr") is not None:
                    grp.create_dataset("yerr", data=np.asarray(curve["yerr"], float), compression="gzip")
                if curve.get("xerr") is not None:
                    grp.create_dataset("xerr", data=np.asarray(curve["xerr"], float), compression="gzip")
        return True
    except Exception as exc:
        QMessageBox.critical(gw, "Save failed", str(exc))
        return False


# ── ITX (Igor Pro text format) ─────────────────────────────────────────────

def save_itx(gw: "GraphWindow", filepath: str | None = None) -> bool:
    """Export curves as Igor Pro Text (.itx) file with full formatting."""
    curves = gw.get_curves()
    if not curves:
        QMessageBox.warning(gw, "No data", "No curves to export.")
        return False

    if filepath is None:
        filepath, _ = QFileDialog.getSaveFileName(
            gw, "Save as Igor Pro ITX", _default_save_path(gw.get_title(), ".itx"),
            "Igor Pro Text (*.itx);;All files (*)",
        )
    if not filepath:
        return False
    if not filepath.lower().endswith(".itx"):
        filepath += ".itx"

    def _safe_name(label: str) -> str:
        """Convert label to a valid Igor wave name (max 31 chars, alphanumeric+_)."""
        import re
        name = re.sub(r"[^A-Za-z0-9_]", "_", label)
        if name and name[0].isdigit():
            name = "w_" + name
        return name[:31] or "wave"

    def _hex_to_igor(hex_color: str) -> tuple[int, int, int]:
        """Convert #rrggbb hex to Igor Pro 0-65535 RGB tuple."""
        h = hex_color.lstrip("#")
        if len(h) == 6:
            r = int(h[0:2], 16)
            g = int(h[2:4], 16)
            b = int(h[4:6], 16)
        else:
            r = g = b = 0
        return r * 257, g * 257, b * 257

    lines = ["IGOR"]
    # (x_name, y_name, label, color) tuples for formatting commands
    wave_names = []

    for i, curve in enumerate(curves):
        lbl    = curve["label"]
        x_arr  = np.asarray(curve["x"], float)
        y_arr  = np.asarray(curve["y"], float)
        color  = curve.get("style", {}).get("color", "#000000")
        suffix = f"_{i+1:02d}" if len(curves) > 1 else ""
        x_name = _safe_name(f"X_{lbl}{suffix}")
        y_name = _safe_name(f"Y_{lbl}{suffix}")
        wave_names.append((x_name, y_name, lbl, color))

        # X wave
        lines.append(f"WAVES/D  {x_name}")
        lines.append("BEGIN")
        for v in x_arr:
            lines.append(f"  {v:.10g}")
        lines.append("END")

        # Y wave
        lines.append(f"WAVES/D  {y_name}")
        lines.append("BEGIN")
        for v in y_arr:
            lines.append(f"  {v:.10g}")
        lines.append("END")

        # Optional Yerr wave
        if curve.get("yerr") is not None:
            e_name = _safe_name(f"Yerr_{lbl}{suffix}")
            lines.append(f"WAVES/D  {e_name}")
            lines.append("BEGIN")
            for v in np.asarray(curve["yerr"], float):
                lines.append(f"  {v:.10g}")
            lines.append("END")

    # ── Display / graph commands ───────────────────────────────────────────
    lines.append("")
    for j, (xn, yn, lbl, _color) in enumerate(wave_names):
        if j == 0:
            lines.append(f'X Display {yn} vs {xn} as "{lbl}"')
        else:
            lines.append(f"X AppendToGraph {yn} vs {xn}")

    # ── Log / linear axes ─────────────────────────────────────────────────
    if gw.is_log_x():
        lines.append("X ModifyGraph log(bottom)=1")
    if gw.is_log_y():
        lines.append("X ModifyGraph log(left)=1")

    # ── Curve colors ──────────────────────────────────────────────────────
    for _xn, yn, _lbl, color in wave_names:
        r, g, b = _hex_to_igor(color)
        lines.append(f"X ModifyGraph rgb({yn})=({r},{g},{b})")

    # ── Axis labels ───────────────────────────────────────────────────────
    x_label = gw.get_x_label()
    y_label = gw.get_y_label()
    if x_label:
        lines.append(f'X Label bottom "{x_label}"')
    if y_label:
        lines.append(f'X Label left "{y_label}"')

    # ── Window title ──────────────────────────────────────────────────────
    title = gw.get_title()
    if title:
        lines.append(f'X TextBox/C/N=title0/A=MC/X=0/Y=5 "{title}"')

    # ── Legend ────────────────────────────────────────────────────────────
    legend_parts = []
    for _xn, yn, lbl, _color in wave_names:
        legend_parts.append(f"\\\\s({yn}) {lbl}")
    if legend_parts:
        legend_text = "\\r".join(legend_parts)
        lines.append(f'X Legend/C/N=text0 "{legend_text}"')

    try:
        with open(filepath, "w", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")
        return True
    except Exception as exc:
        QMessageBox.critical(gw, "Save failed", str(exc))
        return False


# ── Matplotlib ─────────────────────────────────────────────────────────────

def open_matplotlib(gw: "GraphWindow") -> None:
    """Reproduce the graph as a matplotlib figure in a new window."""
    curves = gw.get_curves()
    if not curves:
        QMessageBox.warning(gw, "No data", "No curves to display.")
        return

    try:
        import matplotlib
        matplotlib.use("TkAgg")
    except Exception:
        try:
            import matplotlib
        except ImportError:
            QMessageBox.warning(gw, "matplotlib not found",
                                "Install matplotlib to use this feature.")
            return

    import matplotlib.pyplot as plt

    log_x = gw.is_log_x()
    log_y = gw.is_log_y()
    title = gw.get_title()
    x_label = gw.get_x_label()
    y_label = gw.get_y_label()

    fig, ax = plt.subplots(figsize=(8, 6))

    for curve in curves:
        x = np.asarray(curve["x"], float)
        y = np.asarray(curve["y"], float)
        style = curve.get("style", {})
        color = style.get("color", None)
        width = style.get("width", 1.5)
        symbol = style.get("symbol", None)
        line_style = "-" if symbol is None else "none"
        marker = _pg_symbol_to_mpl(symbol)

        yerr = curve.get("yerr")
        if yerr is not None and not log_y:
            ax.errorbar(x, y, yerr=np.asarray(yerr, float),
                        label=curve["label"], color=color,
                        linewidth=width, marker=marker, fmt="-")
        else:
            ax.plot(x, y, label=curve["label"], color=color,
                    linewidth=width, marker=marker, linestyle=line_style)

    if log_x:
        ax.set_xscale("log")
    if log_y:
        ax.set_yscale("log")

    ax.set_xlabel(x_label or "X")
    ax.set_ylabel(y_label or "Y")
    if title:
        ax.set_title(title)
    ax.legend(loc="best")
    ax.grid(True, which="both", linestyle="--", alpha=0.4)
    fig.tight_layout()
    plt.show()


def _pg_symbol_to_mpl(symbol: str | None) -> str:
    """Convert a pyqtgraph symbol string to a matplotlib marker."""
    _map = {
        "o": "o", "s": "s", "t": "^", "d": "D", "+": "+",
        "x": "x", "star": "*", "p": "p", "h": "h",
    }
    return _map.get(symbol or "", "None")
