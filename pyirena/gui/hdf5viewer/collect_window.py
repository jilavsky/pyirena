"""
CollectWindow — floating window that shows collected 0D values from multiple
HDF5 files as both a table and a pyqtgraph scatter plot.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pyqtgraph as pg

try:
    from PySide6.QtWidgets import (
        QWidget, QVBoxLayout, QSplitter, QTableWidget, QTableWidgetItem,
        QToolBar, QLabel, QPushButton, QAbstractItemView, QHeaderView,
        QFileDialog, QMessageBox,
    )
    from PySide6.QtCore import Qt
except ImportError:
    from PyQt6.QtWidgets import (  # type: ignore[no-redef]
        QWidget, QVBoxLayout, QSplitter, QTableWidget, QTableWidgetItem,
        QToolBar, QLabel, QPushButton, QAbstractItemView, QHeaderView,
        QFileDialog, QMessageBox,
    )
    from PyQt6.QtCore import Qt  # type: ignore[no-redef]

from . import export as _export


class CollectWindow(QWidget):
    """
    Displays a table and plot of values collected from a set of HDF5 files.

    Parameters
    ----------
    rows : list[dict]
        Each dict has keys: "file", "x_value", "y_value", "y_error" (or None).
    x_label : str
        Label for the X axis.
    y_label : str
        Label for the Y axis / collected quantity.
    title : str
        Window title.
    """

    def __init__(
        self,
        rows: list[dict],
        x_label: str = "X",
        y_label: str = "Value",
        title: str = "Collected Values",
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent, Qt.WindowType.Window)
        self.setWindowTitle(title)
        self.resize(700, 500)
        self._rows = rows
        self._x_label = x_label
        self._y_label = y_label
        self._window_title = title

        self._build_ui()
        self._populate()

    # ── UI construction ────────────────────────────────────────────────────

    def _build_ui(self) -> None:
        layout = QVBoxLayout(self)
        layout.setContentsMargins(2, 2, 2, 2)
        layout.setSpacing(0)

        # Toolbar
        tb = QToolBar()
        tb.setMovable(False)

        csv_btn = QPushButton("Save CSV")
        csv_btn.setFixedWidth(75)
        csv_btn.clicked.connect(self._save_csv)
        tb.addWidget(csv_btn)

        jpeg_btn = QPushButton("Save JPEG")
        jpeg_btn.setFixedWidth(80)
        jpeg_btn.clicked.connect(self._save_jpeg)
        tb.addWidget(jpeg_btn)

        layout.addWidget(tb)

        # Splitter: table (top) + plot (bottom)
        splitter = QSplitter(Qt.Orientation.Vertical)

        # Table
        self._table = QTableWidget(0, 4)
        self._table.setHorizontalHeaderLabels(
            ["File", self._x_label, self._y_label, "±Error"]
        )
        self._table.setSelectionBehavior(
            QAbstractItemView.SelectionBehavior.SelectRows
        )
        self._table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self._table.horizontalHeader().setSectionResizeMode(
            0, QHeaderView.ResizeMode.Stretch
        )
        splitter.addWidget(self._table)

        # Plot
        glw = pg.GraphicsLayoutWidget()
        glw.setBackground("w")
        self._plot: pg.PlotItem = glw.addPlot()
        self._plot.setLabel("left",   self._y_label)
        self._plot.setLabel("bottom", self._x_label)
        self._plot.showGrid(x=True, y=True, alpha=0.3)
        self._plot.getAxis('bottom').enableAutoSIPrefix(False)
        self._plot.getAxis('left').enableAutoSIPrefix(False)
        splitter.addWidget(glw)

        splitter.setSizes([180, 280])
        layout.addWidget(splitter, 1)

        # Status
        self._status = QLabel("")
        self._status.setStyleSheet("font-size:9pt; color:#555; padding:1px 4px;")
        layout.addWidget(self._status)

    # ── Data population ────────────────────────────────────────────────────

    def _populate(self) -> None:
        rows = self._rows
        self._table.setRowCount(0)

        xs, ys, yes = [], [], []
        for r in rows:
            row_i = self._table.rowCount()
            self._table.insertRow(row_i)

            fname = Path(r.get("file", "?")).name
            xv    = r.get("x_value")
            yv    = r.get("y_value")
            ye    = r.get("y_error")

            def _cell(val):
                s = f"{val:.6g}" if isinstance(val, float) else str(val or "")
                item = QTableWidgetItem(s)
                item.setTextAlignment(
                    Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
                )
                return item

            self._table.setItem(row_i, 0, QTableWidgetItem(fname))
            self._table.setItem(row_i, 1, _cell(xv))
            self._table.setItem(row_i, 2, _cell(yv))
            self._table.setItem(row_i, 3, _cell(ye))

            if xv is not None and yv is not None:
                xs.append(float(xv))
                ys.append(float(yv))
                yes.append(float(ye) if ye is not None else None)

        self._status.setText(f"{len(rows)} file(s) collected")

        if not xs:
            return

        xs_arr = np.array(xs)
        ys_arr = np.array(ys)

        # Scatter plot
        scatter = pg.ScatterPlotItem(
            xs_arr, ys_arr,
            pen=pg.mkPen("#2980b9"),
            brush=pg.mkBrush("#2980b9"),
            size=8,
        )
        self._plot.addItem(scatter)

        # Error bars (NaN-separated vertical lines)
        if any(e is not None for e in yes):
            ex, ey = [], []
            for xi, yi, ei in zip(xs_arr, ys_arr, yes):
                if ei is None:
                    continue
                ey += [yi - ei, yi + ei, float("nan")]
                ex += [xi, xi, float("nan")]
            if ex:
                self._plot.plot(
                    list(ex), list(ey),
                    pen=pg.mkPen("#2980b9", width=1.5),
                    connect="finite",
                )

        # Auto-range
        self._plot.autoRange()

    # ── Export ─────────────────────────────────────────────────────────────

    def _default_path(self, ext: str) -> str:
        safe = re.sub(r"[^\w\s-]", "", self._window_title).strip().replace(" ", "_")
        return str(Path.cwd() / ((safe or "collected_values") + ext))

    def _save_csv(self) -> None:
        filepath, _ = QFileDialog.getSaveFileName(
            self, "Save CSV", self._default_path(".csv"),
            "CSV files (*.csv);;All files (*)",
        )
        if not filepath:
            return
        if not filepath.lower().endswith(".csv"):
            filepath += ".csv"

        lines = [f"File,{self._x_label},{self._y_label},Error"]
        for r in self._rows:
            fname = Path(r.get("file", "?")).name
            xv = r.get("x_value", "")
            yv = r.get("y_value", "")
            ye = r.get("y_error", "")
            lines.append(f"{fname},{xv},{yv},{ye}")

        try:
            with open(filepath, "w", encoding="utf-8") as f:
                f.write("\n".join(lines))
        except Exception as exc:
            QMessageBox.critical(self, "Save failed", str(exc))

    def _save_jpeg(self) -> None:
        filepath, _ = QFileDialog.getSaveFileName(
            self, "Save as JPEG", self._default_path(".jpg"),
            "JPEG images (*.jpg);;All files (*)",
        )
        if not filepath:
            return
        if not filepath.lower().endswith((".jpg", ".jpeg")):
            filepath += ".jpg"
        self.grab().save(filepath, "JPEG", 95)
