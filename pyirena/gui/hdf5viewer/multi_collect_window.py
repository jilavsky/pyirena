"""
MultiCollectWindow — floating window that shows a table of multiple values
collected from each of a set of HDF5 files.

Rows  : one per file.
Columns: File | X-value | item1 | item2 | … | itemN
"""

from __future__ import annotations

import re
from pathlib import Path

try:
    from PySide6.QtWidgets import (
        QWidget, QVBoxLayout, QTableWidget, QTableWidgetItem,
        QToolBar, QLabel, QPushButton, QAbstractItemView, QHeaderView,
        QFileDialog, QMessageBox,
    )
    from PySide6.QtCore import Qt
except ImportError:
    from PyQt6.QtWidgets import (  # type: ignore[no-redef]
        QWidget, QVBoxLayout, QTableWidget, QTableWidgetItem,
        QToolBar, QLabel, QPushButton, QAbstractItemView, QHeaderView,
        QFileDialog, QMessageBox,
    )
    from PyQt6.QtCore import Qt  # type: ignore[no-redef]


class MultiCollectWindow(QWidget):
    """
    Displays a table of multiple values collected from a set of HDF5 files.

    Parameters
    ----------
    rows : list[dict]
        Each dict has:
            "file"    : str             — filename stem
            "x_value" : float | None
            "values"  : list[object]   — one entry per item (float, str, or None)
    x_label : str
        Column header for the X column.
    item_labels : list[str]
        Column headers for each collected item.
    title : str
        Window title.
    """

    def __init__(
        self,
        rows: list[dict],
        x_label: str = "X",
        item_labels: list[str] | None = None,
        title: str = "Multi-Collect",
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent, Qt.WindowType.Window)
        self.setWindowTitle(title)
        n_items = len(item_labels) if item_labels else 0
        self.resize(max(700, 180 + 110 * n_items), 400)
        self._rows        = rows
        self._x_label     = x_label
        self._item_labels = item_labels or []
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

        # Table: File | X | col1 | col2 | ...
        n_cols = 2 + len(self._item_labels)
        self._table = QTableWidget(0, n_cols)
        headers = ["File", self._x_label] + list(self._item_labels)
        self._table.setHorizontalHeaderLabels(headers)
        self._table.setSelectionBehavior(
            QAbstractItemView.SelectionBehavior.SelectRows
        )
        self._table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        # File column stretches; all others resize to content
        self._table.horizontalHeader().setSectionResizeMode(
            0, QHeaderView.ResizeMode.Stretch
        )
        for c in range(1, n_cols):
            self._table.horizontalHeader().setSectionResizeMode(
                c, QHeaderView.ResizeMode.ResizeToContents
            )
        layout.addWidget(self._table, 1)

        # Status bar
        self._status = QLabel("")
        self._status.setStyleSheet("font-size:9pt; color:#555; padding:1px 4px;")
        layout.addWidget(self._status)

    # ── Data population ────────────────────────────────────────────────────

    def _populate(self) -> None:
        self._table.setRowCount(0)
        n_ok = 0

        for r in self._rows:
            row_i = self._table.rowCount()
            self._table.insertRow(row_i)

            fname  = str(r.get("file", "?"))
            x_val  = r.get("x_value")
            values = r.get("values", [])

            def _num_cell(val):
                s = f"{val:.6g}" if isinstance(val, float) else str(val if val is not None else "")
                it = QTableWidgetItem(s)
                it.setTextAlignment(
                    Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
                )
                return it

            def _str_cell(val):
                return QTableWidgetItem(str(val) if val is not None else "")

            self._table.setItem(row_i, 0, _str_cell(fname))
            self._table.setItem(row_i, 1, _num_cell(x_val))

            for col_i, v in enumerate(values):
                if isinstance(v, str):
                    self._table.setItem(row_i, 2 + col_i, _str_cell(v))
                else:
                    self._table.setItem(row_i, 2 + col_i, _num_cell(v))

            if any(v is not None for v in values):
                n_ok += 1

        n_files = len(self._rows)
        n_items = len(self._item_labels)
        self._status.setText(
            f"{n_ok}/{n_files} file(s) collected,  {n_items} item(s)"
        )

    # ── Export ─────────────────────────────────────────────────────────────

    def _default_path(self, ext: str) -> str:
        safe = re.sub(r"[^\w\s-]", "", self._window_title).strip().replace(" ", "_")
        return str(Path.cwd() / ((safe or "multi_collected") + ext))

    def _save_csv(self) -> None:
        filepath, _ = QFileDialog.getSaveFileName(
            self, "Save CSV", self._default_path(".csv"),
            "CSV files (*.csv);;All files (*)",
        )
        if not filepath:
            return
        if not filepath.lower().endswith(".csv"):
            filepath += ".csv"

        hdr = ["File", self._x_label] + list(self._item_labels)
        lines = [",".join(hdr)]
        for r in self._rows:
            fname  = str(r.get("file", ""))
            x_val  = r.get("x_value")
            values = r.get("values", [])
            row_data = [fname, str(x_val if x_val is not None else "")]
            for v in values:
                row_data.append(str(v) if v is not None else "")
            lines.append(",".join(row_data))

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
