"""
MultiCollectWindow — floating window that shows a table of multiple values
collected from each of a set of HDF5 files.

Rows  : one per file.
Columns: File | X-value | item1 | item2 | … | itemN
"""

from __future__ import annotations

import re
from pathlib import Path

from pyirena.gui._qt import (
    QAbstractItemView,
    QHeaderView,
    QLabel,
    QPushButton,
    Qt,
    QTableWidget,
    QToolBar,
    QVBoxLayout,
    QWidget,
)
from pyirena.gui.plot_export import save_widget_image
from pyirena.gui.table_utils import (
    attach_table_copy,
    enable_table_sorting,
    make_numeric_item,
    make_text_item,
    populating,
    save_rows_as_csv,
)


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
        # Clipboard copy + numeric column sorting, shared with every other table.
        attach_table_copy(self._table, on_save_csv=self._save_csv)
        enable_table_sorting(self._table)
        layout.addWidget(self._table, 1)

        # Status bar
        self._status = QLabel("")
        self._status.setStyleSheet("font-size:9pt; color:#555; padding:1px 4px;")
        layout.addWidget(self._status)

    # ── Data population ────────────────────────────────────────────────────

    def _populate(self) -> None:
        self._table.setRowCount(0)
        n_ok = 0

        # Sorting must be off while filling (Qt re-sorts on every setItem).
        with populating(self._table):
            for r in self._rows:
                row_i = self._table.rowCount()
                self._table.insertRow(row_i)

                fname  = str(r.get("file", "?"))
                x_val  = r.get("x_value")
                values = r.get("values", [])

                self._table.setItem(row_i, 0, make_text_item(fname))
                self._table.setItem(row_i, 1, make_numeric_item(x_val))

                for col_i, v in enumerate(values):
                    if isinstance(v, str):
                        self._table.setItem(row_i, 2 + col_i, make_text_item(v))
                    else:
                        self._table.setItem(row_i, 2 + col_i, make_numeric_item(v))

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
        # Shared writer: item labels containing a comma are quoted rather than
        # silently splitting the row, and floats keep full precision.
        headers = ["File", self._x_label] + list(self._item_labels)
        rows = [
            [str(r.get("file", "")), r.get("x_value"), *r.get("values", [])]
            for r in self._rows
        ]
        save_rows_as_csv(self, headers, rows, self._default_path(".csv"), "Save CSV")

    def _save_jpeg(self) -> None:
        """Toolbar button: save the whole window as a PNG or JPEG image."""
        stem = re.sub(r"[^\w\s-]", "", self._window_title).strip().replace(" ", "_")
        save_widget_image(self, self, stem or "multi_collected")
