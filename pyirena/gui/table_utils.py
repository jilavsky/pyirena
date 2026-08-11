"""
Shared table behaviour for every ``QTableWidget`` in the GUI.

Users arrive from Igor Pro, Origin, Excel and Matlab, where *any* table can be
selected and copied to the clipboard, sorted by clicking a column header, and
written out as text.  In pyIrena those features used to exist on two tables
only, with two different bespoke implementations, so seven other tables were
dead ends — the "collected values" window in particular could only export CSV,
forcing a file round-trip for three numbers.

This module is the single implementation, following the same pattern as
:mod:`pyirena.gui.file_filter`: one module, one behaviour, one tooltip, one
test file.  Panels call :func:`attach_table_copy` (and, where row order is not
meaningful, :func:`enable_table_sorting`) instead of writing their own.

Semantics
---------
- **Ctrl+C** copies the current selection as tab-separated text.  Tabs paste
  cleanly into Excel, Igor Pro, Origin and any text editor.
- **Ctrl+Shift+C** copies the same block with the column headers prepended.
- Right-click offers Copy / Copy with Column Headers / Copy Whole Table, plus
  "Save as CSV…" when the panel supplies a save callback.
- With nothing selected, Copy falls back to the whole table rather than
  silently putting nothing on the clipboard.
- Non-contiguous selections (ctrl-click) copy as a dense grid over just the
  selected rows and columns, so two non-adjacent columns paste as a clean
  two-column block instead of everything in between.
- Numeric columns must sort numerically, not lexically ("10" before "9" is its
  own bug report).  Build those cells with :func:`make_numeric_item` — it keeps
  the formatted display text but sorts on the underlying float.

CSV output goes through :func:`rows_to_csv_text`, which quotes with the stdlib
``csv`` module (labels containing commas used to corrupt the exported file) and
writes floats with ``%.10g`` — the same precision the ITX exporters use.
"""

from __future__ import annotations

import csv
import io
import math
from contextlib import contextmanager
from typing import Any, Callable, Iterable, List, Optional, Sequence, Tuple

from pyirena.gui._qt import (
    QAbstractItemView,
    QApplication,
    QFileDialog,
    QKeySequence,
    QMenu,
    QMessageBox,
    QShortcut,
    Qt,
    QTableWidget,
    QTableWidgetItem,
)

__all__ = [
    "COPY_TOOLTIP",
    "NumericTableWidgetItem",
    "attach_table_copy",
    "cell_text",
    "copy_table_selection",
    "enable_table_sorting",
    "grid_to_tsv",
    "make_numeric_item",
    "make_text_item",
    "numeric_sort_key",
    "populating",
    "rows_to_csv_text",
    "save_rows_as_csv",
    "save_table_as_csv",
    "table_headers",
    "table_rows",
]


COPY_TOOLTIP = (
    "Ctrl+C copies the selected cells as tab-separated text (pastes into\n"
    "Excel, Igor Pro, Origin).  Ctrl+Shift+C includes the column headers.\n"
    "Right-click for more copy options.  With nothing selected the whole\n"
    "table is copied."
)

SORT_TOOLTIP = "Click a column header to sort; click again to reverse."


# ──────────────────────────────────────────────────────────────────────────
#  Numeric-aware items
# ──────────────────────────────────────────────────────────────────────────

def numeric_sort_key(text: str) -> Optional[float]:
    """Parse a table cell's text as a number, or return None.

    Handles the decorations tables use: leading ``±``, trailing ``%``, thousands
    separators, and the em-dash / blank placeholders used for "not applicable".

    Args:
        text: Cell text as displayed.

    Returns:
        The float value, or None when the text is not numeric (including NaN,
        which must sort with the blanks rather than at an arbitrary end).
    """
    if text is None:
        return None
    s = str(text).strip()
    if not s:
        return None
    s = s.replace("±", "").replace(",", "").replace("%", "").strip()
    if not s:
        return None
    try:
        value = float(s)
    except (TypeError, ValueError):
        return None
    if math.isnan(value):
        return None
    return value


class NumericTableWidgetItem(QTableWidgetItem):
    """Table item that displays formatted text but sorts on a numeric value.

    ``QTableWidgetItem`` compares display strings, so a column of "9", "10",
    "100" sorts as 10, 100, 9.  This subclass overrides the comparison to use
    an explicit sort value (or the number parsed out of the text).  Cells whose
    text is not numeric — blanks, "—", "(ref)" — always sort after the numeric
    ones in ascending order, which is where users expect missing data.
    """

    def __init__(self, text: str = "", sort_value: Optional[float] = None) -> None:
        super().__init__(text)
        if sort_value is None:
            sort_value = numeric_sort_key(text)
        elif isinstance(sort_value, float) and math.isnan(sort_value):
            sort_value = None
        self._sort_value = sort_value

    @property
    def sort_value(self) -> Optional[float]:
        """The float this item sorts on, or None when it is not numeric."""
        return self._sort_value

    def __lt__(self, other: Any) -> bool:  # noqa: D105 - Qt operator<
        mine = self._sort_value
        theirs = getattr(other, "sort_value", None)
        if theirs is None and not isinstance(other, NumericTableWidgetItem):
            theirs = numeric_sort_key(other.text()) if hasattr(other, "text") else None
        if mine is None and theirs is None:
            return self.text() < (other.text() if hasattr(other, "text") else "")
        if mine is None:
            return False          # blanks sort last (ascending)
        if theirs is None:
            return True
        return mine < theirs


def make_numeric_item(
    value: Any,
    fmt: str = "{:.6g}",
    empty: str = "",
    align_right: bool = True,
) -> NumericTableWidgetItem:
    """Build a right-aligned, numerically sortable cell.

    Args:
        value: Number to display.  Non-numeric values (None, strings) are shown
            as-is and sort with the blanks.
        fmt: ``str.format`` spec used for floats.  Integers are shown verbatim.
        empty: Text used for None.
        align_right: Right-align the cell, as numbers are shown in every other
            scientific package.

    Returns:
        A :class:`NumericTableWidgetItem`.
    """
    if value is None:
        text: str = empty
        sort_value: Optional[float] = None
    elif isinstance(value, bool):
        text, sort_value = str(value), float(value)
    elif isinstance(value, int):
        text, sort_value = str(value), float(value)
    elif isinstance(value, float):
        if math.isnan(value):
            text, sort_value = empty, None
        else:
            text, sort_value = fmt.format(value), value
    else:
        text = str(value)
        sort_value = numeric_sort_key(text)

    item = NumericTableWidgetItem(text, sort_value)
    if align_right:
        item.setTextAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
    return item


def make_text_item(value: Any, center: bool = False) -> QTableWidgetItem:
    """Build a plain, lexically sorted cell (filenames, labels, status text)."""
    item = QTableWidgetItem("" if value is None else str(value))
    if center:
        item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
    return item


# ──────────────────────────────────────────────────────────────────────────
#  Reading a table back out
# ──────────────────────────────────────────────────────────────────────────

def table_headers(table: QTableWidget) -> List[str]:
    """Return the horizontal header labels, blank where a header is unset."""
    out = []
    for c in range(table.columnCount()):
        hdr = table.horizontalHeaderItem(c)
        out.append(hdr.text() if hdr is not None else "")
    return out


def cell_text(table: QTableWidget, row: int, col: int) -> str:
    """Text of one cell, including cells that hold a widget.

    Some tables put a combo box or check box in a column instead of an item
    (the Contrast isotope selector, for one).  Those cells have no
    ``QTableWidgetItem``, so a naive copy would export a blank column; read the
    widget's current text instead.
    """
    item = table.item(row, col)
    if item is not None:
        return item.text()
    widget = table.cellWidget(row, col)
    if widget is None:
        return ""
    for getter in ("currentText", "text", "toPlainText"):
        fn = getattr(widget, getter, None)
        if callable(fn):
            try:
                return str(fn())
            except TypeError:      # e.g. QAbstractButton.text() overloads
                continue
    return ""


def table_rows(table: QTableWidget) -> List[List[str]]:
    """Return every cell's display text as a list of rows (visual order)."""
    return [
        [cell_text(table, r, c) for c in range(table.columnCount())]
        for r in range(table.rowCount())
    ]


def grid_to_tsv(
    grid: Sequence[Sequence[str]],
    headers: Optional[Sequence[str]] = None,
) -> str:
    """Join a grid of strings into tab-separated clipboard text.

    Tabs (not commas) because that is what Excel, Igor Pro and Origin paste
    into cells without an import dialog.  Embedded tabs and newlines in a cell
    are replaced with spaces so the row/column structure survives the paste.
    """
    def clean(cell: Any) -> str:
        text = "" if cell is None else str(cell)
        return text.replace("\t", " ").replace("\r", " ").replace("\n", " ")

    lines = []
    if headers is not None:
        lines.append("\t".join(clean(h) for h in headers))
    lines.extend("\t".join(clean(cell) for cell in row) for row in grid)
    return "\n".join(lines)


def _selection_grid(table: QTableWidget) -> Tuple[List[int], List[int], List[List[str]]]:
    """Build a dense grid over the selected rows x selected columns.

    Non-contiguous selections (ctrl-clicked cells, several header-selected
    columns) copy as a compact block; unselected intersections stay blank.
    """
    items = table.selectedItems()
    if not items:
        return [], [], []

    rows = sorted({it.row() for it in items})
    cols = sorted({it.column() for it in items})
    row_index = {r: i for i, r in enumerate(rows)}
    col_index = {c: i for i, c in enumerate(cols)}

    grid = [["" for _ in cols] for _ in rows]
    for it in items:
        grid[row_index[it.row()]][col_index[it.column()]] = it.text()
    return rows, cols, grid


def copy_table_selection(
    table: QTableWidget,
    include_headers: bool = False,
    whole_table: bool = False,
) -> str:
    """Put the selection (or the whole table) on the clipboard as TSV.

    Args:
        table: Source table.
        include_headers: Prepend the header labels of the copied columns.
        whole_table: Ignore the selection and copy every row and column.

    Returns:
        The text placed on the clipboard (returned so tests can assert on it
        without touching the clipboard).
    """
    headers = table_headers(table)

    if whole_table or not table.selectedItems():
        grid = table_rows(table)
        used_headers = headers
    else:
        _rows, cols, grid = _selection_grid(table)
        used_headers = [headers[c] if c < len(headers) else "" for c in cols]

    text = grid_to_tsv(grid, used_headers if include_headers else None)
    clipboard = QApplication.clipboard()
    if clipboard is not None:
        clipboard.setText(text)
    return text


# ──────────────────────────────────────────────────────────────────────────
#  CSV
# ──────────────────────────────────────────────────────────────────────────

def _csv_value(value: Any) -> str:
    """Format one value for CSV: floats with %.10g, None as blank."""
    if value is None:
        return ""
    if isinstance(value, float):
        if math.isnan(value):
            return ""
        return f"{value:.10g}"
    return str(value)


def rows_to_csv_text(
    headers: Optional[Sequence[Any]],
    rows: Iterable[Sequence[Any]],
) -> str:
    """Render headers + rows as CSV text.

    Uses the stdlib ``csv`` writer, so a label containing a comma or a quote is
    quoted instead of splitting the row — the hand-rolled ``",".join(...)``
    writers this replaces corrupted such files silently.
    """
    buf = io.StringIO()
    writer = csv.writer(buf, lineterminator="\n")
    if headers is not None:
        writer.writerow([_csv_value(h) for h in headers])
    for row in rows:
        writer.writerow([_csv_value(v) for v in row])
    return buf.getvalue()


def save_rows_as_csv(
    parent,
    headers: Optional[Sequence[Any]],
    rows: Iterable[Sequence[Any]],
    default_path: str = "",
    title: str = "Save as CSV",
) -> Optional[str]:
    """Ask for a filename and write headers + rows as CSV.

    Returns the path written, or None if the user cancelled or the write
    failed (a message box reports the failure).
    """
    path, _ = QFileDialog.getSaveFileName(
        parent, title, default_path, "CSV files (*.csv);;All files (*)"
    )
    if not path:
        return None
    if not path.lower().endswith(".csv"):
        path += ".csv"
    try:
        with open(path, "w", newline="", encoding="utf-8") as fh:
            fh.write(rows_to_csv_text(headers, rows))
    except OSError as exc:
        QMessageBox.critical(parent, "Save failed", str(exc))
        return None
    return path


def save_table_as_csv(
    table: QTableWidget,
    parent=None,
    default_path: str = "",
    title: str = "Save as CSV",
) -> Optional[str]:
    """Write a table's displayed contents to CSV (visual row order)."""
    return save_rows_as_csv(
        parent if parent is not None else table,
        table_headers(table),
        table_rows(table),
        default_path,
        title,
    )


# ──────────────────────────────────────────────────────────────────────────
#  Attaching behaviour
# ──────────────────────────────────────────────────────────────────────────

@contextmanager
def populating(table: QTableWidget):
    """Fill a table without letting sorting reshuffle rows mid-insert.

    ``setSortingEnabled(True)`` makes Qt re-sort after every ``setItem`` call,
    which scrambles a row that is still being filled in.  Wrap population::

        with populating(self._table):
            ...  # insertRow / setItem
    """
    was_sorting = table.isSortingEnabled()
    table.setSortingEnabled(False)
    try:
        yield table
    finally:
        table.setSortingEnabled(was_sorting)


def enable_table_sorting(
    table: QTableWidget,
    default_column: Optional[int] = None,
    ascending: bool = True,
) -> QTableWidget:
    """Turn on click-a-header sorting.

    Only for tables whose row order carries no meaning.  Do **not** call it on
    tables with grouping/section rows or where the row order is itself the
    user's configuration.

    Args:
        table: Table to make sortable.
        default_column: Column to sort by immediately, if any.
        ascending: Direction for ``default_column``.
    """
    header = table.horizontalHeader()
    order = Qt.SortOrder.AscendingOrder if ascending else Qt.SortOrder.DescendingOrder

    if header is not None:
        header.setSortIndicatorShown(True)
        if not header.toolTip():
            header.setToolTip(SORT_TOOLTIP)
        if default_column is None:
            # No indicator = keep the rows in the order the panel filled them
            # (file order, usually the user's chosen file sort).  Without this
            # Qt would sort by column 0 the moment sorting is switched on, and
            # again after every repopulate.  Qt >= 6.1 also lets the user click
            # a third time to come back to that natural order.
            setter = getattr(header, "setSortIndicatorClearable", None)
            if callable(setter):
                setter(True)
            header.setSortIndicator(-1, order)

    table.setSortingEnabled(True)
    if default_column is not None:
        table.sortItems(default_column, order)
    return table


def attach_table_copy(
    table: QTableWidget,
    headers: bool = True,
    on_save_csv: Optional[Callable[[], Any]] = None,
    tooltip: bool = True,
) -> QTableWidget:
    """Give a table the standard clipboard behaviour.

    Installs the right-click menu, Ctrl+C and Ctrl+Shift+C.  Safe to call on
    any ``QTableWidget``; if the table was created with ``NoSelection`` it is
    switched to extended cell selection, because a table you cannot select is a
    table you cannot copy.

    Args:
        table: The table to equip.
        headers: Offer the copy-with-column-headers action and shortcut.
        on_save_csv: Optional callback; when given, a "Save as CSV…" entry is
            added to the context menu so the same menu covers both exits.
        tooltip: Set :data:`COPY_TOOLTIP` on the table when it has no tooltip.

    Returns:
        The same table, so calls can be chained onto the constructor.
    """
    if table.selectionMode() == QAbstractItemView.SelectionMode.NoSelection:
        table.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectItems)

    table.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)

    def show_menu(pos) -> None:
        # Right-clicking a cell outside the selection acts on that cell, which
        # is what every spreadsheet does.
        item_at = table.itemAt(pos)
        if item_at is not None and item_at not in table.selectedItems():
            table.setCurrentItem(item_at)

        menu = QMenu(table)
        act_copy = menu.addAction("Copy")
        act_copy_hdr = menu.addAction("Copy with Column Headers") if headers else None
        act_copy_all = menu.addAction("Copy Whole Table")
        act_csv = None
        if on_save_csv is not None:
            menu.addSeparator()
            act_csv = menu.addAction("Save as CSV…")

        chosen = menu.exec(table.viewport().mapToGlobal(pos))
        if chosen is None:
            return
        if chosen is act_copy:
            copy_table_selection(table, include_headers=False)
        elif act_copy_hdr is not None and chosen is act_copy_hdr:
            copy_table_selection(table, include_headers=True)
        elif chosen is act_copy_all:
            copy_table_selection(table, include_headers=True, whole_table=True)
        elif act_csv is not None and chosen is act_csv:
            on_save_csv()

    # Attaching twice must not show the menu twice; only our own previous
    # handler is disconnected, so a panel's unrelated connections survive.
    previous = getattr(table, "_pyirena_copy_menu", None)
    if previous is not None:
        try:
            table.customContextMenuRequested.disconnect(previous)
        except (RuntimeError, TypeError):
            pass
    table.customContextMenuRequested.connect(show_menu)
    table._pyirena_copy_menu = show_menu

    shortcuts = []
    sc_copy = QShortcut(QKeySequence.StandardKey.Copy, table)
    sc_copy.setContext(Qt.ShortcutContext.WidgetWithChildrenShortcut)
    sc_copy.activated.connect(lambda: copy_table_selection(table, include_headers=False))
    shortcuts.append(sc_copy)

    if headers:
        sc_hdr = QShortcut(QKeySequence("Ctrl+Shift+C"), table)
        sc_hdr.setContext(Qt.ShortcutContext.WidgetWithChildrenShortcut)
        sc_hdr.activated.connect(lambda: copy_table_selection(table, include_headers=True))
        shortcuts.append(sc_hdr)

    # Keep the shortcuts alive for the table's lifetime and let a second
    # attach_table_copy() call replace rather than duplicate them.
    for old in getattr(table, "_pyirena_copy_shortcuts", []):
        old.setEnabled(False)
        old.setParent(None)
    table._pyirena_copy_shortcuts = shortcuts

    if tooltip and not table.toolTip():
        table.setToolTip(COPY_TOOLTIP)

    return table
