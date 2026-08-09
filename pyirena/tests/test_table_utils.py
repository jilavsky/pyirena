"""
Tests for :mod:`pyirena.gui.table_utils` — the shared table behaviour used by
every ``QTableWidget`` in the GUI.

The whole module is skipped when no Qt binding is installed (the plain `test`
CI job); the `test-gui` job runs it with ``QT_QPA_PLATFORM=offscreen``, which
is enough for widgets, selections, sorting and the clipboard.

What is checked here is what users report when it breaks:
  * numeric columns sort 9 < 10 < 100, not "10" < "9";
  * blanks/placeholders sort last instead of in the middle of the numbers;
  * copying a non-contiguous selection produces a dense block;
  * copy with nothing selected still puts the whole table on the clipboard;
  * CSV quoting survives labels containing commas and quotes.
"""

import pytest


def _require_qt() -> None:
    """Skip this module when no Qt binding is installed.

    ``pytest.importorskip`` is not enough: ``pyirena.gui._qt`` raises a plain
    ImportError carrying an installation hint, and pytest >= 8.2 treats that as
    a *broken* module rather than a missing one — so the plain (no-GUI) CI job
    would report a collection error instead of a skip.
    """
    try:
        import pyirena.gui._qt  # noqa: F401
    except ImportError:
        pytest.skip("Qt (PySide6/PyQt6) not available", allow_module_level=True)


_require_qt()

from pyirena.gui._qt import QApplication, Qt, QTableWidget, QTableWidgetItem  # noqa: E402
from pyirena.gui.table_utils import (  # noqa: E402
    NumericTableWidgetItem,
    attach_table_copy,
    cell_text,
    copy_table_selection,
    enable_table_sorting,
    grid_to_tsv,
    make_numeric_item,
    make_text_item,
    numeric_sort_key,
    populating,
    rows_to_csv_text,
    table_headers,
    table_rows,
)


@pytest.fixture(scope="module")
def qapp():
    """A single QApplication for the module (Qt allows only one)."""
    app = QApplication.instance() or QApplication([])
    yield app


@pytest.fixture
def table(qapp):
    """A 3x3 table: text column, numeric column, numeric column with a blank."""
    t = QTableWidget(3, 3)
    t.setHorizontalHeaderLabels(["File", "Rg", "Power"])
    values = [
        ("sample_100C", 100.0, 4.0),
        ("sample_9C", 9.0, None),
        ("sample_10C", 10.0, 3.5),
    ]
    for r, (name, rg, power) in enumerate(values):
        t.setItem(r, 0, make_text_item(name))
        t.setItem(r, 1, make_numeric_item(rg))
        t.setItem(r, 2, make_numeric_item(power))
    return t


# ── Pure formatting logic ────────────────────────────────────────────────

def test_numeric_sort_key_parses_decorated_numbers():
    assert numeric_sort_key("12.5") == 12.5
    assert numeric_sort_key("  1e-3 ") == pytest.approx(1e-3)
    assert numeric_sort_key("± 0.25") == 0.25
    assert numeric_sort_key("1,234") == 1234.0
    assert numeric_sort_key("45 %") == 45.0


def test_numeric_sort_key_rejects_non_numbers():
    for text in ("", "   ", "—", "— (ref)", "sample_10C", None, "nan"):
        assert numeric_sort_key(text) is None


def test_numeric_item_sorts_numerically_not_lexically():
    items = [NumericTableWidgetItem(t) for t in ("100", "9", "10")]
    assert [it.text() for it in sorted(items)] == ["9", "10", "100"]


def test_blank_items_sort_after_numbers():
    numeric = NumericTableWidgetItem("1.0")
    blank = NumericTableWidgetItem("")
    assert numeric < blank
    assert not (blank < numeric)


def test_make_numeric_item_keeps_display_text_and_sort_value():
    item = make_numeric_item(0.000123456789)
    assert item.text() == "0.000123457"          # %.6g display
    assert item.sort_value == pytest.approx(0.000123456789)   # full precision sort

    empty = make_numeric_item(None)
    assert empty.text() == ""
    assert empty.sort_value is None

    nan = make_numeric_item(float("nan"))
    assert nan.sort_value is None


def test_grid_to_tsv_flattens_embedded_tabs_and_newlines():
    text = grid_to_tsv([["a\tb", "c\nd"]], headers=["h1", "h2"])
    assert text == "h1\th2\na b\tc d"


def test_rows_to_csv_text_quotes_commas_and_quotes():
    text = rows_to_csv_text(["label", "value"], [['Rg, nm', 1.5], ['say "hi"', None]])
    assert text == 'label,value\n"Rg, nm",1.5\n"say ""hi""",\n'


def test_rows_to_csv_text_uses_10_significant_digits():
    text = rows_to_csv_text(None, [[1 / 3]])
    assert text.strip() == "0.3333333333"


# ── Qt behaviour ─────────────────────────────────────────────────────────

def test_table_headers_and_rows(table):
    assert table_headers(table) == ["File", "Rg", "Power"]
    assert table_rows(table)[1] == ["sample_9C", "9", ""]


def test_cell_text_reads_cell_widgets(qapp):
    from pyirena.gui._qt import QComboBox

    t = QTableWidget(1, 2)
    t.setItem(0, 0, make_text_item("H"))
    combo = QComboBox()
    combo.addItem("2  [b_c = 6.671 fm]")
    t.setCellWidget(0, 1, combo)
    assert cell_text(t, 0, 0) == "H"
    assert cell_text(t, 0, 1).startswith("2 ")


def test_copy_whole_table_when_nothing_selected(table):
    attach_table_copy(table)
    text = copy_table_selection(table, include_headers=True)
    lines = text.split("\n")
    assert lines[0] == "File\tRg\tPower"
    assert lines[1] == "sample_100C\t100\t4"
    assert len(lines) == 4
    assert QApplication.clipboard().text() == text


def test_copy_non_contiguous_selection_is_a_dense_block(table):
    attach_table_copy(table)
    # Columns 0 and 2 of rows 0 and 2 — the classic ctrl-click selection.
    for r in (0, 2):
        for c in (0, 2):
            table.item(r, c).setSelected(True)
    text = copy_table_selection(table, include_headers=True)
    assert text.split("\n") == [
        "File\tPower",
        "sample_100C\t4",
        "sample_10C\t3.5",
    ]


def test_attach_table_copy_makes_unselectable_table_selectable(qapp):
    from pyirena.gui._qt import QAbstractItemView

    t = QTableWidget(1, 1)
    t.setSelectionMode(QAbstractItemView.SelectionMode.NoSelection)
    attach_table_copy(t)
    assert t.selectionMode() != QAbstractItemView.SelectionMode.NoSelection
    assert t.contextMenuPolicy() == Qt.ContextMenuPolicy.CustomContextMenu
    # Shortcuts are parented to the table so they survive as long as it does.
    assert len(t._pyirena_copy_shortcuts) == 2


def test_attach_table_copy_twice_does_not_duplicate_shortcuts(table):
    attach_table_copy(table)
    attach_table_copy(table)
    assert len(table._pyirena_copy_shortcuts) == 2


def test_enable_table_sorting_sorts_numeric_column(table):
    enable_table_sorting(table)
    table.sortItems(1, Qt.SortOrder.AscendingOrder)
    assert [table.item(r, 1).text() for r in range(3)] == ["9", "10", "100"]
    table.sortItems(1, Qt.SortOrder.DescendingOrder)
    assert [table.item(r, 1).text() for r in range(3)] == ["100", "10", "9"]


def test_sorting_puts_blanks_last(table):
    enable_table_sorting(table)
    table.sortItems(2, Qt.SortOrder.AscendingOrder)
    assert [table.item(r, 2).text() for r in range(3)] == ["3.5", "4", ""]


def test_text_column_still_sorts_lexically(table):
    enable_table_sorting(table)
    table.sortItems(0, Qt.SortOrder.AscendingOrder)
    assert [table.item(r, 0).text() for r in range(3)] == [
        "sample_100C", "sample_10C", "sample_9C",
    ]


def test_enable_table_sorting_does_not_reorder_rows_by_itself(table):
    """Switching sorting on must leave the panel's own row order alone.

    Collect windows list files in the order the user sorted them in the file
    browser; sorting by column 0 the moment the feature is enabled would look
    like the table shuffling itself.
    """
    before = [table.item(r, 0).text() for r in range(3)]
    enable_table_sorting(table)
    assert [table.item(r, 0).text() for r in range(3)] == before
    with populating(table):
        pass
    assert [table.item(r, 0).text() for r in range(3)] == before


def test_enable_table_sorting_with_default_column_sorts_immediately(table):
    enable_table_sorting(table, default_column=1)
    assert [table.item(r, 1).text() for r in range(3)] == ["9", "10", "100"]


def test_populating_restores_sorting_state(table):
    enable_table_sorting(table)
    with populating(table):
        assert not table.isSortingEnabled()
        table.insertRow(0)
        table.setItem(0, 1, make_numeric_item(1.0))
    assert table.isSortingEnabled()
    assert table.rowCount() == 4


def test_populating_keeps_each_row_intact_when_sorting_is_active(qapp):
    """Filling a sorted table must not tear rows apart mid-insert.

    Without ``populating()`` Qt re-sorts after every ``setItem`` call, so a row
    whose name cell is set but whose value cell is not yet set moves, and the
    value then lands on a different row.  Afterwards the user's chosen sort is
    re-applied to the finished rows, which is what they expect.
    """
    t = QTableWidget(0, 2)
    enable_table_sorting(t)
    t.sortItems(1, Qt.SortOrder.AscendingOrder)
    with populating(t):
        for name, value in [("a", 3.0), ("b", 1.0), ("c", 2.0)]:
            r = t.rowCount()
            t.insertRow(r)
            t.setItem(r, 0, make_text_item(name))
            t.setItem(r, 1, make_numeric_item(value))

    pairs = {t.item(r, 0).text(): t.item(r, 1).text() for r in range(t.rowCount())}
    assert pairs == {"a": "3", "b": "1", "c": "2"}
    # Sorting stays applied to the completed table.
    assert [t.item(r, 1).text() for r in range(3)] == ["1", "2", "3"]


def test_copy_ignores_plain_items_without_numeric_data(qapp):
    """A table built with plain QTableWidgetItems still copies correctly."""
    t = QTableWidget(1, 2)
    t.setItem(0, 0, QTableWidgetItem("x"))
    t.setItem(0, 1, QTableWidgetItem("y"))
    assert copy_table_selection(t) == "x\ty"


# ── Consumers ────────────────────────────────────────────────────────────
#  The windows that motivated this module.  Constructing them offscreen is
#  cheap and catches a panel that stops calling the helper.

COLLECT_ROWS = [
    {"file": "/data/s_100C.h5", "x_value": 100.0, "y_value": 9.0, "y_error": 0.1},
    {"file": "/data/s_9C.h5", "x_value": 9.0, "y_value": 100.0, "y_error": None},
    {"file": "/data/s_10C.h5", "x_value": 10.0, "y_value": 10.0, "y_error": 0.2},
]


def test_collect_window_table_is_copyable_and_sortable(qapp):
    pytest.importorskip("pyqtgraph")
    from pyirena.gui.hdf5viewer.collect_window import CollectWindow

    win = CollectWindow(COLLECT_ROWS, "Temperature, C", "Rg, A", "Collected")
    tbl = win._table

    # Rows keep the file order the panel supplied…
    assert [tbl.item(r, 0).text() for r in range(3)] == [
        "s_100C.h5", "s_9C.h5", "s_10C.h5",
    ]
    # …until the user sorts, which must be numeric.
    tbl.sortItems(1, Qt.SortOrder.AscendingOrder)
    assert [tbl.item(r, 1).text() for r in range(3)] == ["9", "10", "100"]

    text = copy_table_selection(tbl, include_headers=True)
    assert text.startswith("File\tTemperature, C\tRg, A")
    assert len(text.split("\n")) == 4
    win.close()


def test_multi_collect_window_csv_quotes_labels_with_commas(qapp):
    from pyirena.gui.hdf5viewer.multi_collect_window import MultiCollectWindow

    rows = [{"file": "s1", "x_value": 1.0, "values": [1.0, "a,b"]}]
    win = MultiCollectWindow(rows, "T", ["Rg, A", "note"], "Multi")
    assert [win._table.item(0, c).text() for c in range(4)] == ["s1", "1", "1", "a,b"]

    csv_text = rows_to_csv_text(
        ["File", "T", "Rg, A", "note"],
        [[r["file"], r["x_value"], *r["values"]] for r in rows],
    )
    assert csv_text.splitlines()[0] == 'File,T,"Rg, A",note'
    assert csv_text.splitlines()[1] == 's1,1,1,"a,b"'
    win.close()


def test_tabulate_results_window_table_is_copyable(qapp):
    pytest.importorskip("pyqtgraph")
    from pyirena.gui.data_selector.results_windows import TabulateResultsWindow

    win = TabulateResultsWindow()
    win.set_data(["File", "Rg"], [["a.h5", 100.0], ["b.h5", 9.0]])
    win.table.sortItems(1, Qt.SortOrder.AscendingOrder)
    assert [win.table.item(r, 1).text() for r in range(2)] == ["9.0", "100.0"]
    assert copy_table_selection(win.table).splitlines()[0] == "b.h5\t9.0"
    win.close()
