"""
pyirena.gui.data_loading — shared data-loading helpers and the DataFileLoaderRow widget.

This module provides:

``load_data_file(parent, file_path, error_fraction)``
    The canonical "load one SAS file" function used by every tool panel.
    Text files are routed through ``ensure_nxcansas_sibling`` (clean + convert)
    so all callers receive a valid HDF5 path and ``is_nxcansas=True``.
    HDF5 files with multiple datasets show a picker dialog.
    Returns ``(data_dict, hdf5_path_str, display_name)`` or ``None``.

``DataFileLoaderRow``
    A thin reusable QWidget — one read-only filename field + "Open…" button.
    Emits ``data_loaded(data, hdf5_path, display_name)`` after a successful
    file open via its own dialog.  Also keeps the filename display in sync
    when data is pushed in externally (via the Data Selector).

All tool panels use ``DataFileLoaderRow`` at the top of their left panel so
a file can be opened either from the Data Selector or directly from the tool.
"""

from __future__ import annotations
import logging

log = logging.getLogger(__name__)


from pathlib import Path
from typing import Callable, Optional

from pyirena.gui._qt import (
    QHBoxLayout, QInputDialog, QLineEdit, QMessageBox, QPushButton, QWidget, Signal,
)


# ── Standalone helpers (no Qt parent required) ────────────────────────────────

def prompt_dataset_choice(parent: QWidget, filename: str, datasets: list) -> Optional[str]:
    """Ask the user which SAS dataset to load from a multi-dataset HDF5 file.

    Parameters
    ----------
    parent    : QWidget to use as dialog parent.
    filename  : File name for the dialog title.
    datasets  : List of ``{'path', 'name', ...}`` dicts from
                ``list_nxcansas_datasets``.

    Returns
    -------
    HDF5 data_path string, or ``None`` if the user cancelled.
    """
    items = [f"{d['name']}   [{d['path']}]" for d in datasets]
    item, ok = QInputDialog.getItem(
        parent,
        "Multiple data sets",
        f"'{filename}' contains {len(datasets)} SAS data sets.\n"
        "Select the one to load:",
        items,
        0,      # default — first (the file's @default dataset)
        False,  # not editable
    )
    if not ok:
        return None
    return datasets[items.index(item)]['path']


def read_nxcansas_with_picker(
    parent: QWidget,
    path: str,
    filename: str,
    status_cb: Optional[Callable[[str], None]] = None,
) -> Optional[dict]:
    """Load NXcanSAS data, prompting when the file holds several datasets.

    Parameters
    ----------
    parent    : QWidget to use as dialog parent.
    path      : Directory containing the file.
    filename  : File name.
    status_cb : Optional callable to receive status-text updates.

    Returns
    -------
    Data dict on success, or ``None`` on failure / user cancel.
    """
    from pyirena.io.hdf5 import readGenericNXcanSAS, list_nxcansas_datasets, _filter_smr

    try:
        datasets = list_nxcansas_datasets(path, filename)
    except Exception as e:
        QMessageBox.critical(parent, "Load Error", f"Could not read {filename}:\n{e}")
        return None

    selectable = _filter_smr(datasets)
    data_path = None
    if len(selectable) > 1:
        data_path = prompt_dataset_choice(parent, filename, selectable)
        if data_path is None:
            if status_cb:
                status_cb("Load cancelled — no data set selected.")
            return None
    elif selectable:
        data_path = selectable[0]['path']

    data = readGenericNXcanSAS(path, filename, data_path=data_path)
    if data is None:
        QMessageBox.critical(parent, "Load Error", f"Could not load data from {filename}")
    return data


def load_data_file(
    parent: QWidget,
    file_path: str,
    error_fraction: float = 0.05,
) -> Optional[tuple]:
    """Load one SAS file for use by a fitting tool.

    Text files (.txt/.dat) are converted to a cleaned NXcanSAS HDF5 sibling
    on first use (see ``pyirena.io.text_import``).  All callers receive
    ``is_nxcansas=True`` and a valid HDF5 filepath.

    Parameters
    ----------
    parent         : QWidget to use as dialog parent for error messages.
    file_path      : Path to the data file (text or HDF5).
    error_fraction : Fractional uncertainty to synthesize when absent.

    Returns
    -------
    ``(data_dict, hdf5_path_str, display_name)`` on success, or ``None``.

    ``data_dict``     has keys Q, Intensity, Error, … (from readGenericNXcanSAS).
    ``hdf5_path_str`` is the path tools should use for result saving.
    ``display_name``  is a human-readable label (original filename/stem).
    """
    from pyirena.io.text_import import ensure_nxcansas_sibling

    fp = Path(file_path)
    try:
        if fp.suffix.lower() in ('.txt', '.dat'):
            h5_path = ensure_nxcansas_sibling(fp, error_fraction=error_fraction)
            data = read_nxcansas_with_picker(
                parent, str(h5_path.parent), h5_path.name
            )
            if data is None:
                return None
            return data, str(h5_path), fp.stem + h5_path.suffix
        else:
            data = read_nxcansas_with_picker(parent, str(fp.parent), fp.name)
            if data is None:
                return None
            return data, str(fp), fp.name
    except Exception as exc:
        QMessageBox.critical(parent, "Load Error", f"Could not load '{fp.name}':\n{exc}")
        return None


# ── DataFileLoaderRow widget ──────────────────────────────────────────────────

class DataFileLoaderRow(QWidget):
    """A one-row widget: read-only filename field + 'Open…' button.

    Emits ``data_loaded(data, hdf5_path, display_name)`` after the user picks
    a file through the dialog and it loads successfully.

    Also call ``set_filename(name)`` when data arrives via ``set_data`` from
    the Data Selector so the display stays in sync.

    Usage::

        self.data_loader = DataFileLoaderRow(state_manager=self.state_manager)
        self.data_loader.data_loaded.connect(self._on_data_loaded)
        layout.addWidget(self.data_loader)

        # In set_data():
        self.data_loader.set_filename(display_name)

        # Slot:
        def _on_data_loaded(self, data, hdf5_path, display_name):
            self.set_data(data['Q'], data['Intensity'], data.get('Error'),
                          display_name, filepath=hdf5_path, is_nxcansas=True)
    """

    data_loaded = Signal(object, str, str)   # data_dict, hdf5_path, display_name

    def __init__(self, state_manager=None, parent: Optional[QWidget] = None):
        super().__init__(parent)
        self._state_manager = state_manager
        self._error_fraction: float = 0.05

        row = QHBoxLayout(self)
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(4)

        self._edit = QLineEdit()
        self._edit.setReadOnly(True)
        self._edit.setPlaceholderText("(no file selected)")
        self._edit.setToolTip(
            "Currently loaded data file.\n"
            "Click 'Open…' to load a different file, or select one in the Data Selector."
        )
        row.addWidget(self._edit, 1)

        btn = QPushButton("Open…")
        btn.setFixedWidth(60)
        btn.setToolTip(
            "Open an NXcanSAS/HDF5 or ASCII text (.dat/.txt) file.\n"
            "Text files are automatically cleaned and converted to NXcanSAS."
        )
        btn.clicked.connect(self._on_open_clicked)
        row.addWidget(btn)

    # ── Public API ────────────────────────────────────────────────────────

    def set_filename(self, name: str) -> None:
        """Update the display field (call from set_data when data arrives externally)."""
        self._edit.setText(name)

    def set_error_fraction(self, value: float) -> None:
        """Override the default 5 % fractional uncertainty used for text files."""
        self._error_fraction = value

    # ── Internal ──────────────────────────────────────────────────────────

    def _last_folder(self) -> str:
        if self._state_manager is not None:
            return (self._state_manager.get('data_selector') or {}).get('last_folder', '')
        return ''

    def _save_last_folder(self, folder: str) -> None:
        if self._state_manager is not None:
            try:
                self._state_manager.set('data_selector', 'last_folder', folder)
            except Exception:
                log.debug("suppressed exception", exc_info=True)

    def _on_open_clicked(self) -> None:
        try:
            from PySide6.QtWidgets import QFileDialog
        except ImportError:
            from PyQt6.QtWidgets import QFileDialog

        path, _ = QFileDialog.getOpenFileName(
            self,
            "Open SAS data file",
            self._last_folder(),
            "SAS data (*.h5 *.hdf5 *.hdf *.nxs *.dat *.txt);;All files (*)",
        )
        if not path:
            return

        self._save_last_folder(str(Path(path).parent))

        # Sync error_fraction from state manager if available
        if self._state_manager is not None:
            self._error_fraction = self._state_manager.get(
                'data_selector', 'error_fraction', 0.05
            )

        res = load_data_file(self, path, error_fraction=self._error_fraction)
        if res is None:
            return
        data, hdf5_path, display_name = res
        self._edit.setText(display_name)
        self.data_loaded.emit(data, hdf5_path, display_name)
