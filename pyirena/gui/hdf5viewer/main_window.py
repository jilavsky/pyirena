"""
HDF5ViewerWindow — main window for the HDF5 Viewer / Data Extractor tool.

Three-panel horizontal splitter:
  left   — FileTreeWidget  (folder/subfolder/file browser)
  centre — HDF5BrowserWidget  (lazy HDF5 content tree)
  right  — PlotControlsPanel  (tabs: 1D Graph / Collect Values)

Manages floating GraphWindow and CollectWindow instances.
State is saved/loaded via StateManager.
"""

from __future__ import annotations

from pathlib import Path

try:
    from PySide6.QtWidgets import (
        QMainWindow, QWidget, QSplitter, QStatusBar, QHBoxLayout,
        QVBoxLayout, QLabel,
    )
    from PySide6.QtCore import Qt
    from PySide6.QtGui import QAction, QCloseEvent
except ImportError:
    from PyQt6.QtWidgets import (  # type: ignore[no-redef]
        QMainWindow, QWidget, QSplitter, QStatusBar, QHBoxLayout,
        QVBoxLayout, QLabel,
    )
    from PyQt6.QtCore import Qt  # type: ignore[no-redef]
    from PyQt6.QtGui import QAction, QCloseEvent  # type: ignore[no-redef]

from .file_tree import FileTreeWidget
from .hdf5_browser import HDF5BrowserWidget
from .plot_controls import PlotControlsPanel
from .graph_window import GraphWindow
from .collect_window import CollectWindow
from . import pyirena_readers as _readers


class HDF5ViewerWindow(QMainWindow):
    """
    Main window for the HDF5 Viewer / Data Extractor tool.

    Parameters
    ----------
    initial_folder : str, optional
        If given, the FileTreeWidget is pre-loaded with this folder.
    state_manager : StateManager, optional
        Shared state manager.  If None, a fresh one is created.
    parent : QWidget, optional
    """

    def __init__(
        self,
        initial_folder: str | None = None,
        state_manager=None,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("HDF5 Viewer / Data Extractor")
        self.resize(1150, 700)

        # State
        self._state_manager = state_manager
        if self._state_manager is None:
            try:
                from pyirena.state.state_manager import StateManager
                self._state_manager = StateManager()
            except Exception:
                self._state_manager = None

        # Active graph window for "Add to active graph"
        self._active_graph: GraphWindow | None = None

        # Keep references to all open windows so they aren't GC'd
        self._graph_windows: list[GraphWindow] = []
        self._collect_windows: list[CollectWindow] = []

        self._build_ui()
        self._wire_signals()
        self._restore_state()

        if initial_folder:
            self._file_tree.set_folder(initial_folder)

    # ── UI construction ────────────────────────────────────────────────────

    def _build_ui(self) -> None:
        central = QWidget()
        self.setCentralWidget(central)

        vl = QVBoxLayout(central)
        vl.setContentsMargins(0, 0, 0, 0)
        vl.setSpacing(0)

        # Application title
        title_lbl = QLabel("pyIrena HDF5 Viewer / Data Extractor")
        title_lbl.setStyleSheet(
            "font-size:13pt; font-weight:bold; color:#2c3e50;"
            "padding:6px 10px 5px 10px; background:#ecf0f1;"
            "border-bottom:1px solid #bdc3c7;"
        )
        vl.addWidget(title_lbl)

        splitter = QSplitter(Qt.Orientation.Horizontal)

        # Left panel — file tree
        self._file_tree = FileTreeWidget()
        self._file_tree.setMinimumWidth(220)
        splitter.addWidget(self._file_tree)

        # Centre panel — HDF5 content browser
        self._hdf5_browser = HDF5BrowserWidget()
        self._hdf5_browser.setMinimumWidth(260)
        splitter.addWidget(self._hdf5_browser)

        # Right panel — plot controls
        self._plot_controls = PlotControlsPanel()
        self._plot_controls.setMinimumWidth(280)
        splitter.addWidget(self._plot_controls)

        # Default widths
        splitter.setSizes([260, 320, 380])
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 0)
        splitter.setStretchFactor(2, 1)

        vl.addWidget(splitter, 1)

        # Status bar
        self._status_bar = QStatusBar()
        self.setStatusBar(self._status_bar)
        self._status_bar.showMessage("Ready")

        # Menu bar
        self._build_menu()

    def _build_menu(self) -> None:
        menu_bar = self.menuBar()

        file_menu = menu_bar.addMenu("&File")
        close_act = QAction("Close", self)
        close_act.triggered.connect(self.close)
        file_menu.addAction(close_act)

        graph_menu = menu_bar.addMenu("&Graph")
        new_graph_act = QAction("New Graph", self)
        new_graph_act.triggered.connect(
            lambda: self._plot_controls._on_new_graph()
        )
        graph_menu.addAction(new_graph_act)

    # ── Signal wiring ──────────────────────────────────────────────────────

    def _wire_signals(self) -> None:
        # File tree → HDF5 browser
        self._file_tree.selection_changed.connect(self._on_file_selection_changed)

        # HDF5 browser → plot controls
        self._hdf5_browser.add_dataset_requested.connect(
            self._plot_controls.set_dataset_role
        )
        self._hdf5_browser.plot_known_type_requested.connect(
            self._on_plot_known_type
        )
        self._hdf5_browser.collect_value_requested.connect(
            self._plot_controls.set_collect_custom_path
        )
        self._hdf5_browser.set_x_axis_path_requested.connect(
            self._plot_controls.set_x_axis_path
        )

        # Plot controls → main window
        self._plot_controls.new_graph_requested.connect(self._open_new_graph)
        self._plot_controls.add_to_active_graph_requested.connect(
            self._add_to_active_graph
        )
        self._plot_controls.collect_requested.connect(self._open_collect_window)
        self._plot_controls.status_message.connect(self._status_bar.showMessage)

    # ── File selection ─────────────────────────────────────────────────────

    def _on_file_selection_changed(self, paths: list[str]) -> None:
        # Show first selected file in HDF5 browser
        first = paths[0] if paths else None
        self._hdf5_browser.load_file(first)
        # Tell plot controls which files are selected
        self._plot_controls.set_selected_files(paths)
        n = len(paths)
        if n == 0:
            self._status_bar.showMessage("No files selected")
        elif n == 1:
            self._status_bar.showMessage(Path(paths[0]).name)
        else:
            self._status_bar.showMessage(f"{n} files selected")

    # ── GraphWindow management ─────────────────────────────────────────────

    def _open_new_graph(self, curves: list[dict]) -> None:
        """Create a new GraphWindow and populate it with *curves*."""
        gw = GraphWindow(parent=None)
        self._graph_windows.append(gw)
        gw.became_active.connect(lambda w=gw: self._set_active_graph(w))
        gw.closed.connect(lambda w=gw: self._on_graph_closed(w))
        gw.show()
        self._set_active_graph(gw)

        for c in curves:
            gw.add_curve(
                x=c["x"], y=c["y"],
                label=c.get("label", ""),
                yerr=c.get("yerr"),
                xerr=c.get("xerr"),
                suggest_log_x=c.get("suggest_log_x", False),
                suggest_log_y=c.get("suggest_log_y", False),
            )

    def _add_to_active_graph(self, curves: list[dict]) -> None:
        """Add *curves* to the currently active GraphWindow."""
        if self._active_graph is None or not self._active_graph.isVisible():
            # No active graph — open a new one
            self._open_new_graph(curves)
            return

        for c in curves:
            self._active_graph.add_curve(
                x=c["x"], y=c["y"],
                label=c.get("label", ""),
                yerr=c.get("yerr"),
                xerr=c.get("xerr"),
                suggest_log_x=c.get("suggest_log_x", False),
                suggest_log_y=c.get("suggest_log_y", False),
            )

    def _set_active_graph(self, gw: GraphWindow) -> None:
        self._active_graph = gw

    def _on_graph_closed(self, gw: GraphWindow) -> None:
        if gw in self._graph_windows:
            self._graph_windows.remove(gw)
        if self._active_graph is gw:
            self._active_graph = (
                self._graph_windows[-1] if self._graph_windows else None
            )

    # ── Known-type plotting (from HDF5 browser right-click) ───────────────

    def _on_plot_known_type(self, type_key: str, group_path: str) -> None:
        """
        Immediately read data for *type_key* from the currently browsed file
        and open a new GraphWindow.
        """
        filepath = self._hdf5_browser._filepath
        if not filepath:
            self._status_bar.showMessage("No file loaded in HDF5 browser.")
            return

        stem = Path(filepath).stem
        curves = []

        try:
            if type_key == "nxcansas":
                result = _readers.read_nxcansas(filepath)
                if result:
                    curves.append({
                        "label": stem,
                        "x": result["Q"], "y": result["I"],
                        "yerr": result.get("dI"),
                        "xerr": None,
                        "suggest_log_x": True,
                        "suggest_log_y": True,
                    })

            elif type_key == "unified_fit":
                result = _readers.read_unified_fit(filepath)
                if result:
                    # Data + model
                    if result.get("I_data") is not None:
                        curves.append({
                            "label": f"{stem}  data",
                            "x": result["Q"], "y": result["I_data"],
                            "yerr": result.get("dI_data"),
                            "xerr": None,
                            "suggest_log_x": True,
                            "suggest_log_y": True,
                        })
                    curves.append({
                        "label": f"{stem}  UF model",
                        "x": result["Q"], "y": result["I_model"],
                        "yerr": None, "xerr": None,
                        "suggest_log_x": True,
                        "suggest_log_y": True,
                    })

            elif type_key == "sizes":
                result = _readers.read_sizes(filepath)
                if result:
                    curves.append({
                        "label": f"{stem}  Sizes I(Q)",
                        "x": result["Q"], "y": result["I_model"],
                        "yerr": None, "xerr": None,
                        "suggest_log_x": True,
                        "suggest_log_y": True,
                    })

            elif type_key == "waxs":
                result = _readers.read_waxs(filepath)
                if result and len(result["Q"]) > 0:
                    curves.append({
                        "label": f"{stem}  WAXS fit",
                        "x": result["Q"], "y": result["I_fit"],
                        "yerr": None, "xerr": None,
                        "suggest_log_x": False,
                        "suggest_log_y": False,
                    })

            elif type_key == "simple_fit":
                result = _readers.read_simple_fit(filepath)
                if result:
                    curves.append({
                        "label": f"{stem}  {result.get('model_name', 'Simple Fit')}",
                        "x": result["Q"], "y": result["I_model"],
                        "yerr": None, "xerr": None,
                        "suggest_log_x": True,
                        "suggest_log_y": True,
                    })

        except Exception as exc:
            self._status_bar.showMessage(f"Error reading {type_key}: {exc}")
            return

        if curves:
            self._open_new_graph(curves)
            self._status_bar.showMessage(
                f"Plotted {len(curves)} curve(s) from {stem}"
            )
        else:
            self._status_bar.showMessage(
                f"No {type_key} data found in {stem}"
            )

    # ── Collect window ─────────────────────────────────────────────────────

    def _open_collect_window(
        self, spec: dict, x_spec: dict, label_spec: dict
    ) -> None:
        """
        Collect 0D values from all selected files and open a CollectWindow.
        """
        selected = self._file_tree.get_selected_paths()
        if not selected:
            self._status_bar.showMessage("No files selected.")
            return

        rows = []
        for i, filepath in enumerate(selected):
            stem = Path(filepath).stem

            # Collect Y value
            y_val = None
            y_err = None
            try:
                collected = _readers.collect_value(filepath, spec)
                if collected is not None:
                    if isinstance(collected, (list, tuple)) and len(collected) == 2:
                        y_val, y_err = float(collected[0]), float(collected[1])
                    else:
                        y_val = float(collected)
            except Exception:
                pass

            # Collect X value
            x_val: float | None = None
            try:
                if x_spec.get("type") == "order":
                    x_val = float(i + 1)
                elif x_spec.get("type") == "sortkey":
                    x_val = _readers.extract_sort_key_value(
                        filepath, x_spec.get("sort_index", 6)
                    )
                    if x_val is None:
                        x_val = float(i + 1)
                elif x_spec.get("type") == "path":
                    x_val = _readers.read_metadata_value(
                        filepath, x_spec.get("path", "")
                    )
                    if x_val is None:
                        x_val = float(i + 1)
            except Exception:
                x_val = float(i + 1)

            rows.append({
                "file": stem,
                "x_value": x_val if x_val is not None else float(i + 1),
                "y_value": y_val,
                "y_error": y_err,
            })

        if not any(r["y_value"] is not None for r in rows):
            self._status_bar.showMessage("No values could be collected.")
            return

        cw = CollectWindow(
            rows=rows,
            x_label=label_spec.get("x_label", ""),
            y_label=label_spec.get("y_label", ""),
            title=label_spec.get("title", ""),
            parent=None,
        )
        self._collect_windows.append(cw)
        cw.destroyed.connect(lambda obj=None, w=cw: self._on_collect_closed(w))
        cw.show()
        self._status_bar.showMessage(
            f"Collected {sum(r['y_value'] is not None for r in rows)}/{len(rows)} values."
        )

    def _on_collect_closed(self, cw: CollectWindow) -> None:
        if cw in self._collect_windows:
            self._collect_windows.remove(cw)

    # ── State save / restore ───────────────────────────────────────────────

    def _restore_state(self) -> None:
        if self._state_manager is None:
            return
        try:
            st = self._state_manager.get("hdf5_viewer") or {}
            folder = st.get("last_folder", "")
            if folder:
                self._file_tree.set_folder(folder)
            sort_idx = st.get("sort_index", 6)
            self._file_tree.set_sort_index(sort_idx)
            flt = st.get("filter_text", "")
            if flt:
                self._file_tree.set_filter(flt)
        except Exception:
            pass

    def _save_state(self) -> None:
        if self._state_manager is None:
            return
        try:
            st = {
                "last_folder": self._file_tree.root_folder,
                "sort_index": self._file_tree._sort_index,
                "filter_text": self._file_tree._filter_text,
            }
            self._state_manager.update("hdf5_viewer", st)
            self._state_manager.save()
        except Exception:
            pass

    # ── Window lifecycle ───────────────────────────────────────────────────

    def closeEvent(self, event: QCloseEvent) -> None:
        self._save_state()
        # Close child windows
        for gw in list(self._graph_windows):
            gw.close()
        for cw in list(self._collect_windows):
            cw.close()
        self._hdf5_browser.clear()
        super().closeEvent(event)
