"""
export_to_igor_tab.py — "Export to Igor" tab for the Data Explorer.

Lets users select one or more pyirena HDF5 files from the viewer file list,
choose which analysis tools to export, specify an output .h5xp path, and run
the batch export in a background thread.

Used by PlotControlsPanel which owns the QTabWidget.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

try:
    from PySide6.QtWidgets import (
        QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
        QPushButton, QLabel, QLineEdit, QCheckBox, QRadioButton,
        QGroupBox, QPlainTextEdit, QProgressBar, QFileDialog,
        QSizePolicy,
    )
    from PySide6.QtCore import Qt, QThread, Signal
    from PySide6.QtGui import QFont
except ImportError:
    from PyQt6.QtWidgets import (  # type: ignore[no-redef]
        QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
        QPushButton, QLabel, QLineEdit, QCheckBox, QRadioButton,
        QGroupBox, QPlainTextEdit, QProgressBar, QFileDialog,
        QSizePolicy,
    )
    from PyQt6.QtCore import Qt, QThread, pyqtSignal as Signal  # type: ignore[no-redef]
    from PyQt6.QtGui import QFont  # type: ignore[no-redef]


# ---------------------------------------------------------------------------
# Style constants
# ---------------------------------------------------------------------------

_BTN_GREEN = ("QPushButton { background:#27ae60; color:white; font-weight:bold; "
              "border-radius:4px; padding:6px 10px; }"
              "QPushButton:hover { background:#2ecc71; }"
              "QPushButton:disabled { background:#95a5a6; }")

_BTN_GREY  = ("QPushButton { background:#7f8c8d; color:white; "
              "border-radius:4px; padding:4px 8px; }"
              "QPushButton:hover { background:#95a5a6; }")

_TOOL_ROWS = [
    ("unified_fit",       "Unified Fit"),
    ("size_distribution", "Size Distribution"),
    ("waxs_peakfit",      "WAXS Peak Fit"),
    ("simple_fits",       "Simple Fits"),
    ("modeling",          "Modeling"),
    ("saxs_morph",        "SAXS Morph"),
    ("fractals",          "Fractals"),
]


# ---------------------------------------------------------------------------
# Background worker thread
# ---------------------------------------------------------------------------

class _ExportWorker(QThread):
    progress = Signal(str)
    finished = Signal(bool)

    def __init__(self, sources: list[Path], h5xp_path: Path,
                 tools: list[str] | None, category: str | None,
                 overwrite: bool, build_results_table: bool):
        super().__init__()
        self._sources = sources
        self._h5xp    = h5xp_path
        self._tools   = tools
        self._category = category
        self._overwrite = overwrite
        self._build_rt  = build_results_table

    def run(self):
        try:
            from pyirena.io.h5xp_extractor import batch_extract_to_h5xp
            self.progress.emit(f"Writing → {self._h5xp.name}\n")
            summaries = batch_extract_to_h5xp(
                self._sources, self._h5xp,
                tools=self._tools,
                category=self._category,
                overwrite=self._overwrite,
                build_results_table=self._build_rt,
            )
            for s in summaries:
                folder  = s.get("folder", "?")
                ok      = s.get("tools_written", [])
                skipped = s.get("tools_skipped", [])
                absent  = s.get("tools_absent", [])
                iq_ok   = "IQ ok" if s.get("iq_written") else "IQ MISSING"
                lines   = [f"  {Path(folder).name}: {iq_ok}"]
                if ok:
                    lines.append(f"    wrote: {', '.join(ok)}")
                if skipped:
                    lines.append(f"    skipped (no data): {', '.join(skipped)}")
                if absent:
                    lines.append(f"    absent: {', '.join(absent)}")
                self.progress.emit("\n".join(lines) + "\n")
            self.progress.emit(
                f"\nDone.  {len(summaries)} file(s) → {self._h5xp}\n"
            )
            self.finished.emit(True)
        except Exception as exc:
            import traceback
            self.progress.emit(f"\nERROR: {exc}\n{traceback.format_exc()}\n")
            self.finished.emit(False)


# ---------------------------------------------------------------------------
# Export to Igor tab widget
# ---------------------------------------------------------------------------

class ExportToIgorTab(QWidget):
    """
    Self-contained 'Export to Igor' panel.

    Call :meth:`set_selected_files` whenever the file selection changes in the
    parent FileTreeWidget so the export button knows what to process.
    """

    status_message = Signal(str)

    def __init__(self, parent: QWidget | None = None):
        super().__init__(parent)
        self._selected_files: list[Path] = []
        self._worker: Optional[_ExportWorker] = None
        self._setup_ui()

    # ── Public API ───────────────────────────────────────────────────────────

    def set_selected_files(self, paths: list[str]) -> None:
        """Called by PlotControlsPanel whenever the file selection changes."""
        self._selected_files = [Path(p) for p in paths]
        n = len(paths)
        self._sel_lbl.setText(
            f"{n} file(s) selected" if n else "No files selected"
        )

    # ── UI construction ──────────────────────────────────────────────────────

    def _setup_ui(self):
        lay = QVBoxLayout(self)
        lay.setContentsMargins(6, 6, 6, 6)
        lay.setSpacing(6)

        # Selection indicator
        self._sel_lbl = QLabel("No files selected")
        self._sel_lbl.setStyleSheet("color:#555; font-style:italic;")
        lay.addWidget(self._sel_lbl)

        # ── Tools ────────────────────────────────────────────────────────────
        grp_tools = QGroupBox("Tools to export")
        tools_lay = QVBoxLayout(grp_tools)
        tools_lay.setSpacing(2)

        self._tool_checks: dict[str, QCheckBox] = {}
        for key, label in _TOOL_ROWS:
            cb = QCheckBox(label)
            cb.setChecked(True)
            tools_lay.addWidget(cb)
            self._tool_checks[key] = cb

        row_sel = QHBoxLayout()
        btn_all = QPushButton("All")
        btn_all.setStyleSheet(_BTN_GREY)
        btn_all.setFixedWidth(50)
        btn_all.clicked.connect(lambda: self._set_all_tools(True))
        btn_none = QPushButton("None")
        btn_none.setStyleSheet(_BTN_GREY)
        btn_none.setFixedWidth(50)
        btn_none.clicked.connect(lambda: self._set_all_tools(False))
        row_sel.addWidget(btn_all)
        row_sel.addWidget(btn_none)
        row_sel.addStretch()
        tools_lay.addLayout(row_sel)
        lay.addWidget(grp_tools)

        # ── Scattering geometry ───────────────────────────────────────────────
        grp_cat = QGroupBox("Scattering geometry")
        cat_lay = QHBoxLayout(grp_cat)
        self._cat_auto = QRadioButton("Auto-detect")
        self._cat_saxs = QRadioButton("SAXS")
        self._cat_waxs = QRadioButton("WAXS")
        self._cat_auto.setChecked(True)
        for rb in (self._cat_auto, self._cat_saxs, self._cat_waxs):
            cat_lay.addWidget(rb)
        cat_lay.addStretch()
        lay.addWidget(grp_cat)

        # ── Output file ───────────────────────────────────────────────────────
        grp_out = QGroupBox("Output h5xp file")
        out_lay = QGridLayout(grp_out)
        self._out_edit = QLineEdit()
        self._out_edit.setPlaceholderText("output.h5xp")
        out_lay.addWidget(self._out_edit, 0, 0)
        btn_browse = QPushButton("Browse…")
        btn_browse.clicked.connect(self._browse_output)
        out_lay.addWidget(btn_browse, 0, 1)

        self._overwrite_cb = QCheckBox("Overwrite if file exists")
        self._overwrite_cb.setChecked(True)
        out_lay.addWidget(self._overwrite_cb, 1, 0, 1, 2)

        self._results_cb = QCheckBox(
            "Build Results table  (scalar parameters vs. sample)"
        )
        self._results_cb.setChecked(True)
        out_lay.addWidget(self._results_cb, 2, 0, 1, 2)
        lay.addWidget(grp_out)

        # ── Export button ─────────────────────────────────────────────────────
        self._btn_export = QPushButton("Export selected files to Igor h5xp")
        self._btn_export.setStyleSheet(_BTN_GREEN)
        self._btn_export.clicked.connect(self._run_export)
        lay.addWidget(self._btn_export)

        # ── Progress bar ──────────────────────────────────────────────────────
        self._progress = QProgressBar()
        self._progress.setRange(0, 0)
        self._progress.setVisible(False)
        lay.addWidget(self._progress)

        # ── Log ───────────────────────────────────────────────────────────────
        lay.addWidget(QLabel("Export log:"))
        self._log = QPlainTextEdit()
        self._log.setReadOnly(True)
        self._log.setMaximumBlockCount(2000)
        self._log.setFont(QFont("Courier New", 9))
        lay.addWidget(self._log, 1)

        # ── Open folder ───────────────────────────────────────────────────────
        self._btn_open = QPushButton("Open output folder in Explorer")
        self._btn_open.setStyleSheet(_BTN_GREY)
        self._btn_open.setEnabled(False)
        self._btn_open.clicked.connect(self._open_output_folder)
        lay.addWidget(self._btn_open)

    # ── Helpers ──────────────────────────────────────────────────────────────

    def _set_all_tools(self, state: bool) -> None:
        for cb in self._tool_checks.values():
            cb.setChecked(state)

    def _browse_output(self) -> None:
        # Default to parent folder of first selected file, or CWD
        start = ""
        if self._selected_files:
            start = str(self._selected_files[0].parent)
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Igor h5xp", start,
            "Igor packed experiment (*.h5xp)"
        )
        if path:
            if not path.lower().endswith(".h5xp"):
                path += ".h5xp"
            self._out_edit.setText(path)

    def _open_output_folder(self) -> None:
        out = self._out_edit.text().strip()
        if out:
            try:
                from PySide6.QtGui import QDesktopServices
                from PySide6.QtCore import QUrl
            except ImportError:
                from PyQt6.QtGui import QDesktopServices  # type: ignore
                from PyQt6.QtCore import QUrl  # type: ignore
            QDesktopServices.openUrl(
                QUrl.fromLocalFile(str(Path(out).parent))
            )

    def _selected_tools(self) -> list[str]:
        return [k for k, cb in self._tool_checks.items() if cb.isChecked()]

    def _category(self) -> str | None:
        if self._cat_saxs.isChecked():
            return "SAXS"
        if self._cat_waxs.isChecked():
            return "WAXS"
        return None

    # ── Export ────────────────────────────────────────────────────────────────

    def _run_export(self) -> None:
        if not self._selected_files:
            self._log.appendPlainText(
                "No files selected — use the File Select panel on the left.\n"
            )
            return

        out_str = self._out_edit.text().strip()
        if not out_str:
            self._log.appendPlainText("No output file specified.\n")
            return

        out_path = Path(out_str)

        self._log.clear()
        self._log.appendPlainText(
            f"Exporting {len(self._selected_files)} file(s)…\n"
            f"Output: {out_path}\n"
        )
        self._btn_export.setEnabled(False)
        self._btn_open.setEnabled(False)
        self._progress.setVisible(True)
        self.status_message.emit(
            f"Exporting {len(self._selected_files)} file(s) to h5xp…"
        )

        self._worker = _ExportWorker(
            sources=self._selected_files,
            h5xp_path=out_path,
            tools=self._selected_tools() or None,
            category=self._category(),
            overwrite=self._overwrite_cb.isChecked(),
            build_results_table=self._results_cb.isChecked(),
        )
        self._worker.progress.connect(self._log.appendPlainText)
        self._worker.finished.connect(self._on_finished)
        self._worker.start()

    def _on_finished(self, ok: bool) -> None:
        self._progress.setVisible(False)
        self._btn_export.setEnabled(True)
        self._btn_open.setEnabled(ok)
        msg = "Export complete." if ok else "Export FAILED — see log above."
        self._log.appendPlainText(msg)
        self.status_message.emit(msg)
