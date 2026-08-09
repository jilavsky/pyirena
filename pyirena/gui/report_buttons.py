"""
"Copy results" / "Save report…" buttons for the fit panels.

Igor Irena wrote every fit to the notebook, so users arrive expecting a text
trail they can paste into an e-mail, a logbook or a paper draft.  In pyIrena
the only ways out of a fit panel used to be save-to-HDF5 (then tabulate in the
Data Selector) or reading numbers off the widgets one at a time — while the
api/control layer had ``export_fit_report`` for AI agents and not for people.

This module is the GUI half of the fix.  The text itself comes from
:mod:`pyirena.core.reporting`, the same builder the Data Selector's *Create
Report* and the MCP ``export_fit_report`` use, so a value can never be
formatted one way in the panel and another way in the report.

A panel wires both buttons in one call::

    from pyirena.gui.report_buttons import make_report_buttons

    row = make_report_buttons(self, self.results_for_report,
                              tool_key='fit_results',
                              default_stem='unified_fit')
    layout.addLayout(row)

where ``results_for_report()`` returns the dict shape
``pyirena.io.load_<tool>_results()`` produces.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Callable, Optional

from pyirena.core.reporting import build_report
from pyirena.gui._qt import (
    QApplication,
    QFileDialog,
    QHBoxLayout,
    QMessageBox,
    QPushButton,
)

log = logging.getLogger(__name__)

__all__ = ["build_panel_report", "make_report_buttons", "copy_report", "save_report"]

_COPY_STYLE = (
    "QPushButton{background:#16a085;color:white;font-weight:bold;"
    "border-radius:4px;padding:4px 10px;}"
    "QPushButton:hover{background:#138d75;}"
    "QPushButton:disabled{background:#bdc3c7;}"
)
_SAVE_STYLE = (
    "QPushButton{background:#2980b9;color:white;font-weight:bold;"
    "border-radius:4px;padding:4px 10px;}"
    "QPushButton:hover{background:#2471a3;}"
    "QPushButton:disabled{background:#bdc3c7;}"
)

COPY_TOOLTIP = (
    "Copy the current results as Markdown text — parameters, fit quality and\n"
    "setup — ready to paste into an e-mail, logbook or manuscript.\n"
    "Identical to what Create Report writes for the saved file."
)
SAVE_TOOLTIP = (
    "Save the current results as a Markdown (.md) report file."
)


def build_panel_report(
    results: Optional[dict],
    tool_key: str,
    file_path: str = "",
    data_info: Optional[dict] = None,
) -> str:
    """Render one tool's results as a Markdown report.

    Args:
        results: Results in ``load_<tool>_results()`` shape, or None.
        tool_key: Which section of :func:`pyirena.core.reporting.build_report`
            the results belong to — ``fit_results`` (Unified), ``sizes_results``,
            ``simple_fit_results``, ``waxs_peakfit_results``,
            ``modeling_results``.
        file_path: Source file, for the report header.
        data_info: Optional ``{'Q','I','I_error'}`` for the data summary.

    Returns:
        Markdown text, or an empty string when there is nothing to report.
    """
    if not results:
        return ""
    sections = {tool_key: results}
    if data_info is not None:
        sections["data_info"] = data_info
    return build_report(file_path or "results", **sections)


def copy_report(parent, text: str, status_setter: Optional[Callable] = None) -> bool:
    """Put report *text* on the clipboard, telling the user what happened."""
    if not text:
        QMessageBox.information(
            parent, "Nothing to report",
            "There are no results yet — run a fit first.",
        )
        return False
    clipboard = QApplication.clipboard()
    if clipboard is None:
        return False
    clipboard.setText(text)
    n_lines = len(text.splitlines())
    if status_setter is not None:
        try:
            status_setter(f"Results copied to clipboard ({n_lines} lines of Markdown).")
        except Exception:
            log.debug("suppressed exception setting status", exc_info=True)
    return True


def save_report(parent, text: str, default_stem: str = "pyirena_report",
                folder=None, status_setter: Optional[Callable] = None) -> Optional[str]:
    """Write report *text* to a ``.md`` file chosen by the user."""
    if not text:
        QMessageBox.information(
            parent, "Nothing to report",
            "There are no results yet — run a fit first.",
        )
        return None

    # Reuse the shared export-folder memory so the report lands where the last
    # graph or table did.
    from pyirena.gui.plot_export import export_folder, remember_export_folder

    default = str(Path(export_folder(folder)) / f"{default_stem}_report.md")
    path, _ = QFileDialog.getSaveFileName(
        parent, "Save results report",
        default,
        "Markdown files (*.md);;Text files (*.txt);;All files (*)",
    )
    if not path:
        return None
    if Path(path).suffix.lower() not in (".md", ".txt"):
        path += ".md"

    try:
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(text)
    except OSError as exc:
        QMessageBox.warning(parent, "Save failed", f"Could not save report:\n{exc}")
        return None

    remember_export_folder(path)
    if status_setter is not None:
        try:
            status_setter(f"Report saved to {Path(path).name}")
        except Exception:
            log.debug("suppressed exception setting status", exc_info=True)
    return path


def make_report_buttons(
    parent,
    results_provider: Callable[[], Optional[dict]],
    tool_key: str,
    default_stem: str = "pyirena",
    file_path_provider: Optional[Callable[[], str]] = None,
    data_info_provider: Optional[Callable[[], Optional[dict]]] = None,
    status_setter: Optional[Callable] = None,
    folder_provider: Optional[Callable] = None,
) -> QHBoxLayout:
    """Build the standard "Copy results" / "Save report…" button row.

    Args:
        parent: Panel the buttons belong to (dialog parent).
        results_provider: Callable returning the results dict, evaluated on
            click so the report always reflects the current fit.
        tool_key: Report section name — see :func:`build_panel_report`.
        default_stem: File stem offered in the save dialog.
        file_path_provider: Callable returning the source file path, for the
            report header.
        data_info_provider: Callable returning ``{'Q','I','I_error'}`` so the
            report can include the data summary.
        status_setter: Callable taking a message, to confirm in the panel's
            status line.
        folder_provider: Callable returning the folder the save dialog should
            prefer (usually the data folder).

    Returns:
        A ``QHBoxLayout`` holding the two buttons.
    """
    def _text() -> str:
        try:
            results = results_provider()
        except Exception as exc:
            log.debug("report provider failed", exc_info=True)
            QMessageBox.warning(parent, "Report failed",
                                f"Could not collect results:\n{exc}")
            return ""
        file_path = ""
        if file_path_provider is not None:
            try:
                file_path = file_path_provider() or ""
            except Exception:
                log.debug("suppressed exception reading file path", exc_info=True)
        data_info = None
        if data_info_provider is not None:
            try:
                data_info = data_info_provider()
            except Exception:
                log.debug("suppressed exception reading data info", exc_info=True)
        return build_panel_report(results, tool_key, file_path, data_info)

    copy_btn = QPushButton("Copy results")
    copy_btn.setStyleSheet(_COPY_STYLE)
    copy_btn.setToolTip(COPY_TOOLTIP)
    copy_btn.clicked.connect(lambda: copy_report(parent, _text(), status_setter))

    save_btn = QPushButton("Save report…")
    save_btn.setStyleSheet(_SAVE_STYLE)
    save_btn.setToolTip(SAVE_TOOLTIP)
    save_btn.clicked.connect(
        lambda: save_report(
            parent, _text(), default_stem,
            folder_provider() if folder_provider is not None else None,
            status_setter,
        )
    )

    row = QHBoxLayout()
    row.setSpacing(6)
    row.addWidget(copy_btn)
    row.addWidget(save_btn)
    row.addStretch()
    # Kept reachable for tests and for panels that want to enable/disable them.
    row.copy_button = copy_btn
    row.save_button = save_btn
    row.report_text = _text
    return row
