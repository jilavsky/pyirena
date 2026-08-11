"""
pyirena.gui.data_selector.report — Markdown fit-report builder.

The builder itself moved to :mod:`pyirena.core.reporting` so that the fit
panels' "Copy results" / "Save report…" buttons and the api/control layer's
``export_fit_report`` (used by AI agents over MCP) can share it — ``api`` and
``core`` may not import from ``gui``.  The generated text is unchanged.

This module re-exports the names the Data Selector already imported.
"""

from pyirena.core.reporting import _build_report, _quality_report_rows, build_report

__all__ = ["_build_report", "build_report", "_quality_report_rows"]
