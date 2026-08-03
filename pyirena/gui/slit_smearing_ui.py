"""
Shared "Slit smeared data" GUI row + wiring for pyIrena fitting panels.

Every fitting tool that can fit slit-smeared USAXS data shows the same compact
control: a "Slit smeared data" checkbox, an editable slit-length field, and a
status label.  The behaviour (SasView-style) is uniform:

* Slit-length (``dQl``) presence in the loaded NXcanSAS data drives smearing
  automatically; the checkbox reflects/toggles it.
* When a file carries BOTH a desmeared and a slit-smeared copy (Matilda),
  toggling the checkbox reloads the appropriate dataset from file.
* The slit length is file-derived but user-editable (rarely changed).

Panels mix in :class:`SlitSmearingMixin` and:

1. call ``self._build_slit_row(layout)`` while constructing controls;
2. call ``self._refresh_slit_ui_from_data(slit_length, is_slit_smeared,
   filepath)`` from ``set_data``;
3. read ``self.slit_active()`` / ``self.current_slit_length()`` when building
   the fit/model object;
4. optionally implement ``_reload_data_with_smearing(prefer)`` (to reload the
   smeared/desmeared sibling) and ``_replot_after_slit_change()`` (to redraw).

Qt import indirection keeps the module usable under both PySide6 and PyQt6.
"""

from __future__ import annotations

from pathlib import Path

try:  # PySide6 first, fall back to PyQt6 (project-wide convention)
    from PySide6.QtWidgets import (
        QCheckBox, QHBoxLayout, QLabel, QLineEdit, QVBoxLayout, QWidget,
    )
except ImportError:  # pragma: no cover - exercised only on PyQt6 installs
    from PyQt6.QtWidgets import (
        QCheckBox, QHBoxLayout, QLabel, QLineEdit, QVBoxLayout, QWidget,
    )


class SlitSmearingMixin:
    """Provides the shared "Slit smeared data" control row and its logic."""

    # ── construction ────────────────────────────────────────────────────────
    def _build_slit_row(self, layout) -> None:
        """Create the slit-smearing control row and add it to ``layout``."""
        self._slit_length = 0.0
        self._has_smr_sibling = False
        self._updating_slit_ui = False

        # Compact controls on one row; the (potentially long) status text goes
        # on its own word-wrapped line below so it never forces the control
        # panel wider than the fitting widgets need — narrow panels (Sizes /
        # Modeling) were being stretched by a single wide slit row.
        outer = QVBoxLayout()
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(1)

        row = QHBoxLayout()
        row.setContentsMargins(0, 0, 0, 0)
        self.slit_smear_check = QCheckBox("Slit smeared")
        self.slit_smear_check.setToolTip(
            "When checked, the model is slit smeared (Lake infinite-slit) before\n"
            "comparison with the data, so fitted parameters are ideal-space.\n"
            "For files with both desmeared and slit-smeared copies, this selects\n"
            "which dataset pyIrena uses.  Auto-enabled for slit-smeared data."
        )
        self.slit_smear_check.stateChanged.connect(self._on_slit_smear_changed)
        row.addWidget(self.slit_smear_check)
        row.addSpacing(6)
        row.addWidget(QLabel("SL (1/Å):"))
        self.slit_length_edit = QLineEdit("0.0")
        self.slit_length_edit.setFixedWidth(80)
        self.slit_length_edit.setToolTip(
            "Slit (half-)length dQl in 1/Å.  File-derived; edit to override.\n"
            "Rarely changed — the file value is the most trusted."
        )
        self.slit_length_edit.editingFinished.connect(self._on_slit_length_edited)
        row.addWidget(self.slit_length_edit)
        row.addStretch()
        outer.addLayout(row)

        self.slit_status_label = QLabel("")
        self.slit_status_label.setStyleSheet("color:#666;font-style:italic;font-size:10px;")
        self.slit_status_label.setWordWrap(True)
        outer.addWidget(self.slit_status_label)

        self.slit_row_widget = QWidget()
        self.slit_row_widget.setLayout(outer)
        # Don't let this row dictate the panel's minimum width.
        self.slit_row_widget.setMaximumWidth(16777215)
        self.slit_status_label.setMinimumWidth(0)
        self.slit_row_widget.setVisible(False)   # shown once data is loaded
        layout.addWidget(self.slit_row_widget)

    # ── public accessors panels read when building their fit object ──────────
    def slit_active(self) -> bool:
        """True when smearing should be applied to the model."""
        if not hasattr(self, "slit_smear_check"):
            return False
        return bool(self.slit_smear_check.isChecked()) and self.current_slit_length() > 0

    def current_slit_length(self) -> float:
        """The slit length to use (0.0 => pinhole), read from the edit field."""
        try:
            return float(self.slit_length_edit.text())
        except (ValueError, TypeError, AttributeError):
            return float(getattr(self, "_slit_length", 0.0) or 0.0)

    # ── sync from freshly-loaded data ───────────────────────────────────────
    def _refresh_slit_ui_from_data(self, slit_length, is_slit_smeared, filepath) -> None:
        """Sync the controls to freshly loaded data (call from ``set_data``)."""
        if not hasattr(self, "slit_row_widget"):
            return
        self._slit_length = float(slit_length or 0.0)
        self._has_smr_sibling = False
        if filepath:
            try:
                from pyirena.io.hdf5 import file_has_smr_entry
                p = Path(filepath)
                self._has_smr_sibling = file_has_smr_entry(str(p.parent), p.name)
            except Exception:
                self._has_smr_sibling = False

        self._updating_slit_ui = True
        try:
            self.slit_row_widget.setVisible(True)
            self.slit_smear_check.setChecked(bool(is_slit_smeared))
            self.slit_length_edit.setText(f"{self._slit_length:.5g}")
            if self._has_smr_sibling:
                self.slit_status_label.setText("desmeared + smeared in file")
                self.slit_status_label.setToolTip(
                    "File has both a desmeared and a slit-smeared copy; the "
                    "checkbox selects which pyIrena loads.")
            elif is_slit_smeared:
                self.slit_status_label.setText("slit-smeared data")
                self.slit_status_label.setToolTip("")
            else:
                self.slit_status_label.setText("pinhole data")
                self.slit_status_label.setToolTip("")
        finally:
            self._updating_slit_ui = False
        self._sync_smearing_hook()

    # ── handlers ────────────────────────────────────────────────────────────
    def _on_slit_smear_changed(self, _state) -> None:
        if getattr(self, "_updating_slit_ui", False):
            return
        checked = self.slit_smear_check.isChecked()
        # File carries both copies -> the checkbox selects which to load.
        if self._has_smr_sibling and getattr(self, "_reload_data_with_smearing", None):
            try:
                self._reload_data_with_smearing(checked)
            except Exception:
                # Reload failed — revert the checkbox so the UI matches the
                # still-loaded (unchanged) dataset rather than lying about it.
                self._updating_slit_ui = True
                try:
                    self.slit_smear_check.setChecked(not checked)
                finally:
                    self._updating_slit_ui = False
            return
        self.slit_length_edit.setEnabled(checked)
        self._sync_smearing_hook()
        self._replot_after_slit_change()

    def _on_slit_length_edited(self) -> None:
        if getattr(self, "_updating_slit_ui", False):
            return
        try:
            self._slit_length = float(self.slit_length_edit.text())
        except (ValueError, TypeError):
            pass
        self._sync_smearing_hook()
        if self.slit_active():
            self._replot_after_slit_change()

    # ── optional model-sync hook ─────────────────────────────────────────────
    def _sync_smearing_hook(self) -> None:
        """Call the panel's ``_sync_smearing_to_model()`` if it defines one.

        Panels (e.g. Unified Fit) that hold a persistent model object push the
        current slit state onto it here.  Panels that read ``slit_active()`` /
        ``current_slit_length()`` at fit time (Sizes/Simple/Modeling) don't
        define this and the hook is a no-op.
        """
        fn = getattr(self, "_sync_smearing_to_model", None)
        if callable(fn):
            try:
                fn()
            except Exception:
                pass

    # ── default hooks (panels may override) ─────────────────────────────────
    def _replot_after_slit_change(self) -> None:
        """Redraw after a slit-setting change.  Best-effort default."""
        for name in ("graph_data", "update_plot", "_replot", "replot", "graph"):
            fn = getattr(self, name, None)
            if callable(fn):
                try:
                    fn()
                except Exception:
                    pass
                return

    # panels that support the desmeared/smeared toggle implement:
    #   def _reload_data_with_smearing(self, prefer_slit_smeared: bool) -> None
