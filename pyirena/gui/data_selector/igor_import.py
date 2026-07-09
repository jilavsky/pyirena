"""
pyirena.gui.data_selector.igor_import — Igor pxp/h5xp import dialog.

Split from the original monolithic data_selector.py (no behavior change).
"""

from pathlib import Path
from typing import Dict, List


from pyirena.gui.data_selector._qt import (
    QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QLabel, QLineEdit, QFileDialog, QDialog, QFormLayout, QDialogButtonBox, QGroupBox, QCheckBox, Qt,
)



class _IgorImportDialog(QDialog):
    """Modal dialog asking the user where to send extracted NeXus files and
    which techniques to keep, before the actual Igor → NeXus import runs.

    Used only by :meth:`DataSelectorPanel.launch_igor_import`. The dialog is
    short and purely a settings prompt; the real work happens in
    :func:`pyirena.batch.igor_to_nexus`, which dispatches to either the
    .pxp or .h5xp reader based on the input extension.
    """

    def __init__(self, pxp_path, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Import Igor Experiment")
        self.pxp_path = Path(pxp_path)
        self.setMinimumWidth(540)

        layout = QVBoxLayout(self)

        info = QLabel(
            f"<b>Source:</b> {self.pxp_path.name}<br>"
            f"<b>Folder:</b> {self.pxp_path.parent}"
        )
        info.setTextFormat(Qt.TextFormat.RichText)
        layout.addWidget(info)

        layout.addSpacing(8)

        # Output folder row
        form = QFormLayout()
        default_out = self.pxp_path.with_name(f"{self.pxp_path.stem}_data")
        self._out_edit = QLineEdit(str(default_out))
        out_row = QHBoxLayout()
        out_row.addWidget(self._out_edit, stretch=1)
        out_browse = QPushButton("Browse…")
        out_browse.clicked.connect(self._browse_output)
        out_row.addWidget(out_browse)
        out_w = QWidget(); out_w.setLayout(out_row)
        form.addRow("Output folder:", out_w)
        layout.addLayout(form)

        layout.addSpacing(4)

        # Techniques group
        grp = QGroupBox("Techniques to export")
        gv = QVBoxLayout(grp)
        self._cb_usaxs = QCheckBox("USAXS")
        self._cb_saxs  = QCheckBox("SAXS")
        self._cb_waxs  = QCheckBox("WAXS")
        for cb in (self._cb_usaxs, self._cb_saxs, self._cb_waxs):
            cb.setChecked(True)
            gv.addWidget(cb)
        layout.addWidget(grp)

        # Overwrite option
        self._cb_overwrite = QCheckBox(
            "Overwrite existing files  (off = append _2, _3, … to keep both)"
        )
        layout.addWidget(self._cb_overwrite)

        layout.addStretch()

        # Buttons
        btn_box = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok
            | QDialogButtonBox.StandardButton.Cancel
        )
        btn_box.button(QDialogButtonBox.StandardButton.Ok).setText("Import")
        btn_box.accepted.connect(self.accept)
        btn_box.rejected.connect(self.reject)
        layout.addWidget(btn_box)

    def _browse_output(self):
        start = self._out_edit.text() or str(self.pxp_path.parent)
        folder = QFileDialog.getExistingDirectory(
            self, "Output folder", start, QFileDialog.Option.ShowDirsOnly
        )
        if folder:
            self._out_edit.setText(folder)

    def options(self) -> Dict:
        """Return a dict of the user's choices."""
        techs: List[str] = []
        if self._cb_usaxs.isChecked():
            techs.append("USAXS")
        if self._cb_saxs.isChecked():
            techs.append("SAXS")
        if self._cb_waxs.isChecked():
            techs.append("WAXS")
        return {
            'output_folder': self._out_edit.text().strip() or None,
            'techniques':    techs if techs else None,
            'overwrite':     self._cb_overwrite.isChecked(),
        }
