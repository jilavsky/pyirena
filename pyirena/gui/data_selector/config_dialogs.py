"""
pyirena.gui.data_selector.config_dialogs — Configure… and Manage Config… dialogs.

Split from the original monolithic data_selector.py (no behavior change).
"""



from pyirena.gui.data_selector._qt import (
    QVBoxLayout, QHBoxLayout, QGridLayout, QPushButton, QLabel, QLineEdit, QComboBox, QMessageBox, QFrame, QDialog, QFormLayout, QDialogButtonBox, QGroupBox, QCheckBox, QColorDialog, QDoubleValidator,
)



class ConfigManagerDialog(QDialog):
    """Dialog for inspecting and pruning sections from a pyirena_config.json file.

    Shows each tool section as a labelled checkbox.  Unchecking a section and
    clicking Save removes that section from the file.  No values are exposed for
    editing, avoiding accidental corruption.
    """

    # Known tool keys in the order we want to display them
    _KNOWN_TOOLS = [
        'unified_fit', 'sizes', 'modeling', 'simple_fits', 'waxs_peakfit',
        'saxs_morph', 'data_merge',
    ]

    def __init__(self, config_path: str, parent=None):
        super().__init__(parent)
        self.config_path = config_path
        self.setWindowTitle("Manage Config Sections")
        self.setMinimumWidth(420)

        import json
        try:
            with open(config_path, 'r') as f:
                self._config = json.load(f)
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Cannot read config file:\n{e}")
            self._config = {}

        self._checkboxes: dict = {}  # key → QCheckBox

        layout = QVBoxLayout(self)

        # File path label
        path_lbl = QLabel(f"<small>{config_path}</small>")
        path_lbl.setWordWrap(True)
        layout.addWidget(path_lbl)

        # Written-by / modified info from _pyirena_config header
        meta = self._config.get('_pyirena_config', {})
        if meta:
            modified = meta.get('modified', '')
            written_by = meta.get('written_by', '')
            info = []
            if written_by:
                info.append(f"Written by: {written_by}")
            if modified:
                info.append(f"Last modified: {modified[:19]}")
            if info:
                meta_lbl = QLabel("<small>" + " &nbsp;|&nbsp; ".join(info) + "</small>")
                meta_lbl.setStyleSheet("color: #7f8c8d;")
                layout.addWidget(meta_lbl)

        sep = QFrame()
        sep.setFrameShape(QFrame.Shape.HLine)
        sep.setFrameShadow(QFrame.Shadow.Sunken)
        layout.addWidget(sep)

        # Instruction
        instr = QLabel(
            "Uncheck sections you want to remove, then click <b>Save</b>.\n"
            "Removed sections will no longer be run by fit_pyirena()."
        )
        instr.setWordWrap(True)
        layout.addWidget(instr)

        # One checkbox per tool section found in the file
        tool_keys_in_file = [k for k in self._KNOWN_TOOLS if k in self._config]
        # Also show any unknown tool keys (future-proofing)
        extra_keys = [k for k in self._config
                      if k not in self._KNOWN_TOOLS and k != '_pyirena_config']
        all_keys = tool_keys_in_file + extra_keys

        if not all_keys:
            layout.addWidget(QLabel("No tool sections found in this file."))
        else:
            grid = QGridLayout()
            grid.setColumnStretch(0, 1)
            for row, key in enumerate(all_keys):
                cb = QCheckBox(key)
                cb.setChecked(True)
                cb.setStyleSheet("font-weight: bold;")
                self._checkboxes[key] = cb
                grid.addWidget(cb, row, 0)
                # Show a timestamp hint if stored per-section (not currently written,
                # but future-safe) or fall back to file-level modified date
                ts = ''
                sec = self._config.get(key, {})
                if isinstance(sec, dict):
                    ts = sec.get('_modified', meta.get('modified', ''))
                if ts:
                    ts_lbl = QLabel(f"<small>{ts[:19]}</small>")
                    ts_lbl.setStyleSheet("color: #95a5a6;")
                    grid.addWidget(ts_lbl, row, 1)
            layout.addLayout(grid)

        sep2 = QFrame()
        sep2.setFrameShape(QFrame.Shape.HLine)
        sep2.setFrameShadow(QFrame.Shadow.Sunken)
        layout.addWidget(sep2)

        btn_row = QHBoxLayout()
        save_btn = QPushButton("Save")
        save_btn.setDefault(True)
        save_btn.clicked.connect(self._save)
        cancel_btn = QPushButton("Cancel")
        cancel_btn.clicked.connect(self.reject)
        btn_row.addStretch()
        btn_row.addWidget(save_btn)
        btn_row.addWidget(cancel_btn)
        layout.addLayout(btn_row)

    def _save(self):
        import json
        import datetime
        to_remove = [key for key, cb in self._checkboxes.items() if not cb.isChecked()]
        if not to_remove:
            self.accept()
            return
        for key in to_remove:
            self._config.pop(key, None)
        # Update modification timestamp
        if '_pyirena_config' in self._config:
            self._config['_pyirena_config']['modified'] = (
                datetime.datetime.now().isoformat(timespec='seconds')
            )
        try:
            with open(self.config_path, 'w') as f:
                json.dump(self._config, f, indent=2)
            QMessageBox.information(
                self, "Saved",
                f"Removed {len(to_remove)} section(s): {', '.join(to_remove)}"
            )
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Could not save config:\n{e}")
            return
        self.accept()


class DataSelectorConfigDialog(QDialog):
    """
    Extensible configuration dialog for the Data Selector.

    Settings are defined as a list of field specifications (FIELD_SPECS).
    Adding a new configurable parameter only requires adding one entry to that list.

    Supported field types
    ---------------------
    'float'  — QLineEdit with QDoubleValidator
    'int'    — QLineEdit with integer validation
    'str'    — plain QLineEdit
    'bool'   — QCheckBox
    'color'  — QPushButton that opens QColorDialog (stores hex color string)
    'choice' — QComboBox of (label, value) pairs; stores the value
    """

    # ---------------------------------------------------------------------------
    # Field specifications — add new settings here
    # ---------------------------------------------------------------------------
    FIELD_SPECS = [
        {
            'group':    'Text File Options',
            'key':      'error_fraction',
            'label':    'Generated uncertainty fraction',
            'tooltip':  (
                'When a text file has only two columns (Q and I) and no uncertainty\n'
                'column, the uncertainty is generated as:  σ = I × this_value\n'
                'Default: 0.05  (5 % of intensity)'
            ),
            'type':     'float',
            'default':  0.05,
            'min':      0.0,
            'max':      100.0,
            'decimals': 4,
        },
        {
            'group':   'Graph Options',
            'key':     'max_legend_items',
            'label':   'Maximum items in legend',
            'tooltip': (
                'When more datasets than this limit are plotted, the legend shows\n'
                'only the first, last, and evenly spaced items in between.\n'
                'Default: 12'
            ),
            'type':    'int',
            'default': 12,
            'min':     2,
            'max':     200,
            'decimals': 0,
        },
        {
            'group':   'Batch Script Options',
            'key':     'batch_mc_uncertainty',
            'label':   'Run with MC uncertainty',
            'tooltip': (
                'When enabled, batch scripts for Unified Fit, Size Distribution\n'
                'and Simple Fits will calculate parameter uncertainties via Monte\n'
                'Carlo after each fit.  This significantly increases run time.\n'
                '(WAXS Peak Fit uses least-squares uncertainties — not affected.)\n'
                'Default: off'
            ),
            'type':    'bool',
            'default': False,
        },
        {
            'group':   'Batch Script Options',
            'key':     'batch_mc_n_runs',
            'label':   'MC passes',
            'tooltip': (
                'Number of noise-perturbed fits used for MC uncertainty estimation.\n'
                'Default: 10'
            ),
            'type':    'int',
            'default': 10,
            'min':     1,
            'max':     500,
            'decimals': 0,
        },
        {
            'group':   'ASCII Export Options',
            'key':     'ascii_delimiter',
            'label':   'Column delimiter',
            'tooltip': (
                'Character separating columns in exported .dat files.\n'
                'Space (default): widest legacy-tool compatibility (LoadWave/J,\n'
                'awk, np.loadtxt, gnuplot).\n'
                'Comma: CSV-style; readable by Excel.'
            ),
            'type':    'choice',
            'choices': [('Space', ' '), ('Comma', ',')],
            'default': ' ',
        },
        {
            'group':   'ASCII Export Options',
            'key':     'ascii_precision',
            'label':   'Significant figures',
            'tooltip': (
                'Number of significant figures used for numeric columns and\n'
                'header parameters.\n'
                '7 (default): single-precision-safe; old Fortran codes parse\n'
                'these reliably.\n'
                '12: double-precision; for tools needing full accuracy.'
            ),
            'type':    'choice',
            'choices': [('7 (single-prec safe)', 7), ('12 (double-prec)', 12)],
            'default': 7,
        },
        {
            'group':   'ASCII Export Options',
            'key':     'ascii_include_header',
            'label':   'Write metadata header',
            'tooltip': (
                'When checked, each .dat file starts with up to 25 lines of\n'
                "'# key = value' metadata (sample, instrument, wavelength,\n"
                'energy, fit parameters, etc.).\n'
                'Uncheck for bare numeric data (some old tools refuse # lines).'
            ),
            'type':    'bool',
            'default': True,
        },
        {
            'group':   'ASCII Export Options',
            'key':     'ascii_include_models',
            'label':   'Also write model curves',
            'tooltip': (
                "When checked, for every fit-result checkbox that's enabled,\n"
                'an additional 4-column file ({stem}_<acronym>.dat) is written\n'
                'with Q, I_model, I_data, dI columns and model parameters in\n'
                'the header.  Acronyms: _unif, _simp, _mod, _sd, _waxs.\n'
                'Files without that result type are silently skipped.'
            ),
            'type':    'bool',
            'default': True,
        },
        # -----------------------------------------------------------------------
        # Future settings — just append a dict here, no other code changes needed
        # -----------------------------------------------------------------------
        # {
        #     'group':   'Text File Options',
        #     'key':     'q_units_scale',
        #     'label':   'Q unit scale factor',
        #     'tooltip': 'Multiply Q by this factor on load (1.0 = no change)',
        #     'type':    'float',
        #     'default': 1.0,
        #     'min':     0.0,
        #     'max':     1000.0,
        #     'decimals': 6,
        # },
        # {
        #     'group':   'Display',
        #     'key':     'plot_color',
        #     'label':   'Default plot color',
        #     'tooltip': 'Color used for single-file plots',
        #     'type':    'color',
        #     'default': '#3498db',
        # },
    ]

    def __init__(self, current_values: dict, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Data Selector — Configure")
        self.setMinimumWidth(420)
        self._widgets = {}   # key -> (field_type, widget)
        self._init_ui(current_values)

    def _init_ui(self, current_values: dict):
        outer = QVBoxLayout()
        outer.setSpacing(12)

        # Group fields by their 'group' key
        groups = {}
        for spec in self.FIELD_SPECS:
            g = spec.get('group', 'General')
            groups.setdefault(g, []).append(spec)

        for group_name, specs in groups.items():
            box = QGroupBox(group_name)
            form = QFormLayout()
            form.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapLongRows)

            for spec in specs:
                key       = spec['key']
                label_txt = spec['label']
                ftype     = spec['type']
                value     = current_values.get(key, spec.get('default', ''))
                tooltip   = spec.get('tooltip', '')

                if ftype in ('float', 'int'):
                    widget = QLineEdit(str(value))
                    validator = QDoubleValidator(
                        float(spec.get('min', -1e300)),
                        float(spec.get('max',  1e300)),
                        int(spec.get('decimals', 6)),
                    )
                    widget.setValidator(validator)
                    widget.setMaximumWidth(120)

                elif ftype == 'bool':
                    widget = QCheckBox()
                    widget.setChecked(bool(value))

                elif ftype == 'color':
                    widget = QPushButton()
                    widget._color = str(value)
                    widget.setStyleSheet(f"background-color: {value};")
                    widget.setFixedSize(60, 24)
                    widget.clicked.connect(
                        lambda checked, btn=widget: self._pick_color(btn)
                    )

                elif ftype == 'choice':
                    widget = QComboBox()
                    choices = spec.get('choices', [])
                    widget._values = [v for _, v in choices]
                    for label, _val in choices:
                        widget.addItem(label)
                    # Select the entry whose value matches *value*
                    selected_idx = 0
                    for i, v in enumerate(widget._values):
                        if v == value:
                            selected_idx = i
                            break
                    widget.setCurrentIndex(selected_idx)

                else:   # 'str'
                    widget = QLineEdit(str(value))

                if tooltip:
                    widget.setToolTip(tooltip)

                lbl = QLabel(label_txt)
                if tooltip:
                    lbl.setToolTip(tooltip)

                form.addRow(lbl, widget)
                self._widgets[key] = (ftype, widget)

            box.setLayout(form)
            outer.addWidget(box)

        # OK / Cancel buttons
        btn_box = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        btn_box.accepted.connect(self.accept)
        btn_box.rejected.connect(self.reject)
        outer.addWidget(btn_box)

        self.setLayout(outer)

    def _pick_color(self, btn):
        color = QColorDialog.getColor(parent=self)
        if color.isValid():
            btn._color = color.name()
            btn.setStyleSheet(f"background-color: {color.name()};")

    def get_values(self) -> dict:
        """Return validated values from all widgets keyed by field key."""
        result = {}
        for spec in self.FIELD_SPECS:
            key   = spec['key']
            ftype = spec['type']
            _, widget = self._widgets[key]

            if ftype in ('float', 'int'):
                try:
                    result[key] = float(widget.text()) if ftype == 'float' else int(widget.text())
                except ValueError:
                    result[key] = spec.get('default', 0)
            elif ftype == 'bool':
                result[key] = widget.isChecked()
            elif ftype == 'color':
                result[key] = widget._color
            elif ftype == 'choice':
                idx = widget.currentIndex()
                values = getattr(widget, '_values', [])
                if 0 <= idx < len(values):
                    result[key] = values[idx]
                else:
                    result[key] = spec.get('default', '')
            else:
                result[key] = widget.text()
        return result
