"""
Regression tests: Modeling form-factor / structure-factor parameters must
propagate into the Markdown report ("Create Report") and the CSV table
("Tabulate Results"), not only the HDF5 file. Covers the core-shell-shell
sphere whose parameters were previously dropped by both outputs.
"""

import os

import numpy as np
import pytest


# A synthetic loaded-modeling-results dict in the shape produced by
# pyirena.io.nxcansas_modeling.load_modeling_results().
def _css_modeling_results():
    return {
        'chi_squared': 1.23,
        'background': 0.001,
        'q_min': 0.005,
        'q_max': 0.35,
        'populations': [{
            'enabled': True,
            'pop_type': 'size_dist',
            'population_index': 1,
            'label': 'shell test',
            'dist_type': 'gauss',
            'form_factor': 'css_sphere_by_core',
            'structure_factor': 'none',
            'dist_params': {'mean_size': 80.0, 'width': 4.0},
            'ff_params': {
                'sld_core': 10.0, 'sld_shell1': 1.0, 'sld_shell2': 5.0,
                'sld_solvent': 9.46, 't_shell1': 21.0, 't_shell2': 22.0,
            },
            'sf_params': {},
            'scale': 0.01,
            'contrast': 1.0,
            'derived': {'volume_fraction': 0.01, 'vol_mean_r': 80.0},
        }],
    }


class TestReportIncludesFormFactorParams:
    def test_report_lists_ff_and_dist_params(self):
        from pyirena.gui.data_selector import _build_report

        md = _build_report('sample.h5', modeling_results=_css_modeling_results())

        # Form-factor parameters must appear (the core complaint).
        for key in ('sld_core', 'sld_shell1', 'sld_shell2', 'sld_solvent',
                    't_shell1', 't_shell2'):
            assert key in md, f"{key} missing from report"
        # A couple of the actual values should be present.
        assert '21' in md and '22' in md
        # Distribution parameters too.
        assert 'mean_size' in md


def _make_css_modeling_result():
    """A real ModelingResult for a core-shell-shell population (nothing fit,
    so fit() just evaluates and returns quickly)."""
    from pyirena.core.modeling import (
        ModelingEngine, ModelingConfig, SizeDistPopulation,
    )
    q = np.logspace(np.log10(0.005), np.log10(0.35), 120)
    pop = SizeDistPopulation()
    pop.enabled = True
    pop.dist_type = 'gauss'
    pop.dist_params = {'mean_size': 80.0, 'width': 4.0}
    pop.dist_params_fit = {'mean_size': False, 'width': False}
    pop.form_factor = 'css_sphere_by_core'
    pop.ff_params = {
        'sld_core': 10.0, 'sld_shell1': 1.0, 'sld_shell2': 5.0,
        'sld_solvent': 9.46, 't_shell1': 21.0, 't_shell2': 22.0,
    }
    pop.ff_params_fit = {k: False for k in pop.ff_params}
    pop.scale = 0.01
    pop.fit_scale = False
    pop.contrast = 1.0
    pop.fit_contrast = False
    cfg = ModelingConfig(populations=[pop], background=0.0,
                         fit_background=False, q_min=q.min(), q_max=q.max())
    eng = ModelingEngine()
    I = eng.total_intensity(cfg, q, use_cache=False)[0]
    dI = 0.02 * I + 1e-7
    return eng.fit(cfg, q, I, dI)


class TestCsvIncludesFormFactorParams:
    def test_tabulate_emits_ff_columns(self, tmp_path, monkeypatch):
        pytest.importorskip("pyirena.gui.data_selector")
        try:
            try:
                from PySide6.QtWidgets import QApplication, QListWidgetItem
            except ImportError:
                from PyQt6.QtWidgets import QApplication, QListWidgetItem
        except Exception:
            pytest.skip("Qt not available")
        os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

        from pyirena.io.nxcansas_modeling import save_modeling_results
        from pyirena.gui import data_selector as ds

        # Write a real HDF5 file with a core-shell-shell modeling result.
        h5 = tmp_path / 'sample.h5'
        save_modeling_results(h5, _make_css_modeling_result())

        app = QApplication.instance() or QApplication([])
        panel = ds.DataSelectorPanel()
        panel.current_folder = str(tmp_path)
        panel.file_list.clear()
        panel.file_list.addItem(QListWidgetItem('sample.h5'))
        panel.file_list.item(0).setSelected(True)
        # Only tabulate modeling results.
        for cb_name in ('unified_fit_result_checkbox', 'size_dist_checkbox',
                        'simple_fits_checkbox', 'waxs_peakfit_checkbox',
                        'saxs_morph_checkbox'):
            getattr(panel, cb_name).setChecked(False)
        panel.modeling_checkbox.setChecked(True)

        captured = {}

        def _capture(_self, headers, rows, default_path):
            captured['headers'] = headers
            captured['rows'] = rows

        monkeypatch.setattr(ds.TabulateResultsWindow, 'set_data', _capture)
        panel.tabulate_results()

        headers = captured.get('headers', [])
        # Dynamic ff_ columns for the css keys must be present.
        for key in ('sld_core', 'sld_shell1', 'sld_shell2',
                    't_shell1', 't_shell2'):
            assert f'MOD_P1_ff_{key}' in headers, f'MOD_P1_ff_{key} missing'
        # And the value must be populated in the single data row.
        col = headers.index('MOD_P1_ff_t_shell1')
        assert captured['rows'][0][col] == pytest.approx(21.0)
