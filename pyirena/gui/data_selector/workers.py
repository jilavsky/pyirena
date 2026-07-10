"""
pyirena.gui.data_selector.workers — background QThread workers for batch fitting.

Split from the original monolithic data_selector.py (no behavior change).
"""

import os


from pyirena.gui.data_selector._qt import (
    QThread, Signal,
)
from pyirena.batch import (
    fit_unified, fit_sizes, fit_simple_from_config,
    fit_waxs_peaks_from_config, fit_modeling, fit_saxs_morph,
)



class BatchWorker(QThread):
    """
    Background thread for running batch fitting (Unified Fit or Size Distribution)
    on a list of files with a shared pyirena_config.json.
    """
    progress = Signal(str)              # emitted before each file: "Working: N/M — name"
    finished = Signal(int, int, list)   # n_ok, n_fail, per-file messages

    def __init__(self, tool: str, file_paths: list, config_file: str,
                 with_uncertainty: bool = False, n_mc_runs: int = 10, parent=None):
        super().__init__(parent)
        self.tool = tool                # 'unified', 'sizes', or 'simple_fits'
        self.file_paths = file_paths
        self.config_file = config_file
        self.with_uncertainty = with_uncertainty
        self.n_mc_runs = n_mc_runs

    def run(self):
        n_ok, n_fail = 0, 0
        messages = []
        if self.tool == 'unified':
            fit_fn = fit_unified
        elif self.tool == 'modeling':
            fit_fn = fit_modeling
        elif self.tool == 'simple_fits':
            fit_fn = fit_simple_from_config
        elif self.tool == 'waxs_peakfit':
            fit_fn = fit_waxs_peaks_from_config
        elif self.tool == 'saxs_morph':
            fit_fn = fit_saxs_morph
        else:
            fit_fn = fit_sizes
        total = len(self.file_paths)

        # MC uncertainty is supported by unified/modeling/sizes/simple_fits/saxs_morph, not waxs
        mc_kwargs = {}
        if self.with_uncertainty and self.tool != 'waxs_peakfit':
            mc_kwargs = {'with_uncertainty': True, 'n_mc_runs': self.n_mc_runs}

        for i, fp in enumerate(self.file_paths):
            fname = os.path.basename(fp)
            self.progress.emit(f"Working: {i + 1}/{total} — {fname}")
            try:
                result = fit_fn(fp, self.config_file, save_to_nexus=True, **mc_kwargs)
                if result is None:
                    result = {'success': False, 'message': 'fit returned None'}
                if result.get('success', False):
                    n_ok += 1
                    messages.append(f"✓ {fname}")
                else:
                    n_fail += 1
                    messages.append(f"✗ {fname}: {result.get('message', 'fit failed')}")
            except Exception as exc:
                n_fail += 1
                messages.append(f"✗ {fname}: {exc}")

        self.finished.emit(n_ok, n_fail, messages)
