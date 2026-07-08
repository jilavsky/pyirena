"""Tests for pyirena.io.text_import — clean_sas_arrays and ensure_nxcansas_sibling."""

import time
from pathlib import Path

import numpy as np
import pytest

from pyirena.io.text_import import (
    clean_sas_arrays,
    converted_sibling_path,
    ensure_nxcansas_sibling,
    _PROVENANCE_ATTR,
)


# ── clean_sas_arrays ──────────────────────────────────────────────────────────

class TestCleanSasArrays:

    def test_all_clean_passthrough(self):
        Q = np.array([0.01, 0.02, 0.03])
        I = np.array([100.0, 50.0, 10.0])
        E = np.array([5.0, 2.5, 0.5])
        Qc, Ic, Ec, dQc, report = clean_sas_arrays(Q, I, E)
        np.testing.assert_array_equal(Qc, Q)
        np.testing.assert_array_equal(Ic, I)
        np.testing.assert_array_equal(Ec, E)
        assert dQc is None
        assert report['n_kept'] == 3
        assert report['n_q_removed'] == 0
        assert report['n_i_removed'] == 0
        assert report['n_e_synthesized'] == 0

    def test_q_zero_removed(self):
        Q = np.array([0.0, 0.01, 0.02])
        I = np.array([1.0, 100.0, 50.0])
        E = np.array([0.05, 5.0, 2.5])
        Qc, Ic, Ec, _, report = clean_sas_arrays(Q, I, E)
        assert len(Qc) == 2
        assert report['n_q_removed'] == 1
        assert report['n_i_removed'] == 0

    def test_q_negative_removed(self):
        Q = np.array([-0.01, 0.01, 0.02])
        I = np.array([1.0, 100.0, 50.0])
        Qc, Ic, Ec, _, report = clean_sas_arrays(Q, I, None)
        assert len(Qc) == 2
        assert report['n_q_removed'] == 1

    def test_beamstop_zeros_removed(self):
        Q = np.array([0.01, 0.02, 0.03, 0.04])
        I = np.array([0.0, 0.0, 50.0, 10.0])  # first two are beamstop zeros
        E = np.array([0.0, 0.0, 2.5, 0.5])
        Qc, Ic, Ec, _, report = clean_sas_arrays(Q, I, E)
        assert len(Qc) == 2
        assert report['n_i_removed'] == 2
        assert np.all(Ic > 0)

    def test_negative_intensity_removed(self):
        Q = np.array([0.01, 0.02, 0.03])
        I = np.array([-1.0, 50.0, 10.0])
        Qc, Ic, _, _, report = clean_sas_arrays(Q, I, None)
        assert len(Qc) == 2
        assert report['n_i_removed'] == 1

    def test_nan_intensity_removed(self):
        Q = np.array([0.01, 0.02, 0.03])
        I = np.array([np.nan, 50.0, 10.0])
        Qc, Ic, _, _, report = clean_sas_arrays(Q, I, None)
        assert len(Qc) == 2

    def test_zero_error_synthesized(self):
        Q = np.array([0.01, 0.02, 0.03])
        I = np.array([100.0, 50.0, 10.0])
        E = np.array([0.0, 2.5, 0.5])  # first error is zero
        Qc, Ic, Ec, _, report = clean_sas_arrays(Q, I, E, error_fraction=0.05)
        assert report['n_e_synthesized'] == 1
        assert Ec[0] == pytest.approx(100.0 * 0.05)
        assert Ec[1] == pytest.approx(2.5)

    def test_no_error_column_all_synthesized(self):
        Q = np.array([0.01, 0.02, 0.03])
        I = np.array([100.0, 50.0, 10.0])
        Qc, Ic, Ec, _, report = clean_sas_arrays(Q, I, None, error_fraction=0.1)
        assert report['n_e_synthesized'] == 3
        np.testing.assert_allclose(Ec, Ic * 0.1)

    def test_dq_aligned(self):
        Q  = np.array([0.0, 0.01, 0.02, 0.03])  # first point removed (Q=0)
        I  = np.array([1.0, 100.0, 50.0, 10.0])
        E  = np.array([0.05, 5.0, 2.5, 0.5])
        dQ = np.array([0.001, 0.001, 0.002, 0.003])
        Qc, Ic, Ec, dQc, report = clean_sas_arrays(Q, I, E, dQ)
        assert len(dQc) == 3
        np.testing.assert_array_equal(dQc, dQ[1:])

    def test_dq_none_when_length_mismatch(self):
        Q  = np.array([0.01, 0.02])
        I  = np.array([100.0, 50.0])
        dQ = np.array([0.001])  # wrong length
        _, _, _, dQc, _ = clean_sas_arrays(Q, I, None, dQ)
        assert dQc is None

    def test_report_keys(self):
        Q = np.array([0.01, 0.02])
        I = np.array([100.0, 50.0])
        _, _, _, _, report = clean_sas_arrays(Q, I, None)
        for key in ('n_original', 'n_kept', 'n_q_removed', 'n_i_removed',
                    'n_e_synthesized', 'error_fraction'):
            assert key in report


# ── converted_sibling_path ────────────────────────────────────────────────────

class TestConvertedSiblingPath:

    def test_dat_to_h5(self):
        p = Path('/tmp/mydata.dat')
        assert converted_sibling_path(p) == Path('/tmp/mydata.h5')

    def test_txt_to_h5(self):
        p = Path('/some/folder/sample_123.txt')
        assert converted_sibling_path(p) == Path('/some/folder/sample_123.h5')


# ── ensure_nxcansas_sibling ───────────────────────────────────────────────────

def _write_text_file(path: Path, rows=None):
    """Write a minimal 3-column SAS text file."""
    if rows is None:
        rows = [
            (0.01, 100.0, 5.0),
            (0.02, 80.0,  4.0),
            (0.05, 30.0,  1.5),
            (0.10, 10.0,  0.5),
        ]
    with open(path, 'w') as f:
        f.write("# Q I dI\n")
        for q, i, e in rows:
            f.write(f"{q}  {i}  {e}\n")


class TestEnsureNxcanSASSibling:

    def test_creates_valid_h5(self, tmp_path):
        txt = tmp_path / 'sample.dat'
        _write_text_file(txt)
        h5 = ensure_nxcansas_sibling(txt)
        assert h5.exists()
        assert h5.suffix == '.h5'
        assert h5.stem == 'sample'

    def test_sibling_is_full_nxcansas(self, tmp_path):
        """The created file must have both sasdata and the provenance marker."""
        import h5py
        txt = tmp_path / 'sample.dat'
        _write_text_file(txt)
        h5 = ensure_nxcansas_sibling(txt)
        with h5py.File(h5, 'r') as f:
            # Provenance marker
            assert _PROVENANCE_ATTR in f.attrs
            # Must contain a Q array somewhere (sasdata or subentry)
            # Walk to find 'I' dataset
            found_intensity = []
            def _visit(name, obj):
                if hasattr(obj, 'shape') and name.endswith('/I'):
                    found_intensity.append(name)
            f.visititems(_visit)
            assert found_intensity, "No 'I' dataset found in HDF5"

    def test_cleans_beamstop_zeros(self, tmp_path):
        """I=0 rows must be removed from the converted file."""
        import h5py
        txt = tmp_path / 'sample.dat'
        _write_text_file(txt, rows=[
            (0.0,  1.0,  0.05),   # Q=0  → removed
            (0.01, 0.0,  0.0),    # I=0  → removed
            (0.02, 80.0, 4.0),
            (0.05, 30.0, 1.5),
        ])
        h5 = ensure_nxcansas_sibling(txt)
        with h5py.File(h5, 'r') as f:
            # Find the Q dataset
            q_data = None
            def _find_q(name, obj):
                nonlocal q_data
                if hasattr(obj, 'shape') and name.endswith('/Q'):
                    q_data = obj[()]
            f.visititems(_find_q)
        assert q_data is not None
        assert len(q_data) == 2          # only the 2 clean rows
        assert np.all(q_data > 0)

    def test_mtime_cache(self, tmp_path):
        """Second call reuses the sibling without rewriting it."""
        txt = tmp_path / 'sample.dat'
        _write_text_file(txt)
        h5 = ensure_nxcansas_sibling(txt)
        mtime1 = h5.stat().st_mtime
        time.sleep(0.05)
        h5b = ensure_nxcansas_sibling(txt)
        mtime2 = h5b.stat().st_mtime
        assert mtime1 == mtime2   # not re-written

    def test_force_reconverts(self, tmp_path):
        """force=True must rewrite even when the sibling is newer."""
        txt = tmp_path / 'sample.dat'
        _write_text_file(txt)
        h5 = ensure_nxcansas_sibling(txt)
        mtime1 = h5.stat().st_mtime
        time.sleep(0.05)
        h5b = ensure_nxcansas_sibling(txt, force=True)
        mtime2 = h5b.stat().st_mtime
        assert mtime2 > mtime1

    def test_collision_guard_uses_fallback(self, tmp_path):
        """If <stem>.h5 is a foreign file, fallback to <stem>_NX.h5."""
        txt = tmp_path / 'sample.dat'
        _write_text_file(txt)
        # Write a foreign HDF5 with the preferred name (no provenance marker)
        import h5py
        foreign = tmp_path / 'sample.h5'
        with h5py.File(foreign, 'w') as f:
            f.attrs['creator'] = 'someone_else'
        h5 = ensure_nxcansas_sibling(txt)
        assert h5.name == 'sample_NX.h5', f"Expected fallback, got {h5.name}"
        assert foreign.exists()   # original untouched
        import h5py
        with h5py.File(foreign, 'r') as f:
            assert f.attrs.get('creator') == 'someone_else'

    def test_raises_on_missing_file(self, tmp_path):
        with pytest.raises(RuntimeError, match="not found"):
            ensure_nxcansas_sibling(tmp_path / 'nonexistent.dat')

    def test_raises_when_no_valid_points(self, tmp_path):
        """All-zero intensities → all removed → RuntimeError."""
        txt = tmp_path / 'bad.dat'
        _write_text_file(txt, rows=[
            (0.01, 0.0, 0.0),
            (0.02, 0.0, 0.0),
        ])
        with pytest.raises(RuntimeError, match="No valid data points"):
            ensure_nxcansas_sibling(txt)

    def test_regression_sibling_has_sasdata_group(self, tmp_path):
        """Regression: the produced file must contain a sasdata group so that
        save_sizes_results (which appends to the file) finds it valid.
        """
        import h5py
        txt = tmp_path / 'sample.dat'
        _write_text_file(txt)
        h5 = ensure_nxcansas_sibling(txt)
        with h5py.File(h5, 'r') as f:
            # Walk the file for any group named 'sasdata'
            sasdata_groups = []
            def _find_sasdata(name, obj):
                if isinstance(obj, h5py.Group) and name.endswith('sasdata'):
                    sasdata_groups.append(name)
            f.visititems(_find_sasdata)
        assert sasdata_groups, (
            "No 'sasdata' group found — the converted file is not a valid "
            "NXcanSAS file and tools will fail to save results into it."
        )
