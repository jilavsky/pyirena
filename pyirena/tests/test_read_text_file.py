"""Tests for pyirena.io.hdf5.readTextFile — delimiter auto-detection and
Q-unit conversion."""

import numpy as np
import pytest

from pyirena.io.hdf5 import readTextFile, Q_UNIT_TO_ANGSTROM


ROWS = [
    (0.01, 100.0, 5.0),
    (0.02, 80.0,  4.0),
    (0.05, 30.0,  1.5),
    (0.10, 10.0,  0.5),
]


def _write(path, text):
    path.write_text(text)


class TestReadTextFile:

    def test_whitespace_delimited(self, tmp_path):
        _write(tmp_path / 'a.dat', "# Q I dI\n" + "\n".join(f"{q} {i} {e}" for q, i, e in ROWS))
        data = readTextFile(str(tmp_path), 'a.dat')
        np.testing.assert_allclose(data['Q'], [r[0] for r in ROWS])
        np.testing.assert_allclose(data['Intensity'], [r[1] for r in ROWS])
        np.testing.assert_allclose(data['Error'], [r[2] for r in ROWS])

    def test_comma_delimited_no_header(self, tmp_path):
        _write(tmp_path / 'b.csv', "\n".join(f"{q},{i},{e}" for q, i, e in ROWS))
        data = readTextFile(str(tmp_path), 'b.csv')
        np.testing.assert_allclose(data['Q'], [r[0] for r in ROWS])
        np.testing.assert_allclose(data['Intensity'], [r[1] for r in ROWS])
        np.testing.assert_allclose(data['Error'], [r[2] for r in ROWS])

    def test_comma_delimited_with_text_header(self, tmp_path):
        text = "Q,I,dI\n" + "\n".join(f"{q},{i},{e}" for q, i, e in ROWS)
        _write(tmp_path / 'c.csv', text)
        data = readTextFile(str(tmp_path), 'c.csv')
        assert len(data['Q']) == len(ROWS)
        np.testing.assert_allclose(data['Q'], [r[0] for r in ROWS])

    def test_comma_delimited_two_columns_synthesizes_error(self, tmp_path):
        text = "\n".join(f"{q},{i}" for q, i, _ in ROWS)
        _write(tmp_path / 'd.csv', text)
        data = readTextFile(str(tmp_path), 'd.csv', error_fraction=0.1)
        np.testing.assert_allclose(data['Error'], np.array([r[1] for r in ROWS]) * 0.1)

    def test_missing_file_returns_none(self, tmp_path):
        assert readTextFile(str(tmp_path), 'nope.csv') is None


class TestQUnitConversion:

    @pytest.mark.parametrize('q_unit', ['1/nm', '1/pm', '1/um', '1/mm'])
    def test_q_and_dq_scaled_to_angstrom(self, tmp_path, q_unit):
        rows = [(1.0, 100.0, 5.0, 0.01), (2.0, 50.0, 2.5, 0.02)]
        text = "\n".join(f"{q} {i} {e} {dq}" for q, i, e, dq in rows)
        _write(tmp_path / 'a.dat', text)
        data = readTextFile(str(tmp_path), 'a.dat', q_unit=q_unit)
        scale = Q_UNIT_TO_ANGSTROM[q_unit]
        np.testing.assert_allclose(data['Q'], [r[0] * scale for r in rows])
        np.testing.assert_allclose(data['dQ'], [r[3] * scale for r in rows])
        # Intensity/Error are untouched by Q-unit conversion
        np.testing.assert_allclose(data['Intensity'], [r[1] for r in rows])
        np.testing.assert_allclose(data['Error'], [r[2] for r in rows])

    def test_default_unit_is_no_op(self, tmp_path):
        rows = [(1.0, 100.0, 5.0), (2.0, 50.0, 2.5)]
        text = "\n".join(f"{q} {i} {e}" for q, i, e in rows)
        _write(tmp_path / 'a.dat', text)
        data = readTextFile(str(tmp_path), 'a.dat')
        np.testing.assert_allclose(data['Q'], [r[0] for r in rows])

    def test_unknown_q_unit_raises(self, tmp_path):
        _write(tmp_path / 'a.dat', "1.0 100.0 5.0\n")
        with pytest.raises(ValueError, match="Unknown q_unit"):
            readTextFile(str(tmp_path), 'a.dat', q_unit='1/parsec')
