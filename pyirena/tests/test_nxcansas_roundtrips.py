"""
Round-trip tests for the io/nxcansas_* save/load pairs.

These guard the on-disk NXcanSAS file format: every writer is exercised
against its reader (or h5py directly) on a tmp_path file, and the values
that come back must equal the values that went in.
"""

import numpy as np
import pytest
import h5py

from pyirena.io.hdf5 import find_matching_groups
from pyirena.io.nxcansas_unified import (
    create_nxcansas_file,
    save_unified_fit_results,
    load_unified_fit_results,
)


def _sasdata(f):
    """Locate the sasdata group the same way the readers do."""
    paths = find_matching_groups(
        f, required_attributes={"canSAS_class": "SASdata"}, required_items={}
    )
    assert paths, "no SASdata group found"
    return f[paths[0]]


# ---------------------------------------------------------------------------
# Shared synthetic data
# ---------------------------------------------------------------------------

@pytest.fixture
def q_i_err():
    q = np.logspace(-3, -0.5, 200)
    intensity = 100.0 * np.exp(-q**2 * 50.0**2 / 3.0) + 1e-4 * q**-4 + 0.01
    error = 0.02 * intensity
    return q, intensity, error


@pytest.fixture
def nx_file(tmp_path, q_i_err):
    """A fresh NXcanSAS data file to attach results to."""
    q, intensity, error = q_i_err
    path = tmp_path / "sample.h5"
    create_nxcansas_file(path, q, intensity, error=error, sample_name="sample")
    return path


# ---------------------------------------------------------------------------
# create_nxcansas_file
# ---------------------------------------------------------------------------

class TestCreateNxcansasFile:
    def test_arrays_round_trip(self, nx_file, q_i_err):
        q, intensity, error = q_i_err
        with h5py.File(nx_file, "r") as f:
            grp = _sasdata(f)
            np.testing.assert_allclose(grp["Q"][:], q)
            np.testing.assert_allclose(grp["I"][:], intensity)
            np.testing.assert_allclose(grp["Idev"][:], error)

    def test_canSAS_attributes(self, nx_file):
        with h5py.File(nx_file, "r") as f:
            assert f["entry"].attrs["NX_class"] == "NXentry"
            assert _sasdata(f).attrs["canSAS_class"] == "SASdata"

    def test_dq_written_when_given(self, tmp_path, q_i_err):
        q, intensity, error = q_i_err
        dq = 0.01 * q
        path = tmp_path / "with_dq.h5"
        create_nxcansas_file(path, q, intensity, error=error, dq=dq)
        with h5py.File(path, "r") as f:
            np.testing.assert_allclose(_sasdata(f)["Qdev"][:], dq)

    def test_core_loader_reads_it_back(self, nx_file, q_i_err):
        from pyirena.core.unified import load_data_from_nxcansas
        q, intensity, error = q_i_err
        data = load_data_from_nxcansas(str(nx_file))
        np.testing.assert_allclose(data["Q"], q)
        np.testing.assert_allclose(data["Intensity"], intensity)


# ---------------------------------------------------------------------------
# Unified Fit results
# ---------------------------------------------------------------------------

def _unified_levels():
    return [
        {"G": 100.0, "Rg": 50.0, "B": 1e-4, "P": 4.0, "RgCutoff": 0.0,
         "ETA": 0.0, "PACK": 0.0, "correlated": False},
        {"G": 5000.0, "Rg": 400.0, "B": 1e-8, "P": 3.2, "RgCutoff": 50.0,
         "ETA": 0.0, "PACK": 0.0, "correlated": False},
    ]


class TestUnifiedFitRoundTrip:
    def test_round_trip(self, nx_file, q_i_err):
        q, intensity, error = q_i_err
        model = intensity * 1.01
        save_unified_fit_results(
            nx_file, q, intensity, model, intensity - model,
            levels=_unified_levels(), background=0.03,
            chi_squared=1.7, num_levels=2, error=error,
        )
        loaded = load_unified_fit_results(nx_file)
        assert loaded["num_levels"] == 2
        assert loaded["background"] == pytest.approx(0.03)
        assert loaded["chi_squared"] == pytest.approx(1.7)
        np.testing.assert_allclose(loaded["Q"], q)
        np.testing.assert_allclose(loaded["intensity_model"], model)
        for sent, back in zip(_unified_levels(), loaded["levels"]):
            for key in ("G", "Rg", "B", "P"):
                assert back[key] == pytest.approx(sent[key]), key

    def test_second_save_overwrites(self, nx_file, q_i_err):
        q, intensity, error = q_i_err
        model = intensity * 1.01
        for chi2 in (5.0, 1.2):
            save_unified_fit_results(
                nx_file, q, intensity, model, intensity - model,
                levels=_unified_levels()[:1], background=0.0,
                chi_squared=chi2, num_levels=1, error=error,
            )
        loaded = load_unified_fit_results(nx_file)
        assert loaded["chi_squared"] == pytest.approx(1.2)
        assert loaded["num_levels"] == 1

    def test_load_result_dispatcher(self, nx_file, q_i_err):
        from pyirena.io.results import load_result
        q, intensity, error = q_i_err
        model = intensity * 1.01
        save_unified_fit_results(
            nx_file, q, intensity, model, intensity - model,
            levels=_unified_levels(), background=0.03,
            chi_squared=1.7, num_levels=2, error=error,
        )
        r = load_result(nx_file, "unified_fit")
        assert r["found"] is True
        assert r["chi_squared"] == pytest.approx(1.7)
        assert len(r["levels"]) == 2

    def test_load_result_not_found(self, nx_file):
        from pyirena.io.results import load_result
        r = load_result(nx_file, "unified_fit")
        assert r["found"] is False


# ---------------------------------------------------------------------------
# Size distribution results
# ---------------------------------------------------------------------------

class TestSizesRoundTrip:
    def test_round_trip(self, nx_file, q_i_err):
        from pyirena.io.nxcansas_sizes import save_sizes_results, load_sizes_results
        q, intensity, error = q_i_err
        r_grid = np.logspace(1, 3, 60)
        dist = np.exp(-((np.log(r_grid) - np.log(100.0)) ** 2))
        params = {
            "method": "regularization",
            "shape": "spheroid",
            "volume_fraction": 0.0123,
            "chi_squared": 0.98,
            "aspect_ratio": 1.0,
        }
        save_sizes_results(
            nx_file, q, intensity, intensity * 0.99, intensity * 0.01,
            r_grid, dist, params, intensity_error=error,
        )
        loaded = load_sizes_results(nx_file)
        np.testing.assert_allclose(loaded["r_grid"], r_grid)
        np.testing.assert_allclose(loaded["distribution"], dist)
        # loader returns scalar metadata flat, not under 'params'
        assert loaded["method"] == "regularization"
        assert loaded["volume_fraction"] == pytest.approx(0.0123)

    def test_load_result_dispatcher(self, nx_file, q_i_err):
        from pyirena.io.nxcansas_sizes import save_sizes_results
        from pyirena.io.results import load_result
        q, intensity, error = q_i_err
        r_grid = np.logspace(1, 3, 60)
        dist = np.exp(-((np.log(r_grid) - np.log(100.0)) ** 2))
        save_sizes_results(
            nx_file, q, intensity, intensity * 0.99, intensity * 0.01,
            r_grid, dist,
            {"method": "regularization", "volume_fraction": 0.0123},
            intensity_error=error,
        )
        r = load_result(nx_file, "size_distribution")
        assert r["found"] is True
        np.testing.assert_allclose(r["r_grid"], r_grid)


# ---------------------------------------------------------------------------
# Simple fits results (via a real Guinier fit — fast)
# ---------------------------------------------------------------------------

class TestSimpleFitRoundTrip:
    def test_round_trip(self, nx_file):
        from pyirena.core.simple_fits import SimpleFitModel
        from pyirena.io.nxcansas_simple_fits import (
            save_simple_fit_results, load_simple_fit_results,
        )
        rng = np.random.default_rng(0)
        q = np.linspace(0.002, 0.02, 80)
        true_G, true_Rg = 500.0, 40.0
        intensity = true_G * np.exp(-q**2 * true_Rg**2 / 3.0)
        intensity *= 1.0 + 0.01 * rng.standard_normal(len(q))
        err = 0.01 * intensity

        m = SimpleFitModel()          # defaults to the Guinier model
        result = m.fit(q, intensity, err)
        assert result.get("success", True)

        save_simple_fit_results(nx_file, result, model_obj=m,
                                intensity_data=intensity, intensity_error=err)
        loaded = load_simple_fit_results(nx_file)
        assert loaded["model"] == "Guinier"
        # Fitted Rg must survive the round trip and be near the truth
        rg = loaded["params"].get("Rg")
        assert rg == pytest.approx(true_Rg, rel=0.05)


# ---------------------------------------------------------------------------
# WAXS peak fit results
# ---------------------------------------------------------------------------

class TestWaxsPeakFitRoundTrip:
    def test_round_trip(self, nx_file):
        from pyirena.io.nxcansas_waxs_peakfit import (
            save_waxs_peakfit_results, load_waxs_peakfit_results,
        )
        q = np.linspace(1.0, 4.0, 300)
        # save() expects Igor-style parameter dicts: {name: {"value":..., "lo":..., "hi":...}}
        peaks = [
            {"shape": "Gauss",
             "A": {"value": 50.0}, "Q0": {"value": 2.0}, "FWHM": {"value": 0.05}},
            {"shape": "Lorentz",
             "A": {"value": 20.0}, "Q0": {"value": 3.1}, "FWHM": {"value": 0.08}},
        ]
        result = {
            "peaks": peaks,
            "peaks_std": [{"Q0": 0.001}, {"Q0": 0.002}],
            "bg_shape": "Constant",
            "bg_params": {"bg0": {"value": 1.5}},
            "bg_params_std": {"bg0": 0.1},
            "chi2": 2.5,
            "reduced_chi2": 1.1,
            "dof": 290,
            "I_model": np.ones_like(q),
        }
        save_waxs_peakfit_results(nx_file, result, q)
        loaded = load_waxs_peakfit_results(nx_file)
        assert len(loaded["peaks"]) == 2
        p0 = loaded["peaks"][0]
        assert p0["shape"] == "Gauss"
        # loader returns Igor-style {name: {"value": ...}} dicts, like save input
        assert p0["Q0"]["value"] == pytest.approx(2.0)
        assert p0["A"]["value"] == pytest.approx(50.0)
        assert loaded["bg_params"]["bg0"]["value"] == pytest.approx(1.5)
        # a derived peak area must always be present
        assert p0["area"] > 0


# ---------------------------------------------------------------------------
# Fractal aggregates (tiny aggregate — fast)
# ---------------------------------------------------------------------------

class TestFractalAggregateRoundTrip:
    def test_save_list_load(self, nx_file):
        from pyirena.core.fractals import grow_aggregate, GrowthConfig
        from pyirena.io.nxcansas_fractals import (
            save_fractal_aggregate, list_fractal_aggregates, load_fractal_aggregate,
        )
        agg = grow_aggregate(GrowthConfig(z=30, seed=42))
        group_path = save_fractal_aggregate(nx_file, agg)

        listed = list_fractal_aggregates(nx_file)
        assert any(entry["group_path"] == group_path for entry in listed)

        loaded = load_fractal_aggregate(nx_file, group_path)
        np.testing.assert_array_equal(loaded.positions, agg.positions)
        assert loaded.uuid == agg.uuid


# ---------------------------------------------------------------------------
# Data manipulation / merge writers
# ---------------------------------------------------------------------------

class TestSaveManipulatedData:
    def test_nxcansas_source_round_trip(self, tmp_path, nx_file, q_i_err):
        from pyirena.io.nxcansas_data_manipulation import save_manipulated_data
        q, intensity, error = q_i_err
        out_dir = tmp_path / "out"
        scaled = intensity * 2.0
        out = save_manipulated_data(
            out_dir, nx_file, True, q, scaled, error, None,
            operation="scaled", provenance={"factor": 2.0},
        )
        assert out.exists()
        assert out.name == "sample_scaled.h5"
        with h5py.File(out, "r") as f:
            np.testing.assert_allclose(_sasdata(f)["I"][:], scaled)
            prov = f["entry/data_manipulation_results"]
            assert prov.attrs["NX_class"] == "NXprocess"
            assert float(prov["factor"][()]) == pytest.approx(2.0)

    def test_strips_nonpositive_points(self, tmp_path, nx_file, q_i_err):
        from pyirena.io.nxcansas_data_manipulation import save_manipulated_data
        q, intensity, error = q_i_err
        bad = intensity.copy()
        bad[:5] = -1.0          # subtraction artifacts
        bad[10] = np.nan
        out = save_manipulated_data(
            tmp_path / "out", nx_file, True, q, bad, error, None,
            operation="sub", provenance={},
        )
        with h5py.File(out, "r") as f:
            I = _sasdata(f)["I"][:]
            assert len(I) == len(q) - 6
            assert (I > 0).all()

    def test_refuses_degenerate_output(self, tmp_path, nx_file, q_i_err):
        from pyirena.io.nxcansas_data_manipulation import save_manipulated_data
        q, intensity, error = q_i_err
        with pytest.raises(ValueError):
            save_manipulated_data(
                tmp_path / "out", nx_file, True, q, -np.abs(intensity), error,
                None, operation="sub", provenance={},
            )


class TestSaveMergedData:
    def test_round_trip_with_provenance(self, tmp_path, nx_file, q_i_err):
        from pyirena.io.nxcansas_data_merge import save_merged_data
        q, intensity, error = q_i_err
        merged_I = intensity * 1.5
        out = save_merged_data(
            tmp_path / "out", nx_file, True, q, merged_I, error, None,
            merge_result_dict={"scale": 1.5, "q_shift": 0.0, "background": 0.0,
                               "chi_squared": 1.0, "q_overlap_min": 0.01,
                               "q_overlap_max": 0.02, "scale_dataset": "DS2",
                               "fit_scale": True, "fit_qshift": False,
                               "split_at_left_cursor": False},
        )
        assert out.name == "sample_merged.h5"
        with h5py.File(out, "r") as f:
            np.testing.assert_allclose(_sasdata(f)["I"][:], merged_I)
            prov = f["entry/data_merge_results"]
            assert float(prov["scale"][()]) == pytest.approx(1.5)

    def test_existing_results_stripped_from_copy(self, tmp_path, nx_file, q_i_err):
        """Old fit results must NOT survive into the merged output file."""
        from pyirena.io.nxcansas_data_merge import save_merged_data
        q, intensity, error = q_i_err
        model = intensity * 1.01
        save_unified_fit_results(
            nx_file, q, intensity, model, intensity - model,
            levels=_unified_levels()[:1], background=0.0,
            chi_squared=1.0, num_levels=1, error=error,
        )
        out = save_merged_data(
            tmp_path / "out", nx_file, True, q, intensity, error, None,
            merge_result_dict={"scale": 1.0},
        )
        with h5py.File(out, "r") as f:
            assert "entry/unified_fit_results" not in f
