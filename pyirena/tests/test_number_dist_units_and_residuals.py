"""Number distributions are stored in Irena units, Modeling residuals are saved.

Number distributions are particles per cm³ of sample per Å of radius
[1/(cm³·Å)], as in Irena.  Files written before v1.2.0 held the same quantity
per Å³ of sample (1e24 smaller); the loaders must convert them so every reader
sees one unit regardless of file age.
"""
import h5py
import numpy as np

from pyirena.core.modeling import ModelingConfig, ModelingEngine, SizeDistPopulation
from pyirena.core.sizes import number_distribution
from pyirena.io._nxcansas_common import NUMBER_DIST_UNITS
from pyirena.io.nxcansas_modeling import load_modeling_results, save_modeling_results
from pyirena.io.nxcansas_sizes import load_sizes_results, save_sizes_results


def _modeling_fit():
    q = np.geomspace(1e-3, 0.3, 120)
    pop = SizeDistPopulation()
    pop.enabled = True
    cfg = ModelingConfig(populations=[pop], background=0.01)
    eng = ModelingEngine()
    I_true, _, _, _ = eng.total_intensity_maybe_smeared(cfg, q, use_cache=True)
    rng = np.random.default_rng(0)
    dI = I_true * 0.02
    I_data = I_true + rng.normal(0.0, 1.0, q.size) * dI
    return eng.fit(cfg, q, I_data, dI), q, I_data, dI


def test_modeling_residuals_roundtrip(tmp_path):
    res, q, I_data, dI = _modeling_fit()
    assert res.residuals is not None
    np.testing.assert_allclose(res.residuals, (I_data - res.model_I) / dI)

    out = tmp_path / "mod.h5"
    save_modeling_results(out, res)
    back = load_modeling_results(out)
    np.testing.assert_allclose(back["residuals"], res.residuals)
    with h5py.File(out, "r") as f:
        assert f["entry/modeling_results/residuals"].attrs["units"] == "dimensionless"



def test_modeling_residuals_reach_the_agent_reader(tmp_path):
    """pyirena.api.read_modeling (and pyirena_read_modeling over MCP) returns them."""
    from pyirena.api import read_modeling

    res, *_ = _modeling_fit()
    out = tmp_path / "mod.h5"
    save_modeling_results(out, res)

    assert read_modeling(str(out))["residuals"] is None   # arrays are opt-in
    got = read_modeling(str(out), include_arrays=True)
    np.testing.assert_allclose(got["residuals"], res.residuals)
    assert len(got["residuals"]) == len(got["model_q"])


def test_modeling_number_dist_is_per_cm3(tmp_path):
    res, *_ = _modeling_fit()
    r = res.radius_grids[0]
    V_cm3 = (4.0 / 3.0) * np.pi * r**3 * 1e-24
    np.testing.assert_allclose(res.number_dists[0], res.volume_dists[0] / V_cm3, rtol=1e-12)

    out = tmp_path / "mod.h5"
    save_modeling_results(out, res)
    with h5py.File(out, "r") as f:
        ds = f["entry/modeling_results/pop_01/number_dist"]
        assert ds.attrs["units"] == NUMBER_DIST_UNITS
    back = load_modeling_results(out)
    np.testing.assert_allclose(back["populations"][0]["number_dist"], res.number_dists[0])


def test_legacy_modeling_number_dist_is_rescaled(tmp_path):
    """Pre-1.2.0 files: number_dist per Å³ of sample, no units attribute."""
    res, *_ = _modeling_fit()
    out = tmp_path / "mod.h5"
    save_modeling_results(out, res)
    with h5py.File(out, "a") as f:
        pg = f["entry/modeling_results/pop_01"]
        del pg["number_dist"]
        pg.create_dataset("number_dist", data=res.number_dists[0] * 1e-24)
        del f["entry/modeling_results/residuals"]
    back = load_modeling_results(out)
    np.testing.assert_allclose(back["populations"][0]["number_dist"], res.number_dists[0])
    assert back["residuals"] is None


def _sizes_file(tmp_path, r, dist):
    q = np.geomspace(1e-3, 0.3, 50)
    out = tmp_path / "sz.h5"
    with h5py.File(out, "w") as f:
        f.create_group("entry")
    save_sizes_results(out, q, q, q, np.zeros_like(q), r, dist,
                       params={"shape": "sphere", "aspect_ratio": 1.0})
    return out


def test_sizes_number_dist_units_and_legacy(tmp_path):
    r = np.linspace(10.0, 500.0, 60)
    dist = np.exp(-(((r - 150.0) / 40.0) ** 2)) * 1e-4
    nd, cnd = number_distribution(dist, r)
    np.testing.assert_allclose(nd, dist / ((4.0 / 3.0) * np.pi * r**3 * 1e-24), rtol=1e-12)

    out = _sizes_file(tmp_path, r, dist)
    back = load_sizes_results(out)
    np.testing.assert_allclose(back["number_dist"], nd)
    np.testing.assert_allclose(back["cumul_num_dist"], cnd)
    with h5py.File(out, "r") as f:
        g = f["entry/sizes_results"]
        assert g["surface_dist"].attrs["units"] == "1/angstrom^2"
        assert g["cumul_surf_dist"].attrs["units"] == "1/angstrom"

    # Rewrite as a pre-1.2.0 file would have it: per Å³, mislabelled units.
    with h5py.File(out, "a") as f:
        g = f["entry/sizes_results"]
        for name, data, units in (("number_dist", nd * 1e-24, "1/angstrom"),
                                  ("cumul_num_dist", cnd * 1e-24, "dimensionless")):
            del g[name]
            g.create_dataset(name, data=data)
            g[name].attrs["units"] = units
    back = load_sizes_results(out)
    np.testing.assert_allclose(back["number_dist"], nd)
    np.testing.assert_allclose(back["cumul_num_dist"], cnd)
