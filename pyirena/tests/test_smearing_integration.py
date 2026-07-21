"""
Integration tests for slit-smearing across pyIrena tools (plan §7 / F9).

These lock in the contracts that the code review (F1–F10) and Jan's testing
exposed:

* Simple Fits **display** smears the model (compute()), not just the fit
  objective — the sphere/spheroid oscillations must visibly damp (A1).
* Background prefits are ideal-space (no double-smearing) when smearing is on
  (F1).
* Data Merge propagates slit provenance: a smeared+pinhole merge reads back as
  slit smeared, and mismatched nonzero slit lengths warn (F3).
* Data Manipulation refuses to mix incompatible smearing status, and its output
  carries no stale ``_SMR`` twin (F2).
* Modeling / Simple Fits / Sizes save→load preserve ``slit_length`` and the
  ideal curve (F4 + parity).
"""

import os
from pathlib import Path

import numpy as np
import pytest

from pyirena.core.smearing import SlitSmearer


# --------------------------------------------------------------------------- #
# A1 — Simple Fits display smears the model
# --------------------------------------------------------------------------- #
def test_simple_fits_compute_smears_sphere():
    from pyirena.core.simple_fits import SimpleFitModel

    q = np.linspace(0.002, 0.03, 400)
    m = SimpleFitModel()
    m.set_model("Sphere")
    m.params["Scale"] = 1.0
    m.params["R"] = 600.0

    m.use_slit_smearing = False
    I_ideal = m.compute(q)

    m.use_slit_smearing = True
    m.slit_length = 0.018
    I_smeared = m.compute(q)

    # Oscillation contrast (std of log I) must drop once the model is smeared.
    def _contrast(I):
        I = I[np.isfinite(I)]
        return float(np.nanstd(np.log(I)))

    assert _contrast(I_smeared) < _contrast(I_ideal)

    # The first form-factor minimum must fill in (dip becomes shallower).
    dip = (q > 0.006) & (q < 0.009)
    assert np.nanmin(I_smeared[dip]) > np.nanmin(I_ideal[dip])


def test_simple_fits_compute_noop_when_pinhole():
    from pyirena.core.simple_fits import SimpleFitModel

    q = np.linspace(0.002, 0.03, 200)
    m = SimpleFitModel()
    m.set_model("Sphere")
    m.params["R"] = 300.0
    m.use_slit_smearing = False
    a = m.compute(q)
    m.use_slit_smearing = True
    m.slit_length = 0.0            # SL<=0 => strict no-op even if flag set
    b = m.compute(q)
    np.testing.assert_array_equal(a[np.isfinite(a)], b[np.isfinite(b)])


# --------------------------------------------------------------------------- #
# F1 — background prefit is ideal-space (no double-smearing)
# --------------------------------------------------------------------------- #
def test_prefit_power_law_recovers_ideal_bp_on_smeared_data():
    from pyirena.core.saxs_morph import fit_power_law_bg, fit_power_law_bg_fixed_p

    B_true, P_true, SL = 5e-4, 3.6, 0.018
    q = np.geomspace(0.02, 0.2, 200)
    sm = SlitSmearer(q, SL)
    I_smeared = sm.smear_model(lambda x: B_true * x ** -P_true)

    # Pinhole (SL=0) fit to smeared data is biased.
    B_pin, P_pin = fit_power_law_bg(q, I_smeared, q.min(), q.max(), slit_length=0.0)
    assert abs(P_pin - P_true) > 0.05

    # Smear-aware prefit recovers the ideal power law.
    B_sm, P_sm = fit_power_law_bg(q, I_smeared, q.min(), q.max(), slit_length=SL)
    assert abs(P_sm - P_true) / P_true < 0.01
    assert abs(B_sm - B_true) / B_true < 0.02

    # Fixed-P prefactor recovery is exact (smearing is linear in B).
    B_fx = fit_power_law_bg_fixed_p(q, I_smeared, q.min(), q.max(), P_true, slit_length=SL)
    assert abs(B_fx - B_true) / B_true < 0.01


# --------------------------------------------------------------------------- #
# F3 — Data Merge slit provenance
# --------------------------------------------------------------------------- #
def test_merge_merged_slit_length_rules():
    from pyirena.core.data_merge import DataMerge

    assert DataMerge.merged_slit_length(0.018, 0.0)[:2] == (0.018, True)
    assert DataMerge.merged_slit_length(0.0, 0.02)[:2] == (0.02, True)
    assert DataMerge.merged_slit_length(0.0, 0.0)[:2] == (0.0, False)
    sl, smeared, warn = DataMerge.merged_slit_length(0.018, 0.02)
    assert smeared and sl == 0.02 and warn        # mismatch keeps larger + warns


def test_merge_output_reads_back_slit_smeared(tmp_path):
    """A slit-smeared + pinhole merge writes dQl so the output auto-detects."""
    from pyirena.io.nxcansas_unified import create_nxcansas_file
    from pyirena.io.nxcansas_data_merge import save_merged_data
    from pyirena.io.hdf5 import readGenericNXcanSAS

    # DS1: slit-smeared USAXS (low Q).  Create a plain NXcanSAS file; the merge
    # writes the merged dQl regardless of DS1's own resolutions.
    q = np.geomspace(1e-3, 0.3, 120)
    I = q ** -3.0 + 0.01
    ds1 = tmp_path / "ds1.h5"
    create_nxcansas_file(ds1, q, I, error=I * 0.03, sample_name="ds1")

    prov = {
        "scale": 1.0, "background": 0.0, "chi_squared": 1.0,
        "slit_length_ds1": 0.018, "slit_length_ds2": 0.0,
        "slit_length_merged": 0.018, "is_slit_smeared_merged": True,
    }
    out = save_merged_data(
        output_folder=tmp_path, ds1_path=ds1, ds1_is_nxcansas=True,
        q=q, I=I, dI=I * 0.03, dQ=None, merge_result_dict=prov,
        ds2_path=None,
    )
    back = readGenericNXcanSAS(str(out.parent), out.name)
    assert back["is_slit_smeared"] is True
    assert abs(back["slit_length"] - 0.018) < 1e-9


# --------------------------------------------------------------------------- #
# F2 — Data Manipulation guard + stale _SMR twin
# --------------------------------------------------------------------------- #
def test_manipulation_subtract_guard():
    from pyirena.core.data_manipulation import DataManipulation, SubtractConfig

    q = np.geomspace(1e-3, 0.3, 80)
    I = q ** -3.0 + 1.0
    dI = I * 0.03

    # Matching smeared operands: OK, result carries slit provenance.
    r = DataManipulation.subtract(
        q, I, dI, None, q, 0.5 * I, 0.5 * dI, SubtractConfig(1.0),
        is_slit_smeared_sample=True, slit_length_sample=0.018,
        is_slit_smeared_buffer=True, slit_length_buffer=0.018,
    )
    assert r.metadata["is_slit_smeared"] is True
    assert abs(r.metadata["slit_length"] - 0.018) < 1e-9

    # Mixed smearing status: refuse.
    with pytest.raises(ValueError):
        DataManipulation.subtract(
            q, I, dI, None, q, 0.5 * I, 0.5 * dI, SubtractConfig(1.0),
            is_slit_smeared_sample=True, slit_length_sample=0.018,
            is_slit_smeared_buffer=False, slit_length_buffer=0.0,
        )

    # Different slit lengths: refuse.
    with pytest.raises(ValueError):
        DataManipulation.divide(
            q, I, dI, None, q, 0.5 * I, 0.5 * dI,
            __import__("pyirena.core.data_manipulation", fromlist=["DivideConfig"]).DivideConfig(1.0),
            is_slit_smeared_num=True, slit_length_num=0.018,
            is_slit_smeared_den=True, slit_length_den=0.03,
        )


def test_drop_smr_entries(tmp_path):
    """A copied Matilda-style file loses its _SMR twin on manipulation save."""
    import h5py
    from pyirena.io.nxcansas_unified import create_nxcansas_file
    from pyirena.io._nxcansas_common import drop_smr_entries

    p = tmp_path / "dual.h5"
    q = np.geomspace(1e-3, 0.3, 60)
    create_nxcansas_file(p, q, q ** -3.0, error=q ** -3.0 * 0.03, sample_name="s")
    # Add a fake _SMR sibling group under entry.
    with h5py.File(p, "a") as f:
        entry = f["entry"] if "entry" in f else f
        g = entry.create_group("s_SMR")
        g.attrs["NX_class"] = "NXsubentry"
    with h5py.File(p, "r") as f:
        assert any(k.endswith("_SMR") for k in f["entry"].keys())

    n = drop_smr_entries(p)
    assert n == 1
    with h5py.File(p, "r") as f:
        assert not any(k.endswith("_SMR") for k in f["entry"].keys())


# --------------------------------------------------------------------------- #
# F4 + parity — Modeling / Simple Fits save→load preserve slit provenance
# --------------------------------------------------------------------------- #
def test_simple_fits_roundtrip_slit(tmp_path):
    from pyirena.core.simple_fits import SimpleFitModel
    from pyirena.io.nxcansas_simple_fits import (
        save_simple_fit_results, load_simple_fit_results,
    )

    SL = 0.018
    q = np.geomspace(0.005, 0.2, 150)
    m = SimpleFitModel()
    m.set_model("Guinier")
    m.params.update({"I0": 100.0, "Rg": 40.0})
    m.use_slit_smearing = True
    m.slit_length = SL
    # Self-consistent smeared data so the fit succeeds.
    sm = SlitSmearer(q, SL)
    I = sm.smear_model(lambda x: 100.0 * np.exp(-x ** 2 * 40.0 ** 2 / 3.0))
    res = m.fit(q, I, I * 0.02)
    assert res["success"]
    assert res["slit_length"] == SL
    assert res["I_model_ideal"] is not None

    out = tmp_path / "sf.h5"
    save_simple_fit_results(out, res)
    back = load_simple_fit_results(out)
    assert abs(back["slit_length"] - SL) < 1e-9
    assert back["data_is_slit_smeared"] is True


def test_modeling_saves_ideal_curve(tmp_path):
    """Modeling result exposes model_I_ideal when smearing is on; it round-trips."""
    from pyirena.core.modeling import (
        ModelingEngine, ModelingConfig, SizeDistPopulation,
    )
    from pyirena.io.nxcansas_modeling import (
        save_modeling_results, load_modeling_results,
    )

    q = np.geomspace(1e-3, 0.3, 120)
    # One sphere size-distribution population, evaluate-only (nothing to fit).
    pop = SizeDistPopulation()
    pop.enabled = True
    cfg = ModelingConfig(populations=[pop], background=0.01)
    cfg.use_slit_smearing = True
    cfg.slit_length = 0.018

    eng = ModelingEngine()
    # Self-consistent data: smeared model + noise-free so the eval path runs.
    I_data, _, _, _ = eng.total_intensity_maybe_smeared(cfg, q, use_cache=True)
    res = eng.fit(cfg, q, I_data, I_data * 0.02)
    assert res.model_I_ideal is not None
    assert res.model_I_ideal.shape == res.model_I.shape

    out = tmp_path / "mod.h5"
    save_modeling_results(out, res)
    back = load_modeling_results(out)
    assert back is not None
    assert back["model_I_ideal"] is not None
    assert abs(back["slit_length"] - 0.018) < 1e-9
    assert back["data_is_slit_smeared"] is True
