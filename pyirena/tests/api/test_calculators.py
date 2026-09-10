"""Tests for pyirena.api.calculators — the stateless calculators group.

No fixtures: these functions take no dataset and touch no user paths. The
optional ``periodictable`` / ``xraydb`` dependencies (the ``contrast``
extra) are skipped per-test so the suite still runs on a base install.
"""
from __future__ import annotations

import json

import pytest

from pyirena.api import (
    calc_compound,
    calc_contrast,
    calc_contrast_energy_scan,
    list_compound_library,
    load_compound,
    lookup_element,
)
from pyirena.api.calculator_schemas import (
    CALCULATOR_SCHEMA_BY_NAME,
    CALCULATOR_TOOL_SCHEMAS,
)
from pyirena.api.calculators import COMPOSITION_MODES

pt = pytest.importorskip("periodictable")


def _assert_json_safe(result):
    """Every calculator return must survive strict JSON (the MCP wire)."""
    json.dumps(result)


# ---------------------------------------------------------------------------
# calc_contrast
# ---------------------------------------------------------------------------

def test_calc_contrast_tio_ti2o3():
    """The originating question: X-ray contrast between TiO and Ti2O3."""
    r = calc_contrast("TiO", 4.95, "Ti2O3", 4.49)
    assert "error" not in r
    # SLDs of the two titanium oxides, 10^10 cm^-2
    assert r["compound_1"]["xray_sld"] == pytest.approx(39.46, rel=0.01)
    assert r["compound_2"]["xray_sld"] == pytest.approx(36.05, rel=0.01)
    # (39.46 - 36.05)^2 ~ 11.6, in 10^20 cm^-4 — directly usable as the
    # 'contrast' parameter of a Sizes or Modeling fit.
    assert r["xray_contrast"] == pytest.approx(11.63, rel=0.01)
    assert r["neutron_contrast"] == pytest.approx(0.76, rel=0.02)
    assert r["units"]["xray_contrast"] == "10^20 cm^-4"
    _assert_json_safe(r)


def test_calc_contrast_silica_water_matches_core():
    """Same value the core-level test locks in (SiO2 vs H2O ~ 87)."""
    r = calc_contrast("SiO2", 2.2, "H2O", 1.0)
    assert r["xray_contrast"] == pytest.approx(87.0, rel=0.05)


def test_calc_contrast_is_symmetric():
    a = calc_contrast("SiO2", 2.2, "H2O", 1.0)
    b = calc_contrast("H2O", 1.0, "SiO2", 2.2)
    assert a["xray_contrast"] == pytest.approx(b["xray_contrast"])
    assert a["neutron_contrast"] == pytest.approx(b["neutron_contrast"])


@pytest.mark.parametrize("formula,density", [("", 0.0), ("", 1.0), ("SiO2", 0.0)])
def test_vacuum_from_empty_formula_or_zero_density(formula, density):
    """Panel rule: empty formula OR zero density means vacuum."""
    r = calc_contrast("SiO2", 2.2, formula, density)
    assert r["compound_2"]["formula_str"] == "vacuum"
    assert r["compound_2"]["xray_sld"] == 0.0
    assert r["compound_2"]["density"] == 0.0


def test_vacuum_does_not_leak_the_core_singleton():
    """A returned vacuum must not be the shared mutable VACUUM object."""
    from pyirena.core.scattering_contrast import VACUUM

    before = VACUUM.name
    calc_contrast("SiO2", 2.2, "", 0.0, name_2="my vacuum")
    assert VACUUM.name == before


def test_isotope_override_gives_d2o_neutron_contrast():
    light = calc_contrast("SiO2", 2.2, "H2O", 1.0)
    heavy = calc_contrast("SiO2", 2.2, "H2O", 1.107, isotopes_2={"H": "2"})
    # D2O neutron SLD is +6.37 vs -0.56 for H2O, so the neutron contrast
    # against silica changes markedly while the X-ray contrast barely moves.
    assert heavy["compound_2"]["neutron_sld"] == pytest.approx(6.37, rel=0.02)
    assert heavy["neutron_contrast"] != pytest.approx(light["neutron_contrast"])


def test_calc_contrast_anomalous_adds_absorption_fields():
    pytest.importorskip("xraydb")
    r = calc_contrast("TiO", 4.95, "Ti2O3", 4.49, energy_keV=12.0, thickness_mm=0.5)
    for key in (
        "xray_contrast_anom", "xray_sld_anom_1", "xray_sld_anom_2",
        "mu_1", "mu_2", "transmission_1", "transmission_2",
        "transmission_sample",
    ):
        assert key in r, key
    assert r["energy_keV"] == 12.0
    assert r["thickness_mm"] == 0.5
    assert 0.0 <= r["transmission_sample"] <= 1.0
    _assert_json_safe(r)


def test_calc_contrast_without_energy_omits_anomalous():
    r = calc_contrast("TiO", 4.95, "Ti2O3", 4.49)
    assert "xray_contrast_anom" not in r
    assert "energy_keV" not in r


# ---------------------------------------------------------------------------
# calc_compound
# ---------------------------------------------------------------------------

def test_calc_compound_water():
    r = calc_compound("H2O", 1.0, name="Water")
    assert r["name"] == "Water"
    assert r["mol_weight"] == pytest.approx(18.015, abs=0.02)
    assert r["xray_sld"] == pytest.approx(9.47, rel=0.01)
    assert r["neutron_sld"] == pytest.approx(-0.56, abs=0.03)
    assert "_element_counts" not in r  # core-private plumbing stays internal
    _assert_json_safe(r)


def test_calc_compound_anomalous_block():
    pytest.importorskip("xraydb")
    r = calc_compound("Fe", 7.87, energy_keV=7.10)
    assert r["anomalous"]["energy_keV"] == 7.10
    assert r["anomalous"]["mu_linear"] > 0
    # f' is strongly negative near the Fe K edge (7.112 keV)
    far = calc_compound("Fe", 7.87, energy_keV=12.0)
    assert r["anomalous"]["xray_sld_anom"] < far["anomalous"]["xray_sld_anom"]


@pytest.mark.parametrize(
    "formula,mode",
    [
        ("Au0.35Ag0.65", "weight_fraction_elements"),
        ("Y2O3:0.10 ZrO2:0.90", "weight_fraction_compounds"),
    ],
)
def test_weight_fraction_modes_warn_about_basis(formula, mode):
    r = calc_compound(formula, 6.0, mode=mode)
    assert "error" not in r
    assert r["xray_sld"] > 0
    assert "basis_warning" in r


def test_atomic_ratio_mode_has_no_basis_warning():
    assert "basis_warning" not in calc_compound("SiO2", 2.2)


# ---------------------------------------------------------------------------
# lookup_element
# ---------------------------------------------------------------------------

def test_lookup_element_titanium():
    r = lookup_element("Ti")
    assert r["Z"] == 22
    assert r["mass"] == pytest.approx(47.867, abs=0.01)
    assert r["neutron_b_c"] is not None
    labels = [iso["label"] for iso in r["isotopes"]]
    assert "natural" in labels
    _assert_json_safe(r)


def test_lookup_element_unknown_symbol_returns_error():
    """Guards the core's get_element_info NameError on an unknown symbol."""
    r = lookup_element("Zz")
    assert r["code"] == "BAD_ELEMENT"
    assert "Zz" in r["error"]


def test_lookup_element_empty_symbol_returns_error():
    assert lookup_element("")["code"] == "BAD_ELEMENT"


# ---------------------------------------------------------------------------
# calc_contrast_energy_scan
# ---------------------------------------------------------------------------

def test_energy_scan_returns_arrays_and_best():
    pytest.importorskip("xraydb")
    r = calc_contrast_energy_scan(
        "TiO", 4.95, "Ti2O3", 4.49, e_start_keV=4.5, e_end_keV=5.5, n_points=20
    )
    assert "error" not in r
    assert len(r["energy"]) == 20
    assert len(r["xray_contrast_anom"]) == 20
    assert 4.5 <= r["best"]["energy_keV"] <= 5.5
    _assert_json_safe(r)


def test_energy_scan_decimates_but_keeps_best_full_resolution():
    pytest.importorskip("xraydb")
    r = calc_contrast_energy_scan(
        "TiO", 4.95, "Ti2O3", 4.49,
        e_start_keV=4.5, e_end_keV=5.5, n_points=60, max_points=10,
    )
    assert len(r["energy"]) == 10
    assert r["n_points"] == 60
    assert 4.5 <= r["best"]["energy_keV"] <= 5.5


def test_energy_scan_rejects_reversed_range():
    """Validation that previously existed only in the Qt panel."""
    r = calc_contrast_energy_scan("TiO", 4.95, "Ti2O3", 4.49, 5.5, 4.5)
    assert r["code"] == "BAD_ENERGY_RANGE"


@pytest.mark.parametrize("n_points", [1, 5000])
def test_energy_scan_rejects_bad_n_points(n_points):
    r = calc_contrast_energy_scan(
        "TiO", 4.95, "Ti2O3", 4.49, 4.5, 5.5, n_points=n_points
    )
    assert r["code"] == "BAD_N_POINTS"


# ---------------------------------------------------------------------------
# Error paths
# ---------------------------------------------------------------------------

def test_bad_mode_lists_valid_modes():
    r = calc_contrast("TiO", 4.95, "Ti2O3", 4.49, mode="bogus")
    assert r["code"] == "BAD_MODE"
    for mode in COMPOSITION_MODES:
        assert mode in r["suggestion"]


def test_bad_formula_suggestion_is_mode_specific():
    r = calc_compound("H:0.1 O:0.9", 1.0, mode="weight_fraction_elements")
    assert r["code"] == "BAD_FORMULA"
    assert "Au0.35Ag0.65" in r["suggestion"]


def test_negative_density_returns_error():
    assert calc_compound("SiO2", -2.2)["code"] == "BAD_DENSITY"


def test_errors_are_json_safe():
    _assert_json_safe(calc_contrast("Zzz9", 1.0, "H2O", 1.0))


# ---------------------------------------------------------------------------
# Compound library (read-only)
# ---------------------------------------------------------------------------

def test_list_compound_library_shape():
    r = list_compound_library()
    assert isinstance(r["compounds"], list)
    assert r["library_path"].endswith(".h5")
    _assert_json_safe(r)


def test_load_missing_compound_returns_error():
    r = load_compound("definitely-not-a-saved-compound")
    assert r["code"] in ("UNKNOWN_COMPOUND", "NO_LIBRARY")


def test_library_is_read_only():
    """No save/delete is exposed through the api surface."""
    import pyirena.api.calculators as calculators

    for forbidden in ("save_compound", "delete_compound"):
        assert not hasattr(calculators, forbidden)


# ---------------------------------------------------------------------------
# Schema parity (mirrors test_control_schemas.py for the calculators registry)
# ---------------------------------------------------------------------------

def _sig_params(fn) -> dict[str, bool]:
    import inspect

    return {
        name: p.default is not inspect.Parameter.empty
        for name, p in inspect.signature(fn).parameters.items()
    }


@pytest.mark.parametrize(
    "schema", CALCULATOR_TOOL_SCHEMAS, ids=lambda s: s["name"]
)
def test_calculator_schema_matches_callable_signature(schema):
    import pyirena.api.calculators as calculators

    fn = getattr(calculators, schema["name"], None)
    assert callable(fn), f"schema '{schema['name']}' has no calculators callable"

    params = _sig_params(fn)
    mandatory = {n for n, has_default in params.items() if not has_default}

    ins = schema["input_schema"]
    props = set(ins.get("properties", {}))
    required = set(ins.get("required", []))

    assert ins["type"] == "object"
    assert schema["description"].strip()
    assert props <= set(params), \
        f"{schema['name']}: schema exposes non-parameters {props - set(params)}"
    assert required == mandatory, \
        f"{schema['name']}: required {required} != mandatory params {mandatory}"
    assert set(params) <= props, \
        f"{schema['name']}: parameters missing from schema {set(params) - props}"


def test_every_calculator_is_exported_from_the_api_facade():
    import pyirena.api as papi

    for name in CALCULATOR_SCHEMA_BY_NAME:
        assert name in papi.__all__, f"{name} missing from pyirena.api.__all__"
        assert callable(getattr(papi, name))
