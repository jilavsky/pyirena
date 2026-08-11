"""
``to_dict`` / ``from_dict`` on the core model objects.

The golden rule is that all state lives in the core model object and travels
``dict → JSON → HDF5`` unchanged.  Two things have to hold for that to be true,
and both are tested here for every tool that has the pair:

* **A round trip changes nothing.**  Through ``json.dumps``/``loads``, not just
  in memory — that is where tuples silently become lists and numpy scalars
  refuse to serialise at all.
* **An old file still loads.**  ``from_dict`` supplies a default for every
  field, so a dict written before a field existed opens with that field at its
  default rather than raising.

The Unified tests additionally pin the *panel* vocabulary (``RgCutoff``,
``correlated``, ``estimate_B``, ``Rg_low``…), which is baked into saved setups
and batch config files and therefore cannot be renamed, only translated.
"""

from __future__ import annotations

import json

import numpy as np
import pytest

from pyirena.core.unified import UnifiedFitModel, UnifiedLevel


def _through_json(d: dict) -> dict:
    """What the dict looks like after a trip through a file."""
    return json.loads(json.dumps(d))


# ── Unified: the model ───────────────────────────────────────────────────

def _populated_model() -> UnifiedFitModel:
    model = UnifiedFitModel(num_levels=2)
    model.levels[0].Rg = 55.0
    model.levels[0].G = 12.5
    model.levels[0].fit_P = True
    model.levels[0].P_limits = (2.0, 4.5)
    model.levels[1].correlations = True
    model.levels[1].PACK = 3.0
    model.levels[1].link_B = True
    model.background = 0.017
    model.fit_background = False
    model.use_slit_smearing = True
    model.slit_length = 0.031
    return model


def test_a_unified_model_survives_a_json_round_trip():
    model = _populated_model()
    restored = UnifiedFitModel.from_dict(_through_json(model.to_dict()))
    assert restored.to_dict() == model.to_dict()


def test_the_round_trip_keeps_the_fit_it_describes():
    """Equal dicts are necessary but not sufficient — the physics must match."""
    model = _populated_model()
    restored = UnifiedFitModel.from_dict(_through_json(model.to_dict()))
    q = np.logspace(-3, -1, 60)
    assert np.allclose(model.calculate_intensity(q), restored.calculate_intensity(q))


def test_bounds_come_back_as_tuples_not_lists():
    """The fitter unpacks bounds; a list from JSON must not leak through."""
    restored = UnifiedFitModel.from_dict(_through_json(_populated_model().to_dict()))
    assert isinstance(restored.levels[0].P_limits, tuple)
    assert restored.levels[0].P_limits == (2.0, 4.5)


def test_data_and_results_are_not_part_of_the_state():
    """to_dict is settings, not outcomes — the same convention as Sizes."""
    model = _populated_model()
    model.q_data = np.array([1.0, 2.0])
    model.chi_squared = 1.23
    d = model.to_dict()
    assert "q_data" not in d and "chi_squared" not in d and "fit_result" not in d
    assert json.dumps(d)          # and so it is serialisable at all


# ── Unified: files written by older versions ─────────────────────────────

def test_a_dict_with_only_levels_still_loads():
    model = UnifiedFitModel.from_dict({"levels": [{"Rg": 7.0}, {"Rg": 70.0}]})
    assert model.num_levels == 2
    assert [lv.Rg for lv in model.levels] == [7.0, 70.0]
    assert model.background == 0.0            # default, not a crash


def test_fields_a_previous_version_never_wrote_take_their_defaults():
    old = {"Rg": 40.0, "G": 2.0, "P": 4.0, "B": 1e-5}    # no flags, no bounds
    level = UnifiedLevel.from_dict(old)
    fresh = UnifiedLevel()
    assert level.Rg == 40.0
    assert level.mass_fractal is fresh.mass_fractal
    assert level.fit_RgCO is fresh.fit_RgCO
    assert level.P_limits == fresh.P_limits


@pytest.mark.parametrize("junk", [None, {}, {"levels": []}])
def test_an_empty_or_missing_dict_gives_a_usable_model(junk):
    model = UnifiedFitModel.from_dict(junk)
    assert model.num_levels == 1 and len(model.levels) == 1


def test_a_level_count_that_disagrees_with_the_levels_follows_the_levels():
    model = UnifiedFitModel.from_dict({"num_levels": 5, "levels": [{"Rg": 3.0}]})
    assert model.num_levels == 1 and len(model.levels) == 1


# ── Unified: the panel vocabulary ────────────────────────────────────────

PANEL_PARAMS = {
    "G": 12.0, "Rg": 45.0, "B": 2e-5, "P": 3.8, "RgCutoff": 5.0,
    "ETA": 20.0, "PACK": 2.0,
    "fit_G": True, "fit_Rg": True, "fit_B": False, "fit_P": True,
    "fit_ETA": True, "fit_PACK": False,
    "estimate_B": True, "correlated": True, "link_rgco": True,
    "G_low": 1.0, "G_high": 100.0, "Rg_low": 10.0, "Rg_high": 90.0,
    "B_low": 1e-9, "B_high": 1e-1, "P_low": 3.0, "P_high": 4.2,
    "ETA_low": 1.0, "ETA_high": 60.0, "PACK_low": 0.0, "PACK_high": 8.0,
}


def test_the_panel_key_names_map_onto_the_model():
    level = UnifiedLevel.from_panel_params(PANEL_PARAMS)
    assert level.RgCO == 5.0                  # RgCutoff → RgCO
    assert level.correlations is True         # correlated → correlations
    assert level.link_B is True               # estimate_B → link_B
    assert level.link_RGCO is True            # link_rgco → link_RGCO
    assert level.Rg_limits == (10.0, 90.0)    # Rg_low/Rg_high → Rg_limits
    assert level.fit_B is False


def test_a_caller_can_refuse_the_links_because_they_change_the_intensity():
    """Drawing the curve the user typed is not the same as fitting it.

    ``link_B`` recomputes B from G, Rg and P at every evaluation, so the graph
    and undo paths pass ``with_links=False`` to plot the typed value.
    """
    linked = UnifiedLevel.from_panel_params(PANEL_PARAMS)
    plain = UnifiedLevel.from_panel_params(PANEL_PARAMS, with_links=False)
    assert linked.link_B is True and plain.link_B is False
    assert plain.B == PANEL_PARAMS["B"]


def test_no_limits_mode_leaves_the_wide_defaults_in_place():
    wide = UnifiedLevel.from_panel_params(PANEL_PARAMS, with_limits=False)
    assert wide.Rg_limits == UnifiedLevel().Rg_limits
    assert wide.G_limits == UnifiedLevel().G_limits
    assert wide.Rg == 45.0                     # values still come across


def test_from_panel_params_matches_the_construction_it_replaced():
    """Locks the mapping against the code that was written out by hand.

    Six call sites in the GUI and two in batch each built a UnifiedLevel from
    this dict; this is the constructor call they used, kept here so a change to
    the mapping has to be deliberate.
    """
    p = PANEL_PARAMS
    expected = UnifiedLevel(
        Rg=p["Rg"], G=p["G"], P=p["P"], B=p["B"],
        RgCO=p["RgCutoff"], ETA=p["ETA"], PACK=p["PACK"],
        correlations=p["correlated"],
        fit_Rg=p["fit_Rg"], fit_G=p["fit_G"], fit_P=p["fit_P"],
        fit_B=p["fit_B"], fit_ETA=p["fit_ETA"], fit_PACK=p["fit_PACK"],
        link_B=p["estimate_B"], link_RGCO=p["link_rgco"],
        Rg_limits=(p["Rg_low"], p["Rg_high"]),
        G_limits=(p["G_low"], p["G_high"]),
        P_limits=(p["P_low"], p["P_high"]),
        B_limits=(p["B_low"], p["B_high"]),
        ETA_limits=(p["ETA_low"], p["ETA_high"]),
        PACK_limits=(p["PACK_low"], p["PACK_high"]),
    )
    assert UnifiedLevel.from_panel_params(p).to_dict() == expected.to_dict()


def test_a_partial_panel_dict_does_not_raise():
    """Older saved setups are missing keys added since; they must still open."""
    level = UnifiedLevel.from_panel_params({"G": 1.0, "Rg": 20.0})
    assert level.Rg == 20.0 and level.RgCO == 0.0


def test_unreadable_numbers_fall_back_rather_than_crash():
    level = UnifiedLevel.from_panel_params({**PANEL_PARAMS, "Rg": "not a number"})
    assert level.Rg == 0.0 and level.G == 12.0


# ── Unified: the batch config shape ──────────────────────────────────────

def test_the_batch_config_flattens_to_the_panel_names():
    from pyirena.batch.unified import _flatten_level_config

    flat = _flatten_level_config({
        "G": {"value": 3.0, "fit": True, "low_limit": 1.0, "high_limit": 9.0},
        "Rg": {"value": 30.0, "fit": False},
        "RgCutoff": 2.0,
        "correlated": True,
        "estimate_B": True,
    })
    assert flat["G"] == 3.0 and flat["fit_G"] is True
    assert flat["G_low"] == 1.0 and flat["G_high"] == 9.0
    assert flat["Rg"] == 30.0 and flat["fit_Rg"] is False
    # A bound the config omits falls back to the GUI's wide default.
    assert flat["Rg_low"] == 0.1 and flat["Rg_high"] == 1e6
    assert flat["RgCutoff"] == 2.0 and flat["correlated"] is True


def test_a_hand_written_config_may_use_bare_numbers():
    from pyirena.batch.unified import _flatten_level_config

    flat = _flatten_level_config({"G": 5.0, "Rg": 25.0})
    assert flat["G"] == 5.0 and flat["Rg"] == 25.0


def test_the_batch_path_builds_the_model_the_config_describes():
    from pyirena.batch.unified import _state_to_model

    model = _state_to_model({
        "num_levels": 2,
        "background": {"value": 0.05, "fit": True},
        "levels": [
            {"G": {"value": 10.0, "fit": True, "low_limit": 1.0, "high_limit": 99.0},
             "Rg": {"value": 50.0, "fit": True}, "correlated": False},
            {"G": {"value": 1.0}, "Rg": {"value": 5.0}, "estimate_B": True},
        ],
    })
    assert model.num_levels == 2
    assert model.background == 0.05 and model.fit_background is True
    assert model.levels[0].G == 10.0 and model.levels[0].G_limits == (1.0, 99.0)
    assert model.levels[1].link_B is True


# ── Modeling: six population types, one mechanism ────────────────────────

def _populated_config():
    from pyirena.core.modeling import (
        DiffractionPeakPopulation,
        GuinierPorodPopulation,
        MassFractalPopulation,
        ModelingConfig,
        SizeDistPopulation,
        SurfaceFractalPopulation,
        UnifiedLevelPopulation,
    )

    sd = SizeDistPopulation()
    sd.dist_params['mean_size'] = 123.0
    sd.dist_params_limits['mean_size'] = (10.0, 900.0)
    sd.structure_factor = 'hard_sphere'
    sd.sf_params = {'radius': 50.0, 'volume_fraction': 0.2}
    cfg = ModelingConfig(
        populations=[
            sd,
            UnifiedLevelPopulation(Rg=44.0, correlations=True),
            DiffractionPeakPopulation(peak_type='voigt', position=0.35),
            GuinierPorodPopulation(),
            MassFractalPopulation(),
            SurfaceFractalPopulation(use_porod_transition=True),
        ],
        background=0.004,
        q_min=0.002, q_max=0.8,
        fit_method='global', de_workers=4,
        use_slit_smearing=True, slit_length=0.02,
    )
    return cfg


def test_a_modeling_setup_survives_a_json_round_trip():
    from pyirena.core.modeling import ModelingConfig

    cfg = _populated_config()
    restored = ModelingConfig.from_dict(_through_json(cfg.to_dict()))
    assert restored.to_dict() == cfg.to_dict()


def test_every_population_type_comes_back_as_its_own_class():
    from pyirena.core.modeling import POPULATION_CLASSES, ModelingConfig

    cfg = _populated_config()
    restored = ModelingConfig.from_dict(_through_json(cfg.to_dict()))
    assert [type(p).__name__ for p in restored.populations] == \
           [type(p).__name__ for p in cfg.populations]
    assert len(POPULATION_CLASSES) == 6


def test_nested_bounds_are_tuples_again_after_json():
    """dist_params_limits is a dict of tuples; JSON flattens both levels."""
    from pyirena.core.modeling import ModelingConfig

    restored = ModelingConfig.from_dict(_through_json(_populated_config().to_dict()))
    assert restored.populations[0].dist_params_limits['mean_size'] == (10.0, 900.0)
    assert isinstance(restored.populations[0].dist_params_limits['mean_size'], tuple)
    assert isinstance(restored.populations[1].Rg_limits, tuple)


def test_the_round_trip_keeps_the_modelled_intensity():
    from pyirena.core.modeling import ModelingConfig, ModelingEngine

    cfg = _populated_config()
    restored = ModelingConfig.from_dict(_through_json(cfg.to_dict()))
    q = np.logspace(-2.5, -0.5, 40)
    engine = ModelingEngine()
    before = engine.total_intensity(cfg, q, use_cache=False)[0]
    after = engine.total_intensity(restored, q, use_cache=False)[0]
    assert np.allclose(before, after)


def test_a_population_type_from_a_newer_pyirena_does_not_break_the_load():
    """Better to open the rest of the fit than to refuse the file."""
    from pyirena.core.modeling import SizeDistPopulation, population_from_dict

    pop = population_from_dict({'pop_type': 'something_new', 'enabled': True})
    assert isinstance(pop, SizeDistPopulation)


def test_a_modeling_field_added_since_the_file_was_written_takes_its_default():
    from pyirena.core.modeling import ModelingConfig

    cfg = ModelingConfig.from_dict({'populations': [{'pop_type': 'unified_level',
                                                     'Rg': 12.0}]})
    fresh = ModelingConfig()
    assert cfg.populations[0].Rg == 12.0
    assert cfg.fit_method == fresh.fit_method          # never written back then
    assert cfg.mc_workers == fresh.mc_workers
    assert cfg.populations[0].PACK_limits == (0.0, 16.0)


def test_the_gui_serializer_is_the_dataclass_serializer():
    """The panel used to keep its own copy; they must not diverge again."""
    pytest.importorskip("pyqtgraph")
    try:
        from pyirena.gui.modeling_panel import _pop_from_dict, _pop_to_dict
    except ImportError:
        pytest.skip("Qt not available")

    from pyirena.core.modeling import POPULATION_CLASSES

    for pop_type, cls in POPULATION_CLASSES.items():
        pop = cls()
        assert _pop_to_dict(pop) == pop.to_dict(), pop_type
        # …and the panel's read path understands what its write path produces.
        assert type(_pop_from_dict(_pop_to_dict(pop))) is cls


def test_the_batch_and_gui_paths_agree_on_a_population():
    """Two entry points, one config file, one result."""
    pytest.importorskip("h5py")
    from pyirena.core.modeling import UnifiedLevelPopulation, population_from_dict

    saved = UnifiedLevelPopulation(Rg=33.0, fit_P=True, label='shell').to_dict()
    assert population_from_dict(saved).to_dict() == saved


# ── WAXS Peak Fit ────────────────────────────────────────────────────────

def _waxs_model():
    from pyirena.core.waxs_peakfit import (
        BG_SHAPES,
        PEAK_SHAPES,
        WAXSPeakFitModel,
        default_bg_params,
        default_peak_params,
    )

    model = WAXSPeakFitModel(
        bg_shape=BG_SHAPES[1],
        peaks=[default_peak_params(PEAK_SHAPES[0], 1.5, 2.0, 0.05),
               default_peak_params(PEAK_SHAPES[2], 2.1, 0.5, 0.03)],
    )
    model.bg_params = default_bg_params(BG_SHAPES[1])
    return model


def test_a_waxs_setup_survives_a_json_round_trip():
    from pyirena.core.waxs_peakfit import WAXSPeakFitModel

    model = _waxs_model()
    restored = WAXSPeakFitModel.from_dict(_through_json(model.to_dict()))
    assert restored.to_dict() == model.to_dict()


def test_the_waxs_round_trip_keeps_the_pattern():
    from pyirena.core.waxs_peakfit import WAXSPeakFitModel

    model = _waxs_model()
    restored = WAXSPeakFitModel.from_dict(_through_json(model.to_dict()))
    q = np.linspace(1.0, 2.5, 80)
    assert np.allclose(model.predict(q, model.bg_params, model.peaks),
                       restored.predict(q, restored.bg_params, restored.peaks))


def test_a_waxs_dict_without_a_background_still_opens():
    from pyirena.core.waxs_peakfit import BG_SHAPES, WAXSPeakFitModel

    model = WAXSPeakFitModel.from_dict({'peaks': []})
    assert model.bg_shape in BG_SHAPES
    assert model.bg_params                      # defaults, not None


def test_an_unknown_waxs_background_falls_back_rather_than_raising():
    from pyirena.core.waxs_peakfit import BG_SHAPES, WAXSPeakFitModel

    model = WAXSPeakFitModel.from_dict({'bg_shape': 'Chebyshev-9000'})
    assert model.bg_shape in BG_SHAPES


def test_waxs_to_dict_does_not_hand_out_the_live_peak_list():
    """A caller editing the dict must not edit the model behind its back."""
    model = _waxs_model()
    d = model.to_dict()
    d['peaks'][0]['A']['value'] = 999.0
    assert model.peaks[0]['A']['value'] != 999.0


# ── The rule itself ──────────────────────────────────────────────────────

def test_the_tools_that_claim_the_pair_really_have_it():
    """One import, so a half-added tool is caught here rather than in a panel."""
    from pyirena.core.modeling import ModelingConfig
    from pyirena.core.simple_fits import SimpleFitModel
    from pyirena.core.sizes import SizesDistribution
    from pyirena.core.waxs_peakfit import WAXSPeakFitModel

    for cls in (UnifiedFitModel, ModelingConfig, WAXSPeakFitModel,
                SizesDistribution, SimpleFitModel):
        assert callable(getattr(cls, "to_dict", None)), cls.__name__
        assert callable(getattr(cls, "from_dict", None)), cls.__name__


# ── The embedded panel setup (U10) ───────────────────────────────────────

def _morph_result():
    from pyirena.core.saxs_morph import SaxsMorphConfig, SaxsMorphResult

    q = np.logspace(-2, -1, 20)
    ones = np.ones_like(q)
    return SaxsMorphResult(
        config=SaxsMorphConfig(), chi_squared=1.0, reduced_chi_squared=1.0,
        dof=10, timestamp='2026-08-10 12:00:00',
        data_q=q, data_I=ones, data_dI=ones * 0.1, data_I_corr=ones,
        model_q=q, model_I=ones,
        r_grid=q, gamma_r=ones, spectral_k=q, spectral_F=ones,
        voxelgram=np.zeros((4, 4, 4), dtype=np.uint8), voxel_size=4,
        box_size_A=1000.0, voxel_pitch_A=250.0, phi_actual=0.3, rng_seed_used=1,
    )


def test_saxs_morph_embeds_the_panel_setup(tmp_path):
    """Reopening a SAXS Morph result must restore the controls, as elsewhere.

    The config scalars saved beside the result describe the calculation; they
    do not say which input mode the user was in or where the cursors were.
    """
    pytest.importorskip("h5py")
    from pyirena.io.nxcansas_saxs_morph import save_saxs_morph_results
    from pyirena.io.setup_config import read_setup_config

    path = tmp_path / "morph.h5"
    state = {"box_size_A": 1234.0, "input_mode": "both",
             "power_law_P": 3.5, "q_min": 0.01, "rng_seed": None}
    save_saxs_morph_results(path, _morph_result(), setup_state=state)

    restored = read_setup_config(path, "entry/saxs_morph_results", "saxs_morph")
    assert restored == state


def test_a_result_saved_without_a_setup_still_reads_back(tmp_path):
    """Batch runs and older files carry no setup; that is not an error."""
    pytest.importorskip("h5py")
    from pyirena.io.nxcansas_saxs_morph import (
        load_saxs_morph_results,
        save_saxs_morph_results,
    )
    from pyirena.io.setup_config import read_setup_config

    path = tmp_path / "nosetup.h5"
    save_saxs_morph_results(path, _morph_result())
    assert read_setup_config(path, "entry/saxs_morph_results", "saxs_morph") is None
    assert load_saxs_morph_results(path) is not None      # the result itself is fine


def test_the_shared_setup_dialog_knows_every_tool_that_embeds_one():
    """The two lists have to agree or a tool gets a button that cannot work."""
    pytest.importorskip("pyqtgraph")
    try:
        from pyirena.gui.setup_loader import TOOL_GROUP_PATH, TOOL_LABEL
    except ImportError:
        pytest.skip("Qt not available")

    assert set(TOOL_GROUP_PATH) == set(TOOL_LABEL)
    assert "saxs_morph" in TOOL_GROUP_PATH
    assert TOOL_GROUP_PATH["saxs_morph"] == "entry/saxs_morph_results"


def test_the_documented_policy_lists_every_tool():
    """docs/HDF5_NxcanSAS_structure.md is where the per-tool decision lives.

    A tool missing from that table is the state this step existed to fix: an
    asymmetry that reads as an oversight because nobody wrote down the reason.
    """
    from pathlib import Path

    import pyirena

    doc = (Path(pyirena.__file__).parent.parent / "docs"
           / "HDF5_NxcanSAS_structure.md")
    if not doc.exists():                     # installed package, not a checkout
        pytest.skip("docs are not part of the installed package")
    text = doc.read_text(encoding="utf-8")
    for tool in ("Unified Fit", "Size Distribution", "Modeling", "Simple Fits",
                 "WAXS Peak Fit", "SAXS Morph", "Fractals", "Data Merge",
                 "Data Manipulation"):
        assert tool in text, f"{tool} is missing from the _pyirena_config table"
