"""The completeness contract for a tool's configuration, across all six tools.

A config that loses a setting does not raise — it returns a plausible wrong
number.  Four such bugs were found in one week in one small area, none of them
caught by the suite, because every test that existed compared *values* and the
things being lost were the decisions *about* the values: which parameters the
scientist pinned, and the bounds that kept the fit physical.

So these are loops over the tool table, not one test per tool, in the spirit of
``test_tool_registration.py``: adding a tool means adding a row, and adding a
field to a model means the contract either covers it or makes you say why not.

1. :func:`test_every_setting_is_serialised_or_explicitly_excluded` — a scalar
   setting is in ``to_dict()`` or named in the class's ``_NOT_SERIALISED``
   allowlist with a reason.
2. :func:`test_a_model_rebuilt_from_its_own_dict_can_run_a_fit` — the rebuilt
   object *fits*, rather than merely comparing equal.  ``SimpleFitModel``
   round-tripped into an object that raised ``AttributeError`` on its first
   fit for exactly as long as the only test on it compared dicts.
3. :func:`test_flags_and_bounds_survive_the_config_round_trip` — **the one
   that matters.**  Every fit flag inverted from its default and every bound
   narrowed, through ``to_dict() → build_setup() → to_dict()``, in each of the
   four envelopes pyIrena writes.  Every bug in ``planning/config-dialects/``
   §2.1–2.3 fails this and passes a values-only comparison.
4. :func:`test_every_config_pyirena_has_written_still_replays` — the real
   exported configs and result files in ``testData/Scripting/``, which is the
   only check that the shapes on users' disks are the shapes we read.
"""

from __future__ import annotations

import copy
import json
from pathlib import Path

import numpy as np
import pytest

from pyirena.core.tool_config import TOOL_SECTIONS, build_setup

SCRIPTING = Path(__file__).resolve().parents[2] / "testData" / "Scripting"


# ── The tool table ───────────────────────────────────────────────────────
#
# One row per tool: how to build a populated model, and anything build_setup
# is *allowed* to change on the way back — each with the reason, because an
# unexplained exemption is how a real loss gets waved through.

def _unified():
    from pyirena.core.unified import UnifiedFitModel
    return UnifiedFitModel(num_levels=3)


def _sizes():
    from pyirena.core.sizes import SizesDistribution
    s = SizesDistribution()
    s.shape = "spheroid"
    s.shape_params = {"aspect_ratio": 2.5}
    s.use_slit_smearing = True
    s.slit_length = 0.0217
    return s


def _simple():
    from pyirena.core.simple_fits import SimpleFitModel
    m = SimpleFitModel()
    m.set_model("Porod")
    return m


def _modeling():
    from pyirena.core.modeling import (
        DiffractionPeakPopulation,
        GuinierPorodPopulation,
        MassFractalPopulation,
        ModelingConfig,
        SizeDistPopulation,
        SurfaceFractalPopulation,
        UnifiedLevelPopulation,
    )
    cfg = ModelingConfig()
    # Every population type, because five of the six were the ones being lost.
    cfg.populations = [
        SizeDistPopulation(), UnifiedLevelPopulation(), DiffractionPeakPopulation(),
        GuinierPorodPopulation(), MassFractalPopulation(), SurfaceFractalPopulation(),
    ]
    return cfg


def _waxs():
    from pyirena.core.waxs_peakfit import WAXSPeakFitModel, default_peak_params
    return WAXSPeakFitModel("Linear", [default_peak_params("Gauss"),
                                       default_peak_params("Voigt")])


def _carbon():
    from pyirena.core.carbon_fit import CarbonFitModel
    return CarbonFitModel()


#: ``tool -> (factory, {key: why build_setup may change it})``.
TOOL_TABLE = {
    "unified_fit":  (_unified, {}),
    # Size Distribution is an inversion, not a least-squares fit, so it has
    # no per-parameter fit flags or bounds for the mutator below to invert.
    # Its two Step 0 losses are caught by the other two tests instead: the
    # slit settings by the completeness check (they were in vars() and not in
    # to_dict()), the spheroid aspect ratio by the value round trip (hence
    # the non-default shape on the factory above).
    "sizes":        (_sizes, {
        "montecarlo_n_repetitions":
            "The main fit always uses a single MC run, matching the GUI; the "
            "repetition count belongs to the separate uncertainty pass.",
    }),
    "simple_fits":  (_simple, {}),
    "modeling":     (_modeling, {}),
    "waxs_peakfit": (_waxs, {}),
    "carbon_fit":   (_carbon, {}),
}


def test_the_table_covers_every_tool_build_setup_knows():
    """A tool added to ``TOOL_SECTIONS`` must gain a row here, or it is untested."""
    assert set(TOOL_TABLE) == set(TOOL_SECTIONS)


# ── 1. Completeness: serialised, or excluded on purpose ──────────────────

@pytest.mark.parametrize("tool", sorted(TOOL_TABLE))
def test_every_setting_is_serialised_or_explicitly_excluded(tool):
    model = TOOL_TABLE[tool][0]()
    serialised = model.to_dict()
    excluded = getattr(type(model), "_NOT_SERIALISED", {})

    unaccounted = {
        name: value for name, value in vars(model).items()
        if not name.startswith("_")
        and isinstance(value, (bool, int, float, str))
        and name not in serialised
        and name not in excluded
    }
    assert not unaccounted, (
        f"{type(model).__name__} has settings that to_dict() drops and "
        f"_NOT_SERIALISED does not mention: {sorted(unaccounted)}. Either "
        f"serialise them or add them to _NOT_SERIALISED with the reason — a "
        f"setting that is silently dropped comes back as a different fit."
    )
    for name, reason in excluded.items():
        assert isinstance(reason, str) and len(reason) > 20, (
            f"{type(model).__name__}._NOT_SERIALISED['{name}'] needs a reason "
            f"a reader can act on, not {reason!r}."
        )


# ── 2. The rebuilt object has to fit, not just compare equal ─────────────

def _tiny_curve(q_min=0.002, q_max=0.3, n=120):
    q = np.logspace(np.log10(q_min), np.log10(q_max), n)
    intensity = 1e-3 * q ** -3.6 + 0.01
    return q, intensity, intensity * 0.03


@pytest.mark.parametrize("tool", sorted(TOOL_TABLE))
def test_a_model_rebuilt_from_its_own_dict_can_run_a_fit(tool):
    """Equal dicts are necessary and not sufficient — the object must work.

    ``SimpleFitModel.from_dict`` rebuilt with ``cls.__new__`` and hand-assigned
    fields, so a field set only in ``__init__`` went missing and the first fit
    raised ``AttributeError``. A dict comparison never saw it.
    """
    model = TOOL_TABLE[tool][0]()
    rebuilt = build_setup({tool: json.loads(json.dumps(model.to_dict()))}).model

    q, intensity, error = _tiny_curve()
    # Not "the fit converges" — these are synthetic curves and some tools are
    # being handed a model of the wrong physics. What must not happen is the
    # object falling over because a field never made it back.
    try:
        _run_one_fit(tool, rebuilt, q, intensity, error)
    except AttributeError as exc:
        pytest.fail(
            f"{tool}: a model rebuilt from its own dict is missing an "
            f"attribute its fit needs: {exc}"
        )


def _run_one_fit(tool, model, q, intensity, error):
    """Drive one short fit per tool, however that tool spells it.

    Each tool is capped to a handful of iterations: the question is whether
    the rebuilt object *works*, not whether a synthetic curve converges.
    """
    if tool == "unified_fit":
        return model.fit(q, intensity, error, max_iterations=3)
    if tool == "sizes":
        return model.fit(q, intensity, error)
    if tool == "simple_fits":
        return model.fit(q, intensity, error)
    if tool == "modeling":
        from pyirena.core.modeling import ModelingEngine
        model.q_min, model.q_max = float(q[0]), float(q[-1])
        return ModelingEngine().fit(model, q, intensity, error)
    if tool == "waxs_peakfit":
        return model.fit(q, intensity, error, model.bg_params, model.peaks)
    if tool == "carbon_fit":
        model.max_iterations = 3
        return model.fit(q, intensity, error)
    raise AssertionError(f"no fit driver for {tool}")


# ── 3. Flags and bounds, through every envelope ──────────────────────────

def _number(value):
    """The value as a float, or None when it is not a plain number."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    return float(value)


def _narrow(pair, value):
    """Narrow a bound the way a scientist does: inward, value still inside.

    An *open* bound (``None``, which WAXS and Simple Fits both write) becomes
    a finite one, because typing a limit where there was none is the same
    journey through the config and is just as easy to drop on the way.
    """
    lo, hi, v = pair[0], pair[1], _number(value)
    if lo is None:
        lo = (v - abs(v) - 1.0) if v is not None else -1e3
    if hi is None:
        hi = (v + abs(v) + 1.0) if v is not None else 1e3
    lo, hi = float(lo), float(hi)
    if v is not None and lo <= v <= hi:
        return [lo + (v - lo) * 0.5, hi - (hi - v) * 0.5]
    return [lo + (hi - lo) * 0.25, hi - (hi - lo) * 0.25]


#: The two spellings of a bound inside a ``{"value": …, "fit": …}`` block.
_BOUND_SPELLINGS = (("low_limit", "high_limit"), ("lo", "hi"))


def _is_param_block(value):
    """True for the ``{value, fit, lo/low_limit, hi/high_limit}`` encoding."""
    return (isinstance(value, dict) and "value" in value
            and ("fit" in value or any(k in value for pair in _BOUND_SPELLINGS
                                       for k in pair)))


def _mutate_param_block(block):
    out = dict(block)
    if isinstance(out.get("fit"), bool):
        out["fit"] = not out["fit"]
    for lo_key, hi_key in _BOUND_SPELLINGS:
        if lo_key in out or hi_key in out:
            lo, hi = _narrow([out.get(lo_key), out.get(hi_key)], out.get("value"))
            out[lo_key], out[hi_key] = lo, hi
    return out


def _invert_flags_and_narrow_bounds(node):
    """Flip every fit flag and pull in every bound, at any depth.

    Walks the serialised form rather than the object, so it reaches all five
    encodings of "a parameter with a fit flag and bounds" that the six tools
    use between them (``planning/config-dialects/`` §1.5) without having to
    know which tool spells it which way:

    * ``X`` / ``fit_X`` / ``X_limits``            — Unified, Carbon, Modeling
    * ``X: {value, fit, low_limit, high_limit}``  — the Unified panel
    * ``X: {value, fit, lo, hi}``                 — WAXS peaks and background
    * ``params`` / ``limits`` / ``param_fixed``   — Simple Fits
    * ``dist_params`` / ``dist_params_fit`` / ``dist_params_limits`` — Modeling
    """
    if isinstance(node, list):
        return [_invert_flags_and_narrow_bounds(item) for item in node]
    if not isinstance(node, dict):
        return node

    out = {}
    for key, value in node.items():
        if key.startswith("fit_") and isinstance(value, bool):
            out[key] = not value
        elif key in ("param_fixed",) and isinstance(value, dict):
            # Simple Fits inverts the sense of the flag as well as the name.
            out[key] = {k: (not v) if isinstance(v, bool) else v
                        for k, v in value.items()}
        elif key.endswith("_fit") and isinstance(value, dict):
            # Modeling's {param: bool} maps (dist_params_fit, ff_params_fit, …)
            out[key] = {k: (not v) if isinstance(v, bool) else v
                        for k, v in value.items()}
        elif _is_param_block(value):
            out[key] = _mutate_param_block(value)
        elif key.endswith("_limits") and isinstance(value, list) and len(value) == 2:
            out[key] = _narrow(value, node.get(key[:-len("_limits")]))
        elif key in ("limits", "param_limits") and isinstance(value, dict):
            values = node.get("params") or {}
            out[key] = {
                k: _narrow(v, values.get(k))
                if isinstance(v, (list, tuple)) and len(v) == 2 else v
                for k, v in value.items()
            }
        elif key.endswith("_limits") and isinstance(value, dict):
            values = node.get(key[:-len("_limits")]) or {}
            out[key] = {
                k: _narrow(v, values.get(k))
                if isinstance(v, (list, tuple)) and len(v) == 2 else v
                for k, v in value.items()
            }
        else:
            out[key] = _invert_flags_and_narrow_bounds(value)
    return out


def _envelopes(tool, section):
    """The same settings inside each wrapper pyIrena writes.

    All four have to read back identically; for most of pyIrena's life only
    the first one did, and the two a remote caller actually holds — the one
    inside a result file and the one the service just handed back — were the
    two that were rejected.
    """
    return {
        "sidecar": {"_pyirena_config": {"tool": tool}, tool: copy.deepcopy(section)},
        "hdf5_state": {"_pyirena_config": {"tool": tool}, "state": copy.deepcopy(section)},
        "export_results": {"ok": True, "tool": tool, "config": copy.deepcopy(section),
                           "fit_q_range": {"q_min": None, "q_max": None}},
        "tool_model": {"tool": tool, "model": copy.deepcopy(section)},
        "bare_section": copy.deepcopy(section),
    }


@pytest.mark.parametrize("tool", sorted(TOOL_TABLE))
@pytest.mark.parametrize("envelope", ["sidecar", "hdf5_state", "export_results",
                                      "tool_model", "bare_section"])
def test_flags_and_bounds_survive_the_config_round_trip(tool, envelope):
    """Not the values — the decisions about them.

    A replay that keeps every number and loses every fit flag floats the
    parameters the scientist pinned and drops the bounds that kept the fit
    physical. It does not error, and it returns a plausible number.
    """
    factory, allowed_changes = TOOL_TABLE[tool]
    want = _invert_flags_and_narrow_bounds(
        json.loads(json.dumps(factory().to_dict())))

    config = _envelopes(tool, want)[envelope]
    got = json.loads(json.dumps(
        build_setup(config, tool if envelope == "bare_section" else None).model.to_dict()))

    for key in allowed_changes:
        want.pop(key, None)
        got.pop(key, None)

    differing = sorted(k for k in set(want) | set(got) if want.get(k) != got.get(k))
    assert not differing, (
        f"{tool} in the {envelope} envelope did not replay unchanged: "
        f"{differing}.\n"
        + "\n".join(f"  {k}:\n    wrote {want.get(k)!r}\n    read  {got.get(k)!r}"
                    for k in differing[:3])
    )


def test_simple_fits_held_parameters_survive_the_config_round_trip():
    """The one flag the round trip above cannot see, because it is not in the model.

    Simple Fits keeps "Fit?" outside ``to_dict()``, in the config section, as
    ``param_fixed = {name: True when *held*}`` — a different name from every
    other tool's and the opposite sense. An inverted boolean under a different
    name reads as correct in review and produces a fit that refits everything
    the scientist pinned.
    """
    model = _simple()
    section = model.to_dict()
    held = sorted(section["params"])[:1]
    section["param_fixed"] = {name: (name in held) for name in section["params"]}

    setup = build_setup({"simple_fits": section})
    assert sorted(setup.fixed_params) == held, (
        f"held parameters {held} did not survive; got {sorted(setup.fixed_params)}"
    )
    for name in held:
        assert setup.fixed_params[name] == section["params"][name]

    # And the sense is not inverted: a parameter left free stays free.
    free = [n for n in section["params"] if n not in held]
    for name in free:
        assert name not in setup.fixed_params


# ── 4. The shapes that are actually on users' disks ──────────────────────

#: ``filename -> tool``, the *Export Parameters* sidecars in testData.
SIDECAR_FIXTURES = {
    "UnifiedFit":            "unified_fit",
    "SizeDis.json":          "sizes",
    "SimpleFits_Porod.json": "simple_fits",
    "modeling.json":         "modeling",
    "WAXS":                  "waxs_peakfit",
    "carbon_CE_950.json":    "carbon_fit",
    "carbon_CE_1400.json":   "carbon_fit",
}

#: ``filename -> HDF5 group holding a ``_pyirena_config`` attribute``.
RESULT_FILE_FIXTURES = {
    "UnifiedFit_PP15.h5": "entry/unified_fit_results",
    "SizeDis_PP15.h5":    "entry/sizes_results",
    "Modeling_PP15.h5":   "entry/modeling_results",
}


@pytest.mark.parametrize("filename", sorted(SIDECAR_FIXTURES))
def test_every_config_pyirena_has_written_still_replays(filename):
    """Every real exported config builds a model, with the tool it claims."""
    path = SCRIPTING / filename
    if not path.exists():
        pytest.skip(f"{filename} not present")
    setup = build_setup(json.loads(path.read_text()))
    assert setup.tool == SIDECAR_FIXTURES[filename]
    assert setup.model is not None


@pytest.mark.parametrize("filename", sorted(RESULT_FILE_FIXTURES))
def test_the_setup_embedded_in_a_result_file_replays(filename):
    """The ``{_pyirena_config, state}`` wrapper on every saved ``.h5``.

    This is the envelope a remote caller most obviously holds — it read the
    result file — and it was rejected outright with ``Config has no '<tool>'
    section`` until the reader was widened.
    """
    import h5py

    path = SCRIPTING / filename
    if not path.exists():
        pytest.skip(f"{filename} not present")
    with h5py.File(path, "r") as f:
        group = RESULT_FILE_FIXTURES[filename]
        if group not in f:
            pytest.skip(f"{filename} has no {group}")
        raw = f[group].attrs.get("_pyirena_config")
    if raw is None:
        pytest.skip(f"{filename} carries no embedded setup")

    setup = build_setup(json.loads(raw.decode() if isinstance(raw, bytes) else raw))
    assert setup.model is not None


def test_a_saved_modeling_setup_keeps_every_population_it_described():
    """The panel nests all but size-distribution populations; readers must look.

    Read flat, a ``unified_level`` population comes back at its dataclass
    defaults — Rg 10 Å instead of the 1e10 Å that switches the Guinier term
    off — so the fit runs and models something else entirely.
    """
    import h5py

    path = SCRIPTING / "Modeling_PP15.h5"
    if not path.exists():
        pytest.skip("Modeling_PP15.h5 not present")
    with h5py.File(path, "r") as f:
        raw = f["entry/modeling_results"].attrs["_pyirena_config"]
    envelope = json.loads(raw.decode() if isinstance(raw, bytes) else raw)

    populations = build_setup(envelope).model.populations
    on_disk = envelope["state"]["populations"]
    assert len(populations) == len(on_disk)

    from pyirena.core.modeling import POPULATION_PANEL_BLOCKS

    checked = 0
    for saved, pop in zip(on_disk, populations):
        block = POPULATION_PANEL_BLOCKS.get(saved.get("pop_type"))
        if not block or not isinstance(saved.get(block), dict):
            continue
        for key, value in saved[block].items():
            if isinstance(value, (int, float, bool, str)) and hasattr(pop, key):
                assert getattr(pop, key) == value, (
                    f"{saved['pop_type']}.{key} replayed as {getattr(pop, key)!r}, "
                    f"but the file says {value!r}"
                )
                checked += 1
    assert checked, "the fixture no longer contains a nested population block"


# ── 5. The documentation describes a config that actually works ──────────

def test_the_config_examples_in_the_docs_still_build():
    """``docs/batch_api.md`` shows a config; it has to be one.

    A worked example is the first thing someone copies, and the only part of
    the documentation that can be checked mechanically. Both dialects appear
    there — the current one in the annotated example, the legacy one in the
    collapsed block beside it — and both must still read.
    """
    import re

    doc = Path(__file__).resolve().parents[2] / "docs" / "batch_api.md"
    if not doc.exists():
        pytest.skip("docs/batch_api.md not present")

    blocks = re.findall(r"```json\n(.*?)```", doc.read_text(), re.DOTALL)
    assert len(blocks) >= 2, "the config-format section lost its examples"

    setup = build_setup(json.loads(blocks[0]))
    assert setup.tool == "unified_fit"
    assert setup.model.num_levels == 2
    assert (setup.fit_q_min, setup.fit_q_max) == (0.003, 0.45)
    # The documented flags and bounds have to arrive, not just the values.
    second = setup.model.levels[1]
    assert second.Rg_limits == (1.0, 100.0)
    assert second.link_RGCO is True
    assert second.fit_P is False

    legacy_level = json.loads(blocks[1])
    legacy = build_setup({"unified_fit": {"num_levels": 1, "levels": [legacy_level]}})
    assert legacy.model.levels[0].G_limits == (1e8, 1e12)
    assert legacy.model.levels[0].fit_P is False
