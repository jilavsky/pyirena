"""
Source-level guard: every analysis tool is wired into every registry.

pyIrena has grown enough surfaces that adding a tool means touching ~45 files,
and about a dozen of those touches are a **key in a hand-maintained registry**
— a dict or list somewhere that enumerates the tools.  Miss one and the tool
still works in the GUI while silently vanishing from Igor export, the batch
pipeline, setup restore, or the agent surface.  That is not hypothetical: when
this test was written it immediately found two live gaps (SAXS Morph missing
from ``batch.pipeline``, and SAXS Morph + Fractals missing from
``PYIRENA_RESULT_GROUPS``, so their results leaked into derived data files).

The table below is therefore the **canonical tool registry of pyIrena** — the
one place that knows all the tools and, per surface, either the key that
registers them or a written reason why they are deliberately absent.  Adding a
tool means adding a row here and then making this file green; a gap you did not
consider fails the build, a gap you did consider is one line of documented
reasoning.

The *keys differ between registries* — ``size_distribution`` in the HDF5 schema
is ``sizes`` to the batch pipeline, ``simple_fit`` to the Data Explorer and
``simple`` to the MCP dispatcher.  That is baked into saved files and config
sections and cannot be renamed, so the table carries the alias per surface
rather than pretending there is one spelling.

Source text and import-only: no Qt, no display, so this runs in every CI job.
See ``docs/developer_adding_features.md`` § "Adding a whole new tool".
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from pyirena.io._nxcansas_common import PYIRENA_RESULT_GROUPS
from pyirena.io.igor_names import TOOL_CROSS_REF
from pyirena.io.schema import TOOL_REGISTRY

REPO = Path(__file__).resolve().parents[2]


def _src(rel: str) -> str:
    return (REPO / rel).read_text(encoding="utf-8")


# ── The canonical table ─────────────────────────────────────────────────────
#
# One row per tool.  A surface value is either the key that registers the tool
# there, or ``None`` — in which case ABSENT_REASON must explain it.

#: schema key → per-surface key (None = not registered on that surface)
TOOLS: dict[str, dict[str, str | None]] = {
    "unified_fit": {
        "igor":     "unified_fit",           # io/igor_names.py TOOL_CROSS_REF
        "setup":    "unified_fit",           # gui/setup_loader.py TOOL_GROUP_PATH
        "batch":    "unified_fit",           # batch/pipeline.py _TOOL_REGISTRY
        "state":    "unified_fit",           # state/state_manager.py defaults
        "report":   "fit_results",           # core/reporting.py _build_report kwarg
        "api":      "read_unified_fit",      # api/results.py reader
        "control":  "unified_fit",           # api/control/<module>.py
        "dispatch": "unified",               # mcp/dispatch.py _CATEGORY_BY_MODULE
        "viewer":   "unified_fit",           # gui/hdf5viewer/pyirena_readers.py
    },
    "size_distribution": {
        "igor": "size_distribution", "setup": "sizes", "batch": "sizes",
        "state": "sizes", "report": "sizes_results", "api": "read_size_distribution",
        "control": "sizes", "dispatch": "sizes", "viewer": "sizes",
    },
    "simple_fits": {
        "igor": "simple_fits", "setup": "simple_fits", "batch": "simple_fits",
        "state": "simple_fits", "report": "simple_fit_results", "api": "read_simple_fit",
        "control": "simple_fits", "dispatch": "simple", "viewer": "simple_fit",
    },
    "modeling": {
        "igor": "modeling", "setup": "modeling", "batch": "modeling",
        "state": "modeling", "report": "modeling_results", "api": "read_modeling",
        "control": "modeling", "dispatch": "modeling", "viewer": "modeling",
    },
    "waxs_peakfit": {
        "igor": "waxs_peakfit", "setup": "waxs_peakfit", "batch": "waxs_peakfit",
        "state": "waxs_peakfit", "report": "waxs_peakfit_results",
        "api": "read_waxs_peakfit", "control": "waxs_peakfit", "dispatch": "waxs",
        "viewer": "waxs",
    },
    "saxs_morph": {
        "igor": "saxs_morph", "setup": "saxs_morph", "batch": "saxs_morph",
        "state": "saxs_morph", "report": "saxs_morph_results", "api": "read_saxs_morph",
        "control": None, "dispatch": None, "viewer": None,
    },
    "fractals": {
        "igor": "fractals", "setup": None, "batch": None, "state": "fractals",
        "report": None, "api": "read_fractals", "control": None, "dispatch": None,
        "viewer": None,
    },
    "data_merge": {
        "igor": "data_merge", "setup": None, "batch": None, "state": "data_merge",
        "report": None, "api": "read_merge_provenance", "control": None,
        "dispatch": "data", "viewer": None,
    },
    "data_manipulation": {
        "igor": "data_manipulation", "setup": None, "batch": None,
        "state": "data_manipulation", "report": None,
        "api": "read_manipulation_provenance", "control": None, "dispatch": "data",
        "viewer": None,
    },
}

#: (tool, surface) → why it is deliberately not registered there.
#: Every ``None`` above needs an entry; every entry needs a ``None``.
ABSENT_REASON: dict[tuple[str, str], str] = {
    ("saxs_morph", "control"):
        "Visualization/morphology, not a least-squares fit — no agent control "
        "surface. Decided per tool in planning/ai-agent/README.md.",
    ("saxs_morph", "dispatch"):
        "Follows from having no control surface.",
    ("saxs_morph", "viewer"):
        "No per-file scalar worth trending yet. Open gap rather than a "
        "decision — wire it if a user asks for SAXS Morph trend plots.",

    ("fractals", "setup"):
        "Visualization tool with no JSON config by design (PLAN.md).",
    ("fractals", "batch"):
        "Visualization tool, not an analysis technique — no batch use case "
        "(PLAN.md, 'A batch path for Fractals').",
    ("fractals", "report"):
        "Produces aggregates and images, not fitted parameters.",
    ("fractals", "control"): "No batch path, so no agent control surface.",
    ("fractals", "dispatch"): "Follows from having no control surface.",
    ("fractals", "viewer"): "No fitted scalars to trend.",

    ("data_merge", "setup"):
        "Writes its own config file next to its pipeline stage, not the shared "
        "analysis setup (PLAN.md, 'A shared batch config loader').",
    ("data_merge", "batch"):
        "Has batch/merge.py, but runs from its own config at a different "
        "pipeline step — deliberately not in fit_pyirena's registry.",
    ("data_merge", "report"): "Produces a data file and provenance, not a fit.",
    ("data_merge", "control"):
        "Reachable through the dispatcher's stateless 'data' category instead "
        "of a session-based control module.",
    ("data_merge", "viewer"): "Provenance, not trendable parameters.",

    ("data_manipulation", "setup"): "Same as data_merge — own config file.",
    ("data_manipulation", "batch"): "Same as data_merge — own config file.",
    ("data_manipulation", "report"): "Produces a data file and provenance, not a fit.",
    ("data_manipulation", "control"): "Reachable through the 'data' category.",
    ("data_manipulation", "viewer"): "Provenance, not trendable parameters.",
}

SURFACES = ("igor", "setup", "batch", "state", "report", "api", "control",
            "dispatch", "viewer")


# ── The table itself must stay honest ───────────────────────────────────────

def test_table_covers_every_tool_in_the_schema():
    """TOOL_REGISTRY is the source of truth for what a 'tool' is."""
    assert set(TOOLS) == set(TOOL_REGISTRY), (
        "pyirena/tests/test_tool_registration.py::TOOLS is out of step with "
        "io/schema.py::TOOL_REGISTRY. Add the new tool's row (and a reason for "
        "every surface it is deliberately absent from)."
    )


def test_every_row_declares_every_surface():
    for tool, row in TOOLS.items():
        assert set(row) == set(SURFACES), f"{tool}: expected keys {SURFACES}, got {sorted(row)}"


def test_every_gap_has_a_recorded_reason():
    gaps = {(t, s) for t, row in TOOLS.items() for s in SURFACES if row[s] is None}
    missing = gaps - set(ABSENT_REASON)
    assert not missing, (
        f"Unexplained gaps: {sorted(missing)}. Either wire the tool into that "
        f"surface or add a line to ABSENT_REASON saying why not."
    )
    stale = set(ABSENT_REASON) - gaps
    assert not stale, f"ABSENT_REASON explains a gap that no longer exists: {sorted(stale)}"


# ── io layer ────────────────────────────────────────────────────────────────

def test_igor_cross_ref_covers_every_tool():
    """Missing here = the tool's curves cannot be exported to an Igor experiment."""
    for tool, row in TOOLS.items():
        if row["igor"] is None:
            continue
        assert row["igor"] in TOOL_CROSS_REF, (
            f"{tool}: no entry in io/igor_names.py::TOOL_CROSS_REF"
        )


def test_result_groups_cover_every_tool():
    """Missing here = stale results leak into files seeded from a source file."""
    for tool, schema in TOOL_REGISTRY.items():
        assert schema["group"] in PYIRENA_RESULT_GROUPS, (
            f"{tool}: {schema['group']} is not stripped by copy_and_strip_results"
        )


def test_schema_group_names_match_the_io_module():
    """entry/<x>_results must have a matching io/nxcansas_*.py that writes it."""
    io_dir = REPO / "pyirena" / "io"
    for tool, schema in TOOL_REGISTRY.items():
        group = schema["group"]
        hits = [p.name for p in io_dir.glob("nxcansas_*.py") if group in p.read_text()]
        hits += [p.name for p in io_dir.glob("*_io.py") if group in p.read_text()]
        assert hits, f"{tool}: no io module mentions {group!r}"


# ── batch ───────────────────────────────────────────────────────────────────

def test_pipeline_registry_matches_the_table():
    """fit_pyirena silently skips config sections it does not know about."""
    src = _src("pyirena/batch/pipeline.py")
    body = src.split("_TOOL_REGISTRY", 1)[1].split("\n    }", 1)[0]
    registered = set(re.findall(r"^\s+'([a-z_]+)':", body, re.M))
    expected = {row["batch"] for row in TOOLS.values() if row["batch"]}
    assert registered == expected, (
        f"batch/pipeline.py::_TOOL_REGISTRY has {sorted(registered)}, "
        f"table expects {sorted(expected)}"
    )


# ── state and setup ─────────────────────────────────────────────────────────

def test_state_manager_has_a_defaults_block():
    """No defaults block = the panel's settings do not survive a restart."""
    src = _src("pyirena/state/state_manager.py")
    for tool, row in TOOLS.items():
        if row["state"] is None:
            continue
        assert re.search(rf'^\s+"{re.escape(row["state"])}":\s*\{{', src, re.M), (
            f"{tool}: no \"{row['state']}\" section in state_manager defaults"
        )


def test_setup_loader_maps_agree():
    """Missing here = 'Load Setup from File…' cannot restore the tool."""
    src = _src("pyirena/gui/setup_loader.py")
    for name in ("TOOL_GROUP_PATH", "TOOL_LABEL"):
        body = src.split(name, 1)[1].split("}", 1)[0]
        got = set(re.findall(r'"([a-z_]+)":', body))
        expected = {row["setup"] for row in TOOLS.values() if row["setup"]}
        assert got == expected, f"setup_loader.{name} has {sorted(got)}, expected {sorted(expected)}"


def test_setup_loader_group_paths_match_the_schema():
    src = _src("pyirena/gui/setup_loader.py")
    body = src.split("TOOL_GROUP_PATH", 1)[1].split("}", 1)[0]
    paths = dict(re.findall(r'"([a-z_]+)":\s*"([^"]+)"', body))
    for tool, row in TOOLS.items():
        if row["setup"] is None:
            continue
        assert paths[row["setup"]] == TOOL_REGISTRY[tool]["group"], (
            f"{tool}: setup_loader group path disagrees with io/schema.py"
        )


# ── reporting ───────────────────────────────────────────────────────────────

def test_report_builder_accepts_every_fitting_tool():
    """Missing here = results never reach the Markdown report or export_fit_report."""
    import inspect

    from pyirena.core import reporting

    params = set(inspect.signature(reporting._build_report).parameters)
    for tool, row in TOOLS.items():
        if row["report"] is None:
            continue
        assert row["report"] in params, (
            f"{tool}: core/reporting.py::_build_report has no {row['report']!r} argument"
        )


# ── api / mcp ───────────────────────────────────────────────────────────────

def test_api_exposes_a_reader_for_every_tool():
    import pyirena.api as api

    for tool, row in TOOLS.items():
        if row["api"] is None:
            continue
        assert hasattr(api, row["api"]), f"{tool}: pyirena.api has no {row['api']}()"


def test_control_modules_exist_and_are_dispatched():
    from pyirena.mcp import dispatch

    for tool, row in TOOLS.items():
        if row["control"] is None:
            continue
        mod = REPO / "pyirena" / "api" / "control" / f"{row['control']}.py"
        assert mod.exists(), f"{tool}: missing api/control/{row['control']}.py"
        assert dispatch._CATEGORY_BY_MODULE.get(row["control"]) == row["dispatch"], (
            f"{tool}: mcp/dispatch.py maps {row['control']!r} to "
            f"{dispatch._CATEGORY_BY_MODULE.get(row['control'])!r}, "
            f"table expects {row['dispatch']!r}"
        )


def test_every_dispatch_category_has_a_blurb():
    """A category with no blurb is invisible to the agent choosing a tool."""
    from pyirena.mcp import dispatch

    for category in set(dispatch._CATEGORY_BY_MODULE.values()):
        assert category in dispatch._CATEGORY_BLURBS, f"no blurb for category {category!r}"


def test_control_tools_all_have_schemas():
    """A control function with no schema cannot be called by an agent."""
    from pyirena.api.control import schemas as control_schemas
    from pyirena.mcp import dispatch

    missing = [n for n in dispatch._REGISTRY
               if n not in control_schemas.TOOL_SCHEMA_BY_NAME
               and n not in dispatch.CALCULATOR_SCHEMA_BY_NAME
               and n not in dispatch.DATA_OP_SCHEMA_BY_NAME]
    assert not missing, f"control functions with no JSON schema: {sorted(missing)}"


# ── Data Explorer ───────────────────────────────────────────────────────────

def test_hdf5_viewer_detects_every_wired_tool():
    """Missing here = the tool's results are invisible to trend plots."""
    src = _src("pyirena/gui/hdf5viewer/pyirena_readers.py")
    for tool, row in TOOLS.items():
        if row["viewer"] is None:
            continue
        group = TOOL_REGISTRY[tool]["group"]
        assert group in src, f"{tool}: pyirena_readers.py never looks for {group!r}"
        assert f'"{row["viewer"]}"' in src, (
            f"{tool}: pyirena_readers.py does not use the key {row['viewer']!r}"
        )


# ── GUI panels ──────────────────────────────────────────────────────────────

def test_every_tool_has_a_gui_panel_and_docs():
    """A tool nobody can reach from the GUI, or read about, is not finished."""
    panels = {p.name for p in (REPO / "pyirena" / "gui").glob("*.py")}
    docs = {p.name for p in (REPO / "docs").glob("*.md")}
    # tool key → (panel file, doc file); unified_fit predates the _panel suffix
    known = {
        "unified_fit":       ("unified_fit.py", "unified_fit_gui.md"),
        "size_distribution": ("sizes_panel.py", "sizes_methods.md"),
        "simple_fits":       ("simple_fits_panel.py", "simple_fits_gui.md"),
        "modeling":          ("modeling_panel.py", "modeling_gui.md"),
        "waxs_peakfit":      ("waxs_peakfit_panel.py", "waxs_peakfit_gui.md"),
        "saxs_morph":        ("saxs_morph_panel.py", "saxs_morph_gui.md"),
        "fractals":          ("fractals_panel.py", "fractals_gui.md"),
        "data_merge":        ("data_merge_panel.py", "data_merge_gui.md"),
        "data_manipulation": ("data_manipulation_panel.py", "data_manipulation_gui.md"),
    }
    assert set(known) == set(TOOL_REGISTRY), "add the new tool's panel/doc pair here"
    for tool, (panel, doc) in known.items():
        assert panel in panels, f"{tool}: missing pyirena/gui/{panel}"
        assert doc in docs, f"{tool}: missing docs/{doc}"


@pytest.mark.parametrize("tool", sorted(TOOLS))
def test_tool_is_launchable_from_the_data_selector(tool):
    """Every tool is reached from the Data Browser — the app's entry point."""
    src = _src("pyirena/gui/data_selector/panel.py")
    stem = TOOL_REGISTRY[tool]["group"].removeprefix("entry/").removesuffix("_results")
    assert stem in src or tool in src, (
        f"{tool}: gui/data_selector/panel.py never mentions it — the tool has "
        f"no launcher in the Data Browser"
    )
