# CLAUDE.md

Orientation file for AI agents working in this repository. Read this first; it
tells you *where* things are and *what rules apply*, not what every function
does. For details, follow the pointers into `docs/`.

Last update date: 05-08-2026 ; version:1.1.0b3

pyIrena is a Python port of the Igor Pro **Irena** small-angle scattering
package (SAXS/SANS/USAXS analysis). Coded almost entirely by Claude; planned,
specified, debugged and validated by Jan Ilavsky. Scientific correctness
outranks everything else in this codebase.

---

## 1. Commands

```bash
pip install -e ".[all]"      # dev install (gui + gui3d + mcp + plotting + dev)
pytest                       # full suite (testpaths = pyirena/tests)
pytest pyirena/tests/api     # api/mcp layer only
ruff check pyirena/          # lint (line-length 100, config in pyproject.toml)
pyirena-gui                  # launch the main GUI
pyirena-mcp                  # run the MCP server (stdio)
```

Console entry points are declared in `pyproject.toml` under
`[project.scripts]` — that is the authoritative list, don't duplicate it here.

Test data lives in `testData/`. Tests must not require a display; Qt-dependent
tests mock dialogs (see `pyirena/tests/api/conftest.py`).

---

## 2. Architecture — the layer stack

**This is the most important section. Every analysis tool in pyIrena is built
from the same seven layers.** Learn the stack once and you can find anything.

| Layer | Path | Responsibility | May import |
|---|---|---|---|
| **core** | `pyirena/core/<tool>.py` | All math + a serialisable model object | numpy, scipy only |
| **io** | `pyirena/io/nxcansas_<tool>.py` | Save/load results into NXcanSAS HDF5 at `entry/<tool>_results`, embedding GUI state as `_pyirena_config` | h5py, core |
| **gui** | `pyirena/gui/<tool>_panel.py` | Thin Qt shell over the core object; full control state round-trips through `StateManager` | Qt (via `gui/_qt.py`), core, io |
| **batch** | `pyirena/batch/<tool>.py` | Headless execution of the same core from a dict or JSON config section | core, io |
| **api** | `pyirena/api/` | Stable, JSON-serialisable facade for AI/scripting; read access plus `api/control/` for interactive fitting sessions | core, io |
| **mcp** | `pyirena/mcp/server.py` | Protocol wrapper exposing `pyirena.api` to MCP clients | api |
| **plotting** | `pyirena/plotting/` | Headless matplotlib rendering | core, matplotlib |

**The golden rule:** all state lives in the core model object and flows
`dict → JSON → HDF5` unchanged. A control that exists only as a Qt widget will
not survive scripting, batch runs, or setup restore — that is a bug, not a
limitation.

The layer stack is the *target* architecture and older tools do not all reach
it yet — serialisation helpers are named `to_dict`/`from_dict` in newer core
modules but not universally, and panel state methods appear as
`_collect_state`, `collect_state`, `_save_state` and `save_state` depending on
vintage. Follow the convention in the module you are editing; follow
`docs/developer_adding_features.md` for anything new. HDF5 group names are
mostly `entry/<tool>_results` but not mechanically derived — Simple Fits writes
`entry/simple_fit_results`. Check `pyirena/io/nxcansas_<tool>.py` rather than
assuming.

Two cross-cutting consumers read the saved HDF5 groups rather than the core
objects:

- `pyirena/gui/data_selector/` — graphing, reports, CSV tabulation
- `pyirena/gui/hdf5viewer/` — trend plots, Igor-experiment export
  (via `io/h5xp_extractor.py` + `io/igor_names.py`)

Support modules that are not tools: `pyirena/state/` (StateManager, setup
save/restore), `pyirena/core/form_factors.py`, `distributions.py`,
`fit_metrics.py`, `smearing.py`, `feature_detect.py`, `similarity.py`.

### Layering invariants — do not break these

1. `core/`, `io/`, `api/`, `batch/` **never import Qt.** Verify with
   `grep -rl "PySide6\|PyQt" pyirena/core pyirena/api pyirena/io pyirena/batch`
   — it must return nothing. `core/`, `io/` and `batch/` are also
   matplotlib-free; `api/plotting.py` and `api/control/` are the only
   non-GUI modules allowed to use matplotlib, and they import it lazily inside
   functions and force the `Agg` backend.
2. All Qt imports in `gui/` go through `pyirena/gui/_qt.py` — never
   `from PySide6 import ...` directly. It normalises the PySide6/PyQt6 fallback.
3. `pyirena.api` returns JSON-serialisable dicts only. No numpy scalars, no
   model objects, no file handles.
4. `from_dict()` must supply defaults for every field so that files written by
   older versions still load. Backwards compatibility of saved data is not
   optional.
5. Optional dependencies (GUI, 3D, MCP, plotting) must degrade gracefully —
   importing `pyirena` with no extras installed must work. See
   `pyirena/tests/test_optional_dep_modules.py`.

---

## 3. Tool map

Each analysis tool occupies one row. **When a tool is added, add a row — do not
restructure this section.** A blank cell means that surface does not exist yet
for that tool; that is often a legitimate gap worth fixing rather than a
deliberate choice.

| Tool | core | io | gui | batch | docs |
|---|---|---|---|---|---|
| Unified Fit | `unified.py` | `nxcansas_unified.py` | `unified_fit.py` | `unified.py` | `unified_fit_gui.md`, `unified_fit_features.md` |
| Size Distribution | `sizes.py` | `nxcansas_sizes.py` | `sizes_panel.py` | `sizes.py` | `sizes_methods.md` |
| Modeling | `modeling.py` | `nxcansas_modeling.py` | `modeling_panel.py` | `modeling.py` | `modeling_gui.md` |
| Simple Fits | `simple_fits.py` | `nxcansas_simple_fits.py` | `simple_fits_panel.py` | `simple.py` | `simple_fits_gui.md` |
| Fractals | `fractals.py` | `nxcansas_fractals.py` | `fractals_panel.py` | — | `fractals_gui.md` |
| SAXS Morph | `saxs_morph.py` | `nxcansas_saxs_morph.py` | `saxs_morph_panel.py` | `saxs_morph.py` | `saxs_morph_gui.md`, `saxs_morph_method_comparison.md` |
| WAXS Peak Fit | `waxs_peakfit.py` | `nxcansas_waxs_peakfit.py` | `waxs_peakfit_panel.py` | `waxs.py` | `waxs_peakfit_gui.md` |
| Data Merge | `data_merge.py` | `nxcansas_data_merge.py` | `data_merge_panel.py` | `merge.py` | `data_merge_gui.md` |
| Data Manipulation | `data_manipulation.py` | `nxcansas_data_manipulation.py` | `data_manipulation_panel.py` | `manipulate.py` | `data_manipulation_gui.md` |
| Scattering Contrast | `scattering_contrast.py` | `contrast_io.py` | `contrast_panel.py` | — | `scattering_contrast_gui.md` |
| Slit Smearing | `smearing.py` | — | `slit_smearing_ui.py` | — | `slit_smearing.md` |
| Feature Identifier | `feature_detect.py` | — | `feature_identifier.py` | — | `feature_identifier.md` |
| Fit Quality Metrics | `fit_metrics.py` | `nxcansas_fit_quality.py` | `quality_display.py` | — | `fit_quality_metrics.md` |

Paths are relative to `pyirena/core/`, `pyirena/io/`, `pyirena/gui/`,
`pyirena/batch/` and `docs/` respectively. Note that batch module names are
occasionally shortened (`simple.py`, `waxs.py`, `merge.py`).

---

## 4. Where to look

Start here rather than grepping. `docs/` has ~33 files; these are the ones that
change how you should work.

| If you are… | Read |
|---|---|
| Adding *any* feature — start here, always | `docs/developer_adding_features.md` (master checklist of every integration surface) |
| Adding a form factor | `docs/developer_adding_form_factors.md` |
| Adding a structure factor | `docs/developer_adding_structure_factors.md` |
| Touching HDF5 read/write | `docs/HDF5_NxcanSAS_structure.md` |
| Working on the api/MCP layer | `pyirena/api/README.md`, `docs/ai_tools_reference.md`, `docs/ai_integration.md` |
| Writing or fixing tests | `docs/testing.md` |
| Working on batch/scripting | `docs/batch_api.md` |
| Importing Igor `.pxp` files | `docs/igor_pxp_import.md` |
| Loading or cleaning input data | `docs/data_import_and_cleaning.md` |
| Packaging or release work | `docs/distribution.md`, `.github/workflows/publish.yml` |
| Understanding a specific GUI panel | `docs/<tool>_gui.md` (see tool map above) |
| Looking for design intent on unfinished work | `planning/` and `IMPROVEMENT_PLAN.md` |

`docs/GUI_README.md`, `docs/gui_quickstart.md`, `docs/QUICK_START.md` and
`docs/usage_guide.md` are user-facing; consult them to keep terminology
consistent with what users see.

---

## 5. Conventions

**Scientific.** Q is in Å⁻¹ throughout unless a function's docstring says
otherwise. Intensity is cm⁻¹ (absolute) where calibration exists. Single-letter
names `I`, `l`, `Q`, `R` are idiomatic here and ruff's `E741` is disabled for
that reason — keep using the physics notation rather than renaming to
`intensity_array`. Preserve the Irena/Igor terminology users already know; when
an algorithm comes from Irena, say so in the docstring and keep the same
parameter names where practical.

**Code.** Line length 100. Google-style docstrings on public functions. Ruff
with `E741`, `E701`, `E702`, `E402` intentionally disabled (see the comments in
`pyproject.toml` — read them before "fixing" a lint suppression). Type hints
where they help; `py.typed` is shipped. Python 3.9 is the floor, so no
`match`, no PEP 604 `X | Y` in runtime annotations.

**GUI.** One panel per tool, one file per panel. Panels are thin — if you are
writing a loop with numpy in a `*_panel.py`, it belongs in `core/`. Every panel
must expose state collection and a `load_state()` so setup save/restore works.

**Testing.** New math needs a test in `pyirena/tests/`. New HDF5 fields need a
round-trip test in `test_nxcansas_roundtrips.py`. New api surface needs a test
in `pyirena/tests/api/`. Tests must be headless and finish well inside the
300 s per-test timeout.

**Versioning.** `pyproject.toml` version is authoritative and the publish
workflow refuses to release if the git tag disagrees. Update `CHANGELOG.md` for
user-visible changes.

---

## 6. Working style in this repo

- Prefer editing the existing tool that already does 80% of the job over adding
  a parallel implementation. Duplication across the seven layers is expensive.
- When a change spans layers, work bottom-up: core → io → batch → gui → api.
  The core model object's `to_dict()` shape is the contract everything else
  depends on.
- Do not add dependencies without asking. The dependency set is deliberately
  small and split across extras.
- `temp/`, `testData/`, `scripts/` and `planning/` are excluded from ruff and
  are not part of the shipped package. `scripts/` holds diagnostic one-offs, not
  supported tooling.
- Ask before large refactors. The codebase is validated against Igor Irena
  results; a "cleaner" rewrite that shifts a numerical result is a regression.

---

## 7. Maintaining this file

This file is a map, not documentation. It is deliberately written at a level of
abstraction that survives ordinary development — most commits should require no
change here.

**Update it when, and only when:**

- a new analysis tool is added → add one row to §3
- a new top-level package appears under `pyirena/` → add one row to §2
- a layering invariant in §2 changes → fix the invariant list
- a command in §1 changes (test runner, lint tool, entry points)
- a new `docs/developer_*.md` guide is written → add one row to §4

**Do not add here:** function signatures, parameter lists, model formulas,
per-release notes, or anything already written in `docs/`. Those rot within
weeks and belong in the code or in `docs/`.

If §3 or §4 has drifted, regenerate just that table by listing the directories
rather than rewriting the whole file.
