# pyirena — Code Quality Improvement Plan

Code review performed July 2026. All planned work landed on `cleanup/code-quality`
and is documented in CHANGELOG.md `[Unreleased]`. Branch is ready to merge.

## Completed (branch `cleanup/code-quality`)

- Lint & packaging: removed ~130 unused imports, bare `except:`, shadow re-imports;
  added `[tool.ruff]` config; rewrote `requirements.txt`.
- IO dedup: 4 copy-pasted NXcanSAS write helpers → `io/_nxcansas_common.py`.
- CI re-enabled: Python 3.9/3.11/3.13 matrix + ruff lint job + headless GUI job.
- Unified Fit math: `core.unified.compute_invariant_sv()` shared by GUI and batch;
  fixed missing 1e4 Å⁻¹→m²/cm³ factor on Sv; local Guinier/power-law unified.
- Monolith splits: `batch.py` (2,616 lines) → `pyirena/batch/` package;
  `gui/data_selector.py` (5,206 lines) → `pyirena/gui/data_selector/` package.
- Logging system: `logging_setup.py`, rotating logs in `~/.pyirena/logs/`,
  `batch.py` / `io/` / `state/` prints → loggers, 137 silent GUI excepts now log.
- Qt shim: `gui/_qt.py`; ~30 GUI modules no longer repeat the PySide6/PyQt6 block.
- Test suite: ~180 new tests (370 total) covering io round-trips, form factors,
  distributions, CorMap, all 14 simple-fit models, fractals, morphology, WAXS;
  fixed LSW normalization, CorMap recursion overflow, and CorMap float overflow.
- Single-source version via `importlib.metadata`; CHANGELOG archived pre-0.8 entries.
- Version field: `__version__` from `importlib.metadata.version("pyirena")`.

## Future / optional

These are large-file comfort refactors — no bugs, no user impact. Revisit if the
files keep growing or become hard to navigate.

- `gui/data_selector/panel.py` (~2,800 lines) — could split by mixin/topic
  (file list, plotting, batch actions, menus).
- `gui/unified_fit.py` (~4,300 lines) — split along panel/widget lines.
- `gui/modeling_panel.py` (~3,800 lines) and `gui/sizes_panel.py` (~3,000 lines) —
  same recipe as `data_selector`.

## Test suite TODO (address off the slit-smearing branch)

- **Hanging test:** `pyirena/tests/test_modeling_global_fit.py::TestExportJsonCarriesFitMethod::test_export_includes_selected_fit_method`
  hangs indefinitely and must be deselected to run the suite. Confirmed
  **pre-existing** — it hangs on a clean checkout of `main`/HEAD too, unrelated
  to the slit-smearing work. Likely a `scipy.optimize` global-fit path or a
  multiprocessing/`differential_evolution` worker that never returns in the test
  environment. Investigate: add a timeout guard, force serial workers in the
  test, or mock the export path. Until fixed, run:
  `pytest --deselect "pyirena/tests/test_modeling_global_fit.py::TestExportJsonCarriesFitMethod::test_export_includes_selected_fit_method"`
  (rest of the suite is green: 439 passed, 5 skipped as of 2026-07-21).
