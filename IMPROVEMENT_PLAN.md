# pyirena — Code Quality Improvement Plan

Status of a code review performed July 2026, updated as work lands on the
`cleanup/code-quality` branch. Remaining items are listed in suggested
priority order; each should ship as its own reviewable commit with the test
suite green.

## Done (branch `cleanup/code-quality`)

**Lint & packaging cleanup**
- Fixed undefined forward references (`SimpleFitModel` in `batch.py`,
  `h5py.*` in `io/pxp_to_nexus.py`); replaced bare `except:`; removed
  shadowed re-imports and ~130 unused imports (all compile/import-verified).
- Extracted the 4 copy-pasted NXcanSAS write helpers into
  `io/_nxcansas_common.py`.
- Added `[tool.ruff]` config (E741/E701/E702/E402 intentionally ignored);
  `ruff` added to the `dev` extra.
- Re-enabled CI (`tests.yml`: Python 3.9/3.11/3.13 + ruff lint job).
  **Watch the first GitHub run** — it may need dependency tweaks.
- Rewrote stale `requirements.txt` (pyproject.toml is the source of truth).

**Unified Fit math unification**
- Single Igor-faithful invariant/Sv implementation:
  `core.unified.compute_invariant_sv()`; GUI panel and batch delegate to it.
  Pinned regression tests guarantee the numbers did not change.
- Fixed `calculate_invariant()`: missing 1e4 Å⁻¹→m²/cm³ factor on Sv,
  Porod-tail gating aligned with Igor.

**Logging (Phases A–C)**
- `pyirena/logging_setup.py`: rotating logs in `~/.pyirena/logs/`
  (2 MB × 5 backups per tool), file=DEBUG, console=INFO, session header,
  uncaught-exception hook. Wired into all GUI entry points (`gui.log`) and
  the MCP server (`mcp.log`). "Locate Logs…" button in the data selector.
- `batch.py` (114 prints), `io/` progress messages, and `state_manager.py`
  converted to logging; unconfigured scripts keep seeing console output via
  `ensure_console_output()`. Pretty-printers (`print_*_results`) and CLI
  summaries intentionally keep `print()`.
- All broad silent `except Exception: pass` in `io/`, `state/`, `core/` now
  log tracebacks; one promoted to WARNING (scattering_contrast silently
  omitted an element's f1 contribution when absent from Chantler tables).

## Remaining — high priority

### 1. Close test gaps
No direct tests for `core/`: `data_manipulation`, `data_merge`,
`distributions`, `form_factors`, `fractals`, `simple_fits`, `waxs_peakfit`,
`scattering_contrast`, `similarity`, `morphology`, `diffraction_lines`.
Also untested: the `io/nxcansas_*` save/load round-trips — these guard the
on-disk file format and are cheap to test with `tmp_path`. Suggested order:
io round-trips first (highest value per line), then `form_factors` /
`distributions` (pure functions, easy), then the rest.

Also: 2 pre-existing failures in `test_modeling_report_csv.py` (present on
`main` before the cleanup branch) need a fix or an updated expectation.

### 2. Finish the model-math unification
- The inline local Guinier / power-law `curve_fit` models in
  `gui/unified_fit.py` (~lines 3080–3330) overlap with the local-fit tools
  in `api/control/unified_fit.py` — compare outputs, then unify into core
  the same way as the invariant (pin numbers first).
- `core/modeling.py:_sphere_amplitude` duplicates
  `core/unified.py:sphere_amplitude` — keep one.

## Remaining — medium priority

### 3. Split monolith modules
- `gui/data_selector.py` (5,206 lines, 12 classes) → one module per class
  under `gui/data_selector/`.
- `batch.py` (2,616 lines, one function per analysis tool) → `batch/`
  package, one module per tool, re-exporting the public names from
  `batch/__init__.py` so `from pyirena.batch import fit_unified` keeps
  working.
- `gui/unified_fit.py` (4,340 lines) — shrinks further after item 2.

Pure moves, no behavior change; one commit per split.

### 4. GUI logging follow-through
- ~133 silent `except ...: pass` in `gui/` — many are legitimate
  widget-lifetime races, but each should get
  `log.debug(..., exc_info=True)` opportunistically when the file is next
  touched.
- ~32 `print()` in `gui/` panels and ~14 in `core/` (mostly `__main__`
  demos) → logger.
- A `pyirena/gui/_qt.py` shim for the repeated PySide6/PyQt6 dual-import
  blocks (~15 files) would remove duplication *and* most of the remaining
  lint noise (see item 5) in one move — do these together.

### 5. Reduce remaining ruff findings (63)
- 46 F401: unused names inside the PySide6/PyQt6 `try/except ImportError`
  fallback blocks — solved for free by the `_qt.py` shim in item 4.
- 17 F841 unused local variables — each needs a human decision (dead code
  vs. intentionally discarded value).

## Remaining — low priority

- **Single-source the version**: `__version__` is duplicated in
  `pyirena/__init__.py` and `pyproject.toml`. Use
  `importlib.metadata.version("pyirena")` in `__init__.py`.
- **Trim CHANGELOG.md** (141 KB): archive old entries to
  `docs/CHANGELOG_archive.md`.
- **Repo clutter**: `pyIrena_icon.png` (426 KB) in root; `codeFragments/`,
  `IgorCodeFragments/`, `planning/` are development scratch — consider
  moving to a separate branch or excluding from the sdist (verify with
  `python -m build` + inspect the tarball).
