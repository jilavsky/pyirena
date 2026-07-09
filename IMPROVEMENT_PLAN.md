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

**Test gaps — closed** (~180 new tests; suite now ~370)
- io round-trips for every `nxcansas_*` writer/reader pair, plus
  `create_nxcansas_file`, data merge/manipulation writers and the
  `load_result` dispatcher.
- `form_factors` and `distributions` against physical limits and known
  values (caught and fixed the missing LSW (3/2)^(11/3) normalization,
  which made `lsw_cdf` discontinuous at the cutoff).
- `data_manipulation`, `data_merge`, `morphology`, `waxs_peakfit` peak
  math, `fractals` growth + `compute_fractal_params` (df = dmin·c) +
  `intensity_unified` (plateau = z, Porod slope −4).
- `similarity` / CorMap — caught and fixed two real bugs: recursion
  overflow for n ≳ 300 points, and a float overflow returning p = 0 for
  n ≳ 1024 (silently flagging long datasets as dissimilar). Pinned to an
  exact analytic case.
- All 14 simple-fit models: finite evaluation, fit self-consistency from
  perturbed starts (parameter recovery for the 11 identifiable models,
  curve reproduction for the 3 degenerate ones), dict round-trips.
- `scattering_contrast` against textbook SLDs and `diffraction_lines`
  against the Si powder pattern (new `testData/Si.cif` fixture; these
  skip without the optional GUI extras).
- The 2 old `test_modeling_report_csv` "failures" were pytest ≥ 8.2
  skip-behavior, not bugs; they now skip cleanly without Qt.

## Done — high priority

### 1. Model-math unification — done
- The inline local Guinier / power-law `curve_fit` models in
  `gui/unified_fit.py` and the local-fit tools in
  `api/control/unified_fit.py` now delegate to a single core
  implementation: `core.unified.fit_local_guinier` and
  `core.unified.fit_local_power_law`. Numbers pinned first by
  `tests/test_local_fits.py` (GUI-unweighted and API-weighted paths both
  reproduce their prior results exactly). The GUI power-law fit now also
  drops non-positive points before fitting, matching the API.
- `core/modeling.py:_sphere_amplitude` now imports the single module-level
  `core.unified.sphere_amplitude`; the duplicate was removed.

**Monolith splits — done** (one commit each, pure moves, no behavior change)
- `batch.py` (2,616 lines) → `pyirena/batch/` package: `_common`,
  `unified`, `sizes`, `simple`, `waxs`, `merge`, `modeling`,
  `saxs_morph`, `manipulate`, `convert`, `pipeline`. Largest module now
  485 lines. All previous imports (`from pyirena.batch import
  fit_unified`, tests' `_compute_invariant_sv`) re-exported unchanged.
- `gui/data_selector.py` (5,206 lines) → `pyirena/gui/data_selector/`
  package: `_qt` (single PySide6/PyQt6 import shim), `plot_utils`,
  `sorting`, `report`, `config_dialogs`, `results_windows`, `workers`,
  `igor_import`, `panel`. All legacy names re-exported. Verified by
  static analysis + stubbed-Qt import execution; **needs one manual GUI
  launch check** (sandbox could not load Qt system libraries).

## Done — medium priority

### 3. GUI logging follow-through — done
- 137 silent `except ...: pass` in `gui/` now emit
  `log.debug(..., exc_info=True)` (file-only DEBUG, no console noise), each
  file given a module logger.
- `print()` in `gui/` panels and the non-demo `core/` diagnostics converted
  to loggers; intentional CLI echoes and `gui/launch.py`'s
  dependency-missing message were intentionally kept.
- `pyirena/gui/_qt.py` shim added; ~30 GUI modules now import Qt names from
  the single shim instead of repeating the PySide6/PyQt6 dual-import block.

### 4. Ruff findings — done (63 → 0)
- 45 F401: removed by the `_qt.py` shim migration (+ two stray non-Qt
  unused imports).
- 18 F841: each reviewed individually — dead assignments removed,
  side-effect-only calls stripped of their unused target, QApplication refs
  in tests renamed `_app`, and one genuine ruff false positive
  (`pxp_to_nexus.py` reassigned-exception var) suppressed with a documented
  `# noqa`.

## Done — low priority

- **CI job with GUI extras** — added a `test-gui` job installing `.[gui]`
  and running the suite headless (`QT_QPA_PLATFORM=offscreen`).
- **Single-source the version** — `__version__` now comes from
  `importlib.metadata.version("pyirena")`.
- **Trim CHANGELOG.md** — entries 0.7.2 and older archived to
  `docs/CHANGELOG_archive.md` (2540 → ~830 lines).
- **Repo clutter** — `codeFragments/`, `IgorCodeFragments/`, `planning/`,
  and root `pyIrena_icon.png` excluded from the sdist (verified with
  `python -m build`; the scratch dirs are kept in the working tree).

## Remaining — medium priority

### 2. Further monolith reduction (optional)
- `gui/data_selector/panel.py` is still ~2,800 lines (the
  DataSelectorPanel class itself) — could be split by mixin/topic
  (file list, plotting, batch actions, menus) if it keeps growing.
- `gui/unified_fit.py` (~4,300 lines) — split along panel/widgets lines
  if desired.
- `gui/modeling_panel.py` (~3,800) and `gui/sizes_panel.py` (~3,000) are
  next in line by the same recipe as data_selector.
