# pyirena — Code Quality Improvement Plan

Status of a code review performed July 2026. The items in "Done" were fixed on
the `cleanup/code-quality` branch. The remaining items are larger refactors
that need care (and in some cases scientific validation), listed in suggested
priority order.

## Done (branch `cleanup/code-quality`)

- Fixed undefined forward references (`SimpleFitModel` in `batch.py`,
  `h5py.*` in `io/pxp_to_nexus.py`) via `TYPE_CHECKING` imports.
- Replaced 4 bare `except:` with `except Exception:` in `gui/unified_fit.py`.
- Removed 2 shadowed re-imports (F811) in `gui/data_selector.py` and
  `gui/waxs_peakfit_panel.py`.
- Removed ~130 unused imports and no-placeholder f-strings (ruff `--fix`,
  every change compile-checked and import-checked).
- Extracted 4 copy-pasted helpers shared by `io/nxcansas_data_merge.py` and
  `io/nxcansas_data_manipulation.py` into `io/_nxcansas_common.py`
  (both files shrank from ~285 to ~160 lines).
- Added `[tool.ruff]` config to `pyproject.toml` (E741/E701/E702/E402
  intentionally ignored — see comments there); added `ruff` to the `dev` extra.
- Re-enabled CI: `.github/workflows/tests.yml` (was `.disabled`), modernized
  to actions/checkout@v4, setup-python@v5, Python 3.9/3.11/3.13 on ubuntu,
  plus a ruff lint job. **Watch the first run** — it has not executed on GitHub
  in a while and may need dependency tweaks.
- Rewrote stale `requirements.txt` (dropped `six`, matched pyproject minimums,
  pointed readers at `pyproject.toml` as the source of truth).

## Remaining — high priority

### 1. Unify GUI model math with `core/` — MOSTLY DONE
Fitting/plotting always used core math. The duplicated invariant/Sv
computation (GUI `_calculate_unified_intensity`/`_sphere_amplitude`, batch
`_compute_invariant_sv`) is now unified into
`core.unified.compute_invariant_sv()` with pinned regression tests, and the
`calculate_invariant()` Sv-units bug (missing 1e4) is fixed.

Still open: the inline local Guinier / power-law `curve_fit` models in
`gui/unified_fit.py` (~lines 3080–3330) overlap with the local-fit tools in
`api/control/unified_fit.py` — compare and unify the same way. Also
`core/modeling.py:_sphere_amplitude` duplicates `unified.sphere_amplitude`.

### 2. Exception-handling audit
461 `except Exception` blocks; ~100 are silent `pass`. In `io/` and `state/`
a swallowed exception can hide file corruption. Suggested policy:

- `io/`, `state/`, `core/`: never swallow silently — log via `logging` and
  re-raise or return an explicit error value.
- `gui/`: swallowing is sometimes legitimate (e.g. removing plot items that
  may already be gone) but should still `logger.debug()` the exception.

Do this module-by-module, with the test suite green after each module.

### 3. Adopt `logging` in place of `print()`
~250 `print()` calls (117 in `batch.py`, 85 in `io/`). Introduce
`logging.getLogger("pyirena.<module>")` per module, keep console output by
default via a root handler in the CLI entry points, and give `batch.py`
functions a `verbose`/`quiet` control. Do together with item 2.

## Remaining — medium priority

### 4. Split monolith modules
- `gui/data_selector.py` (5,193 lines, 12 classes) → one module per class
  under `gui/data_selector/`.
- `batch.py` (2,630 lines, one function per analysis tool) → `batch/`
  package, one module per tool, re-exporting the public names from
  `batch/__init__.py` so `from pyirena.batch import fit_unified` keeps working.
- `gui/unified_fit.py` (4,444 lines) — shrinks naturally after item 1.

Pure moves, no behavior change; do each as a single commit so diffs stay
reviewable.

### 5. Close test gaps
No direct tests for `core/`: `data_manipulation`, `data_merge`,
`distributions`, `form_factors`, `fractals`, `simple_fits`, `waxs_peakfit`,
`scattering_contrast`, `similarity`, `morphology`, `diffraction_lines`.
Also untested: the `io/nxcansas_*` save/load round-trips — these guard the
on-disk file format and are cheap to test with `tmp_path`. Suggested order:
io round-trips first (highest value per line), then `form_factors` /
`distributions` (pure functions, easy), then the rest.

### 6. Reduce remaining ruff findings
`ruff check pyirena` currently reports ~63 leftovers:
- 46 F401 unused imports inside PySide6/PyQt6 `try/except ImportError`
  fallback blocks (ruff won't auto-fix the `try` branch; clean up manually
  and keep both branches importing the same names).
- 17 F841 unused local variables — each needs a human decision (dead code
  vs. intentionally discarded value).

## Remaining — low priority

- **Single-source the version**: `__version__` is duplicated in
  `pyirena/__init__.py` and `pyproject.toml`. Use
  `importlib.metadata.version("pyirena")` in `__init__.py`.
- **Trim CHANGELOG.md** (141 KB): archive pre-1.0 entries to
  `docs/CHANGELOG_archive.md`.
- **Repo clutter**: `pyIrena_icon.png` (426 KB) in root; `codeFragments/`,
  `IgorCodeFragments/`, `planning/` are development scratch — consider moving
  to a separate branch or excluding from the sdist (verify with
  `python -m build` + inspect the tarball).
- **GUI Qt import fallbacks**: the PySide6/PyQt6 dual-import blocks are
  repeated in ~15 files; a single `pyirena/gui/_qt.py` shim would remove the
  duplication and the associated lint noise.
