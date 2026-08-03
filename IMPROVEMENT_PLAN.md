# pyirena — Code Quality Improvement Plan

Living record of code-review findings and their disposition. Fixed items are
documented in `CHANGELOG.md` and removed from the active plan below.

## History — `cleanup/code-quality` (landed pre-1.1.0)

- Lint & packaging: removed ~130 unused imports, bare `except:`, shadow re-imports;
  added `[tool.ruff]` config; rewrote `requirements.txt`.
- IO dedup: 4 copy-pasted NXcanSAS write helpers → `io/_nxcansas_common.py`.
- CI re-enabled: Python 3.9/3.11/3.13 matrix + ruff lint job + headless GUI job.
- Unified Fit math: `core.unified.compute_invariant_sv()` shared by GUI and batch.
- Monolith splits: `batch.py` → `pyirena/batch/`; `gui/data_selector.py` →
  `pyirena/gui/data_selector/`.
- Logging system: `logging_setup.py`, rotating logs, prints → loggers.
- Qt shim: `gui/_qt.py`; ~30 GUI modules no longer repeat the PySide6/PyQt6 block.
- Test suite: ~180 new tests; single-source version via `importlib.metadata`.

## Resolved in 1.1.0b2 (independent review, 2026-07-23)

The 2026-07-23 review of `feature/slit-smearing` @ `756dc1f` raised P0–P3
findings. All valid ones were fixed in **1.1.0b2** — see `CHANGELOG.md` for
details. Summary of disposition:

- **P0 — Enforce `PYIRENA_DATA_ROOT` on the control/MCP write surface** — Fixed.
  `open_dataset` / `save_fit` / `save_sizes_fit` now resolve reads and writes
  through `resolve_safe*`; tests added; stale `_paths.py` note corrected.
- **P1 — NumPy floor vs `numpy.trapezoid`** — Fixed. Raised floor to NumPy ≥ 2.0
  in `pyproject.toml` and the conda recipe (matches the code).
- **P1 — Slit-smearing control-API parity** — Fixed. `use_slit_smeared` added to
  the `open_dataset` schema; both save adapters now write `slit_length`, the
  ideal (`*_ideal`) curve, and (Sizes) `data_is_slit_smeared`. Tests added.
- **P1 — Green, bounded CI** — Fixed. 12 ruff findings cleared; `pytest-timeout`
  added (300 s watchdog); the long-hanging modeling-export test fixed (it
  blocked on an unmocked modal `QMessageBox`).
- **P1 — `output_path` produces a complete data file** — Fixed. New output paths
  are seeded from the source via `copy_and_strip_results`; original preserved.
- **P2 — Validate the public smearing input contract** — Fixed. Central
  validation in `core/smearing.py` (only when slit length > 0). Tests added.
- **P2 — Packaging / conda / optional-deps** — Fixed. `plotting` extra added;
  matplotlib import made lazy; beta classifier; SPDX license; conda recipe
  aligned (+ `igor2`, documented `sha256`); tests excluded from the wheel/sdist.
- **P2 — MCP / schema contract tests** — Fixed (proportionate). Full 68-tool MCP
  surface + all 51 control schemas checked for structure and schema↔signature
  parity. (See "Deferred" for the exhaustive default-value variant.)
- **P2 — Release-workflow guardrails** — Fixed (core). `publish.yml` now verifies
  the release tag equals the `pyproject.toml` version. (See "Deferred" for the
  heavier test-gating / artifact-reuse variant.)
- **P3 — Docs / distribution housekeeping** — Fixed. `docs/distribution.md`
  updated to the single-source-version reality; `scratch_sizes_diagnosis/`
  untracked; wheel no longer ships tests.

## Reviewed — deferred / not worth fixing now

These were part of the review's recommendations but are disproportionate for a
single-maintainer scientific package at this stage. They are recorded here (not
in the active plan) so the decision is explicit. Revisit if the motivating
problem recurs.

- **Dedicated lowest-supported-dependency CI job.** The NumPy floor now matches
  the code (≥ 2.0), so there is no 1.x/2.x split left to guard. A lowest-pins
  resolver job would add CI cost for little benefit unless lower floors return.
- **Exhaustive schema default-value parity.** The parity test already enforces
  `properties ⊆ parameters`, `required == mandatory parameters`, and that every
  parameter is exposed. Asserting each JSON `default` equals the Python default
  is brittle (several defaults are descriptive, e.g. "scipy default") for little
  extra safety.
- **Full install smoke matrix (core / plotting / gui / mcp / all).** CI already
  runs a core `test` job and a `[gui]` `test-gui` job. A 5-way extras × OS
  install matrix is more machinery than this project needs today.
- **Heavier release pipeline** (gate publish on a re-run of the Tests workflow;
  reuse exact previously-validated artifacts). The cheap, high-value
  tag==version guard is in place; full workflow-gating/attestation is deferred.

## Future / optional — large-file comfort refactors

No bugs, no user impact. Revisit if the files keep growing or become hard to
navigate.

- `gui/data_selector/panel.py` (~2,800 lines) — split by mixin/topic.
- `gui/unified_fit.py` (~4,300 lines) — split along panel/widget lines.
- `gui/modeling_panel.py` (~3,800 lines) and `gui/sizes_panel.py` (~3,000 lines).
