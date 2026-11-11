# One fit runner per tool — decisions and the remaining GUI step

Started 10-11-2026 for issue #38 (review §3.2); the batch/API part was merged
11-11-2026 (PR #39). What shipped is in `CHANGELOG.md`.

**Done:** `core/fit_runner.py` (one runner per fitting tool) and
`core/fit_data.py` (one point/σ rule). `pyirena.batch`, `api/control` and
`analyze()` call them; the six GUI panels apply the same data rule. A snapshot
of every batch/`analyze` number on the six `testData/Scripting` fixtures was
identical before and after each tool's change.

**Kept here** because the code cites the decisions (D1–D4) below.

## Decisions (Jan, 10-11-2026)

| # | Question | Decision |
|---|---|---|
| D1 | Unified: walk limits? | **Walk** (`fit_with_limit_walking`) for bounded fits on every path — GUI limits are auto-windows reset by *Fix limits* after each fit. With `no_limits`: no walking (control's behaviour). The core docstring saying batch should not walk is corrected. |
| D2 | Points with σ ≤ 0 | **Repair** to σ = I × the Data Selector "Generated uncertainty fraction" (the text-import rule), on every path. Today: clamped to 1e-30 (one point dominates χ²) in GUI/control, dropped in Sizes/Simple batch, repaired at 5 % in `analyze`. |
| D3 | Points with I ≤ 0 | **Keep** in fits on every path. Text import's own cleaning is unchanged. `analyze` (array sessions) and the Simple Fits GUI stop dropping them. (Core math that needs I > 0, e.g. Sizes' internal mask, is unchanged.) |
| D4 | Sizes background pre-fits on replay | **Re-run when `fit_power_law_B/P` are checked** (as batch does); an empty power-law Q window falls back to the **cursor range** like the GUI, not the full curve. |
| — | Data without uncertainties | (Decided earlier today.) GUI/batch: σ = I × generated fraction for MC; API refuses MC (no setting on a server). |
| — | Re-conversion of text files | (Decided earlier today.) Back up + recreate; GUI asks. |

Defaults taken without a separate question (say if any is wrong):

- **WAXS peaks outside the fitted Q range** are removed before the fit on
  every path (GUI behaviour), and the runner reports which.
- **WAXS weight mode** missing from a config stays `standard` (the documented
  API default); only old configs are affected.
- **GUI writes `load_slit_smeared`** into exported/embedded state when the
  slit-smeared dataset is loaded, so batch fits the same data (today batch
  loads the desmeared dataset and smears the model against it). New optional
  key; old files still load.
- **Every bounded tool honours `no_limits`** from the config in batch
  (Simple Fits and Carbon batch ignored it).
- **Envelope `fit_q_range` / `unwrap_config`** used by batch too (AGENTS.md
  invariant 6); batch Unified/Simple/WAXS/Carbon bypassed it.
- **Session settings stored on the model win over argument defaults** in
  control (Carbon `weighting`/`n_mc_runs`/Q range). WAXS `no_limits` keeps
  its `False` default — `test_no_limits.py` makes "honour the bounds unless
  asked" a rule for every bounded tool, and `analyze` passes a config's own
  value explicitly.

## Remaining — the GUI step (`FitPanelBase`)

The panels still orchestrate their own fits (load → Q mask → slit settings →
fit → MC → store), next to the runners. Open work, in order:

1. **`FitPanelBase`** (review §4.2). The six fit panels each re-define the
   same ~20 methods — `set_data`, `plot_data`, `plot_residuals`,
   `save_state`/`load_state`, `_load_setup_from_file`,
   `_data_info_for_report`, `get/set_cursor_range`,
   `_reload_data_with_smearing`, `_get_data_folder`, `closeEvent`. Those are
   logic, not widget assembly, so they drift (the reason PLAN.md gives for
   *not* sharing the file browsers cuts the other way here). The base owns
   data loading, Q range, state and setup I/O, and calls the tool's runner,
   so the panel only maps widgets ↔ model. `gui/fit_worker.py`
   (`FitJobRunner`) is already shared and is the threading half of it.
2. **WAXS Peak Fit GUI: stored results recompute the background over the
   full Q range** instead of the fitted range, so a re-displayed stored
   result does not match the fit. A GUI-only bug found by the audit; fix it
   in the move (or before, if a user hits it).
3. **Large-panel split** (`PLAN.md`, "Large-file comfort refactors") —
   `FitPanelBase` is its natural seam; do it after, tool by tool.

**Guard:** the snapshot script covers batch/`analyze` only. For the GUI step,
compare stored HDF5 results from each panel before/after on the same fixtures
(offscreen), and keep `validationData/run_validation_report.py` at 195/195.
