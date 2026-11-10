# One fit runner per tool — audit and decisions (issue #38, review §3.2)

Started 10-11-2026 on `fix/review-2026-11-issue38`.

## Why

Fit orchestration — data preparation → Q mask → slit settings → fit →
Monte-Carlo → result — existed three times per tool: the GUI panel,
`pyirena/batch/<tool>.py` and `pyirena/api/control/<tool>.py` (plus the
`analyze()` replay). A read-only audit of all six fitting tools (10-11-2026)
found ~35 places where the copies do something different. Most are plain bugs;
some were genuine choices, decided below.

**Scope of this step:** batch and api/control (incl. `analyze`) call one
runner per tool. The GUI keeps its own orchestration for now (a later
`FitPanelBase` step), but the *data rules* below apply to the GUI too.

**Guard:** `validationData/run_validation_report.py` calls the core fits
directly, so it does not cover orchestration. A snapshot of every number batch
`fit_*` and `analyze()` produce on the six `testData/Scripting` fixtures is
compared before/after each tool's change; it must be identical except where a
decision below changes a number deliberately.

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
  control (Carbon `weighting`/`n_mc_runs`/Q range, WAXS `no_limits`).

## Audit — real differences, by tool

Unified Fit: batch never walks (D1); batch slit-length rule differs from
`apply_data_slit_settings` (config length wins even when smearing is off);
SMR dataset not reloaded in batch (GUI never writes `load_slit_smeared`); MC
refit windows differ by path; batch needs both cursors, `analyze` accepts one.

Simple Fits: batch ignores `no_limits`; points dropped differ on every path
(D2/D3); control skips the Invariant background prefit; batch MC run count
ignores the config.

Size Distribution: control never slit-smears; control MC inversion runs 10
repetitions (core default) where GUI/batch/analyze use 1; prefit trigger and
empty-window fallback (D4); σ ≤ 0 dropped vs repaired (D2); `tool_config`
ignores `montecarlo_convergence` / `montecarlo_max_iter`, maps no `'mcsas'`
alias, and defaults max iterations to 300 vs 1000; batch saved quality
metrics add the unsmeared background.

Modeling: batch builds the config by hand (`_build_pop`) instead of
`modeling_model_from_config` (missing `enabled` → disabled; missing q range
→ data range); σ for data without errors is 5 % in control/analyze vs the
setting in GUI/batch; σ ≤ 0 floored vs repaired (D2); Export Parameters omits
the slit keys; batch MC run count ignores the config; DE unseeded everywhere.

WAXS Peak Fit: control never runs the Q0 pre-search; out-of-range peaks pruned
only in the GUI; batch crashes on files without uncertainties; batch ignores
envelope/panel Q range; control ignores the model's `no_limits`; with no peaks
GUI fits background, batch auto-finds, control/analyze refuse; GUI *stored*
results recompute the background over the full range (GUI bug — for the
`FitPanelBase` step).

Carbon model: batch ignores `no_limits`; control overwrites the model's
`weighting`/`n_mc_runs`/Q range with argument defaults; batch reads Q range
from the model dict only; a one-sided range from a bare model dict can be
lost in `analyze` (q_max = 0 meaning "none"); batch runs MC whenever the
config's `n_mc_runs > 1`.

## Status

- [ ] Data rules D2/D3 in one shared helper, used by batch, control, analyze and the GUI panels
- [ ] Runner: ~~Unified~~ · ~~Sizes~~ (done) · Simple Fits · Modeling · WAXS · Carbon
- [ ] GUI writes `load_slit_smeared`
- [ ] Snapshot comparison clean (or every change explained) per tool
