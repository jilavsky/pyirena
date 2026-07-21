# Code Review — Slit Smearing Implementation (2026-07-20)

Reviewed: uncommitted working tree on `feature/slit-smearing` (30 modified
files, 4 new; +1183/−114 lines). Focus: logic and implementation correctness;
tests confirmed green locally (test_smearing 16/16, plus unified/local-fits/
sizes/simple-fits/roundtrip suites all pass).

## Verdict

The implementation is faithful to the plan where it matters most: the core
engine is mathematically correct (independently re-verified below), the
SL≤0-is-a-no-op contract is applied uniformly, provenance attrs + dual model
curves are written consistently by Unified/Sizes/Simple Fits, detection is
attribute-driven (dQl / `resolutions="dQw,dQl"`) rather than name-driven, and
batch enforcement fails loudly. One deliberate deviation from the plan
(smearing local-fit models on the *sub-range* grid instead of the full-grid
Igor pattern) is actually an improvement — with an analytic model plus
`build_extended_q`, sub-range smearing is exact, and E2 degrades to a worded
warning as intended.

Independent checks performed during review:

* `W` row sums are exactly 1 (constants preserved regardless of interpolation).
* Smeared power-law slopes: −(P−1) for q≪SL, −P for q≫SL (test 3 asserts this).
* Real-file check: `smear_curve(DSM)` reproduces the SMR entry, median
  |rel.dev.| ≈ 0.1% mid-range (test 7 asserts <1%).
* Synthetic Unified fit on smeared data recovers ideal G/Rg/P and is
  non-circular (data made via `smear_curve`, fit uses the matrix path).

## Findings (ordered by priority)

### F1 — Background prefits are not smeared (correctness inconsistency)

`sizes.py::fit_power_law` / `fit_background_term` and
`simple_fits.py::prefit_background` fit the **ideal** `B·q^-P` directly to
slit-smeared data. The main fit then smears `compute_complex_background`
(sizes) — so the background actually subtracted is `smear(B'·q^-P')` where
`B',P'` already describe the *smeared* data: effectively double-smeared.
Simple Fits has the mirror problem (prefit unsmeared, main model smeared).
Magnitude is O((SL/q)²) — small when the bg window sits at q≫SL, but for
USAXS SL≈0.018 with windows at 0.05–0.1 it is percent-level and systematic.
**Fix:** wrap the prefit model functions with `smear_model` when
`use_slit_smearing` is active (same one-line pattern used elsewhere), so
prefit returns ideal-space B/P consistent with the main fit.

### F2 — Data Manipulation output can carry a stale SMR twin (data-integrity trap)

`save_manipulated_data` copies the input file and replaces `Q/I/Idev/Qdev` in
the default sasdata group. For a Matilda file the copied output still contains
the **original, unmanipulated `_SMR` entry** next to the manipulated desmeared
curve. Anything that later loads that output with `prefer_slit_smeared=True`
(GUI checkbox, batch `load_slit_smeared`) gets data inconsistent with the
file's default entry — silently. The PLAN status note ("provenance dQl is
preserved through the copy-based save") is true but is precisely the problem:
it preserves a now-wrong sibling.
**Fix options:** (a) drop `_SMR` entries from manipulation outputs, or
(b) apply the same operation to both entries where meaningful (trim/scale),
dropping otherwise. (a) is safe and one line.
Related: the new subtract/divide guard lives only in the GUI panel
(`_slit_compatible`); `core/data_manipulation.py` and `batch/manipulate.py`
have no equivalent, and since the manipulation panel never loads SMR data
today, the guard is currently unreachable. Move the check into the core
`subtract`/`divide` so every caller inherits it.

### F3 — Data Merge scope was deferred (needs Jan's sign-off)

Plan §5.6 (merge slit-smeared USAXS with pinhole SAXS, provenance-labelled
output) is not implemented — no changes in `core/data_merge.py`,
`gui/data_merge_panel.py`, `batch/merge.py`, or `io/nxcansas_data_merge.py`.
The PLAN status records this as "deferred (utility-tier)". That was an
explicit user-requested capability, so the deferral is a decision to confirm,
not a bug — flagging so it doesn't get lost.

### F4 — Modeling doesn't save the ideal model curve

`save_modeling_results` writes only the `slit_length` /
`data_is_slit_smeared` attrs; `model_I` (and per-population curves) are the
smeared ones, with no `*_ideal` datasets — inconsistent with the
save-both-curves convention honoured by Unified, Sizes, and Simple Fits.
`ModelingEngine.total_intensity_maybe_smeared` already computes `I_ext`
(ideal on `q_ext`); truncating its data-grid nodes to get the ideal curve is
nearly free.

### F5 — UnifiedFitPanel duplicates the slit UI instead of using the mixin

`slit_smearing_ui.SlitSmearingMixin` is used by Sizes/Simple Fits/Modeling,
but `UnifiedFitPanel` carries its own ~110-line copy of the same row,
handlers, and reload logic (only real differences: `_sync_smearing_to_model`
and state-restore). Divergence risk on the next UI change; worth folding into
the mixin with a `_sync` hook.

### F6 — Per-iteration W rebuild in Simple Fits (performance)

The smeared `base_func` closure in `SimpleFitModel.fit` calls `smear_model`
on every objective evaluation, rebuilding the extended grid and sparse W each
time (q is constant across the fit). With numeric Jacobians this is hundreds
of rebuilds — likely seconds of avoidable cost. Hold a `SlitSmearer` (as
Unified and Modeling do). Same theme in Modeling: size-distribution
populations now build their G matrices on the ~8×-refined `q_ext` every
iteration under `use_cache=False`; watch DE runtimes with size-dist
populations, and consider a lower `refine` for that path if it bites.

### F7 — Fractals smearing uses the tabulated path everywhere + silent except

`_smear_for_display` applies `smear_curve` (power-law extrapolation) to both
the MC curve (correct — it is tabulated) and the analytic unified curve
(where `smear_model` would be exact, no extrapolation). It also swallows all
exceptions and returns the **unsmeared** curve — a silent wrong-overlay
failure mode; at minimum log it. The plan's suggestion to extend the MC q
grid to √(qmax²+SL²) (avoiding extrapolation entirely) was not taken up.

### F8 — E3 not implemented

`check_level_feasibility` and `detect_features` have no smearing awareness:
no warning when fitted Rg ≫ π/SL, and feature detection still assumes pinhole
peak/knee shapes on smeared data. Plan §5.1 called for feasibility notes in
both. Small, self-contained follow-up.

### F9 — Integration-test gaps vs plan §7

Engine tests (1–7) are strong. Missing: (8) SMR-entry vs DSM-entry fit of the
real test file agreeing on Rg/G/B/P; (9) round-trip asserting
`intensity_model_ideal` + `slit_length` survive save→load for each tool;
(11/12) batch `load_slit_smeared` and control-API
`open_dataset(use_slit_smeared=True)` smoke tests. These are the tests that
would have caught F1/F2-class regressions.

### F10 — Nits

* `smear_curve` drops non-positive intensity points before log-log
  interpolation — biases the result upward near the noise floor of subtracted
  data; fine for its current uses, deserves a docstring warning.
* Mixin: if `_reload_data_with_smearing` fails, the checkbox stays toggled
  while the old dataset remains loaded (UI/state mismatch).
* `modeling_panel._reload_data_with_smearing` pokes
  `self.data_loader._edit.text()` (private attr).
* CHANGELOG.md not updated (plan phase 8 item).

## Suggested order of fixes

1. F1 (prefit smearing) + F2(a) (drop stale _SMR from manipulation output) +
   core-level guard — correctness, ~1 short session.
2. F9 round-trip/integration tests — locks in the contract before refactors.
3. F4, F7, F8 — small consistency items.
4. F5, F6 — refactor/perf, no behavior change.
5. F3 — separate decision with Jan (implement merge support or formally
   descope it in the plan).

## Resolution (2026-07-21)

All findings addressed on `feature/slit-smearing` (see CHANGELOG "Unreleased"):

- **F1** — Sizes (`fit_power_law`/`fit_background_term`) and Simple Fits
  (`prefit_background` + the `saxs_morph` bg helpers, via a new `slit_length`
  kwarg) now smear the prefit model; verified to recover ideal B/P on smeared
  synthetic data.
- **F2** — `save_manipulated_data`/`save_merged_data` drop stale `_SMR` twins;
  `replace_nxcansas_data` clears an orphaned `dQl`; the compatibility guard now
  lives in `DataManipulation.check_slit_compatible` (core), inherited by batch.
- **F3** — Data Merge implemented (propagate-only per Q3): core `MergeConfig`/
  `MergeResult` + `merged_slit_length`, GUI/batch wiring, `append_dql` writes
  `dQl` to the merged output.
- **F4** — Modeling saves `model_I_ideal` (+ read-back); a q-agnostic G-cache
  bug exposed by this was fixed (`total_intensity_ideal` uses the q-keyed cache).
- **F5** — `UnifiedFitPanel` folded into `SlitSmearingMixin` via an optional
  `_sync_smearing_to_model` hook.
- **F6** — Simple Fits holds one `SlitSmearer` across the fit (no per-iteration
  `W` rebuild); `compute()` reuses the same one-instance pattern.
- **F7** — Fractals `_smear_for_display` logs + flags failures instead of
  silently returning the unsmeared curve.
- **F8/E3** — `UnifiedLevel.slit_smearing_note` and `detect_features(slit_length=)`
  add soft warnings (GUI feasibility label + feature-identifier summary).
- **F9** — `pyirena/tests/test_smearing_integration.py` covers display smearing,
  prefit recovery, merge round-trip, manipulation guard, `_SMR` drop, and
  Simple-Fits/Modeling save→load of slit provenance.
- **F10** — `smear_curve` docstring warns about the non-positive-point drop;
  mixin reverts the checkbox on a failed reload; `modeling_panel` no longer
  pokes `data_loader._edit`; CHANGELOG updated.

**Jan's testing bugs:** the Simple Fits display bug (root cause: `compute()` did
not smear) is fixed; Guinier/Porod linearization is now labelled best-effort for
smeared data.

## additional issues observed by Jan in testing
Simple fits tool:
Sphere model does not show Bessel function oscillations smeared out in display when lit smearing is on. That is impossible - for example - for R=600 the first oscillation is around Q~0.007, slit length is 0.018, first oscillation should be smeared out strongly. Looks perfectly sharp. Is the model smeared? 

Spheroid, Troubler Strey look also same with and without smearing

Guinier fit - how is this handled? Linearization is probably wrong for smeared data, is there linearization for slit smeared data? Probably not, so this may just stay problematic and display incorrectly. This was one of teh reasons not all simple fits in Igor supported slit smeared data, w shoudl cal this best effort.  

Porod - the fit looks good. Linearization is same challenge as Guinier - theoretically Porod linearization is I*Q^3 vs Q^3, but only for infinite slit length. Our data have Porod slope transition below my tested Q range.   
