# Fit Quality Metrics — Implementation Plan

**Branch:** `feature/fit-quality-metrics`
**Spec:** [`04-fit-quality-metrics.md`](04-fit-quality-metrics.md) (the "what/why")
**This doc:** the "how" — phased, backward-compatible, with decision gates.
**Status:** Phases 1–3 implemented (2026-06-18) — at **Decision Gate 1**
(upstream-AI MCP evaluation). Phases 4–5 pending gate outcome.
**Date:** 2026-06-18

> **Progress (2026-06-18):**
> - Phase 1 ✅ `pyirena/core/fit_metrics.py` + `tests/test_fit_metrics.py` (10 tests pass).
> - Phase 2 ✅ `fit_quality_metrics` exported from `pyirena`.
> - Phase 3 ✅ MCP/control: new `get_fit_quality` tool (control fn + `TOOL_SCHEMAS`
>   + `pyirena_ctrl_get_fit_quality`); `run_fit` result gains a `quality` block;
>   `get_residuals` gains `rescaled_residual`, `frac_misfit_percent`, and
>   `summary.robust_scale_s`. Docs updated in `docs/ai_tools_reference.md`.
>   All additive — no existing keys changed. Full suite: 171 pass (2 pre-existing
>   unrelated failures: unified invariant NaN, saxs_morph Gibbs).
> - **→ Hand to upstream AI test bed.** Then Decision Gate 1 before Phase 4.

---

## 0. Executive summary

The spec in `04-fit-quality-metrics.md` is technically sound and the design
(separate σ *scale* from σ *shape*; use MAD for a robust noise scale; use the
σ-independent fractional misfit as a gross-misfit backstop; band by Q; analyze
sign structure) is the right approach for SAXS data where reported σ are
frequently mis-scaled.

The implementation risk is **low** because the spec is deliberately additive:
one new pure-array function, new result keys, new GUI views — nothing existing is
removed or changed. The work is sequenced so we can **prove value on the MCP path
the upstream AI uses, and on Unified Fit, before committing effort to the other
five fit tools.**

The plan has a hard **decision gate after Phase 3** (MCP) and a soft one after
Phase 4 (Unified Fit GUI). Propagation to the other tools (Phase 5) is explicitly
conditional on those results.

---

## 1. What I verified in the current code (integration surface)

| Concern | Finding | File:line |
|---|---|---|
| Core fit returns | `chi_squared`, `reduced_chi_squared`, `fit_intensity`, raw `residuals` (I−M) | `core/unified.py:472-481` |
| Weighted residuals used internally | `(I−M)/σ` with weights=1 when σ missing | `core/unified.py:394-411` |
| Sign-run / CorMap already exists | `_longest_run()`, `P_longest_run()`, full CorMap | `core/similarity.py:50,75,150` |
| Public API exports | `__all__` list, no metrics fn yet | `pyirena/__init__.py:51-71` |
| MCP control: chi² tool | `get_chi_squared()` | `api/control/unified_fit.py:1430` |
| MCP control: residuals tool | `get_residuals()` returns normalized resid + rms/max_abs/mean | `api/control/unified_fit.py:1443-1471` |
| MCP run_fit result | dict with chi², reduced chi², iterations, params | `api/control/unified_fit.py:1415-1423` |
| Array serialization | `_decimated()` / `array_to_list()` → list, NaN/inf→None, decimate to `PYIRENA_MAX_ARRAY_POINTS` (500) | `api/control/unified_fit.py:122`, `api/control/schemas.py:60-98` |
| TOOL_SCHEMAS | list of `{name, description, input_schema}` | `api/control/schemas.py:13-622` |
| MCP server | FastMCP thin pass-through, auto-serializes dicts | `mcp/server.py` |
| Matplotlib residual panel | normalized if σ else raw | `plotting/unified_plots.py:78-92` |
| GUI residual plot | `plot_residuals(q, resid)`, generic "Residuals" label | `gui/unified_fit.py:1008-1024,586-622` |
| GUI residual computation | `(I−M)/σ`, falls back to `(I−M)/I` | `gui/unified_fit.py:2747-2753` |
| Tests | pytest, `pyirena/tests/test_*.py`, synthetic-data pattern | `tests/test_unified.py:86-123` |

**Other fit tools (Phase 5 targets), all return chi²/residuals in their own way:**
`core/sizes.py:222-226` (regularized — see caveat §3.3), `core/simple_fits.py:633`,
`core/modeling.py` (`ModelingResult` dataclass), `core/waxs_peakfit.py:1002-1006`.

---

## 2. Critical evaluation of the spec (before we build)

The spec is good. These are refinements / risks to fold into the implementation —
none change the core design.

1. **`(I−M)/I` is fragile when I → 0 or I < 0** (post background-subtraction,
   or noisy high-Q). The fractional misfit and `frac_uncertainty` can explode or
   flip sign. → **Mask/guard:** compute fractional arrays only where `I` is finite
   and `|I|` exceeds a small floor relative to its own scale; otherwise emit `NaN`
   for that point and exclude it from `max_abs_frac_misfit`. Document this.

2. **Regularized fits (`sizes`, MaxEnt/TNNLS) have no clean `n_params`.** The spec
   takes `n_params` as input, so `reduced_chi2` is the caller's responsibility.
   Worth noting that the *robust* metrics (`robust_scale_s`, fractional misfit,
   sign runs) **do not need `n_params`** and are arguably *more* meaningful than
   reduced χ² for those fits. This strengthens the case for the robust metrics,
   and means Phase 5 for sizes is valuable even though its reduced χ² is murky.

3. **Reuse `similarity.py`, don't reinvent.** `_longest_run()` already computes the
   longest same-sign run, and `P_longest_run()` gives the CorMap probability the
   Unified GUI already displays. `fit_metrics.py` should call these (or a shared
   helper) for `longest_same_sign_run` rather than duplicating logic. Add only the
   lag-1 sign autocorrelation, which is new.

4. **Banding needs a minimum-points guard.** With few valid points, 4 log-spaced
   bands can have 0–1 points → meaningless per-band χ²/MAD. → Adaptively reduce
   `n_bands` so each band has ≥ a small minimum (e.g. 5); report the actual band
   count used.

5. **MAD center.** Use `median(r_i)` as the center (spec already does). For a fit
   with a systematic offset this keeps `s` honest; the offset itself shows up in
   the sign-run / mean diagnostics.

6. **Determinism.** All operations (median, MAD, sums) are deterministic; no RNG.
   Good — matches "cheap, deterministic, computed once."

7. **Naming for the consumer.** Keep the spec's field names verbatim — the
   pyirena-ai decision tree (spec §6) is written against them.

**Verdict:** proceed. Build the function exactly as specified, with the guards in
points 1 and 4, and reuse `similarity.py` per point 3.

---

## 3. The function contract (Phase 1 deliverable)

`pyirena/core/fit_metrics.py`, one public function, pure (arrays in → dict out),
**no thresholds, no verdicts** (per spec §3.1).

```python
def fit_quality_metrics(
    q, intensity, model, sigma, n_params, n_bands=4
) -> dict
```

### 3.1 Returned dict (stable schema — this is the contract consumers depend on)

```
sigma_available: bool
n_valid: int
n_params: int
dof: int | None

# per-point arrays (aligned to the *valid-point* q, which is also returned)
q_valid:            np.ndarray
norm_residual:      np.ndarray | None    # (I−M)/σ ; None if no σ
frac_residual:      np.ndarray           # (I−M)/I, NaN where |I| too small
frac_uncertainty:   np.ndarray | None    # σ/I

# global scalars
reduced_chi2:                  float | None
robust_scale_s:                float | None   # 1.4826·MAD(norm_residual)
sigma_misscale_factor:         float | None   # alias of robust_scale_s
realistic_reduced_chi2_floor:  float | None   # s²
median_frac_uncertainty:       float
max_abs_frac_misfit:           float
q_at_max_frac_misfit:          float
n_outliers_3s:                 int | None     # |norm_residual| > 3s
frac_outliers_3s:              float | None

# structure scalars
longest_same_sign_run:  int          # via similarity._longest_run on (I−M)
sign_autocorr_lag1:     float        # lag-1 autocorr of sign(I−M)

# per-band (adaptive count; each band ≥ min points)
n_bands_used: int
bands: list[dict]   # {q_lo,q_hi,n,reduced_chi2,robust_scale_s,max_abs_frac_misfit}
```

### 3.2 Behavior contract (locks the edge cases for tests)
- σ None / all ≤0 / non-finite → `sigma_available=False`; `norm_residual`,
  `robust_scale_s`, `reduced_chi2`, outliers, `frac_uncertainty` are `None`;
  fractional arrays still populated; **no crash, no invented σ**.
- `dof ≤ 0` → `reduced_chi2=None` (don't divide).
- Points with σ≤0 or non-finite I/M/q/σ are masked out before any stat.
- `frac_residual` NaN where `|I|` below floor; excluded from `max_abs_frac_misfit`.
- Adaptive `n_bands_used ≤ n_bands`.

### 3.3 Caller responsibility note
`n_params` is supplied by each fit. For regularized fits where it is ill-defined,
the caller passes its best estimate (or the value it already uses for its own
reduced χ²); the robust metrics remain valid regardless.

---

## 4. Phased plan

Each phase is independently shippable and backward-compatible. **Existing keys
(`chi_squared`, `reduced_chi_squared`, `residuals`) are never removed or changed.**

### Phase 1 — Core function + unit tests  *(foundation, zero integration risk)*
- Add `pyirena/core/fit_metrics.py` with `fit_quality_metrics()` per §3.
- Reuse `core/similarity._longest_run`; add lag-1 sign autocorr helper.
- Add `pyirena/tests/test_fit_metrics.py` covering the spec §5 cases:
  well-calibrated σ (`s≈1`), σ underestimated ×3 (`s≈3`, floor≈9, **no** outliers),
  localized misfit (outliers>0, hot band, `q_at_max_frac_misfit` inside band),
  systematic wrong-slope (long run, high autocorr), σ missing (`sigma_available
  False`, no crash), dof≤0 guard, and the I→0 fractional guard.
- **Risk:** none — nothing imports it yet.
- **Done when:** `pytest pyirena/tests/test_fit_metrics.py` green.

### Phase 2 — Public API export  *(make it callable everywhere)*
- Export `fit_quality_metrics` from `pyirena/__init__.py` (`__all__` + import).
- One line in the package docstring/example.
- **Risk:** none (pure addition to `__all__`).

### Phase 3 — MCP / control layer  ★ primary test bed (upstream AI) ★
Additive only. Three sub-steps, each shippable:
- **3a.** New control fn `get_fit_quality(session_id)` in
  `api/control/unified_fit.py` → calls `fit_quality_metrics` on the session's last
  fit, serializes via `_decimated()` for arrays. New MCP tool
  `pyirena_ctrl_get_fit_quality` + new `TOOL_SCHEMAS` entry (`schemas.py`).
- **3b.** Enrich `run_fit` result dict with the decision-relevant scalars
  (`robust_scale_s`, `realistic_reduced_chi2_floor`, `max_abs_frac_misfit`,
  `q_at_max_frac_misfit`, `n_outliers_3s`, `longest_same_sign_run`) — additive keys.
- **3c.** Enrich `get_residuals` summary with `robust_scale_s` and offer the
  rescaled `r' = r/s` and fractional-% arrays **as additional keys** (keep
  `residuals` exactly as-is).
- Update `docs/ai_tools_reference.md` (run_fit / get_residuals / new tool +
  interpretation guidance) and `docs/ai_integration.md` tool count.
- **Risk:** low — only new keys/tool. Existing agent calls unaffected.
- **Verify:** unit-test the control fns against a synthetic session; smoke-test
  with MCP inspector. Then hand to upstream AI for real-world evaluation.

> **★ DECISION GATE 1 (after Phase 3).** Upstream AI evaluates whether these
> metrics improve fit-quality judgments on real data. **If no clear benefit →
> stop here; revert nothing (it's additive) but don't propagate.** If beneficial →
> proceed to Phase 4 and consider Phase 5.

### Phase 4 — Unified Fit (GUI + matplotlib)  ★ human test bed ★
- **Summary line:** alongside reduced χ², show `robust_scale_s` and
  `max_abs_frac_misfit` (% with its Q) in the fit success message / status
  (`gui/unified_fit.py:2928-2946`).
- **Selectable residual view:** add a small selector (normalized `r` /
  rescaled `r'=r/s` / fractional misfit %) to the residuals panel. **Default stays
  the current normalized view** initially so nothing changes for existing users;
  switching the default to `r'` is a follow-up decision after we see it in use.
  Touches `gui/unified_fit.py:1008-1024` (plot) and `:586-622` (panel) + the
  matplotlib `plotting/unified_plots.py:78-92`.
- Annotate the plot with `robust_scale_s` and `realistic_reduced_chi2_floor`
  ("σ look ~3× low; χ²≈9 is the real floor").
- **Risk:** low-moderate (GUI only). The view is a *toggle* — the underlying data
  and fit are untouched; worst case a display bug, not a data bug.
- **Verify:** `/run` Unified Fit on a test dataset, switch views, confirm
  annotations and that the default view is unchanged.

> **DECISION GATE 2 (after Phase 4).** Does the rescaled/fractional view help
> users judge fits better than normalized residuals? Decide (a) whether to make
> `r'` the default view, and (b) whether to propagate to other tools.

### Phase 5 — Propagate to other fit tools  *(conditional on gates 1 & 2)*
Apply the same additive pattern, **one tool per commit**, in value order:
1. **sizes** (`core/sizes.py`) — robust metrics especially valuable (no clean dof);
   + `gui/sizes_panel.py` residual view.
2. **simple_fits** (`core/simple_fits.py`) + `gui/simple_fits_panel.py`.
3. **modeling** (`core/modeling.py` `ModelingResult`).
4. **waxs_peakfit** (`core/waxs_peakfit.py`).
Each: add `results['quality'] = fit_quality_metrics(...)` (or dataclass field),
mirror the GUI summary line + selectable residual view, extend the read-side API
schemas (`api/.../results.py`, `schemas.py`) and docs. Add per-tool tests.
- **Risk:** low per tool, but multiplied across four tools + GUIs — hence gated.

---

## 5. Backward-compatibility & rollback

- **Additive contract:** no existing dict key, function signature, or default GUI
  behavior changes through Phase 4. Old code/agents keep working untouched.
- **Default-preserving GUI:** new residual views are opt-in selectors; the default
  remains today's normalized residual until an explicit gate-2 decision.
- **Rollback:** each phase is its own commit(s). Reverting Phase N never affects
  Phase <N. Phase 1–3 carry essentially no user-visible risk.
- **Schema stability:** the §3.1 dict is the consumer contract — once Phase 3
  ships to the agent, add fields but don't rename/remove (mirrors the
  STATE_SCHEMA_VERSION discipline used elsewhere in the project).

---

## 6. Open questions for you

1. **Default residual view (gate 2):** keep normalized as default and offer `r'`
   as a choice, or eventually flip the default to `r'`? (Plan assumes keep-default,
   decide later.)
2. **Phase 3c scope:** do you want the rescaled/fractional *arrays* in
   `get_residuals` from the start, or just the `robust_scale_s` scalar first?
3. **Phase 5 priority:** is sizes the right first propagation target (I argue yes —
   robust metrics matter most where dof is ill-defined), or a different tool?

---

## 7. Suggested commit sequence

```
1. docs(planning): fit-quality-metrics implementation plan        [this file]
2. feat(core): fit_quality_metrics + tests                        [Phase 1]
3. feat(api): export fit_quality_metrics                          [Phase 2]
4. feat(api/control): get_fit_quality tool + enrich run_fit/get_residuals + docs   [Phase 3]
   --- DECISION GATE 1: upstream AI evaluation ---
5. feat(gui): unified-fit quality summary + selectable residual view               [Phase 4]
   --- DECISION GATE 2: keep/flip default, decide propagation ---
6..N. feat(core+gui): propagate to sizes / simple_fits / modeling / waxs_peakfit   [Phase 5]
```
