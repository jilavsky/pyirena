# Fit Quality Metrics — robust residual diagnostics

**Status:** proposed (requested by pyirena-ai side)
**Owner:** pyirena core
**Consumers:** pyirena-ai agent, in-GUI AI advisor, MCP server, and end users

## 1. Motivation / the problem

Fit quality is currently judged by `reduced_chi_squared` plus raw `residuals`
(`I − M`). That breaks down in practice because **reported uncertainties σ are
not a reliable measure of the actual noise**:

- σ scale is often wrong (over- or under-estimated), so a target of
  reduced χ² ≈ 1 may be physically unreachable.
- Consequence A (over-fitting): an agent or user chases χ² < 1 forever, wasting
  effort on a target the data can't support.
- Consequence B (missing gross misfits): when told "σ are unreliable, don't
  worry about exact χ²", the metric gets ignored and genuine misfits slip
  through. We have seen reasoning like *"a normalized residual of 20–50 is fine
  because σ are unreliable"* — which is wrong, see §2.

We want diagnostics that are **robust to a mis-scaled σ** yet still **catch real
misfits**, and that are cheap, deterministic, and computed once per fit so the
consumer (agent/advisor/user) just reads numbers instead of re-deriving them
from the raw arrays.

## 2. The key idea — separate the *scale* of σ from its *shape*

Reported σ carry two independent pieces of information that fail differently:

1. **Relative shape** (which points are noisier than which — high-Q noisier,
   low-Q depends on signal/background): usually trustworthy.
2. **Absolute scale** (overall noise level): often unreliable; this is what
   breaks the χ² ≈ 1 target.

So: keep σ for *relative weighting*, but let the residuals themselves recover
the *true scale*, then judge against that.

### 2a. Fractional misfit — the σ-independent anchor

Definitions per point i:

```
r_i = (I_i − M_i) / σ_i          normalized residual
f_i = σ_i / I_i                  fractional uncertainty ("a few %")

(I_i − M_i) / I_i  =  r_i · f_i   fractional misfit
```

So **fractional misfit = normalized residual × fractional uncertainty.**

This is why "r = 20–50 is fine" is wrong: with f ≈ 5 %, r = 20 ⇒ the model is
off by 100 % of the data value (M ≈ 0). With f ≈ 2 %, r = 50 ⇒ also 100 % off.
A point with |(I−M)/I| ≳ 30–50 % **cannot** be within few-% noise no matter how
unreliable σ is. The fractional misfit is a σ-independent backstop for gross
misfits.

### 2b. Robust noise scale — self-calibration (guards against over-fitting)

Compute a robust scale of the normalized residuals using the median absolute
deviation:

```
s = 1.4826 · median( | r_i − median(r_i) | )      # MAD → Gaussian-σ estimate
```

Interpretation:

- `s ≈ 1` → σ are well-calibrated; reduced χ² ≈ 1 is an honest target.
- `s ≈ 3` → σ underestimated ~3×; the *realistic* reduced χ² floor is ~9.
  Chasing χ² < 1 is futile — stop.
- `s ≈ 0.3` → σ overestimated; over-fitting risk.

**Why MAD and not plain reduced χ²:** MAD reflects the *bulk* of points and
ignores a handful of outliers. That is exactly what separates the two failure
modes:

- σ merely mis-scaled → the *whole* distribution of r_i is wide → `s` is large
  and **no points sit far above s**.
- Real misfit → the bulk r_i stay small (`s ≈ 1–3`) but a few points sit at
  r = 20–50 ≫ s → flagged as genuine outliers.

Then rescale: `r'_i = r_i / s`. A residual of 20–50 in a dataset whose bulk
scatter is s ≈ 2 gives r' = 10–25 — unambiguously a gross misfit. The classic
mistake is comparing residuals to reported σ (or to aggregate χ²) instead of to
this robust scale.

### 2c. Structure — scale problem vs model problem

Magnitude alone is not enough. **Structure** in the residual sequence
distinguishes a bad-σ-scale (leave it alone) from an inadequate model (fix it):

- **Sign runs / autocorrelation:** long stretches of same-sign residuals ⇒ wrong
  functional form (systematic), even at modest magnitude.
- **Localized clusters** of large r'_i over a Q-feature ⇒ missing/mis-shaped
  component.
- **Random scatter, no runs, just wide** ⇒ scale problem only; do not "fix" it.

### 2d. Respect Q-dependence — band it

Noise genuinely varies with Q, so a global MAD can be dominated by one regime.
Compute the scale and χ² in 3–4 Q-windows (≈ per decade of Q) as well as
globally. Uneven per-decade χ² is itself a misfit signal.

## 3. What to implement

A single, **model-agnostic, standalone** function — it operates only on arrays,
so it serves unified, sizes, simple_fits, modeling, waxs_peakfit, etc.

### 3.1 New module `pyirena/core/fit_metrics.py`

```python
def fit_quality_metrics(
    q: np.ndarray,
    intensity: np.ndarray,        # I (data)
    model: np.ndarray,            # M (fit_intensity)
    sigma: np.ndarray | None,     # error_data; may be None or contain <=0
    n_params: int,
    n_bands: int = 4,             # Q-decade-ish windows
) -> dict:
    ...
```

Computation steps:

1. **Validate / mask.** Drop points where `sigma <= 0` or non-finite, and where
   `I/M/q` are non-finite. If `sigma` is `None` or all invalid, set a flag
   `sigma_available = False` and fall back to using `I` as the scale for the
   normalized residual (i.e. report fractional residuals only; `s`, χ² become
   `None`). Do **not** silently invent σ.
2. **Per-point arrays** (over valid points, aligned to q):
   - `norm_residual` = (I − M) / σ
   - `frac_residual` = (I − M) / I
   - `frac_uncertainty` = σ / I
3. **Global scalars:**
   - `reduced_chi2` = Σ r_i² / (N_valid − n_params)   (guard dof > 0)
   - `robust_scale_s` = 1.4826 · MAD(norm_residual)
   - `sigma_misscale_factor` = `robust_scale_s` (alias; this is the factor by
     which actual scatter exceeds reported σ)
   - `realistic_reduced_chi2_floor` = `robust_scale_s ** 2`
   - `median_frac_uncertainty` = median(|σ/I|)
   - `max_abs_frac_misfit` and the q at which it occurs (`q_at_max_frac_misfit`)
   - `n_outliers_3s` = count(|norm_residual| > 3·robust_scale_s)
   - `frac_outliers_3s` = n_outliers_3s / N_valid
4. **Structure scalars:**
   - `longest_same_sign_run` of (I − M)
   - `sign_autocorr_lag1` (lag-1 autocorrelation of sign of residuals)
5. **Per-band** (split valid points into `n_bands` log-spaced Q windows; return a
   list of dicts): `q_lo, q_hi, n, reduced_chi2, robust_scale_s,
   max_abs_frac_misfit`.

Return everything in one flat-ish dict (scalars + arrays + `bands` list +
`sigma_available`). Keep arrays as numpy; callers serialize as needed.

> **Do NOT bake decision thresholds into this function.** No "good/bad" verdict,
> no "stop fitting" flag. Thresholds and the decision tree live in the consumer
> (the pyirena-ai agent strategy/skill). This function returns *facts only*.
> That separation is deliberate: different workflows interpret the same numbers
> differently.

### 3.2 Wire into existing fits (additive, backward-compatible)

For each `fit()` that already returns a results dict (start with
`core/unified.py`, then sizes / simple_fits / modeling), call
`fit_quality_metrics(...)` after the fit converges and add the result under a new
key, e.g. `results['quality'] = fit_quality_metrics(...)`. Do **not** remove or
change existing keys (`chi_squared`, `reduced_chi_squared`, `residuals`).

Note in `unified.py`: the internal `_residuals()` are weighted ((I−M)/σ) but the
returned `residuals` key is raw (I−M). The new function takes raw arrays (q, I,
M, σ) and derives everything itself, so there is no ambiguity to worry about.

## 4. User-facing exposure (requested)

The requester's intent: **expose these to end users too, and prefer the robust /
fractional residual over the plain normalized residual where a single residual
view is shown** — but this will be tested in pyirena first before the pyirena-ai
side relies on it.

Suggested, in priority order:

1. **`fit_quality_metrics` is public API** (export from `pyirena/__init__.py` or
   the documented core API) so users and the MCP server can call it on any fit.
2. **Plotting:** where a normalized-residual panel is currently shown
   (`plotting/unified_plots.py`, and the GUI residual displays), offer the
   **rescaled residual `r' = r/s`** and/or the **fractional misfit `(I−M)/I` in
   %** as the default or as a selectable view. Annotate the plot with
   `robust_scale_s` and `realistic_reduced_chi2_floor` so the user immediately
   sees "your σ look ~3× underestimated; χ²≈9 is the real floor."
3. **GUI / advisor summary line:** alongside reduced χ², show `robust_scale_s`
   and `max_abs_frac_misfit` (%, with its Q). These two numbers are the most
   decision-relevant.
4. **MCP server:** include the `quality` dict in fit-result responses.

Keep the raw normalized residual available too — this adds a view, it does not
remove one.

## 5. Testing

In `pyirena/tests/` add `test_fit_metrics.py` covering:

- **Well-calibrated σ:** synthetic data = model + Gaussian noise at exactly σ ⇒
  `robust_scale_s ≈ 1`, `reduced_chi2 ≈ 1`, no outliers, short sign runs.
- **Under-estimated σ (×3):** same data, σ divided by 3 ⇒ `robust_scale_s ≈ 3`,
  `realistic_reduced_chi2_floor ≈ 9`, still **no** 3s outliers (proves the
  scale problem is not mistaken for misfit).
- **Localized misfit:** inject a bump the model can't fit in one Q-band ⇒ bulk
  `robust_scale_s ≈ 1`, `n_outliers_3s > 0`, `max_abs_frac_misfit` large and its
  q inside the injected band; that band's per-band χ² much larger than others.
- **Systematic misfit (wrong slope):** ⇒ long `longest_same_sign_run`, high
  `sign_autocorr_lag1`, even if magnitudes are modest.
- **σ missing / all ≤ 0:** ⇒ `sigma_available == False`, `s`/χ² are `None`,
  fractional arrays still populated, no crash.
- **dof ≤ 0** guard.

## 6. Decision logic stays in pyirena-ai (FYI, not to implement here)

For context only — the agent side will consume the dict roughly like this, so
the field names above are chosen to make this easy:

- Bulk `|r'| < ~3`, no sign-runs, `max_abs_frac_misfit < ~15 %`
  → fit is as good as the data allows; **stop** (χ² > 1 is just mis-scaled σ).
- Any `n_outliers_3s > 0` with large `max_abs_frac_misfit` (≳ 30 %), or a hot
  per-band χ², or long sign-runs → **real misfit; investigate that Q-region.**
- `robust_scale_s ≪ 1` while still adjusting → **over-fitting; back off.**

That logic is intentionally **not** part of `fit_metrics.py`.
