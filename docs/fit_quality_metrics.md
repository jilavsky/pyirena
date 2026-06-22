# Understanding Fit Quality Metrics

This guide explains the robust fit-quality diagnostics available in pyirena's Unified Fit and API tools. These metrics help you judge whether a fit is truly good, or whether a high χ² is merely a consequence of mis-scaled reported uncertainties σ.

**TL;DR:** Don't chase reduced χ² ≈ 1. Instead, read `robust_scale_s` (how many times actual scatter exceeds reported σ) and `max_abs_frac_misfit` (the largest |(I−M)/I|, a σ-independent gross-misfit backstop). If `robust_scale_s` is large, the realistic χ² floor is `robust_scale_s²`, not 1.

---

## The Problem

Reported uncertainties σ on scattering data are frequently **mis-scaled** — they are often 2–5× too small or too large. This breaks the standard χ² ≈ 1 target:

- **σ too small (common):** reduced χ² shoots to 10–100, and you waste effort chasing a target the data can't support.
- **σ too large (rare):** reduced χ² stays near 1, but you may miss real misfits because the metric is too forgiving.

**The solution:** separate the *absolute scale* of σ (which fails) from the *relative shape* (which is usually trustworthy). Use a robust estimate of the true noise scale to re-calibrate, and add a σ-independent backstop to catch real misfits.

---

## The Metrics

### Core Metric: `robust_scale_s`

**What it is:** A MAD-based (median absolute deviation) estimate of how many times the *actual* scatter exceeds the *reported* σ.

**Formula:** `s = 1.4826 × median(|r_i − median(r_i)|)` where `r_i = (I_i − M_i)/σ_i` are the normalized residuals.

**How to read it:**

| `robust_scale_s` | Meaning | What to do |
|---|---|---|
| ≈ 1.0 | σ are honest; χ² ≈ 1 is a real target. | Standard interpretation: χ² too high → real misfit, fix the model. χ² too low → over-fitting, be careful. |
| ≈ 2–3 | σ are ~2–3× too small. | The *realistic* reduced χ² floor is `robust_scale_s²` ≈ 4–9. A fit at that floor with no outliers is **done — stop tightening.** |
| ≈ 5–10 | σ are ~5–10× too small (extreme mis-scaling, but seen in SAXS). | Realistic floor is 25–100. Don't chase χ² at all; use structure + fractional misfit instead. |
| ≪ 1 (e.g., 0.3) | σ are over-estimated; potential over-fitting risk. | Monitor other metrics. If `longest_same_sign_run` is short and scatter is random, you're OK. If it grows, back off on parameters. |

**Why MAD?** The median absolute deviation ignores outliers, so it reflects the *bulk* noise distribution, not isolated misfits. This is exactly what separates "σ is mis-scaled" (the whole distribution is wide) from "the model is wrong at one Q-region" (a few outliers stick out).

### The Realistic χ² Floor: `realistic_reduced_chi2_floor`

**What it is:** `robust_scale_s²` — the lowest reduced χ² the data can physically support given the actual scatter.

**How to read it:**
- If `robust_scale_s = 3`, then `realistic_reduced_chi2_floor ≈ 9`.
- A fit that reaches χ² = 9 with no outliers, no sign-structure, and no hot Q-bands is **as good as the data allows**. Do not keep tightening.
- If your fit is at χ² = 12 but the floor is 9, you're only 3 points above the limit — likely acceptable.

---

### Sigma-Independent Backstop: `max_abs_frac_misfit` (+ `q_at_max_frac_misfit`)

**What it is:** The largest |(I − M)/I| — the fractional error at the worst point, expressed as a percentage.

**Why it matters:**
- **σ-independent:** whether or not σ are reliable, a point where the model is 50 % off from the data is a real misfit.
- **Gross-misfit detector:** rules out the argument *"normalized residuals are 20–50, but σ are unreliable so that's fine"* — **wrong**. If σ ≈ 5 % and normalized residual = 20, then |(I−M)/I| ≈ 100 % — the model is completely wrong at that point.

**How to read it:**

| `max_abs_frac_misfit` (%) | Meaning |
|---|---|
| < 5 % | Excellent; the model stays within ~5 % of the data everywhere. |
| 5–15 % | Reasonable; one or two points deviate by 5–15 %. Typical for noisy SAXS. |
| 15–30 % | Noticeable misfit in a Q-region; may be acceptable if it's an edge (very low or high Q) or if structure diagnostics are clean. |
| > 30 % | **Likely a real misfit.** The model is off by more than 30 % at some point — investigate that Q-region. May need more parameters or a different model. |

**Example:** You see `reduced_chi2 = 25` and `robust_scale_s = 5` (so floor ≈ 25 is OK), but `max_abs_frac_misfit = 45 %`. The high χ² is not just mis-scaled σ — the model is genuinely 45 % off somewhere. Fix that region.

---

### Structure Diagnostics: `longest_same_sign_run` and `sign_autocorr_lag1`

**What they are:**
- **`longest_same_sign_run`:** the longest consecutive stretch of residuals with the same sign (both positive or both negative).
- **`sign_autocorr_lag1`:** lag-1 autocorrelation of the residual signs; near +1 means stretches of same-sign; near 0 means random scatter.

**Why they matter:** distinguish a *wrong functional form* (systematic) from a *mis-scaled σ* (random).

**How to read them:**

| Pattern | `longest_same_sign_run` | `sign_autocorr_lag1` | Diagnosis |
|---|---|---|---|
| Random scatter (σ issue) | short (< N/4) | near 0 | σ mis-scaled; fit is OK structurally. |
| Systematic misfit (model wrong) | long (> N/2) | > 0.5 | Wrong slope, missing level, or wrong functional form. Fix the model. |
| Clean fit | very short | near 0 | Excellent; bulk distribution is random. |

**Example 1:** `longest_same_sign_run = 150` (out of 250 points), `sign_autocorr_lag1 = 0.8` → the model systematically over- or under-predicts over long Q-stretches. The slope may be wrong, or a level is missing. Not a σ-scale issue.

**Example 2:** `longest_same_sign_run = 40`, `sign_autocorr_lag1 = 0.1` → random scatter; no systematic pattern. If χ² is high, it's because σ are small, not because the model is bad.

---

### Per-Band Breakdown: `bands` (Q-decade windows)

**What it is:** the same metrics (reduced χ², robust scale, max fractional misfit) computed in 3–4 log-spaced Q windows.

**Why it matters:** uneven per-band χ² is a red flag for localized misfit.

**How to read it:**

```
bands:
  - {q_lo: 0.001, q_hi: 0.01, n: 60, reduced_chi2: 2.1, robust_scale_s: 1.2, max_abs_frac_misfit: 0.08}
  - {q_lo: 0.01,  q_hi: 0.1,  n: 95, reduced_chi2: 12.5, robust_scale_s: 3.5, max_abs_frac_misfit: 0.35}
  - {q_lo: 0.1,   q_hi: 1.0,  n: 95, reduced_chi2: 1.8, robust_scale_s: 1.1, max_abs_frac_misfit: 0.05}
```

**Interpretation:** Band 2 (Q = 0.01–0.1) is hot — χ² = 12.5, `robust_scale_s = 3.5` (σ ~3.5× underestimated there), and max misfit = 35 %. Either:
1. σ are especially mis-scaled in that Q-range (instrument issue, beamstop edge, etc.), or
2. The model fails in that Q-range (missing feature, wrong level parameters, etc.).

Check the Q-region visually (plot residuals), and adjust model or re-examine error bars.

---

## Interpretation Decision Tree

Use this flowchart to decide what to do next:

```
1. Compute fit_quality_metrics() after fitting.

2. Check: robust_scale_s large?
   ├─ YES (s > 3)
   │  ├─ Check: max_abs_frac_misfit > 30 %?
   │  │  ├─ YES  → real misfit + mis-scaled σ. Investigate that Q-region & fix model.
   │  │  └─ NO   → pure σ-scale issue. χ² is at realistic floor; fit is done.
   │  └─ Check: longest_same_sign_run long?
   │     ├─ YES  → systematic misfit even after accounting for σ scale. Fix model.
   │     └─ NO   → random scatter; fit is OK.
   │
   └─ NO (s ≈ 1)
      ├─ Check: reduced_chi2 > 5?
      │  ├─ YES  → σ honest but fit is poor. Real misfit.
      │  └─ NO   → good fit.
      └─ Check: longest_same_sign_run long?
         ├─ YES  → systematic misfit.
         └─ NO   → random scatter; fit is OK.

3. Check: per-band χ² uneven?
   ├─ YES (one band hot, others cool) → localized misfit in that Q-band.
   └─ NO  → problem is global or absent.
```

---

## Common Scenarios

### Scenario 1: High χ² but looks fine visually

```
reduced_chi2 = 15
robust_scale_s = 4.0
max_abs_frac_misfit = 8 %
longest_same_sign_run = 35 (out of 250)
```

**Diagnosis:** σ are ~4× too small; the realistic floor is ~16. Your fit is at that floor with no outliers and random scatter. **The fit is done.** Do not keep tightening parameters.

---

### Scenario 2: High χ² and visually wrong

```
reduced_chi2 = 20
robust_scale_s = 1.5
max_abs_frac_misfit = 55 %
longest_same_sign_run = 180
```

**Diagnosis:** σ are slightly mis-scaled (s = 1.5), but that's not the main issue. Max fractional misfit is 55 % (model is off by more than half) and there's a long systematic run of same-sign residuals. **Real misfit.** The model is wrong — fix it (add a level, change P, etc.).

---

### Scenario 3: χ² low but suspicious

```
reduced_chi2 = 0.3
robust_scale_s = 0.2
```

**Diagnosis:** σ are ~5× over-estimated; you may be over-fitting. Check `longest_same_sign_run` and whether the number of free parameters is reasonable. If parameters are being driven to implausible values, back off.

---

### Scenario 4: One Q-band is hot

```
overall reduced_chi2 = 3.5
bands:
  - q: 0.001–0.01, chi2 = 1.2
  - q: 0.01–0.1,  chi2 = 25.0   ← hot
  - q: 0.1–1.0,   chi2 = 1.5
max_abs_frac_misfit = 0.38 (in the hot band)
```

**Diagnosis:** the model works well at low and high Q, but fails at intermediate Q (0.01–0.1). This Q-region may have:
- A missing structural feature (e.g., an aggregated level).
- Mis-estimated uncertainty (beamstop edge, slits, etc.).
- True data anomaly (crystal, contamination).

Focus your effort on that Q-band.

---

## Where to Find These Metrics

### In the Unified Fit GUI
- **Success message popup** after a fit displays: `robust_scale_s`, `realistic_reduced_chi2_floor`, and `max_abs_frac_misfit`.
- **Residual plot** now shows rescaled residuals `r' = r/s` (more informative when σ are mis-scaled).

### In the API / MCP
Call `pyirena.api.control.get_fit_quality(session_id)` to get the full dict (all scalars + arrays + per-band).

Call `pyirena.api.control.run_fit(session_id)` and read the `quality` block in the result (decision-relevant scalars only).

Call `pyirena.api.control.get_residuals(session_id)` to inspect `summary.robust_scale_s` and the rescaled / fractional residual arrays.

### Programmatically
```python
from pyirena.core.fit_metrics import fit_quality_metrics
import numpy as np

# After fitting your model:
metrics = fit_quality_metrics(
    q=q_data,
    intensity=I_data,
    model=I_model,  # fitted intensity
    sigma=error_data,  # may be None
    n_params=number_of_free_parameters,
)

print(f"robust_scale_s: {metrics['robust_scale_s']:.2f}")
print(f"max_abs_frac_misfit: {metrics['max_abs_frac_misfit']*100:.1f}%")
print(f"longest_same_sign_run: {metrics['longest_same_sign_run']}")
```

---

## References

- **Full specification:** `planning/ai-agent/04-fit-quality-metrics.md`
- **Implementation details:** `planning/ai-agent/05-fit-quality-metrics-implementation-plan.md`
- **AI agent tools:** `docs/ai_tools_reference.md` (section "Quality assessment")
- **Source code:** `pyirena/core/fit_metrics.py`

---

## FAQ

**Q: Can I still use the old normalized residuals?**

A: Yes. The rescaled residuals are the default in the GUI and matplotlib plots, but you can switch back by setting `_USE_RESCALED_RESIDUALS = False` in the code. The raw residuals are always available via the API.

**Q: What if σ are not provided?**

A: The σ-dependent metrics (`robust_scale_s`, `reduced_chi2` for dof) become `None`. The σ-independent metrics (fractional misfit, longest same-sign run) are still computed and useful.

**Q: Which metric is most important?**

A: In order: (1) `robust_scale_s` — tells you the σ situation immediately. (2) `max_abs_frac_misfit` — σ-independent gross-misfit detector. (3) `longest_same_sign_run` — structure diagnostic. (4) Per-band χ² — localization. Use all four together.

**Q: Do I need to understand the MAD formula?**

A: No. Just remember: `robust_scale_s` = how many times actual scatter exceeds reported σ. The MAD is robust to outliers, so it's a good measure of the *bulk* noise, not driven by one bad point.

**Q: What if I disagree with the diagnostics?**

A: Check your data (bad points, instrument glitches), your σ (are they really 3× too small?), and your model (is the functional form right?). The metrics reflect what's in the data; they don't lie, but they need good input to give good answers.
