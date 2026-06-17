# Feature Identifier — User Guide

The Feature Identifier analyses your I(Q) data in log-log space and segments
the curve into regions where the power-law slope is approximately constant.
This tells you where to look for structural levels in a Unified Fit and what
Q ranges to hand to the local Guinier / Porod fitting buttons.

---

## How to open it

In the Unified Fit panel, load any dataset, then click the **Identify
Features…** button (blue, in the top control row).  A separate non-modal
window opens.  Click **Detect segments** to run the analysis on the
currently loaded data.  The dialog stays open while you work; close it to
clear all markers from the graph.

**Visualisation only** — the detector never changes any level parameters or
the fit model.  It is safe to run at any point without undoing work.

---

## What is detected

### Segments

A segment is a Q range over which the local log-log slope
`d(log I) / d(log Q)` is approximately constant.  Three kinds:

| Kind | Colour | What it means |
|---|---|---|
| **Power-law (PLS)** | Orange | Slope is steep and stable (Porod, fractal, or any constant power-law exponent) |
| **Guinier plateau (GP)** | Green | Slope is near zero — intensity is approximately flat in log-log |
| **Background** | Grey | Slope is near zero **and** the segment touches the high-Q end of the data |

Segments are listed in the summary **from high-Q to low-Q** — this is the
order used in Unified Fit, where Level 1 corresponds to the smallest
structures at the high-Q end.

### Guinier knees

A Guinier knee is the transition between two adjacent segments where the
slope changes from steep (high-Q Porod tail) to shallow (low-Q Guinier
region).  Physically this means:

- At high Q: intensity falls steeply, dominated by the Porod power-law
- At low Q: intensity flattens toward the Guinier plateau of that level

A knee is only reported when **the low-Q segment has a shallower slope than
the high-Q segment** (|slope_low_Q| < |slope_high_Q|).  If the slope at
low Q is *steeper* than at high Q, there is no Guinier knee — that kind of
transition is more likely a mass-fractal crossover, aggregate structure, or
instrument artefact; it does not imply a level boundary.

### Recommended Guinier windows

For each detected Guinier knee, the dialog suggests a Q window suitable for
the "Fit Rg/G btwn cursors" button:

- **q_min_guinier**: left edge of the Guinier plateau (low-Q side of the knee)
- **q_max_guinier**: beginning of the Porod tail × 1.2 (a small margin into the tail)
- **q_min_powerlaw**: where the Porod tail begins (use this as the left cursor
  for "Fit P/B btwn cursors")

---

## Advanced parameters

Click **Show advanced params** to expose the segmentation controls.  The
values are saved automatically between sessions and across dataset changes —
you do not need to re-enter them when you load a new file into Unified Fit.

### Smoothing

| Parameter | Default | Description |
|---|---|---|
| **Stability window (dec)** | 0.40 | Width in log(Q) decades of the sliding window used to evaluate whether the slope is locally stable.  Larger values capture wider slope plateaus but may merge two distinct levels. |

The slope is computed after Gaussian smoothing of log(I) with σ = 0.15 decades
(not user-adjustable).  This suppresses point-to-point noise without hiding
real slope transitions, which always span ≥ 0.3 decades in SAS data.

### Stability threshold

| Parameter | Default | Description |
|---|---|---|
| **Stability std max** | 0.40 | A point is marked stable if the standard deviation of slope within the ± stability-window is below this value.  The default (0.40) is calibrated against 31 hand-labelled SAXS/USAXS curves.  Decrease to 0.20-0.25 if you see too many spurious segments; increase toward 0.60 if real Porod regions are missed. |

**Why 0.40 and not something smaller?** Real SAS data has smooth Guinier-knee
transitions that look like gradually changing slope over ≈ 0.5-1 decade.
Requiring slope stability too tightly (std < 0.10-0.15) only detects the
flat top of the Porod tail and misses the Guinier region on the other side
of the knee.

### Minimum segment width

| Parameter | Default | Description |
|---|---|---|
| **Min segment width (dec)** | 0.20 | Stable regions narrower than this are discarded.  Prevents noise spikes or very narrow diffraction artefacts from appearing as structural features.  Should be at least 2 × σ_smooth ≈ 0.30 decades in practice; the default 0.20 is intentionally a little looser to catch short segments near data edges. |

### Merging adjacent segments

| Parameter | Default | Description |
|---|---|---|
| **Merge slope tolerance** | 0.40 | Two adjacent stable runs are merged into one segment if their mean slopes differ by less than this value.  This handles the common case where a single power-law region is interrupted by a brief unstable point (noise spike or a few bad data points).  Set to 0 to disable merging. |

The gap between two runs must also be narrower than **0.20 decades** (not
user-adjustable) for merging to be considered.

### Guinier knee sensitivity

| Parameter | Default | Description |
|---|---|---|
| **Min knee Δslope** | 0.50 | A knee is reported between two adjacent segments only if their mean slopes differ by at least this amount.  Larger value → only report obvious knees; smaller value → report gentle transitions between regions of similar slope. |

Note: even if Δslope ≥ this threshold, a knee is **not** reported if the
low-Q side is steeper than the high-Q side (see the physical constraint
described above).

### Q clipping

| Parameter | Default | Description |
|---|---|---|
| **Q max clip (Å⁻¹, 0=off)** | 0.60 | Data above this Q value is excluded before segmentation.  The small-angle approximation breaks down above Q ≈ 0.6 Å⁻¹ for most systems; beyond that the curve typically shows amorphous diffraction peaks, Bessel function oscillations, or other features that are not SAS power-law structure and should not be modelled with Unified Fit.  Set to 0 to analyse the full Q range. |

---

## Interpreting the markers

After clicking **Detect segments**, the main I(Q) graph shows:

| Marker | Meaning |
|---|---|
| Orange region | Power-law slope (PLS) — stable Porod/fractal/diffraction region |
| Green region | Guinier plateau (GP) — roughly flat region at low Q |
| Grey region | Background — flat region at the high-Q data end |
| Red dashed band | Guinier knee — transition zone between two slope segments |

The label on each power-law region shows `PLS P=X.X` where X is the mean
absolute slope (the Porod exponent P).  The label on a knee shows
`GK Δ=X.X` where X is the magnitude of the slope change.

**What is NOT shown:** the transition zones themselves (where slope changes
from one stable value to another) are between the coloured regions.  These
are the actual Guinier knees; they are marked separately with the red band.

---

## Practical workflow

1. Load a dataset into Unified Fit.
2. Click **Identify Features…** → **Detect segments**.
3. Look at the summary panel (right side): it lists segments from high-Q to
   low-Q with their Q ranges and slopes.  The entry at the top ("L1") should
   correspond to Level 1 in Unified Fit (the smallest structures, highest Q).
4. **Choose number of levels**: each orange (PLS) segment that has a Guinier
   knee to its low-Q side is a candidate level.  The "Suggested Unified Fit
   levels" count is a starting point; you may need more if the algorithm missed
   a subtle knee.
5. **Place cursors**: drag the Unified Fit cursors to bracket a suggested
   Guinier window, then click "Fit Rg/G btwn cursors" for each level in turn.
6. Click **Clear markers** when done, then close the dialog.

---

## Limitations and known difficult cases

- **Bessel function oscillations** (dense monodisperse spheres): the oscillation
  envelope looks like a series of alternating steep/flat regions and will be
  over-segmented.  These cannot be modelled by Unified Fit; set the Unified Fit
  Q range to avoid this region.
- **Structure-factor peaks** (concentrated suspensions): the intensity dip and
  peak near the structure-factor maximum creates a false slope transition.
  The detector cannot distinguish this from a true Guinier knee.
- **Very noisy data** at low Q (USAXS, extremely dilute samples): segments may
  be missed if most points in a region have σ_I/I > 0.5 (the noise-filter
  threshold).
- **Overlapping levels** with similar Rg: two levels whose Guinier knees are
  less than ≈ 0.5-1 decade apart may be detected as a single broad segment.

---

## Technical notes for developers

The detector lives in `pyirena/core/feature_detect.py`.  The main entry point
is `detect_features(q, I, sigma_I=None, config=None) → FeatureDetectResult`.

The result schema:
```python
FeatureDetectResult.segments        # list of segment dicts (low-Q → high-Q)
FeatureDetectResult.guinier_knees   # list of knee dicts
FeatureDetectResult.background_q_min
FeatureDetectResult.recommended_guinier_windows
FeatureDetectResult.recommended_nlevels
FeatureDetectResult.log_decades
FeatureDetectResult.n_points
FeatureDetectResult.q_min_analysed
FeatureDetectResult.q_max_analysed
FeatureDetectResult.n_segments_found
```

Each knee dict has keys: `q_min`, `q_max`, `q_center`, `slope_low_q`,
`slope_high_q`, `delta_slope`.  A knee is only emitted when
`|slope_low_q| < |slope_high_q|` (slope shallower at low Q).

The algorithm is also available to the AI agent via the MCP tool
`pyirena_ctrl_detect_features(session_id, q_max_clip=0.6, config_overrides={})`.
See `docs/ai_tools_reference.md` for the full tool list.
