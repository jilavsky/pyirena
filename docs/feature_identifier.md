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

## How the algorithm works

The detector takes four conceptual steps:

1. **Smooth log(I) in log(Q)** with a Gaussian kernel of width
   `sigma_smooth = 0.15 decades` (not user-adjustable).  This suppresses
   point-to-point noise without hiding real slope transitions, which in SAS
   data always span at least 0.3 decades.

2. **Compute the local slope** d(log I)/d(log Q) at every Q value using a
   central difference over ±0.15 decades around each point.

3. **Detect change-points** in the slope profile.  For each candidate
   boundary point, the algorithm compares the *mean* slope in a window to
   its left against the *mean* slope in a window to its right.  Local
   maxima of |mean_left − mean_right| that exceed a threshold are
   change-points; segments are the intervals between them.  This runs in
   **two passes**:
   - **Pass 1 (coarse)**: wider window, larger threshold — finds the major
     transitions cleanly.
   - **Pass 2 (refinement)**: any segment wider than
     `wide_region_decades` (default 1.0) is re-scanned with tighter
     parameters to find hidden sub-structure.  This catches the case where
     a single broad region actually contains multiple distinct power-law
     slopes that drift smoothly between each other — the Unified Fit model
     says smooth drift implies a blend of constant-slope levels.

4. **Merge and filter** the resulting segments: adjacent segments whose
   mean slopes are within `merge_slope_tol` are merged (suppressing
   spurious splits at noise spikes); segments narrower than
   `min_segment_decades` are dropped — except for segments at the very
   lowest or highest Q values, which use the looser
   `edge_min_segment_decades` so that narrow low-Q Guinier plateaus and
   high-Q backgrounds are not lost.

Then each segment is classified by its mean slope as `background`,
`guinier_plateau`, or `power_law` (see below), and adjacent segments are
checked to see if a Guinier knee can be inferred between them.

The change-point approach replaces an earlier variance-based detector that
could not distinguish "slope is constant" from "slope is slowly drifting" —
which led to missed knees at data edges (sample15) and lumping smoothly-
drifting multi-level data into one segment (sample25).

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

Click **Show advanced params** to expose the segmentation controls.  Values
are saved automatically between sessions (under the `feature_detect` key in
`~/.pyirena/state.json`); you do not need to re-enter them when you load a
new file into Unified Fit.

### Change-point detection — Pass 1 (coarse)

| Parameter | Default | Description |
|---|---|---|
| **Pass-1 window (dec)** | 0.30 | Half-width of the left/right windows used to compute the slope-mean difference.  A larger window averages over more points so the detector is less twitchy at noise spikes but cannot resolve transitions narrower than the window. |
| **Pass-1 threshold** | 0.40 | A change-point is declared at a local maximum of |mean_L − mean_R| only if the difference exceeds this value.  Larger value → fewer change-points (only the biggest transitions); smaller value → more change-points (catches subtler boundaries). |

Pass 1 is meant to find the *major* boundaries reliably.  If you want only
the most obvious transitions, raise the threshold.

### Change-point detection — Pass 2 (wide-region refinement)

| Parameter | Default | Description |
|---|---|---|
| **Pass-2 window (dec)** | 0.20 | Tighter window for the second pass, which re-scans regions wider than `wide_region_decades` to find hidden sub-structure. |
| **Pass-2 threshold** | 0.20 | Lower than Pass-1 threshold, so subtler slope changes inside a wide region are caught.  Wider regions are more likely to contain multiple slope drifts that should be separated. |
| **Wide-region threshold (dec)** | 1.0 | A segment found in Pass 1 wider than this gets a Pass-2 re-scan.  Smaller value → Pass 2 runs on more segments (catches more sub-structure but risks over-segmentation). |

The two-pass approach handles the common SAS case where two adjacent
structural levels have similar Porod slopes that connect smoothly — the
first pass treats the whole region as one segment, the second pass detects
the subtle mean-slope difference inside it.

### Segment filtering

| Parameter | Default | Description |
|---|---|---|
| **Min segment width (dec)** | 0.10 | Interior segments narrower than this are discarded.  Prevents noise spikes or one-or-two-point artefacts from appearing as structural features. |
| **Edge min width (dec)** | 0.05 | Looser width threshold for segments that touch the lowest-Q or highest-Q data point.  Recognises that real low-Q Guinier plateaus and high-Q backgrounds may genuinely be narrow when the data extent does not extend further in those directions. |
| **Merge slope tolerance** | 0.15 | Two adjacent segments are merged into one if their mean slopes differ by less than this value.  Suppresses spurious splits when a noise spike triggers a change-point inside an otherwise-uniform region.  Set to 0 to disable merging. |

### Guinier knee sensitivity

| Parameter | Default | Description |
|---|---|---|
| **Min knee Δslope** | 0.10 | A Guinier knee is reported between two adjacent segments only if their mean slopes differ by at least this amount.  Larger value → only report obvious knees; smaller value → report gentle transitions between segments of similar slope. |

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

## Tuning guidance

### "I am getting too many segments"

The default thresholds are set to be sensitive — they detect every slope
transition the data supports.  To get fewer, larger segments:

- Raise **Pass-1 threshold** (e.g. 0.60 instead of 0.40) — coarser
  segmentation, only big transitions detected
- Raise **Merge slope tolerance** (e.g. 0.30 instead of 0.15) — adjacent
  similar segments get merged
- Raise **Wide-region threshold** (e.g. 2.0 instead of 1.0) — Pass 2 runs
  less often

### "I know there is a knee here but the algorithm misses it"

If your suspected knee is at the very low-Q or high-Q end of the data:

- The default detection requires the boundary to be at least 0.05 decades
  inside the data extreme.  Try lowering **Edge min width**.

If the knee is in the interior and the slope transition is gradual:

- The transition may not produce a sharp enough mean-slope difference to
  cross **Pass-1 threshold**.  Lower it (e.g. 0.25 instead of 0.40).
- Or lower **Pass-2 threshold** (e.g. 0.10) so subtler transitions inside
  wide regions are caught.

### "Two segments have nearly identical slopes — why are they separate?"

Lower **Merge slope tolerance** is supposed to keep them separate; raise it
above their slope difference to merge them.

### "I want to see the bare-bones segmentation without the refinement pass"

Set **Wide-region threshold** to a very large value (e.g. 10.0) so no
segment is ever wide enough to trigger Pass 2.

---

## Practical workflow

1. Load a dataset into Unified Fit.
2. Click **Identify Features…** → **Detect segments**.
3. Look at the summary panel (right side): it lists segments from high-Q to
   low-Q with their Q ranges and slopes.  The entry at the top ("L1") should
   correspond to Level 1 in Unified Fit (the smallest structures, highest Q).
4. **Choose number of levels**: each orange (PLS) segment that has a Guinier
   knee to its low-Q side is a candidate level.  The "Suggested Unified Fit
   levels" count is a starting point; you may need fewer if some PLS regions
   are noise.
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
  less than ≈ 0.5-1 decade apart may be detected as a single broad segment
  unless Pass 2 catches the slight slope drift.

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

### Version history

- **v0.8.5** (current): change-point segmentation with two-pass refinement
  and edge-aware width filter.  Config schema changed; old saved state is
  silently ignored via a schema-version check.
- **v0.8.4**: variance-based stability detection.  Could not distinguish
  smooth slope drift from constant slope; missed edge knees because the
  stability window straddled the transition zone.
- **v0.8.3**: threshold-based plateau detector (initial release).
