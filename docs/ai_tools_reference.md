# pyirena MCP Tools — Reference for AI Agents

This document is written for AI assistants and the humans who craft
system prompts for them. It explains every tool the `pyirena-mcp` server
exposes, when to use which, and recommended call patterns.

If you are an AI assistant with access to pyirena tools, read this
document — or have it inlined into your system prompt — before asking
questions about a user's analysis data.

> **Tool name prefixes:**
> - `pyirena_` — read-only tools that query existing fit results
> - `pyirena_ctrl_` — control tools that drive fitting interactively
>
> Both sets are globally unambiguous when the client connects to multiple
> MCP servers. The underlying Python library uses `pyirena.api` (read)
> and `pyirena.api.control` (control) with no prefix.

---

## What pyirena gives you

pyirena works in two modes. **Read mode** retrieves existing fit results
stored in NXcanSAS HDF5 files. **Control mode** lets you actually run fits —
configure models, set parameters, execute the fitting algorithm, evaluate
results, and save back to HDF5.

### Read-mode analysis tools (results stored in NXcanSAS files)

- **Reduced data** — the raw I(Q) curve, plus measurement uncertainties.
- **Unified Fit** — Beaucage multi-level model with parameters G, Rg, B,
  P, ETA, PACK, RgCutoff (up to 5 levels per file).
- **Simple Fits** — single-model fits (Guinier, Porod, Debye-Bueche…).
- **Size Distribution** — particle size distribution P(r), with bulk
  volume fraction.
- **Modeling** — parametric forward model with multiple populations:
  size-distribution, unified-level, diffraction-peak, Guinier-Porod,
  mass-fractal, surface-fractal.
- **SAXS Morph** — voxelgram-based forward modeling.
- **WAXS Peak Fit** — per-peak position, FWHM, amplitude, area.
- **Fractals** — grown mass-fractal aggregates (Z, df, dmin, c, Rg).
- **Data Merge / Data Manipulation** — provenance for combined or
  transformed datasets.

Your job is to help the user understand what's in their data, and — when
asked — to fit new datasets autonomously using the control tools.

### Control-mode fitting (Phase 1: Unified Fit)

The `pyirena_ctrl_*` tools let you drive pyirena's Unified Fit model
end-to-end. Sessions are in-memory for the lifetime of the MCP server
process and persist across tool calls within a conversation.

---

## Recommended discovery workflow

Almost every interaction should start with cheap orientation calls:

1. **`pyirena_summarize_folder(folder)`** — get a one-shot snapshot:
   how many files, which samples, which analyses are present per file.
   Use this first whenever the user mentions a folder.
2. **`pyirena_list_files(folder, sort="mtime_desc", limit=20)`** — when
   you need to identify specific files. The newest are usually most
   relevant.
3. **`pyirena_inspect_file(path)`** — drill into one file: sample name,
   exact list of analyses present, Q range, number of points.
4. **`pyirena_read_<tool>(path)`** — only now load the actual results.

Do not call `pyirena_read_*` functions until you know which analyses
exist — otherwise you'll get `{"found": false}` and waste a turn.

---

## Common conventions

- **Every call returns a dict.** Never a raw array, never an object.
  All values are JSON-safe (NaN and inf are converted to `null`).
- **Missing data returns `{"found": false}`** — never an exception. You
  can defensively check `result["found"]` before using fields.
- **Arrays are decimated by default** to about 500 points so responses
  stay small. To get the full curve, pass `include_arrays=True` on
  per-tool readers or `include_full=True` on `pyirena_read_reduced_data`.
  Only do this when you specifically need the high-resolution shape.
- **File paths must be absolute** (or relative to `PYIRENA_DATA_ROOT` if
  that env var is set by the operator).
- **`PYIRENA_DATA_ROOT` may restrict access.** If a file path is rejected
  with a security error, tell the user to widen `PYIRENA_DATA_ROOT` or
  put the file inside the configured root.

---

## Tool reference (grouped)

### Discovery — call these first

#### `pyirena_summarize_folder(folder, sample_filter=None)`
Aggregate snapshot. Returns total file count, distinct sample names, per-
analysis file counts (e.g. `{"unified_fit": 30, "modeling": 12}`), and
mtime range. Optional `sample_filter` is a case-insensitive substring.

When to use: first call when the user mentions a folder or "all my data".

#### `pyirena_list_files(folder, pattern="*.h5,...", sort="mtime_desc", limit=100, deep=True)`
One row per file. Sort keys: `name_asc/desc`, `mtime_asc/desc`,
`size_asc/desc`. Use `deep=False` for a faster shallow listing if you
only need paths/sizes/mtime.

When to use: you need to identify specific files (latest 3, oldest,
matching a name pattern).

#### `pyirena_inspect_file(path)`
Single-file deep look. Returns sample, scan number, mtime, list of
`analyses_present`, Q range, number of points.

When to use: before any `pyirena_read_*` call, to confirm what's in the
file.

---

### Reading reduced data

#### `pyirena_read_reduced_data(path, decimate=500, include_full=False)`
The raw measured I(Q). Default decimation keeps response small. Set
`include_full=True` only if you need the actual point density (rare for
question-answering).

#### `pyirena_read_metadata(path)`
Sample name, label, thickness, blank file, instrument, timestamp.

---

### Reading per-tool results

All take `(path, include_arrays=False, max_points=500)` unless noted.
Arrays are **omitted by default** to keep responses compact — pass
`include_arrays=True` only when you specifically need the curves.

#### `pyirena_read_simple_fit(path)`
Returns `{model, success, chi_squared, reduced_chi_squared, dof, q_min,
q_max, params{}, params_std{}, derived{}, ...}`. Single model per file.

#### `pyirena_read_unified_fit(path)`
Returns `{num_levels, background, chi_squared, levels: [...]}` where each
level has `{G, Rg, B, P, RgCutoff, ETA, PACK, correlations,
G_err, Rg_err, ...}`. Use this for the most common question:
hierarchical structure with multiple Rg's.

#### `pyirena_read_size_distribution(path)`
Returns `{method, shape, volume_fraction, rg, chi_squared, n_iterations,
q_power, r_min, r_max, n_bins, ...}`. With `include_arrays=True`:
`r_grid`, `distribution` (volume-weighted P(r)), `number_dist`,
`cumul_vol_dist`.

#### `pyirena_read_modeling(path)`
Returns `{chi_squared, background, populations: [...]}` where each
population has `{pop_type, label, parameters{}, derived{}}`. `pop_type`
is one of: `size_dist`, `unified_level`, `diffraction_peak`,
`guinier_porod`, `mass_fractal`, `surface_fractal`. Population
parameters vary by type — refer to the `parameters` dict for what's
present.

#### `pyirena_read_saxs_morph(path)`
Returns `{chi_squared, volume_fraction, contrast, rg_A, phi_actual,
voxel_size, box_size_A, morphology_metrics: {...}, ...}`. The 3-D
voxelgram itself is intentionally not exposed (too large). Morphology
metrics include pore-size statistics, Euler number, open/closed
porosity.

#### `pyirena_read_waxs_peakfit(path)`
Returns `{n_peaks, bg_shape, chi_squared, bg_params{}, peaks: [...]}`
where each peak has `{shape, params{Q0, A, FWHM, eta, ...},
params_std{}, area, area_std}`.

#### `pyirena_read_fractals(path)`
Returns metadata for each grown fractal aggregate in the file:
`{aggregates: [{Z, df, dmin, c, RgPrimary, RgAggregate, ...}]}`. Full
intensity arrays per aggregate are not exposed via this tool — they
require direct HDF5 access.

#### `pyirena_read_merge_provenance(path)`
Returns `{scale, q_shift, background, chi_squared, q_overlap_min,
q_overlap_max, ds1_file, ds2_file}`. Tells you how a merged file was
constructed.

#### `pyirena_read_manipulation_provenance(path)`
Returns `{operation, source_file, parameters{...}}`. Tells you whether a
file was scaled, trimmed, rebinned, subtracted, divided, or averaged.

---

### Cross-file aggregation

#### `pyirena_tabulate_parameter(folder, tool, parameter, x_axis="scan_number", subgroup_index=None, sample_filter=None)`
The workhorse for trend questions. Extracts ONE scalar across every file
in `folder` that contains the requested tool's results.

- `tool` is one of: `unified_fit, simple_fits, size_distribution,
  modeling, saxs_morph, waxs_peakfit, fractals, data_merge,
  data_manipulation`.
- `parameter` is a scalar key from that tool's schema. Common ones:
  - unified_fit: `background`, `chi_squared`, `Rg` (per-level →
    `subgroup_index=1` for level 1)
  - simple_fits: `chi_squared`, `param_Rg`, `param_I0`, `param_B`
  - size_distribution: `volume_fraction`, `rg`, `chi_squared`, `q_power`
  - modeling: `chi_squared`, `background`
  - waxs_peakfit: `peak_Q0`, `peak_FWHM`, `peak_A`, `peak_Area` (each
    requires `subgroup_index` = peak number, 1-based)
- `x_axis` controls row ordering: `scan_number` (from filename),
  `mtime`, or `name`.

Returns `{n_rows, rows: [{path, name, sample, scan_number, mtime, value,
stddev}], units, label}`. Use `value` for the plot and `stddev` for
error bars when available.

#### `pyirena_summarize_sample(folder, sample)`
Condenses one sample's history across the folder: how many files,
analyses run, and min/max/n of every top-level scalar parameter. Use
when the user asks "tell me everything about sample_X".

---

### Plotting (returns mixed text + image content)

Both plotting tools return a **two-item content list**, not a single
value. You must handle both items by content type:

| Index | `type` | What it contains |
|-------|--------|-----------------|
| 0 | `text` | `"Plot saved to: /abs/path/to/file.png"` — the on-disk location |
| 1 | `image` | The PNG encoded as base64; `mimeType` is `"image/png"` |

**Critical: never print or forward item 1 as a string.** It is raw
base64 binary data. If your client or agent loop iterates over content
items without checking `type`, it will dump thousands of characters of
garbled text to the user. Always branch on `content[i].type`:

```python
for item in tool_result.content:
    if item.type == "text":
        print(item.text)          # file path — safe to show
    elif item.type == "image":
        display_image(item.data)  # base64 PNG — render, don't print
```

If the AI client renders content items by type natively (Claude Desktop,
Claude Code, MCP Inspector), the image appears inline automatically and
you do not need to handle it manually.

If the client does not render images (AnythingLLM in some modes, custom
pipelines), skip the image item entirely and tell the user the file path
from item 0 so they can open the PNG manually.

#### `pyirena_plot_iq(paths, overlay=True, log_x=True, log_y=True, output_path=None)`
Plots I(Q) for one or more files. Saves PNG to `PYIRENA_PLOT_CACHE`
(or `output_path`) and returns both items described above.

- `overlay=True` (default): one axes, all curves overlaid (colour-coded
  by file).
- `overlay=False`: a grid, one subplot per file. Use when there are
  many curves and overlay would be unreadable.
- WAXS data: pass `log_x=False, log_y=False` for the standard
  linear-linear WAXS view.

#### `pyirena_plot_parameter_trend(folder, tool, parameter, x_axis="scan_number", subgroup_index=None, sample_filter=None, output_path=None)`
Internally calls `pyirena_tabulate_parameter` then renders the result as
a line+marker plot with error bars (when `stddev` is available). Saves
PNG and returns both items described above.

Use for the very common request: "plot how Rg evolved across all my
scans for sample X."

---

## Control tools reference (`pyirena_ctrl_` prefix)

These tools use sessions. Always start with `pyirena_ctrl_open_dataset()`,
capture the returned `session_id`, and pass it to every subsequent call.

### Recommended fitting workflow

Feature detection is **step 0 before selecting a model** — it tells you how
many levels the curve needs and where to place Q-windows, so you don't
start fitting blind.

```
pyirena_ctrl_open_dataset(file_path)           → session_id + data summary

# Step 0: understand the curve structure FIRST
pyirena_ctrl_detect_features(session_id)       → segments, knees, recommended_nlevels
# → use recommended_nlevels for nlevels below
# → use guinier_knees[i].q_min..q_max for fit_local_guinier Q windows
# → use segments for initial P estimate per level

pyirena_ctrl_select_model(session_id, nlevels=N)   # N from detect_features
pyirena_ctrl_get_model_description(session_id) → read before fitting
pyirena_ctrl_get_data_q_range(session_id)      → know your Q range
pyirena_ctrl_set_fit_q_range(session_id, q_min=…, q_max=…)  # if needed

# Optional: get local starting values from the detected regions
pyirena_ctrl_fit_local_guinier(session_id, q_min=knee.q_min, q_max=knee.q_max)
pyirena_ctrl_fit_local_power_law(session_id, q_min=seg.q_min, q_max=seg.q_max)
pyirena_ctrl_set_parameter_value(session_id, "Rg_1", rg_from_local)
pyirena_ctrl_set_parameter_value(session_id, "P_1", p_from_local)

pyirena_ctrl_fix_all_except(session_id, ["Rg_1","G_1","background"])
pyirena_ctrl_run_fit(session_id)               → chi_squared, params
pyirena_ctrl_get_fit_image(session_id)         → inspect visually
# if chi_squared > ~5, free more params or add a level, run again
pyirena_ctrl_free_parameter(session_id, "P_1")
pyirena_ctrl_run_fit(session_id)               → improved chi_squared
pyirena_ctrl_save_fit(session_id)              → write back to HDF5
pyirena_ctrl_export_fit_report(session_id)     → markdown summary
```

### Session lifecycle

| Tool | Returns |
|------|---------|
| `pyirena_ctrl_open_dataset(file_path)` | `session_id`, data summary |
| `pyirena_ctrl_list_open_sessions()` | all open sessions |
| `pyirena_ctrl_get_session_summary(session_id)` | file, model, Q range, χ² |
| `pyirena_ctrl_close_session(session_id)` | frees memory |

### Model selection

| Tool | Returns |
|------|---------|
| `pyirena_ctrl_list_available_models()` | `["unified_fit"]` (Phase 1) |
| `pyirena_ctrl_select_model(session_id, model_name="unified_fit", nlevels=1)` | full parameter table |
| `pyirena_ctrl_get_model_parameters(session_id)` | current parameter table |
| `pyirena_ctrl_get_model_description(session_id)` | physical meaning of each param + tips |

### Parameter control

All parameter names follow the convention `<param>_<level>` for
level-specific params (e.g. `Rg_1`, `G_2`, `P_1`) and plain names for
model-wide params (`background`).

| Tool | Effect |
|------|--------|
| `pyirena_ctrl_set_parameter_value(session_id, param_name, value)` | set starting value |
| `pyirena_ctrl_set_parameter_bounds(session_id, param_name, lo, hi)` | constrain range |
| `pyirena_ctrl_fix_parameter(session_id, param_name)` | hold fixed |
| `pyirena_ctrl_free_parameter(session_id, param_name)` | release for fitting |
| `pyirena_ctrl_fix_all_except(session_id, ["Rg_1", "G_1", …])` | staged fitting setup |
| `pyirena_ctrl_reset_parameters_to_defaults(session_id)` | factory defaults |

**Staged fitting strategy** (recommended):
1. `fix_all_except(["background"])` → fit background first
2. `fix_all_except(["Rg_1","G_1","background"])` → add Rg and G
3. If χ²ᵣ > 5, `free_parameter("P_1")` → add power-law slope
4. If residuals show structure at low-Q, add a level

### Level management (Unified Fit)

| Tool | Effect |
|------|--------|
| `pyirena_ctrl_add_unified_level(session_id, position=-1)` | add level (−1 = append) |
| `pyirena_ctrl_remove_unified_level(session_id, level)` | remove 1-based level |

After adding/removing a level, all parameters are renumbered and prior fit results are cleared.

### Feature detection — call before selecting a model

#### `pyirena_ctrl_detect_features(session_id, q_min=None, q_max=None, q_max_clip=0.6, config_overrides=None)`

Analyses the loaded I(Q) curve in log-log space and segments it into regions
where the power-law slope `d(log I)/d(log Q)` is approximately constant.
Returns a structured description of the curve's features without modifying
the model.

**When to call:** always call this *before* `pyirena_ctrl_select_model`.  The
return value tells you how many Unified Fit levels are needed and where to
find them, so you don't pick the wrong model complexity.

**Parameters:**
- `q_min`, `q_max` — restrict analysis to a Q sub-range (default: full data range)
- `q_max_clip` — silently drop data above this Q (default: 0.6 Å⁻¹, the SAS
  approximation limit; diffraction features above 0.6 should not be treated as
  SAS structure). Pass `null` to disable.
- `config_overrides` — dict of `FeatureDetectConfig` field names to adjust
  sensitivity.  Useful fields:
  - `"change_threshold_1"` (default 0.40) — raise to get fewer coarser segments;
    lower to catch subtler slope transitions
  - `"change_threshold_2"` (default 0.20) — threshold for the refinement pass
    inside wide (> 1.0 decade) segments
  - `"min_segment_decades"` (default 0.10) — minimum width of an interior segment
  - `"merge_slope_tol"` (default 0.15) — merge adjacent segments with similar slopes

**Returns dict with:**

| Key | Type | Description |
|-----|------|-------------|
| `ok` | bool | `true` on success |
| `segments` | list | Power-law slope segments, sorted **low-Q to high-Q** |
| `guinier_knees` | list | Inferred Guinier knees between adjacent segments |
| `recommended_nlevels` | int | Suggested number of Unified Fit levels |
| `recommended_guinier_windows` | list | Q windows for `fit_local_guinier` per knee |
| `background_q_min` | float\|null | Q where high-Q flat background begins |
| `n_segments_found` | int | Total segments detected |
| `log_decades` | float | Total Q range analysed in log10(Q) decades |
| `q_min_analysed`, `q_max_analysed` | float | Actual Q range after clipping |
| `n_points` | int | Data points used |

**Each segment has:**

| Field | Description |
|-------|-------------|
| `q_min`, `q_max` | Q range of the segment (Å⁻¹) |
| `slope` | Mean `d(log I)/d(log Q)` — negative for power-law, near zero for plateaus |
| `slope_std` | Std dev of slope within the segment (lower = more constant) |
| `kind` | `"power_law"` \| `"guinier_plateau"` \| `"background"` |
| `intensity_mid` | I at the geometric centre Q — useful initial G estimate |
| `width_decades` | Segment width in log10(Q) decades |

**Each Guinier knee has:**

| Field | Description |
|-------|-------------|
| `q_min`, `q_max`, `q_center` | Q location of the transition |
| `slope_low_q` | Mean slope of the lower-Q (shallower) segment |
| `slope_high_q` | Mean slope of the higher-Q (steeper) segment |
| `delta_slope` | `\|slope_high_q − slope_low_q\|` |

A knee is only reported when `|slope_low_q| < |slope_high_q|` — the slope
gets shallower going high-Q to low-Q, the physical Guinier-knee signature.

**Each recommended Guinier window has:**

| Field | Description |
|-------|-------------|
| `q_min_guinier` | Use as `q_min` for `fit_local_guinier` |
| `q_max_guinier` | Use as `q_max` for `fit_local_guinier` |
| `q_min_powerlaw` | Use as `q_min` for `fit_local_power_law` |

**Interpreting the result:**

1. **How many levels?** Use `recommended_nlevels`.  Each `"power_law"` or
   `"guinier_plateau"` segment with a Guinier knee on its low-Q side is one
   Unified Fit level; `"background"` segments are not levels.

2. **Starting P?** Each segment has a `P` field (positive Porod exponent,
   I ∝ Q⁻ᴾ).  Use it directly as the starting value for `set_parameter_value`.

3. **Starting Rg and G?** Call `fit_local_guinier` with the Q window from
   `recommended_guinier_windows[i]`; apply the returned Rg and G with
   `set_parameter_value`.

4. **Background starting value?** If `background_q_min` is non-null, read
   `intensity_mid` of the background segment and use it as the background
   starting value.

5. **Level ordering:** `segments` and `recommended_guinier_windows` are already
   ordered **high-Q → low-Q**, matching Unified Fit level numbering.
   `segments[0]` corresponds to Level 1 (smallest structure, highest Q);
   `segments[1]` to Level 2, and so on.  No reversal needed.

---

### Q range

| Tool | Effect |
|------|--------|
| `pyirena_ctrl_get_data_q_range(session_id)` | full data Q range |
| `pyirena_ctrl_get_fit_q_range(session_id)` | current fit Q range |
| `pyirena_ctrl_set_fit_q_range(session_id, q_min=…, q_max=…)` | restrict fit range |
| `pyirena_ctrl_reset_fit_q_range(session_id)` | restore full data range |

Use `set_fit_q_range` to exclude beam-stop artefacts at low-Q or noisy
high-Q tails before fitting. Either end can be `null` to leave it unchanged.

### Fit execution

#### `pyirena_ctrl_run_fit(session_id, max_iter=None)`
Runs the fitting algorithm synchronously. Returns:
- `success` (bool)
- `chi_squared`, `reduced_chi_squared` (float)
- `iterations` (int)
- `message` (fit status from scipy)
- `parameters_updated` (list of {name, value})
- `quality` (dict) — robust fit-quality scalars: `robust_scale_s`,
  `realistic_reduced_chi2_floor`, `max_abs_frac_misfit`, `q_at_max_frac_misfit`,
  `median_frac_uncertainty`, `n_outliers_3s`, `longest_same_sign_run`,
  `sign_autocorr_lag1`, `sigma_available`. See **Quality assessment** below for
  how to read these. (Full per-point arrays + per-band breakdown:
  `pyirena_ctrl_get_fit_quality`.)

**Interpreting reduced_chi_squared — read this carefully:**
Reported uncertainties σ in SAXS are *frequently mis-scaled*, so reduced χ² alone
is unreliable. A reduced χ² of 9 does **not** necessarily mean a bad fit — it may
just mean σ are ~3× too small. **Do not chase reduced χ² ≈ 1** and do not dismiss
a fit just because reduced χ² is large. Instead, combine it with `quality`:
- If `robust_scale_s` ≈ 1 → σ are honest; the usual reading applies (~1 excellent,
  2–10 reasonable, >10 poor, <0.5 over-fitting).
- If `robust_scale_s` ≈ s > 1 → σ are ~s× too small; the *realistic* reduced-χ²
  floor is `realistic_reduced_chi2_floor` ≈ s². A reduced χ² near that floor with
  no outliers and no sign-structure is **as good as the data allows — stop.**
- Regardless of σ scale, a large `max_abs_frac_misfit` (≳ 0.3) or a long
  `longest_same_sign_run` indicates a **real** misfit to investigate.

Re-calling `run_fit` after a partial convergence continues from current
parameter values — this is intentional and useful.

### Quality assessment

| Tool | Returns |
|------|---------|
| `pyirena_ctrl_get_chi_squared(session_id)` | `chi_squared`, `reduced_chi_squared` |
| `pyirena_ctrl_get_residuals(session_id)` | `residuals` (normalised), `rescaled_residual`, `frac_misfit_percent`, `summary` (rms / max_abs / mean / `robust_scale_s`) |
| `pyirena_ctrl_get_fit_quality(session_id, n_bands=4)` | full robust diagnostics (scalars + per-point arrays + per-band) |
| `pyirena_ctrl_get_fit_image(session_id, width=1024, height=768)` | inline PNG (data + model + residuals subplot) |
| `pyirena_ctrl_get_residuals_image(session_id)` | same image; requires completed fit |

#### `pyirena_ctrl_get_fit_quality(session_id, n_bands=4)` — robust diagnostics

The recommended way to judge a fit when σ may be mis-scaled. Returns **facts
only** (no good/bad verdict — you apply the thresholds). Fields:

**Global scalars**
- `sigma_available` (bool) — `False` if no usable σ; then σ-dependent fields are
  `null` and only the fractional/structure fields are meaningful.
- `reduced_chi2`, `dof`, `n_valid`, `n_params`.
- `robust_scale_s` — MAD-based estimate of how many × the *actual* scatter
  exceeds the reported σ. **≈ 1** σ honest; **≈ 3** σ ~3× too small; **≪ 1** σ
  too large / over-fitting risk. This is the single most useful number.
- `sigma_misscale_factor` — alias of `robust_scale_s`.
- `realistic_reduced_chi2_floor` — `robust_scale_s²`; the lowest reduced χ² the
  data can physically support. Don't try to beat it.
- `max_abs_frac_misfit` + `q_at_max_frac_misfit` — the largest |(I−M)/I| and the
  Q where it occurs. **σ-independent backstop**: ≳ 0.3 (30 %) means a gross local
  misfit no matter how unreliable σ is.
- `median_frac_uncertainty` — typical σ/I ("a few %").
- `n_outliers_3s`, `frac_outliers_3s` — points beyond 3·`robust_scale_s`. These
  are genuine outliers *even after* accounting for a mis-scaled σ.

**Structure scalars** (distinguish a wrong σ-scale from a wrong model)
- `longest_same_sign_run` — long run of same-sign (I−M) ⇒ wrong functional form.
- `sign_autocorr_lag1` — near +1 ⇒ systematic; near 0 ⇒ random scatter.

**Per-band** (`bands`, count in `n_bands_used`): each `{q_lo, q_hi, n,
reduced_chi2, robust_scale_s, max_abs_frac_misfit}`. One hot band points at the
Q-region to fix; uneven per-band χ² is itself a misfit signal.

**Decision sketch** (yours to tune):
- bulk `robust_scale_s` small-ish, no sign-runs, `max_abs_frac_misfit` < ~0.15 →
  fit is as good as the data allows; **stop** (a high reduced χ² is just σ-scale).
- `n_outliers_3s` > 0 with `max_abs_frac_misfit` ≳ 0.3, or a hot band, or a long
  sign-run → **real misfit; investigate that Q-region.**
- `robust_scale_s` ≪ 1 while still tightening → **over-fitting; back off.**

`get_fit_image` works **before and after** a fit:
- Before: shows data + model at current starting parameter values (useful for sanity-checking starting conditions)
- After: adds a normalised residuals subplot below the main plot

Handle the image return value the same way as `pyirena_plot_iq`:
check `content.type`, render the `image` item, show the `text` item as a label.

### Persistence

| Tool | Effect |
|------|--------|
| `pyirena_ctrl_save_fit(session_id, output_path=None)` | write fit to HDF5 (default: overwrites source) |
| `pyirena_ctrl_export_fit_report(session_id, format="markdown")` | returns report text |

---

## Example interactions

### "What's in this folder?"
```
1. pyirena_summarize_folder("/data/run42")
   → {"n_files": 47, "samples": ["catalyst_A", "catalyst_B"],
       "analyses_count": {"unified_fit": 47, "size_distribution": 47,
                          "modeling": 12}}
2. Tell the user in plain language.
```

### "Show me the last three I(Q) curves for catalyst_A."
```
1. pyirena_list_files("/data/run42", sort="mtime_desc", limit=3)
   → filter rows where sample == "catalyst_A"
2. pyirena_plot_iq([row["path"] for row in rows], overlay=True)
   → content[0]: text — "Plot saved to: /tmp/.../iq_xxxxx.png"
     content[1]: image — base64 PNG (render inline; do NOT print as text)
3. Display the image; mention the saved path.
```

### "Is Rg trending across this batch?"
```
1. pyirena_summarize_folder(...) — confirm unified_fit results are present.
2. pyirena_plot_parameter_trend(folder, tool="unified_fit",
                                 parameter="Rg", subgroup_index=1)
   → content[0]: text — file path
     content[1]: image — base64 PNG (render inline; do NOT print as text)
3. Describe the trend in words (slope direction, scatter, any plateau).
```

### "Compare the two most recent fits in detail."
```
1. pyirena_list_files(..., limit=2)
2. For each of the two paths:
   - pyirena_inspect_file(path) to see which analyses are present
   - pyirena_read_unified_fit(path) — get parameters
3. Present the parameters side by side and comment on differences.
```

### "Fit this dataset with Unified Fit."
```
1. pyirena_ctrl_open_dataset("/data/scan_042.h5")
   → session_id = "a1b2c3d4"

2. pyirena_ctrl_detect_features("a1b2c3d4")
   → segments=[{q_min, q_max, P, kind, intensity_mid, ...}, ...],
     guinier_knees=[{q_min, q_max, P_low_q, P_high_q, delta_P, ...}, ...],
     recommended_nlevels=2,
     recommended_guinier_windows=[
       {q_min_guinier=0.01,  q_max_guinier=0.06,  q_min_powerlaw=0.05},  # Level 1 (high Q)
       {q_min_guinier=0.001, q_max_guinier=0.008, q_min_powerlaw=0.007}, # Level 2 (low Q)
     ]
   → Note: segments and windows are ordered HIGH-Q → LOW-Q (Level 1 first).
   → recommended_nlevels=2, so use nlevels=2 below.

3. pyirena_ctrl_select_model("a1b2c3d4", nlevels=2)
4. pyirena_ctrl_get_model_description("a1b2c3d4")
   → read fitting tips

5. # segments[0] = Level 1 (high Q, small structure)
   # segments[1] = Level 2 (low Q, large structure)
   # recommended_guinier_windows follow the same order.
   local_L1 = pyirena_ctrl_fit_local_guinier("a1b2c3d4", q_min=0.01, q_max=0.06)
   → Rg=28, G=1e3
   local_L2 = pyirena_ctrl_fit_local_guinier("a1b2c3d4", q_min=0.001, q_max=0.008)
   → Rg=320, G=5e5
   pyirena_ctrl_set_parameter_value("a1b2c3d4", "Rg_1", 28)
   pyirena_ctrl_set_parameter_value("a1b2c3d4",  "G_1", 1e3)
   pyirena_ctrl_set_parameter_value("a1b2c3d4", "Rg_2", 320)
   pyirena_ctrl_set_parameter_value("a1b2c3d4",  "G_2", 5e5)

6. pyirena_ctrl_fix_all_except("a1b2c3d4", ["Rg_1","G_1","Rg_2","G_2","background"])
7. pyirena_ctrl_run_fit("a1b2c3d4")
   → reduced_chi_squared=3.8 (reasonable; Rg and G now near-optimal)
8. pyirena_ctrl_free_parameter("a1b2c3d4", "P_1")
   pyirena_ctrl_free_parameter("a1b2c3d4", "P_2")
   pyirena_ctrl_run_fit("a1b2c3d4")
   → reduced_chi_squared=9.1, but quality.robust_scale_s=3.0
   → σ are ~3× underestimated; realistic_reduced_chi2_floor≈9 — this IS converged,
     not a bad fit. Don't keep tightening.
9. pyirena_ctrl_get_fit_quality("a1b2c3d4")
   → max_abs_frac_misfit=0.08 (8%), n_outliers_3s=0, longest_same_sign_run short,
     bands all similar → no real misfit; the high χ² is purely a σ-scale artefact.
10. pyirena_ctrl_get_fit_image("a1b2c3d4")
    → confirm residuals visually (random scatter, just wide)
11. pyirena_ctrl_save_fit("a1b2c3d4")
12. pyirena_ctrl_export_fit_report("a1b2c3d4", format="markdown")
    → summarise for user
```

### "Why might Rg jump at scan 17?"
```
1. pyirena_tabulate_parameter(..., tool="unified_fit", parameter="Rg",
                                subgroup_index=1)
2. Identify the jump in the rows.
3. For scan 17 specifically:
   - pyirena_inspect_file → confirm what data is there
   - pyirena_read_metadata → sample notes, thickness, blank
   - pyirena_read_reduced_data(include_full=True) — get the actual curve
     to comment on changes in slope / Guinier region
4. Hypothesise (changed sample? exposure? beam intensity?
   contamination?). Always frame as hypotheses, not certainties.
```

---

## Things to avoid

- **Don't call `pyirena_read_*` blindly.** Use `pyirena_inspect_file`
  first to confirm the analysis is present. Otherwise you'll get
  `{"found": false}` and burn a turn.
- **Don't request full arrays unless necessary.** `include_arrays=True`
  / `include_full=True` make responses much larger and consume your
  context budget. For most questions, the scalar summary is enough.
- **Don't fabricate parameter values.** If `value` is `null`, say so.
  Do not interpolate or invent.
- **Don't speculate beyond the data.** Physical interpretations (e.g.
  "this Rg means …") should be framed as "consistent with…" and tied to
  what the user has told you about their sample.
- **Don't chase reduced χ² ≈ 1, and don't dismiss a fit for a large reduced χ².**
  Reported σ are often mis-scaled. Call `pyirena_ctrl_get_fit_quality` (or read
  the `quality` block from `run_fit`): if `robust_scale_s` ≈ 3 then σ are ~3× too
  small and reduced χ² ≈ 9 is the *realistic floor* — a fit sitting there with no
  outliers and no sign-structure is done. Conversely, a normalised residual of
  20–50 is **not** "fine because σ are unreliable": check `max_abs_frac_misfit` —
  if it is ≳ 0.3 the model is off by ≳ 30 % of the data, a real misfit.
- **Don't skip `get_fit_image` after a fit.** Visual inspection of the
  residuals subplot is the most reliable way to spot systematic
  deviations that χ² alone misses.
- **Don't free too many parameters at once.** The standard staged
  approach (background first, then Rg+G, then P, then add levels) is
  more robust than releasing everything simultaneously.
- **Don't mix read tools and control tools for the same file.** The
  read tools (`pyirena_read_unified_fit`) see what was previously saved
  in the HDF5 file. The control tools (`pyirena_ctrl_*`) operate on the
  in-memory session. Call `pyirena_ctrl_save_fit` first, then the read
  tool, if you want to compare.

---

## When pyirena tools aren't enough

- **Fitting models beyond Unified Fit** (Modeling, Size Distribution,
  Simple Fits, WAXS) — Phase 1 only covers Unified Fit in the control
  API. Use the pyirena GUI tools (`pyirena-gui`, etc.) or the headless
  batch API (`pyirena.batch.*`) for those until Phase 2 ships.
- **Instrument / beamline control** — use the user's instrument-control
  agent (e.g. EPICS / pyepics), not pyirena.
- **Data reduction** (2D → 1D, sector integration, masking) — that's
  upstream of pyirena; use AreaDetector / pyFAI / etc.
