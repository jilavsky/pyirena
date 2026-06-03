# pyirena MCP Tools — Reference for AI Agents

This document is written for AI assistants and the humans who craft
system prompts for them. It explains every tool the `pyirena-mcp` server
exposes, when to use which, and recommended call patterns.

If you are an AI assistant with access to pyirena tools, read this
document — or have it inlined into your system prompt — before asking
questions about a user's analysis data.

> **All MCP tool names are prefixed `pyirena_`** (e.g.
> `pyirena_summarize_folder`, `pyirena_list_files`). The prefix makes
> tools globally unambiguous when a client connects to multiple MCP
> servers, and avoids confusing small models with the `server-tool`
> dash convention some clients render. The underlying Python library
> functions in `pyirena.api` are unprefixed (`api.summarize_folder`,
> `api.list_files`) — only the MCP wrappers carry the prefix.

---

## What pyirena gives you

pyirena reads small-angle X-ray scattering (SAXS / USAXS) analysis
results stored in NXcanSAS HDF5 files. The user has typically run one or
more analyses on each file:

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

Your job is to help the user understand what's in their data and what it
means physically.

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
  first to confirm the analysis is there. Otherwise you'll get
  `{"found": false}` and burn a turn.
- **Don't request full arrays unless necessary.** `include_arrays=True`
  / `include_full=True` make responses much larger and consume your
  context budget. For most questions, the scalar summary is enough.
- **Don't fabricate parameter values.** If `value` is `null`, say so. Do
  not interpolate or invent.
- **Don't speculate beyond the data.** Physical interpretations (e.g.
  "this Rg means …") should be framed as "consistent with…" and tied to
  what the user has told you about their sample.
- **Don't try to write or modify files.** v0.7 is read-only. There are
  no `pyirena_save_*`, `pyirena_update_*`, or `pyirena_run_*` tools.

---

## When pyirena tools aren't enough

If a user asks a question that pyirena's surface cannot answer (running
a new fit, modifying a file, doing instrument control), say so
explicitly. Suggest:

- For new fits / data manipulation: the pyirena GUI tools (
  `pyirena-gui`, `pyirena-modeling`, `pyirena-datamerge`, etc.) or the
  headless batch API (`pyirena.batch.fit_unified`,
  `pyirena.batch.merge_data`, etc.).
- For instrument / beamline control: the user's instrument-control
  agent (e.g. the EPICS / pyepics tool), not pyirena.
