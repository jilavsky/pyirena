# pyirena.api

Stable, AI-friendly facade over pyirena's HDF5 readers. Returns
JSON-serializable dicts; no Qt/pyqtgraph imports.

## Quick reference

```python
from pyirena import api

# Discovery
api.summarize_folder("/data/run42")        # counts of samples, analyses
api.list_files("/data/run42", limit=20)    # one row per file
api.inspect_file("/data/run42/scan_007.h5")

# Reading
api.read_reduced_data("scan_007.h5")
api.read_unified_fit("scan_007.h5")        # arrays decimated by default
api.read_modeling("scan_007.h5")

# Aggregation
api.tabulate_parameter("/data/run42",
                       tool="unified_fit", parameter="Rg",
                       subgroup_index=1)
api.summarize_sample("/data/run42", sample="catalyst_3")

# Plotting (headless matplotlib)
api.plot_iq(["scan_007.h5", "scan_008.h5"], output_path="/tmp/iq.png")
api.plot_parameter_trend("/data/run42",
                          tool="unified_fit", parameter="Rg",
                          subgroup_index=1)

# Calculators (stateless: no dataset, no file access)
api.calc_contrast("TiO", 4.95, "Ti2O3", 4.49)   # -> xray_contrast 11.63
api.calc_compound("SiO2", 2.2)                  # -> xray_sld 18.8
api.calc_contrast_energy_scan("TiO", 4.95, "Ti2O3", 4.49,
                              e_start_keV=4.5, e_end_keV=5.5)
api.lookup_element("Ti")
api.list_compound_library(); api.load_compound("Alumina")

# Data operations (the only group that WRITES data files)
api.average_data(["s1.h5", "s2.h5", "s3.h5"])      # -> run_manip/s1_avg.h5
api.subtract_data("sample.h5", "buffer.h5")        # -> run_manip/sample_sub.h5
api.merge_datasets("usaxs/s_001.h5", "saxs/s_001.h5")
api.scale_data("s1.h5", scale_I=2.0)
api.trim_data("s1.h5", q_min=0.01, q_max=0.3)
api.rebin_data("s1.h5", n_points=200)
api.divide_data("s1.h5", "reference.h5")
api.match_merge_files("usaxs/", "saxs/")           # read-only pairing helper
```

## Environment overrides

| Variable | Purpose | Default |
|----------|---------|---------|
| `PYIRENA_DATA_ROOT` | Restrict all file access to this subtree | none (any abs path OK) |
| `PYIRENA_MAX_ARRAY_POINTS` | Default decimation cap for returned arrays | 500 |
| `PYIRENA_PLOT_CACHE` | Where plot PNGs are saved when no path given | `<tempdir>/pyirena-mcp` |

## Tool keys

`tabulate_parameter()` and `summarize_sample()` use the keys from
`pyirena.io.schema.TOOL_REGISTRY`:

- `simple_fits`, `unified_fit`, `size_distribution`
- `modeling`, `saxs_morph`, `waxs_peakfit`
- `fractals`, `data_merge`, `data_manipulation`

Each tool's available scalar parameters are listed in
`TOOL_REGISTRY[tool]["scalars"]`. Parameters marked `per_subgroup: True`
require a `subgroup_index` argument (1-based).

## Design notes

- All functions return `dict` (via `dataclasses.asdict`). Missing data
  returns `{"found": False, ...}` rather than raising.
- Numpy arrays in returned dicts are decimated to bounded length and
  NaN/inf are replaced with `None` so output is strict-JSON-safe.
- Set `include_arrays=True` on result readers to keep arrays (still
  decimated to `max_points`). Set `include_full=True` on
  `read_reduced_data` for full-fidelity I(Q).
- **Calculators are the exception to the two rules above.** They answer
  experiment-planning questions from first principles rather than reading a
  file, so they take no path, never return `found`, and are the only api
  module outside the `PYIRENA_DATA_ROOT` sandbox (there is nothing to
  sandbox — the sole file touched is the user's own compound library, read
  only). They also *return* `{"error", "suggestion", "code"}` dicts instead
  of raising, matching `pyirena.api.control`: they are reached through the
  same MCP dispatcher, and their failure modes (a typo in a formula, an
  unknown element, a missing optional dependency) are agent-recoverable.
  They need the `pyirena[contrast]` extra.
- `calc_contrast()` returns `xray_contrast` as (Δρ)² in 10²⁰ cm⁻⁴ — the same
  units and convention as the `contrast` parameter of a Sizes `set_shape()`
  or a Modeling population, so the result can be fed straight into a fit.
- **Data operations are the only api group that creates data files.** Output
  goes to a *sibling* of the source folder (`/data/run42` →
  `/data/run42_manip`, or `_merged`) with a per-operation filename suffix
  (`_avg`, `_sub`, `_div`, `_scaled`, `_trimmed`, `_rebinned`, `_merged`),
  matching what the GUI and `pyirena.batch` already do. Repeating an
  operation overwrites its own previous result; different operations do not
  collide. Pass `output_folder` to override — required when
  `PYIRENA_DATA_ROOT` is the data folder itself, since the sibling then
  falls outside it (you get `PATH_NOT_ALLOWED`).
- Every data operation returns `n_points_in` / `n_points_result` /
  `n_points_written` / `n_dropped_nonpositive`. Points vanish for two
  legitimate reasons — log-log interpolation does not extrapolate, and the
  saver strips non-positive intensities — so check these rather than
  assuming the output has the same length as the input. A subtraction that
  drops many points is over-subtracted.
- These wrap `core` + `io` directly rather than `pyirena.batch`, which
  attaches a stdout log handler (fatal for MCP's stdio transport) and
  returns bare `None` on every failure.
