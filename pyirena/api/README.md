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
