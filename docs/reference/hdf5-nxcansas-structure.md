# HDF5 NXcanSAS Structure Reference

This document provides instructions for agents to locate reduced Intensity ($I$), $Q$, and uncertainties ($dI$) within a pyirena-generated HDF5 file following the NXcanSAS standard.

## 1. Locating Primary Data (Reduced Intensity)

The primary scattering data is stored within the `entry/` group. To find the correct sample data, you must first identify the default sample name.

### Step 1: Identify Sample Name
Read the `default` attribute of the `entry/` group.
- **Path:** `entry/`
- **Attribute:** `@default`
- **Value:** This will be the name of the primary (desmeared) sample.

### Step 2: Access the `sasdata` Group
Once you have the sample name (e.g., `<sample_name>`), navigate to its `sasdata` group.
- **Path:** `entry/<sample_name>/sasdata/`

### Step 3: Read Datasets
The following datasets contain the reduced scattering data:

| Dataset | Description | Units | Attributes to Note |
|---------|-------------|--------|---------------------|
| `Q`     | Scattering vector | $1/\text{\AA}$ | `@long_name`: "Q (A^-1)", `@resolutions`: see below |
| `Qdev` / `dQw` | Per-point $Q$ resolution (one value per $Q$ point) | $1/\text{\AA}$ | Optional |
| `dQl`   | Slit length (half-height) — a **single scalar**, not per-point | $1/\text{\AA}$ | Present ⇒ the curve is slit smeared |
| `I`     | Reduced Intensity | $1/\text{cm}$ | `@long_name`: "Intensity", `@uncertainties`: "Idev" |
| `Idev`  | Uncertainties ($dI$) | $1/\text{cm}$ or $\text{cm}^2/\text{cm}^3$ | `@long_name`: "Uncertainties" |

**Note on $dI$:** If the `Idev` dataset is missing from the `sasdata` group, it means uncertainties were not recorded or are unavailable for that specific dataset.

**Note on `Q@resolutions`:** this attribute lists *every* resolution
contribution, comma separated — `"Qdev"`, `"dQw,dQl"`, or just `"dQl"`.
Only the per-point entries (`Qdev` / `dQw`) are arrays parallel to `Q`;
`dQl` is a single scalar describing the slit geometry. **Never read the
first token blindly** — slit-smeared data with no per-point width
legitimately declares `resolutions="dQl"` alone, and treating that scalar
as a per-point resolution corrupts everything downstream. pyIrena reads the
first non-`dQl` token into `dQ` (validating that it is parallel to `Q`) and
`dQl` separately into `slit_length` / `is_slit_smeared`.

## 2. Handling USAXS (SMR vs. Desmeared)
For USAXS files, two variants of the sample group may exist:
1. **Desmeared (Primary):** `entry/<sample_name>/`
2. **Slit-Smeared (SMR):** `entry/<sample_name>_SMR/`

**Rule:** Always prefer the desmeared variant unless specifically instructed otherwise. The `entry/@default` attribute automatically points to the desmeared variant.

## 3. Summary of Paths for Agents
To extract the core data, use these relative paths from the root:

- **Sample Name:** `entry/@default`
- **Q Array:** `entry/<sample_name>/sasdata/Q`
- **I Array:** `entry/<sample_name>/sasdata/I`
- **dI Array:** `entry/<sample_name>/sasdata/Idev`

## 4. Example Python Snippet (h5py)
```python
import h5py

with h5py.File("file.h5", "r") as f:
    # 1. Get the default sample name
    sample_name = f["entry"].attrs["default"]
    
    # 2. Access the sasdata group
    sd = f[f"entry/{sample_name}/sasdata"]
    
    # 3. Extract data
    Q  = sd["Q"][:]
    I  = sd["I"][:]
    dI = sd["Idev"][:] if "Idev" in sd else None
```

## 5. Results groups and the embedded setup (`_pyirena_config`)

Each analysis tool writes its results to `entry/<tool>_results` (check
`pyirena/io/nxcansas_<tool>.py` — the names are conventional, not derived;
Simple Fits writes `entry/simple_fit_results`).

Most tools additionally embed the **panel setup** as a JSON string in the
`_pyirena_config` attribute on that group. That is what makes "reopen the
result file and get every control back" work, including the things the physics
parameters do not capture — which check boxes were ticked, the cursor Q range,
bounds, the fit method. `pyirena/io/setup_config.py` writes and reads it, and
`pyirena/gui/setup_loader.py` provides the shared *Load Setup from File…*
dialog.

Not every tool does this, and the differences are deliberate:

| Tool | Setup embedded? | Why |
|---|---|---|
| Unified Fit | yes | fit flags, bounds, links and cursors are not recoverable from the fitted values |
| Size Distribution | yes | same |
| Modeling | yes | same |
| Simple Fits | yes | same |
| WAXS Peak Fit | yes | same |
| SAXS Morph | yes | added in this release; the config scalars stored beside the result describe the calculation but not the input mode or the cursor Q range |
| Fractals | **no** | a visualization tool rather than an analysis one. Every growth parameter that defines an aggregate is written as an explicit dataset (`StickingProbability`, `NumberOfTestPaths`, `Seed`, …) and `load_fractal_aggregate()` rebuilds the `GrowthConfig` from them, so the aggregate *is* fully restorable — there is simply no panel state worth preserving on top of it |
| Data Merge | **no — provenance instead** | merging runs inside data-reduction pipelines, driven by its own JSON config. What matters afterwards is what was done to the data, so the tool appends an NXprocess **provenance** group naming both inputs and every merge parameter. A setup blob would describe a GUI session that never existed in a pipeline run |
| Data Manipulation | **no — provenance instead** | same reasoning: the operation and its parameters are recorded as NXprocess provenance, read back by `pyirena_read_manipulation_provenance` |

**Rule for a new tool:** if a user can change a control that changes the result
but is not itself written out as a number, embed the setup. If the tool is a
step in a pipeline whose inputs are a config file, write provenance instead.
Whichever you choose, add a row here — the asymmetry above looked like an
oversight until it was written down.

### 5.1 The embedded setup is meant to be read by hand

Users open result files in HDFView, in Igor and with `h5dump`, and read
`_pyirena_config` as text. That is a constraint on its design, not an
accident of it:

* **One JSON string attribute on the results group.** Not a binary blob, not
  a sub-group per control. Anything that needs a tool to decode is not
  readable at the beamline.
* **The key names are the ones in the GUI and in the maths.** `Rg`, `G`, `P`,
  `RgCO`, `phi` — a scientist reading the attribute should recognise the
  Irena names. This is why the core field names were kept rather than renamed
  to the panel's `RgCutoff` / `correlated` / `estimate_B`.
* **It says what it is.** The header carries a `_note` — *"pyIrena tool
  setup. Replay it with pyirena.batch, or with
  pyirena.api.control.analyze(data, config)."* — alongside `tool`,
  `saved_by`, `saved_at` and `schema_version`, so whoever finds the attribute
  does not have to guess what to do with it.

The envelope is `{"_pyirena_config": {header}, "state": {the tool's
settings}}`, and `pyirena.core.tool_config.build_setup()` reads it directly —
the same reader that takes an *Export Parameters* sidecar and an
`export_results()` reply. So the setup inside a result file can be replayed
without unwrapping it by hand:

```python
import h5py, json
from pyirena.core.tool_config import build_setup

with h5py.File("sample.h5") as f:
    envelope = json.loads(f["entry/unified_fit_results"].attrs["_pyirena_config"])
setup = build_setup(envelope)     # -> model, fitted Q range, slit settings
```

Settings are written in the **core** dialect (`model.to_dict()`); the older
panel dialect, in every file written before 1.2, is still read. See
`docs/batch_api.md` § *Configuration file format* for both, and
`planning/config-dialects/` for why there were two.
