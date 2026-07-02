# Text File Import and Automatic Cleaning

pyIrena's native container format is **NXcanSAS HDF5** (`.h5`).  All fitting
tools, result saving, and the Data Selector's result-viewer panels require
HDF5.  Text files (`.dat`, `.txt`) are a common input format from many
instruments, but cannot store fit results.

## How it works

When you select a `.dat` or `.txt` file in the Data Selector — whether to
plot it, open it in a fitting tool, or export it to ASCII — pyIrena
automatically:

1. **Reads** the text file (Q, I, and optionally dI and dQ columns).
2. **Cleans** the data (see rules below).
3. **Writes** a cleaned NXcanSAS HDF5 sibling file next to the original:
   `mydata.dat` → `mydata.h5`
4. **Uses the sibling** for all subsequent operations.

This happens silently on first use.  The sibling is **cached by modification
time**: if it already exists and is newer than the text file, it is reused
without re-reading the source.  Deleting `mydata.h5` or touching `mydata.dat`
forces reconversion.

## Cleaning rules

These rules are applied silently.  The counts are recorded inside the HDF5
file under `entry/notes/` so the provenance is never lost.

| Rule | Action |
|------|--------|
| Q ≤ 0 or non-finite Q | Point removed |
| I ≤ 0 or non-finite I | Point removed (beamstop zeros) |
| dI ≤ 0 or non-finite dI (on surviving points) | dI replaced by `I × error_fraction` |
| No dI column at all | dI synthesized as `I × error_fraction` |

**Why remove I ≤ 0?**  Instruments routinely leave zero-intensity points in
reduced data for channels behind the beam-stop.  These are invisible on a
log-intensity plot but they inject large, numerically exact zeros that corrupt
all fitting methods — the minimizer drives the model intensity to zero at those
Q values, pulling the solution completely off course.  The cleaning step
removes them before they can do damage.

**Why Q ≤ 0?**  Occasional Q = 0 rows appear in some instrument outputs (the
direct-beam channel).  They make log-Q axes undefined and must be excluded.

The `error_fraction` default is **0.05** (5 % of intensity).  You can change
it in *Data Selector → Configure → Text File Options → Generated uncertainty
fraction*.

## Naming and collision guard

The converted file is placed **next to the original text file** with the same
stem and a `.h5` extension:

```
/data/sample42.dat  →  /data/sample42.h5
```

**Collision guard**: if `sample42.h5` already exists and was *not* created by
pyIrena (it has no internal provenance marker), pyIrena falls back to
`sample42_NX.h5` and prints a one-line notice.  Your existing HDF5 file is
never overwritten.

## What you see

In the Data Selector file list you still see `sample42.dat` — the file type
selector shows text files normally.  Internally, once you click "Create
Graph", a fitting tool, or "ASCII Export", pyIrena is actually working with
the cleaned `sample42.h5` next to it.  Fit results are saved into that file
and are visible in the result-viewer panels when `sample42.h5` is selected.

## Fitting and saving results

Because tools always receive a valid NXcanSAS file path, result saving works
exactly as it does for native HDF5 data:

- **Sizes / Unified Fit / Simple Fits / Modeling** — results are appended to
  `sample42.h5` under the appropriate NXcanSAS group.
- **WAXS Peak Fit / SAXS Morph** — same.
- **ASCII Export** — converts the text file first (creating `sample42.h5`),
  then exports the raw data (and any stored fit results) as `.dat` files in
  an `ascii_export/` subfolder.

## Batch API

`pyirena.batch` functions (`fit_sizes`, `fit_unified`, etc.) apply the same
conversion automatically:

```python
from pyirena.batch import fit_sizes

# Works with text files — creates sample42.h5 and saves results into it
result = fit_sizes("sample42.dat", "pyirena_config.json")
print(result['output_file'])   # → .../sample42.h5
```

## Data Merge and Data Manipulation

These panels read text files directly (they are pre-processing tools that
*produce* NXcanSAS output).  The same cleaning rules (Q ≤ 0 removed,
I ≤ 0 removed, zero errors synthesized) are applied to their text input
in-memory before any further processing.

## Low-level API

```python
from pyirena.io.text_import import (
    clean_sas_arrays,         # clean arrays in-memory, return report
    converted_sibling_path,   # compute <stem>.h5 path (no I/O)
    ensure_nxcansas_sibling,  # full convert-once workflow
)

# Clean arrays without writing a file
Q, I, E, dQ, report = clean_sas_arrays(Q_raw, I_raw, E_raw, dQ_raw,
                                        error_fraction=0.05)
print(report)
# {'n_original': 500, 'n_kept': 493, 'n_q_removed': 1,
#  'n_i_removed': 6, 'n_e_synthesized': 0, 'error_fraction': 0.05}

# Force reconversion even if a sibling already exists
from pathlib import Path
h5_path = ensure_nxcansas_sibling(Path("mydata.dat"), force=True)
```
