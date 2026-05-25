# Importing Igor Pro Packed Experiments (.pxp)

The **Import Igor Experiment** tool lets you bring data from a legacy Igor
Pro packed experiment (`.pxp`) into pyIrena for analysis. Each USAXS, SAXS,
or WAXS sample in the experiment becomes a stand-alone NXcanSAS `.h5`
file that any pyIrena tool (Unified Fit, Size Distribution, Modeling,
…) can open directly.

The tool is available three ways:

| Mode             | How to launch                                                          |
|------------------|------------------------------------------------------------------------|
| Interactive GUI  | **Data Processing & Reference → Import Igor Experiment…** in the Data Selector |
| Python API       | `from pyirena.batch import pxp_to_nexus`                                |
| CLI              | `python -m pyirena.io.pxp_to_nexus legacy.pxp -v`                       |

---

## Contents

1. [When to use this tool](#when-to-use-this-tool)
2. [What gets exported](#what-gets-exported)
3. [GUI walkthrough](#gui-walkthrough)
4. [Output folder layout](#output-folder-layout)
5. [Python API](#python-api)
6. [CLI reference](#cli-reference)
7. [Wave-name and folder-name conventions](#wave-name-and-folder-name-conventions)
8. [Wave-note metadata](#wave-note-metadata)
9. [Limitations](#limitations)
10. [Troubleshooting](#troubleshooting)

---

## When to use this tool

Use it whenever your reduced 1-D scattering data lives in an Igor Pro
experiment (`.pxp` or `.pxt`) and you want to analyse it in pyIrena
without going back to raw detector files.

Typical scenarios:

- **Legacy archives** — years of USAXS/SAXS/WAXS data reduced with the
  original Igor pipeline, stored as one `.pxp` per beamtime.
- **Quick triage** — somebody hands you an Igor experiment with 30
  samples; you want to fit Unified or size distributions to all of them
  in one go without manual export.
- **Round-tripping** — combine with the existing **Export to Igor**
  (h5xp writer) for Python → Igor → Python workflows.

The reverse direction — exporting pyIrena results *into* an Igor `.h5xp`
packed experiment — is documented in
[`hdf5_viewer_gui.md`](hdf5_viewer_gui.md#export-to-igor).

## What gets exported

For each sample folder the tool recognises (see [conventions](#wave-name-and-folder-name-conventions)),
one NXcanSAS `.h5` file is written containing:

| Dataset                        | Source wave                                   | Notes                                  |
|--------------------------------|-----------------------------------------------|----------------------------------------|
| `entry/<sample>/sasdata/Q`     | `DSM_Qvec` (USAXS) or `R_Qvec` (SAXS/WAXS)    | units = `1/angstrom`                  |
| `entry/<sample>/sasdata/I`     | `DSM_Int`  or `R_Int`                          | units = `1/cm`                        |
| `entry/<sample>/sasdata/Idev`  | `DSM_Error` or `R_Error`                       | always written                        |
| `entry/<sample>/sasdata/Qdev`  | `DSM_dQ`                                       | written only when the source wave exists |
| `entry/sample/*`               | parsed from `NXSampleStart…End` in the note    | when present                          |
| `entry/instrument/*`           | parsed from `NXInstrumentStart…End`            | wavelength is hoisted to `entry/instrument/beam/incident_wavelength` |
| `entry/notes/*`                | everything else from the wave note             | preserved verbatim                    |

Samples whose folders **don't** contain the expected wave triple are
skipped silently — typically these are "blank" / "raw" folders that
were never desmeared or reduced.

## GUI walkthrough

1. Open the Data Selector (`pyirena-gui`).
2. Click **Import Igor Experiment…** in the *Data Processing &
   Reference* group (purple button).
3. Pick the `.pxp` file in the file dialog.
4. In the import dialog:
   - **Output folder** — defaults to `<pxp_stem>_data` next to the input.
   - **Techniques to export** — three checkboxes (USAXS / SAXS / WAXS),
     all on by default.
   - **Overwrite existing files** — off by default; an unused suffix
     (`_2`, `_3`, …) is appended to keep both copies if the target
     filename exists.
5. Click **Import**. Extraction runs synchronously (a 16-MB legacy
   experiment with 50+ samples finishes in 1–2 s on a typical laptop).
6. A summary dialog reports the per-technique tally and offers to load
   the output folder as the current data folder, so you can start
   analysing immediately.

## Output folder layout

A 60-sample APS USAXS experiment with all three techniques produces:

```
legacy_2024_data/
├── USAXS/
│   ├── Sample01.h5
│   ├── Sample02.h5
│   └── ...
├── SAXS/
│   ├── Sample01.h5
│   └── ...
└── WAXS/
    ├── Sample01.h5
    └── ...
```

If your Igor experiment uses nested sub-folders for in-situ runs
(e.g. `root:SAXS:heater_run_42:step_03:Sample01`), the importer
**mirrors** that depth in the output:

```
legacy_2024_data/
└── SAXS/
    └── heater_run_42/
        └── step_03/
            └── Sample01.h5
```

This preserves the organisation you already chose in Igor — useful when
the folder names encode experimental conditions you care about.

## Python API

```python
from pyirena.batch import pxp_to_nexus

result = pxp_to_nexus(
    pxp_file="legacy.pxp",
    output_folder=None,         # default: <pxp_stem>_data next to input
    techniques=["USAXS"],       # default: None = all present
    overwrite=False,            # default: append _2, _3, …
    verbose=True,
)

print(f"Wrote {result['n_written']} files to {result['output_folder']}")
for f in result['files']:
    if f['status'] != 'ok':
        print(f"  {f['status']}: {f['source']} — {f['message']}")
```

Returned dict:

| Key                       | Type     | Meaning                                                                 |
|---------------------------|----------|-------------------------------------------------------------------------|
| `success`                 | bool     | always `True` if the function returns                                    |
| `output_folder`           | str      | absolute path of the output folder                                       |
| `n_written`               | int      | number of `.h5` files successfully written                               |
| `n_skipped`               | int      | sample folders that didn't have a recognised wave triple                 |
| `n_errors`                | int      | sample folders that failed during file write                             |
| `n_unparseable_records`   | int      | wave records in the `.pxp` that igor2 couldn't decode (usually 0–1)      |
| `files`                   | list     | per-folder dicts with `source`, `output`, `technique`, `n_points`, `status`, `message` |

Returns `None` only if the input path does not exist.

## CLI reference

```
python -m pyirena.io.pxp_to_nexus PXP [-o OUTPUT] [-t TECHNIQUE] [--overwrite] [-v]
```

Options:

| Flag              | Meaning                                                                |
|-------------------|------------------------------------------------------------------------|
| `-o, --output`    | Output folder. Default: `<pxp_stem>_data` next to the input.           |
| `-t, --technique` | Only export this technique. Repeat for multiple (e.g. `-t USAXS -t SAXS`). |
| `--overwrite`     | Overwrite existing output files instead of appending `_2`, `_3`, …     |
| `-v, --verbose`   | Print a per-file summary table.                                        |

Exit status is 0 if all files wrote successfully, 1 otherwise.

## Wave-name and folder-name conventions

The recognition tables are **data-driven dicts at the top of
`pyirena/io/pxp_to_nexus.py`** — you can extend them in-place if your
group uses different conventions.

### Top-level folder → technique

```python
TECHNIQUE_FOLDERS = {
    "USAXS":         "USAXS",
    "SAXS":          "SAXS",
    "WAXS":          "WAXS",
    "Imported SAXS": "SAXS",
    "SAS":           "SAXS",
    "Imported":      "SAXS",
}
```

Matching is case-insensitive. The folder name **inside** Igor (e.g.
`root:Imported SAXS:...`) is what's checked.

### Per-technique wave triples

```python
WAVE_PICKERS = {
    "USAXS": [("DSM_Qvec", "DSM_Int", "DSM_Error", "DSM_dQ")],
    "SAXS":  [("R_Qvec",   "R_Int",   "R_Error",   None)],
    "WAXS":  [("R_Qvec",   "R_Int",   "R_Error",   None)],
}
```

Each entry is `(Q_name, I_name, Err_name, dQ_name_or_None)`. The
importer tries the entries in order; the first one whose Q, I, and
Err waves all exist in a sample folder wins. dQ is optional — if the
named wave isn't there, `Qdev` is simply omitted from the NeXus file.

To add a new pattern (e.g. exporting slit-smeared USAXS as well), append
a tuple:

```python
WAVE_PICKERS["USAXS"].append(
    ("SMR_Qvec", "SMR_Int", "SMR_Error", "SMR_dQ"),
)
```

## Wave-note metadata

Igor wave notes use the convention `key=value;key=value;`. The APS USAXS
pipeline wraps NeXus-style metadata in sentinel markers:

```
DATAFILE=03_31_run.dat;DATE=2026-03-31 12:22:05;COMMENT=Sample01;
Nexus_attributesStartHere;
NXUserStart;...;NXUserEnd;
NXSampleStart;name=Sample01;thickness=4;temperature=20;NXSampleEnd;
NXInstrumentStart;name=APS USAXS;wavelength=0.5904;NXInstrumentEnd;
NXMetadataStart;...;NXMetadataEnd;
Nexus_attributesEndHere;
SlitLength=0.024;NumberOfSteps=100;...
```

The parser maps these to NXcanSAS as follows:

| Marker pair                                  | Lands in              |
|----------------------------------------------|-----------------------|
| `NXSampleStart` … `NXSampleEnd`              | `entry/sample/<key>`  |
| `NXInstrumentStart` … `NXInstrumentEnd`     | `entry/instrument/<key>` |
| `NXUserStart` … `NXUserEnd`                  | `entry/notes/user_<key>` |
| `NXMetadataStart` … `NXMetadataEnd`         | `entry/notes/<key>`   |
| `DATAFILE=...` (top-level)                   | `entry/title` + `entry/notes/DATAFILE` |
| `DATE=...` (top-level)                       | `entry/start_time` + `entry/notes/DATE` |
| `COMMENT=...` (top-level)                    | `entry/title` (only if no DATAFILE) |
| Any other top-level `key=value`              | `entry/notes/<key>`   |

The `instrument.wavelength` value is also hoisted into the canonical
NXcanSAS location `entry/instrument/beam/incident_wavelength` so that
downstream tools find it at the spec-mandated path.

Bare wave-level annotations (`units=...`, `long_name=...`, `Wname=...`)
are dropped — they describe the wave itself, not the sample or
instrument.

## Limitations

- **Reduced 1-D data only.** 2-D detector frames, calibration tables,
  and Igor procedure files are ignored. The goal is to bring
  analysis-ready I(Q) into pyIrena, not to fully round-trip an entire
  Igor experiment.
- **No analysis-result import.** Fits stored in Igor (Unified, Size
  Distribution, etc.) are *not* parsed. Re-fit them in pyIrena using
  the imported NeXus files — that way the results are stored in the
  pyIrena schema and can be browsed by the HDF5 Viewer.
- **One known igor2 quirk.** The upstream `igor2` package occasionally
  fails to decode a single wave record in older Igor experiments
  (typically a text or dependency-formula wave). The importer is
  defensive — it skips that record and reports it in
  `n_unparseable_records` rather than aborting, so the rest of the
  experiment still extracts cleanly.

## Troubleshooting

**"0 files written" but the experiment definitely has samples** —
your folder layout probably doesn't match `TECHNIQUE_FOLDERS`. Open
the `.pxp` in Igor and check the top-level folder names. If they're
e.g. `MyData_USAXS`, `MyData_SAXS`, add them to the dict.

**Sample folders show up as skipped** — they don't contain the
expected wave triple. Likely causes:

- It's a "blank" folder that holds raw detector counts only (intended
  behaviour — these are correctly skipped).
- The reduction pipeline used a different wave-name convention. Add
  a tuple to `WAVE_PICKERS` for the relevant technique.

**"igor2 is required" ImportError** — install it: `pip install igor2`
(or `pip install -e ".[gui]"` to get the full pyIrena stack).
`igor2` is a hard dependency now and is pulled in automatically by
both `pip install pyirena` and `environment.yml`.

**Imported files don't appear in the Data Selector** — the importer
already calls the same NXcanSAS locator (`find_matching_groups`) that
pyIrena's reader uses, so this shouldn't happen. If it does, file an
issue with one of the offending `.h5` files attached.
