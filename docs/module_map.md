# Module map

One line per module, so you can find a thing without grepping 113,000 lines.
This is a *lookup table*, not documentation — for how the layers fit together
read [AGENTS.md](../AGENTS.md) §2, and for how to wire something new read
[developer_adding_features.md](developer_adding_features.md).

Per-tool modules (`core/<tool>.py`, `io/nxcansas_<tool>.py`,
`gui/<tool>_panel.py`, `batch/<tool>.py`) follow the tool map in AGENTS.md §3
and are listed here only for completeness. **The tables worth reading are the
shared ones** — those are what gets reimplemented by accident.

---

## "Where is the code that…"

| …does this | look here |
|---|---|
| decides how a table copies, sorts, exports CSV | `gui/table_utils.py` |
| exports a plot (clipboard, PNG/SVG, CSV, Igor ITX) | `gui/plot_export.py` |
| picks colours / makes a widget readable on any OS theme | `gui/theme.py` |
| remembers window size, position, splitter widths | `gui/window_state.py` |
| handles files dropped onto a window | `gui/file_drop.py` |
| filters a file list | `gui/file_filter.py` |
| sorts filenames (numeric, date, name) | `core/file_sorting.py` |
| decides which extensions are data files | `core/file_types.py` |
| draws the standard I(Q) log-log plot | `gui/sas_plot.py` |
| provides the editable Q-min/Q-max fit-range control | `gui/q_range_ui.py` |
| builds a Markdown fit report (one builder, all consumers) | `core/reporting.py` |
| adds the "Copy results / Save report…" buttons | `gui/report_buttons.py` |
| shows the fit-quality summary line | `gui/quality_display.py` |
| imports Qt (the only place allowed to) | `gui/_qt.py` |
| persists panel settings between sessions | `state/state_manager.py` |
| embeds GUI setup inside a result file | `io/setup_config.py` |
| reloads that setup ("Load Setup from File…") | `gui/setup_loader.py` |
| describes result groups to generic consumers | `io/schema.py` (`TOOL_REGISTRY`) |
| maps pyIrena results to Igor wave names | `io/igor_names.py` (`TOOL_CROSS_REF`) |
| formats numbers for display | `core/fmt_utils.py` |
| writes the log file | `logging_setup.py` |
| diagnoses a broken install (`pyirena-doctor`) | `diagnostics.py` |

---

## core/ — math and serialisable model objects (numpy/scipy only)

**Per tool**

| Module | Lines | Purpose |
|---|---:|---|
| `unified.py` | 1687 | Unified Fit model (Beaucage levels) |
| `sizes.py` | 1824 | Particle size-distribution fitting (MaxEnt, Reg, TNNLS, MCSaS) |
| `modeling.py` | 2698 | Modeling tool engine — multi-population forward model |
| `simple_fits.py` | 1706 | Simple Fits model library (`MODEL_REGISTRY`) |
| `waxs_peakfit.py` | 1418 | WAXS profiles, backgrounds, peak finder, fit engine |
| `fractals.py` | 1097 | Mass-fractal aggregate generator and analyzer |
| `saxs_morph.py` | 1401 | 3-D two-phase voxelgram model via Gaussian random fields |
| `scattering_contrast.py` | 681 | Scattering-contrast calculator engine |
| `data_merge.py` | 888 | Merge two SAS datasets |
| `data_manipulation.py` | 579 | Headless trim / scale / average / subtract / rebin |

**Shared physics and support**

| Module | Lines | Purpose |
|---|---:|---|
| `form_factors.py` | 952 | Form factors for size-distribution analysis |
| `distributions.py` | 630 | Parametric size distributions for Modeling |
| `smearing.py` | 463 | Slit-smearing engine |
| `fit_metrics.py` | 361 | Robust, σ-scale-independent fit diagnostics |
| `feature_detect.py` | 684 | Feature detection in I(Q) — drives the Feature Identifier |
| `similarity.py` | 307 | Outlier / radiation-damage detection across datasets |
| `morphology.py` | 306 | Morphology metrics for binary 3-D voxelgrams |
| `diffraction_lines.py` | 201 | Theoretical powder-diffraction stick patterns |
| `reporting.py` | 754 | **The** Markdown fit-report builder (GUI, API and batch all use it) |
| `file_sorting.py` | 186 | Filename sort keys for every file browser |
| `file_types.py` | 114 | Data-file type table and folder listing |
| `fmt_utils.py` | 68 | Number formatting |

## io/ — NXcanSAS HDF5 persistence and format conversion

| Module | Lines | Purpose |
|---|---:|---|
| `nxcansas_unified.py` | 573 | Unified Fit results ↔ `entry/unified_fit_results` |
| `nxcansas_sizes.py` | 373 | Size-distribution results |
| `nxcansas_modeling.py` | 404 | Modeling results |
| `nxcansas_simple_fits.py` | 361 | Simple Fits results (group is `entry/simple_fit_results`) |
| `nxcansas_waxs_peakfit.py` | 397 | WAXS peak-fit results |
| `nxcansas_saxs_morph.py` | 367 | SAXS Morph results |
| `nxcansas_fractals.py` | 237 | Fractal aggregate results |
| `nxcansas_data_merge.py` | 188 | Data Merge provenance |
| `nxcansas_data_manipulation.py` | 187 | Data Manipulation provenance |
| `nxcansas_fit_quality.py` | 161 | Shared `fit_quality/` sub-group for every residual-based tool |
| `contrast_io.py` | 373 | Contrast calculator's own compound library |
| `_nxcansas_common.py` | 233 | Shared write helpers; `PYIRENA_RESULT_GROUPS` |
| `schema.py` | 477 | **`TOOL_REGISTRY`** — machine-readable description of every result group |
| `setup_config.py` | 186 | Embed/read the `_pyirena_config` GUI setup blob |
| `results.py` | 351 | High-level "load whatever results this file has" |
| `scattering.py` | 203 | Read-only data discovery and loading API |
| `text_import.py` | 313 | Clean and convert ASCII SAS files to NXcanSAS |
| `ascii_export.py` | 1111 | ASCII export |
| `hdf5.py` | 1321 | Low-level HDF5 support shared with Matilda |
| `pxp_to_nexus.py` | 1382 | Import legacy Igor `.pxp` packed experiments |
| `h5xp_writer.py` | 662 | Write Igor Pro `.h5xp` packed experiments |
| `h5xp_extractor.py` | 940 | Extract pyIrena results into an Igor experiment |
| `igor_names.py` | 376 | **`TOOL_CROSS_REF`** — pyIrena ↔ Igor/Irena wave-name mapping |

## batch/ — headless execution from a dict or JSON config

| Module | Lines | Purpose |
|---|---:|---|
| `pipeline.py` | 157 | `fit_pyirena` — runs every tool present in a `pyirena_config.json`; holds `_TOOL_REGISTRY` |
| `unified.py` | 512 | `fit_unified` |
| `sizes.py` | 424 | `fit_sizes` |
| `simple.py` | 321 | `fit_simple`, `fit_simple_from_config` |
| `modeling.py` | 229 | `fit_modeling` |
| `waxs.py` | 295 | `fit_waxs_peaks` |
| `saxs_morph.py` | 194 | `fit_saxs_morph` |
| `merge.py` | 245 | `merge_data` |
| `manipulate.py` | 288 | `manipulate_data`, `average_data` |
| `convert.py` | 153 | Igor `.pxp` / `.h5xp` → NXcanSAS |
| `_common.py` | 117 | Shared config and data loading |

## api/ — stable JSON-serialisable facade (no Qt, no numpy scalars)

**Read side**

| Module | Lines | Purpose |
|---|---:|---|
| `results.py` | 581 | Per-tool result readers |
| `data.py` | 142 | Raw I(Q) and sample metadata |
| `discovery.py` | 289 | List files, summarise a folder, inspect one file |
| `aggregate.py` | 310 | Parameter trends and sample summaries across files |
| `plotting.py` | 225 | Headless matplotlib (Agg) rendering |
| `schemas.py` | 517 | Result schemas |
| `_paths.py` | 88 | `PYIRENA_DATA_ROOT` sandboxing for every read and write |

**Write / control side** (`api/control/` — interactive fitting sessions)

| Module | Lines | Purpose |
|---|---:|---|
| `unified_fit.py` | 2057 | Agent-drivable Unified Fit |
| `modeling.py` | 1096 | Agent-drivable Modeling |
| `sizes.py` | 893 | Agent-drivable Size Distribution |
| `waxs_peakfit.py` | 898 | Agent-drivable WAXS Peak Fit |
| `simple_fits.py` | 710 | Agent-drivable Simple Fits |
| `schemas.py` | 1800 | **`TOOL_SCHEMA_BY_NAME`** — JSON schema per control function |
| `session.py` | 88 | In-memory session registry |
| `errors.py` | 45 | Structured error helpers |
| `_images.py` | 61 | Shared PNG rendering |

**Stateless helpers exposed to agents**

| Module | Lines | Purpose |
|---|---:|---|
| `calculators.py` | 689 | Scattering contrast and friends (+ `calculator_schemas.py`) |
| `data_ops.py` | 1217 | Merge and manipulation for agents (+ `data_op_schemas.py`) |

## mcp/ — protocol wrapper

| Module | Lines | Purpose |
|---|---:|---|
| `server.py` | 631 | MCP stdio server; the only place that builds MCP content blocks |
| `dispatch.py` | 203 | Collapses ~90 control functions into 4 dispatcher tools (provider tool cap) |

## gui/ — one panel per tool, plus the shared UX contract

**Tool panels**

| Module | Lines |
|---|---:|
| `unified_fit.py` | 4442 |
| `sizes_panel.py` | 3195 |
| `modeling_panel.py` | 3895 |
| `waxs_peakfit_panel.py` | 2706 |
| `simple_fits_panel.py` | 2352 |
| `data_manipulation_panel.py` | 2370 |
| `data_merge_panel.py` | 1999 |
| `saxs_morph_panel.py` | 1896 |
| `fractals_panel.py` | 1569 |
| `contrast_panel.py` | 1558 |

**Shared behaviour — use these, never reimplement**

| Module | Lines | Purpose |
|---|---:|---|
| `_qt.py` | 159 | The single PySide6/PyQt6 import point |
| `theme.py` | 484 | `apply_theme`, colour tokens, background+foreground helpers |
| `table_utils.py` | 557 | Copy, numeric sort, CSV for every `QTableWidget` |
| `plot_export.py` | 883 | Clipboard / image / CSV / ITX export for every plot |
| `sas_plot.py` | 1147 | The standard I(Q) plot (`make_sas_plot`, `plot_iq_data`) |
| `window_state.py` | 779 | Window geometry, splitter widths, Shift-click reset |
| `file_drop.py` | 328 | Drag-and-drop file opening |
| `file_filter.py` | 104 | The shared filter box |
| `q_range_ui.py` | 289 | Editable Q-range fields tied to graph cursors |
| `report_buttons.py` | 347 | "Copy results" / "Save report…" |
| `quality_display.py` | 85 | Uniform fit-quality readout |
| `setup_loader.py` | 145 | "Load Setup from File…" dialog flow |
| `data_loading.py` | 313 | Shared loader row widget |
| `slit_smearing_ui.py` | 207 | Shared "Slit smeared data" row |
| `ai_advisor.py` | 865 | In-GUI LLM advisor panel |
| `launch.py` | 56 | `pyirena-gui` entry point |

**Auxiliary dialogs**

`feature_identifier.py` (457), `sizes_feature_identifier.py` (115),
`diffraction_lines_panel.py` (752, embedded in WAXS), `saxs_morph_3d.py`
(1034), `fractals_workers.py` (272).

**`gui/data_selector/` — the Data Browser, and the main entry point**

| Module | Lines | Purpose |
|---|---:|---|
| `panel.py` | 3009 | `DataSelectorPanel` — launches every tool; `main()` lives here |
| `results_windows.py` | 1035 | Stored-fit viewer windows and the tabulate window |
| `config_dialogs.py` | 472 | Configure… / Manage Config… dialogs |
| `plot_utils.py` | 141 | Shared pyqtgraph palette and styling |
| `igor_import.py` | 119 | Igor import dialog |
| `workers.py` | 89 | Background batch-fitting threads |
| `report.py` | 14 | Thin shim onto `core/reporting.py` |
| `sorting.py` | 36 | Thin shim onto `core/file_sorting.py` |

**`gui/hdf5viewer/` — the Data Explorer (trend plots, Igor export)**

| Module | Lines | Purpose |
|---|---:|---|
| `main_window.py` | 590 | `HDF5ViewerWindow` |
| `plot_controls.py` | 1096 | Tabbed control panel; **hardcoded per-tool "collect" lists** |
| `pyirena_readers.py` | 730 | Knows where each tool stores results in HDF5 |
| `hdf5_browser.py` | 470 | Lazy HDF5 tree browser |
| `file_tree.py` | 393 | Folder/file tree with lazy sub-folder expansion |
| `graph_window.py` | 893 | Floating pyqtgraph plot window |
| `export.py` | 381 | Export utilities for `GraphWindow` |
| `export_to_igor_tab.py` | 327 | "Export to Igor" tab |
| `collect_window.py` | 323 | Collected 0-D values as table + scatter |
| `multi_collect_window.py` | 178 | Multi-value collection table |

## Everything else

| Module | Lines | Purpose |
|---|---:|---|
| `state/state_manager.py` | 930 | Schema-versioned JSON settings; per-tool default blocks |
| `plotting/unified_plots.py` | 409 | Headless matplotlib for Unified Fit |
| `plotting/plot_saxs.py` | 395 | `plot_saxs` multi-file CLI plotting utility |
| `diagnostics.py` | 609 | `pyirena-doctor` |
| `logging_setup.py` | 166 | Rotating log configuration |
| `version_check.py` | 79 | GitHub-based update notification |
| `examples/basic_demo.py` | 390 | Unified Fit demonstration script |
| `gui/ai_skills/*.md` | — | Tool-specific prompt context for the AI advisor |

---

## Maintaining this file

Regenerate rather than hand-edit when it drifts: every line comes from the
module's own first docstring line. A module with no docstring is the reason a
row here reads badly — fix the docstring, not the table.

```bash
python - <<'PY'
import ast, pathlib
for p in sorted(pathlib.Path('pyirena').rglob('*.py')):
    if '__pycache__' in p.parts or 'tests' in p.parts:
        continue
    doc = ast.get_docstring(ast.parse(p.read_text(errors='replace'))) or ''
    first = ' '.join(doc.strip().split('\n\n')[0].split())
    print(f'{str(p):50s} {sum(1 for _ in p.open(errors="replace")):5d}  {first[:110]}')
PY
```
