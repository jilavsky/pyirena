# Feature Parity & Uniformity Review — 2026-08-08

**Purpose.** Find user-expected features that are missing or unevenly implemented across
pyIrena's user-contact surfaces, and backend non-uniformities that will cause the next
"regex filter" incident. Review only — no code was changed. Each item is meant to be
implemented as its own small, debuggable session with an agent.

**Method.** Systematic grep/read audit of `pyirena/gui`, `pyirena/core`, `pyirena/io`,
`pyirena/batch`, `pyirena/api`. Every claim below was verified against the code on this
date; file:line references are included so an implementing agent can re-verify before
touching anything (the codebase will drift).

**Context.** Users come from Igor Pro, Origin, Excel, Matlab. The baseline expectation is:
any table can be copied to the clipboard, any graph can be copied/exported, any list can
be sorted, results can be gotten out as text without a file round-trip, and every panel
that looks the same behaves the same.

**Good news first.** Two feared problems are already solved: the filename filter is
unified in `pyirena/gui/file_filter.py` (regex with substring fallback, shared tooltip,
headless tests in `pyirena/tests/test_file_filter.py`) and is used by all four file
browsers. HDF5 group naming is centralized in `pyirena/io/schema.py`. The layering
invariant (no Qt in core/io/api/batch) holds — verified with grep, zero hits.

---

## Priority key

- **P1** — users will hit this quickly; the "clipboard incident" class. Do first.
- **P2** — expected by Igor/Origin/Excel users; will surface as support email eventually.
- **P3** — polish / consistency; batch up opportunistically.
- **U**  — backend unification; do before adding new tools.

Effort: **S** = under ~50 lines touching 1–2 files; **M** = one focused session;
**L** = multi-session, needs a plan.

---

# Part A — Missing user-facing features

## A1. Clipboard copy for every table (P1, S–M)

Exactly two of the nine table widgets support copy today, with two *different*
implementations:

| Table | File | Copy? | Notes |
|---|---|---|---|
| Data Selector tabulate | `gui/data_selector/results_windows.py:965–1065` | ✅ | Context menu + Ctrl+C + Ctrl+Shift+C (with headers). The reference implementation. |
| Contrast T-dependence table | `gui/contrast_panel.py:833, 904–930` | ✅ | Own bespoke copy; menu-only Ctrl+C shortcut (verify it actually fires when menu is closed — QAction shortcut created inside the context-menu handler typically does not). |
| Contrast isotope table | `gui/contrast_panel.py:662` | ❌ | Sits next to a table that *has* copy. |
| HDF5 Viewer collect window | `gui/hdf5viewer/collect_window.py:100` | ❌ | Has "Save CSV" only (`:78–81`) — the exact scenario from the user report: forced CSV round-trip for a few numbers. |
| HDF5 Viewer multi-collect | `gui/hdf5viewer/multi_collect_window.py:94` | ❌ | Same: CSV only (`:80–83`). |
| Data Manipulation similarity results | `gui/data_manipulation_panel.py:936` | ❌ | No copy **and no CSV export at all** — results are view-only (see A10). |
| HDF5 Viewer multi-collect config table | `gui/hdf5viewer/plot_controls.py:371` | ❌ | Config table, low value; include only if the shared helper makes it free. |

**Recommendation.** Extract one helper (e.g. `pyirena/gui/table_utils.py`) modeled on
`results_windows.py`: `attach_table_copy(table, headers=True)` installing the context
menu, Ctrl+C, and Ctrl+Shift+C-with-headers, tab-separated (pastes cleanly into Excel,
Igor, Origin). Then apply it to every `QTableWidget` above and delete the two bespoke
implementations. This is the same play as `file_filter.py` — one module, one behavior,
one tooltip, one test file.

Verify with: `grep -rn "QTableWidget(" pyirena/gui` and confirm each hit calls the helper.

## A2. Column-click sorting on tables (P2, S)

`setSortingEnabled` appears **nowhere** in the codebase. Every mainstream tool sorts a
table when you click the header; Igor tables, Excel, Origin all do. Most valuable on the
Data Selector tabulate window (many rows of fit results) and the HDF5 Viewer
collect/multi-collect windows. One line per table (`table.setSortingEnabled(True)`) plus
a check that numeric columns sort numerically (use `QTableWidgetItem.setData(Qt.EditRole,
float)` or a numeric-aware item subclass — text items will sort "10" before "9", which
would be its own bug report). Fold into A1's shared helper.

## A3. Copy graph to clipboard + more image formats (P1, M)

Every plot's export menu today is "Save graph as JPEG…" (file dialog only) and, on some
plots, "Save as Igor Pro ITX…". Missing, and expected from Igor (Edit→Copy), Origin, and
Matlab users:

1. **Copy graph to clipboard** as an image. In pyqtgraph:
   `QApplication.clipboard().setImage(ImageExporter(plot).export(toBytes=True))` — a few
   lines, huge daily value (paste straight into PowerPoint/Word/email).
2. **PNG option** (JPEG artifacts on line plots with text are visibly ugly; PNG should
   arguably be the default). SVG is cheap via `pyqtgraph.exporters.SVGExporter` if wanted.
3. The JPEG save dialog defaults to `Path.home()` (`gui/sas_plot.py:777`) instead of the
   remembered last folder that data dialogs use via StateManager (see A8).

**Where:** implement once in the shared export helper (see U2 — there are currently ~5
duplicated implementations) so every plot gains all three at once.

## A4. Export plotted curves as CSV/text, not only ITX (P2, S–M)

`save_itx_from_plot` (`gui/sas_plot.py:856`) exports plot curves as Igor .itx only.
Origin/Excel/Matlab/Python users need CSV/tab-text of the same curves. Add "Save curve
data as CSV…" beside the ITX entry in the same menu, reusing the same curve-collection
code that the ITX exporter already has. Goes into the shared export helper (U2) so it
appears on every plot uniformly.

## A5. Fit results as text: copy/report button on every fit panel (P1, M)

Verified: no fit panel (Unified, Sizes, WAXS, Simple Fits, Modeling, Fractals) can put
its numeric results on the clipboard or into a text report. The only exits are
save-to-HDF5 (then Data Selector tabulation) or reading numbers off the widgets. Igor
Irena writes every fit to the notebook; users expect a text trail. Meanwhile the
API already has exactly this: `pyirena_ctrl_export_fit_report` exists for agents but not
for GUI users.

**Recommendation.** One shared "Results report" mechanism: a small function per tool that
formats the current model into the same text the ctrl-API report uses (put the formatter
in core or api so GUI/batch/MCP all share it), plus two buttons on each fit panel: "Copy
results" (clipboard) and "Save report…" (.txt). Implement the shared plumbing once, then
one small PR per panel.

## A6. Drag-and-drop file opening (P2, M)

`setAcceptDrops` / `dropEvent` appear nowhere. Igor, Origin, and every modern analysis
GUI accept a file dropped onto the window/panel. Highest value: dropping data files onto
the Data Selector file list (behave like Select Folder + select those files) and onto the
HDF5 Viewer. Moderate Qt boilerplate, well-contained per widget. Do it in the shared file
browser if U4 happens first; otherwise start with Data Selector only.

## A7. Params-to-JSON parity + config naming (P2, S each)

"Save params to JSON / Load params from JSON" exists on: Unified (`unified_fit.py:2404`),
Sizes (`sizes_panel.py:1679`), Modeling (`modeling_panel.py:2605`), Simple Fits
(`simple_fits_panel.py:950`), WAXS (`waxs_peakfit_panel.py:1565`), SAXS Morph (config
export, `saxs_morph_panel.py:886`). Missing on:

- **Fractals** (`fractals_panel.py` — zero `json` hits). No JSON save/load, and also no
  batch module (U8), so fractal setups cannot be reproduced or scripted at all.
- **Data Manipulation** (zero `json` hits) — even though `batch/manipulate.py` exists, the
  GUI cannot export the config that batch would consume.
- **Contrast** — no JSON (has its own HDF5 compound library via `contrast_io.py`, so
  arguably fine; decide deliberately rather than by omission).

Also: Data Merge saves `merge_config.json` (`data_merge_panel.py:1449`) while every other
tool appends a section to `pyirena_config.json`. Unify on the `pyirena_config.json`
section convention so one file drives the whole batch pipeline.

## A8. File-dialog folder memory is inconsistent (P3, S)

StateManager-based last-folder memory exists and works for data loading
(`gui/data_loading.py:239–247`), Data Selector (`panel.py:95,1019`), Diffraction Lines
(`diffraction_lines_panel.py:458`). But graph JPEG export defaults to `Path.home()`
(`sas_plot.py:777`), and other save dialogs vary. Users notice when one dialog remembers
and the next doesn't. Sweep all `QFileDialog.get(Save|Open)FileName` calls to use one
shared `last_folder` accessor. Cheap; combine with A3/U2.

## A9. Window geometry not persisted (P3, S–M)

No `QSettings`, no `saveGeometry`/`restoreGeometry` anywhere; windows open at hard-coded
sizes (e.g. `unified_fit.py:406`). Igor and Origin users expect windows to reopen where
they left them. Since StateManager already round-trips setup state, add
geometry (as `x,y,w,h`) per panel to its saved state — or a one-shot `QSettings`-based
mixin. Low urgency, high perceived polish.

## A10. Similarity results are trapped in the GUI (P2, S)

The Data Manipulation similarity results table (`data_manipulation_panel.py:936`) has no
copy, no CSV, no save. A scientist who computes similarity across 50 files can only look
at it. Fix = A1 helper (copy) + a "Save CSV" button matching the collect windows.

## A11. Plots with no export menu at all (P2, S each)

Panels whose plots verifiably lack both JPEG and ITX export: **Data Merge**, **Data
Manipulation** (`DataManipulationGraphWindow`, `data_manipulation_panel.py:338`),
**Fractals**, **Feature Identifier**. Simple Fits has JPEG only — no ITX
(`simple_fits_panel.py:150–199`). Once U2 exists, wiring each is a few lines; users
currently discover the gap when they try to right-click "like in the Unified panel."

---

# Part B — Cross-surface parity matrix

Quick reference; ✅ verified present, ❌ verified absent, ◐ partial.

| Feature | Data Selector | HDF5 Viewer | Manip | Merge | Unified | Sizes | Modeling | Simple | WAXS | Fractals | SAXS Morph | Contrast |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Regex file filter (shared) | ✅ | ✅ | ✅ | ✅ | – | – | – | – | – | – | – | – |
| Named sort modes (temp/time/order…) | ✅ | ✅ | ✅ | ◐ order-only | – | – | – | – | – | – | – | – |
| Table clipboard copy | ✅ | ❌ | ❌ | – | – | – | – | – | – | – | – | ◐ one of two tables |
| Table CSV export | ✅ | ✅ | ❌ | – | – | – | – | – | – | – | – | ❌ |
| Plot JPEG export | ✅ | ✅ (mpl) | ❌ | ❌ | ✅ | ✅ | ✅ | ✅ | ✅ | ❌ | ✅ | ✅ |
| Plot ITX export | ✅ | ✅ (Igor tab) | ❌ | ❌ | ✅ | ✅ | ✅ | ❌ | ✅ | ❌ | ? verify | ? verify |
| Copy graph to clipboard | ❌ everywhere | | | | | | | | | | | |
| Column-click table sorting | ❌ everywhere | | | | | | | | | | | |
| Drag-and-drop open | ❌ everywhere | | | | | | | | | | | |
| Params/config JSON save+load | – | – | ❌ | ◐ own filename | ✅ | ✅ | ✅ | ✅ | ✅ | ❌ | ✅ | ❌ |
| Results text report / copy | ◐ via tabulate | – | – | – | ❌ | ❌ | ❌ | ❌ | ❌ | ❌ | – | – |
| Setup save/restore (StateManager) | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |

---

# Part C — Backend unification (do before adding tools)

## U1. Shared table utilities module (enables A1, A2, A10)

Create `pyirena/gui/table_utils.py` following the `file_filter.py` pattern: Qt-touching
but self-contained, one behavior, one test file (Qt parts smoke-tested, formatting logic
headless). Contents: `attach_table_copy()`, numeric-sort-aware item helper,
`table_to_csv()` (shared by every "Save CSV" button so quoting/precision is identical).
Delete the bespoke copies in `results_windows.py` and `contrast_panel.py` once migrated.

## U2. One plot-export implementation (enables A3, A4, A11)

There are at least **five** parallel implementations of "add JPEG/ITX export to a plot":
`sas_plot.py:750` (the most complete), `waxs_peakfit_panel.py:650`,
`unified_fit.py:706–736`, `data_selector/plot_utils.py:122`, and inline
`ImageExporter` uses in `modeling_panel.py:2088`, `contrast_panel.py:204`,
`sizes_panel.py:443`, `simple_fits_panel.py:194`. Any new export feature (clipboard
copy, PNG, CSV) currently has to be added five+ times — the regex-filter failure mode,
already in progress. Consolidate on the `sas_plot.py` version, extend it with
clipboard/PNG/CSV once, and migrate every panel to call it. This is the single highest-
leverage refactor in this review.

## U3. De-duplicate the filename sort keys

`_sort_key_name/temperature/time/order/pressure` exist verbatim in **three** places:
`gui/data_selector/sorting.py`, `gui/hdf5viewer/file_tree.py:50–98`, and
`gui/data_manipulation_panel.py:89–126`, plus a lone `_sort_key_order` in
`data_merge_panel.py:96`. A new sort key (e.g. `_percent`, `_wt%`, voltage) or a regex
fix will be added to one and silently missed in the others. Promote
`data_selector/sorting.py` to a shared module (e.g. `pyirena/gui/list_sorting.py`, or
Qt-free under `pyirena/` so batch/api can reuse it), import it everywhere, delete copies.
Bonus: give Data Merge the full sort-mode dropdown the other three browsers have.

## U4. One file-browser widget (longer term)

Four independent implementations of "folder picker + filter box + file list":
Data Selector (`data_selector/panel.py`), `_ManipFileBrowser`
(`data_manipulation_panel.py:171`), `_DatasetSelectorWidget` (`data_merge_panel.py:687`),
and `hdf5viewer/file_tree.py`. They already share filter logic; they don't share sort
modes, context menus, double-click behavior, or (future) drag-and-drop. A shared
`FileBrowserWidget` would make A6 and every future browser feature a one-place change.
**L effort** — do after U1–U3, panel by panel, comparing behavior carefully (Data
Selector's browser has the most features; treat it as the spec).

## U5. Standardize panel state-method names

Current spread, all meaning the same two operations:
save: `save_state` (unified, sizes, simple, manip, merge) / `_save_state` (modeling,
waxs, contrast, saxs_morph, fractals, diffraction, hdf5viewer);
load: `load_state` / `_load_state` / `_restore_state` (contrast) / `_apply_state`
(waxs, unified);
collect: `collect_state` / `_collect_state`.
Pick the public pair `save_state()` / `load_state()` (CLAUDE.md §5 already mandates
exposing them) with private `_collect_state()` / `_apply_state(dict)` as the
implementation halves, rename mechanically, and record the convention in
`docs/developer_adding_features.md`. Pure rename, no logic change — good candidate for
one supervised agent pass with the full test suite as the gate.

## U6. `to_dict`/`from_dict` rollout in core

Present: `simple_fits.py`, `sizes.py`, `morphology.py`, `feature_detect.py` (to only).
Absent from core model objects of: `unified.py`, `modeling.py`, `fractals.py`,
`waxs_peakfit.py`, `saxs_morph.py`, `data_merge.py`, `data_manipulation.py`,
`scattering_contrast.py` — their serialization lives ad hoc in gui/io layers instead of
on the model ("all state lives in the core model object" is the stated golden rule).
Not urgent per tool, but each new agent-built feature on a tool without `to_dict` will
invent its own serialization. Suggested order: `unified.py` first (most used, and its
GUI/JSON/HDF5 code already implies the dict shape — formalize it), then modeling, waxs.
Every addition needs a round-trip test and `from_dict` defaults for old files.

## U7. Consolidate the Qt import shims

Invariant says all Qt goes through `gui/_qt.py`, but there are three shims —
`gui/_qt.py`, `gui/data_selector/_qt.py`, and an inline one in `slit_smearing_ui.py:31` —
plus scattered local `try: PySide6 / except: PyQt6` blocks in `sizes_panel.py:750`,
`unified_fit.py:227`, `data_loading.py:253`, `ai_advisor.py`, `fractals_panel.py`,
`saxs_morph_panel.py`, `hdf5viewer/export_to_igor_tab.py`. Add the missing symbols to
`gui/_qt.py`, redirect everything, and add a test/lint check (grep in CI) so it stays
fixed.

## U8. Batch-layer gaps: Fractals (and decide on Contrast)

`pyirena/batch/` covers every tool except **Fractals** (no batch, no JSON — a fractals
analysis is currently unreproducible outside the GUI) and **Contrast** (a calculator; a
thin batch/API "compute contrast for these compounds" would still serve beamline scripts
and agents). Fractals batch is the real gap given the beamline use case.

## U9. api/control coverage for agents

Interactive control sessions (`api/control/`, exposed via MCP) exist for **Unified Fit
and Sizes only**. Read access covers all tools (`api/results.py` has readers for all
ten), but the stated direction — agents helping users or running independently — needs
control for Modeling, Simple Fits, and WAXS Peak Fit at minimum. The session/Q-range
infrastructure in `api/control/session.py` is shared, so each new tool is a bounded **M**
task following the `sizes.py` template. Suggested order by scientific demand: Modeling →
Simple Fits → WAXS.

## U10. `_pyirena_config` GUI-state embedding parity in io

Embedded by: unified, sizes, modeling, simple_fits, waxs_peakfit. Not embedded by:
**fractals, saxs_morph, data_merge, data_manipulation** — so "reopen results file, restore
panel exactly" works for the first five tools and silently doesn't for the rest. For
merge/manipulation the provenance groups may be a deliberate substitute; for fractals and
saxs_morph it looks like plain omission. Decide per tool and document the decision in
`docs/HDF5_NxcanSAS_structure.md`.

---

# Part D — Suggested sequencing

Each line = one agent session with tests. Order chosen so shared helpers land before the
many small consumers.

1. **U1** table_utils (copy + numeric sort + CSV) → then **A1/A2/A10** across all tables
   (one panel per session).
2. **U2** plot-export consolidation → then **A3** (clipboard/PNG), **A4** (CSV curves),
   **A11** (missing menus), **A8** (dialog folder memory) ride along.
3. **A5** results report: shared formatter (reuse ctrl-API report code) → per-panel buttons.
4. **U3** shared sort keys (+ Merge sort dropdown).
5. **A7** JSON parity: Fractals, Data Manipulation; Merge filename convention.
6. **U8** Fractals batch module (unblocks reproducibility).
7. **U5** state-method rename (mechanical, one pass).
8. **U9** api/control for Modeling, then Simple Fits, then WAXS.
9. **A6** drag-and-drop; **A9** geometry persistence; **U7** Qt shims; **U6**/**U10**
   per-tool as those tools are next touched.
10. **U4** shared file browser — last, once its feature set (filter+sort+copy+DnD+menus)
    is fully settled by the items above.

# Part E — Preventing recurrence: the "standard UX contract"

The root cause of both incidents was that "obvious" features live in nobody's spec. Add a
short checklist to `docs/developer_adding_features.md` (and reference it from CLAUDE.md)
that every new GUI surface must satisfy before it is done:

Every **table**: clipboard copy (tab-separated, with/without headers), column sorting,
CSV export if >1 row of results. Every **plot**: right-click export menu from the shared
helper (clipboard, PNG/JPEG, ITX, CSV). Every **file list**: shared filter
(`file_filter.py`), shared sort modes, remembered folder. Every **panel**:
`save_state`/`load_state`, params-to-JSON as a `pyirena_config.json` section, tooltips on
buttons. Every **result**: reachable via HDF5 *and* text report *and* `pyirena.api`.

With the shared helpers from Part C in place, each checklist line is one function call —
which is what makes the checklist realistic for agents to follow.

---

*Generated by code review on 2026-08-08 against version 1.1.0b6. Re-verify file:line
references before implementing; several panels exceed 2000 lines and shift often.*
