# Raised buttons — plan

Branch: `feature/raised-buttons` · Status: **P1 done, P2 waiting on the role
table (§9)** · 2026-11-23

## 0. Decisions (Jan, 2026-11-23)

1. **Option B, raised QSS.** Corner radius **7 px**, bevel depth **0.35**
   (`DEFAULT_RADIUS` / `DEFAULT_BEVEL` in the module).
2. **Vendored file + sync script.** Each repo must install and run on its own.
   The global `~/.claude/CLAUDE.md` tells agents to run the sync check when they
   touch a package's GUI. Minor per-package differences are acceptable; the
   overall look must match.
3. **Merge near-duplicate colours.** Use a small set of colours, and give
   buttons with the same function the same colour in every tool (table in §9).
4. **No keyboard-focus indicator** for now. It would be another feature to
   maintain, and keyboard-only use is not expected.

Files:

- `pyirena/gui/raised_buttons.py` — the shared module (source of truth).
- `scripts/sync_raised_buttons.py` — check / refresh the vendored copies.
- `pyirena/tests/test_raised_buttons.py` — module contract + render test.
- `demo_raised_buttons.py` — side-by-side demo with live tuning and a
  benchmark. Run:
  `/opt/miniconda3/envs/pyirena/bin/python planning/raised-buttons/demo_raised_buttons.py`
  (`--bench` prints the benchmark; `--shot out.png` saves a screenshot).

## 1. Diagnosis — why pyIrena's buttons are flat

pyIrena already runs Fusion (`pyirena/gui/theme.py::apply_theme`). The flat look
comes from our own stylesheets, which overrule it:

1. The app-wide rule in `theme._base_stylesheet()` —
   `QPushButton { background-color: …; border: 1px solid …; }` — flattens
   **every** button that has no inline style.
2. About 110 inline `QPushButton { background: … }` strings in 15 panel files,
   plus 37 calls to `accent_button_css` / `soft_button_css` / `chip_button_css`,
   flatten the coloured ones.
3. 219 `QPushButton(` constructions in total; `fit_worker.py` swaps the style
   at runtime (Fit → Cancel).

Matilda, Bernardyn and AIDA set no Fusion style and no theme. They run the
native style, where any `background-color` on a button already gives a flat
box (AIDA has 8 such QSS strings, Matilda and Bernardyn have 1 each).

## 2. Recommendation — option B, raised QSS via role properties

| | A Today | **B Raised QSS** | C Fusion + palette |
|---|---|---|---|
| Looks like a button | no | **yes, strongest** | mildly |
| Sunken pressed state | colour change only | **yes, inverted bevel + text shift** | yes |
| Keeps our colours exactly | yes | **yes** (gradient is derived from them) | lighter, hue shifts |
| Same on macOS/Win/Linux | yes | **yes** (QSS paints everything) | only under Fusion |
| Works with our existing QSS | — | **yes** | no: `QWidget{color}` overrides palette text (white labels turn dark), disabled needs per-group palette |

Benchmark (300 buttons, macOS, this machine; run `--bench` yourself):

| approach | build + polish | repaint all |
|---|---:|---:|
| A today (inline per button) | 31.9 ms | 5.27 ms |
| **B raised, btnRole + one sheet** | **28.9 ms** | **5.18 ms** |
| B raised, inline per button | 56.7 ms | 5.69 ms |
| C Fusion palette tint | 25.1 ms | 7.81 ms |
| raised + `QGraphicsDropShadowEffect` | 44.2 ms | 15.10 ms |

So: gradients cost nothing when painting. What does cost something is
**per-widget `setStyleSheet`**: each one is parsed and polished separately.
Using role properties is slightly *faster* than what we do today. Shadow
effects triple the repaint cost and are ruled out.

## 3. Design

### 3.1 The shared module `raised_buttons.py`

- **Pure strings, no Qt import.** Duck-typed `set_role()` is the only function
  that touches a widget. So it drops into any package regardless of its
  `_qt.py` shim, and it is testable without a display.
- API:
  - `stylesheet(roles=None, *, bevel=0.35, radius=7)` — the app-level chunk:
    neutral default for `QPushButton`/`QToolButton`, one rule set per role,
    size variants.
  - `set_role(button, role, size=None)` — sets `btnRole` / `btnSize` dynamic
    properties and re-polishes that one widget (works for runtime swaps).
  - `inline_qss(role_or_hex, size=…)` — self-contained fallback for a widget
    in an app that did not install the sheet.
  - `button_qss`, `size_qss`, `lighten/darken/mix` — building blocks.
- Roles keep pyIrena's existing colours and names: `neutral, green, mint,
  blue, teal, purple, orange, red, grey, soft-green, soft-amber`. Sizes:
  `normal, small` (the old chip), `large` (Fit buttons).
- States: raised (gradient + two-tone border), hover (lighter), pressed
  (inverted gradient and border, label 1 px down, outer size unchanged),
  checked (sunken, selection blue), disabled (flat grey, no bevel).
- Host passes its own role table (pyIrena: from `theme.py` tokens), so colour
  stays owned by each package's theme and the bevel logic stays shared.

### 3.2 pyIrena integration

`theme._base_stylesheet()` replaces its "buttons" block with
`raised_buttons.stylesheet(roles=<built from theme tokens>)`. Panels call
`set_role(btn, "green")` instead of `btn.setStyleSheet(...)`.

## 4. Distribution across packages

**Recommended: vendored single file, byte-identical, synced by script.**

- Source of truth: `pyirena/gui/raised_buttons.py` (header says so, carries
  `__version__`).
- `scripts/sync_raised_buttons.py` in pyIrena copies it to the known
  locations (`Matilda/matilda/gui/`, `Bernardyn/bernardyn/gui/`,
  `Aida/src/aida/ui/qt/`, …) and with `--check` reports any copy that has
  drifted.
- Each repo gets one small test: the module imports without Qt, and every
  role defines a text colour.
- No new dependency, no install-order coupling, each package still builds
  alone. It is ~250 lines that change rarely.

Alternative: a tiny shared package on PyPI/conda (`jil-qtstyle` or similar).
That is cleaner if more shared GUI code follows, but it adds a dependency to
every package and a release step. Not chosen (decision 2, §0).

## 5. Implementation phases (pyIrena)

Each phase is a separate commit on `feature/raised-buttons` and leaves the GUI
working.

**P1 — module + global default — DONE 2026-11-23**
- Move `raised_buttons.py` into `pyirena/gui/`.
- `theme.py`: button block → `raised_buttons.stylesheet(...)`; rebuild
  `accent_button_css` / `soft_button_css` / `chip_button_css` on
  `inline_qss`, so their 37 call sites turn raised with no call-site edits.
- Tests: `tests/test_raised_buttons.py` (no Qt import in the module, valid
  QSS colours, every role has a text colour, pressed outer size == normal);
  extend `test_gui_theme_contract.py` to render one raised button off-screen.
- Result: every unstyled button and every helper-styled button is raised.

**P2 — migrate inline strings to roles (mechanical, one file per commit)**
- Replace the ~110 raw `QPushButton{background:…}` strings and the module
  constants (`_BTN_GREEN`, `_BTN_FIT`, …) with `set_role(...)`. Map
  near-duplicate colours (`#3498db` vs `#2980b9`, `#81c784`, `#e74c3c`) to
  the nearest role. This changes button colours very slightly and nothing
  numerical.
- `fit_worker.py`: Fit ↔ Cancel via `set_role(btn, "red")` / restore role.
- Then convert the P1 helpers' callers to `set_role` and deprecate the
  helpers.
- Contract test: no new inline `QPushButton` stylesheet that sets a
  background. Known offenders go in a list that may shrink and never grow,
  the same pattern as `KNOWN_BARE_BACKGROUNDS`. Colour-swatch buttons stay
  exempt.

**P3 — sister packages (after pyIrena has been in use for a while)**
- Copy the file via the sync script. In each app's startup:
  `app.setStyle("Fusion")` (recommended) +
  `app.setStyleSheet(existing + raised_buttons.stylesheet())`.
- Matilda and Bernardyn: almost nothing to migrate. AIDA: 8 strings.
- Docs: one paragraph in each repo's AGENTS.md / developer notes.

## 6. Performance rules (go into the module docstring and AGENTS.md)

1. Install the button stylesheet **once** on the `QApplication`.
2. Choose colours with `set_role` (dynamic property), not per-widget
   `setStyleSheet`.
3. No `QGraphicsDropShadowEffect` / `QGraphicsEffect` on buttons.
4. Never call `setStyleSheet` from a paint path, timer or per-data-point loop.
5. No broad `QWidget { background: … }` rules. They flatten every child button.
   Set panel backgrounds through the palette.

## 7. Risks and things to verify

- **Embedded panels**: Matilda imports pyIrena. If it opens pyIrena panels
  in its own `QApplication` without `apply_theme`, role-only buttons fall
  back to Qt default. Fix: Matilda installs the same sheet (P3), or
  `apply_theme` is called when a panel opens. Needs checking.
- **Keyboard focus**: with QSS, Qt no longer draws a focus ring. Today's look
  has the same issue. A `:focus` border tint is optional; say if you want it.
- **Contrast**: white on `mint` (#52c77a) is about 2.3:1, an existing issue.
  We could darken mint slightly while we are touching these. Your call.
- **High-DPI Windows (125/150 %)**: the 1 px two-tone border may blur. Check
  the demo on a Windows beamline machine before P2.
- Dynamic-property changes need a re-polish. `set_role` does this, but
  hand-written `setProperty("btnRole", …)` calls would not, so the contract
  test should flag them.

## 8. P1 notes

- Normal/large padding stays at the old `3px 8px`. A wider 10 px clipped
  fixed-width buttons ("Open…" in Sizes), so it was reverted.
- Role rules only ever *add* `font-weight: bold`. Writing `normal` would
  override bold fonts that panels set in code.
- Inline-styled buttons (≈110 strings) still look flat until P2 converts them.
- The demo's column A is a frozen copy of the old look, so the comparison
  survives the theme change.

## 9. Button roles — colour by function (proposal for P2)

Built from an audit of 131 styled buttons in `pyirena/gui`. Most functions
already have one colour everywhere; the rows marked ★ change what a user
sees.

| Role | Function | Buttons (examples) |
|---|---|---|
| `green` | Main action: fit, compute, commit | Fit, Fit sizes, Fit all, Grow, Calculate 3D, Optimize merge, Save merged data, Apply & save, Export to Igor h5xp, Add to active graph |
| `mint` | Put something on the graph | Graph model / Graph unified; ★ Results to graph(s) (now `#81c784`) |
| `blue` | Secondary tools and output files | Identify features, Batch run, Calculate (contrast), Detect segments, Check similarity, New graph, Save report, Save to HDF5/NeXus, Export/Import to file; ★ Find peaks, Import CIF, Sort by Q (`#3498db`), ★ Export/Import file (`#2471a3`) |
| `teal` | Uncertainty, Monte Carlo, tabular results | Calc. uncertainty (MC), Compute MC I(q), Find best growth, Copy results, Save CSV, Add all to library; ★ Collect from all selected files (now teal *or* purple) |
| `purple` | Helper sub-fits, special analyses | Fit background, Fit P/B, Ask AI advisor, Energy scan; ★ Fit flat / power-law bckg in SAXS Morph (now mint); ★ Fix limits? (green in Modeling/Unified, purple in WAXS) |
| `orange` | Undo / reset | Revert back, Reset to defaults, Clear all CIFs |
| `red` | Stop / destructive | Cancel fit / MC / job, Remove, ×, Auto-reject; ★ `#e74c3c` merged into `#c0392b` |
| `soft-green` | Settings to/from file | Store in File, Save/Load params to JSON, + Add peak; ★ Store in File in WAXS (now mint); ★ Save/Load JSON config in Merge (now grey), Save config in SAXS Morph (now purple) |
| `soft-amber` | Load a saved setup | Load setup from file… |
| `neutral` | Utility, navigation, dialogs | Open…, Browse, Select/Create folder, Parse ▸, All/None, Open output folder, Clear markers, Pop out; ★ everything now grey `#7f8c8d` with white text (Set vacuum, dialog Cancel, …) |

That leaves **10 colours** (`sage` and `grey` are dropped from pyIrena's
table; the shared module keeps `grey` for other packages). Disabled is the one
shared flat grey; the per-file disabled shades (`#bdc3c7`, `#95a5a6`, …) go.

Open point in this table:

- **"? Help" is red in all 13 panels.** It is consistent, but red otherwise
  means stop/destructive. Options: keep red (familiar), or `blue`.
