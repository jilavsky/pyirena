# P2 brief — migrate inline button stylesheets to function roles

You are converting push-button styling in a few pyIrena GUI files. Other
agents are converting *other* files at the same time in the same working tree:
**edit only the files assigned to you, and do not run any git command.**

## Background (already done, do not change)

- `pyirena/gui/raised_buttons.py` — shared module that draws raised buttons.
- `pyirena/gui/theme.py` — `apply_theme()` installs one app-wide stylesheet
  with a rule per *role*. Panels choose a role with
  `set_button_role(button, role, size=None)` (import from
  `pyirena.gui.theme`). This only sets a dynamic property, so it is cheap and
  can be called again at runtime.
- Read `theme.py` around `BUTTON_FUNCTIONS` and `set_button_role`, and
  `pyirena/gui/fit_worker.py` (swaps Fit → Cancel via roles) for the pattern.

## The task

Replace every `QPushButton` / `QToolButton` stylesheet that sets a
**background colour** with `set_button_role(...)`, in your files only:

- inline strings (`btn.setStyleSheet("QPushButton { background: #27ae60; ... }")`),
- module or class constants (`_BTN_GREEN = "QPushButton {...}"`) — convert
  their users and delete the constant when it is no longer used,
- the old helpers `accent_button_css(...)`, `soft_button_css(...)`,
  `CHIP_BUTTON_CSS`, `chip_button_css(...)` when applied to a button. Remove
  imports that become unused.

## Roles — choose by what the button DOES, not by its current colour

| role | function | examples |
|---|---|---|
| `primary` | main action: fit, compute, commit | Fit, Fit sizes, Fit all, Grow, Grow N, Calculate 3D, Optimize merge, Save merged data, Apply & save, Export to Igor h5xp, Add to active graph |
| `graph` | put something on the graph | Graph model, Graph unified, Results to graph(s) |
| `tool` | secondary tools and output files | Identify features…, Batch run / Batch all selected, Calculate (contrast), Detect segments, Check similarity, New graph, Save report…, Save to HDF5/NeXus, Export to file…, Import from file…, Find peaks, Import CIF…, Sort by Q, Load selected into slot, Fit cmplx. bckg & sizes |
| `analysis` | uncertainty, Monte Carlo, tabular results | Calc. uncertainty (MC), Compute Monte Carlo I(q), Find best growth, Copy results, Save CSV…, Save as CSV…, Add all to local library, Collect from all selected files |
| `helper` | helper sub-fits, special analyses | Fit background, Fit P/B, Fit flat bckg, Fit power-law bckg, Fix limits / Fix limits?, Ask AI advisor, Energy scan |
| `revert` | undo / reset | Revert back, Reset to defaults, Clear all CIFs |
| `stop` | stop or destroy | Cancel fit / Cancel MC / Cancel selected job, Remove, ×, Auto-reject |
| `settings` | settings to/from file | Store in File, Save/Load params to JSON, Save/Load JSON config…, Save config to JSON…, + Add peak |
| `setup` | load a saved setup | Load setup from file… |
| `help` | help | ? Help |
| `None` (neutral) | utility, navigation, dialogs | Open…, Browse, Select existing folder…, Create default folder, Set vacuum, Parse ▸, All / None, Open output folder, Clear markers, Pop out, dialog Cancel, Load aggregate(s) from NeXus…, Load JSON config is `settings` though |

A button not in the table: pick by function using the table's logic, and
list it in your report. If two buttons with the same function would get
different roles, give them the same one and say so.

## Sizes

`size="large"` when the old style had `font-size: 12px` or larger (big
Fit buttons). `size="small"` when it was ≤ 10px or a chip
(`CHIP_BUTTON_CSS`, ÷10/×10 style helper buttons). Otherwise omit `size`.

## What to drop and what to keep

- Drop from the QSS: background, color, font-weight, border, border-radius,
  padding, hover/pressed/disabled colours. The role rule supplies all of it,
  including a uniform disabled look.
- Keep **sizing** that layout depends on: if the QSS had `min-width`,
  `max-width`, `min-height` or `max-height` and the code does not already set
  it another way, use the widget API (`setMinimumWidth(...)`, etc.).
- Keep any explicit `setFixedWidth` / `setFixedHeight` / `setFont` calls.
- Then the button should carry **no** `setStyleSheet` at all. If a non-colour
  QSS property genuinely must stay (rare), keep a stylesheet with *only* that
  property and say why in the report.

## Do not touch

- Colour-swatch / colour-picker buttons where the background *is* the value
  shown (line colour, marker colour). Leave them as they are.
- Stylesheets on non-buttons (labels, frames, line edits, group boxes).
- Checkable toggle buttons that only style `:checked` without a background
  on the base state — the global rule already sinks them. If one sets a
  base background, give it a role and report it.
- Behaviour, texts, tooltips, layout order, signals.

## Runtime colour changes

Search for places that restyle a button later (state changes, enable/disable
with a colour swap, Fit→Cancel done by hand, "highlight that results are
stale" and so on). Convert each to `set_button_role(btn, ...)` with the
appropriate role. If a button previously changed colour to signal a state that
no role expresses, keep the behaviour by switching roles (for example to
`revert` or `stop`) and report it. Never restore a stylesheet that you have
removed.

## Verify before reporting

Use `/opt/miniconda3/envs/pyirena/bin/python` for everything.

1. `python -m ruff check <your files>` — clean.
2. For each file: grep that no `QPushButton`/`QToolButton` stylesheet with a
   background remains (except colour swatches), and that each
   `set_button_role` role is one of: neutral roles above /
   `primary graph tool analysis helper revert stop settings setup help`.
3. Construct each panel/widget you changed off-screen to prove it still
   builds, e.g.
   `QT_QPA_PLATFORM=offscreen python -c "from pyirena.gui._qt import QApplication; from pyirena.gui.theme import apply_theme; app=QApplication([]); apply_theme(app); from pyirena.gui.sizes_panel import SizesFitPanel; w=SizesFitPanel(); w.show(); app.processEvents(); print('ok')"`
   (adapt the class; if a widget needs constructor arguments, find how its
   caller builds it).
4. `python -m pytest -q -p no:cacheprovider pyirena/tests -k "gui or theme or state or panel"` —
   report the result. A failure in a file you did not touch: report it, do
   not fix it.

## Report (your final message)

Per file: a compact list `label → role [size]`, then
- buttons you mapped by judgement (not in the table),
- runtime restyles you converted,
- anything left as a stylesheet, and why,
- verification results.
Keep it short. Do not paste diffs.
