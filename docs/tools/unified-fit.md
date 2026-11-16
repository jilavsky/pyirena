(unified-fit)=
# Unified Fit

:::{admonition} At a glance
:class: tip

**Launch:** `pyirena-gui` → **Models → Unified Fit**, or the **Unified Fit (GUI)** button
(it opens the first selected file). Standalone: `python -m pyirena.gui.unified_fit`. \
**Python:** `pyirena.api.control` (`open_dataset`, `select_model`, `run_fit`, …) or
`pyirena.batch.fit_unified(data_file, config_file)`. \
**MCP:** `pyirena_ctrl_open_dataset`, then the category `"unified"` tools through
`pyirena_call`; `pyirena_read_unified_fit` reads stored results. \
**Igor equivalent:** {doc}`Irena Unified Fit <igor:Irena/UnifiedFit>`
:::

The Unified Fit GUI provides an interactive interface for fitting small-angle scattering data
using the Unified Fit model (Beaucage model). The GUI closely follows the design of the Igor Pro
Irena package. Each level describes one population of scatterers at a characteristic size Rgᵢ.
If the curve has a Bragg peak or crystalline feature, consider whether it requires a different
model.

(unified-fit-window)=
## The window

```{figure} img/unified-fit/unified-fit-01-window.png
:alt: Unified Fit window after a two-level fit of unified_two_level.h5. Controls on the left, the I versus Q graph with data, fitted model and cursors A and B on the right, residuals below it and the fit-quality message at the bottom. Numbered circles mark the parts described in the legend.

The Unified Fit window after a two-level fit of `validationData/unified/unified_two_level.h5`.
```

The panel uses a split-screen layout: controls on the left, graphs on the right.

1. **Data file**: the loaded file; **Open…** picks another (or drop a file on the panel).
2. **Slit smeared** and slit length **SL**: see [Slit smearing](#unified-fit-slit-smearing).
3. **Q range for fit**: the same two values as cursors A and B; type a Q to move a cursor.
4. **Number of levels** (1–5), **No limits?** and **Identify Features…**.
5. **Level tabs**: G, Rg, B, P and the other parameters of each level.
6. **SAS Background** and **Fit Bckg?**.
7. **Graph Unified**, **Fit** and **Fix limits?**.
8. **Results**: annotate the graph, revert, store, save and load setups, reports.
9. **Passes** and **Calc. Uncertainty (MC)**: Monte Carlo uncertainty.
10. **Graph tabs**: **I vs Q** and **Porod (I·Q⁴)**.
11. **Main plot**: data, model, local fits and the draggable cursors A and B.
12. **Residuals plot**.
13. **Message bar**: the result of the last fit, with the fit-quality metrics.

(unified-fit-workflow)=
## Typical workflow

The steps use `validationData/unified/unified_two_level.h5`, a synthetic two-level curve with
known parameters (Rg ≈ 120 Å and Rg ≈ 1200 Å with the level-2 cut-off linked to level 1, see
`validationData/README.md`).

1. **Load data.** Select a data file in the main GUI and open Unified Fit. Data appears
   immediately in the graph.
2. **Judge the number of levels.** Switch to the **Porod (I·Q⁴)** tab. Flat plateaus separated
   by transitions indicate distinct levels. Count the plateaus to estimate how many levels to
   use.

   ```{figure} img/unified-fit/unified-fit-02-porod-tab.png
   :alt: Porod tab showing I times Q to the fourth versus Q for unified_two_level.h5, with two plateau-and-transition regions.

   The **Porod (I·Q⁴)** tab for `unified_two_level.h5` before any fit.
   ```

3. **Initial parameter estimation.** Start with **1 level** (increase later if residuals show
   structure). For each level, working from largest structure (Level 1) to smallest:

   - **G** — read the low-Q plateau height on the I-Q plot
   - **Rg** — position the cursors in the Guinier region and click **Fit Rg/G btwn cursors**
   - **P** — position the cursors in the power-law region and click **Fit P/B btwn cursors**
   - **B** — use **Estimate B from G/Rg/P?** to get a starting value, or read from Fit P/B
     result

   % CHECK(Jan): "largest structure (Level 1)" — see the phase-5 status report. In this file
   % level 1 is the Rg ≈ 120 Å level and level 2 (RgCO linked to Rg of level 1) the 1200 Å one.

   ```{figure} img/unified-fit/unified-fit-03-guinier-fit.png
   :alt: Level 2 tab after Fit Rg/G between cursors; cursors at Q 0.0002 and 0.0012, green dashed local Guinier curve over the low-Q plateau, G and Rg updated in the level tab.

   **Fit Rg/G btwn cursors** on level 2, cursors at Q = 0.0002–0.0012 Å⁻¹: G ≈ 4.0×10³ cm⁻¹,
   Rg ≈ 1.2×10³ Å.
   ```

4. **Fit.** Check **Fit?** for the parameters to optimise and click **Fit**. Three strategies
   (sequential, cursor-based, full simultaneous) are described under
   [Fitting strategy](#unified-fit-strategy).
5. **Add more levels if needed.** If the residuals show systematic structure (humps, dips):
   increase "Number of levels", switch to the Porod tab to identify where the new level sits in
   Q, set initial parameters for the new level, and fit.
6. **Assess quality.**
   - Residuals should be flat and random around 0
   - χ² (shown in status bar after fit) should be close to 1 for well-scaled data
   - Physical checks: Rg values should be ordered (Level 1 largest), P in 1–6, G positive

   ```{figure} img/unified-fit/unified-fit-04-fit-result.png
   :alt: Graph area after the two-level fit with Results to graphs: model on the data, flat random residuals, parameter boxes for both levels, and the fit-quality message.

   After **Fit** and **Results to graphs**: reduced χ² ≈ 1.06, flat residuals.
   ```

7. **Save results.** **Store in File** writes to the HDF5 source file; **Calc. Uncertainty (MC)**
   runs Monte Carlo uncertainty before saving if needed. See
   [Results and outputs](#unified-fit-results).

% SCREENSHOT unified-fit-06-real-data: a good multi-level fit on real USAXS/SAXS data from
% testData/ (Jan picks the sample), I vs Q tab with "Results to graphs" annotations.

(unified-fit-controls)=
## Controls

### Model and fit controls

| Control | Meaning | Default | Tips |
|---|---|---|---|
| **Graph Unified** | Calculate and plot the model with current parameters (no fitting) | | Preview before fitting |
| **Number of levels** | Select 1–5 structural levels | 1 | Only levels up to this number are active |
| **Update automatically?** | Auto-recalculate whenever a parameter changes (throttled to ~150 ms) | off | Watch both graphs update while scrubbing |
| **Show selected level?** | Overlay individual-level Guinier and power-law contributions on the graph | off | The overlay follows the selected level tab |
| **No limits?** | Disable parameter bounds during fitting | off | Hides the Low / High fields |
| **Identify Features…** | Opens the [Feature Identifier](feature-identifier.md): suggests the number of levels and Guinier Q-windows | | Visualisation only; does not modify the fit |
| **SAS Background** | Constant incoherent background [cm⁻¹] | 1e-06 | |
| **Fit Bckg?** | Include background in fitting | off | |
| **Fit** | Run least-squares optimisation (scipy TRF with bounds, or Nelder-Mead when "No limits?" is checked) | | |
| **Revert back** | Restore parameters to state before last fit | | |
| **Reset to Defaults** | Reset all parameters to defaults | | Cannot be undone; **Revert back** undoes only the last fit |
| **Fix limits?** | Lock current low/high values | | Sets every limit to ±5× the current value |
| **Store local fits** | Persist Guinier/Porod level fits for the "Display local fits" overlay | | |

% CHECK(Jan): "Store local fits" was not found in the current panel; "Show selected level?" is
% listed where the old page had "Display local (Porod & Guinier) fits?".

### Level tabs (1–5)

Each level has its own tab. Only levels up to the selected "Number of levels" are active.

| Control | Meaning | Default | Tips |
|---|---|---|---|
| **G** | Guinier prefactor [cm⁻¹]: low-Q amplitude | 100 | |
| **Rg** | Radius of gyration [Å]: characteristic size | 100 | |
| **Fit Rg/G btwn cursors** | Local Guinier fit to data between the A/B cursors | | |
| **Estimate B from G/Rg/P?** | Auto-calculate B using the Hammouda relationship | off | B is then not fitted |
| **B** | Porod constant [cm⁻¹ Å⁻ᴾ]: high-Q amplitude | 0.01 | |
| **P** | Power-law slope (1–6) | 4 | |
| **Fit P/B btwn cursors** | Local power-law fit to data between the cursors | | |
| **Fit?** | Uncheck to hold parameter fixed during fitting | | One per parameter |
| **Low / High limit** | Bounds used by **Fit** | | Hidden when "No limits?" is checked |
| **RgCutoff** | Beaucage transition radius between levels [Å] | 0 | Levels 2–5; **Link RgCutoff** ties it to Rg of the level below |
| **Correlations?** | Enable Born-Green structure-factor correlations (ETA, PACK) | off | ETA/PACK are only fitted when this is on |
| **Sv [m²/cm³] / Invariant** | Read-only surface-area and scattering-invariant estimates | | |
| **Copy / Move / Swap level** | Reorganize level parameters across tabs | | |

Each parameter value field supports mouse-wheel scrubbing: plain = ±1 %, Shift = ±10 %,
Ctrl = ±0.1 %.

### Graphs

The graph window has two tabs you can switch between at any time. Both update together whenever
the model is recalculated.

**Tab 1 — "I vs Q", main plot (top, ~80 % of height)**

- Log-log scale: Intensity [cm⁻¹] vs Q [Å⁻¹]
- Data points (blue) with error bars
- Unified Fit model curve (red)
- Individual level Guinier (green dashed) and Porod (green dotted) contributions, if "Display
  local fits" is checked
- SAS background level (brown dashed horizontal line), if enabled
- Secondary top axis showing feature radius R = π/Q [Å]
- Two draggable vertical cursors (**A** red, **B** blue) — define the Q range for "Fit btwn
  cursors" and for the global fit

**Tab 1, residuals plot (bottom, ~20 % of height)**

- Log(Q) vs (Data − Fit) / σ
- Dashed horizontal line at 0
- A flat, random residuals trace indicates a good fit

**Graph interaction**

- Scroll wheel: zoom in/out
- Click-drag: pan
- Right-click → *Copy graph to clipboard*, *Save graph as image…* (PNG/JPEG/SVG),
  *Save whole window as image…*, *Save curve data as CSV…*, *Save as Igor Pro ITX…*
  (see "Getting graphs and tables out of pyIrena" on the [Data Selector](data-selector.md) page)
- Zoom is preserved across auto-updates; only resets when new data is loaded

**Cursors**

- Drag cursor A (red) or cursor B (blue) to bracket the Q range of interest
- Cursor positions feed into "Fit Rg/G btwn cursors", "Fit P/B btwn cursors", and the global Fit
  Q range
- Cursors clamp to data boundaries when a new dataset is loaded
- The **Q range for fit** box in the control panel shows the same two values and is editable:
  type a Q and press Return and the cursor moves there. Entries are validated (positive,
  distinct), swapped if reversed, and clamped to the data's Q range. Simple Fits and Modeling
  have the identical control.

**Tab 2 — "Porod (I·Q⁴)"**

- Log-log scale: I·Q⁴ [cm⁻¹·Å⁻⁴] vs Q [Å⁻¹]
- Same data (blue) and model (red) as Tab 1, transformed into Porod presentation
- Local Guinier and Porod level fits also shown (if enabled), scaled by Q⁴
- **No cursors** — this tab is for visual inspection only
- The Porod plot makes it much easier to judge how many levels are needed and where each level's
  Q range lies: a flat region indicates a well-described Porod scatterer; shoulders or bumps
  indicate additional levels
- Right-click → the same export menu as Tab 1 (clipboard, image, curve CSV, ITX)
- Zoom is preserved independently of Tab 1

(unified-fit-results)=
## Results and outputs

| Control | What it does |
|---|---|
| **Store in File** | Write fit result to the HDF5 source file |
| **Export ASCII** | Write Q, I(Q), fit, and residuals as a text file |
| **Copy results** | Put the current fit on the clipboard as Markdown text (levels with uncertainties, fit quality, background, data summary), ready to paste into an e-mail, logbook or manuscript. No file round-trip needed. |
| **Save report…** | The same text as a `.md` file. **Ctrl-click** (⌘-click on macOS) also saves the graph as a PNG beside it and embeds it in the report, so the figure and the numbers travel together |
| **Results to graphs** | Annotate the main graph with level parameters |
| **Passes / Calc. Uncertainty (MC)** | Monte Carlo uncertainty analysis: run N forward-model samples with noise to estimate parameter confidence intervals |
| **Load Setup from File…** | Restore every Unified Fit control from a file stored earlier |
| **Save params to JSON / Load params from JSON** | Save the current parameters to a pyIrena JSON file, and restore them |

% CHECK(Jan): "Export ASCII" was not found in the current panel.

**Store in File** appends the parameters and the model curve to the source NXcanSAS file as
`entry/unified_fit_results`, with the complete panel setup embedded as `_pyirena_config`; see
{ref}`hdf5-results-groups`. To reuse a stored fit:

- **Load Setup from File…** restores every control (parameters, fit flags, bounds, links,
  cursors) from such a file, including files written by the pyirena-ai agent;
- the [Data Selector](data-selector.md) graphs, reports and tabulates stored results across many
  files;
- `pyirena_read_unified_fit(path)` returns them to an AI agent.

(unified-fit-scripting)=
## From Python or AIDA

::::{tab-set}
:::{tab-item} Python
The control API drives the same model as the window. Run from the repository root:

```python
from pyirena.api.control import (
    open_dataset, select_model, set_level_option, set_parameter_value,
    free_parameter, run_fit, save_fit,
)

sid = open_dataset("validationData/unified/unified_two_level.h5")["session_id"]
select_model(sid, nlevels=2)
set_level_option(sid, 2, "link_RGCO", True)    # RgCutoff of level 2 = Rg of level 1
for name, value in {"G_1": 10, "Rg_1": 100, "B_1": 1e-7, "P_1": 4,
                    "G_2": 3000, "Rg_2": 1000, "B_2": 1e-6, "P_2": 3}.items():
    set_parameter_value(sid, name, value)
free_parameter(sid, "P_2")
fit = run_fit(sid)
print(fit["reduced_chi_squared"])               # about 1.07
for p in fit["parameters_updated"]:
    print(p["name"], p["value"])
save_fit(sid, output_path="unified_two_level_fit.h5")   # input file stays untouched
```

For many files with one setup, export the panel's setup with **Save params to JSON** and use
`pyirena.batch.fit_unified(data_file, config_file)`; see [Batch API](../guides/batch-api.md).
:::
:::{tab-item} AIDA / MCP
Ask in plain language, for example *"Fit validationData/unified/unified_two_level.h5 with a
two-level Unified model and link RgCO of level 2 to Rg of level 1."* The agent calls
`pyirena_ctrl_open_dataset`, then `pyirena_call` with `detect_features`, `select_model`,
`set_parameter_value`, `run_fit` and `save_fit`. The recommended call sequence is in the
[AI tools reference](../reference/ai-tools-reference.md#control-tools-reference--unified-fit-dispatcher-category-unified).
:::
::::

(unified-fit-method)=
## Method notes

The model sums contributions from N structural levels:

```
I(q) = Σᵢ [ Gᵢ exp(−q² Rgᵢ²/3)
           + exp(−q² Rgᵢ₋₁²/3) · Bᵢ · ([erf(q Rgᵢ/√6)]³/q)^Pᵢ ]
       + Background
```

Each level describes one population of scatterers at a characteristic size Rgᵢ.

| Parameter | Symbol | Units | Typical Range | Physical meaning |
|-----------|--------|-------|---------------|-----------------|
| Guinier prefactor | G | cm⁻¹ | 1–10¹⁰ | Low-Q scattering amplitude |
| Radius of gyration | Rg | Å | 1–10 000 | Characteristic size |
| Power-law slope | P | — | 1–6 | Surface / mass fractal |
| Porod constant | B | cm⁻¹Å⁻ᴾ | 10⁻²⁰–10¹⁰ | High-Q amplitude |
| Cutoff radius | RgCO | Å | 0–10 000 | Level-to-level transition |

| P | Physical interpretation |
|---|------------------------|
| 1 | Rods (1-D objects) |
| 2 | Platelets (2-D objects) |
| 3 | Rough surface (surface fractal) |
| 4 | Smooth / sharp interface (Porod law) |
| 5–6 | Mass fractals |

(unified-fit-quality)=
### Interpreting fit quality

After a fit, the success message includes robust fit-quality metrics:
- **σ-scale**: how many times actual scatter exceeds reported uncertainties (MAD-based estimate)
- **χ² floor**: the realistic lower bound for reduced χ² given your data's noise level
- **max|(I−M)/I|**: the largest fractional misfit (σ-independent backstop for gross misfits)

**Key insight:** if reduced χ² is high but σ-scale is large (e.g., σ-scale ≈ 3), the *realistic*
floor is `σ-scale²` ≈ 9. A fit at that floor with random scatter and no outliers is **done — you
cannot improve it further** without changing the model or error bars.

See the **[Fit Quality Metrics guide](../methods/fit-quality-metrics.md)** for detailed
interpretation of all diagnostics, decision trees, and common scenarios.

Also check the **residuals plot**: it now shows rescaled residuals `r' = r/σ(robust)`, which
compare scatter against the data's own noise floor rather than (possibly unreliable) reported σ.
Structure in the rescaled plot (long stretches of same-sign) indicates a systematic model
problem; random scatter indicates a mis-scaled σ only.

(unified-fit-slit-smearing)=
### Slit smearing (USAXS)

When slit-smeared USAXS data are loaded (an NXcanSAS `dQl`, e.g. from Matilda), a **"Slit
smeared"** control row appears and the model is smeared to match the data — all levels, the
background, the local Guinier/Porod cursor fits, and the "Show selected level" overlay. Fitted
parameters remain ideal-space (pinhole-equivalent); the invariant is computed from the ideal
model. Note that `ETA`/`PACK` are only fitted when a level's **correlations** are enabled — leave
them unchecked otherwise, or the fit will iterate needlessly (slit-smeared fits amplify this).
See **[Slit smearing](../methods/slit-smearing.md)** for the full reference, scripting contract,
and performance notes.

```{figure} img/unified-fit/unified-fit-05-slit-smeared.png
:alt: Unified Fit window with unified_one_level_slitsmeared.h5 loaded: Slit smeared checked with SL 0.03, one-level fit through the smeared data, flat residuals.

`validationData/unified/unified_one_level_slitsmeared.h5` (dQl = 0.03 Å⁻¹ in the file): slit
smearing switches on by itself and the fitted G, Rg and B are the ideal-space values.
```

(unified-fit-strategy)=
## Fitting strategy

**Sequential** (most reliable):
1. Check "Fit?" only for G and Rg of Level 1 → **Fit**
2. Then also enable B and P → **Fit**
3. Enable background if needed → **Fit**

**Cursor-based** (quick local estimates before global fit):
1. Bracket each region with cursors A and B
2. Click **Fit Rg/G btwn cursors** or **Fit P/B btwn cursors**
3. Then run the global **Fit**

**Full simultaneous fit**:
1. Set reasonable initial values for all parameters
2. Enable all "Fit?" checkboxes
3. Click **Fit** — may need tight bounds if starting far from solution

(unified-fit-tips)=
## Tips and troubleshooting

### Reading the Porod tab

The Porod presentation (I·Q⁴ vs Q) is the fastest way to assess structural complexity:
- Flat regions = well-described Porod scattering at that scale
- Peaks or bumps = additional levels needed
- The slope of the tail in Porod space relates to P − 4

### Mouse-wheel scrubbing

All parameter fields support the mouse wheel with modifier keys:
- **No modifier** — ±1 % of current value per notch
- **Shift** — ±10 % (coarse)
- **Ctrl / Cmd** — ±0.1 % (fine)

Watch both graphs update in real time when "Update automatically?" is checked.

### Cursor placement

- Place cursors to bracket a single structural feature
- For Level 1 Guinier: A at low-Q flat region onset, B just past the Guinier knee
- For Level 1 Porod: A just past the Guinier region, B before the next structural feature

### Common issues

**Fit doesn't converge**
- Improve initial parameter guesses
- Reduce the number of simultaneously fitted parameters
- Tighten the bounds (set realistic low/high limits)

**Unphysical P values (< 1 or > 6)**
- Check if another level is needed in that Q range
- Add a constraint via the limit fields

**Large residuals at specific Q**
- Zoom in to that region on the I-Q tab to inspect the data
- Consider whether a Bragg peak or crystalline feature requires a different model
- An additional level may help

**Fit snaps to wrong region**
- Use cursor-based local fitting first to guide the global fit
- Fix the well-determined parameters and only fit the uncertain ones

## References

- Beaucage, G. "Approximations Leading to a Unified Exponential/Power-Law Approach to Small-Angle
  Scattering" *J. Appl. Cryst.* (1995) **28**, 717-728
- Beaucage, G. "Small-Angle Scattering from Polymeric Mass Fractals of Arbitrary Mass-Fractal
  Dimension" *J. Appl. Cryst.* (1996) **29**, 134-146
- Ilavsky, J. & Jemian, P. R. "Irena: tool suite for modeling and analysis of small-angle
  scattering" *J. Appl. Cryst.* (2009) **42**, 347-353

## Support

For issues or questions: https://github.com/jilavsky/pyirena/issues
