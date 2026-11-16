(simple-fits)=
# Simple Fits

:::{admonition} At a glance
:class: tip

**Launch:** `pyirena-gui` → **Models → Simple Fits**, or select one file in the Data Selector
and click **Simple Fits (GUI)**; the panel opens with that file's data loaded and plotted. \
**Python:** `pyirena.batch.fit_simple(data_file, config)`, or `pyirena.api.control`
(`select_simple_model`, `run_simple_fit`, …). \
**MCP:** `pyirena_ctrl_open_dataset`, then the category `"simple"` tools through `pyirena_call`;
`pyirena_read_simple_fit` reads stored results. \
**Igor equivalent:** {doc}`Irena Simple Fits <igor:Irena/SimpleFits>` and
{doc}`System specific models <igor:Irena/SysSpecModels>`
:::

The Simple Fits tool provides direct analytical model fitting for common small-angle
scattering functions.  It is ported from the Igor Pro `IR3_SimpleFits.ipf` and
`IR3_SystemSpecificModels.ipf` modules of the Irena package.

Simple Fits fits a single analytical scattering model to I(Q) data over a
user-selected Q range.  It is most useful when the scattering is well described
by one function — for example, a Guinier region at low Q, a Porod tail at high Q,
or scattering from monodisperse spherical particles.

For complex multi-level structures use the **Unified Fit** tool instead.
For model-independent particle size distributions use the **Size Distribution** tool.

(simple-fits-window)=
## The window

```{figure} img/simple-fits/simple-fits-01-window.png
:alt: Simple Fits window after a Guinier fit of simple_guinier.h5. Controls on the left with the Guinier model, Q range, parameter grid with I0 and Rg, and results; on the right the I versus Q graph with data, model and cursors A and B, the residuals, and the Guinier linearization ln I versus Q squared with its straight-line fit. Numbered circles mark the parts described in the legend.

The Simple Fits window after a Guinier fit of `validationData/simple_fits/simple_guinier.h5`.
```

1. **Data file**: the loaded file; **Open…** picks another (or drop a file on the panel).
2. **Slit smeared** and slit length **SL**: see [Slit smearing](#simple-fits-slit-smearing).
3. **Model**: the analytical model (see [Models](#simple-fits-models)).
4. **Q range for fit**: Q min and Q max, the same values as cursors A and B.
5. **No limits?** and **Complex background (B·Q⁻ᴾ + flat)**.
6. **Parameters**: one row per parameter — **Fit?**, value, lo, hi and ± std.
7. **Graph model** and **Fit** (**Calculate Invariant** for the Invariant).
8. **Fit results**: χ² and reduced χ².
9. **Result buttons**: annotate the graph, revert, store, save and load setups, reports, reset.
10. **Passes** and **Calc. Uncertainty (MC)**.
11. **I(Q) plot**: data, model and the cursors.
12. **Residuals plot** (rescaled, see [Fit quality](#simple-fits-fit-quality)).
13. **Linearization plot** (Guinier and Porod families, see [Linearization plots](#simple-fits-linearization)).
14. **Status line**: the result of the last fit, with the fit-quality metrics.

(simple-fits-workflow)=
## Typical workflow

The steps use `validationData/simple_fits/simple_guinier.h5`, a Guinier curve with I0 = 250 cm⁻¹
and Rg = 45 Å (see `validationData/README.md`).

1. Select a **model** from the dropdown: Guinier.
2. Set the **Q range** (see [Q range selection](#simple-fits-q-range)). For a Guinier fit keep
   `Q·Rg < 1.3`.
3. Enter **starting parameter values** and optionally adjust lower/upper **limits**.
   - Uncheck "Fit?" to hold a parameter fixed at its current value.
   - Check "No limits" to ignore all bounds during fitting (uses very wide defaults).
4. Click **Fit**.  The model curve appears on the I(Q) plot; residuals update; the
   linearization panel updates (if applicable).  Each parameter field shows `± σ`
   (uncertainty from the covariance matrix). Here Rg ≈ 45.8 ± 1.1 Å (see the window above).
5. If convergence is poor:
   - Try better starting values (the linearization plot can guide initial Rg/Rc estimates).
   - Narrow the Q range to exclude non-model contributions.
   - Try "No limits" if the optimiser is hitting a bound unexpectedly.
6. **Estimate uncertainties** with **Calc. Uncertainty (MC)**, then **Store in File**.

The same steps fit a monodisperse sphere. With
`validationData/simple_fits/simple_sphere.h5` (Scale = 8 cm⁻¹, R = 120 Å, flat background
0.01 cm⁻¹) choose **Sphere**, check **Complex background**, hold B = 0 and P fixed and fit
`flat`:

```{figure} img/simple-fits/simple-fits-02-sphere.png
:alt: Graph area after a Sphere fit of simple_sphere.h5 with a flat background: the model follows three form-factor oscillations, residuals are flat; there is no linearization for this model.

**Sphere** with a flat complex background: Scale ≈ 8.0 cm⁻¹, R ≈ 120 Å.
```

% SCREENSHOT simple-fits-04-real-data: a Guinier or Porod fit on real data from testData/
% (Jan picks), with the linearization plot showing the fitted range highlighted.

(simple-fits-controls)=
## Controls

| Control | Meaning | Default | Tips |
|---|---|---|---|
| **Model** | The analytical model to fit (see [Models](#simple-fits-models)) | Guinier | |
| **Q min / Q max** | Q range for fit, the same as cursors A and B (see [Q range selection](#simple-fits-q-range)) | | |
| **No limits?** | Ignore all parameter bounds and fit unconstrained | off | The lo/hi fields are hidden |
| **Complex background (B·Q⁻ᴾ + flat)** | Add a background term `B·Q⁻ᴾ + flat` (see [Background options](#simple-fits-background)) | off | Not available for Porod and Power Law |
| **Extend by Porod tail (Kp/Qmax)** | Invariant only: add the analytic high-Q tail | off | |
| **Refit background from saved ranges** | Invariant only: re-run the background prefits before every Calculate | off | Also applies in scripted/batch runs |
| **Fit?** / value / **lo** / **hi** | Per parameter: fit or hold, starting value, bounds | registry defaults | Uncheck "Fit?" to hold a parameter fixed |
| **± std** | Uncertainty from the covariance matrix, or from Monte Carlo | | |
| **Fit B/P btwn cursors** | Fit the background power-law B·Q⁻ᴾ between the cursors | | With complex background |
| **Fit Flat btwn cursors** | Estimate the flat term between the cursors | | With complex background |
| **Graph model** | Compute and display the current model curve without fitting | | |
| **Fit** | Fit the selected model on the current Q range | | **Calculate Invariant** for the Invariant |
| **Revert back** | Restore the parameter values from before the last fit | | |
| **Reset to defaults** | Reset all parameters for the current model to their registry defaults | | |
| **Passes / Calc. Uncertainty (MC)** | Repeat the fit on noise-perturbed data | 10 | See [Monte Carlo uncertainty](#simple-fits-monte-carlo) |

(simple-fits-q-range)=
### Q range selection

The **Q range for fit** box and the two graph cursors are two views of the same
thing — change either and the other follows.

**Manual entry:** type a value into **Q min** or **Q max** and press Return (or
click away).  The cursors move to match.  Values are checked as you commit
them: both must be positive and different, a pair entered the wrong way round
is swapped for you, and values outside the loaded data are clamped to its Q
range.  A rejected entry reverts and the status line says why.

**Cursor-driven:** drag the two vertical cursor lines (**A** = Q min, **B** =
Q max) on the I(Q) plot.  The fields update as you drag.

Only data points inside the Q range are used for fitting.  The model is evaluated
and plotted over the full data Q range so you can see extrapolation quality.

The same **Q range for fit** control appears in Modeling and Unified Fit and
behaves identically there.

(simple-fits-background)=
### Background options

**Models with explicit Background (Porod, Power Law).** Background is a regular fitted
parameter — enter a starting value and bounds like any other parameter.

**Models with complex background (all others).** Check **"Complex background"** to add a
background term of the form:

```
BG(Q) = B · Q^(−P) + flat
```

Three additional parameters appear in the grid, labelled **B**, **P**, **flat**
(stored internally as `BG_B`, `BG_P`, `BG_flat`):

| Parameter | Meaning | Typical use |
|-----------|---------|-------------|
| `B` | Power-law prefactor | Parasitic scattering |
| `P` | Power-law exponent | Usually 4 (Porod) |
| `flat` | Flat (incoherent) background | Incoherent scattering |

The naming mirrors the Unified Fit convention (B for prefactor, P for exponent).

To add only a flat background: set `B = 0` (fixed) and fit `flat`.
To model a Porod-law background: fix `P = 4` and fit `B`.

**Prefit helpers.**  When complex background is active, two small buttons appear
below the parameter grid:

- **Fit B/P btwn cursors** — fits the power-law `B·Q⁻ᴾ` to the data between the
  two cursors and writes the result as starting values.  If P's "Fit?" box is
  unchecked, only B is fit at the current (model-guided) P.
- **Fit Flat btwn cursors** — estimates the flat term from the median residual
  over the cursor window (using the current B and P).

Place the cursors over the low-Q power-law tail for **Fit B/P**, and over a
high-Q region where the data has flattened for **Fit Flat**.

(simple-fits-results)=
## Results and outputs

(simple-fits-report)=
### Results as text — Copy results / Save report…

**Copy results** puts the current fit on the clipboard as Markdown: the model,
χ² and reduced χ², degrees of freedom, the fitted Q range, every parameter with
its 1σ uncertainty, the derived quantities, and the data summary.  Paste it
straight into an e-mail, a logbook or a manuscript — no HDF5 round-trip.
**Save report…** writes the same text to a `.md` file.  **Ctrl-click**
(⌘-click on macOS) also captures the graph as a PNG next to the report and
embeds it as a figure — one action instead of save-report, save-graph, then
inserting the image by hand.

The text is produced by the same builder as the Data Selector's *Create Report*
and the AI agents' `export_fit_report`, so all three agree to the digit.

(simple-fits-store)=
### Store in File (HDF5 only)

Click **"Store in File"** to write results into `entry/simple_fit_results` in the
NXcanSAS HDF5 file, with the full panel setup embedded so **Load Setup from File…** can restore
every control (see {ref}`hdf5-results-groups`).  The group contains:

```
entry/simple_fit_results/
    attrs:  model, success, chi_squared, reduced_chi_squared, dof,
            q_min, q_max, use_complex_bg, timestamp, program
    Q             [array, Å⁻¹]
    I_model       [array]
    residuals     [array]
    intensity_data  [array]   (copy of fitted I data)
    intensity_error [array]
    params/         one scalar per parameter
    params_std/     one scalar per parameter (MC or covariance uncertainty)
    derived/        model-specific computed values (xi, d, thickness, …)
```

Stored results load back with `load_simple_fit_results` (see the Python tab under
[From Python or AIDA](#simple-fits-scripting)); `pyirena_read_simple_fit(path)` returns them to
an AI agent.

% CHECK(Jan): the old page also showed `from pyirena import load_result;
% load_result("sample.h5", "simple_fits")` as a top-level convenience function. It raises
% ValueError: load_result supports only 'unified_fit' and 'size_distribution'. Removed from the
% runnable example; restore once load_result supports Simple Fits.

(simple-fits-export-import)=
### Export / Import Parameters

**Export Parameters** writes the current model configuration — model name, parameter
values, limits, Q range, complex-background flag — into a `pyirena_config.json` file
under the `simple_fits` key.  Multiple tool configurations can coexist in the same
JSON file.

**Import Parameters** reads the `simple_fits` section back and restores all values.

These files drive batch processing (see [From Python or AIDA](#simple-fits-scripting)).

% CHECK(Jan): the buttons are labelled "Save params to JSON" and "Load params from JSON" in the
% current window.

(simple-fits-data-selector)=
### Data Selector integration

With **"Simple Fits"** checked in the Data Selector:

| Button | Behaviour |
|--------|-----------|
| **Create Graph** | Opens a two-panel window (I(Q) + residuals) showing stored fit results from all selected HDF5 files.  Each file gets a unique colour; the model line uses a darker shade for contrast. |
| **Create Report** | Adds a "## Simple Fits" section to each file's Markdown report: model, chi², reduced chi², DOF, Q range, all parameters ± std, derived quantities. |
| **Tabulate Results** | Adds `SF_model`, `SF_chi2`, `SF_reduced_chi2`, `SF_dof`, `SF_q_min`, `SF_q_max`, `SF_use_complex_bg`, `SF_<param>`, `SF_<param>_std`, and `SF_derived_<name>` columns to the CSV table. |
| **Simple Fits (script)** | Reads `pyirena_config.json`, fits all selected files, saves results. |

To process many files without the GUI, use **"Simple Fits (script)"** — this reads
the `simple_fits` section from `pyirena_config.json` and fits all selected files,
saving results to each HDF5 file automatically.

(simple-fits-scripting)=
## From Python or AIDA

::::{tab-set}
:::{tab-item} Python
Run from the repository root. `fit_simple` stores its results into the data file by default,
so the example works on copies:

```python
import shutil

from pyirena import MODEL_NAMES
from pyirena.batch import fit_simple
from pyirena.io.nxcansas_simple_fits import load_simple_fit_results

# List all available model names
print(MODEL_NAMES)

# Work on copies: fit_simple stores its results into the data file by default
shutil.copy("validationData/simple_fits/simple_guinier.h5", "guinier.h5")
shutil.copy("validationData/simple_fits/simple_sphere.h5", "sphere.h5")
shutil.copy("validationData/sizes/sizes_sphere_flat_background.h5", "invariant.h5")

# Fit one file — minimal config
result = fit_simple("guinier.h5",
                    config={'model': 'Guinier'})   # uses registry defaults

# Fit with explicit starting values and bounds
result = fit_simple("guinier.h5", config={
    'model':  'Guinier',
    'params': {'I0': 100.0, 'Rg': 30.0},
    'limits': {'I0': [1e-3, 1e6], 'Rg': [1.0, 500.0]},
}, q_min=0.002, q_max=0.03)

if result and result['success']:
    print(f"Rg = {result['params']['Rg']:.2f} ± "
          f"{result['params_std']['Rg']:.2f} Å")          # true Rg = 45 Å
    print(f"Reduced chi² = {result['reduced_chi2']:.3f}")

# With Monte Carlo uncertainty; the flat background goes in the complex background
result = fit_simple("sphere.h5",
                    config={'model': 'Sphere',
                            'use_complex_bg': True,
                            'params': {'Scale': 5.0, 'R': 100.0,
                                       'BG_B': 0.0, 'BG_P': 4.0, 'BG_flat': 0.01}},
                    fixed_params={'BG_B': 0.0, 'BG_P': 4.0},
                    with_uncertainty=True, n_mc_runs=100)
if result and result['success']:
    print(result['params'], result['params_std'])        # true Scale 8, R 120 Å

# Invariant (calculation — no fitting): volume fraction from Q*
result = fit_simple("invariant.h5", config={
    'model': 'Invariant',
    'use_complex_bg': True,                  # subtract B·Q⁻ᴾ + flat first
    'params': {'Contrast': 100.0,            # Δρ² in 10²⁰ cm⁻⁴
               'BG_B': 0.0, 'BG_P': 4.0, 'BG_flat': 0.02},
    'invariant_porod_tail': False,           # True → add Kp/Qmax tail
    # Recommended for batch: refit the background per file from saved
    # Q windows (recorded by the GUI prefit buttons) before calculating.
    'bg_prefit': {
        'enabled': True,
        'power_law': {'use': False},
        'flat':      {'use': True, 'q_min': 0.2, 'q_max': 0.35},  # high-Q
    },
}, fixed_params={'BG_B': 0.0}, q_min=0.0008, q_max=0.35)
if result and result['success']:
    d = result['derived']
    print(f"Q* = {d['Invariant']:.4g} cm⁻⁴,  φ = {d['VolumeFraction']:.4g},"
          f"  QmaxUsed = {d['QmaxUsed']:.4g} Å⁻¹")       # true φ = 0.006

# Loading stored results (fit_simple stored the Guinier fit in guinier.h5 above)
stored = load_simple_fit_results("guinier.h5")
print(stored['model'], stored['chi_squared'])
print(stored['params'])        # dict of parameter values
print(stored['params_std'])    # dict of uncertainties
print(stored['derived'])       # model-specific computed values
```

For many files, export the panel's setup with **Save params to JSON** and pass its
`simple_fits` section as `config`, or run every configured tool at once with
`pyirena.batch.fit_pyirena(data_file, config_file)`. For the full batch API reference, see
[Batch fitting API](../guides/batch-api.md).
:::
:::{tab-item} AIDA / MCP
Ask in plain language, for example *"Do a Guinier fit of
validationData/simple_fits/simple_guinier.h5 and give me Rg with its uncertainty."* The agent
calls `pyirena_ctrl_open_dataset`, then `pyirena_call` with `select_simple_model`,
`set_simple_parameter`, `run_simple_fit` and `save_simple_fit`. See the
[AI tools reference](../reference/ai-tools-reference.md#control-tools-reference--simple-fits-dispatcher-category-simple).
:::
::::

% CHECK(Jan): the old Scripting API block used placeholder files (sample.h5, data/,
% pyirena_config.json). It now runs on validationData copies; the Sphere example adds the
% flat complex background the dataset needs, and the Invariant example uses
% sizes_sphere_flat_background.h5 (true φ = 0.006) with a flat-only background prefit. The
% loop over data/*.h5 and the fit_pyirena call, which need a config file, are described in
% the sentence below the code instead.

(simple-fits-method)=
## Method notes

(simple-fits-models)=
### Models

#### Guinier family

| Model | Formula | Linearization |
|-------|---------|--------------|
| **Guinier** | `I(Q) = I0·exp(−Q²Rg²/3)` | ln I vs Q² → slope = −Rg²/3 |
| **Guinier Rod** | `I(Q) = I0·exp(−Q²Rc²/2)/Q` | ln(QI) vs Q² → slope = −Rc²/2 |
| **Guinier Sheet** | `I(Q) = I0·exp(−Q²Rg²)/Q²` | ln(Q²I) vs Q² → slope = −Rg² |

Guinier Sheet also derives thickness `t = 2√3·Rg`.

Valid Q range for Guinier approximation: `Q·Rg < 1.3`.

#### Porod / Power Law

| Model | Formula | Parameters |
|-------|---------|-----------|
| **Porod** | `I(Q) = Kp·Q⁻⁴ + Background` | Kp, Background |
| **Power Law** | `I(Q) = P·Q⁻ⁿ + Background` | P (prefactor), Exponent n, Background |

These models have an explicit flat Background parameter.  The complex background
option is not available for these two models.

#### Debye Polymer Chain

| Model | Formula |
|-------|---------|
| **Debye Polymer Chain** | `I(Q) = Scale · 2(exp(−x) − 1 + x) / x²`  where `x = Q²·Rg²` |

The Debye function (P. Debye, 1947) describes single-chain scattering from an ideal
Gaussian polymer.  It is exact for a random-walk chain with no excluded volume.
Key properties:

- `P(Q) → 1` as `Q → 0` (forward-scattering limit)
- At high Q, `P(Q) ∝ Q⁻²` — characteristic of a chain with fractal dimension 2

Parameters: **Scale** (forward-scattering intensity at Q=0), **Rg** (radius of
gyration of the chain).  Complex background is supported.

#### Sphere and Spheroid

| Model | Formula |
|-------|---------|
| **Sphere** | `I(Q) = Scale·|FF(QR)|²`;  `FF(x) = 3(sin x − x cos x)/x³` |
| **Spheroid** | orientational average of sphere FF with semi-axis ratio β (oblate β<1, prolate β>1) |

The spheroid orientational average is computed with Gauss-Legendre quadrature
(50-point rule) for efficiency.

Both models support complex background.  A flat background can be added via the
complex background option (set B = 0, flat = desired value).

#### Correlation models and Invariant

| Model | Formula / description |
|-------|-----------------------|
| **Debye-Bueche** | `I(Q) = Prefactor·Eta²·ξ³/(1+Q²ξ²)²` — two-phase random media |
| **Teubner-Strey** | `I(Q) = Prefactor/(A + C1·Q² + C2·Q⁴)` — micro-emulsion / lamellar |
| **Benedetti-Ciccariello** | Correlated two-phase model with interfacial layer; derives from SLD contrast |
| **Hermans** | Paracrystalline lamellar model; complex exponentials for d1, d2, σ1, σ2 |
| **Hybrid Hermans** | Hermans + one Guinier term + one Unified Fit level |
| **Unified Born Green** | Two-level Unified model with Born-Green structure factor S(Q); gives correlation length ξ and inter-particle distance |

Teubner-Strey derives correlation length ξ = (A/C2)^¼/√2 and repeat distance d.
These appear in the "Derived" section of the results.

(simple-fits-linearization)=
### Linearization plots

| Model | X axis | Y axis | Slope encodes |
|-------|--------|--------|--------------|
| Guinier | Q² | ln I | −Rg²/3 |
| Guinier Rod | Q² | ln(QI) | −Rc²/2 |
| Guinier Sheet | Q² | ln(Q²I) | −Rg² |
| Porod | Q⁴ | IQ⁴ | — (extrapolates to Kp at Q=0) |

The fitted Q-range points are highlighted; points outside the range are shown in grey.
A linear fit line is overlaid with the intercept labelled.

For all other models the linearization panel shows "No linearization available for
\<model\>".

% CHECK(Jan): in the Sphere shot the linearization panel is hidden rather than showing the
% "No linearization available" text.

(simple-fits-monte-carlo)=
### Monte Carlo uncertainty

The covariance-matrix uncertainties from `scipy.optimize.curve_fit` can underestimate
real parameter uncertainty when the model is non-linear or the data have correlated
errors.  The Monte Carlo option provides an independent estimate:

1. Set **N MC runs** (default 50; increase for publication-quality estimates).
2. Click **Calculate Uncertainty**.  The tool runs N independent fits on
   Gaussian-perturbed copies of the data (perturbation amplitude = dI).
3. The standard deviation of each parameter across the N fits updates the ± labels.

MC uncertainties reflect data noise propagation more faithfully than the covariance
matrix for models with strong parameter correlations (e.g., Spheroid Scale and R).

% CHECK(Jan): the window's controls are "Passes" (default 10) and "Calc. Uncertainty (MC)";
% 50 is the default of fit_simple(n_mc_runs=...) in scripts.

(simple-fits-fit-quality)=
### Fit quality

Residuals are displayed **rescaled** as r' = r / σ(robust): the normalized
residual divided by a robust (MAD-based) estimate of the actual noise scale, so
the scatter is judged against the data's own noise floor rather than the (often
mis-scaled) reported σ. The fit status line also reports **σ-scale** (how many ×
the actual scatter exceeds the reported σ), the realistic reduced-χ² floor, and
the largest fractional misfit **max|(I−M)/I|** — a σ-independent gross-misfit
backstop.

See the **[Fit Quality Metrics guide](../methods/fit-quality-metrics.md)** for full
interpretation of these numbers.

(simple-fits-slit-smearing)=
### Slit smearing (USAXS)

When slit-smeared USAXS data are loaded (an NXcanSAS `dQl`, e.g. from Matilda),
a **"Slit smeared"** control row appears and each analytic model is smeared
before comparison, so the displayed "Graph model" curve matches the data
(e.g. Sphere/Spheroid oscillations are correctly smoothed) and fitted parameters
remain ideal-space. Guinier/Porod **linearization** has no closed form for
smeared data, so it is drawn best-effort (ideal-space) and labelled as such. The
**Invariant** integrates the data directly, which is invalid for smeared data,
so it is **disabled with a message** — load the desmeared dataset, or use Unified
Fit's model-based invariant. See **[Slit smearing](../methods/slit-smearing.md)** for the
full reference.

## References

- Debye, P. (1947), the Debye function for an ideal Gaussian polymer chain (Debye Polymer
  Chain model).
- Irena: Igor Pro `IR3_SimpleFits.ipf` and `IR3_SystemSpecificModels.ipf`; see
  {doc}`Irena Simple Fits <igor:Irena/SimpleFits>`.

% CHECK(Jan): the old page had no References section; the two entries above come from the
% text. Add the papers for Teubner-Strey, Debye-Bueche, Benedetti-Ciccariello and Hermans.

(simple-fits-invariant)=
## Invariant (calculation — no fitting)

The **Invariant** entry in the model selector is not a fit: it directly
integrates the data to obtain the Porod invariant

```
Q* = ∫ q²·I(q) dq  =  2π² · Δρ² · φ(1−φ)
```

which is exact for any two-phase system with sharp interfaces, independent of
morphology.  This makes it a very robust way to obtain volume fractions for
polydisperse single-"particle"-type systems (porosity in rocks, precipitates,
voids) — no model, no fitting, just calibrated data and a contrast value.

Workflow:

1. Load data and select **Invariant** as the model.
2. If a background must be removed, check **Complex background** and use the
   **Fit B/P btwn cursors** / **Fit Flat btwn cursors** prefit buttons on a
   background-dominated Q region — the background `B·Q⁻ᴾ + flat` is
   subtracted from the data before integration.  The Q window used by each
   prefit button is **remembered** (shown below the Invariant options and
   stored with the setup).  Check **"Refit background from saved ranges"**
   to re-run these prefits automatically before every Calculate — in the
   GUI *and* in scripted/batch runs.  This is essential for scripting:
   the invariant is only as good as the background subtraction, so batch
   runs re-determine B/P/flat per file from the saved windows instead of
   trusting exported values.  The background parameters' **"Fit?"
   checkboxes gate the replay**: unchecking `BG_flat` leaves the flat term
   at its entered value, unchecking `BG_P` refits B at the held P, and
   unchecking `BG_B` skips the power-law refit entirely (use this — with
   `BG_B` = 0 — for datasets in a sequence that should not carry a
   low-Q power-law term while the flat background is still re-determined).
   The same gating applies in scripted/batch runs via `param_fixed`.
3. Enter **Contrast** Δρ² in units of 10²⁰ cm⁻⁴ (default 100; use the
   Scattering Contrast calculator if needed).  The Invariant is a direct
   integration, not a least-squares fit, so Contrast is used exactly as
   entered — it has no **Fit?** box and no lo/hi limit fields.  (The
   background terms `B`, `P` and flat *are* determined by the prefit and keep
   theirs.)
4. Position the cursors over the Q range to integrate and press
   **Calculate Invariant**.

```{figure} img/simple-fits/simple-fits-03-invariant.png
:alt: Simple Fits window with the Invariant model on sizes_sphere_flat_background.h5: complex background with the flat term from Fit Flat between cursors, background-corrected data in green, the flat background line, and the running integral of q squared I rising to a plateau on the right axis; the status line reports Q star, phi and QmaxUsed.

**Calculate Invariant** on `validationData/sizes/sizes_sphere_flat_background.h5` after **Fit
Flat btwn cursors** at Q = 0.2–0.35 Å⁻¹: φ = 0.0058 (generated: 0.006), QmaxUsed = 0.265 Å⁻¹.
```

The graph shows the background-corrected data (green), the background curve,
and the running integral `∫q²I dq` on the **right axis** (red) — analogous to
the integrated-volume-fraction curve in the Size Distribution tool.  The
invariant is read at the point where the running integral flattens (smoothed
derivative crosses zero); this Q value is reported as **QmaxUsed** and guards
against a slightly over-subtracted background at high Q.

Details and options:

- **Extrapolation to Q = 0** is always applied (the first integrand value is
  repeated at Q = 0, following Igor Irena).  Its contribution is ~I₁q₁³ —
  negligible for any reasonable q_min — so there is no user control for it.
- **Extend by Porod tail (Kp/Qmax)** (optional, off by default = Igor
  behaviour): adds the analytic high-Q tail `∫_{Qmax}^{∞} Kp·q⁻² dq = Kp/Qmax`
  to compensate for measurement truncation.  Kp is estimated as the median of
  `(I−bg)·q⁴` over the last half-decade before QmaxUsed.
- **Volume fraction**: with data on absolute scale [cm⁻¹] and correct
  contrast, φ is obtained by solving `φ(1−φ) = Q*/(2π²Δρ²)` (lower root,
  φ < 0.5 — the invariant cannot distinguish φ from 1−φ).  If
  `φ(1−φ) > 0.249` the value is unphysical: φ is set to NaN and a warning
  points at calibration/contrast problems.

Results stored in `derived/`: `Invariant` [cm⁻⁴], `VolumeFraction`,
`PhiOneMinusPhi`, `QmaxUsed` [Å⁻¹], and — when the tail option is on —
`PorodTail` [cm⁻⁴] and `PorodKp` [cm⁻¹Å⁻⁴].  The running integral is stored
as `Q_integral` / `running_integral` datasets and is exported to Igor
experiments as the `SimFitInvariantIntegral` wave.

**Caveats**: requires absolute intensity calibration and desmeared data
(slit-smeared USAXS must be desmeared first); q² amplifies background errors
at high Q, so background subtraction matters; systems with more than two
phases or diffuse interfaces bias the result.  χ², residuals, linearization
and Monte Carlo uncertainty do not apply and are disabled for this model.
