(waxs-peak-fit)=
# WAXS Peak Fit

:::{admonition} At a glance
:class: tip

**Launch:** `pyirena-gui` → **Models → WAXS Peak Fit**, or the **WAXS Peaks (GUI)** button (it
opens the first selected file). \
**Python:** `pyirena.api.control` (`select_waxs_model`, `find_waxs_peaks`, `run_waxs_fit`, …) or
`pyirena.batch.fit_waxs(data_file, config_file)` / `fit_waxs_peaks(data_file, config)`. \
**MCP:** `pyirena_ctrl_open_dataset`, then the category `"waxs"` tools through `pyirena_call`;
`pyirena_read_waxs_peakfit` reads stored results. \
**Igor equivalent:** {doc}`Irena WAXS <igor:Irena/WAXS>` (powder diffraction fitting)
:::

The **WAXS Peak Fit** tool fits diffraction peaks (Gaussian, Lorentzian, Pseudo-Voigt, or
log-normal profiles) rising above a smooth background in Wide-Angle X-ray Scattering data. All
plots are **linear/linear** (intensity vs. Q in Å⁻¹). The **Diffraction Lines** tab overlays
theoretical stick patterns from CIF files for phase identification; it is not a Rietveld
refinement.

% CHECK(Jan): the peak-shape list (Gaussian, Lorentzian, Pseudo-Voigt, log-normal) omits "Voigt",
% which the Shape combo and core PEAK_SHAPES offer.

(waxs-peak-fit-window)=
## The window

```{figure} img/waxs-peak-fit/waxs-peak-fit-01-window.png
:alt: WAXS Peak Fit window after fitting waxs_three_gauss.h5 with a linear background and three Gaussian peaks. Controls on the left, the linear I versus Q graph with data, total model, each peak and the background on the right, residuals below. Numbered circles mark the parts described in the legend.

The WAXS Peak Fit window after **Find Peaks**, removing three false detections and **Fit** on
`validationData/waxs/waxs_three_gauss.h5`.
```

The window is split into a **left control panel** (≈430 px wide) and a **right graph area**. The
left panel contains two tabs: **WAXS Peak Fit** — peak finding, background, and fitting controls;
and **Diffraction Lines** — theoretical stick-pattern overlay for phase identification (see
[Diffraction Lines tab](#diffraction-lines-tab)).

1. **Data file**: the loaded file; **Open…** picks another.
2. **Fit Q range** (set by the two cursors) and **No limits?**.
3. **Background** shape and its parameters.
4. **Peak Finding** parameters and **Find Peaks**.
5. **Peaks**: one row per peak, with **Remove**; **Sort by Q** above.
6. **Graph Model**, **Fit** and **Fix Limits**.
7. **Weighting**.
8. **Q₀ pre-search**.
9. **Results** buttons and the status line with the fit quality.
10. **Main graph**: data scatter, total model (red), per-peak overlays (colors), background
    (dashed gray), diffraction line sticks (one color per phase), Q-range cursors.
11. **Residuals**: rescaled residuals r' = r/σ(robust); dashed zero line; x-axis linked to main.

(waxs-peak-fit-workflow)=
## Typical workflow

The steps use `validationData/waxs/waxs_three_gauss.h5`: three Gaussian peaks (Q₀ = 1.85, 2.55,
3.05 Å⁻¹) on a linear background, with known parameters (see `validationData/README.md`).

1. Open **WAXS Peaks (GUI)** from the Data Selector.
2. Drag the Q cursors to bracket the diffraction peaks of interest.
3. Choose a **Background** shape.
   - For XRD / powder data try **SNIP** first.
   - For data with a gently-varying background, start with **Cubic** or **Linear**.
4. Click **Find Peaks** — the tool detects peaks automatically.

   ```{figure} img/waxs-peak-fit/waxs-peak-fit-02-find-peaks.png
   :alt: WAXS Peak Fit window right after Find Peaks on waxs_three_gauss.h5 with a Linear background: six detected peaks, three of them tiny false detections, drawn over the data.

   After **Find Peaks** with the default settings: the three real peaks plus three small false
   detections (2, 3 and 5), to be removed in the next step.
   ```

5. Inspect the overlay: use mouse-wheel on **Q0** and **FWHM** to fine-tune initial positions.
   Right-click → **Add Peak** for any missed peaks; use **Remove** for false detections.
6. Select a **Weighting** mode. If peaks have few points compared to the background, try
   **Equal** or **Relative**.
7. If peaks are narrow or their positions may have drifted, enable one or both **Q₀ pre-search**
   options. For a data series where peaks shift coherently, tick **Cross-corr.** and set the
   **Search window** to the expected maximum drift. For very narrow peaks that may have shifted
   independently, also tick **Per-peak Q₀ scan**.
8. Click **Graph Model** to preview the current parameter set without fitting.
9. Click **Fit** — pre-search (if enabled) runs first, then `curve_fit`. Fitted values replace
   the initial guesses.

   ```{figure} img/waxs-peak-fit/waxs-peak-fit-03-fit-result.png
   :alt: Graph area after Fit and Results to graphs: three Gaussian peaks and a linear background through the data, flat residuals, and an annotation listing bg0, bg1 and Q0, FWHM, A of the three peaks.

   After **Fit** and **Results to graphs**: Q₀ = 1.850, 2.550, 3.049 Å⁻¹ and FWHM = 0.090, 0.120,
   0.151 Å⁻¹, as generated.
   ```

10. If the fit diverged, click **Revert**, adjust limits or initial values, and refit.
11. Click **Store in File** to save results.
12. Click **Export Parameters** to write `pyirena_config.json` for batch use (pre-search settings
    are included).

% CHECK(Jan): steps 10 and 12 name "Revert" and "Export Parameters"; the buttons are now
% "Revert back" and "Save params to JSON".

(waxs-peak-fit-controls)=
## Controls

(waxs-peak-fit-graph)=
### Graph and cursors

**Right-click** on the main graph to:
- **Add Peak at Q = X.XXX Å⁻¹** — inserts a peak row initialized at that Q
- **Show data as line / Show data as points** — toggle scatter vs. connected line
- **Hide / Show error bars** — toggle ±σ vertical bars on data
- **Copy graph to clipboard**, **Save graph as image…** (PNG/JPEG/SVG at 1600 px), **Save whole
  window as image…**, **Save curve data as CSV…**, **Save as Igor Pro ITX…** — the standard
  export menu, described on the [Data Selector](data-selector.md) page

Two vertical cursors (cyan dashed lines) define the **Q fit range**. The read-only Q-min / Q-max
display boxes below the title bar update as you drag the cursors. Only data inside this range is
used for fitting; residuals outside are shown as NaN (blank).

% CHECK(Jan): the cursors are drawn blue (Q_min) and red (Q_max) dashed, and the Q range is shown
% as a "Fit Q range: 1.2 – 4 Å⁻¹" label, not boxes.

(waxs-peak-fit-no-limits)=
### No limits?

When **checked**, the Lo limit / Hi limit columns are hidden for all background and peak
parameters, and the fit runs unconstrained. The limits are remembered and restored when the
checkbox is unchecked again.

(waxs-peak-fit-background)=
### Background

Choose a background shape from the combo box (default **SNIP**).

**Polynomial shapes** — optimised simultaneously with the peaks by `scipy.optimize.curve_fit`:

| Shape | Parameters |
|-------|-----------|
| Constant | bg0 |
| Linear | bg0, bg1 |
| Cubic | bg0 … bg3 |
| 5th Polynomial | bg0 … bg5 |

Each row has **Value**, **Fit?** checkbox, and (when limits are shown) **Lo limit** / **Hi
limit** fields.

**Adaptive (data-driven) shapes** — estimated directly from the data before fitting so the
background is guaranteed to stay at or below the measured intensities. These are particularly
effective for XRD / powder-diffraction spectra with many overlapping peaks. Only one numeric
parameter is shown (fraction of the data length):

| Shape | Parameter | Description |
|-------|-----------|-------------|
| SNIP | Half-width (fraction) | Iterative peak-clipping (Statistics-sensitive Non-linear Iterative Peak-clipping); standard for XRD |
| Rolling Quantile Spline | Window + Quantile | Rolling percentile filter + CubicSpline; `quantile=0` gives a rolling minimum |
| Rolling Ball | Radius (fraction) | Morphological grey-erosion + grey-dilation; equivalent to rolling a ball under the spectrum |

> **Tip:** Try **SNIP** first for XRD data. If the background is very steep or irregular,
> **Rolling Quantile Spline** with a small quantile (0.05–0.15) often gives a cleaner result.

(waxs-peak-fit-peak-finding)=
### Peak Finding

| Control | Meaning | Default | Tips |
|---|---|---|---|
| **Min. Prominence (frac.)** | Minimum peak prominence relative to data range | 0.05 | |
| **Min. FWHM (Å⁻¹)** | Peaks narrower than this are ignored | 0.001 | |
| **Max. FWHM (Å⁻¹)** | Peaks wider than this are ignored | 0.500 | |
| **Min. Distance (Å⁻¹)** | Minimum separation between detected peaks | 0.005 | |
| **BG smooth window (%)** | Savitzky-Golay window as % of data length | 15 | |
| **Find Peaks** | Auto-detect peaks and populate the Peaks scroll area. The model is graphed immediately after detection. | | |

(waxs-peak-fit-peaks)=
### Peaks (scroll area)

Each detected or manually added peak appears as a collapsible row:

| Control | Meaning | Default | Tips |
|---|---|---|---|
| **Shape** | Gauss, Lorentz, Pseudo-Voigt, LogNormal | Gauss | |
| **A** | Amplitude/height | | |
| **Q0** | Position, Å⁻¹ | | Mouse wheel: fixed step of **0.001 Å⁻¹** per notch (trackpad-smooth; hold **Shift** for 0.01 Å⁻¹ coarser step) |
| **FWHM** | Full width at half maximum, Å⁻¹ | | Mouse wheel as for Q0 |
| **eta** | Mixing fraction for Pseudo-Voigt only | | |
| **Fit?** | Per parameter; uncheck to hold fixed | | |
| **Lo / Hi** | Limit fields | | Hidden when "No limits?" is checked |
| **Remove** (red) | Deletes this peak | | |
| **Sort by Q** | Re-orders the peak rows by position | | |

**Click "Add Peak Manually"** to insert a new Gaussian peak at the centre of the current Q view.
Alternatively, **right-click on the graph** and choose "Add Peak at Q = …".

% CHECK(Jan): there is no "Add Peak Manually" button in the current panel (code comment: add and
% remove are available via right-click on the graph). The Shape combo also offers "Voigt".

(waxs-peak-fit-fit-controls)=
### Fit controls

| Control | Meaning | Default | Tips |
|---|---|---|---|
| **Graph Model** (light green) | Evaluate model at current parameters and draw overlays | | |
| **Fit** (dark green) | Run `scipy.optimize.curve_fit`; update fields with fitted values | | |
| **Revert** (orange) | Restore parameters to their pre-Fit values | | Labelled **Revert back**, in the Results group |
| **Fix Limits** | Set all fit limits to ±5× the current peak parameter values | | Convenient starting point before running a constrained fit |
| **Weighting** | σ used in the fit, see below | Equal (σ = 1) | |

| Weighting | σ used in fit | When to use |
|------|--------------|-------------|
| 1/σ² (standard) | measured *dI* (or 1 if unavailable) | Uses actual measurement uncertainties; can be slow when σ range is large |
| Equal (σ = 1) | 1 for all points | **Default.** Prevents background points from dominating; robust when σ range is large |
| Relative (σ = dI/I) | *dI/I* (relative error) | Emphasises narrow peaks; useful when peak-to-background ratio is large |

(waxs-peak-fit-presearch)=
### Q₀ pre-search

| Control | Meaning | Default | Tips |
|---|---|---|---|
| **Cross-corr. (global shift)** | Enable FFT cross-correlation pre-step | off | |
| **Per-peak Q₀ scan** | Enable independent grid scan per peak | off | |
| **Search window (Å⁻¹)** | Half-width used by both methods | 0.050 | Labelled **Window (Å⁻¹)** |
| **Grid steps** | Trial positions per peak in the scan | 50 | |

> **Tip:** Enable **both** steps for maximum robustness: cross-correlation handles a large global
> shift quickly, then the per-peak scan fine-tunes individual positions within the remaining
> window.

The pre-search settings are saved with **Save State**, exported to `pyirena_config.json` via
**Export Parameters**, and honoured by the batch API (`fit_waxs_peaks` / `fit_waxs`). How the two
methods work is described under [Method notes](#waxs-peak-fit-presearch-method).

(diffraction-lines-tab)=
## Diffraction Lines tab

The **Diffraction Lines** tab overlays theoretical powder diffraction stick patterns from CIF
files on the main graph. This is a phase-identification aid — it shows expected peak positions
and relative heights for one or more known crystal structures so you can compare them visually
against the measured data. It is **not** a Rietveld refinement; no profile fitting is performed.

```{figure} img/waxs-peak-fit/waxs-peak-fit-04-diffraction-lines.png
:alt: WAXS Peak Fit window with the Diffraction Lines tab selected: one CIF (Si) in the list with hkl labels switched on, and its sticks with Miller indices drawn over the fitted data.

The Diffraction Lines tab with `testData/Si.cif` imported and **hkl** on. Silicon is not in this
synthetic dataset; the sticks only show what the overlay looks like.
```

% SCREENSHOT waxs-peak-fit-05-phase-id: real powder pattern from testData/ with two candidate
% phases imported, sticks aligned with the measured peaks (Jan picks sample and CIFs).

Patterns are computed using [Dans_Diffraction](https://github.com/DanPorter/Dans_Diffraction)
(Apache-2.0 licence). It must be installed separately:

```
pip install Dans-Diffraction
```

or, if you use the conda environment shipped with pyirena:

```
conda env update -f environment.yml
```

### Wavelength and distance correction

| Control | Description |
|---------|-------------|
| **λ (Å)** spinbox | X-ray wavelength used to convert d-spacings to Q (Q = 4π sin θ / λ). Default 1.5406 Å (Cu Kα). |
| **Auto from file** checkbox | When checked, the wavelength is read automatically from the loaded NXcanSAS file (`/entry/instrument/wavelength`). Editing the spinbox manually unchecks this so the value is not overwritten when the next file loads. |
| **ΔL (mm)** spinbox | **Sample-to-detector distance correction.** Range ±10 mm. If your instrument's actual distance L_true differs from the calibration distance L_cal used during data reduction by ΔL mm, this shift moves the overlay sticks (and hkl labels) to where the peaks *actually* appear in your data without touching the data itself. Uses the exact non-linear geometry `tan(2θ_app) = (L_true/L_cal) · tan(2θ_true)`, which is important at WAXS angles where a uniform Q-shift would be wrong. **Default 0** (no correction). |
| **L_cal (mm)** spinbox | **Calibration sample-to-detector distance** (the distance value used during the original I(Q) reduction). Default 100 mm. When **Auto from file** (checkbox below) is checked, this is read automatically from the NXcanSAS file's sample-to-detector distance metadata and updated live. Manual override is always available. |
| **Auto from file** checkbox (L_cal) | When checked, L_cal is auto-loaded from the NXcanSAS file's instrument metadata (tries several standard NeXus paths; honors units attributes; applies a metres↔mm heuristic if units are missing). Falls back to the manual spinbox value if no metadata is found. Mirrors the wavelength Auto behaviour. |

**Calibrating ΔL at the beamline:** Load a known sample with strong, well-resolved diffraction
peaks. Compare the tabulated peak positions (shown by the overlay sticks) to the measured peaks.
If peaks appear at higher Q than tabulated, increase ΔL (true distance is larger than L_cal); if
lower Q, decrease ΔL (true distance is smaller). When the sticks line up with the peaks, |ΔL| is
your distance error.

### CIF files list

Each imported CIF appears as a row:

| Control | Description |
|---------|-------------|
| Visibility checkbox (leftmost) | Show or hide this phase's sticks without removing it. |
| Color swatch button | Click to open a color-picker dialog and change the phase color. |
| Phase name | Chemical formula or CIF filename stem read from the CIF on import. Hover for full file path. |
| Scale spinbox (0.01 – 10×) | Manual multiplier applied on top of the automatic scaling. The automatic scaling sets the tallest stick to the measured data's peak intensity within the same Q span. Adjust this spinbox to make sticks taller or shorter. |
| **d×** spinbox (0.9000 – 1.1000) | **Per-CIF lattice / d-spacing scale (a/a₀).** `1.0000` = use the tabulated CIF values exactly. Use this to correct alloy lattice deviations from the nominal element values — for example, an Al-based alloy whose real lattice parameter is 0.5 % larger than pure Al would use `d× = 1.0050`. All peak Q's are divided by this factor (since Q = 2π/d, scaling every d by a/a₀ scales every Q by a₀/a). **Exact for cubic systems** (Al, Cu, Ni, FCC/BCC steels, …); **first-order accurate for non-cubic** crystals (it captures the dominant volumetric strain but ignores anisotropic strain). Default `1.0000`. Right-click the row → **Reset d× to 1.0** to clear. |
| **hkl** checkbox | Show Miller-index labels (e.g. `(110)`) above each stick. Off by default; only labels for reflections with intensity ≥ 5 % of the phase maximum are shown to avoid clutter. |
| **×** button | Remove this CIF from the list. |

**Right-click** any row for **Reset d× to 1.0** (when d× ≠ 1) and **Delete** context-menu
entries.

**d× vs ΔL — when to use which:** d× is **per-CIF** and corrects an *intrinsic* property of the
crystal (lattice parameter); ΔL is **global** and corrects an *instrument* property
(sample-to-detector distance). If only one phase is misaligned with the data, use d× on that
phase. If all phases are uniformly off in the same direction across the full Q range, suspect
ΔL. The two compose: d× is applied first, then ΔL.

### Buttons

| Button | Action |
|--------|--------|
| **Import CIF…** (blue) | Open a file picker to select one or more `.cif` files. The picker reopens at the last used folder. Invalid or unparseable CIFs produce an error dialog; all others are added to the list. |
| **AMCSD…** | Open the [American Mineralogist Crystal Structure Database](https://www.rruff.net/amcsd/) in your default browser. |
| **COD…** | Open the [Crystallography Open Database](https://www.crystallography.net/cod/search.html) in your default browser. |
| **Clear all CIFs** (orange) | Prompts for confirmation, then removes all CIF entries from the list and clears the overlays. |

### Persistence

All Diffraction Lines settings — CIF file paths, phase names, colors, visibility, intensity
scales, per-CIF d×, hkl-toggle state, wavelength, distance correction (ΔL and L_cal), and the last
folder used for the file picker — are saved automatically to the pyirena state file whenever a
control changes. They are restored the next time the WAXS Peak Fit tool is opened.

CIF files that have been moved or deleted since the last session are silently skipped on state
restore.

### Phase identification workflow

1. Switch to the **Diffraction Lines** tab.
2. Verify the wavelength matches your data (auto-detected from NXcanSAS files).
3. Click **AMCSD…** or **COD…** to open a database in your browser; search for a candidate phase
   and download its CIF file.
4. Click **Import CIF…**, select the downloaded CIF. Sticks appear on the main graph immediately,
   auto-scaled to the data peak.
5. Compare stick positions to the measured peaks. Use the scale spinbox to adjust stick height if
   the auto-scale is off.
6. Enable **hkl** to verify individual reflection assignments.
7. Add more CIFs for additional candidate phases; each gets a distinct color.
8. Uncheck the visibility box for phases you want to temporarily hide.
9. Phases you do not need can be deleted with the **×** button or right-click → **Delete**.

(waxs-peak-fit-results)=
## Results and outputs

| Button | Action |
|--------|--------|
| Reset to Defaults | Clear all peaks; reset background to SNIP, weighting to Equal |
| Save State | Persist GUI state to pyIrena state file |
| Store in File | Save fit results to `entry/waxs_peakfit_results` in the HDF5 file |
| Copy results | Put the current peaks on the clipboard as Markdown — per-peak A, Q₀, FWHM, η with 1σ uncertainties and the integrated area, plus χ², background shape and Q range. Describes the peak rows on screen, so it works before a fit as well as after one |
| Save report… | The same text as a `.md` file. **Ctrl-click** (⌘-click on macOS) also saves the graph as a PNG beside it and embeds it as a figure |
| Results to graphs | Overlay a text annotation with fitted values on the main plot |
| Export Parameters | Write current parameters to `pyirena_config.json` in the data directory |
| Import Parameters | Load parameters from a previously exported JSON file |
| Load Setup from File… | Restore every WAXS Peak Fit control (peak list, background, q-range, …) from a NXcanSAS file previously saved by pyirena or by the pyirena-ai agent |

% CHECK(Jan): "Save State", "Export Parameters" and "Import Parameters" are not in the current
% panel; the JSON buttons are "Save params to JSON" / "Load params from JSON".

(waxs-peak-fit-hdf5)=
### HDF5 output format

Results are stored in `entry/waxs_peakfit_results` (NXprocess group), with the full panel setup
embedded so **Load Setup from File…** can restore it; see {ref}`hdf5-results-groups`.

| Path | Content |
|------|---------|
| `.attrs['n_peaks']` | Number of fitted peaks |
| `.attrs['bg_shape']` | Background shape string |
| `.attrs['chi_squared']` | χ² |
| `.attrs['reduced_chi_squared']` | χ² / dof |
| `Q` | Q array in fit range (Å⁻¹) |
| `I_fit` | Total model curve |
| `I_bg` | Background-only curve |
| `residuals` | `(I_data − I_fit) / σ` |
| `intensity_data` | Measured I in fit range |
| `intensity_error` | Measurement uncertainty σ |
| `background/bg0` … `bgN` | Background coefficients (float64 scalars) |
| `background_std/bg0` … `bgN` | Uncertainties from covariance matrix |
| `peak_01/params/Q0`, `A`, `FWHM`, (`eta`) | Fitted peak parameters |
| `peak_01/params_std/…` | Parameter uncertainties |
| `peak_01/Q_peak` | Q array ±5·FWHM around Q0 |
| `peak_01/I_peak` | Individual peak curve on that sub-range |

(waxs-peak-fit-scripting)=
## From Python or AIDA

::::{tab-set}
:::{tab-item} Python
The control API drives the same model as the window. Run from the repository root:

```python
from pyirena.api.control import (
    find_waxs_peaks, list_waxs_peaks, open_dataset, remove_waxs_peak, run_waxs_fit,
    save_waxs_fit, select_waxs_model,
)

sid = open_dataset("validationData/waxs/waxs_three_gauss.h5")["session_id"]
select_waxs_model(sid, bg_shape="Linear")
find_waxs_peaks(sid, prominence_frac=0.05)

# Remove false detections: peaks lower than 20 % of the tallest one.
peaks = list_waxs_peaks(sid)["peaks"]
height = {p["index"]: p["parameters"][0]["value"] for p in peaks}      # "A"
for index in sorted(height, reverse=True):
    if height[index] < 0.2 * max(height.values()):
        remove_waxs_peak(sid, index)

fit = run_waxs_fit(sid, weight_mode="equal")
print(fit["reduced_chi_squared"])                  # about 0.84
for p in fit["peaks"]:
    values = {q["name"]: round(q["value"], 4) for q in p["parameters"]}
    print(p["shape"], values)                      # Q0 1.85, 2.55, 3.05
save_waxs_fit(sid, output_path="waxs_three_gauss_fit.h5")   # input file stays untouched
```
:::
:::{tab-item} AIDA / MCP
Ask in plain language, for example *"Fit the peaks in validationData/waxs/waxs_three_gauss.h5 on
a linear background."* The agent calls `pyirena_ctrl_open_dataset`, then `pyirena_call` with
`select_waxs_model`, `find_waxs_peaks`, `list_waxs_peaks`, `run_waxs_fit` and `save_waxs_fit`.
The call sequence is in the
[AI tools reference](../reference/ai-tools-reference.md#control-tools-reference--waxs-peak-fit-dispatcher-category-waxs).
`pyirena_read_waxs_peakfit(path)` returns stored results.
:::
::::

(waxs-peak-fit-batch)=
### Batch fitting after interactive setup

Once `pyirena_config.json` exists in the data directory (created by **Export Parameters**, which
includes the pre-search settings):

```text
from pyirena.batch import fit_waxs
from pathlib import Path

data_files = sorted(Path("data/").glob("*_waxs.h5"))
for f in data_files:
    result = fit_waxs(f, "data/pyirena_config.json")
    if result and result['success']:
        print(f"{f.name}: {result['n_peaks']} peaks, "
              f"reduced-χ²={result['reduced_chi2']:.4g}")
```

Pre-search can also be activated or overridden from the config dict directly:

```text
from pyirena.batch import fit_waxs_peaks
import json, numpy as np

config = json.loads(Path("data/pyirena_config.json").read_text())["waxs_peakfit"]

# Enable both pre-search methods with a ±0.03 Å⁻¹ window
config["presearch"] = {
    "cross_corr":    True,
    "per_peak_scan": True,
    "search_window": 0.030,
    "n_steps":       60,
}

for f in sorted(Path("data/").glob("*_waxs.h5")):
    result = fit_waxs_peaks(f, config=config)
    if result and result['success']:
        print(f"{f.name}: reduced-χ²={result['reduced_chi2']:.4g}")
```

Or, from the Data Selector, select all files and click **WAXS Peaks (script)** to run the batch
fit in the GUI's background worker.

% CHECK(Jan): the two batch examples above are shown as text, not Python: they need a
% data/pyirena_config.json and *_waxs.h5 files that do not exist in the repository, so they cannot
% run as written (phase-5 rule). The second one also uses Path without importing it.

(waxs-peak-fit-method)=
## Method notes

(waxs-peak-fit-presearch-method)=
### Q₀ pre-search

Least-squares fitting of sharp peaks is very sensitive to the initial Q0 estimate: if the
starting position is more than ~1–2 peak widths from the true centre, the optimizer sees a
near-zero gradient and fails to converge. The **Q₀ pre-search** group runs one or both pre-steps
*before* the main `curve_fit` call to find better Q0 starting values automatically.

**Cross-corr. (global shift).** Uses FFT cross-correlation between the current model and the
measured data to detect a single rigid Q-shift common to all peaks. The detected shift (clipped to
±**Search window**) is added to every Q0 before fitting.

*When to use:* peaks shift together due to a lattice parameter change with temperature, pressure,
or composition (e.g. all peaks move by the same ΔQ). Very fast (one FFT).

**Per-peak Q₀ scan.** For each peak independently, Q0 is scanned over
`[Q0 − search_window, Q0 + search_window]` in **Grid steps** trial positions. At each trial
position the total model (all peaks + background) is evaluated and compared with the data; the
position with the lowest sum-of-squares residual becomes the new Q0 starting value. All other
peaks remain at their current Q0 during the scan of each peak.

*When to use:* individual peaks may have shifted by different amounts (e.g. two phases with
different thermal expansion), or when the cross-correlation shift is insufficient for very narrow
peaks.

(waxs-peak-fit-quality)=
### Fit quality

Residuals are displayed **rescaled** as r' = r / σ(robust): the normalized residual divided by a
robust (MAD-based) estimate of the actual noise scale, so the scatter is judged against the data's
own noise floor rather than the (often mis-scaled) reported σ. The fit status line also reports
**σ-scale** (how many × the actual scatter exceeds the reported σ), the realistic reduced-χ²
floor, and the largest fractional misfit **max|(I−M)/I|** — a σ-independent gross-misfit
backstop.

See the **[Fit Quality Metrics guide](../methods/fit-quality-metrics.md)** for full
interpretation of these numbers.

## See also

- [Batch API reference — `fit_waxs`](../guides/batch-api.md#fit_waxs)
- [NXcanSAS file format](../reference/hdf5-nxcansas-structure.md)
- [Data Selector quick start](../getting-started/quickstart.md)
