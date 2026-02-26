# WAXS Peak Fit — GUI Guide

The **WAXS Peak Fit** tool fits diffraction peaks (Gaussian, Lorentzian,
Pseudo-Voigt, or log-normal profiles) rising above a smooth background in
Wide-Angle X-ray Scattering data.  All plots are **linear/linear** (intensity
vs. Q in Å⁻¹).

---

## Opening the tool

In the Data Selector:

1. Select one or more HDF5 data files.
2. Click **WAXS Peaks (GUI)** — opens the `WAXSPeakFitPanel` for the first
   selected file.

---

## Panel layout

The window is split into a **left control panel** (≈430 px wide) and a
**right graph area**.

### Graph area

| Panel | Content |
|-------|---------|
| Top (main) | Data scatter, total model (red), per-peak overlays (colors), background (dashed gray) |
| Bottom (residuals) | Normalized residuals `(I − fit) / σ`; dashed zero line; x-axis linked to main |

**Right-click** on the main graph to:
- **Add Peak at Q = X.XXX Å⁻¹** — inserts a peak row initialized at that Q
- **Show data as line / Show data as points** — toggle scatter vs. connected line
- **Hide / Show error bars** — toggle ±σ vertical bars on data
- **Save graph as JPEG…** — export at 1600 px width

### Vertical cursors

Two vertical cursors (cyan dashed lines) define the **Q fit range**.  The
read-only Q-min / Q-max display boxes below the title bar update as you drag
the cursors.  Only data inside this range is used for fitting; residuals
outside are shown as NaN (blank).

---

## Left panel: controls (top → bottom)

### No limits?

When **checked**, the Lo limit / Hi limit columns are hidden for all
background and peak parameters, and the fit runs unconstrained.  The limits
are remembered and restored when the checkbox is unchecked again.

---

### Background

Choose a background shape from the combo box.

**Polynomial shapes** — optimised simultaneously with the peaks by
`scipy.optimize.curve_fit`:

| Shape | Parameters |
|-------|-----------|
| Constant | bg0 |
| Linear | bg0, bg1 |
| Cubic | bg0 … bg3 |
| 5th Polynomial | bg0 … bg5 |

Each row has **Value**, **Fit?** checkbox, and (when limits are shown) **Lo
limit** / **Hi limit** fields.

**Adaptive (data-driven) shapes** — estimated directly from the data before
fitting so the background is guaranteed to stay at or below the measured
intensities.  These are particularly effective for XRD / powder-diffraction
spectra with many overlapping peaks.  Only one numeric parameter is shown
(fraction of the data length):

| Shape | Parameter | Description |
|-------|-----------|-------------|
| SNIP | Half-width (fraction) | Iterative peak-clipping (Statistics-sensitive Non-linear Iterative Peak-clipping); standard for XRD |
| Rolling Quantile Spline | Window + Quantile | Rolling percentile filter + CubicSpline; `quantile=0` gives a rolling minimum |
| Rolling Ball | Radius (fraction) | Morphological grey-erosion + grey-dilation; equivalent to rolling a ball under the spectrum |

> **Tip:** Try **SNIP** first for XRD data.  If the background is very steep
> or irregular, **Rolling Quantile Spline** with a small quantile (0.05–0.15)
> often gives a cleaner result.

---

### Peak Finding Parameters

| Field | Default | Description |
|-------|---------|-------------|
| Min. Prominence (fraction) | 0.05 | Minimum peak prominence relative to data range |
| Min. FWHM (Å⁻¹) | 0.001 | Peaks narrower than this are ignored |
| Max. FWHM (Å⁻¹) | 0.500 | Peaks wider than this are ignored |
| Min. Distance (Å⁻¹) | 0.005 | Minimum separation between detected peaks |
| BG smoothing window (%) | 15 | Savitzky-Golay window as % of data length |

Click **Find Peaks** (blue) to auto-detect peaks and populate the Peaks
scroll area.  The model is graphed immediately after detection.

---

### Peaks (scroll area)

Each detected or manually added peak appears as a collapsible row showing:

- **Shape** combo box: Gauss, Lorentz, Pseudo-Voigt, LogNormal
- Parameter rows: **A** (amplitude/height), **Q0** (position, Å⁻¹),
  **FWHM** (full width at half maximum, Å⁻¹), **eta** (mixing fraction for
  Pseudo-Voigt only)
- **Fit?** checkbox per parameter (uncheck to hold fixed)
- **Lo / Hi** limit fields (hidden when "No limits?" is checked)
- **Remove** button (red) — deletes this peak

**Mouse-wheel on Q0 or FWHM fields:** fixed step of **0.001 Å⁻¹** per notch
(trackpad-smooth; hold **Shift** for 0.01 Å⁻¹ coarser step).

**Click "Add Peak Manually"** to insert a new Gaussian peak at the centre of
the current Q view.  Alternatively, **right-click on the graph** and choose
"Add Peak at Q = …".

---

### Fit controls

| Button | Color | Action |
|--------|-------|--------|
| Graph Model | light green | Evaluate model at current parameters and draw overlays |
| Fit | dark green | Run `scipy.optimize.curve_fit`; update fields with fitted values |
| Revert | orange | Restore parameters to their pre-Fit values |

**Weighting** combo (below the buttons):

| Mode | σ used in fit | When to use |
|------|--------------|-------------|
| 1/σ² (standard) | measured *dI* (or 1 if unavailable) | Uses actual measurement uncertainties; can be slow when σ range is large |
| Equal (σ = 1) | 1 for all points | **Default.** Prevents background points from dominating; robust when σ range is large |
| Relative (σ = dI/I) | *dI/I* (relative error) | Emphasises narrow peaks; useful when peak-to-background ratio is large |

---

### Q₀ pre-search

Least-squares fitting of sharp peaks is very sensitive to the initial Q0
estimate: if the starting position is more than ~1–2 peak widths from the
true centre, the optimizer sees a near-zero gradient and fails to converge.
The **Q₀ pre-search** group runs one or both pre-steps *before* the main
`curve_fit` call to find better Q0 starting values automatically.

#### Cross-corr. (global shift)

Uses FFT cross-correlation between the current model and the measured data to
detect a single rigid Q-shift common to all peaks.  The detected shift (clipped
to ±**Search window**) is added to every Q0 before fitting.

**When to use:** peaks shift together due to a lattice parameter change with
temperature, pressure, or composition (e.g. all peaks move by the same ΔQ).
Very fast (one FFT).

#### Per-peak Q₀ scan

For each peak independently, Q0 is scanned over
`[Q0 − search_window, Q0 + search_window]` in **Grid steps** trial positions.
At each trial position the total model (all peaks + background) is evaluated and
compared with the data; the position with the lowest sum-of-squares residual
becomes the new Q0 starting value.  All other peaks remain at their current Q0
during the scan of each peak.

**When to use:** individual peaks may have shifted by different amounts (e.g.
two phases with different thermal expansion), or when the cross-correlation
shift is insufficient for very narrow peaks.

| Control | Default | Description |
|---------|---------|-------------|
| Cross-corr. (global shift) | off | Enable FFT cross-correlation pre-step |
| Per-peak Q₀ scan | off | Enable independent grid scan per peak |
| Search window (Å⁻¹) | 0.050 | Half-width used by both methods |
| Grid steps | 50 | Trial positions per peak in the scan |

> **Tip:** Enable **both** steps for maximum robustness: cross-correlation
> handles a large global shift quickly, then the per-peak scan fine-tunes
> individual positions within the remaining window.

The pre-search settings are saved with **Save State**, exported to
`pyirena_config.json` via **Export Parameters**, and honoured by the batch
API (`fit_waxs_peaks` / `fit_waxs`).

---

### Additional / Results

| Button | Action |
|--------|--------|
| Reset to Defaults | Clear all peaks; reset background to SNIP, weighting to Equal |
| Save State | Persist GUI state to pyIrena state file |
| Store in File | Save fit results to `entry/waxs_peakfit_results` in the HDF5 file |
| Results to graphs | Overlay a text annotation with fitted values on the main plot |
| Export Parameters | Write current parameters to `pyirena_config.json` in the data directory |
| Import Parameters | Load parameters from a previously exported JSON file |

---

## Workflow

### Typical interactive session

1. Open **WAXS Peaks (GUI)** from the Data Selector.
2. Drag the Q cursors to bracket the diffraction peaks of interest.
3. Choose a **Background** shape.
   - For XRD / powder data try **SNIP** first.
   - For data with a gently-varying background, start with **Cubic** or **Linear**.
4. Click **Find Peaks** — the tool detects peaks automatically.
5. Inspect the overlay: use mouse-wheel on **Q0** and **FWHM** to fine-tune
   initial positions.  Right-click → **Add Peak** for any missed peaks; use
   **Remove** for false detections.
6. Select a **Weighting** mode.  If peaks have few points compared to the
   background, try **Equal** or **Relative**.
7. If peaks are narrow or their positions may have drifted, enable one or both
   **Q₀ pre-search** options.  For a data series where peaks shift coherently,
   tick **Cross-corr.** and set the **Search window** to the expected maximum
   drift.  For very narrow peaks that may have shifted independently, also tick
   **Per-peak Q₀ scan**.
8. Click **Graph Model** to preview the current parameter set without fitting.
9. Click **Fit** — pre-search (if enabled) runs first, then `curve_fit`.
   Fitted values replace the initial guesses.
10. If the fit diverged, click **Revert**, adjust limits or initial values, and
    refit.
11. Click **Store in File** to save results.
12. Click **Export Parameters** to write `pyirena_config.json` for batch use
    (pre-search settings are included).

### Batch fitting after interactive setup

Once `pyirena_config.json` exists in the data directory (created by
**Export Parameters**, which includes the pre-search settings):

```python
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

```python
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

Or, from the Data Selector, select all files and click **WAXS Peaks
(script)** to run the batch fit in the GUI's background worker.

---

## HDF5 output format

Results are stored in `entry/waxs_peakfit_results` (NXprocess group).

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

---

## See also

- [Batch API reference — `fit_waxs`](batch_api.md#fit_waxs)
- [NXcanSAS file format](NXcanSAS_UnifiedFit_Format.md)
- [Data Selector quick start](gui_quickstart.md)
