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

Choose a background shape from the combo box:

| Shape | Parameters |
|-------|-----------|
| Constant | bg0 |
| Linear | bg0, bg1 |
| Cubic | bg0 … bg3 |
| 5th Polynomial | bg0 … bg5 |

Each row has **Value**, **Fit?** checkbox, and (when limits are shown) **Lo
limit** / **Hi limit** fields.

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

---

### Additional / Results

| Button | Action |
|--------|--------|
| Reset to Defaults | Clear all peaks; reset background to Constant |
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
3. Choose a **Background** shape (start with Cubic or Linear).
4. Click **Find Peaks** — the tool detects peaks automatically.
5. Inspect the overlay: use mouse-wheel on **Q0** and **FWHM** to fine-tune
   initial positions.  Right-click → **Add Peak** for any missed peaks; use
   **Remove** for false detections.
6. Click **Graph Model** to preview the current parameter set without fitting.
7. Click **Fit** — fitted values replace the initial guesses.
8. If the fit diverged, click **Revert**, adjust limits or initial values, and
   refit.
9. Click **Store in File** to save results.
10. Click **Export Parameters** to write `pyirena_config.json` for batch use.

### Batch fitting after interactive setup

Once `pyirena_config.json` exists in the data directory:

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
