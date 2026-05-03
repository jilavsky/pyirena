# SAXS Morph Tool (3D Voxelgram)

The SAXS Morph tool generates a **3D voxelgram of a two-phase porous structure**
from your experimental I(Q) using the Gaussian Random Fields (GRF) method,
then computes the model I(Q) from that voxelgram and refines parameters
(volume fraction, contrast, Power-law + Flat background) by fitting back to
the data.

The result is a 3D phase map you can rotate, slice, and inspect — useful for
intuitive visualisation of the *kind* of microstructure consistent with a SAS
profile.

References:

- Berk, N. F. (1991). *Phys. Rev. E* **51**, 4141.
- Roberts, A. P. (1997). *Phys. Rev. E* **55**, R1286.
- Levitz, P. (2007). *Modelling Simul. Mater. Sci. Eng.* **15**, S2.
- Cherny, A. Yu., et al. (2013). *J. Appl. Cryst.* **46**, 365.
- Igor Pro reference: [`IR3_3DTwoPhaseSolid.ipf`](https://github.com/jilavsky/SAXS_IgorCode/blob/master/User%20Procedures/Irena/IR3_3DTwoPhaseSolid.ipf)

---

## Contents

1. [Overview](#overview)
2. [Launching](#launching)
3. [Panel layout](#panel-layout)
4. [Workflow](#workflow)
5. [Parameter reference](#parameter-reference)
6. [Background tabs (Power-law + Flat)](#background-tabs-power-law--flat)
7. [3D viewer controls](#3d-viewer-controls)
8. [2D slice viewer](#2d-slice-viewer)
9. [Pop-out mode](#pop-out-mode)
10. [Fitting and MC uncertainty](#fitting-and-mc-uncertainty)
11. [Performance and memory notes](#performance-and-memory-notes)
12. [HDF5 output format](#hdf5-output-format)
13. [Limitations (v1)](#limitations-v1)
14. [Scripting API](#scripting-api)

---

## Overview

Algorithm pipeline for one Graph Model / Fit iteration:

1. Subtract user background (Power-law + Flat) from the data → `I_corr(Q)`.
2. Sinc-transform `Q² · I_corr(Q)` → Debye autocorrelation γ(r); normalise
   γ(0) = 1.
3. Sinc-transform γ(r) · r² → spectral density F(k); clip negative bins.
4. Threshold parameter α = √2 · erfinv(1 − 2φ).
5. Generate 3D Gaussian white noise → FFT → multiply by √F(|k₃ᴅ|) →
   inverse FFT → real scalar field.
6. Threshold the scalar field at α·σ → binary uint8 cube.
7. FFT the (mean-subtracted) voxelgram → spherically average |F(k)|² →
   resample onto Q grid → multiply by contrast.
8. Add background back and compute χ² against data.

The fit varies any of `volume_fraction`, `contrast` (when not linked),
`power_law_B`, `power_law_P`, `background` to minimise χ².

---

## Launching

**From the Data Selector:**

1. Select one or more HDF5 files in the file list.
2. Click **SAXS Morph (GUI)** (purple button, row 8).

The tool opens with the selected file pre-loaded and all settings restored
from the previous session.

**From the command line:**

```bash
pyirena-saxsmorph path/to/file.h5
# or
python -m pyirena.gui.saxs_morph_panel path/to/file.h5
```

---

## Panel layout

```
┌─────────────────────────────┬────────────────────────────────────────────┐
│  Left controls (~440 px)     │  Right graphs / 3D                         │
│                              │                                            │
│  SAXS Morph    [? Help]      │  ┌─────────────────────────────────────┐   │
│  Data file:  ………  [Open…]   │  │  I(Q)  log-log                      │   │
│  Q range: 0.001 to 0.3 Å⁻¹   │  │   • Data scatter                     │   │
│  ┌─ Voxel grid ─────────┐    │  │   • Data − background scatter        │   │
│  │ Cube side (fit):  128│    │  │   • Red model line                   │   │
│  │ Cube side (rendr):256│    │  │   • Two cursors (red Qmin / blue)    │   │
│  │ Box size [Å]:    1000│    │  └─────────────────────────────────────┘   │
│  │ RNG seed:  (blank=…) │    │  ┌────────────────┬───────────────────┐   │
│  └──────────────────────┘    │  │ 2D slice       │ 3D PyVista        │   │
│  ┌─ Two-phase ─────────┐     │  │ XY/XZ/YZ combo │ rotate w/ mouse   │   │
│  │ φ:        Fit lo hi │     │  │ slider         │ right-click menu  │   │
│  │ Δρ²:      Fit lo hi │     │  │  [Pop out ⤢]   │  [Pop out ⤢]      │   │
│  │ ☑ Link via invariant│     │  └────────────────┴───────────────────┘   │
│  └─────────────────────┘     │                                            │
│  ┌ Background tabs ─────┐    │  Status bar                                │
│  │ Power-law │  Flat    │    │  Fit done. χ² = 12.34, φ_actual = 0.31    │
│  │  B        Fit lo hi  │    │                                            │
│  │  P        Fit lo hi  │    │                                            │
│  └──────────────────────┘    │                                            │
│  ☐ No limits (Nelder-Mead)   │                                            │
│  [Graph Model] [Fit] [Cancel]│                                            │
│  Passes: 10  [MC] [Revert]   │                                            │
│  Result block (chi², φ, …)   │                                            │
│  [Save Result to HDF5…]      │                                            │
└─────────────────────────────┴────────────────────────────────────────────┘
```

---

## Workflow

1. **Load a data file** (Data Selector → SAXS Morph (GUI), or `Open…` button
   in the panel).
2. **Set background** in the Power-law / Flat tabs. For most porous-glass /
   Vycor samples, fit a Q⁻⁴ Porod tail at high Q first using a separate tool,
   then plug B and P into the Power-law tab here.
3. **Set Q range** by dragging the red and blue cursors on the I(Q) plot.
   The fit only sees data within `[Qmin, Qmax]`.
4. **Set initial parameters**:
   - `φ` (volume fraction): your best guess of the minority-phase fraction.
   - `Δρ²` (contrast): leave the **Link via invariant** checkbox enabled to
     have contrast derived automatically from the Porod invariant — this is
     the recommended default.
   - Voxel size (fit): **128³** is the practical sweet spot during fitting;
     drop to 64³ for very fast iteration on slow machines.
   - Voxel size (render): **256³** is a good default for the saved result;
     384³ / 512³ produce nicer renders at the cost of disk space.
   - Box size (Å): Make sure `box_size_A / voxel_size` is comparable to the
     correlation length you expect to resolve (e.g. 1000 Å / 128 ≈ 8 Å pitch).
   - RNG seed: leave blank for stochastic exploration; set an integer to
     freeze the voxelgram for figures or to reproduce a published result.
5. Click **Graph Model**. After ~5–30 s (depending on voxel size), you see:
   - The model I(Q) as a red curve overlaying the data.
   - A 2D slice in the bottom-left viewer (drag the slider to scrub depth).
   - A 3D isosurface in the bottom-right PyVista viewer (rotate with mouse).
6. Iterate parameters until the model resembles the data, then click **Fit**.
7. After convergence, click **Calc. Uncertainty (MC)** to estimate parameter
   error bars (the spinbox sets the number of passes; 10 is a reasonable
   default).
8. Click **Save Result to HDF5…** to write the compressed voxelgram +
   parameters into the source HDF5 file (or a new file).

---

## Parameter reference

| Parameter | Units | Default | Range | Notes |
|---|---|---|---|---|
| `voxel_size_fit` | voxels | 128 | {64, 128, 256} | Hard-clamped ≤ 256 in fit loop. |
| `voxel_size_render` | voxels | 256 | {64, 128, 256, 384, 512} | One-off post-fit render at this size. |
| `box_size_A` | Å | 1000 | > 0 | Edge length of the cubic simulation box. |
| `volume_fraction` (φ) | — | 0.30 | [0.05, 0.95] | Minority-phase volume fraction. |
| `contrast` (Δρ²) | 10²⁰ cm⁻⁴ | 1.0 | ≥ 0 | Auto-derived when Link is on. |
| `link_phi_contrast` | — | True | — | Derive Δρ² from invariant. |
| `power_law_B` | cm⁻¹ Å^P | 0 | ≥ 0 | Power-law amplitude. |
| `power_law_P` | — | 4.0 | [0, 6] | Power-law exponent. |
| `background` | cm⁻¹ | 0 | ≥ 0 | Flat background. |
| `rng_seed` | — | None | — | Integer freezes randomness. |
| `n_mc_runs` | passes | 10 | [1, 500] | MC uncertainty passes. |

Each fittable parameter has a `Fit?` checkbox and `lo` / `hi` limits.
The **No limits** master toggle hides the lo/hi columns and switches the
optimiser to Nelder-Mead (no bounds).

---

## Background tabs (Power-law + Flat)

The background contributes additively:

```
I_bg(Q) = B · Q⁻ᴾ + flat
```

It is subtracted from the data before the GRF inversion, and re-added when
computing the model I(Q) for χ² evaluation. Each term has its own Fit /
limits row, so you can fix the Porod tail (`P = 4`) while letting `B` and
`flat` float, for example.

The **Power-law** tab is appropriate when your data has a meaningful Q⁻⁴
Porod tail; the **Flat** tab covers electronic / detector noise. Use both
together for the most realistic fit.

---

## 3D viewer controls

The 3D viewer renders the binary voxelgram as an **isosurface (flying-edges)**
— much faster than full volume rendering on integrated GPUs.

| Action | How |
|---|---|
| Rotate | Left-click + drag |
| Pan | Shift + left-click + drag |
| Zoom | Scroll wheel or right-click + drag |
| Reset view | Right-click → **Reset view** |
| Pick color | Right-click → **Pick isosurface color…** |
| Toggle bounding box | Right-click → **Show / Hide bounding box** |
| Save screenshot | Right-click → **Save screenshot…** (PNG/JPEG) |

If PyVista is not installed, the 3D viewer pane shows a yellow install hint
instead of crashing; the rest of the tool (I(Q) plot, 2D slice, fitting)
still works. Install with:

```bash
pip install pyirena[gui3d]
```

or, on macOS Apple Silicon if wheels are missing:

```bash
conda install -c conda-forge vtk pyvista
```

---

## 2D slice viewer

A simple pyqtgraph ImageView with a black/white lookup table (suited to
binary voxels). Controls:

- **Slice plane** combo: choose XY (Z slice), XZ (Y slice), or YZ (X slice).
- **Position slider**: drag to scrub through the 3rd axis. The label below
  shows `Slice 17 / 64 (z = 265.6 Å)`.

---

## Pop-out mode

Each viewer (2D slice and 3D PyVista) has a **Pop out ⤢** button below it.
Click to detach the viewer into a resizable, top-level dialog window — useful
for inspecting the voxelgram in detail or for comparing both viewers
side-by-side at large size. Closing the dialog returns the viewer to its
original slot in the panel.

---

## Fitting and MC uncertainty

- **Graph Model** runs `compute_voxelgram` once with the current parameters
  (no fit). Synchronous — typically 5–30 s at 128³.
- **Fit** runs `SaxsMorphEngine.fit` on a background `QThread`. Status bar
  shows progress; **Cancel** aborts at the next iteration (sub-second
  latency at ≤256³). On success, best-fit values are written back into the
  widgets and the 3D viewer refreshes at the render resolution.
- **Calc. Uncertainty (MC)** perturbs the data by Gaussian noise ~ σ_I and
  refits N times (N = `Passes:` spinbox). Results are reported as ± standard
  deviations in a popup and stored in the result.
- **Revert** restores all widget values to the snapshot taken just before
  the last fit started.

---

## Performance and memory notes

| Voxel size | RAM (uint8) | RAM (FFT, transient) | Mesh build | HDF5 size (gzip) |
|---|---|---|---|---|
| 64³ | 256 KB | ~64 MB | < 1 s | ~10–50 KB |
| 128³ | 2 MB | ~512 MB | ~1 s | ~100 KB–2 MB |
| 256³ | 16 MB | ~4 GB | ~3 s | ~1–10 MB |
| 384³ | 54 MB | ~13 GB | ~10 s | ~3–30 MB |
| 512³ | 128 MB | ~32 GB | ~30 s | ~10–100 MB |

**Recommendation**: keep `voxel_size_fit ≤ 256` (the engine enforces this).
Use 384³ / 512³ only for the final render of a converged model.

---

## HDF5 output format

Saved at `entry/saxs_morph_results/`. See
[HDF5_Structure_Reference.md](HDF5_Structure_Reference.md) for the full
on-disk inventory. The voxelgram is stored as a 3-D `uint8` dataset with
gzip compression (level 4) and chunking `(N, N, 1)` so 2D slices load
cheaply without inflating the whole cube.

---

## Limitations (v1)

- **Cube only** — anisotropic boxes (e.g. 64×256×256) are not yet supported.
- **Two phases only** — ternary thresholding is planned for a follow-up.
- **Fit voxel size hard-clamped to 256³** — the per-iteration FFT memory
  cost grows as N³·complex128, and 384³ / 512³ during fitting will OOM most
  laptops.
- **Cancel latency** — proportional to one iteration of `compute_voxelgram`,
  so cancelling at 256³ takes up to a few seconds.

---

## Scripting API

```python
from pyirena import fit_saxs_morph

result = fit_saxs_morph(
    data_file='sample.h5',
    config_file='pyirena_config.json',  # exported from GUI
    save_to_nexus=True,
    with_uncertainty=True,
    n_mc_runs=10,
)
print(result['message'])
```

Or programmatically:

```python
from pyirena.core.saxs_morph import SaxsMorphEngine, SaxsMorphConfig
from pyirena.io.nxcansas_saxs_morph import save_saxs_morph_results

cfg = SaxsMorphConfig(
    voxel_size_fit=128, voxel_size_render=256,
    box_size_A=1000.0, volume_fraction=0.30,
    link_phi_contrast=True, rng_seed=42,
)
engine = SaxsMorphEngine()
res = engine.fit(cfg, q, I, dI)          # SaxsMorphResult
save_saxs_morph_results('out.h5', res)
print(f'chi^2 = {res.chi_squared:.4g}, phi_actual = {res.phi_actual:.4g}')
```
