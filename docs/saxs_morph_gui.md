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
10. [Calculate 3D](#calculate-3d)
11. [Performance and memory notes](#performance-and-memory-notes)
12. [HDF5 output format](#hdf5-output-format)
13. [Limitations (v1)](#limitations-v1)
14. [Scripting API](#scripting-api)

---

## Overview

The SAXS Morph method is **not an iterative fit** of model parameters.
Given the data invariant Q* and one of (φ, Δρ²), the *other* is uniquely
determined (Q* = 2π² φ(1−φ) Δρ²). Once those are set, the spectral
function F(k) is derived from the data, and the voxelgram is generated
**deterministically** from F(k) plus a random seed. No model parameter
is varied to minimise χ². The only true fits are the two **background
pre-fits**.

Algorithm pipeline for one **Calculate 3D** click:

1. Subtract the user-supplied background (Power-law + Flat) from the
   data → `I_corr(Q)`.
2. Resolve (φ, Δρ²) according to **Input mode**:
   - mode `phi`: user supplies φ; Δρ² is computed from the invariant.
   - mode `contrast`: user supplies Δρ²; φ is computed from the invariant.
   - mode `both`: both are used as supplied (no derivation).
3. Sinc-transform `Q² · I_corr(Q)` → Debye autocorrelation γ(r);
   normalise γ(0) = 1.
4. Sinc-transform γ(r) · r² → spectral density F(k); clip negative bins.
5. Threshold parameter α = √2 · erfinv(1 − 2φ).
6. Generate 3D Gaussian white noise → FFT → multiply by √F(|k₃ᴅ|) →
   inverse FFT → real scalar field.
7. Threshold the scalar field at α·σ → binary uint8 voxelgram.
8. FFT the (mean-subtracted) voxelgram → spherically average |F(k)|² →
   resample onto Q grid → multiply by contrast.
9. Add background back and compute χ² against data (as a quality metric,
   not a fit objective).

The two background pre-fits performed before step 1:

- **Power-law**: linear least-squares fit of log10(I) = log10(B) − P·log10(Q)
  over a user-selected low-Q window where the data is dominated by the
  Porod tail.
- **Flat background**: median of (I − Power-law) over a user-selected
  high-Q window where the structural scattering has decayed.

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
┌─────────────────────────────────┬────────────────────────────────────────────┐
│  Left controls (~440 px)         │  Right graphs / 3D                         │
│                                  │                                            │
│  SAXS Morph     [? Help]         │  ┌─────────────────────────────────────┐  │
│  Data file:  ………   [Open…]      │  │  I(Q)  log-log                      │  │
│  Q range: 0.001 to 0.3 Å⁻¹       │  │   • Data scatter                     │  │
│  ┌─ Voxel grid ──────────┐       │  │   • Data − background scatter        │  │
│  │ Cube side (render): 256│      │  │   • Red model line                   │  │
│  │ Box size [Å]:      1000│      │  │   • Two cursors (red Qmin / blue)    │  │
│  │ RNG seed:    (blank=…) │      │  └─────────────────────────────────────┘  │
│  └────────────────────────┘      │  ┌────────────────┬───────────────────┐  │
│  ┌─ Two-phase ─────────────┐     │  │ 2D slice       │ 3D PyVista        │  │
│  │ Input mode: [φ → derive Δρ²]│ │  │ XY/XZ/YZ combo │ rotate w/ mouse   │  │
│  │ φ:    [0.30]  (input)       │ │  │ slider         │ right-click menu  │  │
│  │ Δρ²:  [3.21]  (derived)     │ │  │  [Pop out ⤢]   │  [Pop out ⤢]      │  │
│  └─────────────────────────────┘ │  └────────────────┴───────────────────┘  │
│  ┌ Background tabs ──────────┐   │  Status bar                                │
│  │ Power-law Bckg │ Flat Bckg│   │  Calculate 3D done. χ² = 12.34            │
│  │  Q range: [0.001..0.005]  │   │                                            │
│  │  [Set from cursors]        │   │                                            │
│  │  [Fit Power-law Bckg]      │   │                                            │
│  │  B = 1.2e-5,  P = 4.0      │   │                                            │
│  └────────────────────────────┘   │                                            │
│  [Calculate 3D]  (big green)     │                                            │
│  Result block (chi², φ, …)       │                                            │
│  [Save Result to HDF5…]          │                                            │
└─────────────────────────────────┴────────────────────────────────────────────┘
```

---

## Workflow

The workflow has **three sequential pre-steps** before generating the 3D voxelgram:

1. **Load a data file** (Data Selector → SAXS Morph (GUI), or `Open…` button
   in the panel).
2. **Pre-fit the Power-law background** (low-Q tail):
   - Drag the red/blue cursors to a low-Q window where the data is dominated
     by the Q⁻⁴ Porod slope.
   - Open the **Power-law Bckg** tab → click **Set from cursors** → the
     Q-range fields fill in.
   - Click **Fit Power-law Bckg** → B and P are populated.
3. **Pre-fit the flat background** (high-Q noise floor):
   - Drag the cursors to a high-Q window where structural scattering has
     decayed and the residual is detector noise.
   - Open the **Flat Bckg** tab → click **Set from cursors** → click
     **Fit Flat Bckg** → background is populated (median of I − Power-law
     in that range).
4. **Set the modelling Q range**: drag the cursors back to the Q range you
   want the GRF method to use (typically the full structural region).
5. **Set parameters** in the Two-phase parameters box:
   - **Input mode** combo:
     - *Input φ → derive Δρ²*  — recommended default. You enter φ, the
       contrast is computed automatically from the data invariant.
     - *Input Δρ² → derive φ* — useful when the contrast is known from
       sample composition.
     - *Use both as-is* — manual override, no derivation.
   - **Voxel cube side (render)**: 256³ is a good default; 384³ / 512³
     produce nicer renders at the cost of disk space and RAM.
   - **Box size (Å)**: ensure `box_size_A / N` is roughly the correlation
     length you expect to resolve (e.g. 1000 Å / 256 ≈ 4 Å pitch).
   - **RNG seed**: leave blank for stochastic exploration; set an integer
     to freeze the voxelgram for reproducible figures.
6. Click **Calculate 3D**. After ~5–30 s (depending on voxel size), you see:
   - The red **model I(Q)** overlaying the data.
   - The black **data − background** trace (what the GRF method actually
     models).
   - A 2D slice in the bottom-left viewer (drag the slider to scrub depth).
   - A 3D isosurface in the bottom-right PyVista viewer (rotate with mouse).
   - The result block updates with χ², φ_actual, contrast, etc.
7. Click **Save Result to HDF5…** to write the compressed voxelgram +
   parameters into the source HDF5 file (or a new file).

If you want to try a different RNG seed, change φ, swap input modes, etc.,
just edit and click **Calculate 3D** again — the voxelgram regenerates in
seconds (no iterative fitting).

---

## Parameter reference

| Parameter | Units | Default | Notes |
|---|---|---|---|
| `voxel_size_render` | voxels | 256 | {64, 128, 256, 384, 512}. One-off render at this size. |
| `box_size_A` | Å | 1000 | Edge length of the cubic simulation box. |
| `input_mode` | — | `phi` | `phi` (derive Δρ²) / `contrast` (derive φ) / `both` (no derivation). |
| `volume_fraction` (φ) | — | 0.30 | Minority-phase volume fraction. |
| `contrast` (Δρ²) | 10²⁰ cm⁻⁴ | 1.0 | Auto-derived in mode `phi`; user-input in mode `contrast`/`both`. |
| `power_law_B` | cm⁻¹ Å^P | 0 | Set by **Fit Power-law Bckg**. |
| `power_law_P` | — | 4.0 | Set by **Fit Power-law Bckg**. |
| `background` | cm⁻¹ | 0 | Set by **Fit Flat Bckg**. |
| `power_law_q_min/q_max` | Å⁻¹ | None | Q range for the Power-law pre-fit. |
| `background_q_min/q_max` | Å⁻¹ | None | Q range for the Flat pre-fit. |
| `q_min/q_max` | Å⁻¹ | (cursor-driven) | Q range used by **Calculate 3D**. |
| `rng_seed` | — | None | Integer freezes randomness across runs. |

There are **no fittable parameters in this workflow** beyond the two
background pre-fits — the GRF method computes the voxelgram
deterministically from the chosen φ and the data spectrum.

The `voxel_size_fit`, `fit_*`, `*_limits`, `no_limits`, `n_mc_runs`
fields in `SaxsMorphConfig` are kept for backward compatibility with
the deprecated `Engine.fit()` method but are not exposed in the GUI.

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

## Calculate 3D

The single big green action button. Runs `SaxsMorphEngine.compute_voxelgram`
synchronously at the render resolution, applying the input-mode rule
(derive φ from invariant + contrast, or vice versa, or use both as-is).

The button blocks the GUI for the duration of the calculation (typically
5–30 s at 256³). For faster iteration on slow machines or large boxes,
either:

- Lower the **Cube side (render)** value temporarily.
- Use the scripting API to run a parameter sweep headlessly.

There is no MC uncertainty in this workflow because the model has no
fittable parameters — variation across RNG seeds reflects only the
stochastic nature of the GRF realisation, not parameter uncertainty.

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
- **Calculate 3D blocks the GUI** — the run is synchronous (typically
  5–30 s at 256³). Use a smaller render size for fast iteration.

---

## Scripting API

The headless `fit_saxs_morph()` function follows the same three-step
workflow as the GUI: Power-law pre-fit → Flat pre-fit → Calculate 3D.

```python
from pyirena import fit_saxs_morph

result = fit_saxs_morph(
    data_file='sample.h5',
    config_file='pyirena_config.json',  # exported from GUI
    save_to_nexus=True,
)
print(result['message'])
```

The `pyirena_config.json` should contain a `'saxs_morph'` section with
the four pre-fit Q ranges plus the modelling Q range. If a pre-fit Q
range is missing, that pre-fit step is skipped and the corresponding
parameter values from the config are used as-is.

Or call the engine + helpers directly:

```python
from pyirena.core.saxs_morph import (
    SaxsMorphEngine, SaxsMorphConfig,
    fit_power_law_bg, fit_flat_bg,
)
from pyirena.io.nxcansas_saxs_morph import save_saxs_morph_results

# Step 1: Power-law pre-fit at low Q.
B, P = fit_power_law_bg(q, I, q_min=1e-3, q_max=5e-3)

# Step 2: Flat-bg pre-fit at high Q (subtracts the power law first).
flat = fit_flat_bg(q, I, q_min=0.4, q_max=0.5,
                   power_law_B=B, power_law_P=P)

# Step 3: Calculate 3D.
cfg = SaxsMorphConfig(
    voxel_size_render=256,
    box_size_A=1000.0,
    input_mode='phi',
    volume_fraction=0.30,
    power_law_B=B, power_law_P=P, background=flat,
    rng_seed=42,
)
engine = SaxsMorphEngine()
res = engine.compute_voxelgram(cfg, q, I, dI,
                               voxel_size_override=cfg.voxel_size_render)
save_saxs_morph_results('out.h5', res)
print(f'phi_actual = {res.phi_actual:.4g}, '
      f'derived contrast = {res.config.contrast:.4g}')
```
