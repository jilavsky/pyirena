(saxs-morph)=
# SAXS Morph Tool (3D Voxelgram)

:::{admonition} At a glance
:class: tip

**Launch:** `pyirena-gui` → **Models → SAXS Morph (3D)**, or the **3D saxsMorph (GUI)** button
(it opens the selected file). Standalone: `pyirena-saxsmorph path/to/file.h5`. \
**Python:** `pyirena.batch.fit_saxs_morph(data_file, config_file)`, or
`pyirena.core.saxs_morph.SaxsMorphEngine` directly. \
**MCP:** no control tools; `pyirena_read_saxs_morph` reads stored results. \
**Igor equivalent:** {doc}`Irena Two Phase Solid model <igor:Irena/TwoPhaseSolid>`
:::

% CHECK(Jan): Igor equivalent marked uncertain (?) in irena-pyirena-map.md; the page's own Igor
% reference is IR3_3DTwoPhaseSolid.ipf, which suggests TwoPhaseSolid is right.

The SAXS Morph tool generates a **3D voxelgram of a two-phase porous structure** from your
experimental I(Q) using the Gaussian Random Fields (GRF) method, then computes the model I(Q)
from that voxelgram and refines parameters (volume fraction, contrast, Power-law + Flat
background) by fitting back to the data.

The result is a 3D phase map you can rotate, slice, and inspect — useful for intuitive
visualisation of the *kind* of microstructure consistent with a SAS profile.

% CHECK(Jan): "refines parameters ... by fitting back to the data" conflicts with Method notes
% ("not an iterative fit"; the only fits are the two background pre-fits).

(saxs-morph-window)=
## The window

```{figure} img/saxs-morph/saxs-morph-01-window.png
:alt: SAXS Morph window after Calculate 3D on simple_debye_bueche.h5. Controls on the left, the I(Q) plot with data, red model and the Qmin and Qmax cursors at top right, the 2D slice of the voxelgram below it, and the result block with chi-squared and morphology metrics. The 3D pane is covered by a grey box. Numbered circles mark the parts described in the legend.

The SAXS Morph window after **Calculate 3D** on `validationData/simple_fits/simple_debye_bueche.h5`.
The 3D view (11) cannot be rendered by the screenshot script and is shown as a grey box; see
the 3D viewer screenshot under [3D viewer controls](#saxs-morph-3d-viewer).
```

1. **Data file**: the loaded file; **Open…** picks another.
2. **Q range**: the modelling Q range, set by the Qmin (red) and Qmax (blue) cursors.
3. **Voxel grid**: cube side, box size, RNG seed and smoothing.
4. **Two-phase parameters**: input mode, volume fraction φ and contrast (Δρ)².
5. **Background tabs**: **Power-law Bckg** and **Flat Bckg** pre-fits.
6. **Calculate 3D**.
7. **Result block**: χ², φ, voxel pitch, S/V and the morphology metrics of the minority phase.
8. **Save/Append to HDF5…**, **Load Setup from File…**, **Save config to JSON…**.
9. **I(Q) plot**: data, data − background, the red model, the cursors and the Q_box / Q_nyq
   markers.
10. **2D slice viewer**.
11. **3D viewer** (PyVista isosurface).
12. **Status bar**.

(saxs-morph-workflow)=
## Typical workflow

The steps use `validationData/simple_fits/simple_debye_bueche.h5`, a Debye-Bueche curve (random
two-phase structure, correlation length 80 Å) with no background. Because the file has no
background, steps 2 and 3 are skipped for it (B = 0 and background = 0 are the defaults); they are
needed for real data.

The workflow has **three sequential pre-steps** before generating the 3D voxelgram:

1. **Load a data file** (Data Selector → SAXS Morph (GUI), or `Open…` button in the panel).
2. **Pre-fit the Power-law background** (low-Q tail):
   - Drag the red/blue cursors to a low-Q window where the data is dominated by the Q⁻⁴ Porod
     slope.
   - Open the **Power-law Bckg** tab → click **Set from cursors** → the Q-range fields fill in.
   - Click **Fit Power-law Bckg** → B and P are populated.
3. **Pre-fit the flat background** (high-Q noise floor):
   - Drag the cursors to a high-Q window where structural scattering has decayed and the residual
     is detector noise.
   - Open the **Flat Bckg** tab → click **Set from cursors** → click **Fit Flat Bckg** →
     background is populated (median of I − Power-law in that range).
4. **Set the modelling Q range**: drag the cursors back to the Q range you want the GRF method to
   use (typically the full structural region). For the example file: Q = 0.005–0.2 Å⁻¹ with a
   4000 Å box (Calculate 3D asks for confirmation when Q_min < 3 × 2π/box size).
5. **Set parameters** in the Two-phase parameters box (see [Controls](#saxs-morph-controls)):
   - **Input mode** combo:
     - *Input φ → derive Δρ²*  — recommended default. You enter φ, the contrast is computed
       automatically from the data invariant.
     - *Input Δρ² → derive φ* — useful when the contrast is known from sample composition.
     - *Use both as-is* — manual override, no derivation.
   - **Voxel cube side (render)**: 256³ is a good default; 384³ / 512³ produce nicer renders at
     the cost of disk space and RAM.
   - **Box size (Å)**: ensure `box_size_A / N` is roughly the correlation length you expect to
     resolve (e.g. 1000 Å / 256 ≈ 4 Å pitch).
   - **RNG seed**: leave blank for stochastic exploration; set an integer to freeze the voxelgram
     for reproducible figures.

   For the example: φ = 0.3, box 4000 Å, cube 256³, seed 42.
6. Click **Calculate 3D**. After ~5–30 s (depending on voxel size), you see:
   - The red **model I(Q)** overlaying the data.
   - The black **data − background** trace (what the GRF method actually models).
   - A 2D slice in the bottom-left viewer (drag the slider to scrub depth).
   - A 3D isosurface in the bottom-right PyVista viewer (rotate with mouse).
   - The result block updates with χ², φ_actual, contrast, etc.

   ```{figure} img/saxs-morph/saxs-morph-02-calculated.png
   :alt: I(Q) plot after Calculate 3D: the red GRF model follows the Debye-Bueche data from Qmin 0.005 to Qmax 0.2, with the Q_box and Q_nyq markers.

   The model I(Q) on `simple_debye_bueche.h5` (φ = 0.3). The contrast derived from the invariant
   is 4.72 × 10²⁰ cm⁻⁴; the Debye-Bueche invariant Prefactor·η²·π/4 gives 4.74.
   ```

   ```{figure} img/saxs-morph/saxs-morph-03-slice.png
   :alt: 2D slice through the middle of the 256-cubed voxelgram, XY plane, showing a bicontinuous pattern of dark minority phase (30 percent) in a white matrix, box 4000 by 4000 angstrom.

   The middle XY slice of the voxelgram (slice 129 of 256).
   ```

   % CHECK(Jan): step 6 says "black data − background trace"; in the screenshot the data −
   % background trace is green.

7. Click **Save Result to HDF5…** to write the compressed voxelgram + parameters into the source
   HDF5 file (or a new file).

   % CHECK(Jan): the button is now labelled "Save/Append to HDF5…".

If you want to try a different RNG seed, change φ, swap input modes, etc., just edit and click
**Calculate 3D** again — the voxelgram regenerates in seconds (no iterative fitting).

(saxs-morph-controls)=
## Controls

### Voxel grid and two-phase parameters

| Control | Meaning | Default | Tips |
|---|---|---|---|
| **Cube side** | `voxel_size_render` [voxels]: {64, 128, 256, 384, 512}. One-off render at this size. | 256³ | Voxel pitch = Box size / Cube side |
| **Box size [Å]** | `box_size_A` [Å]: edge length of the cubic simulation box. | 1000 | Pitch smaller than the correlation length of interest; box at least a few correlation lengths across |
| **RNG seed** | `rng_seed`: integer freezes randomness across runs. | blank (None) | Blank = a fresh random realisation each time |
| **Smoothing σ [vox]** | Gaussian smoothing applied to the voxelgram after thresholding (in voxel units). | 1.0 | 0 = binary output; 1.0 matches Igor's gauss3d step; > 2 may wash out fine features |
| **Input mode** | `input_mode`: `phi` (derive Δρ²) / `contrast` (derive φ) / `both` (no derivation). | `phi` | φ and Δρ² are linked by Q* = 2π² φ(1−φ) Δρ² |
| **Volume fraction φ** | `volume_fraction` (φ): minority-phase volume fraction. | 0.30 | |
| **Contrast (Δρ)²** | `contrast` (Δρ²) [10²⁰ cm⁻⁴]: auto-derived in mode `phi`; user-input in mode `contrast`/`both`. | 1.0 | The derived field is read-only |
| **Q range** (cursors) | `q_min/q_max` [Å⁻¹]: Q range used by **Calculate 3D**. | (cursor-driven) | |
| **Calculate 3D** | Generate the voxelgram from the current parameters. | | Run the two background pre-fits first |

There are **no fittable parameters in this workflow** beyond the two background pre-fits — the GRF
method computes the voxelgram deterministically from the chosen φ and the data spectrum.

`SaxsMorphConfig` used to carry `voxel_size_fit`, `fit_*`, `*_limits`, `no_limits` and
`n_mc_runs` for a deprecated `Engine.fit()`. Both the method and the fields are gone (GitHub issue
#29); `link_phi_contrast` is the one survivor, kept because it is written into saved HDF5 and
older files still read back through it.

(saxs-morph-background)=
### Background tabs (Power-law + Flat)

| Control | Meaning | Default | Tips |
|---|---|---|---|
| **Q min / Q max** (Power-law Bckg) | `power_law_q_min/q_max` [Å⁻¹]: Q range for the Power-law pre-fit. | None | |
| **Set from cursors** | Copy the current main I(Q) cursor positions into the Q-range fields. | | |
| **Fit Power-law Bckg** | Linear least-squares fit of log10(I) vs log10(Q) over the Q range; sets B and P. | | |
| **Amplitude B** | `power_law_B` [cm⁻¹ Å^P]: set by **Fit Power-law Bckg**. | 0 | |
| **Exponent P** | `power_law_P`: set by **Fit Power-law Bckg**. | 4.0 | |
| **Q min / Q max** (Flat Bckg) | `background_q_min/q_max` [Å⁻¹]: Q range for the Flat pre-fit. | None | |
| **Fit Flat Bckg** | Median of (I − Power-law) over the Q range. | | Run the Power-law fit first so it can be subtracted |
| **Background** | `background` [cm⁻¹]: set by **Fit Flat Bckg**. | 0 | |

The background contributes additively:

```
I_bg(Q) = B · Q⁻ᴾ + flat
```

It is subtracted from the data before the GRF inversion, and re-added when computing the model
I(Q) for χ² evaluation. Each term has its own Fit / limits row, so you can fix the Porod tail
(`P = 4`) while letting `B` and `flat` float, for example.

% CHECK(Jan): the current panel has no Fit / limits rows for B, P or background (the code says
% they are not fittable in this workflow).

The **Power-law** tab is appropriate when your data has a meaningful Q⁻⁴ Porod tail; the **Flat**
tab covers electronic / detector noise. Use both together for the most realistic fit.

(saxs-morph-calculate)=
### Calculate 3D

The single big green action button. Runs `SaxsMorphEngine.compute_voxelgram` synchronously at the
render resolution, applying the input-mode rule (derive φ from invariant + contrast, or vice
versa, or use both as-is).

The button blocks the GUI for the duration of the calculation (typically 5–30 s at 256³). For
faster iteration on slow machines or large boxes, either:

- Lower the **Cube side (render)** value temporarily.
- Use the scripting API to run a parameter sweep headlessly.

% CHECK(Jan): the current panel runs Calculate 3D on a background thread (_CalcWorker); the GUI
% stays responsive and only the button is disabled. Same in Limitations below.

There is no MC uncertainty in this workflow because the model has no fittable parameters —
variation across RNG seeds reflects only the stochastic nature of the GRF realisation, not
parameter uncertainty.

(saxs-morph-3d-viewer)=
### 3D viewer controls

The 3D viewer renders the binary voxelgram as an **isosurface (flying-edges)** — much faster than
full volume rendering on integrated GPUs.

% SCREENSHOT saxs-morph-04-3d-view: 3D isosurface of the simple_debye_bueche.h5 voxelgram
% (phi 0.3, box 4000 A, 256^3, seed 42), rotated to show the bicontinuous structure, with the
% bounding box. Manual: PyVista/VTK cannot render under the offscreen screenshot script.

| Action | How |
|---|---|
| Rotate | Left-click + drag |
| Pan | Shift + left-click + drag |
| Zoom | Scroll wheel or right-click + drag |
| Reset view | Right-click → **Reset view** |
| Pick color | Right-click → **Pick isosurface color…** |
| Toggle bounding box | Right-click → **Show / Hide bounding box** |
| Save screenshot | Right-click → **Save screenshot…** (PNG/JPEG) |

If PyVista is not installed, the 3D viewer pane shows a yellow install hint instead of crashing;
the rest of the tool (I(Q) plot, 2D slice, fitting) still works. Install with:

```bash
pip install pyirena[gui3d]
```

or, on macOS Apple Silicon if wheels are missing:

```bash
conda install -c conda-forge vtk pyvista
```

(saxs-morph-slice-viewer)=
### 2D slice viewer

A simple pyqtgraph ImageView with a black/white lookup table (suited to binary voxels). Controls:

- **Slice plane** combo: choose XY (Z slice), XZ (Y slice), or YZ (X slice).
- **Position slider**: drag to scrub through the 3rd axis. The label below shows
  `Slice 17 / 64 (z = 265.6 Å)`.

(saxs-morph-popout)=
### Pop-out mode

Each viewer (2D slice and 3D PyVista) has a **Pop out ⤢** button below it. Click to detach the
viewer into a resizable, top-level dialog window — useful for inspecting the voxelgram in detail or
for comparing both viewers side-by-side at large size. Closing the dialog returns the viewer to
its original slot in the panel.

(saxs-morph-results)=
## Results and outputs

| Control | What it does |
|---|---|
| **Save/Append to HDF5…** | Append the SAXS Morph results (voxelgram + parameters) into an HDF5 file. All other result groups in the file are preserved; only the `saxs_morph_results` group is replaced. The dialog suggests the source file. |
| **Load Setup from File…** | Restore every control from the setup embedded in a SAXS Morph result file (written by Save/Append to HDF5… or by a batch run). |
| **Save config to JSON…** | Append the current SAXS Morph parameters to a `pyirena_config.json` file; other tool sections are preserved. This is the file the Data Selector "SAXS Morph (script)" button reads. |

Saved at `entry/saxs_morph_results/`. See {ref}`hdf5-results-groups` for the results groups and
the embedded setup. The voxelgram is stored as a 3-D `uint8` dataset with gzip compression (level
4) and chunking `(N, N, 1)` so 2D slices load cheaply without inflating the whole cube.

`pyirena_read_saxs_morph(path)` returns the stored parameters and morphology metrics to an AI
agent (not the voxelgram itself).

(saxs-morph-scripting)=
## From Python or AIDA

::::{tab-set}
:::{tab-item} Python
The headless `fit_saxs_morph()` function follows the same three-step workflow as the GUI:
Power-law pre-fit → Flat pre-fit → Calculate 3D. Run from the repository root:

```python
import json

from pyirena import fit_saxs_morph

# A minimal 'saxs_morph' section; "Save config to JSON…" in the GUI writes the full one.
with open("pyirena_config.json", "w") as f:
    json.dump({"_pyirena_config": {"file_type": "pyIrena Configuration File"},
               "saxs_morph": {"q_min": 0.005, "q_max": 0.2, "box_size_A": 4000,
                              "voxel_size_render": 256, "input_mode": "phi",
                              "volume_fraction": 0.3, "rng_seed": 42}}, f, indent=2)

result = fit_saxs_morph(
    data_file='validationData/simple_fits/simple_debye_bueche.h5',
    config_file='pyirena_config.json',  # exported from GUI
    save_to_nexus=False,                # True appends the result to the data file
)
print(result['message'])
```

The `pyirena_config.json` should contain a `'saxs_morph'` section with the four pre-fit Q ranges
plus the modelling Q range. If a pre-fit Q range is missing, that pre-fit step is skipped and the
corresponding parameter values from the config are used as-is.

Or call the engine + helpers directly:

```python
from pyirena.core.saxs_morph import (
    SaxsMorphEngine, SaxsMorphConfig,
    fit_power_law_bg, fit_flat_bg,
)
from pyirena.io import discover_scattering, load_scattering
from pyirena.io.nxcansas_saxs_morph import save_saxs_morph_results

data = load_scattering(discover_scattering(
    'validationData/simple_fits/simple_debye_bueche.h5')[0])
q, I, dI = data.q, data.intensity, data.uncertainty

# This file has no background, so steps 1 and 2 are skipped here (B = 0, flat = 0).
# For data with a low-Q power law and a flat background:
# Step 1: Power-law pre-fit at low Q.
# B, P = fit_power_law_bg(q, I, q_min=1e-3, q_max=5e-3)
# Step 2: Flat-bg pre-fit at high Q (subtracts the power law first).
# flat = fit_flat_bg(q, I, q_min=0.4, q_max=0.5,
#                    power_law_B=B, power_law_P=P)
B, P, flat = 0.0, 4.0, 0.0

# Step 3: Calculate 3D.
cfg = SaxsMorphConfig(
    q_min=0.005, q_max=0.2,
    voxel_size_render=256,
    box_size_A=4000.0,
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
:::
:::{tab-item} AIDA / MCP
SAXS Morph has no control (dispatcher) tools, so an agent cannot run it. It can read a stored
result with `pyirena_read_saxs_morph(path)`: parameters and morphology metrics (pore-size
statistics, Euler number, open/closed porosity), not the voxelgram. See the
[AI tools reference](../reference/ai-tools-reference.md).
:::
::::

(saxs-morph-method)=
## Method notes

The SAXS Morph method is **not an iterative fit** of model parameters. Given the data invariant Q*
and one of (φ, Δρ²), the *other* is uniquely determined (Q* = 2π² φ(1−φ) Δρ²). Once those are
set, the spectral function F(k) is derived from the data, and the voxelgram is generated
**deterministically** from F(k) plus a random seed. No model parameter is varied to minimise χ².
The only true fits are the two **background pre-fits**.

Algorithm pipeline for one **Calculate 3D** click:

1. Subtract the user-supplied background (Power-law + Flat) from the data → `I_corr(Q)`.
2. Resolve (φ, Δρ²) according to **Input mode**:
   - mode `phi`: user supplies φ; Δρ² is computed from the invariant.
   - mode `contrast`: user supplies Δρ²; φ is computed from the invariant.
   - mode `both`: both are used as supplied (no derivation).
3. Sinc-transform `Q² · I_corr(Q)` → Debye autocorrelation γ(r); normalise γ(0) = 1.
4. Sinc-transform γ(r) · r² → spectral density F(k); clip negative bins.
5. Threshold parameter α = √2 · erfinv(1 − 2φ).
6. Generate 3D Gaussian white noise → FFT → multiply by √F(|k₃ᴅ|) → inverse FFT → real scalar
   field.
7. Threshold the scalar field at α·σ → binary uint8 voxelgram.
8. FFT the (mean-subtracted) voxelgram → spherically average |F(k)|² → resample onto Q grid →
   multiply by contrast.
9. Add background back and compute χ² against data (as a quality metric, not a fit objective).

The two background pre-fits performed before step 1:

- **Power-law**: linear least-squares fit of log10(I) = log10(B) − P·log10(Q) over a user-selected
  low-Q window where the data is dominated by the Porod tail.
- **Flat background**: median of (I − Power-law) over a user-selected high-Q window where the
  structural scattering has decayed.

% CHECK(Jan): "low-Q window where the data is dominated by the Porod tail" (also workflow step 2:
% "low-Q window ... dominated by the Q⁻⁴ Porod slope") — a Porod tail is normally a high-Q feature;
% the low-Q power law is usually a different (larger-scale) contribution.

A comparison with other 3D-morphology approaches is in
[SAXS Morph method comparison](../methods/saxs-morph-method-comparison.md).

(saxs-morph-performance)=
## Performance and memory notes

| Voxel size | RAM (uint8) | RAM (FFT, transient) | Mesh build | HDF5 size (gzip) |
|---|---|---|---|---|
| 64³ | 256 KB | ~64 MB | < 1 s | ~10–50 KB |
| 128³ | 2 MB | ~512 MB | ~1 s | ~100 KB–2 MB |
| 256³ | 16 MB | ~4 GB | ~3 s | ~1–10 MB |
| 384³ | 54 MB | ~13 GB | ~10 s | ~3–30 MB |
| 512³ | 128 MB | ~32 GB | ~30 s | ~10–100 MB |

% CHECK(Jan): the Cube side tooltip gives different numbers (e.g. 256³ ~15 s, 125 MB RAM; 512³
% ~120 s, 1 GB). The 256³ Calculate 3D of the example took about 4 s.

**Recommendation**: work at 256³ or below while you are choosing parameters. Use 384³ / 512³ only
for the final render, once the two background pre-fits and φ are settled.

(saxs-morph-limitations)=
## Limitations (v1)

- **Cube only** — anisotropic boxes (e.g. 64×256×256) are not yet supported.
- **Two phases only** — ternary thresholding is planned for a follow-up.
- **Calculate 3D blocks the GUI** — the run is synchronous (typically 5–30 s at 256³). Use a
  smaller render size for fast iteration.

## References

- Berk, N. F. (1991). *Phys. Rev. E* **51**, 4141.
- Roberts, A. P. (1997). *Phys. Rev. E* **55**, R1286.
- Levitz, P. (2007). *Modelling Simul. Mater. Sci. Eng.* **15**, S2.
- Cherny, A. Yu., et al. (2013). *J. Appl. Cryst.* **46**, 365.
- Igor Pro reference:
  [`IR3_3DTwoPhaseSolid.ipf`](https://github.com/jilavsky/SAXS_IgorCode/blob/master/User%20Procedures/Irena/IR3_3DTwoPhaseSolid.ipf)

% CHECK(Jan): "Berk, N. F. (1991). Phys. Rev. E 51, 4141" — Phys. Rev. E volume 51 is 1995
% (Berk's 1991 paper is Phys. Rev. A 44, 5069); year or journal looks inconsistent.
