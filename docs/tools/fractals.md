(fractals)=
# Fractals (Mass-Fractal Aggregate Visualization)

:::{admonition} At a glance
:class: tip

**Launch:** `pyirena-gui` → **Support Tools** → **Fractals**. \
**Python:** `pyirena.core.fractals` (`grow_aggregate`, `intensity_unified`, `optimize_growth`);
no batch entry point. \
**MCP:** `pyirena_read_fractals` lists aggregates stored in a file; no growth tools. \
**Igor equivalent:** {doc}`Irena Fractals <igor:Irena/Fractals>` (mass fractal aggregate model:
{doc}`igor:Irena/3DAggregate`)
:::

The Fractals tool **grows random mass-fractal aggregates by Monte-Carlo random walk on a simple
cubic lattice**, computes their fractal parameters (Z, dmin, c, df, Rg primary, Rg aggregate),
and back-calculates the scattering intensity I(Q) so users can compare the model against an
optional Unified-fit result loaded from a NeXus file.

It is a **visualization / qualitative-modeling tool**. No data fitting — just generate aggregates
whose fractal parameters resemble those measured on real data, and inspect them in 2D and 3D.
GUI-only; no batch / scripting / headless API.

% CHECK(Jan): "GUI-only; no batch / scripting / headless API" — pyirena.core.fractals exposes
% grow_aggregate, intensity_unified and optimize_growth (used in "From Python" below), and
% pyirena_read_fractals reads saved aggregates.

(fractals-window)=
## The window

```{figure} img/fractals/fractals-01-window.png
:alt: Fractals window. Left, the controls: NeXus file box, growth parameters, Grow One / Grow Many / Optimizer tabs, active jobs, stored aggregates, active aggregate parameters, back-calculated intensity and save buttons. Right, the I(Q) plot with data, the loaded Unified fit and the aggregate's analytical curve, the 2D slice viewer and the 3D viewer pane. Numbered circles mark the parts in the legend.

The Fractals window with `modeling_mass_fractal.h5` (with a stored two-level Unified fit, levels
not RgCutoff-linked, so no targets are detected) loaded and one aggregate grown (Z = 250, seed 12345). The 3D pane shows its placeholder here because the
screenshot is made without a display; see the 3D view below.
```

The window is a horizontal splitter: a scrollable controls panel on the left; on the right the
I(Q) plot (top, ~40 %) and the 2D slice and 3D viewers (bottom, ~60 %).

1. **NeXus file (optional)** — file picker plus a coloured target summary.
2. **Growth parameters** — Z, sticking probability, # test paths, Rg primary, allowed neighbor
   distance, multi-particle attraction, random seed.
3. **Mode tabs** — Grow One / Grow Many / Optimizer.
4. **Active jobs** — running and queued growth jobs with progress.
5. **Stored aggregates** (session-only list) — every grown aggregate accumulates here; click to
   make it active.
6. **Active aggregate parameters** — Z, dmin, c, df, R, p, s, true sticking %, Rg primary, Rg
   aggregate, primary diameter, # endpoints. Shown alongside their target values when a target
   is loaded.
7. **I(Q) panel** — Q range and Monte-Carlo button.
8. **Save** — write selected aggregate to NeXus, or load aggregates from an existing NeXus.
9. **I(Q) plot** — log-log I(Q) plot showing optional loaded data, optional loaded Unified-fit
   model, the active aggregate's analytical Unified curve (green), and (if computed) the active
   aggregate's Monte-Carlo curve (blue dashed). Legend shows the four series.
10. **2D slice viewer**, with a `Pop out ⤢` button.
11. **3D PyVista voxelgram viewer**, with a `Pop out ⤢` button.

% SCREENSHOT fractals-04-3d-viewer: 3D PyVista isosurface of a grown aggregate (Z = 250, seed
% 12345) in the right-hand pane; needs a real display (VTK crashes offscreen).

(fractals-workflow)=
## Typical workflow

### Without reference data

1. Set growth parameters (defaults are good starting values).
2. Click **Grow One**. Within ~10 s a new aggregate appears in the "Stored aggregates" list.
3. Click the new entry; the I(Q) plot, 2D slice, and 3D viewer populate.
4. Note the computed fractal parameters in the "Active aggregate parameters" panel.
5. Adjust parameters and grow more aggregates to explore the morphology space. Use **Grow Many**
   to queue a batch.

### With a reference Unified-fit result

The screenshots use `validationData/modeling/modeling_mass_fractal.h5` (a Teixeira mass fractal
of R = 40 Å spheres, Dv = 2.4, ξ = 800 Å; see `validationData/README.md`) after a two-level
Unified fit was stored in it with **Store in File**. The levels in that fit are not
RgCutoff-linked, so the screenshots show no targets.

1. Click **Open NeXus…** and pick an HDF5 file with a Unified-fit result that contains a fractal
   level pair (consecutive levels where the high level's `RgCutoff` ≈ the low level's `Rg`).
2. The target fractal parameters appear in the top NeXus box AND just above the "Active aggregate
   parameters" widget. See [Target parameters from a Unified-fit
   pair](#fractals-targets) for the formulas used.

   % CHECK(Jan): which dataset and Unified fit make a good fractal reference? A linked two-level
   % fit of validationData/modeling/modeling_mass_fractal.h5 gives c = 0.67 < 1 (free B) or
   % c = 1.00 with level 1's B at its lower bound (mass-fractal option), so the panel reports
   % "not a valid fractal". Target detection itself works since 16-11-2026 (RgCutoff fix).
   % SCREENSHOT fractals-05-targets: NeXus box and "Active aggregate parameters" with target
   % values (green/orange/red) after loading a real mass-fractal aggregate with a good linked
   % two-level Unified fit.

3. Set Z in the Growth parameters to the target Z, and use the **Optimizer** tab with
   `target_dmin` and `target_c` set from the targets shown. Click **Find Best Growth**.
4. The optimizer grows several aggregates and selects the one whose computed `dmin` and `c` are
   closest to the targets. All trials accumulate in the Stored aggregates list; the best is
   bolded.

   ```{figure} img/fractals/fractals-03-optimizer.png
   :alt: Controls panel after Find Best Growth with target dmin 1.3 and target c 1.4, two iterations, Z 100, seed 12345: six trial aggregates in the Stored aggregates list, the last marked BEST and selected, its dmin and c shown in green next to the targets.

   Optimizer tab after **Find Best Growth** (targets dmin = 1.3, c = 1.4, 2 iterations, Z = 100,
   seed 12345): the best trial is marked ★ BEST and selected.
   ```

5. Compare the active aggregate's analytical Unified curve (green) to the loaded data and
   Unified-fit model (blue scatter / red line) in the I(Q) plot. Run **Compute Monte Carlo I(Q)**
   to overlay the slow but more direct PDF-based intensity (blue dashed).

   ```{figure} img/fractals/fractals-02-montecarlo-iq.png
   :alt: I(Q) plot with the data, the loaded Unified fit (red), the aggregate's analytical Unified curve (green) and the Monte-Carlo curve (blue dashed).

   The I(Q) plot after **Compute Monte Carlo I(Q)** for the aggregate of the window screenshot.
   ```

(fractals-controls)=
## Controls

(fractals-growth-parameters)=
### Growth parameters

| Control | Meaning | Default | Tips |
|---|---|---|---|
| **Degree of aggregation Z** | Number of primary particles. | 250 | Larger Z → more reliable fractal parameters but slower growth. |
| **Sticking probability [%]** | Probability per contact event that the random-walking particle sticks. | 75 | Low SP → compact (high df). High SP → loose (low df). |
| **Number of test paths** | Maximum unique paths per endpoint enumerated when computing fractal parameters. | 2500 | Affects c (and therefore dmin). |
| **Rg primary [Å]** | Radius of gyration of the primary sphere. | 10 | Sets the physical scale (`primary_diameter = 2·√(5/3)·Rg ≈ 2.58 Rg`). |
| **Allowed neighbor distance** | Which lattice neighbors count as "in contact": Edge (6 nbrs), Face diagonal (18 nbrs), Body diagonal (26 nbrs). | Body diagonal | |
| **Multi-particle attraction** | How sticking probability changes when ≥2 existing particles are within reach: Neutral / Attractive / Repulsive / Not allowed. | Neutral | |
| **Random seed (0 = random)** | Set non-zero for reproducible aggregates. | 0 (random) | |

(fractals-mode-tabs)=
### Mode tabs

**Grow One.** Click **Grow** to grow a single aggregate with the current parameters. The job
runs in a background QThread; the panel stays responsive and you can inspect previously-completed
aggregates while it runs.

**Grow Many.** Set **N aggregates** (default 5) and click **Grow N**. N growths are queued
sequentially. When `seed > 0`, each gets `seed + i` so the batch is reproducible; otherwise each
is fully random. Useful for assessing variance in fractal parameters with the same input.

**Optimizer ("Find Best Growth").** Set `Target dmin`, `Target c`, `Tolerance`, and `Max
iterations`. Click **Find Best Growth**. The optimizer bisects the sticking probability over
`[10, 90]` %, growing 3 trial aggregates per iteration (low / mid / high SP), keeping the best by
objective `(dmin − target_dmin)² + (c − target_c)²`. Stops at `max_iter` or when the objective
drops below `tolerance²`.

| Control | Meaning | Default | Tips |
|---|---|---|---|
| **Target dmin** | Target minimum dimension | 2.0 | From the NeXus targets, or typed |
| **Target c** | Target connectivity dimension | 1.2 | |
| **Tolerance** | Stop when the objective drops below `tolerance²` | 0.05 | |
| **Max iterations** | Bisection iterations, 3 trial aggregates each | 10 | |

Every trial appears in the Stored aggregates list with a label like `opt-iter-3-sp42` so you can
compare them. The best is bolded after the run completes.

(fractals-active-parameters)=
### Active aggregate parameters

When you click an aggregate in the Stored aggregates list, this widget populates with:

| Field | Meaning |
|---|---|
| Z | Actual degree of aggregation (number of particles placed) |
| dmin | Minimum dimension from path statistics |
| c | Connectivity dimension (df / dmin) |
| df | Mass-fractal dimension (log Z / log R_dimensionless) |
| R (lattice) | Weighted endpoint-to-endpoint distance, lattice units |
| p | Weighted average path length |
| s | exp(ln Z / dmin), an Irena-reported aggregate metric |
| True sticking [%] | 100·Z / total random-walk attempts |
| Rg primary [Å] | Primary-particle Rg (echoed from input) |
| Rg aggregate [Å] | sqrt(Rg_centers² + Rg_primary²) — direct from positions, NOT the McGlasson approximation |
| Primary diameter [Å] | 2·√(5/3)·Rg_primary |
| # endpoints | Particles with neighbor count < 2 (chain ends) |

When a target is loaded, each comparable field shows `actual   (target: X)` and is colour-coded:

- Green when within 10 % of the target.
- Orange when within 25 %.
- Red otherwise.

A bold compact target summary line is shown above the table so the relationship between actual
and target stays visible without scrolling back to the top of the panel.

(fractals-intensity-controls)=
### Back-calculated intensity and saving

| Control | Meaning | Default | Tips |
|---|---|---|---|
| **Q min / Q max [Å⁻¹]** | Q range of the back-calculated curves | 0.001 / 1.0 | Auto-filled with the data's Q range when a NeXus file is loaded |
| **# points** | Points in the back-calculated Q grid | 200 | |
| **Compute Monte Carlo I(Q)** | Run the slow but more direct calculation (see [Method notes](#fractals-iq)) | | **Cancel MC** stops it |
| **Save selected aggregate to NeXus…** | Append the selected aggregate to an HDF5/NeXus file | | See [Results and outputs](#fractals-results) |
| **Load aggregate(s) from NeXus…** | Pull every aggregate from an existing file into the session list | | |

(fractals-viewers)=
### 3D and 2D viewers

**3D viewer (right, PyVista isosurface)**

- Renders the voxelgram as a flying-edges isosurface at value 0.5.
- Mesh colour: dark grey by default (right-click → "Pick isosurface colour…" to change).
- White background with black axes and tick labels. Axis units are in Ångström (`X [A]` /
  `Y [A]` / `Z [A]`).
- Right-click context menu: reset view, change colour, toggle bounding box, save screenshot.
- Click `Pop out ⤢` to open the viewer in a separate dialog you can resize / maximize.

The display uses `oversample=10, sphere_voxel_radius=10` (Irena-style "fat spheres" of radius D =
primary_diameter, twice the physically correct R). This is intentional for **visualization only**
— it makes all neighbor types (edge, face diagonal, body diagonal) overlap so the aggregate looks
like a connected blob rather than scattered particles. The MC scattering calculation uses the
physically correct `oversample=20, sphere_voxel_radius=10` independently.

**2D slice viewer (left, pyqtgraph)**

- Shows one slice through the voxelgram.
- White background, dark grey solid (majority phase), white voids (minority phase). Axes black;
  units `[A]`.
- Combobox at top selects the slice plane (XY / XZ / YZ).
- Slider scrubs the slice index; a label shows the current physical coordinate.
- Pop out via the button below.

(fractals-results)=
## Results and outputs

(saving-and-loading-aggregates)=
### Save

Select an aggregate in the list, click **Save selected aggregate to NeXus…**, and pick a file
path. The aggregate is appended as a new `entry/fractals_results/aggregate_{N}` group inside an
HDF5/NeXus file. N auto-increments; multiple aggregates can coexist in the same file. Other groups
(Unified fit, raw data, other aggregates) are preserved. Results groups in general are described
in {ref}`hdf5-results-groups`.

The HDF5 layout:

```
entry/fractals_results/aggregate_{N}/   [NXprocess]
  @NX_class      = "NXprocess"
  @analysis_type = "Mass Fractal Aggregate"
  @program       = "pyirena"
  @timestamp
  positions          (Z×3 int32, gzip-compressed)
  neighbor_list      (Z×26 int32, -1 padded, gzip-compressed)
  neighbor_count     (Z, uint8)
  attempt_value      (scalar int)
  parameters/   [group, scalars]
    Z, dmin, c, df, R_dimensionless, p, s, RgPrimary,
    RgAggregate, PrimaryDiameter, TrueStickingProbability,
    NumEndpoints, NumPathsUsed
  input_params/   [group]
    StickingProbability, NumberOfTestPaths, AllowedNearDistance,
    Seed (scalars), MultiParticleAttraction (string)
  intensity/   [optional group, present if Q + I_unified present]
    Q, I_unified, I_montecarlo
```

### Load

Click **Load aggregate(s) from NeXus…** to pull every aggregate from an existing file into the
session list. Useful for resuming work on saved aggregates or for cross-comparing aggregates from
multiple files.

Aggregates accumulate in the session list in memory only. Use **Save selected aggregate to
NeXus…** to persist any you want to keep across sessions.

(fractals-scripting)=
## From Python or AIDA

::::{tab-set}
:::{tab-item} Python
The core functions behind the window can be called directly. Run from the repository root:

```python
import numpy as np

from pyirena.core.fractals import GrowthConfig, grow_aggregate, intensity_unified
from pyirena.io.nxcansas_fractals import save_fractal_aggregate

agg = grow_aggregate(GrowthConfig(z=100, sticking_prob=75.0, rg_primary=10.0, seed=12345))
p = agg.params
print(f"Z={p.z} df={p.df:.3f} dmin={p.dmin:.3f} c={p.c:.3f} Rg_agg={p.rg_aggregate:.1f} A")

agg.q = np.logspace(-3, 0, 200)
agg.i_unified = intensity_unified(p, agg.q)          # analytical Unified I(Q)
print(save_fractal_aggregate("my_aggregates.h5", agg))   # entry/fractals_results/aggregate_1
```
:::
:::{tab-item} AIDA / MCP
There are no growth tools over MCP. `pyirena_read_fractals(path)` lists the aggregates stored in a
file (one row per aggregate, with its fractal parameters), so an agent can tabulate or compare
aggregates saved from this window.
:::
::::

(fractals-method)=
## Method notes

(fractals-overview)=
### Scattering of a mass-fractal aggregate

A mass-fractal aggregate is an irregular particle made of `Z` primary spheres connected on a
simple cubic lattice. Its scattering profile shows three regimes:

1. **Low Q** (`Q · Rg_aggregate << 1`): Guinier plateau — the aggregate acts as one large object
   of size ~Rg_aggregate.
2. **Intermediate Q** (`1/Rg_aggregate << Q << 1/Rg_primary`): mass-fractal power-law with slope
   `−df`, where `df` is the mass-fractal dimension (typically 1.5–2.5 for diffusion-limited
   cluster aggregation).
3. **High Q** (`Q · Rg_primary >> 1`): primary-particle Porod region with slope `−4` for smooth
   spheres.

Fitting real data with Unified Fit gives back `(Rg_primary, Rg_aggregate, df)` and three other
fractal descriptors (`Z`, `dmin`, `c`). This tool lets you grow a candidate aggregate that
reproduces those numbers and inspect its 3D morphology — answering "what does an aggregate with
these fractal parameters actually look like?"

(target-parameters-from-a-unified-fit-pair)=
(fractals-targets)=
### Target parameters from a Unified-fit pair

The Fractals tool can read both the experimental I(Q) and the stored Unified-fit result from the
same NeXus file. When you click **Open NeXus…**, the tool:

1. Loads the I(Q) data (`Q`, `I`, `Idev`) via the standard NXcanSAS reader.
2. Loads the Unified-fit results (`G_i`, `Rg_i`, `B_i`, `P_i`, `RgCutoff_i` for every level).
3. Scans levels for the first **consecutive pair** `(low, high)` where `levels[high].RgCutoff` is
   within ±25 % of `levels[low].Rg`. This pair is the fractal representation: low = primary
   particles, high = aggregate.
4. Computes target fractal parameters (Igor formulas — see below).
5. Auto-fills the I(Q) plot's Q range with the data's Q range.

If no qualifying level pair is found, the tool falls back to "no targets" and growth still works
without comparison.

The last loaded NeXus path is saved to the application state so the next session starts with the
file pre-loaded (silently cleared if the file has been moved or deleted).

Given a fractal level pair (`level_low` = primary, `level_high` = aggregate):

```
df      = P2
gamma   = exp( lgamma(P2 / 2) )         # Γ(P2 / 2)
dmin    = B2 · Rg2^P2 / (gamma · G2)
c       = P2 / dmin
Z       = G2 / G1 + 1
```

These are exactly the formulas used by Irena (`IR3A_*` in `IR3_3DModels.ipf`). Notes:

- `df` (fractal dimension) is just the Porod exponent of the aggregate level.
- `dmin` (minimum dimension) describes the topology of paths through the aggregate. `dmin = 1` is
  a linear chain; `dmin = df` is a compact sphere. Typical: 1.1–1.4.
- `c = df / dmin` is the connectivity dimension. Typical: 1.1–1.4.
- `Z` (degree of aggregation) = number of primary particles. The `+1` accounts for the central
  particle in Beaucage's level convention.

When no Unified-fit pair is loaded, you can still set `target_dmin` and `target_c` manually in the
Optimizer tab.

(iq-back-calculation)=
(fractals-iq)=
### I(Q) back-calculation

Two paths produce I(Q) for the active aggregate:

**Analytical Unified (always).** A two-level Beaucage Unified-fit closed-form intensity:

- **Level 1** (primary sphere): `G=1, Rg=Rg_primary, P=4, B=4π/Rg_primary⁴`
- **Level 2** (aggregate): `G=Z, Rg=Rg_aggregate, P=df, B = (G·P / Rg^P) · Γ(P/2), RgCutoff =
  Rg_primary`

`Rg_aggregate` here is the **direct measurement from the grown particle positions** via the
parallel-axis (Steiner) theorem:

```
Rg_centers (Å) = sqrt(mean(||p_i − centroid||²)) × primary_diameter
Rg_aggregate   = sqrt(Rg_centers² + Rg_primary²)
```

(The McGlasson approximation `Rg_primary · Z^((1/c−1)/(dmin−df))` is not used because it
underestimates the true Rg by 30–100 % depending on the morphology, which makes the analytical
curve disagree with the MC curve.)

The analytical curve is shown in **green**.

**Monte-Carlo PDF (on demand).** Click **Compute Monte Carlo I(Q)** to run the slow but more
direct calculation:

1. Voxelize the aggregate: place a sphere of physical radius `R = primary_diameter / 2` around
   every lattice position on a 20×-oversampled cubic grid (sphere kernel = 10 voxels, voxel pitch
   = R / 10 in physical units).
2. Sample random pairs of solid voxels uniformly with replacement and build a histogram of their
   Euclidean distances. This histogram IS the Glatter pair-distance distribution
   `p(r) = 4π·r²·γ(r)`.
3. Apply the Debye / Glatter sine transform:
   ```
   I(Q) ∝ ∫₀^∞ p(r) · sinc(Qr) dr
   ```
4. Truncate above `Q_voxel = π / pitch_A` (the voxel-grid Nyquist Q, above which the discrete
   voxel approximation is unreliable).

The MC curve is shown in **blue dashed**. Costs: ~3–5 s per call, ~40 MB voxelgram for typical
Z=80 aggregates.

% CHECK(Jan): the Monte-Carlo steps above describe a voxel-pair sampling method; the code
% (core/fractals.py intensity_montecarlo) now places uniform points in every primary sphere
% (optional 10 % polydispersity) and histograms ALL pair distances (Shape2SAS-style).

**Invariant rescaling.** When a Unified-fit reference is loaded, both model curves (green and
blue) are scaled to the data using the integral `∫ I_data dQ / ∫ I_model dQ` over the **fractal
regime** `Q ∈ [0.5π/Rg_aggregate, 1.5π/Rg_primary]`. This window deliberately EXCLUDES very low Q
(sample-level power-law contamination) and very high Q (instrumental flat background) that the
single-aggregate model never tries to reproduce.

(fractals-slit-smearing)=
### Slit smearing (USAXS)

When the loaded I(Q) data are slit smeared (an NXcanSAS `dQl`), the calculated intensity is slit
smeared for the overlay/comparison so it can be judged against the data on equal footing. If
smearing the comparison curve fails it is logged and flagged rather than silently shown
unsmeared. See **[Slit smearing](../methods/slit-smearing.md)** for the full reference.

(algorithm-details-and-igor-compatibility)=
## Algorithm details and Igor compatibility

The following details are reproduced verbatim from Irena's `IR3A_*` functions to guarantee that a
pyirena aggregate matches an Irena aggregate grown with identical inputs:

- **Sticking probability table** (per contact count `chcnt`):
  - `chcnt == 1`: base SP
  - `chcnt ≥ 2` (Attractive): `(SP+100)/2` for chcnt=2, `(SP+300)/4` for ≥3
  - `chcnt ≥ 2` (Repulsive): `(SP+10)/2` for chcnt=2, `(SP+30)/4` for ≥3
  - `chcnt ≥ 2` (Not allowed): 0
  - `chcnt ≥ 2` (Neutral): base SP
- **Neighbor distance thresholds**: `1.1`, `1.05·√2`, `1.05·√3` (intentionally lenient).
- **Path-walk junction rule** (computing dmin, c): at junctions with more than 3 outgoing
  neighbors, only the first 3 are explored (matches `IR3A_MT_NextPathStep`). Removing this rule
  changes c and dmin and breaks comparison against historical Irena results.
- **Wall launch**: 6 wall faces normally; 2 corner-launch fallbacks for compact aggregates with
  Z < 6.
- **Voxel oversample 10×, sphere kernel 10 voxels** (Irena-style display) — see
  [3D and 2D viewers](#fractals-viewers).

(performance-and-memory-notes)=
## Performance and memory notes

- Growth: ~0.3 s for Z=50, ~3 s for Z=250, ~30 s for Z=1000 (single thread, default parameters).
  Body-diagonal sticking is fastest; edge-only is slowest because there are fewer contact
  opportunities.
- Path enumeration: parallelized across endpoints via `ThreadPoolExecutor` (Python threads, no
  GIL release for pure-Python walks but I/O-bound enough to help).
- Voxelization (display): ~8 MB at Z=80, lighter — uses `oversample=10`.
- Voxelization (MC scattering): ~80 MB at Z=80, ~340 MB at Z=500 — uses `oversample=20` for finer
  grid, doubling the voxel-Nyquist Q.
- Monte-Carlo intensity: 3–5 s for typical aggregates with the default budget of 1e7 pair samples
  / 20 s wall-time.

## Notes and limitations

- This is a visualization tool — there is no fitting of measured data to growth parameters (the
  optimizer just searches for matching `dmin` and `c`, both computed from the grown aggregate, not
  measured from data).
- The 3D viewer requires PyVista; on systems without it, only the 2D slice and the I(Q) plot
  work. Install with `pip install pyirena[gui3d]`.
- No batch / scripting / CLI entry point — intentionally GUI-only given the
  qualitative-visualization purpose.

## References

- Beaucage, G. (1995). *J. Appl. Cryst.* **28**, 717. Unified Fit model.
- Mountain, R. D. and Mulholland, G. W. (1988). *Langmuir* **4**, 1321. Lattice MC growth of
  fractal aggregates.
- Sorensen, C. M. and Roberts, G. C. (1997). *J. Colloid Interface Sci.* **186**, 447. Fractal
  dimension from path statistics (df, c, dmin).
- Igor Pro reference:
  [`IR3_3DModels.ipf`](https://github.com/jilavsky/SAXS_IgorCode/blob/master/User%20Procedures/Irena/IR3_3DModels.ipf)
  and
  [`IR3_3DSupportFunctions.ipf`](https://github.com/jilavsky/SAXS_IgorCode/blob/master/User%20Procedures/Irena/IR3_3DSupportFunctions.ipf).
