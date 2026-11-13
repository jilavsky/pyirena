# Fractals Tool (Mass-Fractal Aggregate Visualization)

The Fractals tool **grows random mass-fractal aggregates by Monte-Carlo
random walk on a simple cubic lattice**, computes their fractal parameters
(Z, dmin, c, df, Rg primary, Rg aggregate), and back-calculates the
scattering intensity I(Q) so users can compare the model against an
optional Unified-fit result loaded from a NeXus file.

It is a **visualization / qualitative-modeling tool**.  No data fitting —
just generate aggregates whose fractal parameters resemble those measured
on real data, and inspect them in 2D and 3D.  GUI-only; no batch /
scripting / headless API.

References:

- Beaucage, G. (1995). *J. Appl. Cryst.* **28**, 717.  Unified Fit model.
- Mountain, R. D. and Mulholland, G. W. (1988). *Langmuir* **4**, 1321.
  Lattice MC growth of fractal aggregates.
- Sorensen, C. M. and Roberts, G. C. (1997). *J. Colloid Interface Sci.*
  **186**, 447.  Fractal dimension from path statistics (df, c, dmin).
- Igor Pro reference: [`IR3_3DModels.ipf`](https://github.com/jilavsky/SAXS_IgorCode/blob/master/User%20Procedures/Irena/IR3_3DModels.ipf)
  and [`IR3_3DSupportFunctions.ipf`](https://github.com/jilavsky/SAXS_IgorCode/blob/master/User%20Procedures/Irena/IR3_3DSupportFunctions.ipf).

---

## Contents

1. [Overview](#overview)
2. [Launching](#launching)
3. [Panel layout](#panel-layout)
4. [Workflow](#workflow)
5. [Loading a NeXus reference file](#loading-a-nexus-reference-file)
6. [Target parameters from a Unified-fit pair](#target-parameters-from-a-unified-fit-pair)
7. [Growth parameter reference](#growth-parameter-reference)
8. [Mode tabs (Grow One / Grow Many / Optimizer)](#mode-tabs)
9. [Active aggregate parameters](#active-aggregate-parameters)
10. [I(Q) back-calculation (Unified analytical + Monte-Carlo)](#iq-back-calculation)
11. [3D and 2D viewers](#3d-and-2d-viewers)
12. [Saving and loading aggregates (NeXus)](#saving-and-loading-aggregates)
13. [Algorithm details and Igor compatibility](#algorithm-details-and-igor-compatibility)
14. [Performance and memory notes](#performance-and-memory-notes)

---

## Overview

A mass-fractal aggregate is an irregular particle made of `Z` primary
spheres connected on a simple cubic lattice.  Its scattering profile
shows three regimes:

1. **Low Q** (`Q · Rg_aggregate << 1`): Guinier plateau — the aggregate
   acts as one large object of size ~Rg_aggregate.
2. **Intermediate Q** (`1/Rg_aggregate << Q << 1/Rg_primary`): mass-fractal
   power-law with slope `−df`, where `df` is the mass-fractal dimension
   (typically 1.5–2.5 for diffusion-limited cluster aggregation).
3. **High Q** (`Q · Rg_primary >> 1`): primary-particle Porod region with
   slope `−4` for smooth spheres.

Fitting real data with Unified Fit gives back `(Rg_primary, Rg_aggregate, df)`
and three other fractal descriptors (`Z`, `dmin`, `c`).  This tool lets
you grow a candidate aggregate that reproduces those numbers and inspect
its 3D morphology — answering "what does an aggregate with these fractal
parameters actually look like?"

---

## Launching

From the Data Selector → **Support Tools** section → click **Fractals**.

The window opens with no aggregate.  Optionally load a NeXus file with
Unified-fit results to display the target fractal parameters at the top
of the panel.  The last loaded NeXus path is remembered between sessions
and re-loaded automatically on the next launch (silently cleared if the
file has been moved or deleted).

---

## Panel layout

The window is a horizontal splitter:

- **Left** — scrollable controls panel:
  1. Title bar with red `? Help` button (this document).
  2. **NeXus file (optional)** — file picker plus a coloured target summary.
  3. **Growth parameters** — Z, sticking probability, # test paths,
     Rg primary, allowed neighbor distance, multi-particle attraction,
     random seed.
  4. **Mode tabs** — Grow One / Grow Many / Optimizer.
  5. **Active jobs** — running and queued growth jobs with progress.
  6. **Stored aggregates** (session-only list) — every grown aggregate
     accumulates here; click to make it active.
  7. **Active aggregate parameters** — Z, dmin, c, df, R, p, s, true
     sticking %, Rg primary, Rg aggregate, primary diameter, # endpoints.
     Shown alongside their target values when a target is loaded.
  8. **I(Q) panel** — Q range and Monte-Carlo button.
  9. **Save** — write selected aggregate to NeXus, or load aggregates
     from an existing NeXus.

- **Right** — vertical splitter:
  - **Top** (~40 %): log-log I(Q) plot showing optional loaded data,
    optional loaded Unified-fit model, the active aggregate's analytical
    Unified curve (green), and (if computed) the active aggregate's
    Monte-Carlo curve (blue dashed).  Legend shows the four series.
  - **Bottom** (~60 %): 2D slice viewer (left) and 3D PyVista voxelgram
    viewer (right), each with a `Pop out ⤢` button.

---

## Workflow

### Without reference data

1. Set growth parameters (defaults are good starting values).
2. Click **Grow One**.  Within ~10 s a new aggregate appears in the
   "Stored aggregates" list.
3. Click the new entry; the I(Q) plot, 2D slice, and 3D viewer populate.
4. Note the computed fractal parameters in the "Active aggregate
   parameters" panel.
5. Adjust parameters and grow more aggregates to explore the morphology
   space.  Use **Grow Many** to queue a batch.

### With a reference Unified-fit result

1. Click **Open NeXus…** and pick an HDF5 file with a Unified-fit result
   that contains a fractal level pair (consecutive levels where the high
   level's `RgCutoff` ≈ the low level's `Rg`).
2. The target fractal parameters appear in the top NeXus box AND just
   above the "Active aggregate parameters" widget.  See [Target
   parameters from a Unified-fit pair](#target-parameters-from-a-unified-fit-pair)
   for the formulas used.
3. Set Z in the Growth parameters to the target Z, and use the
   **Optimizer** tab with `target_dmin` and `target_c` set from the
   targets shown.  Click **Find Best Growth**.
4. The optimizer grows several aggregates and selects the one whose
   computed `dmin` and `c` are closest to the targets.  All trials
   accumulate in the Stored aggregates list; the best is bolded.
5. Compare the active aggregate's analytical Unified curve (green) to
   the loaded data and Unified-fit model (blue scatter / red line) in
   the I(Q) plot.  Run **Compute Monte Carlo I(Q)** to overlay the
   slow but more direct PDF-based intensity (blue dashed).

---

## Loading a NeXus reference file

The Fractals tool can read both the experimental I(Q) and the stored
Unified-fit result from the same NeXus file.  When you click **Open
NeXus…**, the tool:

1. Loads the I(Q) data (`Q`, `I`, `Idev`) via the standard NXcanSAS
   reader.
2. Loads the Unified-fit results (`G_i`, `Rg_i`, `B_i`, `P_i`,
   `RgCutoff_i` for every level).
3. Scans levels for the first **consecutive pair** `(low, high)` where
   `levels[high].RgCutoff` is within ±25 % of `levels[low].Rg`.  This
   pair is the fractal representation: low = primary particles,
   high = aggregate.
4. Computes target fractal parameters (Igor formulas — see next section).
5. Auto-fills the I(Q) plot's Q range with the data's Q range.

If no qualifying level pair is found, the tool falls back to "no targets"
and growth still works without comparison.

The last loaded NeXus path is saved to the application state so the next
session starts with the file pre-loaded.

---

## Target parameters from a Unified-fit pair

Given a fractal level pair (`level_low` = primary, `level_high` = aggregate):

```
df      = P2
gamma   = exp( lgamma(P2 / 2) )         # Γ(P2 / 2)
dmin    = B2 · Rg2^P2 / (gamma · G2)
c       = P2 / dmin
Z       = G2 / G1 + 1
```

These are exactly the formulas used by Irena (`IR3A_*` in
`IR3_3DModels.ipf`).  Notes:

- `df` (fractal dimension) is just the Porod exponent of the aggregate
  level.
- `dmin` (minimum dimension) describes the topology of paths through
  the aggregate.  `dmin = 1` is a linear chain; `dmin = df` is a
  compact sphere.  Typical: 1.1–1.4.
- `c = df / dmin` is the connectivity dimension.  Typical: 1.1–1.4.
- `Z` (degree of aggregation) = number of primary particles.  The `+1`
  accounts for the central particle in Beaucage's level convention.

When no Unified-fit pair is loaded, you can still set `target_dmin` and
`target_c` manually in the Optimizer tab.

---

## Growth parameter reference

| Parameter | Default | Description |
|---|---|---|
| **Degree of aggregation Z** | 250 | Number of primary particles. Larger Z → more reliable fractal parameters but slower growth. |
| **Sticking probability [%]** | 75 | Probability per contact event that the random-walking particle sticks. Low SP → compact (high df). High SP → loose (low df). |
| **Number of test paths** | 2500 | Maximum unique paths per endpoint enumerated when computing fractal parameters. Affects c (and therefore dmin). |
| **Rg primary [Å]** | 10 | Radius of gyration of the primary sphere. Sets the physical scale (`primary_diameter = 2·√(5/3)·Rg ≈ 2.58 Rg`). |
| **Allowed neighbor distance** | Body diagonal | Which lattice neighbors count as "in contact": Edge (6 nbrs), Face diagonal (18 nbrs), Body diagonal (26 nbrs). |
| **Multi-particle attraction** | Neutral | How sticking probability changes when ≥2 existing particles are within reach: Neutral / Attractive / Repulsive / Not allowed. |
| **Random seed** | 0 (random) | Set non-zero for reproducible aggregates. |

---

## Mode tabs

### Grow One

Click **Grow** to grow a single aggregate with the current parameters.
The job runs in a background QThread; the panel stays responsive and
you can inspect previously-completed aggregates while it runs.

### Grow Many

Set **N aggregates** (default 5) and click **Grow N**.  N growths are
queued sequentially.  When `seed > 0`, each gets `seed + i` so the
batch is reproducible; otherwise each is fully random.  Useful for
assessing variance in fractal parameters with the same input.

### Optimizer ("Find Best Growth")

Set `Target dmin`, `Target c`, `Tolerance`, and `Max iterations`.
Click **Find Best Growth**.  The optimizer bisects the sticking
probability over `[10, 90]` %, growing 3 trial aggregates per iteration
(low / mid / high SP), keeping the best by objective
`(dmin − target_dmin)² + (c − target_c)²`.  Stops at `max_iter` or
when the objective drops below `tolerance²`.

Every trial appears in the Stored aggregates list with a label like
`opt-iter-3-sp42` so you can compare them.  The best is bolded after
the run completes.

---

## Active aggregate parameters

When you click an aggregate in the Stored aggregates list, this widget
populates with:

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

When a target is loaded, each comparable field shows
`actual   (target: X)` and is colour-coded:

- Green when within 10 % of the target.
- Orange when within 25 %.
- Red otherwise.

A bold compact target summary line is shown above the table so the
relationship between actual and target stays visible without scrolling
back to the top of the panel.

---

## I(Q) back-calculation

Two paths produce I(Q) for the active aggregate:

### Analytical Unified (always)

A two-level Beaucage Unified-fit closed-form intensity:

- **Level 1** (primary sphere): `G=1, Rg=Rg_primary, P=4, B=4π/Rg_primary⁴`
- **Level 2** (aggregate): `G=Z, Rg=Rg_aggregate, P=df,
  B = (G·P / Rg^P) · Γ(P/2), RgCutoff = Rg_primary`

`Rg_aggregate` here is the **direct measurement from the grown particle
positions** via the parallel-axis (Steiner) theorem:

```
Rg_centers (Å) = sqrt(mean(||p_i − centroid||²)) × primary_diameter
Rg_aggregate   = sqrt(Rg_centers² + Rg_primary²)
```

(The McGlasson approximation `Rg_primary · Z^((1/c−1)/(dmin−df))` is
not used because it underestimates the true Rg by 30–100 % depending
on the morphology, which makes the analytical curve disagree with the
MC curve.)

The analytical curve is shown in **green**.

### Monte-Carlo PDF (on demand)

Click **Compute Monte Carlo I(Q)** to run the slow but more direct
calculation:

1. Voxelize the aggregate: place a sphere of physical radius
   `R = primary_diameter / 2` around every lattice position on a
   20×-oversampled cubic grid (sphere kernel = 10 voxels, voxel pitch
   = R / 10 in physical units).
2. Sample random pairs of solid voxels uniformly with replacement and
   build a histogram of their Euclidean distances.  This histogram IS
   the Glatter pair-distance distribution `p(r) = 4π·r²·γ(r)`.
3. Apply the Debye / Glatter sine transform:
   ```
   I(Q) ∝ ∫₀^∞ p(r) · sinc(Qr) dr
   ```
4. Truncate above `Q_voxel = π / pitch_A` (the voxel-grid Nyquist Q,
   above which the discrete voxel approximation is unreliable).

The MC curve is shown in **blue dashed**.  Costs: ~3–5 s per call,
~40 MB voxelgram for typical Z=80 aggregates.

### Invariant rescaling

When a Unified-fit reference is loaded, both model curves (green and
blue) are scaled to the data using the integral
`∫ I_data dQ / ∫ I_model dQ` over the **fractal regime**
`Q ∈ [0.5π/Rg_aggregate, 1.5π/Rg_primary]`.  This window deliberately
EXCLUDES very low Q (sample-level power-law contamination) and very
high Q (instrumental flat background) that the single-aggregate model
never tries to reproduce.

---

## 3D and 2D viewers

### 3D viewer (right, PyVista isosurface)

- Renders the voxelgram as a flying-edges isosurface at value 0.5.
- Mesh colour: dark grey by default (right-click → "Pick isosurface
  colour…" to change).
- White background with black axes and tick labels.  Axis units are
  in Ångström (`X [A]` / `Y [A]` / `Z [A]`).
- Right-click context menu: reset view, change colour, toggle bounding
  box, save screenshot.
- Click `Pop out ⤢` to open the viewer in a separate dialog you can
  resize / maximize.

The display uses `oversample=10, sphere_voxel_radius=10` (Irena-style
"fat spheres" of radius D = primary_diameter, twice the physically
correct R).  This is intentional for **visualization only** — it makes
all neighbor types (edge, face diagonal, body diagonal) overlap so the
aggregate looks like a connected blob rather than scattered particles.
The MC scattering calculation uses the physically correct
`oversample=20, sphere_voxel_radius=10` independently.

### 2D slice viewer (left, pyqtgraph)

- Shows one slice through the voxelgram.
- White background, dark grey solid (majority phase), white voids
  (minority phase).  Axes black; units `[A]`.
- Combobox at top selects the slice plane (XY / XZ / YZ).
- Slider scrubs the slice index; a label shows the current physical
  coordinate.
- Pop out via the button below.

---

## Saving and loading aggregates

### Save

Select an aggregate in the list, click **Save selected aggregate to
NeXus…**, and pick a file path.  The aggregate is appended as a new
`entry/fractals_results/aggregate_{N}` group inside an HDF5/NeXus
file.  N auto-increments; multiple aggregates can coexist in the same
file.  Other groups (Unified fit, raw data, other aggregates) are
preserved.

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

Click **Load aggregate(s) from NeXus…** to pull every aggregate from
an existing file into the session list.  Useful for resuming work on
saved aggregates or for cross-comparing aggregates from multiple files.

---

## Algorithm details and Igor compatibility

The following details are reproduced verbatim from Irena's `IR3A_*`
functions to guarantee that a pyirena aggregate matches an Irena
aggregate grown with identical inputs:

- **Sticking probability table** (per contact count `chcnt`):
  - `chcnt == 1`: base SP
  - `chcnt ≥ 2` (Attractive): `(SP+100)/2` for chcnt=2, `(SP+300)/4` for ≥3
  - `chcnt ≥ 2` (Repulsive): `(SP+10)/2` for chcnt=2, `(SP+30)/4` for ≥3
  - `chcnt ≥ 2` (Not allowed): 0
  - `chcnt ≥ 2` (Neutral): base SP
- **Neighbor distance thresholds**: `1.1`, `1.05·√2`, `1.05·√3`
  (intentionally lenient).
- **Path-walk junction rule** (computing dmin, c): at junctions with
  more than 3 outgoing neighbors, only the first 3 are explored
  (matches `IR3A_MT_NextPathStep`).  Removing this rule changes c and
  dmin and breaks comparison against historical Irena results.
- **Wall launch**: 6 wall faces normally; 2 corner-launch fallbacks
  for compact aggregates with Z < 6.
- **Voxel oversample 10×, sphere kernel 10 voxels** (Irena-style
  display) — see [3D and 2D viewers](#3d-and-2d-viewers).

---

## Slit smearing (USAXS)

When the loaded I(Q) data are slit smeared (an NXcanSAS `dQl`), the calculated
intensity is slit smeared for the overlay/comparison so it can be judged against
the data on equal footing. If smearing the comparison curve fails it is logged
and flagged rather than silently shown unsmeared. See
**[Slit smearing](slit_smearing.md)** for the full reference.

## Performance and memory notes

- Growth: ~0.3 s for Z=50, ~3 s for Z=250, ~30 s for Z=1000 (single
  thread, default parameters).  Body-diagonal sticking is fastest;
  edge-only is slowest because there are fewer contact opportunities.
- Path enumeration: parallelized across endpoints via
  `ThreadPoolExecutor` (Python threads, no GIL release for pure-Python
  walks but I/O-bound enough to help).
- Voxelization (display): ~8 MB at Z=80, lighter — uses
  `oversample=10`.
- Voxelization (MC scattering): ~80 MB at Z=80, ~340 MB at Z=500 —
  uses `oversample=20` for finer grid, doubling the voxel-Nyquist Q.
- Monte-Carlo intensity: 3–5 s for typical aggregates with the default
  budget of 1e7 pair samples / 20 s wall-time.

---

## Notes and limitations

- This is a visualization tool — there is no fitting of measured data
  to growth parameters (the optimizer just searches for matching
  `dmin` and `c`, both computed from the grown aggregate, not
  measured from data).
- Aggregates accumulate in the session list in memory only.  Use
  **Save selected aggregate to NeXus…** to persist any you want to
  keep across sessions.
- The 3D viewer requires PyVista; on systems without it, only the 2D
  slice and the I(Q) plot work.  Install with `pip install pyirena[gui3d]`.
- No batch / scripting / CLI entry point — intentionally GUI-only
  given the qualitative-visualization purpose.
