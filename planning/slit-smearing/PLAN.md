# Slit Smearing Support — Implementation Plan

Branch: `feature/slit-smearing`
Test data: `testData/BoehNaNO2_10m_34_84min_0923.h5` (Matilda USAXS file with both
desmeared and slit-smeared entries; slit length dQl = 0.017700 1/A, Q range
1.04e-4 – 0.30 1/A for both entries)

Status: PLANNING — written 2026-07-18, execution over multiple future sessions.

---

## 1. Goal and scope

Add the major Irena capability that pyIrena skipped for speed: fitting
**slit-smeared USAXS data** directly. The model (never the data) is smeared
before comparison with data and before display. Slit smearing is a
deterministic, single-pass linear operation, so the core is simple; the work
is in careful wiring and edge-case handling.

**Tools in scope** (per Jan, 2026-07-18):

| Tool | What it gets |
|---|---|
| Unified Fit | Full support: load SMR data, smear model in fit loop, smear local Guinier / power-law fit models |
| Size Distribution | Full support: smear the G matrix columns (linear operator ⇒ smear once, fit fast) |
| Modeling | Full support: smear total model intensity per dataset in fit loop |
| Simple Fits | Full support: smear each analytic model; local-fit edge cases apply |
| Fractals | Option: load SMR data if present, smear calculated intensity for comparison/fitting |
| Data Merge | Option: merge slit-smeared USAXS with pinhole SAXS (valid when slit length ≈ qmin of SAXS) |
| Data Manipulation | Support SMR data as input; warn when mixing smeared + unsmeared curves |

**Out of scope:** WAXS (obviously), saxsMorph (revisit later — smearing a
morphology-derived I(q) is the same one call, but the tool's fitting story is
unclear; leave a hook, no UI), desmearing tool (pyIrena philosophy: fit
smeared data directly, don't desmear).

---

## 2. Physics and algorithm (reference: Irena Igor code)

### 2.1 The smearing integral

For an "infinitely long" slit of half-length `SL` (Lake / APS USAXS
convention, `SL` = dQl):

    I_smeared(q) = (1/SL) * ∫_0^SL  I_ideal( sqrt(q² + l²) )  dl

Igor reference: `IR1B_SmearData` in
`User Procedures/Irena/IR1_Desmearing.ipf` (SAXS_IgorCode repo, line ~604).
Key details of the Igor implementation we will mirror or improve:

1. **l-grid**: points distributed 0 → 2·SL mimicking the spacing of the
   data's own Q grid (`Smear_Q = 2·SL·(Q[p]−Q[0])/(Q[N−1]−Q[0])`),
   integration by trapezoid (`areaXY`) over 0..SL, normalized by 1/SL.
2. **Internal data extension**: evaluation needs I at sqrt(q²+l²) up to
   sqrt(qmax²+SL²) > qmax. Igor extends the tabulated curve by up to 300
   points with a linear ramp to zero at 20·qmax — a crude extrapolation that
   is the source of the "artifacts" concern.
3. **Trapezoidal slit** variant (`IR1B_SmearDataTrapeziod`) exists in Igor
   but is effectively unused at APS; we implement rectangular profile only
   (v1), keeping the function signature open for a slit-width parameter later.

### 2.2 The pyIrena improvement: smear the *model*, not a tabulated curve

Igor smears a tabulated curve because of its wave-based architecture. In
Python we can do better: everywhere we control the model we evaluate it
directly on an **extended q grid**, so no extrapolation of tabulated values
is ever needed and the artifact problem largely disappears:

- Build `q_ext` = data q grid extended to `max(qmax, sqrt(qmax² + SL²))`
  with margin (Irena uses 3·SL as the target; see `IR1A_UnifiedCalculateIntensity`,
  `IR1_UnifiedFitFncts.ipf` line ~196: if `qmax < 3·SL`, append points with
  spacing `2·Δq_last` until 3·SL is reached).
- Evaluate the model on `q_ext` (models are analytic — exact, cheap).
- Smear, then truncate back to the original grid.

Tabulated-curve smearing (`smear_curve`) is still needed for two cases:
validation tests (smear the desmeared data and compare to SMR) and Fractals
Monte-Carlo intensity (tabulated). For those, extrapolation is explicit and
selectable: `power_law` (fit last decade, default), `flat`, `linear_zero`
(Igor-compatible), `none` (raise if extension needed).

### 2.3 Smearing as a matrix (performance)

Smearing is linear in I. For a fixed `(q, q_ext, SL)` build once a sparse
**smearing matrix W** (n_q × n_q_ext), rows = trapezoid weights of the
l-integration after interpolation onto q_ext:

    I_smeared = W @ I_ideal(q_ext)

Benefits: fit loops cost one matvec per iteration; Size Distribution smears
the whole G matrix once (`G_smeared = W @ G(q_ext, r)`) and the inversion
machinery (MaxEnt / regularization / TNNLS / MC) is untouched.

### 2.4 Edge cases to handle explicitly

These are the cases that "bomb and confuse users" in Igor:

- **E1 — qmax < 3·SL** (data don't span the slit): extend model grid as in
  2.2. Warn (log + GUI status) that model beyond qmax is an extrapolation of
  the *model*, which is safe as long as the model is sensible there.
- **E2 — local fit with qmax_fit ≪ SL** (user fits only very low q):
  NEVER evaluate + smear the model only on the fit sub-range. Igor's
  pattern (`IR1_FitLocalGuinier` all-at-once functions,
  `IR1_UnifiedFitFncts.ipf` ~895): evaluate model on the FULL data q grid
  (extended per E1), smear, then extract the fit sub-range by index. We do
  the same. Additionally: if `q_max_fit < SL`, emit a clear warning
  ("fit range is below the slit length; smeared local fits are weakly
  determined") rather than failing cryptically.
- **E3 — Rg vs slit length**: smeared Guinier knee with Rg ≫ π/SL is
  strongly distorted; add a feasibility note in `check_level_feasibility`
  and `detect_features` (they currently assume pinhole shapes).
- **E4 — constant background**: smears to itself (∫bkg/SL = bkg) — no
  special handling needed, but unit-test it.
- **E5 — non-monotonic / unsorted q, NaNs**: sanitize before building W
  (reuse existing data-cleaning path).
- **E6 — SL = 0 or missing**: treat as pinhole; the smeared path must be a
  strict no-op so tools can call it unconditionally.
- **E7 — log-spaced vs linear q**: l-grid construction must not assume
  either; use the Igor "mimic the data spacing" trick which handles both.

---

## 3. Data model: where slit-smeared data live

### 3.1 In files (NXcanSAS — already standard)

Matilda writes both entries; verified in the test file:

```
entry/<sample>/sasdata          Q, I, Idev, Qdev            (desmeared; Q@resolutions="Qdev")
entry/<sample>_SMR/sasdata      Q, I, Idev, dQw, dQl        (slit-smeared; Q@resolutions="dQw,dQl")
```

- `dQl` — scalar dataset = slit (half-)length in 1/A. THE source of SL.
- `dQw` — per-point slit-width resolution; becomes the loaded `dQ`.
- Detection rule already in codebase: group path contains `_SMR`
  (`pyirena/io/hdf5.py::_is_smr`). Keep it, but ALSO detect the NXcanSAS-
  standard way: `Q@resolutions == "dQw,dQl"` — files from other pipelines
  won't use the `_SMR` naming.

### 3.2 In memory — the `slit_length` contract

Every tool currently receives `(q, intensity, error, dq)` in some form.
Add two fields carried alongside, end to end:

- `data_is_slit_smeared: bool`
- `slit_length: float`  (0.0 ⇒ pinhole)

Concretely:

- `io/hdf5.py::_read_one_sasdata` returns dict — add keys
  `'slit_length'` (from `dQl`, else 0.0) and `'is_slit_smeared'`
  (from resolutions attr / `_SMR` path). dQ for SMR = `dQw`.
- `readGenericNXcanSAS` gains `prefer_slit_smeared=False` kwarg:
  when True and the file has an SMR sibling of the default entry, load it
  instead. Implement via `_ordered_sasdata(f, include_smr=...)` — today
  `_filter_smr` unconditionally strips SMR entries; make that behavior
  switchable rather than removing it (the picker default stays desmeared).
- GUI panels: `set_data(...)` in each panel gains
  `slit_length=0.0` keyword (all have compatible signatures:
  `gui/unified_fit.py::set_data` ~2544, `sizes_panel` ~1850,
  `modeling_panel` ~2889, `simple_fits_panel` ~1271).

---

## 4. New core module: `pyirena/core/smearing.py`

Single home for all smearing math. No tool implements its own.

```python
def build_extended_q(q, slit_length, min_span_factor=3.0) -> np.ndarray
    # q extended to >= min_span_factor*SL past nothing-below-qmax,
    # spacing = 2*Δq_last (Igor-compatible); returns q unchanged if long enough

def build_smearing_matrix(q, slit_length, q_ext=None) -> scipy.sparse.csr_matrix
    # W such that I_sm(q) = W @ I(q_ext). Rows: trapezoid weights of
    # (1/SL)∫_0^SL I(sqrt(q²+l²)) dl with linear interpolation onto q_ext.

def smear_model(model_fn, q, slit_length) -> np.ndarray
    # convenience: q_ext = build_extended_q; W = build_smearing_matrix;
    # return W @ model_fn(q_ext)

def smear_curve(q, I, slit_length, extrapolation="power_law") -> np.ndarray
    # for tabulated curves; extrapolation in {"power_law","flat","linear_zero","none"}

class SlitSmearer:
    # caches (q, SL) -> (q_ext, W); the object every fitting loop holds.
    def __init__(self, q, slit_length, min_span_factor=3.0)
    def smear(self, I_on_q_ext) -> np.ndarray
    def smear_model(self, model_fn) -> np.ndarray
    @property def q_ext
    def is_noop  # SL <= 0
```

Design rules:

- SL ≤ 0 ⇒ every call is an exact no-op (identity W); tools call
  unconditionally, no `if use_smearing:` forests.
- Pure NumPy/SciPy, no h5py, no Qt — unit-testable in isolation.
- Validate against Igor: smear a power law q^-4 (analytic slit-smeared
  slope −3), smear a flat curve (identity), smear desmeared test data and
  compare to the file's SMR entry (tolerance a few % in the mid-range,
  looser at extreme ends).

---

## 5. Per-tool wiring

### 5.1 Unified Fit (first tool — proves the whole chain)

Core (`core/unified.py`):
- Replace the placeholder `UnifiedFitModel.apply_slit_smearing` (~line 269,
  currently returns intensity unchanged with a warning) with a `SlitSmearer`
  held on the model; `use_slit_smearing`/`slit_length` attrs already exist.
- `_residuals` / `fit` / `calculate_intensity`: evaluate levels+background on
  `smearer.q_ext`, then `smearer.smear(...)`, truncate. (Mirrors Igor
  `IR1A_UnifiedCalculateIntensity`.)
- `fit_local_guinier` / `fit_local_power_law` (~659/~758): add
  `slit_length=0.0` and `q_full=None` kwargs. When smeared: build model on
  full extended grid, smear, index-select the sub-range (Igor all-at-once
  pattern, E2). Both callers (GUI cursors + control API) updated.
- Invariant (`calculate_invariant`): computed from the IDEAL model —
  document this; no change needed beyond making sure it uses unsmeared
  parameters (it does — it's parameter-based).

GUI (`gui/unified_fit.py`):
- Accept `slit_length` in `set_data`; show read-only indicator
  "Slit-smeared data, SL=0.0177 1/A" + checkbox "Fit slit-smeared data"
  (auto-checked when SMR data loaded; unchecking reverts to treating data
  as pinhole — expert escape hatch, same as Irena).
- Plot: display smeared model over data; add optional second curve
  "ideal (pinhole) model" toggle.
- Local fit buttons: route slit_length through; surface E2 warning in
  status bar instead of exception.

IO (`io/nxcansas_unified.py::save_unified_fit_results`):
- New optional args: `slit_length`, `intensity_model_ideal`.
- Store in `unified_fit_results`: existing `I_model` = **smeared** model
  (matches the data it was fitted to), new dataset `I_model_ideal` =
  pinhole model on same q, attrs `slit_length`, `data_is_slit_smeared`.
  Saving BOTH means downstream consumers (plots, Igor export, HDF5 viewer,
  ASCII export) keep working and the user picks which curve to look at.
- `load_unified_fit_results` reads both back; absent ⇒ pinhole legacy file.

Control API (`api/control/unified_fit.py` + `session.py` + `schemas.py`):
- `open_dataset(file_path, use_slit_smeared=False)`; response reports
  `is_slit_smeared`, `slit_length`, and whether an SMR entry exists.
- `run_fit`, `fit_local_guinier`, `fit_local_power_law`,
  `check_level_feasibility`, `detect_features` become smearing-aware (E2/E3
  warnings in their structured responses so the AI advisor can explain).
- MCP server tool docs updated (`mcp/server.py`).

Batch (`batch/unified.py`, `batch/_common.py`): JSON key
`"use_slit_smeared": true` in data loading block; `_common` loader passes
`prefer_slit_smeared` and returns slit_length.

### 5.2 Size Distribution (`core/sizes.py`, `gui/sizes_panel.py`, `io/nxcansas_sizes.py`)

- `_build_g_matrix` (~402): build on `q_ext`, then `G ← W @ G`. One line of
  wiring once `SlitSmearer` exists; all four methods (MaxEnt, regularization,
  TNNLS, MC) inherit it because they only see G. (Igor reference:
  `IR1_Sizes.ipf` ~3820–3920 does exactly this including 3·SL extension.)
- Background terms: power-law background must be smeared too (it is part of
  the model); flat background is invariant (E4). `fit_power_law` /
  `fit_background_term` on the smeared path use `smear_model`.
- Save: both smeared and ideal fitted-intensity curves + `slit_length` attr
  in `size_distribution_results`. Distribution itself is unaffected (it is
  ideal-space by construction).
- Control API: `sizes_suggest_setup` should read slit_length and set it;
  `pyirena_ctrl_open_dataset` shared change covers loading.

### 5.3 Modeling (`core/modeling.py`, `gui/modeling_panel.py`, `io/nxcansas_modeling.py`)

- Modeling fits N datasets simultaneously; each dataset gets its OWN
  `SlitSmearer` (per-dataset SL — one can be SMR USAXS, another pinhole
  SAXS; this is a headline use case).
- Wire at the point where total model intensity for dataset i is compared
  to data i (the `_DEObjective` / per-dataset residual assembly): evaluate
  populations on that dataset's `q_ext`, apply its W.
- Structure factors are part of the ideal model — apply BEFORE smearing.
- Save both curves per dataset + per-dataset `slit_length`.

### 5.4 Simple Fits (`core/simple_fits.py`, `gui/simple_fits_panel.py`, `io/nxcansas_simple_fits.py`)

- `SimpleFitModel`: hold a `SlitSmearer`; wrap the selected analytic model
  fn with `smear_model` in the residual. All models (Guinier family, Porod,
  power law, Debye, sphere, …) work unchanged underneath.
- E2 applies with force here (users fit narrow ranges by design): model
  evaluation on full extended grid, subselect, plus explicit warning when
  `q_max_fit < SL`.
- **Semantics note for results table**: fitted Rg/Kp/etc. are IDEAL-space
  parameters (because the model is smeared to match data) — this is the
  whole point and must be stated in the tool docs.
- Invariant calculation from SMR data: compute from the smeared *model*'s
  ideal twin, not by integrating smeared data — flag in GUI when data is
  SMR that "extrapolation-based invariant uses the ideal model".
- Save both curves + `slit_length`.

### 5.5 Fractals (`core/fractals.py`, `gui/fractals_panel.py`, `io/nxcansas_fractals.py`)

- Loading measured data for comparison: allow SMR selection (shared loader
  change).
- `intensity_unified` (analytic branch): smear via `smear_model` — exact.
- `intensity_montecarlo` (Debye sum, tabulated on its own q grid): smear via
  `smear_curve` with power-law extrapolation; MC q grid should be built out
  to `sqrt(qmax²+SL²)` when SL known, avoiding extrapolation entirely
  (`mc_q_max` already computes a natural limit — extend request there).
- Optimizer (`optimize_growth`) scores against smeared calculated intensity
  when data is SMR.

### 5.6 Data Merge (`core/data_merge.py`, `gui/data_merge_panel.py`, `io/nxcansas_data_merge.py`, `batch/merge.py`)

Use case: merge slit-smeared USAXS with pinhole SAXS when SL ≈ qmin(SAXS)
(Jan: works in practice on the APS instrument).
- Add option "use slit-smeared USAXS" in dataset-1 selection (loader change
  covers it). NO desmearing, no smearing of the SAXS side — this is a
  deliberate approximation the user opts into.
- Overlap-region scaling math (`_overlap_data_ds1_grid`, `_wls_bg_scale2`)
  works on the curves as given — unchanged.
- **Provenance is the critical part**: output NXcanSAS must declare the
  merged curve honestly. Write `Q@resolutions="dQw,dQl"` with `dQl` = SL
  over the USAXS portion... but NXcanSAS has one resolution declaration per
  curve. DECISION (see §8/Q3): v1 writes the merged curve as slit-smeared
  (`dQl` scalar kept, `dQw` per-point, entry name suffix `_SMR`), with a
  `smeared_portion_qmax` attribute recording where USAXS ends; merge
  provenance note records both inputs. Rationale: for q > qmin(SAXS),
  SL ≲ q so smearing is a small perturbation — labeling the whole curve
  smeared is the conservative choice and downstream pyIrena tools can then
  fit it with smearing enabled (correct at low q, negligible at high q).
- Merge provenance reader (`pyirena_read_merge_provenance`) reports the flag.

### 5.7 Data Manipulation (`core/data_manipulation.py`, `gui/data_manipulation_panel.py`, `io/nxcansas_data_manipulation.py`)

- Allow SMR input (loader change). Scale/trim/rebin/average preserve
  smearing status: propagate `slit_length` + `is_slit_smeared` to output,
  keep `dQw`/`dQl` writing.
- Subtract/divide of two curves: legal only when both operands have the SAME
  smearing status and SL (within tolerance). Mixed operands ⇒ hard error
  with clear message (not a warning — the result would be silently wrong).
- Averaging N SMR curves: require matching SL; output keeps it.

### 5.8 Shared surfaces

- `gui/data_selector/panel.py`: per-file "load slit-smeared" option
  (context-menu or checkbox column shown only when file has an SMR entry —
  `list_nxcansas_datasets` with `include_smr=True` tells us); passes
  `slit_length` into every `tool.set_data(...)` call.
- `io/ascii_export.py`: SMR-aware export exists half-way (reads dQw+dQl);
  make export of SMR data explicit (columns Q, I, dI, dQw + header line
  `# slit_length = ...`).
- HDF5 viewer / Igor export (`gui/hdf5viewer/*`): show both model curves,
  label them; low priority, mostly free because both are plain datasets.
- `api/results.py` / readers (`pyirena_read_unified_fit` etc.): return both
  curves + slit_length when present.
- `io/schema.py` + `docs/HDF5_NxcanSAS_structure.md`: document the new
  datasets/attrs.

---

## 6. Output convention (decided)

Per Jan: save BOTH slit-smeared and pinhole (ideal) model curves so
downstream needs no modification and the user picks.

Standard for every tool's results group:

```
<tool>_results/
    Q                     # data q grid
    I_model               # SMEARED model when fit was smeared (matches data) — name unchanged for back-compat
    I_model_ideal         # pinhole model on same Q (only written when smearing was used)
    @slit_length          # scalar, 1/A; absent or 0 ⇒ pinhole fit
    @data_is_slit_smeared # 1/0
```

Legacy readers see `I_model` and keep working; new readers/plots offer the
toggle. Experimental data saved alongside results keep NXcanSAS SMR
declaration (`dQw`,`dQl`) — "in NXcanSAS that is easy to declare" — so any
canSAS-aware software understands the file.

---

## 7. Testing strategy

New file `pyirena/tests/test_smearing.py` (core) + per-tool additions.

Unit (no files needed):
1. SL=0 ⇒ identity (all APIs, incl. matrix).
2. Flat curve ⇒ unchanged (E4).
3. Power law q^-P ⇒ smeared slope ≈ −(P−1) for q ≪ SL (infinite-slit
   analytic result) and ≈ −P for q ≫ SL; check both log-log regimes.
4. Guinier: smeared curve ≥ has known analytic series; check against direct
   numerical quadrature (`scipy.integrate.quad` per point) to 1e-6.
5. Matrix vs direct quadrature agreement on random smooth curves.
6. Extension: qmax < 3·SL grid building; E2 subrange indexing; E5 NaN input.

Integration (uses `testData/BoehNaNO2_10m_34_84min_0923.h5`):
7. `smear_curve(DSM data, SL=dQl)` ≈ SMR data from same file (the file IS
   the ground truth pair — this is why it was added). **Already validated
   during planning (2026-07-18)**: a 20-line prototype (trapezoid l-integral,
   200 points, power-law extrapolation from the last decade) reproduces the
   file's SMR curve with median |rel. dev.| = 0.11%, 90th pct = 0.20% over
   the mid range (5·qmin – 0.5·qmax); worst point 22% at the extreme ends
   where desmearing itself is ill-conditioned. Test tolerance: median < 1%
   mid-range.
8. Unified Fit on SMR entry vs same fit on DSM entry ⇒ consistent Rg/G/B/P
   within uncertainties (the physics must agree).
9. Round-trip: save smeared unified fit → `load_unified_fit_results` →
   both curves + attrs back (`test_nxcansas_roundtrips.py` pattern).
10. Sizes on SMR vs DSM entry ⇒ similar volume distributions.
11. Batch JSON with `use_slit_smeared` runs headless.
12. Control API: `open_dataset(use_slit_smeared=True)` reports SL; run_fit
    works; E2 warning surfaces in `fit_local_guinier` response.

GUI smoke (manual checklist, `docs/` note): load test file, verify
indicator, fit, toggle ideal/smeared display curve, save, reopen.

---

## 8. Open questions to confirm with Jan before/while coding

- **Q1 — l-grid resolution**: Igor mimics data spacing (≈N points over the
  slit). Proposal: fixed 100–200 log-ish points over [0, SL] for W
  construction — accuracy set by test #4, independent of data gridding. OK?
- **Q2 — GUI default**: when a file has an SMR entry, default remains
  loading the DESMEARED entry (current behavior), and slit smearing
  activates only when the user explicitly picks the SMR dataset. Matches
  "don't surprise existing users". OK?
- **Q3 — merged-curve declaration** (see 5.6): whole merged curve declared
  slit-smeared with `smeared_portion_qmax` attr. Alternative: declare
  desmeared and record "low-q portion smeared" in provenance only.
  Current plan chooses the conservative option; confirm.
- **Q4 — Modeling mixed datasets**: when fitting SMR USAXS + pinhole SAXS
  simultaneously, any need for per-dataset override of SL in the GUI, or is
  file-derived SL always trusted? Plan assumes file-derived + editable field.
- **Q5 — Simple Fits invariant on SMR data**: warn-and-compute-from-ideal-model
  (planned) vs. disable? (Irena precedent unclear.)

---

## 9. Execution phases (multi-session)

Each phase ends green (`pytest pyirena/tests`) and committable.

**Phase 1 — Core smearing engine** (no UI):
`core/smearing.py` + `tests/test_smearing.py` units 1–6 + integration test 7
(the DSM→SMR file check). This validates the math before any wiring.
~1 session.

**Phase 2 — Data plumbing**:
`io/hdf5.py` (SMR-aware `_read_one_sasdata`, `prefer_slit_smeared`,
`include_smr`), `list_nxcansas_datasets`, data selector option, `set_data`
signatures (all panels accept + store slit_length; no behavior yet),
`ascii_export` header. Tests: loader returns SL for test file; picker lists
both entries. ~1 session.

**Phase 3 — Unified Fit end-to-end** (reference implementation):
core model + local fits + GUI + save/load both curves + control API + batch.
Integration tests 8, 9, 12. This is the template PR for every other tool.
~1–2 sessions incl. debugging.

**Phase 4 — Size Distribution**: G-matrix smearing + background terms + GUI
+ save. Test 10. ~1 session.

**Phase 5 — Simple Fits + Fractals**: analytic-model wrapping, E2 warnings,
fractals MC q-grid extension. ~1 session.

**Phase 6 — Modeling**: per-dataset smearers, mixed-dataset test
(SMR USAXS + pinhole SAXS simultaneous fit). ~1 session.

**Phase 7 — Data Merge + Data Manipulation**: options, provenance/labels,
mixed-operand guards. Batch merge flag. ~1 session.

**Phase 8 — Polish**: docs (`docs/slit_smearing.md` user guide +
`developer_adding_features.md` cross-ref + `HDF5_NxcanSAS_structure.md`),
HDF5 viewer labels, MCP tool doc strings, CHANGELOG, GUI smoke checklist.
~1 session.

Dependency order: 1 → 2 → 3 → {4, 5, 6, 7 in any order} → 8.

---

## 10. Risk notes

- **Performance**: W is built once per (q, SL) and cached; Unified fit adds
  one sparse matvec per residual call — negligible. Sizes adds one dense
  multiply at setup — negligible. MC-uncertainty loops in Unified reuse the
  smearer.
- **Back-compat**: `I_model` name reused for the fitted (smeared) curve —
  old readers keep working but will silently plot a smeared curve; the
  `@slit_length` attr is the discriminator. Readers we ship are updated in
  the same phase as each writer.
- **Numerical**: interpolation inside W is linear in log-log space
  (intensities span decades); test 5 guards accuracy.
- **User confusion (the Irena lesson)**: every failure mode surfaces as a
  worded warning at the point of action (GUI status bar / structured API
  warning), never a stack trace. E2 is the one that historically "bombs";
  it gets a dedicated message.
