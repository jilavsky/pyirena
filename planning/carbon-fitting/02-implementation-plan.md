# Implementation plan — wiring the new tool

Follows `docs/developer_adding_features.md` § "Adding a whole new tool"
throughout. That guide is the authority on *how* to wire a tool into
pyIrena's 7-layer stack and ~15 hand-maintained registries; this document is
the carbon-tool-specific application of it, plus a reuse map so we don't
reimplement physics pyIrena already has.

## 0. Proposed name and keys

Working name: **Carbon Fitting** (module/schema key: `carbon_fit`). Follows
the existing naming style (`waxs_peakfit`, `simple_fits`). Per
`test_tool_registration.py`'s own warning, the key differs by surface in
existing tools (e.g. `simple_fit` in Data Explorer vs `simple` in MCP
dispatch) — for a new tool we should **pick one spelling and use it
everywhere it isn't forced otherwise**, to avoid adding to that inconsistency
debt. Proposed: `carbon_fit` everywhere except the MCP dispatch category,
which by existing convention is a short noun (`unified`, `sizes`, `simple`,
`modeling`, `waxs`) — propose `carbon`. Confirm with Jan before coding (see
[03-open-questions.md](03-open-questions.md)).

## 1. What already exists — reuse map (read this before writing new math)

This is the single most important finding of this planning pass: **a large
fraction of the required physics already exists in pyIrena**, mostly inside
the **Modeling** tool's population system and **Simple Fits**' model
registry, because both tools were built generically enough to already cover
several of Saurel et al.'s models under different names. Reimplementing them
inside a new `core/carbon_fit.py` would violate AGENTS.md §6 ("prefer editing
the existing tool that already does 80% of the job... duplication across the
seven layers is expensive") — instead, the new tool's core module should
**import and call** this existing math, not duplicate it.

| Needed physics | Paper eq. | Already in pyIrena as | Notes |
|---|---|---|---|
| Teubner-Strey two-phase model | eq. 8 | `simple_fits.py::_teubner_strey` + `MODEL_REGISTRY['Teubner-Strey']` | Exact formula match (`Prefactor/(A+C1Q²+C2Q⁴)`, A≡1). Already derives `CorrLength`(ξ), `RepeatDist`(d). Missing: `f_a`, and all the ρ/SLD-dependent derived quantities (φ, r, S_mp, w_P, w_C) — those are new, carbon-tool-specific, because Simple Fits has no concept of density/contrast. |
| Fractal aggregate structure factor (Teixeira) | eq. 5, A3.17 | `modeling.py::_mass_fractal_intensity` (`MassFractalPopulation`) | Implements Teixeira 1988 SF combined with a spheroid form factor. Very close to eq. 5–6's fractal+globule combination — check whether the existing spheroid-based form factor is an acceptable substitute for the paper's specific "algebraic globule" (A3.9), or whether the globule needs adding as a new pyIrena form factor (`docs/developer_adding_form_factors.md` gives the exact recipe: `_build_g_your_ff` in `form_factors.py` + `_G_BUILDERS` entry + `FF_LABELS`/`FF_PARAM_LABELS` in `modeling_panel.py`). Likely: add `globule`/`algebraic_globule` as a first-class form factor (short, ~30 lines, per the guide's template) rather than approximating with spheroid. |
| Surface fractal / crumpled-layer structure factor | A3.17 (same function, D≈2.5) | `modeling.py::_surface_fractal_intensity` (`SurfaceFractalPopulation`) | Same Teixeira family; already handles `Ds`, `Ksi`. Compare against eq. 17's `S_3D/Q²` before deciding to reuse vs. write a carbon-specific variant with the extra `/Q²` orientation factor folded in. |
| Beaucage unified Guinier-Porod (algebraic globule building block) | A3.3–A3.9 | `pyirena/core/unified.py` (whole Unified Fit tool) + `form_factors.py` sphere/spheroid builders | The "exp(Guinier) + erf-blended Porod" pattern (A3.3, A3.14) is exactly Beaucage's unified level, already implemented as a general multi-level engine. Consider literally reusing `unified.py`'s single-level evaluator as a library call for the globule/discoid form factors, rather than re-deriving the erf blending by hand. |
| Voigt/Lorentz/Gauss peak shapes, background subtraction, peak finding | eq. 2, A3.30–A3.43 | `pyirena/core/waxs_peakfit.py` (`WAXSPeakFitModel`, `gauss_peak`/`lorentz_peak`/`pseudo_voigt_peak`, `eval_peak`, background methods) | Has Gauss/Lorentz/Pseudo-Voigt (linear mix) but **not a true Voigt (convolution)** — this is genuinely new (see §3). Everything else — the peak table UI pattern, `find_peaks_in_data`, background estimation, `peak_area` — is directly reusable infrastructure and a good model for the WAXS tab's peak list widget. |
| `use_complex_bg` (power-law + flat background, fit-between-cursors, "prefit replay" persistence) | eq. 3/4 (generalized) | `simple_fits.py` `SimpleFitModel.use_complex_bg`, `bg_prefit`, `prefit_background()` | This is **exactly** the "Complex background: low-Q power law slope + flat background" tab Jan described, already built, tested, and battle-hardened (it's the subject of the Invariant case study in `developer_adding_features.md` — read that case study, the "prefit replay" lesson applies directly here: whatever the Background tab does interactively must be replayable from a config, not just from typed values). Reuse this mechanism/pattern directly rather than inventing a new background system. It doesn't have the roughness (`S_rough`/`R_rough`) term — that's the new piece for this tool (see §3). |
| Density → SLD → contrast (formula + density → X-ray SLD, contrast) | SI eq. S1/S2 | `pyirena/core/scattering_contrast.py` (`compute_compound`, `compute_contrast`) | Already does exactly this via NIST/Chantler tables. New: a thin helper that goes `ρ_struc → SLD_C`, `ρ_sample → SLD_sample`, wraps the corrected `ρ_struc` formula (§5 of 01-models.md) which is carbon-crystallography-specific and not in `scattering_contrast.py` today. |
| Diffraction peak population (single Bragg-like peak, Q0/FWHM) | — | `modeling.py::DiffractionPeakPopulation`, `_diffraction_peak_intensity` | Worth reading before designing the WAXS tab's per-peak data model — may already have a compatible parameter shape to imitate, even though it lacks the true-Voigt/Debye-Waller physics needed here. |

**Net effect on scope**: the genuinely *new* physics to write is smaller
than "15 models" suggests. It is, concretely:

1. The Porod+roughness two-term background formula (eq. 3) — new, small.
2. The "algebraic globule" form factor (A3.8/A3.9), if judged not
   adequately covered by the existing spheroid — new, small, follows the
   documented form-factor recipe exactly.
3. A true Voigt peak (Lorentzian⊗Gaussian convolution, not pseudo-Voigt) —
   new, small (closed-form Voigt via the Faddeeva function, `scipy.special.wofz`
   is the standard implementation — no need to numerically convolve).
4. The `exp(−Q²⟨δz²⟩/3)` Debye-Waller-type intensity factor on the WAXS
   peak — new, trivial (one line).
5. The `1/Q²` orientation-averaging factor tying the WAXS peak's 1D
   correlation function to the powder-averaged 3D intensity — new, trivial,
   but easy to silently drop; call this out explicitly in code comments and
   tests (a fit that "looks right" without it will have a wrong `K`).
6. The `ρ_struc` (corrected) → SLD → contrast → φ/r/S_mp/w_P/w_C conversion
   layer — new, this is the "convenience" work, and it's plumbing (calling
   `scattering_contrast.py` + a handful of formulas), not new physics.
7. (v2) The crumpled-layer geometry linking SAXS eq. 5 and WAXS eq. 17
   through shared R/Σ/D — new wiring, deferred.

## 2. Proposed core architecture

`pyirena/core/carbon_fit.py`, following the `_SerialisableDataclass` pattern
from `modeling.py` (inherit it — new fields serialise for free):

```python
@dataclass
class CarbonBackground(_SerialisableDataclass):
    S_macro: float = 5.0
    porod_exponent: float = 4.0        # fittable, default fixed
    fit_porod_exponent: bool = False
    S_rough: float = 0.0
    R_rough: float = 50.0              # Å

@dataclass
class CarbonSaxsRegion(_SerialisableDataclass):
    mode: str = 'fractal'              # 'fractal' | 'teubner_strey'
    # fractal branch
    phi: float = 0.1
    pore_radius: float = 10.0          # Å
    globule_k: float = 1.0
    fractal_D: float = 2.5
    fractal_sigma: float = 100.0       # Å
    # teubner_strey branch
    ts_I0: float = 1.0
    ts_C1: float = -30.0
    ts_C2: float = 5000.0

@dataclass
class CarbonWaxsPeak(_SerialisableDataclass):
    label: str = '002'
    Q_c: float = 1.8                   # Å⁻¹
    w_L: float = 0.1
    w_G: float = 0.1
    K: float = 1.0
    delta_z2: float = 0.0              # ⟨δz²⟩, Å²
    fit_flags: dict = field(default_factory=dict)

@dataclass
class CarbonMaterial(_SerialisableDataclass):
    formula: str = 'C'
    link_geometry: bool = False        # v2: tie WAXS crumpling to SAXS fractal params

class CarbonFitModel:
    background: CarbonBackground
    saxs: CarbonSaxsRegion
    waxs_peaks: list[CarbonWaxsPeak]     # repeatable, like WAXSPeakFitModel's peak list
    material: CarbonMaterial
    def to_dict / from_dict            # composes the above, defaults for every field
    def evaluate(q) -> np.ndarray      # I_Porod + I_mp + I_waxs, calls into
                                        # simple_fits/modeling/waxs_peakfit/scattering_contrast helpers
    def compute_derived(fitted) -> dict  # the ~15-value convenience layer
```

This mirrors Modeling's "list of populations" pattern (here: background +
one SAXS region + N WAXS peaks) rather than Simple Fits' "one model" pattern,
because the tool is inherently a sum of heterogeneous components — closer in
spirit to Modeling than to Simple Fits or WAXS Peak Fit alone. It should
still get its own `gui/carbon_fit_panel.py`, not be bolted onto
`modeling_panel.py`, per Jan's explicit ask for a standalone tool with its
own tab layout and Data Selector entry.

## 3. New physics to write (detail)

- **Porod+roughness** (`_porod_roughness(q, S_macro, n, S_rough, R_rough,
  delta_sld2)`): direct transcription of 01-models.md §2, once the bracket
  power is re-verified (open question).
- **Algebraic globule form factor**: if adding, follow
  `docs/developer_adding_form_factors.md` exactly — `_build_g_globule` in
  `form_factors.py`, register in `_G_BUILDERS`, add to `FF_LABELS`/
  `FF_PARAM_LABELS` in `modeling_panel.py` (benefits Modeling too, not just
  this tool — a good "give back" side effect). Verify against Annex 3
  Fig. A3.1's `b=1.713` sphere-limit case as the ground-truth test.
- **True Voigt peak**: propose adding as a fourth shape (`"Voigt"`) to
  `pyirena/core/waxs_peakfit.py`'s existing `PEAK_SHAPES` /
  `eval_peak`/`eval_peak_derivs`/`peak_area` family, via
  `scipy.special.wofz` (`voigt_profile` is also directly available in
  `scipy.special` since a while — check the pinned scipy version in
  `pyproject.toml` supports it before adding a dependency-version bump).
  This benefits the existing WAXS Peak Fit tool too (real diffraction peaks
  are often better described by a true Voigt than a linear pseudo-Voigt
  mix) — flag this as a nice independent small PR that could land before
  the carbon tool if Jan wants to de-risk it separately.
- **Debye-Waller factor + 1/Q² orientation factor**: one-line multiplications
  onto whatever peak-shape evaluator is used; keep them as clearly separated,
  named factors in the code (not folded into the peak shape function) so
  they're independently testable and match the paper's own factored form
  (eq. 2).
- **Density/SLD conversion helper** (`pyirena/core/carbon_density.py` or a
  section of `carbon_fit.py` — decide based on size): `rho_struc(d002, d100,
  d002_graphite=..., d100_graphite=..., rho_graphite=2.26)`,
  `sld_from_density(formula, density)` (thin wrapper over
  `scattering_contrast.compute_compound`), and the φ/w_P/w_C/S_mp derived
  formulas from 01-models.md §5.

## 4. GUI layout proposal

Mirrors Modeling/Unified Fit's shell (per `docs/developer_adding_features.md`
§ "standard UX contract" — all of the shared helpers listed there apply:
`attach_table_copy`, `attach_plot_export`, `QRangeFields`, `window_state`,
`file_drop`, `report_buttons`, theme tokens):

- Left: tabbed control area — **Background**, **SAXS region**, **WAXS
  region** (peak table + add/remove peak buttons, styled like
  `waxs_peakfit_panel.py`'s own peak list), **Material** (formula + density
  inputs, read-only derived-quantity readout panel — SLD_C, SLD_sample,
  ΔSLD×2, φ, ξ, d, f_a, r, w_P, w_C, S_mp, S_part).
  Global Fit / Fit-between-cursors / MC-uncertainty controls beneath, same
  pattern as Simple Fits/Modeling.
- Right: one or two `pg.PlotItem`s via `make_sas_plot` — recommend **one**
  full-range log-log I(Q) plot with the three component curves overlaid
  (matching the paper's own Fig. 5–7 style: black data, red total fit, blue
  I_Porod, orange I_mp, green I_waxs — reuse `theme.py` tokens for consistent
  colors) plus data/residuals, rather than two separate graphs — the whole
  point of the tool is seeing all three regions on one curve. A second,
  optional zoomed WAXS-only panel (linear Q, to see peak shapes clearly) is
  a reasonable v1.1 addition, not required for v1.
- Editable Q-range via `QRangeFields`, shared cursors.
- Report/export via `report_buttons.make_report_buttons` +
  `core/reporting.py` section (see registry item 7 below).

## 5. Registry checklist (from `docs/developer_adding_features.md`)

One row added to each, key = `carbon_fit` unless noted:

| # | Registry | File | Proposed key |
|---|---|---|---|
| 1 | `TOOL_REGISTRY` | `io/schema.py` | `carbon_fit` |
| 2 | `PYIRENA_RESULT_GROUPS` | `io/_nxcansas_common.py` | derived from (1), no action |
| 3 | `TOOL_CROSS_REF` | `io/igor_names.py` | `carbon_fit`; wave names e.g. `CarbonI`/`CarbonQ`, plus per-component curves (`CarbonPorodI`, `CarbonMpI`, `CarbonWaxsI`) — decide wave-note parameter list to include all ~15 derived values |
| 4 | `TOOL_GROUP_PATH` + `TOOL_LABEL` | `gui/setup_loader.py` | `carbon_fit` |
| 5 | `_TOOL_REGISTRY` | `batch/pipeline.py` | `carbon_fit` |
| 6 | defaults block | `state/state_manager.py` | `carbon_fit` |
| 7 | `_build_report` kwarg + section | `core/reporting.py` | `carbon_fit_results` |
| 8 | `read_carbon_fit()` + result dataclass | `api/results.py`, `api/schemas.py` | `read_carbon_fit` |
| 9 | control module + `TOOL_SCHEMA_BY_NAME` | `api/control/carbon_fit.py`, `api/control/schemas.py`, `api/control/__init__.py` | yes — this should be agent-drivable like the other 5 fitting tools (Jan explicitly asked for MCP wiring) |
| 10 | `_CATEGORY_BY_MODULE` + `_CATEGORY_BLURBS` | `mcp/dispatch.py` | `carbon` |
| 11 | detection + `_collect_carbon_fit` | `gui/hdf5viewer/pyirena_readers.py` | `carbon_fit` |
| 12 | collect-item lists | `gui/hdf5viewer/plot_controls.py` | add all params + derived |
| 13 | wave names + writer | `io/igor_names.py`, `io/h5xp_extractor.py` | see (3) |
| 14 | checkbox, launcher, window, worker, results window, config dialog | `gui/data_selector/` | grep `waxs_peakfit` in `panel.py` per AGENTS.md's own advice, replicate every touch point |
| 15 | window registration | `gui/window_state.py` call site | `carbon_fit` |

Add the row to `pyirena/tests/test_tool_registration.py` **first**, per the
guide's own advice ("add your row first... that converts the list from
something you have to remember into something the build tells you"), before
writing any GUI code.

## 6. Files to create (bottom-up order)

```
pyirena/core/carbon_fit.py            math + CarbonFitModel (to_dict/from_dict)
pyirena/core/carbon_density.py        density/SLD/contrast convenience layer (or fold into above)
pyirena/io/nxcansas_carbon_fit.py     entry/carbon_fit_results
pyirena/batch/carbon_fit.py           headless fit_carbon_fit(data, config)
pyirena/gui/carbon_fit_panel.py       thin Qt panel + graph window
pyirena/api/control/carbon_fit.py     agent control surface
docs/carbon_fit_gui.md                user documentation + scripting example
pyirena/tests/test_carbon_fit.py      math vs analytic ground truth
```

Possible new/modified shared files:
```
pyirena/core/form_factors.py          + globule form factor (if not reusing spheroid)
pyirena/gui/modeling_panel.py         + FF_LABELS entry for globule (if added)
pyirena/core/waxs_peakfit.py          + true Voigt peak shape (benefits WAXS Peak Fit too)
pyirena/core/scattering_contrast.py   possibly + a rho_struc-from-d-spacings helper,
                                       or keep that carbon-specific and only call
                                       compute_compound/compute_contrast from carbon_fit.py
```

## 7. Phased build order

1. **Core + tests** — `carbon_fit.py` math, all three components, unit tests
   against analytic ground truth (synthesize `I(Q)` from known parameters
   using the formulas in 01-models.md directly, assert the fitter recovers
   them — same pattern as `validationData/generate_validation_data.py`).
   No GUI yet. This is where the nm↔Å and nomenclature decisions
   (03-open-questions.md) must be locked in, since they're baked into every
   parameter name and default.
2. **io + registries (1)(3)** — HDF5 round-trip test.
3. **batch + registry (5)** — headless fit from a config dict.
4. **gui + registries (4)(6)(14)(15)** — the panel, obeying the standard UX
   contract; this is the largest single chunk of work (Data Selector wiring
   alone is called out in the guide as "the single biggest wiring job").
5. **api + registries (7)(8)(9)(10)** — reader first, control surface once
   the interactive fitting flow is stable.
6. **HDF5 Data Explorer + Igor export (11)(12)(13)**.
7. **docs** — `docs/carbon_fit_gui.md`, row in `AGENTS.md` §3, row in
   `docs/module_map.md`, `CHANGELOG.md` entry.

At each step, `pytest pyirena/tests/test_tool_registration.py` tells you
what's still unwired — run it continuously, not just at the end.

## 8. Testing strategy

- **Analytic ground truth**: for each of the three components independently,
  generate synthetic `I(Q)` from hand-picked parameters via the exact
  formulas (no noise), fit, assert recovery to tight tolerance. This is
  achievable without Igor and without the original paper's raw data — Jan
  confirmed the Igor implementation isn't worth digging out for this reason.
- **Combined full-range fit**: synthesize all three components summed (as
  eq. 1), fit simultaneously, verify no cross-talk / parameter confusion
  between regions (e.g. background stealing intensity that belongs to the
  WAXS peak tail) — this is the actual hard numerical-fitting question this
  tool raises and is worth a dedicated test with the components' Q-ranges
  deliberately made to overlap.
- **Regression against the paper's own numbers**: Tables 1–3 give refined
  parameters for 7 real samples (graphite, CPC, GC, HC, AC, CDC400, CDC900).
  Digitizing even one or two of these (their Q-range, sample type, and
  published best-fit values) as a `validationData`-style fixture would let
  us assert "pyIrena reproduces the published fit," which is stronger
  evidence than a synthetic self-consistency test. Recommend this as a
  stretch goal for the test suite, not a blocker for v1.
- **Serialization round-trip**: standard `to_dict`/`from_dict` test per
  `test_core_serialization.py` conventions, including the WAXS peak list
  (variable length — test with 0, 1, and 3 peaks).
- **HDF5 round-trip**, **batch config path**: per the master checklist.

## 9. Explicit v1/v2 scope split (record decisions, PLAN.md style)

**v1 (first implementable slice):**
- Background: Porod (fixed exponent 4) + roughness term.
- SAXS region: fractal+globule **or** Teubner-Strey, mutually exclusive mode
  toggle, independently parameterized (no shared geometry with WAXS).
- WAXS region: arbitrary list of Voigt peaks with Debye-Waller factor and
  1/Q² orientation factor, no crumpled-layer S_3D wrapping.
- Material tab: formula + density → SLD/contrast → φ, w_P, w_C, S_mp,
  S_part, per-peak d-spacing. Corrected ρ_struc formula from the corrigendum.
- Full registry wiring per §5 (Data Selector, batch, API/MCP, Igor export).

**v2 (explicitly deferred, record here so it isn't re-litigated):**
- Shared/linked crumpled-layer geometry (R, Σ, D) tying SAXS eq. 5 and WAXS
  eq. 17 together (main text §3.4/§3.6).
- Disk/slit-pore alternative form factor for the SAXS region (paper mentions
  it as an alternative to the globule, not used in their worked examples).
- Pore-pore clustering / Teubner-Strey Fig. 4e style aggregation structure
  factor beyond the single Teixeira fractal already covered.
- Digitized paper-data regression tests (nice-to-have, not required to ship).
