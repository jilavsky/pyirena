# pyirena HDF5 Structure Reference

**Audience:** Igor Pro importer authors (and AI assistants writing such importers).
**Scope:** All HDF5 / NXcanSAS group structures produced by pyirena, with
full dataset, attribute, and sub-group inventories.

This document is the single authoritative reference for the *on-disk*
layout of every pyirena-generated file. Both the raw SAXS/USAXS/WAXS
data layout (NXcanSAS / NXsas conformance) and the per-tool result
groups are documented here.

---

## 1. Conventions used in this document

- HDF5 paths are shown in `monospace`. Group nodes end in `/`.
- An attribute on a group or dataset is written as `@name` (e.g. `@NX_class`).
- Datasets show: `name` *(dtype, shape, units)* — meaning.
- Optional members are tagged **(optional)**. Members marked **(legacy)**
  exist in older files but are no longer written by current pyirena.
- Where Igor Pro loading is non-obvious, an Igor hint is appended.
- Numeric conventions:
  - Q is in inverse angstroms (`1/A` or `1/angstrom`); both spellings
    appear in legacy files.
  - I is in `1/cm` (cm⁻¹) for SAXS/USAXS, `arb` (arbitrary) for WAXS
    when raw 2-D detector counts have not been calibrated.
  - All Q/I/Idev arrays are 1-D, stored as `float64` unless noted.
- A pyirena HDF5 file is permitted to contain *any combination* of the
  result groups in §3 — they coexist freely and are independent.

### Igor Pro loading primer

The HDF5 XOP for Igor Pro 9/10 exposes two main calls:

```
HDF5LoadFile  /R /N=fileID  fileVar         // open a file, get an ID
HDF5LoadGroup /R /Z          rootDF, fileID, "/entry/unified_fit_results"
HDF5LoadData  /Z /N=waveName fileID,        "/entry/<sample>/sasdata/Q"
```

`HDF5LoadGroup` recursively loads an entire subtree as Igor data folders +
waves. `HDF5LoadData` loads a single dataset and creates a wave with the
provided name. Neither call carries HDF5 attributes onto the wave by
default; use `HDF5LoadData /A=...` to fetch named attributes, or
`HDF5ListAttributes` followed by `HDF5LoadAttribute` to enumerate them.

For the recipes in §5, all Igor code assumes the file has been opened
with `HDF5OpenFile /R fileID as path` and closed with `HDF5CloseFile fileID`.

---

## 2. Top-level NXcanSAS layout

### 2.1 Root attributes

```
file.h5
├── @default          = "entry"
├── @file_name        (string)
├── @file_time        (string, ISO timestamp or "YYYY-MM-DD HH:MM:SS")
├── @creator          ("Matilda NeXus writer" | "pyirena")
├── @NeXus_version    "4.3.0"
├── @HDF5_version     (string)
├── @h5py_version     (string)
├── @instrument       "12IDE USAXS"        (USAXS files only)
└── @Matilda_version  (string)             (USAXS files only)
```

USAXS-reduced files written by Matilda also include `@scan_mode`,
`@timestamp`, and `@creator_config_file`.

### 2.2 entry/ group

```
entry/                              [NXentry]
├── @NX_class       = "NXentry"
├── @canSAS_class   = "SASentry"
├── @default        = "<sample_name>"      (points to the desmeared subentry)
├── definition      (string) "NXsas"
├── title           (string)              (optional)
├── experiment_identifier (string)        (optional)
└── ...
```

`entry/@default` is critical for picking the *primary* (desmeared) data
when both SMR and desmeared variants exist — see §2.3 vs §2.4.

### 2.3 entry/{sample_name}/ — desmeared (USAXS) or only (SAXS/WAXS)

```
entry/{sample_name}/                [NXsubentry]
├── @NX_class       = "NXsubentry"
├── @canSAS_class   = "SASentry"
├── @default        = "sasdata"
├── @title          (string)
├── definition      (string) "NXcanSAS"
├── title           (string)
├── run             (string) "run_identifier"
└── sasdata/                        [NXdata]
    ├── @NX_class       = "NXdata"
    ├── @canSAS_class   = "SASdata"
    ├── @signal         = "I"
    ├── @I_axes         = "Q"
    │
    ├── I            (float64, Nq, units="1/cm")
    │   ├── @units            "1/cm"
    │   ├── @uncertainties    "Idev"
    │   ├── @long_name        "Intensity"
    │   ├── @blankname        (string, USAXS only)
    │   ├── @thickness        (float, USAXS only — sample thickness in mm)
    │   ├── @label            (string, USAXS only)
    │   ├── @Kfactor          (float, USAXS only — calibration constant)
    │   └── @OmegaFactor      (float, USAXS only — solid angle factor)
    │
    ├── Q            (float64, Nq, units="1/angstrom")
    │   ├── @units            "1/angstrom"
    │   ├── @long_name        "Q (A^-1)"
    │   └── @resolutions      "Qdev"
    │
    ├── Idev         (float64, Nq)
    │   ├── @units            "1/cm" or "cm2/cm3"
    │   └── @long_name        "Uncertainties"
    │
    └── Qdev         (float64, Nq)         (Q-resolution; desmeared)
        ├── @units            "1/angstrom"
        └── @long_name        "Q (A^-1)"
```

For **SAXS** and **WAXS** files there is normally only this single
subentry — no SMR variant exists. For **USAXS** the desmeared
subentry coexists with the slit-smeared one (§2.4).

**Igor hint:** load the whole sasdata group with
`HDF5LoadGroup /R /Z rootDF, fileID, "/entry/<sample>/sasdata"`
which produces waves `Q`, `I`, `Idev`, `Qdev` in the destination data
folder.

### 2.4 entry/{sample_name}_SMR/ — slit-smeared USAXS variant (USAXS only)

Identical structure to §2.3 except:

```
entry/{sample_name}_SMR/sasdata/
    ├── @signal     = "I"
    ├── @I_axes     = "Q"
    ├── I, Q, Idev as in §2.3 (different numerical values)
    ├── Q.@resolutions = "dQw,dQl"
    ├── dQw          (float64, Nq, units="1/angstrom")
    │                 — point-wise Q widths
    └── dQl          (float64, scalar, units="1/angstrom")
                      — slit length (constant across spectrum)
```

**Distinguishing desmeared from SMR:** the SMR group's name ends in
`_SMR`; its `Q.@resolutions` attribute contains `dQw,dQl` (two
comma-separated names) instead of `Qdev`. `entry/@default` always
points to the desmeared variant.

For Igor importers planning to operate on USAXS data: prefer the
desmeared variant; the slit-smeared one is rarely needed downstream.

### 2.5 entry/QRS_data/ — USAXS reduction debug arrays (USAXS only, optional)

Raw Q/R/S arrays from USAXS reduction, retained for debugging. Not
used by any analysis tool.

```
entry/QRS_data/
├── Intensity   (float64, N_raw, units="arb", "Intensity")
├── Q           (float64, N_raw, units="1/angstrom")
└── Error       (float64, N_raw, units="arb")
```

Length `N_raw` (typically 7999 in current Matilda output) is *not*
the same as the reduced `Nq` in §2.3 / §2.4. There is no NX_class
attribute; treat as plain HDF5 group.

### 2.6 entry/Blank_data/ — blank measurement (USAXS only, optional)

Same shape as §2.5, plus an `Intensity.@blankname` attribute naming
the blank file.

### 2.7 Sample / Instrument / Metadata sub-trees

These groups are populated by the upstream data acquisition pipeline
(Matilda for USAXS at APS 12-ID-E; varies for SAXS/WAXS). They are
*read* by pyirena's metadata loader but never *written* by pyirena
itself. Their schema mirrors the NeXus standard.

```
entry/instrument/                   [NXinstrument]
├── @NX_class      = "NXinstrument"
├── @canSAS_class  = "SASinstrument"
├── name           (string)
├── monochromator/                  [NXmonochromator]
│   ├── @NX_class      = "NXmonochromator"
│   ├── energy         (float, 1, units="keV")
│   └── wavelength     (float, 1, units="A")
└── source/                         [NXsource] (optional)

entry/sample/                       [NXsample]
├── @NX_class      = "NXsample"
├── @canSAS_class  = "SASsample"
├── name           (string)
├── thickness      (float, 1, units="mm")
├── concentration  (float, 1)
├── temperature    (float, 1)
├── description    (uint8 array, often empty)
├── chemical_formula (uint8 array)
└── ...

entry/metadata/                     [NXcollection] (USAXS only)
└── ... 100+ scalars: DCM_energy, DCM_theta, AR_*, upd_bkg*, etc.

entry/Metadata/                     [NXcollection] (SAXS/WAXS only,
                                     note capital M)
└── ... different key set per beamline
```

For an Igor importer aiming only at the **scientific data + analysis
results**, these are useful for header annotation but not required
for plotting.

---

## 3. Per-tool result groups

Every pyirena fitting / analysis tool writes its results into a single
group under `entry/`. All such groups carry `@NX_class = "NXprocess"`
and an `@analysis_type` tag identifying the tool.

### 3.1 entry/unified_fit_results/ — Unified Fit (Beaucage)

Source: `pyirena/io/nxcansas_unified.py`.

Detailed reference document: `docs/NXcanSAS_UnifiedFit_Format.md`.
Summary tree:

```
entry/unified_fit_results/          [NXprocess]
├── @NX_class       = "NXprocess"
├── @analysis_type  = "Unified Fit"
├── @program        = "pyirena"
├── @timestamp      (ISO string)
├── @num_levels     (int)
│
├── background      (float64, scalar)            1/cm
├── chi_squared     (float64, scalar)
├── background_err  (float64, scalar, optional)  MC σ
│
├── Q                 (float64, Nq, units="1/angstrom")
├── intensity_data    (float64, Nq, units="1/cm")
├── intensity_error   (float64, Nq, units="1/cm", optional)
├── intensity_model   (float64, Nq, units="1/cm")
├── residuals         (float64, Nq)              normalized
│
└── level_1/, level_2/, … level_N/   [group]
    ├── @level_number   (int)
    ├── G               (float64, scalar)        Guinier prefactor
    ├── Rg              (float64, scalar)        Å
    ├── B               (float64, scalar)        Porod constant
    ├── P               (float64, scalar)        Porod exponent
    ├── ETA             (float64, scalar)        Å (correlated only)
    ├── PACK            (float64, scalar)        (correlated only)
    ├── RgCutoff        (float64, scalar, optional)
    ├── correlations    (bool, scalar, optional)
    ├── G_err, Rg_err,
    │   B_err, P_err,
    │   ETA_err, PACK_err  (float64 scalars, optional — MC σ, only
    │                       written when > 0)
    ├── @Sv             (float, optional)        m²/cm³ (when P ≈ 4)
    └── @Invariant      (float, optional)        cm⁻⁴
```

**Backward compatibility:** older files (pyirena ≤ 0.4.5) stored
`background`, `chi_squared`, and the per-level parameters as
**attributes** on the group/level rather than scalar datasets. New
loaders try the dataset first then fall back to the attribute. New
writers always emit datasets.

**Igor hint:**
```igor
HDF5LoadGroup /R /Z rootDF, fileID, "/entry/unified_fit_results"
// produces waves Q, intensity_data, intensity_model, residuals
// plus subfolders level_1, level_2, ... each containing scalar waves G, Rg, B, P, ...
```

### 3.2 entry/simple_fit_results/ — Simple Fits

Source: `pyirena/io/nxcansas_simple_fits.py`.

Single-model fits (Guinier, Porod, Sphere, Treubner-Strey, Guinier
Sheet, Unified Born-Green, Debye Polymer Chain, etc.). Only one model
per file at a time — saving a new fit replaces the previous group.

```
entry/simple_fit_results/           [NXprocess]
├── @NX_class            = "NXprocess"
├── @program             = "pyirena.core.simple_fits"
├── @timestamp           (ISO string)
├── @model               (string, e.g. "Guinier", "Sphere", …)
├── @success             (bool)
├── @dof                 (int)
├── @q_min, @q_max       (float)
├── @use_complex_bg      (bool)
├── @n_mc_runs           (int)
│
├── chi_squared          (float64, scalar)
├── reduced_chi_squared  (float64, scalar)
│
├── Q                    (float64, Nq, units="1/angstrom", gzip)
├── I_model              (float64, Nq, units="arb", gzip)
├── residuals            (float64, Nq, gzip)
├── intensity_data       (float64, Nq, units="1/cm", optional, gzip)
├── intensity_error      (float64, Nq, units="1/cm", optional, gzip)
│
├── params/              [group]
│   ├── I0, Rg, B, P, …  (float64 scalars; key set varies by model)
│   │   ├── @limit_low   (optional)
│   │   └── @limit_high  (optional)
│   └── ...
│
├── params_std/          [group]
│   └── <param_name>     (float64 scalar — std from covariance or MC)
│
└── derived/             [group, optional]
    ├── Thickness        Å         (Guinier Sheet)
    ├── CorrLength       Å         (Treubner-Strey)
    ├── RepeatDist       Å         (Treubner-Strey)
    ├── Rad              Å         (Unified Born-Green)
    ├── G1                          (Unified Born-Green)
    ├── Rg2              Å         (Unified Born-Green)
    └── ... (per model)
```

**Backward compatibility:** `chi_squared`, `reduced_chi_squared` may
be group attributes (legacy) instead of datasets.

**Parameter set per model** (informative — enumerated in
`pyirena/core/simple_fits.py::MODEL_REGISTRY`):

| Model               | Free parameters                    |
|---------------------|------------------------------------|
| Guinier             | I0, Rg                             |
| Porod               | B, P, BG_G, BG_P                   |
| Guinier Sheet       | I0, Thickness                      |
| Sphere              | I0, R, BG_G, BG_P                  |
| Treubner-Strey      | I0, ksi (a), CorrLength, BG_G, BG_P|
| Unified Born-Green  | I0, Rg, ETA, PACK, BG_G, BG_P      |
| Debye Polymer Chain | I0, Rg, BG_G, BG_P                 |
| ...                 | (see source)                        |

### 3.3 entry/sizes_results/ — Size Distribution

Source: `pyirena/io/nxcansas_sizes.py`.

Inversion of I(Q) to a particle size distribution P(r). Three
solver methods: MaxEnt, regularization, TNNLS, and Monte Carlo (the
last writes per-bin uncertainties).

```
entry/sizes_results/                [NXprocess]
├── @NX_class               = "NXprocess"
├── @program                = "pyirena.core.sizes"
├── @timestamp              (ISO string)
│
├── Datasets (fit results — new format):
│   ├── chi_squared         (float64, scalar)
│   ├── volume_fraction     (float64, scalar)         total Vf integrated
│   ├── rg                  (float64, scalar)         Å (computed from P(r))
│   ├── n_iterations        (int64,  scalar)
│   └── q_power             (float64, scalar)
│
├── Attributes (model / grid setup):
│   @shape                  (string) "Spheres" | "Cylinders" | …
│   @contrast               (float)  (Δρ)² in 10²⁰ cm⁻⁴
│   @aspect_ratio           (float)
│   @r_min, @r_max          (float)  Å
│   @n_bins                 (int)
│   @log_spacing            (bool)
│   @background             (float)  1/cm
│   @power_law_B, @power_law_P  (float)
│   @method                 (string) "maxent" | "regularization" | "tnnls" | "montecarlo"
│   @error_scale            (float)
│
├── Attributes (method-specific):
│   @maxent_sky_background, @maxent_stability, @maxent_max_iter
│   @regularization_evalue, @regularization_min_ratio
│   @tnnls_approach_param,  @tnnls_max_iter
│   @montecarlo_n_repetitions, @montecarlo_convergence, @montecarlo_max_iter
│
├── Attributes (Q ranges used):
│   @power_law_q_min, @power_law_q_max
│   @background_q_min, @background_q_max
│   @cursor_q_min,    @cursor_q_max
│
├── Arrays (gzip-compressed):
│   ├── Q                   (float64, Nq, units="1/angstrom")
│   ├── intensity_data      (float64, Nq, units="1/cm")
│   ├── intensity_model     (float64, Nq, units="1/cm")
│   ├── residuals           (float64, Nq)
│   ├── intensity_error     (float64, Nq, optional)
│   ├── r_grid              (float64, Nbins, units="angstrom")
│   ├── distribution        (float64, Nbins, units="volume_fraction/angstrom")
│   │                        — volume size distribution P(r)
│   ├── number_dist         (float64, Nbins, units="1/angstrom")
│   │                        — derived: P(r) / [(4/3)πr³]
│   ├── cumul_vol_dist      (float64, Nbins, units="volume_fraction")
│   │                        — derived: cumulative volume integral
│   ├── cumul_num_dist      (float64, Nbins, units="dimensionless")
│   │                        — derived: cumulative number integral
│   └── distribution_std    (float64, Nbins, optional, units="volume_fraction/angstrom")
│                            — per-bin std across MC repetitions (Monte Carlo only)
```

**Backward compatibility:**
- `chi_squared`, `volume_fraction`, `rg`, `n_iterations`, `q_power`
  may be group attributes (legacy) instead of datasets.
- Old files used `mcsas_*` attribute names instead of `montecarlo_*`,
  and `method='mcsas'` instead of `'montecarlo'`. Importers should
  remap.
- `number_dist`, `cumul_vol_dist`, `cumul_num_dist` are derivable
  from `r_grid` and `distribution`; they are written for convenience
  but importers can recompute.

### 3.4 entry/waxs_peakfit_results/ — WAXS Peak Fit

Source: `pyirena/io/nxcansas_waxs_peakfit.py`.

Multi-peak fits (Gauss, Lorentz, pseudo-Voigt) on top of a polynomial
background. Number of peaks is variable per file.

```
entry/waxs_peakfit_results/         [NXprocess]
├── @NX_class               = "NXprocess"
├── @program                = "pyIrena waxs_peakfit"
├── @timestamp              (ISO string)
├── @n_peaks                (int)
├── @bg_shape               (string) "Constant"|"Linear"|"Quadratic"|"Cubic"
├── @dof                    (int)
├── @q_min, @q_max          (float)
│
├── chi_squared             (float64, scalar)
├── reduced_chi_squared     (float64, scalar)
│
├── Q                       (float64, Nq, units="1/angstrom")
├── I_fit                   (float64, Nq, units="arb")    total model
├── I_bg                    (float64, Nq, units="arb")    background only
├── residuals               (float64, Nq)
├── intensity_data          (float64, Nq, units="arb", optional)
├── intensity_error         (float64, Nq, units="arb", optional)
│
├── background/             [group]
│   └── bg0, bg1, … bgM     (float64 scalars — polynomial coefficients)
│       ├── @limit_low
│       └── @limit_high
│
├── background_std/         [group]
│   └── bg0, bg1, … bgM     (float64 scalars — std per coefficient)
│
└── peak_01/, peak_02/, … peak_NN/  [group]
    ├── @shape              (string) "Gauss"|"Lorentz"|"PseudoVoigt"
    ├── Q_peak              (float64, K_p, units="1/angstrom")
    │                        — Q range ±5×FWHM around Q0
    ├── I_peak              (float64, K_p, units="arb")
    │                        — individual peak curve evaluated on Q_peak
    ├── params/             [group]
    │   ├── A               (float64, scalar) amplitude
    │   ├── Q0              (float64, scalar) peak centre, 1/A
    │   ├── FWHM            (float64, scalar) full-width-half-max, 1/A
    │   └── eta             (float64, scalar, PseudoVoigt only)
    │       ├── @limit_low
    │       └── @limit_high
    └── params_std/         [group]
        └── A, Q0, FWHM, [eta]   (float64 scalars — std per param)
```

**Backward compatibility:** `chi_squared`, `reduced_chi_squared` may
be group attributes (legacy) instead of datasets.

The per-peak `Q_peak`/`I_peak` arrays are *additional* convenience
storage — they are not required for re-evaluating the model since
the parameters themselves suffice; importers can ignore them.

### 3.5 entry/modeling_results/ — Modeling (parametric forward)

Source: `pyirena/io/nxcansas_modeling.py`.

Multi-population parametric forward modelling. **Three population
types** can coexist in the same fit, each with a *different* schema
inside its `pop_NN/` group: `size_dist`, `unified_level`,
`diffraction_peak`. Importers must dispatch on the `@pop_type`
attribute when reading per-population fields.

```
entry/modeling_results/             [NXprocess]
├── @NX_class               = "NXprocess"
├── @analysis_type          = "Modeling"
├── @program                = "pyirena"
├── @timestamp              (ISO string)
│
├── chi_squared             (float64, scalar)
├── reduced_chi_squared     (float64, scalar)
├── dof                     (int64,  scalar)
├── background              (float64, scalar)        1/cm
├── q_min, q_max            (float64, scalar)        1/A
├── background_err          (float64, scalar, optional — MC σ)
│
├── model_q                 (float64, Nq, units="1/angstrom")
├── model_I                 (float64, Nq, units="1/cm")    total model
│
└── pop_01/, pop_02/, … pop_NN/  [group]
    ├── @population_index   (int)
    ├── @enabled            (bool)
    ├── @pop_type           (string) one of: "size_dist" | "unified_level"
    │                                       | "diffraction_peak"
    ├── @label              (string)
    │
    ├── model_I             (float64, Nq, units="1/cm")
    │                        — this population's I(Q) contribution
    │
    ├── derived/            [group]
    │   ├── volume_fraction      (size_dist only)
    │   ├── vol_mean_r           (size_dist only)
    │   ├── r_total_mean         (size_dist only)
    │   └── ... (per pop_type)
    │
    └── (type-specific sub-tree — see below)
```

**Per-population schema by `@pop_type`:**

`@pop_type = "unified_level"`:
```
pop_NN/
├── G, Rg, B, P, RgCO       (float64 scalars)
│   ├── @fit                (bool)
│   ├── @limit_lo
│   └── @limit_hi
├── correlations            (bool, scalar)
├── ETA, PACK               (float64 scalars, with same @fit/@limit_*)
└── pop1_<param>_err        (float64 scalars, optional — MC σ)
```

`@pop_type = "diffraction_peak"`:
```
pop_NN/
├── @peak_type              (string) "gaussian"|"lorentzian"|"pseudo_voigt"
├── position                (float64, scalar, 1/A)
├── amplitude               (float64, scalar)
├── width                   (float64, scalar, 1/A — FWHM)
└── eta_voigt               (float64, scalar — for pseudo_voigt)
    ├── @fit                (bool)
    ├── @limit_lo
    └── @limit_hi
```

`@pop_type = "size_dist"`:
```
pop_NN/
├── @dist_type              (string) "gaussian"|"lognormal"|"lsw"|
│                                    "schulz_zimm"|"ardell"
├── @form_factor            (string)
├── @structure_factor       (string) "none"|"interferences"|"hard_sphere"
├── scale, contrast         (float64 scalars; with @fit attribute)
├── use_number_dist         (bool, scalar)
├── n_bins                  (int, scalar)
├── radius_grid             (float64, n_bins, units="angstrom", optional)
├── volume_dist             (float64, n_bins, optional)
├── number_dist             (float64, n_bins, optional)
│
├── dist_params/            [group]
│   └── <param>             (float64 scalar; e.g. mean, sigma)
│       ├── @fit            (bool)
│       ├── @limit_lo
│       └── @limit_hi
├── ff_params/              [group]   form-factor parameters
│   └── <param>             (float64 scalar with @fit)
└── sf_params/              [group]   structure-factor parameters
    └── <param>             (float64 scalar with @fit)
```

**Igor hint:** because the per-population schema is type-dependent,
loading `pop_NN/` with `HDF5LoadGroup /R /Z` is safe (it just
materialises every dataset and group present), but the importer must
read `@pop_type` first to know which fields to expect and how to
interpret them.

### 3.6 entry/saxs_morph_results/ — SAXS Morph (3D voxelgram)

Source: `pyirena/io/nxcansas_saxs_morph.py`.

Generates a 3D binary voxelgram of a two-phase porous structure from
experimental I(Q) using the Gaussian Random Fields method (Berk 1991,
Roberts 1997, Levitz 2007), then refines parameters by fitting the
model I(Q) — recomputed from the voxelgram — back to the data.

The voxelgram itself is the heaviest payload: a 3-D `uint8` cube stored
gzip-compressed and chunked `(N, N, 1)` so a single 2-D slice can be
loaded without inflating the whole array.

```
entry/saxs_morph_results/                [NXprocess]
├── @NX_class               = "NXprocess"
├── @analysis_type          = "SAXS Morph"
├── @program                = "pyirena"
├── @timestamp              (ISO string)
│
├── chi_squared             (float64, scalar)
├── reduced_chi_squared     (float64, scalar)
├── dof                     (int64,  scalar)
│
├── q_min, q_max            (float64, scalar)        1/A
│
├── volume_fraction         (float64, scalar)        target φ
│   ├── @fit                (bool)
│   ├── @limit_lo
│   └── @limit_hi
├── contrast                (float64, scalar)        (Δρ)² in 10²⁰ cm⁻⁴
│   ├── @fit
│   ├── @limit_lo
│   └── @limit_hi
├── link_phi_contrast       (bool, scalar)
├── power_law_B             (float64, scalar)        with @fit / @limit_*
├── power_law_P             (float64, scalar)        with @fit / @limit_*
├── background              (float64, scalar)        with @fit / @limit_*
│
├── voxel_size              (int64,  scalar)         N (cube side)
├── box_size_A              (float64, scalar)        physical box edge [Å]
├── voxel_pitch_A           (float64, scalar)        = box_size_A / N
├── phi_actual              (float64, scalar)        realised φ of voxelgram
├── rng_seed                (int64,  scalar)         seed used (reproducibility)
│
├── data_q                  (float64, Nq, units="1/angstrom")
├── data_I                  (float64, Nq, units="1/cm")        raw data
├── data_dI                 (float64, Nq, units="1/cm")        σ_I
├── data_I_corr             (float64, Nq, units="1/cm")        data − background
├── model_q                 (float64, Nq, units="1/angstrom")
├── model_I                 (float64, Nq, units="1/cm")        full model
│
├── r_grid                  (float64, Nr, units="angstrom")
├── gamma_r                 (float64, Nr, dimensionless)       Debye autocorr
├── spectral_k              (float64, Nk, units="1/angstrom")
├── spectral_F              (float64, Nk, units="angstrom**3") spectral density
│
├── voxelgram               (uint8, (N, N, N))
│   ├── @pitch_A
│   ├── @box_size_A
│   ├── @phi_actual
│   ├── @description        "Binary phase indicator: 0=phase A, 1=phase B"
│   ├── chunks              (N, N, 1)         → cheap 2-D slice loads
│   └── compression         "gzip", level 4
│
└── <param>_err             (float64 scalar, optional — MC σ)
    e.g. volume_fraction_err, contrast_err, power_law_B_err, etc.
```

**Sizing**: a 256³ binary cube compresses to 1–10 MB on porous-media data
(chunked gzip is highly effective on smooth fields); a 512³ cube to
10–100 MB. Memory cost in RAM is the uncompressed size: 16 MB at 256³,
125 MB at 512³.

**Igor hint:** `HDF5LoadData /TYPE=1 /Z /N=voxelgram fileID, "/entry/saxs_morph_results/voxelgram"`
loads the cube as an Igor 3-D wave (the HDF5 XOP transparently
decompresses gzip-chunked datasets). For partial loads — e.g. one
slice — use `HDF5LoadData /SLAB=...` to read only the chunk you need.

### 3.7 entry/data_merge_results/ — Data Merge (provenance only)

Source: `pyirena/io/nxcansas_data_merge.py`.

Provenance for a Data Merge operation. The merged Q/I/Idev/Qdev
arrays are written *into the standard `entry/{sample}/sasdata/` group
of §2.3* — this group only records *how* the merge was done.

```
entry/data_merge_results/           [NXprocess]
├── @NX_class       = "NXprocess"
├── @program        = "pyirena"
├── @version        = "1.0"
├── @timestamp      (ISO string)
├── @ds1_path       (string — legacy attribute)
├── @ds2_path       (string — legacy attribute)
│
├── ds1_file        (string scalar)         input dataset 1 path
├── ds2_file        (string scalar)         input dataset 2 path
│
├── scale                   (float64 scalar)
├── q_shift                 (float64 scalar)
├── background              (float64 scalar)         from DS1 only
├── chi_squared             (float64 scalar)
├── q_overlap_min           (float64 scalar)
├── q_overlap_max           (float64 scalar)
├── scale_dataset           (int scalar)             1 or 2
├── fit_scale               (int scalar — bool 0/1)
├── fit_qshift              (int scalar — bool 0/1)
└── split_at_left_cursor    (int scalar — bool 0/1)
```

Booleans are stored as `int(value)`; None values as `float('nan')`.

**Strip list (informative):** when Data Merge copies the source NXcanSAS
file, it deletes any of the following result groups before writing the
merged data:

```
entry/unified_fit_results, entry/sizes_results, entry/simple_fits_results,
entry/waxs_peakfit_results, entry/data_merge_results, entry/modeling_results,
entry/data_manipulation_results
```

Importers should not assume these groups persist across a merge.

### 3.8 entry/data_manipulation_results/ — Data Manipulation (provenance only)

Source: `pyirena/io/nxcansas_data_manipulation.py`.

Provenance for Scale / Trim / Rebin / Average / Subtract / Divide
operations. Same write-into-sasdata + provenance pattern as §3.6.

```
entry/data_manipulation_results/    [NXprocess]
├── @NX_class       = "NXprocess"
├── @analysis_type  = "Data Manipulation"
├── @program        = "pyirena"
├── @version        = "1.0"
├── @timestamp      (ISO string)
│
├── operation       (string)        one of "scaled" | "trimmed" |
│                                    "rebinned" | "avg" | "sub" | "div"
├── source_file     (string)        input file path
│
└── (operation-specific scalars; varies by tool)
```

The operation-specific keys depend on the tool — for example, `scale`
for Scale, `q_min` / `q_max` for Trim, `n_bins` / `log_spacing` for
Rebin, `denominator_file` for Divide, `buffer_file` for Subtract,
etc. Importers should treat the contents of this group as a free-form
string-keyed map.

Same strip-list as §3.6 applies before re-saving.

---

## 4. How tools coexist in a single file

A typical workflow file may end up with:

```
file.h5
└── entry/
    ├── {sample}/sasdata/             ← raw or merged Q/I/Idev/Qdev
    ├── {sample}_SMR/sasdata/         ← USAXS slit-smeared variant (USAXS only)
    ├── instrument/, sample/, metadata/  ← metadata sub-trees
    │
    ├── unified_fit_results/          ← from Unified Fit tool
    ├── simple_fit_results/           ← from Simple Fits tool
    ├── sizes_results/                ← from Size Distribution tool
    ├── waxs_peakfit_results/         ← from WAXS Peak Fit tool
    ├── modeling_results/             ← from Modeling tool
    ├── data_merge_results/           ← provenance, if file came from Merge
    └── data_manipulation_results/    ← provenance, if file went through Manipulation
```

Tools that mutate the primary data (Data Merge, Data Manipulation)
strip prior result groups before writing — see the strip list in §3.6.

Tools that only compute results (Unified Fit, Simple Fits, Sizes,
WAXS Peak Fit, Modeling) overwrite their *own* group when re-run but
leave other tools' results intact.

An importer that wants to enumerate *all* available analysis results
in a file should simply walk `entry/*` and look for groups with
`@NX_class = "NXprocess"`.

---

## 5. Reading recipes

### 5.1 Python (h5py)

```python
import h5py

with h5py.File("mydata.h5", "r") as f:
    # Primary data — desmeared if both variants exist
    default = f["entry"].attrs["default"]            # "<sample_name>"
    sd = f[f"entry/{default}/sasdata"]
    Q  = sd["Q"][:]                                  # (Nq,)
    I  = sd["I"][:]
    dI = sd["Idev"][:] if "Idev" in sd else None

    # Unified Fit results
    if "entry/unified_fit_results" in f:
        uf = f["entry/unified_fit_results"]
        n  = int(uf.attrs["num_levels"])
        for k in range(1, n + 1):
            lev = uf[f"level_{k}"]
            G, Rg, B, P = (float(lev[p][()]) for p in ("G", "Rg", "B", "P"))
            print(f"L{k}: G={G:.3e}, Rg={Rg:.1f} Å, B={B:.3e}, P={P:.2f}")
```

The pyirena library provides `load_*_results` helpers in
`pyirena.io.nxcansas_*` that wrap each group. See those modules'
docstrings for the returned-dict schema.

### 5.2 Igor Pro 9/10

```igor
Function LoadPyirenaPrimary(path, dfName)
    String path, dfName
    Variable fileID
    HDF5OpenFile /R fileID as path

    // Read the entry/@default attribute to learn the sample name
    String sampleName
    HDF5LoadAttribute /Z /TYPE=2 fileID, "/entry", "default"
    Wave/T attrW = root:default
    sampleName = attrW[0]
    KillWaves attrW

    // Build destination data folder and load the sasdata subtree
    NewDataFolder /O /S $("root:" + dfName)
    HDF5LoadGroup /R /Z :, fileID, "/entry/" + sampleName + "/sasdata"
    // Now this DF contains waves Q, I, Idev, Qdev

    HDF5CloseFile fileID
End

Function LoadUnifiedFitResults(path, dfName)
    String path, dfName
    Variable fileID
    HDF5OpenFile /R fileID as path

    NewDataFolder /O /S $("root:" + dfName)
    HDF5LoadGroup /R /Z :, fileID, "/entry/unified_fit_results"
    // Top-level waves: Q, intensity_data, intensity_model, residuals,
    //                  background, chi_squared
    // Subfolders: level_1, level_2, ... each containing scalar waves
    //             G, Rg, B, P (and optionally ETA, PACK, *_err)

    HDF5CloseFile fileID
End
```

For the other result groups, follow the same pattern: open, walk to
`/entry/<group_name>`, `HDF5LoadGroup` recursively. Then read scalar
waves with `WaveVal = waveName[0]` for fit parameters.

For per-peak (WAXS) and per-population (Modeling) groups, after
loading you will find subfolders `peak_01`, `peak_02`, … or `pop_01`,
`pop_02`, … — iterate over them with `DataFolderRefStatus` /
`GetIndexedObjName`.

---

## 6. Versioning and backward compatibility

Two recurring schema changes affect importers reading old files:

1. **Float scalars: dataset (new) vs attribute (legacy).** Several
   tools moved `chi_squared`, `reduced_chi_squared`,
   `volume_fraction`, etc. from group attributes to scalar datasets so
   that they are visible in HDF5 viewers without expanding "Show
   attributes". Importers should try the dataset first then fall back
   to the attribute — pyirena's own loaders do this. Affected tools:
   Unified Fit (`background`, `chi_squared`), Simple Fits
   (`chi_squared`, `reduced_chi_squared`), Sizes (`chi_squared`,
   `volume_fraction`, `rg`, `n_iterations`, `q_power`), WAXS Peak Fit
   (`chi_squared`, `reduced_chi_squared`).

2. **Sizes: `mcsas_*` → `montecarlo_*` rename.** Old Size Distribution
   files used `method = "mcsas"` and `mcsas_n_repetitions`,
   `mcsas_convergence`, `mcsas_max_iter`. Current files use
   `method = "montecarlo"` and `montecarlo_*`. Importers should remap
   the old names.

The desmeared-vs-SMR detection (§2.3 vs §2.4) has been stable since
the original Matilda format and is not expected to change.

---

## 7. Appendix: Glossary of acronyms

| Acronym  | Tool           | Meaning                                           |
|----------|----------------|---------------------------------------------------|
| G        | Unified Fit    | Guinier prefactor (intensity at Q→0 for level)    |
| Rg       | Unified Fit    | Radius of gyration (Å)                            |
| B        | Unified Fit    | Porod constant                                    |
| P        | Unified Fit    | Porod exponent                                    |
| RgCO     | Modeling/UF    | Cut-off Rg for next-larger structure (Å)          |
| ETA      | Unified/Simple | Correlation hole size (Å)                         |
| PACK     | Unified/Simple | Packing factor (degree of correlation, 0–16)      |
| Sv       | Unified Fit    | Surface to volume ratio (m²/cm³, derived if P≈4)  |
| dist_type| Modeling       | gaussian / lognormal / lsw / schulz_zimm / ardell |
| Vf       | Modeling       | Volume fraction = 0.5(1 − √(1 − 4·scale))         |
| Q0       | WAXS Peak Fit  | Peak centre Q (1/Å)                               |
| FWHM     | WAXS Peak Fit  | Full-width at half-maximum (1/Å)                  |
| eta      | WAXS Peak Fit  | Pseudo-Voigt mixing parameter (0=Gauss, 1=Lorentz)|
| Kfactor  | USAXS metadata | Calibration constant (instrument-specific)        |
| OmegaFactor| USAXS metadata | Solid-angle factor                              |
| dQw      | USAXS SMR      | Point-wise Q widths (1/Å)                         |
| dQl      | USAXS SMR      | Slit length, constant per spectrum (1/Å)          |
| Qdev     | NXcanSAS       | Q-resolution standard deviation (desmeared)       |

---

## See also

- `docs/NXcanSAS_UnifiedFit_Format.md` — extended Unified Fit format
  reference with worked examples.
- `pyirena/io/nxcansas_*.py` — authoritative source code for every
  result group writer.
- `pyirena/io/hdf5.py` — `readGenericNXcanSAS()` and
  `find_matching_groups()` for generic NXcanSAS reading.
- NeXus standard — http://www.nexusformat.org/
- NXcanSAS standard — https://www.cansas.org/formats/canSAS2012/1.1/doc/index.html
