# Validating pyIrena, and comparing it with Igor Pro Irena

pyIrena is a port of the Igor Pro **Irena** package, so two questions follow it
everywhere: *is the mathematics right*, and *does it agree with Irena?*  The
`validationData/` folder answers both with the same material — a set of
synthetic scattering curves whose parameters are known exactly.

## Why synthetic data

Comparing two packages on experimental data tells you whether they agree, not
whether either is right; and when they disagree it gives you no way to decide
which one is wrong.  Data generated from known parameters gives an external
reference both packages are measured against, so "pyIrena and Irena agree" and
"both are correct" become separate, checkable statements.

The generator matters as much as the data.  `validationData/_models.py`
implements every model — sphere and spheroid form factors, the Beaucage Unified
equation, Percus-Yevick hard spheres, Teixeira's mass and surface fractals,
Hammouda's Guinier-Porod form, the Simple-Fit laws, the peak shapes and the
Lake slit-smearing integral — directly from the published expressions, and
**imports nothing from pyIrena**.  Fitting the resulting data therefore tests
pyIrena's mathematics; if pyIrena had a sign error or a wrong prefactor, the
recovered parameters would be wrong.  A generator that called
`pyirena.core` would only ever prove that pyIrena can invert its own forward
model.

## What is in `validationData/`

| | |
|---|---|
| `README.md` | the manifest: every dataset with its complete ground truth, generated from the same source as the data |
| `ground_truth.json` | the same, machine readable |
| `_models.py` | independent model implementations (no pyIrena import) |
| `_spec.py` | the dataset definitions: parameters, Q ranges, derived quantities |
| `generate_validation_data.py` | writes the data files, the manifest and the JSON |
| `run_validation_report.py` | fits every file with pyIrena, writes the results tables |
| `VALIDATION_RESULTS.md` / `.csv` | truth vs pyIrena for every parameter, with empty Irena columns |
| `sizes/ unified/ modeling/ simple_fits/ waxs/ merge/ manipulation/` | the data |

Each dataset ships as three files:

* `<name>.dat` — 3-column ASCII `Q  I  dI`, the noisy "measurement", with the
  full ground truth in the `#` header.  This is the file to load in Irena.
* `<name>_ideal.dat` — 2-column `Q  I`, the exact noise-free model.  Comparing
  a package's *calculated* curve with this needs no fitting at all and isolates
  the forward model from the optimiser.
* `<name>.h5` — NXcanSAS, the noisy data under `entry/sasdata`, plus
  `entry/ground_truth/I_ideal` and `entry/ground_truth/parameters_json`.

## Regenerating and re-running

```bash
python3 validationData/generate_validation_data.py   # rebuild the data + manifest
python3 validationData/run_validation_report.py      # fit everything, write the tables
```

Generation is deterministic: the noise seed for each dataset is
`SHA-256(dataset name)[:8]`, so the files rebuild byte-for-byte and adding a
dataset never perturbs an existing one.  Every fit starts from values offset
from the truth (0.6x, typically) so the tables demonstrate convergence rather
than assume it.

## Comparing against Irena

1. Load `<name>.dat` in Igor Pro Irena (three columns, Q in Å⁻¹, I in cm⁻¹).
2. Use the tool settings quoted for that dataset in `VALIDATION_RESULTS.md` —
   the same size grid, the same number of bins, the same fixed parameters — so
   both packages solve the same problem.
3. Enter the Irena values in the `Irena` column of `VALIDATION_RESULTS.csv`.

The result is a three-way table: known truth, pyIrena, Irena.  That is the form
the comparison should take in a paper or a referee response, because it
separates *agreement between implementations* from *agreement with reality*.

## Reading the tolerances

Some quantities are not individually determined by scattering data, and the
tables say so rather than hiding it:

* **Correlated prefactors.** In the Unified model `B` carries units
  cm⁻¹ Å⁻ᴾ, so a 1 % shift in `P` moves `B` by several percent. `B` is scored
  loosely whenever `P` is free. In Debye-Bueche only the product
  `Prefactor · Eta²` is determined, so that product is what is scored.
* **Polynomial background coefficients.** A cubic background's four
  coefficients are strongly correlated; the *background curve* is what the data
  determine, so `VALIDATION_RESULTS.md` scores the background evaluated at five
  Q values and reports the coefficients for information only.
* **Regularised inversions.** The Size Distribution tool solves an
  ill-conditioned inverse problem and returns a smoothed version of the
  generating distribution. Its integral moments — volume fraction, mean and RMS
  radius — are recovered to a fraction of a percent; the modal radius, which
  depends on the shape, has a wider tolerance.
* **Grid truncation.** Modelled size distributions are integrated over a grid
  that omits 1 % of the cumulative distribution at each tail, so a recovered
  volume fraction can sit a percent or two below the generated value. The
  generator integrates over ±6σ.
* **Nuisance backgrounds.** A small flat background that contributes a few
  percent of the intensity anywhere is weakly determined; it is scored loosely.

Every row also carries a **model-versus-exact-curve** comparison — the median
and maximum point-by-point deviation of the fitted model from `<name>_ideal.dat`
— which is tolerance-free, tool-independent, and directly reproducible in Irena.

## The regression tests

`pyirena/tests/test_validation_data.py` runs on every `pytest` invocation and
takes well under a second. It asserts that

* pyIrena's forward models agree with the independent implementations in
  `_models.py` to machine precision (this is what certifies the generator),
* every generated file exists, loads, and round-trips its ideal curve,
* the slit-smeared dataset declares its `dQl`, and
* a fast subset of the datasets recovers its generating parameters.

The full table across all tools is deliberately *not* a unit test: several of
those fits take tens of seconds. Run `run_validation_report.py` for that.

## What this exercise has already caught

Building the hard-sphere Modeling dataset exposed a real defect: structure-factor
parameters were packed into the fit vector under the same key group the surface
fractal population used for its attributes, so every optimiser step wrote them
onto the population object instead of into `sf_params`. They were "fitted" in
the sense of appearing in the parameter vector, but the model never saw them
change, and `fit()` returned them at their starting values with no warning. See
the entry in `CHANGELOG.md` and `pyirena/tests/test_modeling_structure_factor.py`.
