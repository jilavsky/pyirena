<!--
Hand-maintained commentary appended verbatim to VALIDATION_RESULTS.md by
run_validation_report.py.  The numbers quoted below refer to the Irena values in
irena_values.csv; revisit them if those change.
-->

## Notes on the pyIrena-Irena comparison

### Quantities Irena defines differently

**Modeling diffraction-peak width.**  pyIrena reports the Gaussian standard
deviation sigma; Irena reports the half-width at half-maximum.  This is the only
row in the table where the two packages differ by more than 8 %, and it is
entirely a difference in the reported quantity.  On
`modeling_sizedist_plus_peak` the ratio of the two reported widths is
0.01413 / 0.012 = **1.17750**, against sqrt(2 ln 2) = **1.17741** for a Gaussian
— agreement to five digits.  Converted to sigma the Irena width is 0.012001,
i.e. **+0.008 %** from the true value rather than the +17.8 % the raw numbers
suggest.

**WAXS polynomial background coefficients** (`bg0`, `bg1`, … marked `*`).  The
two packages parameterise the peak-fit background differently, so the
coefficients cannot be compared term by term.  The quantity the data actually
determine is the background *curve*, which is tabulated at five Q values in each
WAXS section and agrees with the truth to better than 2 % in pyIrena.

### Getting a like-for-like comparison

Three settings have to match before the parameters mean the same thing in both
packages.  Each was found the hard way while this table was assembled, and each
produces a large apparent disagreement while leaving the fit quality untouched
— so a good reduced chi-squared on both sides is not on its own evidence that
the two analyses are comparable.

1. **The contrast.** Each dataset's contrast is stated in the `#` header of its
   `.dat` file and in the settings line of its section here; it is not the same
   for every dataset (`modeling_sizedist_plus_peak` uses 80, its neighbours
   100).  Since the fitted scale absorbs the contrast, using the wrong one
   shifts the recovered volume fraction and nothing else.
2. **Cut-off radii in the Unified model.** `unified_two_level` links `RgCO` on
   level 2 to `Rg` of level 1.  Fitting the same data without that cut-off gives
   a model that describes the curve just as well (0.11 % median deviation from
   the exact curve, against 0.14 % with it) but with visibly different `G_1`,
   `B_1`, `B_2` and `P` values.
3. **Width conventions**, as above.

### The remaining rows that deviate by more than a percent or two

Every one of these is a quantity the data determine only weakly.  They are worth
listing because in each case the *pattern* of the disagreement is informative.

**Power-law prefactors.**  `B` carries units cm^-1 A^-P, so it is strongly
correlated with `P`: a 1.6 % shift in `P_1` moves `B_1` by more than 20 %.  In
`unified_two_level` both packages are off the true `B_1` by the same +23.4 %,
and off the true `B_2` by the same +6.0 % — they agree with **each other** to
0.03 % on both.  Two independent implementations landing on the same value of an
ill-determined parameter is a stronger result than either agreeing with the
truth would be.  The `pop2 B` of `modeling_sizedist_plus_unified` is the one
prefactor where the packages differ from each other (+5.3 % against -0.8 %), for
the same reason: their fitted `P` differ by 0.3 %.

**Modal radii from the Size Distribution tool.**  `sizes_sphere_bimodal` mode 1
(+2.3 % against -4.7 %) and `sizes_spheroid_lognormal` (+5.7 % against +2.3 %)
are the two largest.  A regularised inversion returns a smoothed distribution, so
where its peak sits depends on the regularisation, which the two packages do not
tune identically.  The integral moments, which do not depend on the smoothing,
agree to a few tenths of a percent throughout.

**The log-normal shape in `modeling_sizedist_plus_unified`.**  Both packages put
`sdeviation` high (+5.4 % and +7.3 %) and Irena also puts the median high
(+4.5 %).  This population's Guinier region sits on the Porod tail of the Unified
level, so its width and median trade against that level's `B` and `P`; the two
packages are biased in the same direction, which is what a shared degeneracy
looks like rather than an implementation difference.

**Reduced chi-squared** is a property of a fit, not a parameter of the model, and
the two packages weight and count degrees of freedom slightly differently.  The
values are listed for information; a difference of a few percent there carries no
meaning.

### One value read against the table

**`simple_porod`, `Kp`** was entered as `2.5137`, three decades from the true
`2.5e-06`, and has been read as `2.5137e-06` (deviation +0.548 %).  Two entries
in `unified_two_level` (`B_1`, `B_2`) arrived as `2.3206-7` and `1.4765-6`, with
the exponent's `e` lost in transcription, and have been rewritten as `2.3206e-7`
and `1.4765e-6`.

### Provenance of the numbers

The pyIrena column was produced by `validationData/run_validation_report.py`
against pyIrena **1.1.0b9**.  One row depends on the version: the hard-sphere
structure-factor parameters in the Modeling tool were not actually varied during
fitting before 1.1.0b9 (they were packed into the fit vector under a key group
the surface-fractal population already used, so each optimiser step wrote them to
the wrong place).  On `modeling_hardsphere` that defect gave reduced chi-squared
397 with the structure-factor and distribution parameters wrong by 8-40 %; the
row in this table is the fixed behaviour.

### Reproducing this table

```bash
python3 validationData/generate_validation_data.py    # rebuild the data
python3 validationData/run_validation_report.py       # fit it all, write the tables
```

The pyIrena columns come from that second command.  The Irena columns are read
from `validationData/irena_values.csv`, which holds the values obtained with the
Igor package, keyed by dataset and quantity; a quantity with no entry there is
left blank.  Adding more Irena results means adding rows to that CSV and
re-running the report — the numbers are data in the repository, not text inside
a document, so regenerating the report never loses them.

`validationData/fill_irena_deviations.py` fills the `Irena dev %` column of a
hand-edited copy of this table, and with `--export-csv` writes that copy's Irena
column back into `irena_values.csv`.

This prose section is maintained by hand in `validationData/irena_notes.md` and
is appended to the report verbatim.
