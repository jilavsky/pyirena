# Slit Smearing (USAXS)

pyIrena can fit **slit-smeared USAXS data** directly — the model is smeared to
match the data, rather than desmearing the data. This is the transparent,
SasView-style approach: fitted parameters are always **ideal-space**
(pinhole-equivalent), and the data are never modified.

## The physics

For an infinitely long slit of (half-)length `SL` (the APS/USAXS Lake
convention, `SL` = `dQl`):

```
I_smeared(q) = (1/SL) ∫₀^SL  I_ideal( √(q² + l²) )  dl
```

pyIrena evaluates the analytic model on an **extended q grid** and convolves it
with this integral (a fixed sparse operator `W`, so a fit loop costs one matvec
per iteration). For Size Distribution the whole G matrix is smeared once
(`G_smeared = W @ G(q_ext)`), so the inversion machinery is untouched.

Core engine: [`pyirena/core/smearing.py`](../pyirena/core/smearing.py)
(`SlitSmearer`, `smear_model`, `smear_curve`, `build_smearing_matrix`).
`SL ≤ 0` is a strict no-op, so tools call it unconditionally.

## How data is detected (NXcanSAS)

Slit-smeared data declare, on the `Q` dataset:

```
Q@resolutions = "dQw,dQl"
    dQw   per-point Q-width resolution  → becomes the loaded dQ
    dQl   scalar slit (half-)length     → the slit length that drives smearing
```

Presence of `dQl` means the loaded curve **is** slit smeared, and model
smearing is enabled automatically. A file with a plain `Q@resolutions="Qdev"`
is treated as pinhole.

Matilda writes **both** copies of a USAXS measurement:

```
entry/<sample>/sasdata          desmeared   (Q@resolutions="Qdev";  @default)
entry/<sample>_SMR/sasdata      slit smeared (Q@resolutions="dQw,dQl")
```

By default pyIrena loads the `@default` (desmeared) dataset. When a file has
both copies, the fitting panels show a **"Slit smeared data"** checkbox that
selects which one to load; checking it reloads the `_SMR` dataset and enables
smearing. The **Slit length** field is file-derived but editable (rarely
changed — the file value is the most trusted).

## Per-tool behaviour

| Tool | Behaviour with slit-smeared data |
|---|---|
| **Unified Fit** | Model (all levels + background) smeared in the fit loop; local Guinier/Porod cursor fits smeared too; "Show selected level" overlay smeared. Invariant is computed from the **ideal model** (parameter-based) and labelled as such. |
| **Size Distribution** | G matrix + power-law background smeared; recovered distribution is ideal-space. |
| **Simple Fits** | Each analytic model smeared. The **Invariant** integrates the data directly, which is invalid for smeared data — it is **disabled with a message** (load the desmeared data, or use Unified Fit's model-based invariant). |
| **Modeling** | Total model (all populations + structure factors + background) smeared before comparison. |
| **Fractals** | Calculated intensity smeared for overlay/comparison when the loaded data are slit smeared. |
| **Data Merge** | Merging a slit-smeared curve (typically low-Q USAXS) with a pinhole one yields a slit-smeared output: the merged file gets a `dQl` dataset and provenance records the input/merged slit lengths. Two different nonzero slit lengths warn (larger kept). Optimization is unchanged (slit length ≤ SAXS Qmin ⇒ negligible in the overlap). |
| **Data Manipulation** | Subtract/divide refuse to mix a smeared curve with a pinhole one (or different slit lengths); the guard lives in the core engine so batch inherits it. Outputs drop any stale `_SMR` twin and orphaned `dQl`. |

## Performance

Model-side smearing evaluates the model on an **extended, refined q grid**
(≈ 8× the data points) every fit iteration, so a slit-smeared fit is inherently
a few times slower per iteration than a pinhole fit. This is expected and small
in absolute terms.

What is *not* acceptable is a fit that iterates far more than it should. In
Unified Fit, `ETA`/`PACK` only affect the model when a level's **correlations**
are enabled; they are therefore treated as free parameters only in that case.
Leaving them "fit" with correlations off would otherwise make the optimiser
thrash over parameters that cannot change the fit — cheap for pinhole data but
heavily amplified by the extended grid (this previously made some slit-smeared
fits ~50× slower). If a slit-smeared fit feels unexpectedly slow, check that
only parameters that actually affect the model are marked for fitting.

## Saved results

Every fitting tool's result group records:

```
<tool>_results/
    intensity_model         # SMEARED model when smearing was used (matches data)
    intensity_model_ideal   # pinhole model on the same Q (only when smeared)
    @slit_length            # scalar 1/Å; 0 or absent ⇒ pinhole fit
    @data_is_slit_smeared   # 1/0
```

Legacy readers keep working (`intensity_model` name unchanged); new
readers/plots offer a smeared/ideal toggle.

## Scripting (batch API) — the Matilda contract

Canonical JSON keys (used across the batch API, control API, and MCP):

- **`load_slit_smeared: true`** — in a tool's config block, loads the file's
  slit-smeared (`_SMR`) dataset. Smearing is then enabled automatically.
- **`slit_length: <float>`** — optional override of the file-derived slit
  length (1/Å).
- Reported back on load: **`is_slit_smeared`** (bool) and **`slit_length`**.

Enforcement: requesting `load_slit_smeared` on a file that has no `dQl` is a
**hard error** (the batch call returns `None` with a clear message) — never a
silent pinhole fit.

Example (Unified Fit batch config):

```json
{
  "_pyirena_config": {"version": "1.0"},
  "unified_fit": {
    "load_slit_smeared": true,
    "num_levels": 1,
    "levels": [{"level": 1,
                "G":  {"value": 300,  "fit": true},
                "Rg": {"value": 60,   "fit": true},
                "B":  {"value": 1e-3, "fit": true},
                "P":  {"value": 4.0,  "fit": true}}],
    "background": {"value": 0.01, "fit": true}
  }
}
```

## Control API / MCP (AI advisor)

`open_dataset(file_path, use_slit_smeared=True)` loads the `_SMR` dataset and
returns `is_slit_smeared`, `slit_length`, and `has_slit_smeared_entry`. When
the session's data are slit smeared, a Unified Fit model created for it smears
automatically; `run_fit` reports `slit_smearing_applied` and `slit_length` so
the AI advisor can explain that the fitted parameters are ideal-space.
