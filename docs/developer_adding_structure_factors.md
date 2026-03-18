# Developer Guide: Adding a New Structure Factor

This guide explains exactly how to add a new inter-particle structure factor to
the pyIrena Modeling tool.  Follow the steps in order — HDF5 I/O and state
management are handled generically, so only three files ever need to change.

---

## Overview of the architecture

| Layer | File | Responsibility |
|-------|------|---------------|
| Physics | `pyirena/core/modeling.py` | Compute S(Q) and apply it to I_raw |
| GUI | `pyirena/gui/modeling_panel.py` | Expose the new SF in the dropdown + parameter rows |
| Docs | `docs/modeling_gui.md` | Keep the user-facing Structure Factors table current |

`pyirena/io/nxcansas_modeling.py` and `pyirena/state/state_manager.py` store
`sf_params` as arbitrary dicts and require **no changes**.

---

## Step 1 — Physics: `pyirena/core/modeling.py`

### 1a. Add a private S(Q) function

Add the new function near the other structure-factor functions
(`_interferences_sf`, `_hard_sphere_sf`):

```python
def _your_sf(q: np.ndarray, param_a: float = 1.0, param_b: float = 0.1) -> np.ndarray:
    """
    Your structure factor S(Q).

    Returns an array of the same shape as q, with values >= 0.
    S(Q) = 1 means no correlation (dilute limit).

    Args:
        q:       1-D array of Q values [Å⁻¹]
        param_a: Description [units]
        param_b: Description [units]

    Returns:
        S(Q) array, same shape as q
    """
    # Protect against Q=0 if needed
    q_safe = np.maximum(q, 1e-10)
    S = ...   # your formula
    return np.where(q < 1e-10, 1.0, S)
```

**Requirements for S(Q):**
- Same shape as input `q`
- All values ≥ 0 (values > 0 strongly recommended to avoid division issues)
- S(Q→0) should approach a finite value or 1 (not diverge)
- At high Q, S(Q) → 1 (no long-range correlations)

### 1b. Add a branch in `apply_structure_factor()`

In the `apply_structure_factor` function (currently around line 342),
add a branch **before** the final `raise ValueError`:

```python
def apply_structure_factor(
    I_raw: np.ndarray,
    q: np.ndarray,
    pop: SizeDistPopulation,
) -> np.ndarray:
    """Multiply I_raw by S(Q) for the population's structure factor."""
    sf = pop.structure_factor.lower()
    if sf == 'none':
        return I_raw
    if sf == 'interferences':
        eta  = float(pop.sf_params.get('eta',  100.0))
        pack = float(pop.sf_params.get('pack',   0.1))
        return I_raw * _interferences_sf(q, eta, pack)
    if sf == 'hard_sphere':
        rad = float(pop.sf_params.get('radius', 50.0))
        vf  = float(pop.sf_params.get('volume_fraction', 0.1))
        return I_raw * _hard_sphere_sf(q, rad, vf)
    if sf == 'your_sf':                                          # ← add this block
        param_a = float(pop.sf_params.get('param_a', 1.0))
        param_b = float(pop.sf_params.get('param_b', 0.1))
        return I_raw * _your_sf(q, param_a=param_a, param_b=param_b)
    raise ValueError(f"Unknown structure factor: {pop.structure_factor!r}")
```

Use `sf_params.get('key', default)` with sensible physical defaults so the SF
degrades gracefully when parameters have not been explicitly set.

### 1c. Register fittable parameters in `_pack_params()`

In `ModelingEngine._pack_params()`, the structure-factor active-parameter dict
is defined as:

```python
sf_active_keys = {
    'interferences': ['eta', 'pack'],
    'hard_sphere':   ['radius', 'volume_fraction'],
}.get(sf, [])
```

Add your SF's fittable parameter names to this dict:

```python
sf_active_keys = {
    'interferences': ['eta', 'pack'],
    'hard_sphere':   ['radius', 'volume_fraction'],
    'your_sf':       ['param_a', 'param_b'],   # ← add
}.get(sf, [])
```

The rest of the fitting machinery (bounds, "Fit?" checkboxes, MC uncertainty)
picks up the new parameters automatically.

---

## Step 2 — GUI: `pyirena/gui/modeling_panel.py`

### 2a. Add entry to `SF_LABELS`

```python
SF_LABELS = {
    'none':          ('None (dilute)', []),
    'interferences': ('Interferences', ['eta', 'pack']),
    'hard_sphere':   ('Hard Sphere (P-Y)', ['radius', 'volume_fraction']),
    'your_sf':       ('Your SF Display Name', ['param_a', 'param_b']),  # ← add
}
```

The list must exactly match the parameter names used in `apply_structure_factor`
and `_pack_params`.

### 2b. Add display labels to `SF_PARAM_LABELS`

```python
SF_PARAM_LABELS = {
    'eta':             'Correlation dist. η [Å]',
    'pack':            'Packing factor',
    'radius':          'HS radius [Å]',
    'volume_fraction': 'Volume fraction',
    'param_a':         'Your Parameter A [units]',   # ← add
    'param_b':         'Your Parameter B [units]',   # ← add
}
```

The GUI reads `SF_PARAM_LABELS` to build the parameter rows automatically.
No further GUI changes are needed.

---

## Step 3 — State defaults (only if non-trivial defaults are needed)

`pyirena/state/state_manager.py` initialises `sf_params` as an empty dict by
default.  If your SF requires a non-zero/non-one default **and** you want it to
appear when the user first selects the shape, add it to `DEFAULT_STATE`:

```python
# In DEFAULT_STATE, inside each population template:
'sf_params': {
    'eta': 100.0,           # existing
    'pack': 0.1,            # existing
    'radius': 50.0,         # existing
    'volume_fraction': 0.1, # existing
    'param_a': 1.0,         # ← add if needed
    'param_b': 0.1,         # ← add if needed
},
```

This step is **optional** — omit it if the engine handles a missing key safely
via `sf_params.get('param_a', default_value)` in Step 1b.

---

## Step 4 — Documentation: `docs/modeling_gui.md`

Add a row to the **Structure Factors** table in Section 7:

```markdown
| Your SF Name | `your_sf` | param_a [units], param_b [units] | One-sentence description |
```

---

## Verification

Run this quick check after all edits:

```python
python -c "
import numpy as np
from pyirena.core.modeling import (
    SizeDistPopulation, apply_structure_factor, ModelingEngine, ModelingConfig
)

# 1. S(Q) application
pop = SizeDistPopulation(structure_factor='your_sf',
                         sf_params={'param_a': 1.0, 'param_b': 0.1})
q   = np.linspace(0.01, 0.5, 50)
I_sf = apply_structure_factor(np.ones_like(q), q, pop)
assert I_sf.shape == q.shape,          'Wrong shape'
assert np.all(np.isfinite(I_sf)),      'Non-finite S(Q)'
assert np.all(I_sf >= 0),              'Negative S(Q)'
print('S(Q) application OK, mean S =', I_sf.mean())

# 2. Fitting round-trip (parameter packing)
config = ModelingConfig(
    populations=[pop],
    background=0.0, fit_background=False, q_min=0.01, q_max=0.5,
)
pop.sf_params_fit = {'param_a': True, 'param_b': False}
eng = ModelingEngine()
x0, lo, hi, keys = eng._pack_params(config)
sf_keys = [k for k in keys if len(k) >= 4 and k[2] == 'sf']
assert len(sf_keys) == 1, f'Expected 1 fittable SF param, got: {sf_keys}'
print('_pack_params OK:', sf_keys)
"
```

Also verify the GUI loads without error:

```bash
python -c "import pyirena.gui.modeling_panel; print('modeling_panel import OK')"
```

---

## Checklist

- [ ] `_your_sf(q, ...)` function added to `modeling.py`
- [ ] Branch added in `apply_structure_factor()` for `'your_sf'`
- [ ] `'your_sf': ['param_a', ...]` added to `sf_active_keys` dict in `_pack_params()`
- [ ] `SF_LABELS['your_sf']` entry added to `modeling_panel.py`
- [ ] `SF_PARAM_LABELS` updated for all new parameter keys
- [ ] State defaults updated (if non-trivial defaults needed)
- [ ] Row added to Structure Factors table in `docs/modeling_gui.md`
- [ ] Verification script passes
