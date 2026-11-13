# Developer Guide: Adding a New Form Factor

This guide explains exactly how to add a new particle-shape form factor to the
pyIrena Modeling tool.  It is written for both human developers and AI coding
agents — follow the steps in order and the form factor will work in the GUI,
batch scripting, HDF5 output, and the Collect Values tool automatically.

---

## Overview of the architecture

The form-factor system has three layers:

| Layer | File | Responsibility |
|-------|------|---------------|
| Physics | `pyirena/core/form_factors.py` | Compute G-matrix column `V·F²(Q,r)` |
| GUI | `pyirena/gui/modeling_panel.py` | Expose the new shape in the dropdown + parameter rows |
| Docs | `docs/modeling_gui.md` | Keep the user-facing Form Factors table current |

HDF5 I/O (`pyirena/io/nxcansas_modeling.py`) and state management
(`pyirena/state/state_manager.py`) are **generic** — they store `ff_params` as
an arbitrary dict, so no changes are required there.

---

## Step 1 — Physics: `pyirena/core/form_factors.py`

### 1a. Add an amplitude helper (if needed)

If your form factor has a new analytical amplitude function (i.e. it is not
just the sphere amplitude `_sphere_amplitude(qr)` reused with a different
geometry), add it as a module-level private function:

```python
def _cylinder_amplitude(qr_perp: np.ndarray, qr_par: np.ndarray) -> np.ndarray:
    """
    Normalised cylinder amplitude F(Q, r, L, θ).
    ...
    """
    ...
```

Small-argument Taylor branches are strongly recommended for any function that
involves `sin(x)/x` or `(sin x − x·cos x)/x³` to avoid catastrophic
cancellation at Q≈0.  See `_sphere_amplitude` for the pattern.

### 1b. Add a scalar-r form factor function (optional but useful for testing)

This function computes `V(r) × F²(Q,r)` for a single radius `r` (a scalar),
returning a 1-D array over Q.  It is used for documentation and unit tests
but is **not** called by the engine (the G-matrix builder is).

```python
def your_ff(q: np.ndarray, r: float, param_a: float = 1.0) -> np.ndarray:
    """
    Your form factor per unit volume fraction.

    Returns V(r) × F_norm²(Q, r)   [Å³]

    Args:
        q:       1-D array of Q values [Å⁻¹]
        r:       Characteristic radius [Å]
        param_a: Shape-specific parameter

    Returns:
        1-D array of length len(q)  [Å³]
    """
    q = np.asarray(q, dtype=float)
    V = ...                           # particle volume in Å³
    F = ...                           # normalised amplitude F(Q,r,param_a)
    return V * F ** 2
```

### 1c. Add a G-matrix builder

This is the function the engine actually calls.  It must return a 2-D array
of shape `(len(q), len(r_grid))` in units of cm⁻¹ per unit volume fraction.

**Simple form factor** (no orientation averaging — follow `_build_g_sphere`):

```python
def _build_g_your_ff(
    q: np.ndarray,
    r_grid: np.ndarray,
    contrast: float,
    param_a: float = 1.0,
) -> np.ndarray:
    """Vectorized G matrix for Your form factor.

    G[i,j] = V(r[j]) * F²(Q_i, r[j], param_a) * contrast * 1e-4
    """
    # Broadcast: shape (M, N)
    qr = q[:, np.newaxis] * r_grid[np.newaxis, :]
    V  = ...  # volumes, shape (N,)
    F  = _your_ff_amplitude(qr, param_a)   # shape (M, N)
    return V[np.newaxis, :] * F ** 2 * (contrast * 1e-4)
```

**Orientation-averaged form factor** (follow `_build_g_spheroid` — loop over r
bins to keep each slice L2-cache-friendly):

```python
def _build_g_your_ff(
    q: np.ndarray,
    r_grid: np.ndarray,
    contrast: float,
    param_a: float = 1.0,
) -> np.ndarray:
    """G matrix for Your form factor via per-r-bin Gauss-Legendre integration."""
    # Precompute quadrature nodes/weights (cache at module level for performance)
    cos_t   = 0.5 * (_GL_NODES + 1.0)   # orientation nodes, shape (K,)
    weights = 0.5 * _GL_WEIGHTS          # quadrature weights, shape (K,)

    M, N = len(q), len(r_grid)
    G = np.empty((M, N), dtype=float)

    for j, r in enumerate(r_grid):
        V = ...   # volume of particle with semi-axis r
        # orientation-dependent effective radius or parameter
        r_eff = r * np.sqrt(1.0 + (param_a ** 2 - 1.0) * cos_t ** 2)
        qr_eff = q[:, np.newaxis] * r_eff[np.newaxis, :]   # (M, K)
        F_eff  = _sphere_amplitude(qr_eff)                  # (M, K)
        G[:, j] = V * (F_eff ** 2 @ weights)

    return G * (contrast * 1e-4)
```

> **Note:** `_GL_NODES, _GL_WEIGHTS = leggauss(50)` is already defined at
> module level in `form_factors.py`.  Import `leggauss` from
> `numpy.polynomial.legendre` and add the cached pair if it does not yet exist.

### 1d. Register the builder

Add your builder to the `_G_BUILDERS` dict near the end of the module,
before `build_g_matrix`:

```python
_G_BUILDERS: dict[str, callable] = {
    'sphere':   _build_g_sphere,
    'spheroid': _build_g_spheroid,
    'your_ff':  _build_g_your_ff,   # ← add this line
}
```

---

## Step 2 — GUI: `pyirena/gui/modeling_panel.py`

### 2a. Add entry to `FF_LABELS`

```python
FF_LABELS = {
    'sphere':   ('Sphere',   []),
    'spheroid': ('Spheroid', ['aspect_ratio']),
    'your_ff':  ('Your FF Name', ['param_a', 'param_b']),  # ← add
}
```

The list `['param_a', 'param_b']` must exactly match the `**shape_params`
keyword argument names in `_build_g_your_ff`.  Empty list `[]` means the FF
has no extra parameters beyond the universal radius.

### 2b. Add display labels to `FF_PARAM_LABELS`

```python
FF_PARAM_LABELS = {
    'aspect_ratio': 'Aspect ratio',
    'param_a': 'Your parameter A description',  # ← add
    'param_b': 'Your parameter B description',  # ← add
}
```

The GUI reads `FF_PARAM_LABELS` to build the parameter rows automatically.
No further GUI changes are needed.

---

## Step 3 — State defaults (only if non-trivial defaults are needed)

`pyirena/state/state_manager.py` initialises `ff_params` as an empty dict
by default, which is fine for parameters that are set to `0` or `1` when
absent.  If your FF requires a non-zero / non-one default **and** you want
that default to appear when the user first selects the shape, add it to
`DEFAULT_STATE`:

```python
# In DEFAULT_STATE, inside each population template:
'ff_params': {
    'aspect_ratio': 1.0,   # existing
    'param_a': 2.5,        # ← add if needed
},
```

This step is **optional** — omit it if the engine safely handles a missing key
(e.g. via `ff_params.get('param_a', default_value)`).

---

## Step 4 — Documentation: `docs/modeling_gui.md`

Add a row to the **Form Factors** table in Section 6:

```markdown
| Your FF Name | `your_ff` | param_a [units], param_b [units] | One-sentence description |
```

---

## Verification

Run this quick check after all edits to confirm the G-matrix builds correctly:

```python
python -c "
from pyirena.core.form_factors import build_g_matrix
import numpy as np

q      = np.linspace(0.01, 0.5, 50)
r_grid = np.linspace(10, 500, 100)
G = build_g_matrix(q, r_grid, shape='your_ff', contrast=1.0, param_a=2.5)
assert G.shape == (50, 100), f'Wrong shape: {G.shape}'
assert np.all(np.isfinite(G)), 'G contains non-finite values'
assert np.all(G >= 0),        'G contains negative values'
print('G matrix OK:', G.shape, '  min/max:', G.min(), G.max())
"
```

Also verify the GUI loads without error:

```bash
python -c "import pyirena.gui.modeling_panel; print('modeling_panel import OK')"
```

---

## Checklist

- [ ] `_build_g_your_ff` function added to `form_factors.py`
- [ ] Builder registered in `_G_BUILDERS`
- [ ] `FF_LABELS['your_ff']` entry added to `modeling_panel.py`
- [ ] `FF_PARAM_LABELS` updated for all new parameter keys
- [ ] State defaults updated (if non-trivial defaults needed)
- [ ] Row added to Form Factors table in `docs/modeling_gui.md`
- [ ] Verification script passes
