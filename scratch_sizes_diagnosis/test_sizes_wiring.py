"""
Ground-truth wiring test for pyIrena SizesDistribution inversion methods.

We build a KNOWN wide volume distribution of spheres on the solver's own
radius grid, forward-model I(q) = G @ x_true exactly (no model mismatch),
add small noise, then invert with each method and compare the recovered
VOLUME distribution to ground truth.

If a method is correctly wired, result['distribution'] should match Dv_true.
"""
import numpy as np
from pyirena.core.sizes import SizesDistribution
from pyirena.core.form_factors import build_g_matrix, make_r_grid, bin_widths

np.random.seed(0)

# ── Solver configuration (shared by all methods) ──────────────────────────
r_min, r_max, n_bins, log_spacing = 20.0, 6000.0, 120, True
shape, contrast = 'sphere', 1.0

r_grid = make_r_grid(r_min, r_max, n_bins, log_spacing)
dw = bin_widths(r_grid)

# ── Ground-truth VOLUME distribution Dv(r): wide log-normal in radius ──────
# center radius ~500 A (diameter ~1000 A), spanning ~100-5000 A diameter
r0, sigma = 500.0, 0.5           # geometric mean radius, log-width
Dv_true = (1.0 / (r_grid * sigma * np.sqrt(2 * np.pi))
           * np.exp(-(np.log(r_grid / r0)) ** 2 / (2 * sigma ** 2)))
Dv_true *= 0.01 / np.trapezoid(Dv_true, r_grid)   # total vol fraction = 0.01

# per-bin volume fraction that the solver's x_raw represents
x_true = Dv_true * dw

# ── Forward model on the SAME grid: I = G @ x_true ────────────────────────
q = np.logspace(np.log10(3e-4), np.log10(0.3), 250)
G = build_g_matrix(q, r_grid, shape, contrast)
I_clean = G @ x_true
err = 0.01 * I_clean
I_noisy = I_clean + err * np.random.randn(len(q))

Vf_true = float(np.trapezoid(Dv_true, r_grid))
rg_true = float(np.sqrt(np.trapezoid(r_grid ** 2 * Dv_true, r_grid) / Vf_true))
mean_r_true = float(np.trapezoid(r_grid * Dv_true, r_grid) / Vf_true)
print(f"GROUND TRUTH: Vf={Vf_true:.4e}  Rg={rg_true:.1f} A  "
      f"vol-mean r={mean_r_true:.1f} A  peak r={r_grid[np.argmax(Dv_true)]:.1f} A")
print("=" * 78)


def run(method, **kw):
    s = SizesDistribution()
    s.r_min, s.r_max, s.n_bins, s.log_spacing = r_min, r_max, n_bins, log_spacing
    s.shape, s.contrast = shape, contrast
    s.method = method
    s.background = 0.0
    s.fractional_error = False
    for k, v in kw.items():
        setattr(s, k, v)
    res = s.fit(q, I_noisy, err)
    if not res.get('success'):
        print(f"{method:16s} FAILED: {res.get('message')}")
        return
    Dv = res['distribution']
    rg = res['rg']
    Vf = res['volume_fraction']
    peak_r = r_grid[np.argmax(Dv)]
    vmean = np.trapezoid(r_grid * Dv, r_grid) / Vf if Vf > 0 else 0
    # shape agreement with truth (cosine similarity of normalized Dv)
    a = Dv / (np.linalg.norm(Dv) + 1e-300)
    b = Dv_true / np.linalg.norm(Dv_true)
    cos = float(np.dot(a, b))
    # fraction of total volume sitting in the single largest bin
    frac_top = float(np.max(Dv * dw) / np.sum(Dv * dw)) if np.sum(Dv*dw) > 0 else 0
    print(f"{method:16s} chi2={res['chi_squared']/len(q):7.2f}  "
          f"Vf={Vf:.3e}  Rg={rg:7.1f}  vmean_r={vmean:7.1f}  "
          f"peak_r={peak_r:7.1f}  shape_cos={cos:.3f}  top_bin_frac={frac_top:.2f}")


print("method            chi2/N     Vf         Rg      vmean_r   peak_r    shape   topbin")
print("-" * 78)
run('maxent')
run('tnnls')
run('regularization')
run('montecarlo', montecarlo_n_repetitions=6, montecarlo_max_iter=40000)
print("=" * 78)
print("Ideal recovery: peak_r~500, vmean_r~%.0f, Rg~%.0f, shape_cos~1.0, "
      "topbin small." % (mean_r_true, rg_true))
