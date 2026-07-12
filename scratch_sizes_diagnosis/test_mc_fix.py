"""
Confirm the Monte Carlo fix: x_raw from _fit_montecarlo is ALREADY volume
fraction per bin. Dividing by bin width (like every other method) should
recover the true volume distribution. The extra *V(r) in fit() is the bug.
"""
import numpy as np
from pyirena.core.sizes import SizesDistribution
from pyirena.core.form_factors import build_g_matrix, make_r_grid, bin_widths

np.random.seed(0)
r_min, r_max, n_bins = 20.0, 6000.0, 120
r_grid = make_r_grid(r_min, r_max, n_bins, True)
dw = bin_widths(r_grid)
r0, sigma = 500.0, 0.5
Dv_true = (1/(r_grid*sigma*np.sqrt(2*np.pi))*np.exp(-(np.log(r_grid/r0))**2/(2*sigma**2)))
Dv_true *= 0.01/np.trapezoid(Dv_true, r_grid)
x_true = Dv_true*dw
q = np.logspace(np.log10(3e-4), np.log10(0.3), 250)
G = build_g_matrix(q, r_grid, 'sphere', 1.0)
I = G @ x_true
err = 0.01*I
I = I + err*np.random.randn(len(q))

s = SizesDistribution()
s.r_min, s.r_max, s.n_bins, s.log_spacing = r_min, r_max, n_bins, True
s.shape, s.contrast = 'sphere', 1.0
s.montecarlo_n_repetitions = 6
s.montecarlo_max_iter = 40000

# call the raw MC fitter directly
x_raw, n_it, x_std = s._fit_montecarlo(G, I, err, r_grid)

Vf_true = np.trapezoid(Dv_true, r_grid)
vmean_true = np.trapezoid(r_grid*Dv_true, r_grid)/Vf_true

# CURRENT (buggy) pipeline: distribution = x_raw/dw, then *V(r), renormalized
dist = x_raw/np.maximum(dw, 1e-300)
V_r = (4/3)*np.pi*r_grid**3
dist_buggy = dist*V_r
dist_buggy *= np.trapezoid(dist, r_grid)/np.trapezoid(dist_buggy, r_grid)

# FIXED pipeline: distribution = x_raw/dw  (same as maxent/tnnls/regularization)
dist_fixed = dist

def stats(name, D):
    Vf = np.trapezoid(D, r_grid)
    vmean = np.trapezoid(r_grid*D, r_grid)/Vf
    cos = float(np.dot(D/np.linalg.norm(D), Dv_true/np.linalg.norm(Dv_true)))
    print(f"  {name:18s} Vf={Vf:.3e} vmean_r={vmean:7.1f} peak_r={r_grid[np.argmax(D)]:7.1f} cos={cos:.3f}")

print(f"TRUTH              Vf={Vf_true:.3e} vmean_r={vmean_true:7.1f} peak_r={r_grid[np.argmax(Dv_true)]:7.1f}")
stats("MC current (*V_r)", dist_buggy)
stats("MC fixed (x/dw)", dist_fixed)
