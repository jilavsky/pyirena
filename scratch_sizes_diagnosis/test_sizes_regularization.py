"""
Probe the regularization + montecarlo behaviour under more realistic,
harder conditions: flat background, truncated Q range, linear vs log grid.
"""
import numpy as np
from pyirena.core.sizes import SizesDistribution
from pyirena.core.form_factors import build_g_matrix, make_r_grid, bin_widths

np.random.seed(1)


def make_data(qmin, qmax, bkg, r0=500.0, sigma=0.5, vf=0.01, n_q=250,
              r_min=20.0, r_max=6000.0, n_bins=120, log_spacing=True):
    r_grid = make_r_grid(r_min, r_max, n_bins, log_spacing)
    dw = bin_widths(r_grid)
    Dv = (1.0 / (r_grid * sigma * np.sqrt(2*np.pi))
          * np.exp(-(np.log(r_grid/r0))**2 / (2*sigma**2)))
    Dv *= vf / np.trapezoid(Dv, r_grid)
    x_true = Dv * dw
    q = np.logspace(np.log10(qmin), np.log10(qmax), n_q)
    G = build_g_matrix(q, r_grid, 'sphere', 1.0)
    I = G @ x_true + bkg
    err = 0.01 * I
    I += err * np.random.randn(len(q))
    return q, I, err, r_grid, Dv, dw


def run(tag, method, q, I, err, r_grid, Dv_true, dw, **kw):
    s = SizesDistribution()
    s.r_min, s.r_max = r_grid[0], r_grid[-1]
    s.n_bins = len(r_grid)
    s.log_spacing = not np.allclose(np.diff(r_grid), np.diff(r_grid)[0])
    s.shape, s.contrast = 'sphere', 1.0
    s.method = method
    s.background = kw.pop('background', 0.0)
    s.fractional_error = False
    for k, v in kw.items():
        setattr(s, k, v)
    res = s.fit(q, I, err)
    if not res.get('success'):
        print(f"  {tag:28s} {method:14s} FAILED")
        return
    Dv = res['distribution']
    Vf = res['volume_fraction']
    vmean = np.trapezoid(r_grid*Dv, r_grid)/Vf if Vf > 0 else 0
    b = Dv_true/np.linalg.norm(Dv_true)
    a = Dv/(np.linalg.norm(Dv)+1e-300)
    cos = float(np.dot(a, b))
    frac_top = float(np.max(Dv*dw)/np.sum(Dv*dw)) if np.sum(Dv*dw) > 0 else 0
    peak_r = r_grid[np.argmax(Dv)]
    print(f"  {tag:28s} {method:14s} chi2={res['chi_squared']/len(q):6.2f} "
          f"Vf={Vf:.2e} vmean_r={vmean:7.1f} peak_r={peak_r:7.1f} "
          f"cos={cos:.3f} topbin={frac_top:.2f} nit={res['n_iterations']}")


scenarios = [
    ("full-q log", dict(qmin=3e-4, qmax=0.3, bkg=0.0, log_spacing=True)),
    ("full-q linear", dict(qmin=3e-4, qmax=0.3, bkg=0.0, log_spacing=False)),
    ("hi-q-truncated log", dict(qmin=3e-4, qmax=0.02, bkg=0.0, log_spacing=True)),
    ("lo-q-truncated log", dict(qmin=2e-3, qmax=0.3, bkg=0.0, log_spacing=True)),
    ("with-flat-bkg log", dict(qmin=3e-4, qmax=0.3, bkg=0.05, log_spacing=True)),
]

for tag, cfg in scenarios:
    ls = cfg.pop('log_spacing')
    q, I, err, r_grid, Dv_true, dw = make_data(log_spacing=ls, **cfg)
    print(tag, f"(log_spacing={ls})")
    for m in ('maxent', 'regularization', 'montecarlo'):
        extra = dict(montecarlo_n_repetitions=4, montecarlo_max_iter=30000) if m == 'montecarlo' else {}
        run(tag, m, q, I, err, r_grid, Dv_true, dw, **extra)
    print()
