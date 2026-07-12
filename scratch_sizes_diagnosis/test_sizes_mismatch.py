"""
Model-mismatch scenario closer to real data:
 - ground truth generated on an independent FINE grid (continuous dist)
 - optional low-q power-law upturn (large-scale scattering not in sphere model)
 - solver inverts on its own coarser grid
Watch whether regularization dumps a spike into the lowest-r bin while
maxent stays smooth.
"""
import numpy as np
from pyirena.core.sizes import SizesDistribution
from pyirena.core.form_factors import build_g_matrix, make_r_grid, bin_widths

np.random.seed(2)


def truth_intensity(q, r0=500.0, sigma=0.5, vf=0.01):
    rf = np.logspace(np.log10(5), np.log10(1e4), 2000)   # fine independent grid
    dwf = bin_widths(rf)
    Dv = (1.0/(rf*sigma*np.sqrt(2*np.pi))*np.exp(-(np.log(rf/r0))**2/(2*sigma**2)))
    Dv *= vf/np.trapezoid(Dv, rf)
    Gf = build_g_matrix(q, rf, 'sphere', 1.0)
    return Gf @ (Dv*dwf), rf, Dv


def run(tag, method, q, I, err, **kw):
    s = SizesDistribution()
    s.r_min, s.r_max, s.n_bins, s.log_spacing = 20.0, 6000.0, 120, True
    s.shape, s.contrast = 'sphere', 1.0
    s.method = method
    s.fractional_error = False
    for k, v in kw.items():
        setattr(s, k, v)
    res = s.fit(q, I, err)
    if not res.get('success'):
        print(f"  {method:14s} FAILED: {res.get('message')}"); return
    Dv = res['distribution']; rg_grid = res['r_grid']; dw = bin_widths(rg_grid)
    Vf = res['volume_fraction']
    peak_r = rg_grid[np.argmax(Dv)]
    vmean = np.trapezoid(rg_grid*Dv, rg_grid)/Vf if Vf > 0 else 0
    frac_top = float(np.max(Dv*dw)/np.sum(Dv*dw)) if np.sum(Dv*dw) > 0 else 0
    frac_lowest = float((Dv*dw)[0]/np.sum(Dv*dw)) if np.sum(Dv*dw) > 0 else 0
    print(f"  {method:14s} chi2={res['chi_squared']/len(q):7.2f} Vf={Vf:.2e} "
          f"vmean_r={vmean:7.1f} peak_r={peak_r:7.1f} topbin={frac_top:.2f} "
          f"lowest_bin={frac_lowest:.2f}")


q = np.logspace(np.log10(3e-4), np.log10(0.3), 300)
I0, rf, Dv_true = truth_intensity(q)
Vf_true = np.trapezoid(Dv_true, rf)
vmean_true = np.trapezoid(rf*Dv_true, rf)/Vf_true
print(f"TRUTH vmean_r={vmean_true:.1f}  Vf={Vf_true:.2e}\n")

for plaw in (0.0, 1e-9, 1e-8):
    I = I0 + plaw*q**(-3.5)          # low-q power-law upturn (large structures)
    err = 0.01*I
    Id = I + err*np.random.randn(len(q))
    print(f"power-law amp={plaw:g}")
    for m in ('maxent', 'regularization'):
        run(m, m, q, Id, err)
    print()
