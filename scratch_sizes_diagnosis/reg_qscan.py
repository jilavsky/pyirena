"""Regularization sensitivity to the high-Q inversion limit, real USAXS data."""
import json, numpy as np, h5py
from pyirena.core.sizes import SizesDistribution

cfg = json.load(open('/sessions/bold-tender-mccarthy/mnt/uploads/pyirena_config.json'))['sizes']
f = h5py.File('/sessions/bold-tender-mccarthy/mnt/uploads/Al1Si__40C_3min_1516.h5', 'r')
g = f['/entry/Al1Si_-40C_3min/sasdata']
Q = g['Q'][:]; I = g['I'][:]; Idev = g['Idev'][:]
m = np.isfinite(Q) & np.isfinite(I) & (Q > 0)
Q, I, Idev = Q[m], I[m], Idev[m]
print(f"data: {len(Q)} pts, Q {Q.min():.2e}..{Q.max():.2e}")

qlo = cfg['cursor_q_min']
print(f"config inversion window: {qlo:.2e} .. {cfg['cursor_q_max']:.2e}\n")

def build():
    s = SizesDistribution()
    s.r_min, s.r_max, s.n_bins = cfg['r_min'], cfg['r_max'], cfg['n_bins']
    s.log_spacing, s.shape, s.contrast = cfg['log_spacing'], cfg['shape'], cfg['contrast']
    s.background = cfg['background']
    s.power_law_B, s.power_law_P = cfg['power_law_B'], cfg['power_law_P']
    s.method = 'regularization'
    s.regularization_evalue = cfg['regularization_evalue']
    s.regularization_min_ratio = cfg['regularization_min_ratio']
    s.fractional_error = cfg['fractional_error']
    s.fractional_error_value = cfg['fractional_error_value']
    s.error_scale = cfg['error_scale']
    return s

print(f"{'q_hi':>10} {'Npts':>5} {'chi2/N':>8} {'nit':>4} {'Vf':>10} "
      f"{'Rg':>8} {'vmean_r':>8} {'peak_r':>8} {'topbin':>7}")
for qhi in (0.003, 0.005, 0.008, cfg['cursor_q_max'], 0.02, 0.03, 0.05, 0.1):
    sel = (Q >= qlo) & (Q <= qhi)
    s = build()
    res = s.fit(Q[sel], I[sel], Idev[sel])
    if not res.get('success'):
        print(f"{qhi:10.4f} FAILED {res.get('message')}"); continue
    r = res['r_grid']; D = res['distribution']
    from pyirena.core.form_factors import bin_widths
    dw = bin_widths(r); Vf = res['volume_fraction']
    vmean = np.trapezoid(r*D, r)/Vf if Vf > 0 else 0
    topbin = float(np.max(D*dw)/np.sum(D*dw)) if np.sum(D*dw) > 0 else 0
    print(f"{qhi:10.4f} {sel.sum():5d} {res['chi_squared']/sel.sum():8.2f} "
          f"{res['n_iterations']:4d} {Vf:10.3e} {res['rg']:8.1f} {vmean:8.1f} "
          f"{r[np.argmax(D)]:8.1f} {topbin:7.2f}")
