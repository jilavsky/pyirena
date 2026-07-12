import json, numpy as np, h5py
from pyirena.core.sizes import SizesDistribution
from pyirena.core.form_factors import bin_widths

cfg = json.load(open('/sessions/bold-tender-mccarthy/mnt/uploads/pyirena_config.json'))['sizes']
f = h5py.File('/sessions/bold-tender-mccarthy/mnt/uploads/Al1Si__40C_3min_1516.h5','r')
g=f['/entry/Al1Si_-40C_3min/sasdata']; Q=g['Q'][:]; I=g['I'][:]; E=g['Idev'][:]
m=np.isfinite(Q)&np.isfinite(I)&(Q>0); Q,I,E=Q[m],I[m],E[m]
sel=(Q>=cfg['cursor_q_min'])&(Q<=cfg['cursor_q_max'])

s=SizesDistribution()
for k in ('r_min','r_max','n_bins','log_spacing','shape','contrast','background',
          'power_law_B','power_law_P','fractional_error','fractional_error_value','error_scale'):
    setattr(s,k,cfg[k])
s.method='montecarlo'; s.montecarlo_n_repetitions=cfg['montecarlo_n_repetitions']
s.montecarlo_max_iter=cfg['montecarlo_max_iter']
np.random.seed(0)
res=s.fit(Q[sel], I[sel], E[sel])

r=res['r_grid']; D=res['distribution']; dw=bin_widths(r)
xraw = D*dw                      # per-bin volume fraction
Vf=res['volume_fraction']
print(f"chi2/N={res['chi_squared']/sel.sum():.2f}  Vf={Vf:.3e}  Rg(reported)={res['rg']:.1f}")
print(f"peak of density D(r) at r={r[np.argmax(D)]:.1f}")
print(f"peak of per-bin volfrac x_raw at r={r[np.argmax(xraw)]:.1f}")

# where does the VOLUME actually live?
cum = np.cumsum(xraw)/np.sum(xraw)
for frac in (0.5, 0.9):
    idx=np.searchsorted(cum, frac)
    print(f"  {int(frac*100)}% of total volume is below r={r[idx]:.0f} A")
# where does the SCATTERED INTENSITY (invariant-ish, low Q) come from?
# contribution to I at lowest fitted Q from each bin:
from pyirena.core.form_factors import build_g_matrix
G=build_g_matrix(Q[sel], r, cfg['shape'], cfg['contrast'])
contrib_lowQ = G[0,:]*xraw            # intensity at Q_min from each bin
c2=np.cumsum(contrib_lowQ)/np.sum(contrib_lowQ)
idx=np.searchsorted(c2,0.5)
print(f"  50% of the LOW-Q intensity comes from bins below r={r[idx]:.0f} A")

# consistency check: Rg recomputed two ways
Rg_from_D = np.sqrt(np.trapezoid(r**2*D, r)/np.trapezoid(D, r))
print(f"Rg recomputed from displayed D(r) = {Rg_from_D:.1f}")
# fraction of volume above 100 A
above100 = np.sum(xraw[r>100])/np.sum(xraw)
print(f"fraction of volume ABOVE 100 A = {above100:.3f}")
print(f"fraction of volume ABOVE 1000 A = {np.sum(xraw[r>1000])/np.sum(xraw):.3f}")
