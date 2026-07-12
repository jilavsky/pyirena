import json, numpy as np, h5py
from pyirena.core.sizes import SizesDistribution
from pyirena.core.form_factors import bin_widths

cfg = json.load(open('/sessions/bold-tender-mccarthy/mnt/uploads/pyirena_config.json'))['sizes']
f = h5py.File('/sessions/bold-tender-mccarthy/mnt/uploads/Al1Si__40C_3min_1516.h5','r')
g=f['/entry/Al1Si_-40C_3min/sasdata']; Q=g['Q'][:]; I=g['I'][:]; E=g['Idev'][:]
m=np.isfinite(Q)&np.isfinite(I)&(Q>0); Q,I,E=Q[m],I[m],E[m]
sel=(Q>=cfg['cursor_q_min'])&(Q<=cfg['cursor_q_max'])

def fit(meth):
    s=SizesDistribution()
    for k in ('r_min','r_max','n_bins','log_spacing','shape','contrast','background',
              'power_law_B','power_law_P','fractional_error','fractional_error_value','error_scale'):
        setattr(s,k,cfg[k])
    s.method=meth
    if meth=='montecarlo': s.montecarlo_n_repetitions=10; s.montecarlo_max_iter=100000
    np.random.seed(1)
    return s.fit(Q[sel],I[sel],E[sel])

print(f"{'method':14s} {'Rg':>8} {'peak D(r)=dV/dr':>16} {'peak r*D=dV/dlnr':>17} {'%vol>100A':>9}")
for meth in ('maxent','tnnls','regularization','montecarlo'):
    r=fit(meth); rg=r['r_grid']; D=r['distribution']; dw=bin_widths(rg)
    x=D*dw
    vperlnr = rg*D                      # dV/dln(r)  (correct density for log-r axis)
    peakD = rg[np.argmax(D)]
    peakL = rg[np.argmax(vperlnr)]
    above100 = np.sum(x[rg>100])/np.sum(x)
    print(f"{meth:14s} {r['rg']:8.0f} {peakD:16.1f} {peakL:17.1f} {above100:9.2f}")
