import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pyirena.core.sizes import SizesDistribution
from pyirena.core.form_factors import build_g_matrix, make_r_grid, bin_widths

np.random.seed(0)
r_min, r_max, n_bins = 20.0, 6000.0, 120
r_grid = make_r_grid(r_min, r_max, n_bins, True)
dw = bin_widths(r_grid)
Dv_true = (1/(r_grid*0.5*np.sqrt(2*np.pi))*np.exp(-(np.log(r_grid/500.0))**2/(2*0.25)))
Dv_true *= 0.01/np.trapezoid(Dv_true, r_grid)
q = np.logspace(np.log10(3e-4), np.log10(0.3), 250)
G = build_g_matrix(q, r_grid, 'sphere', 1.0)
I = G @ (Dv_true*dw); err = 0.01*I; I = I + err*np.random.randn(len(q))

def fit(method, **kw):
    s = SizesDistribution()
    s.r_min, s.r_max, s.n_bins, s.log_spacing = r_min, r_max, n_bins, True
    s.shape, s.contrast, s.method, s.fractional_error = 'sphere', 1.0, method, False
    for k, v in kw.items(): setattr(s, k, v)
    return s.fit(q, I, err)

res = {m: fit(m) for m in ('maxent', 'tnnls', 'regularization')}
mc = fit('montecarlo', montecarlo_n_repetitions=6, montecarlo_max_iter=40000)
# reconstruct the corrected MC distribution (x_raw/dw, no *V_r)
s = SizesDistribution(); s.r_min,s.r_max,s.n_bins,s.log_spacing=r_min,r_max,n_bins,True
s.shape,s.contrast='sphere',1.0; s.montecarlo_n_repetitions=6; s.montecarlo_max_iter=40000
x_raw,_,_ = s._fit_montecarlo(G, I, err, r_grid)
mc_fixed = x_raw/np.maximum(dw,1e-300)

fig, ax = plt.subplots(1, 2, figsize=(13, 5.2))
for a in ax:
    a.plot(r_grid, Dv_true, 'k-', lw=2.5, label='ground truth', zorder=5)
    a.set_xscale('log'); a.set_xlabel('radius  [Å]')
    a.set_ylabel('volume distribution  D$_v$(r)  [vol-frac / Å]')
    a.grid(alpha=0.3, which='both')

ax[0].plot(r_grid, res['maxent']['distribution'], label='MaxEnt')
ax[0].plot(r_grid, res['tnnls']['distribution'], '--', label='TNNLS')
ax[0].plot(r_grid, res['regularization']['distribution'], ':', lw=2, label='Regularization')
ax[0].plot(r_grid, mc['distribution'], color='red', lw=2, label='Monte Carlo (current)')
ax[0].set_title('As shipped: MaxEnt / TNNLS / Regularization agree;\nMonte Carlo shifted ~2× to large r')
ax[0].legend(fontsize=9)

ax[1].plot(r_grid, mc['distribution'], color='red', lw=2, label='Monte Carlo (current, ×V(r) bug)')
ax[1].plot(r_grid, mc_fixed, color='green', lw=2, label='Monte Carlo (fixed, x/Δr)')
ax[1].set_title('Monte Carlo: removing the extra ×V(r)=(4/3)πr³\nre-weighting restores the correct distribution')
ax[1].legend(fontsize=9)

plt.tight_layout()
plt.savefig('/sessions/bold-tender-mccarthy/mnt/outputs/sizes_method_diagnosis.png', dpi=130)
print('saved')
