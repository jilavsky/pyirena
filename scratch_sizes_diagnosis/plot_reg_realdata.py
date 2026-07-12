import json, numpy as np, h5py
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.optimize import nnls
from pyirena.core.sizes import SizesDistribution
from pyirena.core.form_factors import build_g_matrix, make_r_grid, bin_widths

cfg = json.load(open('/sessions/bold-tender-mccarthy/mnt/uploads/pyirena_config.json'))['sizes']
f = h5py.File('/sessions/bold-tender-mccarthy/mnt/uploads/Al1Si__40C_3min_1516.h5','r')
g=f['/entry/Al1Si_-40C_3min/sasdata']; Q=g['Q'][:]; I=g['I'][:]; Idev=g['Idev'][:]
mm=np.isfinite(Q)&np.isfinite(I)&(Q>0); Q,I,Idev=Q[mm],I[mm],Idev[mm]
qlo=cfg['cursor_q_min']
r_grid=make_r_grid(cfg['r_min'],cfg['r_max'],cfg['n_bins'],cfg['log_spacing']); dw=bin_widths(r_grid)

def fit(qhi):
    s=SizesDistribution()
    for k in ('r_min','r_max','n_bins','log_spacing','shape','contrast','background',
              'power_law_B','power_law_P','fractional_error','fractional_error_value','error_scale'):
        setattr(s,k,cfg[k])
    s.method='regularization'; s.regularization_evalue=cfg['regularization_evalue']
    s.regularization_min_ratio=cfg['regularization_min_ratio']
    sel=(Q>=qlo)&(Q<=qhi)
    return s.fit(Q[sel],I[sel],Idev[sel])

def old_fallback(qhi):
    """Reproduce the previous behaviour: argmin-chi2 (smallest alpha)."""
    s=SizesDistribution()
    for k in ('shape','contrast','background','power_law_B','power_law_P'): setattr(s,k,cfg[k])
    sel=(Q>=qlo)&(Q<=qhi); q=Q[sel]
    Isub=I[sel]-s.compute_complex_background(q)
    err=np.maximum(np.abs(I[sel])*cfg['fractional_error_value'],1e-300)
    G=build_g_matrix(q,r_grid,cfg['shape'],cfg['contrast']); N=G.shape[1]
    L=np.zeros((N-2,N))
    for i in range(N-2): L[i,i]=1;L[i,i+1]=-2;L[i,i+2]=1
    W=1/err; sa=np.sqrt(10.0**-5.0)
    A=np.vstack([G*W[:,None],sa*L]); b=np.concatenate([Isub*W,np.zeros(N-2)])
    x,_=nnls(A,b,maxiter=10*N)
    x=np.maximum(x,cfg['regularization_min_ratio']*x.max())
    return x/np.maximum(dw,1e-300)

fig,ax=plt.subplots(1,2,figsize=(13.5,5.3))

# Left: data + background + windows
s0=SizesDistribution()
for k in ('shape','contrast','background','power_law_B','power_law_P'): setattr(s0,k,cfg[k])
ax[0].loglog(Q,I,'.',ms=3,color='0.5',label='data I(Q)')
ax[0].loglog(Q,s0.compute_complex_background(Q),'r-',lw=1,label='complex background (B·Q$^{-4}$+flat)')
for qv,c,lab in [(0.0108,'green','good window end (0.0108)'),(0.02,'orange','0.02'),(0.05,'red','0.05')]:
    ax[0].axvline(qv,color=c,ls='--',lw=1.2,label=lab)
ax[0].axvline(qlo,color='k',ls=':',lw=1)
ax[0].set_xlabel('Q [Å$^{-1}$]'); ax[0].set_ylabel('I [cm$^{-1}$]')
ax[0].set_title('Al1Si USAXS data — background dominates beyond ~0.01 Å$^{-1}$')
ax[0].legend(fontsize=8)

# Right: recovered distributions (log-log)
res=fit(0.0108); D=res['distribution']
floor=1e-9
ax[1].loglog(r_grid,np.maximum(D,floor),'k-',lw=2.5,label='q_hi=0.0108 (good window)')
ax[1].loglog(r_grid,np.maximum(old_fallback(0.02),floor),color='red',lw=1.4,ls='--',
             label='q_hi=0.02 OLD fallback (spike at r=10)')
for qhi,c in [(0.02,'tab:orange'),(0.03,'tab:blue'),(0.05,'tab:green')]:
    r=fit(qhi); ax[1].loglog(r_grid,np.maximum(r['distribution'],floor),color=c,lw=1.6,
                             label=f'q_hi={qhi} NEW fallback')
ax[1].set_xlim(10,1e4); ax[1].set_ylim(1e-7,1e-1)
ax[1].set_xlabel('radius [Å]'); ax[1].set_ylabel('D$_v$(r) [vol-frac/Å]')
ax[1].set_title('Regularization vs high-Q cutoff:\nOLD explodes into a spike; NEW degrades smoothly')
ax[1].legend(fontsize=8,loc='lower left')
plt.tight_layout()
plt.savefig('/sessions/bold-tender-mccarthy/mnt/outputs/reg_realdata_fallback.png',dpi=130)
print('saved')
