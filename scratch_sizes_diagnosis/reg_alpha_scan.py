import json, numpy as np, h5py
from scipy.optimize import nnls
from pyirena.core.sizes import SizesDistribution
from pyirena.core.form_factors import build_g_matrix, make_r_grid, bin_widths

cfg = json.load(open('/sessions/bold-tender-mccarthy/mnt/uploads/pyirena_config.json'))['sizes']
f = h5py.File('/sessions/bold-tender-mccarthy/mnt/uploads/Al1Si__40C_3min_1516.h5','r')
g=f['/entry/Al1Si_-40C_3min/sasdata']; Q=g['Q'][:]; I=g['I'][:]; Idev=g['Idev'][:]
mm=np.isfinite(Q)&np.isfinite(I)&(Q>0); Q,I,Idev=Q[mm],I[mm],Idev[mm]
qlo=cfg['cursor_q_min']
r_grid=make_r_grid(cfg['r_min'],cfg['r_max'],cfg['n_bins'],cfg['log_spacing']); dw=bin_widths(r_grid)
s=SizesDistribution()
for k in ('shape','contrast','background','power_law_B','power_law_P'): setattr(s,k,cfg[k])

def make_L(N):
    L=np.zeros((N-2,N))
    for i in range(N-2): L[i,i]=1;L[i,i+1]=-2;L[i,i+2]=1
    return L

def prep(qhi):
    sel=(Q>=qlo)&(Q<=qhi); q=Q[sel]
    Isub=I[sel]-s.compute_complex_background(q)
    err=np.maximum(np.abs(I[sel])*cfg['fractional_error_value'],1e-300)
    return build_g_matrix(q,r_grid,cfg['shape'],cfg['contrast']),Isub,err

def solve(G,Iv,err,la):
    N=G.shape[1];W=1/err;L=make_L(N);sa=np.sqrt(10.0**la)
    A=np.vstack([G*W[:,None],sa*L]);b=np.concatenate([Iv*W,np.zeros(N-2)])
    x,_=nnls(A,b,maxiter=10*N);res=(Iv-G@x)/err
    return x,float(np.dot(res,res))

for qhi in (0.02,0.05):
    G,Iv,err=prep(qhi);N=len(Iv)
    c2min=solve(G,Iv,err,-5)[1]
    print(f"q_hi={qhi} N={N} chi2_min(la=-5)={c2min:.0f}")
    print(f"   {'log_a':>6} {'chi2/N':>9} {'chi2/c2min':>10} {'Vf':>10} {'vmean_r':>8} {'peak_r':>8} {'topbin':>6}")
    for la in np.arange(-5,13.001,1.0):
        x,c2=solve(G,Iv,err,la); D=x/np.maximum(dw,1e-300)
        Vf=np.trapezoid(D,r_grid); vmean=np.trapezoid(r_grid*D,r_grid)/Vf if Vf>0 else 0
        topbin=np.max(x)/np.sum(x) if np.sum(x)>0 else 0
        print(f"   {la:6.1f} {c2/N:9.1f} {c2/c2min:10.2f} {Vf:10.3e} {vmean:8.1f} {r_grid[np.argmax(D)]:8.1f} {topbin:6.2f}")
    print()
# reference: stable low-q result
Gr,Ir,er=prep(0.0108)
xr,c2r=solve(Gr,Ir,er,None) if False else (None,None)
s2=SizesDistribution()
for k in ('r_min','r_max','n_bins','log_spacing','shape','contrast','background','power_law_B','power_law_P','fractional_error','fractional_error_value','error_scale'): setattr(s2,k,cfg[k])
s2.method='regularization';s2.regularization_evalue=cfg['regularization_evalue'];s2.regularization_min_ratio=cfg['regularization_min_ratio']
sel=(Q>=qlo)&(Q<=0.0108)
res=s2.fit(Q[sel],I[sel],Idev[sel]); D=res['distribution']
print(f"REFERENCE (q_hi=0.0108, full method): Vf={res['volume_fraction']:.3e} vmean_r={np.trapezoid(r_grid*D,r_grid)/res['volume_fraction']:.1f} peak_r={r_grid[np.argmax(D)]:.1f}")
