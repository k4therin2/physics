"""test_third.py (A3): third-order consistency.  For k = zhat, |q| = 1, a few mu: phi-averaged coef_abc from the SHEAR-form L3
vs the DIRECT-form L3 (they differ by conservation-equation identities; the third-order fields enter through C1's contracted
solution, so agreement tests C1's third-order fields + C2's L3 + the engine's functionals/observer/emission)."""
import numpy as np, time
from thirteen import *
nu, wnu, nv = nu_grid(8)
phis = 2 * np.pi * (np.arange(16) + 0.5) / 16
t0 = time.time()
for mu in (0.31, -0.62):
    Xp = X3c(G, 1.0, 1.0, np.array([mu, -mu]))
    th = {nm: Xp[nm][0] + Xp[nm][1] for nm in ('Phi', 'dPhi', 'Psi', 'dPsi', 'V', 'dV')}
    th['dT3lin'] = (Xp['delta'][0] + Xp['delta'][1]) / 4
    res = {}
    for form in ('L_split', 'L_direct_split'):
        ch, ctx = coef_abc(1.0, 1.0, mu, phis, nv, form=form, third=th)
        tot = sum(ch.values()).mean(axis=0)
        res[form] = (tot, ch)
        print(form, 'mu=%g' % mu, 'phi-avg coef_abc(nu):', np.round(tot[:4], 6), '(%.0fs)' % (time.time() - t0))
    d = res['L_split'][0] - res['L_direct_split'][0]
    print('   max|shear - direct| =', np.max(np.abs(d)), '  scale', np.max(np.abs(res['L_split'][0])))
