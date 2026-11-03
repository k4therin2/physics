"""localisation: the direct = shear identity must hold for ARBITRARY observer-event / emission-time shift symbols
(C2 T2 treats them as free).  Perturb second-order symbols by random amounts and watch direct - shear."""
import numpy as np, sys
import engine
from thirteen import *
which = sys.argv[1:]  # e.g. deta_o2 dxo2 deta_e2
orig = engine.Ctx.symbol
rng = np.random.default_rng(7)
PERT = {}
def sym(self, name, order, comp):
    res = orig(self, name, order, comp)
    tag = '%s%d' % (name, order)
    if tag in which:
        key = (tag, comp)
        if key not in PERT: PERT[key] = {m: rng.normal() + 1j * rng.normal() for m in res}
        res = {m: v + PERT[key][m] for m, v in res.items()}
    return res
engine.Ctx.symbol = sym
nu, wnu, nv = nu_grid(6)
phis = 2 * np.pi * (np.arange(12) + 0.5) / 12
mu = 0.31
Xp = X3c(G, 1.0, 1.0, np.array([mu, -mu]))
th = {nm: Xp[nm][0] + Xp[nm][1] for nm in ('Phi', 'dPhi', 'Psi', 'dPsi', 'V', 'dV')}
th['dT3lin'] = (Xp['delta'][0] + Xp['delta'][1]) / 4
r = {}
for form in ('L_split', 'L_direct_split'):
    ch, ctx = coef_abc(1.0, 1.0, mu, phis, nv, form=form, third=th)
    r[form] = sum(ch.values()).mean(axis=0)
print(which, 'direct - shear:', np.round(r['L_direct_split'] - r['L_split'], 4))
