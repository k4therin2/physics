"""test_sign.py (A3): C2's `definitions` dict writes k1(eta) = k1_init + int_eta^eta_o G1 (and same sign in k2/2).
Flip the sign of the G-integral in the engine and show L2(shear) != L2(direct): the correct form is k1 = k1_init - int G1."""
import numpy as np
import engine
from engine import Ctx, LT, aadd, amap, ascale
orig = Ctx.functional
def flipped(self, kind, order, comp):
    key = ('fun', kind, order, comp)
    if kind == 'k' and order == 1 and key not in self.cache:
        r = self.ray
        k1i = self.evaluate(engine.RC['k_init'][comp].c[1].subs(engine.NZ), 'obs')
        G1 = self.evaluate(engine.RC['G1'][comp].subs(engine.NZ), 'ray')
        self.cache[key] = aadd(amap(k1i, lambda v: v[..., None] * np.ones(r.M)), amap(G1, r.cum), +1.0)
    return orig(self, kind, order, comp)
Ctx.functional = flipped
rng = np.random.default_rng(3)
nv = rng.normal(size=(6, 3)); nv /= np.linalg.norm(nv, axis=1)[:, None]
ctx = Ctx([[[0.7, -0.2, 0.4]], [[-0.3, 0.9, 0.5]], [[0.2, 0.1, -1.1]]], nv)
L2s, _ = ctx.termlist(LT['L_split'][2], maxorder=2); L2d, _ = ctx.termlist(LT['L_direct_split'][2], maxorder=2)
print('C2-definitions sign (k1 = k1_init + int G1): max|shear - direct| =', max(np.max(np.abs(L2s[m] - L2d[m])) for m in (3, 5, 6)))
