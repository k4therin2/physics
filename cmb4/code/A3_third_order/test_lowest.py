"""test_lowest.py (A3): validate the engine at orders 1 and 2.
(1) K1 from C2's L1 lists (shear and direct) vs the type-A expansion with Delta_l (standard SW+ISW+Doppler incl. observer terms).
(2) L2 shear form vs L2 direct form for a generic momentum pair (equal modulo conservation equations => numerically equal):
    tests k1, Dx1, deta_o1, dxo1, deta_e1, second-order fields, and the sign of the ray-functional integrals."""
import numpy as np
from scipy.special import spherical_jn, eval_legendre
from engine import Ctx, LT, aadd
import tf

rng = np.random.default_rng(3)
nv = rng.normal(size=(6, 3)); nv /= np.linalg.norm(nv, axis=1)[:, None]
k1 = np.array([[0.7, -0.2, 0.4]]); k2 = np.array([[-0.3, 0.9, 0.5]]); k3 = np.array([[0.2, 0.1, -1.1]])
ctx = Ctx([k1, k2, k3], nv)
r = ctx.ray
def Delta_l(l, k):
    chs = r.eta_o - r.eta_s
    v = (tf.tdt(k, r.eta_s) + tf.tphi(k, r.eta_s)) * spherical_jn(l, k * chs) - k * tf.tv(k, r.eta_s) * spherical_jn(l, k * chs, derivative=True) \
        + r.integ(2 * tf.tphi(k, r.t, 1) * spherical_jn(l, k * r.chi))
    if l == 0: v -= tf.tdt(k, r.eta_o) + tf.tphi(k, r.eta_o)
    if l == 1: v += k * tf.tv(k, r.eta_o) / 3
    return v
for form in ('L_split', 'L_direct_split'):
    K1, _ = ctx.termlist(LT[form][1])
    kk = np.linalg.norm(k1); nu = nv @ k1[0] / kk
    ref = sum((2 * l + 1) * 1j**l * Delta_l(l, kk) * eval_legendre(l, nu) for l in range(60))
    print(form, 'K1 engine vs type-A: max|diff| =', np.max(np.abs(K1[1][0] - ref)), ' |K1| ~', np.max(np.abs(ref)))
L2s, ps = ctx.termlist(LT['L_split'][2], maxorder=2)
L2d, pd = ctx.termlist(LT['L_direct_split'][2], maxorder=2)
for m in (3, 5, 6):
    print('pair mask', m, ' shear:', np.round(L2s[m][0][:3], 6), ' direct:', np.round(L2d[m][0][:3], 6), ' max|diff| =', np.max(np.abs(L2s[m] - L2d[m])))
