"""e3test.py (W4): numerical test of formula (E3) against the '-3<m1>m2' channel of A3's monochromatic kappa (k0=1, eta_*=2, eta_o=8).
(E3): T_<1>2,l = 8 pi sigma^4 Delta_0 Delta_l conj(Dbar_l),  Dbar_l = (1/2) int dmu K2_{l0}(k zhat, q) / (i^l sqrt(4 pi (2l+1))),
so its contribution to Delta^kappa_l is -3 T/(4 pi Delta_l) = -6 Delta_0 conj(Dbar_l).  K2 = kernel of m2 (= AB coefficient / 2)."""
import sys, os, pickle, numpy as np
from numpy.polynomial import legendre as Lg
from scipy.special import eval_legendre
d = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'A3_third_order'); sys.path.insert(0, d); os.chdir(d)
from thirteen import G, RAY, Delta_l, project, sky_grid
from engine import Ctx, LT, aadd, ascale, amul
N = int(sys.argv[1]) if len(sys.argv) > 1 else 24
x, w = Lg.leggauss(N)
nv, wv = sky_grid(24, 48)
th = np.arccos(np.clip(nv[:, 2], -1, 1)); ph = np.arctan2(nv[:, 1], nv[:, 0])
kv = np.repeat([[0, 0, 1.0]], N, 0)
qv = np.stack([np.sqrt(1 - x**2), 0 * x, x], 1)
ctx = Ctx([kv, qv], nv, g=G, ray=RAY)
L1, _ = ctx.termlist(LT['L_split'][1], maxorder=1)
L2, _ = ctx.termlist(LT['L_split'][2], maxorder=2)
m2 = aadd(ascale(L2, 2.0), amul(L1, L1))[3]            # (N, Npix) AB coefficient
r = pickle.load(open('mono_C.pkl', 'rb'))
Dch = project(r['kappa']['-3<m1>m2'], r['nu'], r['wnu'], 6)
D0 = Delta_l(0, 1.0)
for l in range(1, 7):
    Y = np.sqrt((2 * l + 1) / (4 * np.pi)) * eval_legendre(l, nv[:, 2])
    Kl0 = 4 * np.pi * (m2 / 2) @ (wv * np.conj(Y))   # (N,)
    Dbar = 0.5 * np.sum(w * Kl0) / (1j**l * np.sqrt(4 * np.pi * (2 * l + 1)))
    print('l=%d  E3-pred %+.6f%+.1ei   channel %+.6f%+.1ei' % (l, (-6 * D0 * np.conj(Dbar)).real, (-6 * D0 * np.conj(Dbar)).imag, Dch[l].real, Dch[l].imag))
