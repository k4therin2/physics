"""e2test.py (W4): numerical test of formula (E2) against the '-3<m2>m1' channel of A3's monochromatic kappa (k0=1, eta_*=2, eta_o=8).
(E2): T_<2>1,l = <<m2>> C11_l + 4 pi sigma^4 Delta_l^2 int dmu P_l(mu) Kbar2(k0 khat, -k0 khat'),  mu = khat.khat'.
Its contribution to Delta^kappa_l is -3 T/(4 pi Delta_l) = -3 Delta_l [<<m2>> + int dmu P_l Kbar2]."""
import sys, os, pickle, numpy as np
from numpy.polynomial import legendre as Lg
from scipy.special import eval_legendre
d = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'A3_third_order'); sys.path.insert(0, d); os.chdir(d)
from thirteen import M2_sky, M2_tad, Delta_l, project
N = int(sys.argv[1]) if len(sys.argv) > 1 else 24
x, w = Lg.leggauss(N)
kv = np.repeat([[0, 0, 1.0]], N, 0)
qv = np.stack([np.sqrt(1 - x**2), 0 * x, x], 1)
AB = M2_sky(kv, -qv)            # AB coefficient of m2 for the pair (k, -q)
tad = M2_tad(1.0)
r = pickle.load(open('mono_C.pkl', 'rb'))
ch = r['kappa']['-3<m2>m1']
Dch = project(ch, r['nu'], r['wnu'], 6)
for l in range(1, 7):
    Dl = Delta_l(l, 1.0)
    I = np.sum(w * eval_legendre(l, x) * AB)
    for f, lab in ((0.5, 'Kbar=AB/2'), (1.0, 'Kbar=AB')):
        print('l=%d %-10s E2-pred %+.6f   channel %+.6f' % (l, lab, (-3 * Dl * f * (tad + I)).real, Dch[l].real))
