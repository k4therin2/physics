"""check_planewave.py -- CHECK (iv): single plane wave zeta = cos(k (x + x_o)) (observer at the origin), k = 1, along +x.
Th1(n) = sum_s c_s K_1(s k; n),  Th2(n) = sum_{s1 s2} c_s1 c_s2 K_2(s1 k, s2 k; n),  c_+ = e^{i k x_o}/2, c_- = conj.
Compared with the nonperturbative planar numerics of N1 (a_l0, polar axis +x, real Y_l0) for 4 configurations,
using both C2 forms (shear and direct) of L2.  The (k,-k) pairs have K = 0 exactly (zero_mode prescription in Pair)."""
import sys, os, json
import numpy as np
from numpy.polynomial import legendre as Lg
from scipy.special import eval_legendre
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kernel2 import Geometry, MapKernels
N1 = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'N1_planar', 'output')
mu, wmu = Lg.leggauss(80)
nv = np.stack([mu, np.sqrt(1 - mu**2), 0 * mu], axis=1)
kx = np.array([1.0, 0, 0])
alm = lambda f, l: np.sqrt(4 * np.pi * (2 * l + 1)) * 0.5 * np.sum(wmu * f * eval_legendre(l, mu))
worst = {}
for eo, xo, tag in ((8.0, 0.0, 'to8_xo0'), (4.0, 0.0, 'to4_xo0'), (8.0, np.pi / 2, 'to8_xo0.5pi'), (4.0, np.pi / 2, 'to4_xo0.5pi')):
    g = Geometry(2.0, eo)
    ref = json.load(open(os.path.join(N1, 'coeff_%s.json' % tag)))
    c = {1: np.exp(1j * xo) / 2, -1: np.exp(-1j * xo) / 2}
    for form in ('shear', 'direct'):
        Th1 = 0; Th2 = 0
        for s1 in (1, -1):
            for s2 in (1, -1):
                mk = MapKernels(g, s1 * kx, s2 * kx, nv, form=form)
                K2, ch = mk.K2()
                Th2 = Th2 + c[s1] * c[s2] * K2
                if s2 == -s1:
                    Th1 = Th1 + c[s1] * mk.K1(mk.sa)
        print('== %s  form=%s   max|Im Th1|=%.1e max|Im Th2|=%.1e' % (tag, form, np.max(np.abs(Th1.imag)), np.max(np.abs(Th2.imag))))
        for l in range(1, 13):
            a1, a2 = alm(Th1.real, l), alm(Th2.real, l)
            r1, r2 = ref['Th1 = c1']['a_l0'][l], ref['Th2 = 2 c2']['a_l0'][l]
            print('  l=%2d  Th1 %+.10f (N1 %+.10f)   Th2 %+.10f (N1 %+.10f)  diff %.1e' % (l, a1, r1, a2, r2, abs(a2 - r2)))
            worst[(tag, form)] = max(worst.get((tag, form), 0), abs(a2 - r2), abs(a1 - r1))
    print('  N1 Th2 uncertainty (max over map):', ref['Th2 = 2 c2']['spread'])
print('WORST |a_l0 - N1| (l=1..12):', worst)
