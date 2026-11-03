"""typeB_check.py -- type-B (finite 3j sum) evaluation of the channel  Th2_B = m1^2 - 2<m1> m1  (kernel K1(k1)K1(k2) - Delta_0(k1)K1(k2)
- Delta_0(k2)K1(k1)), monochromatic, and comparison with the diagonal G_BB of the brute-force mu-integral (cl22.py).
C3 (B1) + monopole rules:   C^{BB}_l = C^{m1^2}_l - 4 <<m1>^2> C^(11)_l ,
   C^{m1^2}_l = 2 (4 pi sigma^2)^2 sum_{l1 l2} (2l1+1)(2l2+1)/(4 pi) (l l1 l2;000)^2 Delta_l1^2 Delta_l2^2,  <<m1>^2> = sigma^2 Delta_0^2,
   C^(11)_l = 4 pi sigma^2 Delta_l^2   (Delta_l incl. observer monopole (l=0) and observer dipole (l=1))."""
import sys, os, pickle
import numpy as np
from math import lgamma, exp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kernel2 import Geometry, Delta_l

def w3j000(a, b, c):
    L = a + b + c
    if L % 2 or c > a + b or c < abs(a - b): return 0.0
    gg = L // 2
    lg = 0.5 * (lgamma(L - 2 * a + 1) + lgamma(L - 2 * b + 1) + lgamma(L - 2 * c + 1) - lgamma(L + 2)) + lgamma(gg + 1) - lgamma(gg - a + 1) - lgamma(gg - b + 1) - lgamma(gg - c + 1)
    return (-1)**gg * exp(lg)

def CBB(k0, eta_s, eta_o, lmax, LT=50):
    g = Geometry(eta_s, eta_o)
    D = np.array([Delta_l(l, k0, g) for l in range(LT + 1)])
    out = []
    for l in range(lmax + 1):
        s = 0.0
        for l1 in range(LT + 1):
            for l2 in range(abs(l - l1), min(LT, l + l1) + 1):
                w = w3j000(l, l1, l2)
                if w: s += (2 * l1 + 1) * (2 * l2 + 1) / (4 * np.pi) * w * w * D[l1]**2 * D[l2]**2
        out.append(2 * (4 * np.pi)**2 * s - 4 * D[0]**2 * 4 * np.pi * D[l]**2)
    return np.array(out)

if __name__ == '__main__':
    tag = sys.argv[1]
    r = pickle.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cl22_%s.pkl' % tag), 'rb'))
    i = r['groups'].index('B: L1^2 - 2<m1>m1')
    cb = CBB(r['k0'], r['eta_s'], r['eta_o'], r['lmax'])
    print('tag', tag, ' type-B 3j sum vs mu-integral diagonal G_BB')
    for l in range(1, r['lmax'] + 1):
        print('%3d  3j: %.10e   mu-int: %.10e   rel %.1e' % (l, cb[l], r['G'][i, i, l], abs(cb[l] / r['G'][i, i, l] - 1)))
