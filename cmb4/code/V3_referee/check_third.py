"""check_third.py -- V3: third-order contracted kernels K3c (C1 third_order.pkl) vs the independent E^{(3)} realised on
three discrete modes zeta1 e^{iq.x}, zeta2 e^{-iq.x}, zeta3 e^{ik.x} with explicit second-order fields for the pairs.
Identity tested:  6 * [coefficient of zeta1 zeta2 zeta3 in [eps^3]E_ab]  ==  K3c(k,q,mu;p=|k-q|) + K3c(k,q,-mu;p=|k+q|)
(sum over the two images of q; the tadpole pair (q,-q) is represented by abstract tadPhi2/tadQ2/taddelta2(eta)).
Abstract KPhi2(a,b,c,eta), Kh(a,b,c,eta), tadpoles and their eta-derivatives get random numerical values (same on both sides).
Verified second-order kernels (K_Q, K_S00, K_S0z, K_om, K_vT) are taken from C1 (see check_second.py / check_vector.py)."""
import random, sympy as sp
from pw_tools import *
from third_tools import *
E = E_mine(); K3 = K['third']
z1, z2, z3 = sp.symbols('zeta1 zeta2 zeta3')
worst = 0
for seed in range(3):
    rng = random.Random(100 + seed)
    pt = rand_point(seed); km = pt['k']
    f = three_mode_fields(pt, z1, z2, z3, amp=sp.Rational(1, 2))
    mine = {}
    for ab in [(0, 0), (0, 3), (0, 1), (1, 1), (2, 2), (3, 3), (1, 3)]:
        mine[ab] = 6 * coeff_at_origin(substitute(E[('E',) + ab + (3,)], f), (z1, z2, z3))
    print("point", seed, pt)
    checks = {'S00': (mine[(0, 0)], c1_image_sum(K3['S00'], pt, 1)), 'S0z': (mine[(0, 3)], c1_image_sum(K3['S0z'], pt, 1)),
              'Q (=S_TL/k^2)': ((mine[(3, 3)] - (mine[(1, 1)] + mine[(2, 2)]) / 2) / km**2, c1_image_sum(K3['Q'], pt, 1)),
              'Sxx-Syy': (mine[(1, 1)] - mine[(2, 2)], c1_image_sum(K3['Sxxyy'], pt, 1)),
              'S0x (image s->-s)': (mine[(0, 1)], c1_image_sum(K3['S0x'], pt, -1)), 'Sxz (image s->-s)': (mine[(1, 3)], c1_image_sum(K3['Sxz'], pt, -1))}
    for nm, (a, b) in checks.items():
        va, vb = ev(a, rng, pt['eta']), ev(b, rng, pt['eta']); d = abs(va - vb); worst = max(worst, d)
        print("   %-20s mine=%s  C1=%s  |diff|=%.2e" % (nm, sp.N(va, 12), sp.N(vb, 12), d))
print("WORST third-order contracted-kernel discrepancy:", worst)
