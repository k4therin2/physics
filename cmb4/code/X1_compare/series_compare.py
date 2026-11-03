"""X1: small-k0 series of our Delta_2 (C5 small_k_series, symbolic eta_*, eta_o) vs Adam's superhorizon C_2^(11)
   (his Sec. 9.2): C2^(11) = 4pi/2025 s^2 k^4 D^2 - 4pi/637875 s^2 k^6 (o-e)^2 (o+e)(59o^3+59o^2e-106oe^2+114e^3) + O(k^8), D = o^2-e^2,
   and his leading T_l coefficient (80): T_l = -k^l (o-e)^(l-1) [o + (l-1) e] / (3 (2l+1)!!) + O(k^(l+2))."""
import sys, os
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'C5_linear_Cl'))
import linear_cl as L
k, es, eo, sig2 = L.k, L.eta_s, L.eta_o, L.sigma2
o, e = eo, es
D2 = sp.expand(L.small_k_series(2, 7))
print('our Delta_2 small-k series:', sp.factor(D2))
C11 = sp.expand(4 * sp.pi * sig2 * D2**2)
C11 = sum(C11.coeff(k, n) * k**n for n in (4, 6))
adam = 4 * sp.pi / 2025 * sig2 * k**4 * (o**2 - e**2)**2 \
     - 4 * sp.pi / 637875 * sig2 * k**6 * (o - e)**2 * (o + e) * (59 * o**3 + 59 * o**2 * e - 106 * o * e**2 + 114 * e**3)
print('C11 (ours) - C11 (Adam) through k^6:', sp.simplify(C11 - sp.expand(adam)))
for l in (2, 3, 4, 5):
    Dl = sp.expand(L.small_k_series(l, l + 1))
    lead = Dl.coeff(k, l)
    adam_lead = -(o - e)**(l - 1) * (o + (l - 1) * e) / (3 * sp.factorial2(2 * l + 1))
    print('l=%d leading coefficient ours = %s ; ours - Adam(80) = %s' % (l, sp.factor(lead), sp.simplify(lead - adam_lead)))
