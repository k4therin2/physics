"""check_ic2.py -- V3: claim (4). Independent derivation of Phi2(k, eta->0) from Gaussian zeta:
 zeta_nl = (1/4) ln(1+delta) + (1/2) ln(1-2Phi) expanded by sympy; 00 constraint at eta->0 from MY E^{(2)};
 Q2(0) from MY traceless source; all with symbolic (k, q, mu) and Taylor-expanded transfer functions."""
import sys, sympy as sp
from pw_tools import *
C1 = '/home/user/physics/cmb4/code/C1_einstein'; sys.path.insert(0, C1)
import kernels_api as KA; K = KA.load()
from fourier import eta, k, q, p, mu
E = E_mine()
eps = sp.Symbol('eps')
# 1. zeta expansion
d1, d2, P1, P2 = sp.symbols('delta1 delta2 Phi1 Phi2')
dl = eps * d1 + eps**2 * d2 / 2; Ph = eps * P1 + eps**2 * P2 / 2
zeta = sp.log(1 + dl) / 4 + sp.log(1 - 2 * Ph) / 2
z1 = sp.series(zeta, eps, 0, 3).removeO().coeff(eps, 1); z2 = 2 * sp.series(zeta, eps, 0, 3).removeO().coeff(eps, 2)
print("zeta1 =", sp.expand(z1), ";  zeta^(2) =", sp.expand(z2), "  [C1: delta2/4 - delta1^2/4 - Phi2 - 2 Phi1^2]")
assert sp.expand(z2 - (d2 / 4 - d1**2 / 4 - P2 - 2 * P1**2)) == 0
# 2. symbolic plane waves with Taylor-expanded transfer functions (exact small-eta limits)
kk, qq, m = sp.symbols('k q mu', positive=True)
NS = 10
def Tser(kmag):
    xv = kmag * t / sp.sqrt(3)
    return sp.series(-2 * (sp.sin(xv) - xv * sp.cos(xv)) / xv**3, t, 0, NS).removeO()
import pw_tools
pw_tools.T_Phi = Tser          # monkey-patch: mode() uses pw_tools.T_Phi/T_delta/T_V (module globals)
zq, zp = sp.symbols('zeta_q zeta_p')
kv, qv, pv = vec3(kk, qq, m)
pm = sp.sqrt(kk**2 + qq**2 - 2 * kk * qq * m)
f = {}
for nm, val in pw_tools.mode(zq, qv, qq).items(): f[nm] = f.get(nm, 0) + val
for nm, val in pw_tools.mode(zp, pv, pm).items(): f[nm] = f.get(nm, 0) + val
def src(a, b):
    return sp.expand(coeff_at_origin(substitute(E[('E', a, b, 2)], f), (zq, zp)))
S00 = src(0, 0); STL = src(3, 3) - (src(1, 1) + src(2, 2)) / 2
lim_eta2S00 = sp.simplify(sp.series(sp.expand(t**2 * S00), t, 0, 1).removeO())
Q0 = sp.simplify(sp.series(sp.expand(STL / kk**2), t, 0, 1).removeO())
print("lim eta^2 S00 =", lim_eta2S00, "   [C1: 32/3]")
print("Q2(0) =", sp.factor(Q0))
# zeta^(2) quadratic kernel: discrete coefficient of zeta_q zeta_p in [eps^2](-delta1^2/4 - 2 Phi1^2)... careful:
# [eps^2] of the quadratic terms of zeta_nl = (1/2) * (-delta1^2/4 - 2Phi1^2) -> continuous kernel = discrete coeff of [eps^2]
quad = sp.expand(coeff_at_origin(sp.expand((-f['delta1']**2 / 4 - 2 * f['Phi1']**2) / 2), (zq, zp)))
quad0 = sp.simplify(sp.series(quad, t, 0, 1).removeO())
print("zeta2 quadratic kernel at eta=0:", quad0, "  [expected from Phi1=-2/3, delta1=4/3 per zeta: 2*(-(4/3)^2/4 - 2*(4/9))/2 = ", sp.Rational(-4, 9) - sp.Rational(8, 9), "]")
# 3. delta2(0) from 00: delta2 = (eta^2/3) S00 - (2/3) eta^2 k^2 Phi2 - 2 eta Phi2' - 2 Psi2 ; regular Phi2 -> C at eta->0
C = sp.Symbol('C')
delta2_0 = lim_eta2S00 / 3 - 2 * (C - Q0)
zeta2_0 = delta2_0 / 4 - C + quad0
Csol = sp.factor(sp.solve(zeta2_0, C)[0])
print("Phi2(k, eta->0) kernel (mine):", Csol)
C1K = K['second']['K_Phi2_0'].subs(p, pm).subs({k: kk, q: qq})
print("difference from C1 K_Phi2_0:", sp.simplify(Csol - C1K))
C1Q = K['second']['K_Q0'].subs(p, pm).subs({k: kk, q: qq})
print("difference Q0 from C1 K_Q0:", sp.simplify(Q0 - C1Q))
# 4. PUB (3.12b,c) in our variables (PUB Psi = our Phi; R_PUB^(1) = -zeta; R^(2) = -2 R^(1)^2):
#  Phi2(0) = (2/3)(-2 zeta^2) + Phi1^2 - Delta^-1[dPhi dPhi] + 3 Delta^-2 dd[dPhi dPhi],  zeta^2 -> (9/4) Phi1^2
#  Fourier kernels per zeta_q zeta_p with Phi1 = -2/3 zeta:  Phi1^2 -> 4/9;  Delta^-1[dPhi dPhi] -> (4/9)(-(q.p))/(-k^2) = (4/9) q.p/k^2 ... sign: d_i Phi(q) d_i Phi(p) -> (i q)(i p) = -q.p ; Delta^-1 -> -1/k^2
qdp = qq * pm * ((kk**2 - qq**2 - pm**2) / (2 * qq * pm))        # q.p = (k^2 - q^2 - p^2)/2
kdq = kk * qq * m; kdp = kk**2 - kk * qq * m
inv1 = sp.Rational(4, 9) * (-qdp) / (-kk**2)                    # Delta^-1[d Phi d Phi]
inv2 = sp.Rational(4, 9) * (I * kk) ** 2 * 0                       # placeholder
# Delta^-2 d_J d_I [d_J Phi d_I Phi] -> (1/k^4) (i k_J)(i k_I) (i q_J)(i p_I) (4/9) = (4/9) (k.q)(k.p)/k^4
inv2 = sp.Rational(4, 9) * kdq * kdp / kk**4
PUB_Phi2 = sp.Rational(2, 3) * (-2) * sp.Rational(9, 4) * sp.Rational(4, 9) + sp.Rational(4, 9) - inv1 + 3 * inv2
PUB_Psi2mPhi2 = 4 * sp.Rational(4, 9) + 3 * inv1 - 9 * inv2
print("PUB (3.12b) Phi2(0) - mine:", sp.simplify(PUB_Phi2 - Csol))
print("PUB (3.12c) Psi2-Phi2 (= -Q2(0)) - mine:", sp.simplify(PUB_Psi2mPhi2 + Q0))
# 5. consistency of the k->0 (constant zeta) limit: angle-averaged Phi2(0) should be -8/9 (C1), Q(0) -> -16/9
print("angle average of Phi2(0) at k->0:", sp.limit(sp.integrate(Csol, (m, -1, 1)) / 2, kk, 0), "; of Q2(0):", sp.limit(sp.integrate(Q0, (m, -1, 1)) / 2, kk, 0))
