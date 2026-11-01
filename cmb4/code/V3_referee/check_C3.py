"""check_C3.py -- V3: C1's small-eta model KPhi2_small uses Phi2 = C + (3/10) S_Phi(0) eta^2 - ...; the correct
coefficient from Phi'' + 4Phi'/eta = S0 (or int_0^eta (eta'/3)(1-eta'^3/eta^3) S0) is 1/10.  Re-derive C3 with a
symbolic coefficient alpha and test whether C3 (and the limits Q3c(0), lim eta^2 S00c) depend on it."""
import sys, pickle, sympy as sp
C1 = '/home/user/physics/cmb4/code/C1_einstein'; sys.path.insert(0, C1)
import kernels_api as KA; K = KA.load()
from fourier import eta, k, q, p, mu, TF
S2 = K['second']; K3 = K['third']; IC = K['third_ic']
alpha = sp.Symbol('alpha')
a_, b_, c_ = sp.symbols('a b c', positive=True)
def at(expr, A, B, C):
    return expr.subs({q: a_, p: b_, k: c_}, simultaneous=True).subs({a_: A, b_: B, c_: C}, simultaneous=True)
A, B, C = q, k, p
KPhi2_0 = S2['K_Phi2_0']; SPhi0 = IC['SPhi0_2nd']
# check Green's-function coefficient: int_0^eta (eta'/3)(1 - eta'^3/eta^3) S0 deta' = S0 eta^2/10
etap = sp.Symbol('etap', positive=True)
print("int_0^eta G_Phi(k->0) deta' =", sp.simplify(sp.integrate(etap / 3 * (1 - etap**3 / eta**3), (etap, 0, eta))), " => coefficient 1/10, C1 used 3/10")
KPhi2_small = at(KPhi2_0, A, B, C) + eta**2 * (alpha * at(SPhi0, A, B, C) - C**2 * at(KPhi2_0, A, B, C) / 30)
print("C1's KPhi2_small reproduced with alpha=3/10:", sp.simplify(KPhi2_small.subs(alpha, sp.Rational(3, 10)) - IC['KPhi2_small']) == 0)
Kh_small = sp.Rational(8, 9) * eta**2
TP = lambda kk: TF['Phi'](kk, eta)
Xs = lambda kk: sp.series(KA.T_Phi_closed(kk), eta, 0, 8).removeO()
tads = IC['tads']
KPhi2f = sp.Function('KPhi2'); Khf = sp.Function('Kh')
def small_eta(expr, order=6):
    expr = expr.subs(KPhi2f(A, B, C, eta), KPhi2_small).subs(Khf(A, B, C, eta), Kh_small)
    for nm, fn in tads.items():
        expr = expr.subs(sp.Function(nm)(eta), sp.series(fn, eta, 0, order).removeO())
    expr = expr.doit()
    for kk in (q, p, k):
        expr = expr.subs(sp.Derivative(TP(kk), eta), sp.diff(Xs(kk), eta)).subs(TP(kk), Xs(kk))
    return expr.doit()
psub = {p: sp.sqrt(k**2 + q**2 - 2 * k * q * mu)}
def lim0(expr, pref=1):
    e = sp.expand(pref * small_eta(expr)).subs(psub)
    return sp.simplify(sp.series(e, eta, 0, 1).removeO())
Q3_0 = lim0(K3['Q']); S00_0 = lim0(K3['S00'], eta**2)
print("Q3c(0) depends on alpha:", Q3_0.has(alpha), "; equals C1:", sp.simplify(Q3_0 - IC['Q3_0']) == 0)
print("lim eta^2 S00c depends on alpha:", S00_0.has(alpha), "; equals C1:", sp.simplify(S00_0 - IC['S00_0']) == 0)
kz = lim0(IC['zeta3_contracted'])
print("zeta3 contracted kernel limit depends on alpha:", kz.has(alpha))
C3s = sp.Symbol('C3'); Phi3 = sp.Function('Phi3')(eta); de3 = sp.Function('delta3')(eta)
z3_0 = IC['zeta3_lin'].subs(de3, S00_0 / 3 - 2 * (C3s - Q3_0)).subs(Phi3, C3s)
C3_sol = sp.factor(sp.solve(sp.expand(z3_0 + kz), C3s)[0])
print("C3 depends on alpha:", C3_sol.has(alpha), "; equals C1's C3:", sp.simplify(C3_sol - IC['C3']) == 0)
print("C3 =", C3_sol)
print("squeezed limit q->0:", sp.limit(C3_sol, q, 0))
# also: does S_Phi^(3)c(eta->0) (the finite part) depend on alpha?  (not used for C3 but reported by C1 as 'finite')
SPhi3_0 = lim0(K3['SPhi'])
print("S_Phi^(3)c(eta=0) depends on alpha:", SPhi3_0.has(alpha), "(this O(1) value feeds the O(eta^2) term of Phi3c only)")
