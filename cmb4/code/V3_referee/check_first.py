"""check_first.py -- V3: claim (1), linear order, from the independent E_{mu nu}."""
import sympy as sp
from pw_tools import *
E = E_mine()
kk = sp.Symbol('k', positive=True)
ph = sp.exp(I * kk * z)
names = ['Phi1', 'Psi1', 'delta1', 'v11', 'v21', 'v31']
amps = {nm: sp.Function('A_' + nm)(t) * ph for nm in names}
A = lambda nm: sp.Function('A_' + nm)(t)
L = {ab: sp.expand(substitute(E[('E',) + ab + (1,)], amps).subs({x: 0, y: 0, z: 0})) for ab in [(0, 0), (0, 1), (0, 2), (0, 3), (1, 1), (2, 2), (3, 3), (1, 2), (1, 3), (2, 3)]}
TL = sp.simplify(L[(3, 3)] - (L[(1, 1)] + L[(2, 2)]) / 2)
print("traceless scalar:", sp.factor(TL), " => Psi1 = Phi1")
print("0x, 0y:", L[(0, 1)], L[(0, 2)], " => v_x = v_y = 0")
vz = sp.solve(L[(0, 3)], A('v31'))[0]
V1 = sp.simplify(vz / (I * kk))
print("V1 =", V1, "  claim -(eta^2/2)(Phi'+Psi/eta):", sp.simplify(V1 + t**2 / 2 * (A('Phi1').diff(t) + A('Psi1') / t)) == 0)
de = sp.solve(L[(0, 0)], A('delta1'))[0]
print("delta1 =", de, "  claim:", sp.simplify(de - (-sp.Rational(2, 3) * t**2 * kk**2 * A('Phi1') - 2 * t * A('Phi1').diff(t) - 2 * A('Psi1'))) == 0)
tr = (L[(1, 1)] + L[(2, 2)] + L[(3, 3)]).subs(A('delta1'), de).subs(A('Psi1'), A('Phi1')).doit()
tr = sp.expand(tr / tr.coeff(A('Phi1').diff(t, 2)))
print("ODE:", tr, " claim:", sp.simplify(tr - (A('Phi1').diff(t, 2) + 4 / t * A('Phi1').diff(t) + kk**2 / 3 * A('Phi1'))) == 0)
TP = T_Phi(kk)
print("T_Phi solves ODE:", sp.simplify(tr.subs(A('Phi1'), TP).doit()) == 0, "; T_Phi(0) =", sp.limit(TP, t, 0))
# zeta normalisation: zeta_1 = lim (delta_T1 - Phi1) = lim (delta1/4 - Phi1)
print("lim (T_delta/4 - T_Phi) =", sp.limit(T_delta(kk) / 4 - TP, t, 0), "(must be 1)")
# full consistency: all ten linear equations vanish on the solution
sol = {A('Phi1'): TP, A('Psi1'): TP, A('delta1'): T_delta(kk), A('v31'): I * kk * T_V(kk), A('v11'): 0, A('v21'): 0}
print("all E^(1) vanish on the closed-form solution:", all(sp.simplify(v.subs(sol).doit()) == 0 for v in L.values()))
