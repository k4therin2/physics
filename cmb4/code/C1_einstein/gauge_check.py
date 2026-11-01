"""gauge_check.py -- control check (c): constant zeta is pure gauge.
Exact solution: eta = A etat, x^i = B xt^i (constant rescalings) of the RD background, re-expressed with a(etat)=etat:
  ghat_00 = -A^4, ghat_ij = A^2 B^2 delta_ij, delta = A^-4 - 1, v=0  =>  Psi=(A^4-1)/2, Phi=(1-A^2B^2)/2.
  zeta_nl = ln(1+delta_T) + (1/2) ln(1-2Phi) = ln(1/A) + ln(AB) = ln B.
Verify with the engine that this family solves E_{mu nu}=0 to O(eps^3); then analyse the k=0 second-order equations."""
import sympy as sp
from tpoly import *
d = load_polys('E.pkl'); E = d['E']
eps = sp.Symbol('eps'); eta = sp.Symbol('eta', positive=True)
a1, a2, a3, b1, b2, b3, c = sp.symbols('a1 a2 a3 b1 b2 b3 c')
A = 1 + eps * a1 + eps**2 * a2 / 2 + eps**3 * a3 / 6
B = 1 + eps * b1 + eps**2 * b2 / 2 + eps**3 * b3 / 6
exact = {'Psi': (A**4 - 1) / 2, 'Phi': (1 - A**2 * B**2) / 2, 'delta': A**(-4) - 1}
def coeff(expr, n):   # X_n = n! * [eps^n] X
    return sp.factorial(n) * sp.series(expr, eps, 0, 4).removeO().coeff(eps, n)
vals = {(f, n): sp.expand(coeff(ex, n)) for f, ex in exact.items() for n in (1, 2, 3)}
for kk, vv in sorted(vals.items()):
    print(kk, '=', vv)
def subs_const(f, od, m):
    if any(m):
        return sp.Integer(0)       # constants: all derivatives vanish
    if f in exact:
        return vals[(f, od)]
    return sp.Integer(0)           # om, h, v vanish
for a in range(4):
    for b in range(a, 4):
        ex = sp.expand(E[a][b].to_sympy(eta=eta, funcs=subs_const))
        assert ex == 0, (a, b, ex)
print("VERIFIED: the exact dilation/time-rescaling family (A,B) solves all E_{mu nu} to O(eps^3) for arbitrary a_n, b_n.")
# zeta_nl = ln B exactly:
zeta_series = sp.series(sp.log(B), eps, 0, 4).removeO()
print("zeta_nl = ln B =", sp.expand(zeta_series))
# --- k=0 second-order analysis: first-order regular mode has Psi1=Phi1 => 2a1 = -(a1+b1) => b1 = -3a1; zeta1 = b1 = c
first = {a1: -c / 3, b1: c}
print("first order: Phi1 = Psi1 =", sp.expand(vals[('Phi', 1)].subs(first)), ", delta1 =", vals[('delta', 1)].subs(first))
# generic constant second-order fields: substitute Phi2, Psi2, delta2 symbols, first order fixed
P2, S2, D2 = sp.symbols('Phi2 Psi2 delta2')
def subs2(f, od, m):
    if any(m): return sp.Integer(0)
    if od == 1: return vals[(f, 1)].subs(first) if f in exact else sp.Integer(0)
    if od == 2: return {'Phi': P2, 'Psi': S2, 'delta': D2}.get(f, sp.Integer(0))
    return sp.Integer(0)
eqs = []
for a in range(4):
    for b in range(a, 4):
        ex = sp.expand(2 * E[a][b].order_part(2).to_sympy(eta=eta, funcs=subs2))
        if ex != 0:
            print("E_%d%d^(2) at k=0:" % (a, b), ex)
            eqs.append(ex)
sol = sp.solve(eqs, [D2, S2], dict=True)
print("k=0 second-order equations solved for (delta2, Psi2) in terms of free Phi2:", sol)
# the exact family at second order with b1=-3a1: 
fam = {kk: sp.expand(vals[(f, 2)].subs(first)) for kk, f in (('Phi2', 'Phi'), ('Psi2', 'Psi'), ('delta2', 'delta'))}
print("exact family at 2nd order (free a2, b2):", fam)
# check family satisfies the k=0 equations
for e in eqs:
    assert sp.expand(e.subs({P2: fam['Phi2'], S2: fam['Psi2'], D2: fam['delta2']})) == 0
print("VERIFIED: exact family lies on the k=0 solution set.")
# zeta^(2) = 0 condition: zeta_nl = ln B => zeta2 = 2! [eps^2] ln B = b2 - b1^2 = 0 => b2 = c^2
print("zeta_nl^(2) = 2*[eps^2] ln B =", sp.expand(2 * zeta_series.coeff(eps, 2)), " => b2 = c^2 for zeta_nl = eps c exactly")
# remaining freedom a2: residual gauge at k=0.  The member with Psi2 = Phi2 (k->0 limit if the traceless source vanished):
fam_z = {kk: vv.subs(b2, c**2) for kk, vv in fam.items()}
print("with b2=c^2:", fam_z)
sol_a2 = sp.solve(fam_z['Psi2'] - fam_z['Phi2'], a2)
print("Psi2 = Phi2 selects a2 =", sol_a2)
final = {kk: sp.expand(vv.subs(a2, sol_a2[0])) for kk, vv in fam_z.items()}
print("==> constant-zeta second-order constants with Psi2=Phi2, zeta2=0:", final)
print("    delta_T2 := 2[eps^2](1+delta)^{1/4}:", sp.expand(2 * sp.series((1 + exact['delta'])**sp.Rational(1, 4), eps, 0, 3).removeO().coeff(eps, 2).subs(first).subs(b2, c**2).subs(a2, sol_a2[0])))
import pickle
pickle.dump({'family2': {kk: sp.srepr(vv) for kk, vv in fam.items()}, 'final': {kk: sp.srepr(vv) for kk, vv in final.items()}}, open('gauge_check.pkl', 'wb'))
