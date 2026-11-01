"""first_order.py -- SVT decomposition at first order, ODE, solutions, checks."""
import sympy as sp, pickle
from tpoly import *
from fourier import *

d = load_polys('E.pkl'); E = d['E']; cons = d['cons']
L = lambda a, b: linear_part(E[a][b], 1)
Phi1, Psi1, de1 = amp('Phi', 1), amp('Psi', 1), amp('delta', 1)
vx, vy, vz = [amp('v%d' % i, 1) for i in (1, 2, 3)]
E00, E0x, E0y, E0z = L(0, 0), L(0, 1), L(0, 2), L(0, 3)
Exx, Eyy, Ezz, Exy, Exz, Eyz = L(1, 1), L(2, 2), L(3, 3), L(1, 2), L(1, 3), L(2, 3)
print("E00 =", E00); print("E0z =", E0z); print("E0x =", E0x)
print("Ezz =", Ezz); print("Exx =", Exx); print("Exy =", Exy, " Exz =", Exz)
# traceless scalar part: Ezz - (Exx+Eyy)/2 ∝ (ik)^2 (Phi-Psi)
TL = sp.simplify(Ezz - (Exx + Eyy) / 2)
print("traceless scalar:", sp.factor(TL))
solPsi = sp.solve(TL, Psi1)[0]
print("=> Psi1 =", solPsi)
# vector part of 0i: E0x, E0y
print("0i transverse:", E0x, E0y, " => v_x = v_y = 0 (no transverse velocity at first order)")
solvz = sp.solve(E0z, vz)[0]
print("=> v_z =", solvz, "  i.e. v_i = d_i V with V = -(eta^2/2)(Phi' + Psi/eta)")
V1 = sp.simplify(solvz / (I * k))
soldel = sp.solve(E00, de1)[0]
print("=> delta1 =", soldel)
# Phi equation from the trace
trace = (Exx + Eyy + Ezz)
trace = trace.subs(de1, soldel).subs(Psi1, solPsi).doit()
ode = sp.simplify(trace)
print("trace eq after eliminating delta, Psi:", sp.collect(sp.expand(ode), Phi1))
# normalise to Phi'' + ...
c2 = ode.coeff(sp.Derivative(Phi1, (eta, 2)))
odeN = sp.expand(ode / c2)
print("ODE: ", odeN, " = 0")
assert sp.simplify(odeN - (sp.Derivative(Phi1, (eta, 2)) + 4 / eta * sp.Derivative(Phi1, eta) + k**2 / 3 * Phi1)) == 0
print("VERIFIED: Phi1'' + (4/eta) Phi1' + (k^2/3) Phi1 = 0")

# conservation equations implied?
C = [linear_part(cons[nu], 1) for nu in range(4)]
subsd = {de1: soldel, Psi1: solPsi, vz: solvz, vx: 0, vy: 0}
for nu in range(4):
    c = C[nu].subs({de1: soldel, vz: solvz, vx: 0, vy: 0}).subs(Psi1, solPsi).doit()
    # eliminate Phi'' and Phi''' using the ODE
    phi2 = sp.solve(odeN, sp.Derivative(Phi1, (eta, 2)))[0]
    c = c.subs(sp.Derivative(Phi1, (eta, 3)), sp.diff(phi2, eta)).subs(sp.Derivative(Phi1, (eta, 2)), phi2).doit()
    c = sp.simplify(c.subs(sp.Derivative(Phi1, (eta, 2)), phi2))
    print("conservation nu=%d after substitution:" % nu, c)
    assert c == 0
print("VERIFIED: conservation equations implied by Einstein equations at first order")

# ---- solutions. x = k eta / sqrt3
x = k * eta / sp.sqrt(3)
zeta = sp.Symbol('zeta')
TPhi = -2 * (sp.sin(x) - x * sp.cos(x)) / x**3
assert sp.simplify(odeN.subs(Phi1, zeta * TPhi).doit()) == 0
print("VERIFIED: T_Phi = -2(sin x - x cos x)/x^3 solves the ODE (regular at 0)")
# second (singular) solution
TPhi_sing = (sp.cos(x) + x * sp.sin(x)) / x**3
assert sp.simplify(odeN.subs(Phi1, TPhi_sing).doit()) == 0
print("second solution (cos x + x sin x)/x^3 verified (singular ~ eta^-3)")
Tdel = sp.simplify(soldel.subs(Psi1, solPsi).subs(Phi1, TPhi).doit())
TV = sp.simplify(V1.subs(Psi1, solPsi).subs(Phi1, TPhi).doit())
print("T_delta =", sp.simplify(Tdel)); print("T_V =", sp.simplify(TV))
# delta_T = (1+delta)^{1/4}-1 -> delta/4 at first order
TdT = Tdel / 4
lim = sp.limit(TdT - TPhi, eta, 0)
print("lim_{eta->0} (delta_T1 - Phi1)/zeta =", lim)
assert lim == 1
print("VERIFIED: zeta_1 = lim (delta_T1 - Phi1);  Phi1(0) =", sp.limit(TPhi, eta, 0))
print("T_deltaT(0) =", sp.limit(TdT, eta, 0), "   T_V ~", sp.series(TV, eta, 0, 3))
# express in x
X = sp.Symbol('x', positive=True)
tox = lambda f: sp.simplify(f.subs(eta, sp.sqrt(3) * X / k))
res = {'T_Phi': tox(TPhi), 'T_delta': tox(Tdel), 'T_deltaT': tox(TdT), 'T_V': tox(TV),
       'T_Phi_eta': TPhi, 'T_delta_eta': Tdel, 'T_V_eta': TV, 'ode': odeN,
       'Psi1': 'Phi1', 'delta1_from_00': soldel, 'V1_from_0i': V1}
for kk in ('T_Phi', 'T_delta', 'T_deltaT', 'T_V'):
    print(kk, "(x) =", res[kk])
with open('first_order.pkl', 'wb') as fh:
    pickle.dump({kk: sp.srepr(vv) if isinstance(vv, sp.Basic) else vv for kk, vv in res.items()}, fh)
with open('first_order.tex', 'w') as fh:
    for kk in ('T_Phi', 'T_delta', 'T_deltaT', 'T_V'):
        fh.write("%s(x) = %s\\\\\n" % (kk, sp.latex(res[kk])))
    fh.write("ODE: %s = 0\\\\\n" % sp.latex(odeN))
# vorticity: v1_i = d_i V -> curl = 0 identically (verified structurally: v_x=v_y=0 in k-frame, v_z = i k V)
print("vorticity: v1 = grad V1 (transverse components vanish) -> curl v1 = 0")
