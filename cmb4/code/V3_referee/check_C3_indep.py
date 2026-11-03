"""check_C3_indep.py -- V3 claim (5): independent derivation of the third-order initial amplitude.
Discrete three-mode realisation (zeta1 e^{iq.x}, zeta2 e^{-iq.x}, zeta3 e^{ik.x}) in MY E^{(3)}, with all first- and
second-order fields replaced by their exact small-eta Taylor series (transfer functions expanded; Phi2 = C(0) +
eta^2 [alpha S_Phi(0) - c^2 C(0)/30] with the CORRECT alpha = 1/10 (also run with C1's 3/10); h2 = (8/9) eta^2 Lambda qq;
tadpoles = angle-averaged k->0 limits computed HERE from the (verified) finite-k kernels).  Then
   zeta_nl^(3)(eta->0) = 0  (zeta_nl = 1/4 ln(1+delta) + 1/2 ln(1-2 Phi), expanded HERE)  =>  C3sym,
which must equal C1's C3(k,q,mu) + C3(k,q,-mu) [p -> |k-q|, |k+q|] (image sum: discrete realisation contains both)."""
import sys, random, sympy as sp
import pw_tools
from pw_tools import *
import third_tools as TT
from third_tools import K, KA, eta, k, q, p, mu
E = E_mine()
NS = 9
alpha = sp.Rational(1, 10) if len(sys.argv) < 2 else sp.Rational(sys.argv[1])
# optional zero-mode (tadpole) shifts to test convention dependence: argv[3], argv[4] = shift of tadPhi2, tadPsi2 (constants)
sPhi = sp.Rational(sys.argv[3]) if len(sys.argv) > 3 else 0; sPsi = sp.Rational(sys.argv[4]) if len(sys.argv) > 4 else 0
def Tser_expr(kk):
    xv = kk * t / sp.sqrt(3)
    return sp.series(-2 * (sp.sin(xv) - xv * sp.cos(xv)) / xv**3, t, 0, NS).removeO()
pw_tools.T_Phi = Tser_expr
from fourier import TF
def ser_form(expr, syms=(q, p, k)):
    expr = expr.doit()
    for kk in syms:
        expr = expr.subs(TF['Phi'](kk, eta), Tser_expr(kk).subs(t, eta))
    return expr.doit()
def at(expr, A, B, C):
    return ser_form(expr).subs({q: A, p: B, k: C}, simultaneous=True).subs(eta, t)
TT.at = at
# ---- tadpoles: angle-averaged k->0 limits of the finite-k small-eta kernels
kk_, qq_, mm_ = sp.symbols('kk qq mm', positive=True)
def tad(expr_kqp):
    e = sp.expand(ser_form(expr_kqp).subs(eta, t).subs(p, sp.sqrt(kk_**2 + qq_**2 - 2 * kk_ * qq_ * mm_)).subs({k: kk_, q: qq_}))
    e = sp.series(e, t, 0, 4).removeO()
    out = 0
    for n in range(0, 4):
        c = e.coeff(t, n)
        if c == 0: continue
        cav = sp.integrate(sp.series(c, kk_, 0, 1).removeO().expand(), (mm_, -1, 1)) / 2
        out += sp.simplify(cav) * t**n
    return sp.expand(out)
S2 = K['second']
C0 = S2['K_Phi2_0']; SPhi0 = K['third_ic']['SPhi0_2nd']
tQ = tad(S2['K_Q'])
tC0 = tad(C0); tS0 = tad(SPhi0)
tPhi = tC0 + t**2 * alpha * tS0          # ghat(k eta) -> 1 as k -> 0
tS00 = tad(eta**2 * S2['K_S00']) / t**2
tPhi = tPhi + sPhi; tQ = tQ + sPhi - sPsi   # Psi = Phi - Q
tdelta = sp.expand(t**2 * tS00 / 3 - 2 * t * sp.diff(tPhi, t) - 2 * (tPhi - tQ))   # k^2 Phi term: angle avg of k^2 C0 -> 0
print("tadpoles (per P(q), q-symbol qq): tadQ2 =", tQ, "; tadPhi2 =", tPhi, "; taddelta2 =", tdelta)
# ---- second-order Phi2 series for generic pair
a_, b_, c_ = sp.symbols('a_ b_ c_', positive=True)
def Phi2_ser(A, B, C):
    C0v = C0.subs({q: A, p: B, k: C}, simultaneous=True); S0v = SPhi0.subs({q: A, p: B, k: C}, simultaneous=True)
    return C0v + t**2 * (alpha * S0v - C**2 * C0v / 30)
def subs_abstract(expr, qmag):
    reps = {}
    for f in expr.atoms(sp.core.function.AppliedUndef):
        nm = f.func.__name__
        if nm == 'KPhi2': reps[f] = Phi2_ser(*f.args[:3])
        elif nm == 'Kh': reps[f] = sp.Rational(8, 9) * t**2
        elif nm == 'tadPhi2': reps[f] = tPhi.subs(qq_, qmag)
        elif nm == 'tadQ2': reps[f] = tQ.subs(qq_, qmag)
        elif nm == 'taddelta2': reps[f] = tdelta.subs(qq_, qmag)
    return expr.subs(reps).doit()
z1, z2, z3 = sp.symbols('zeta1 zeta2 zeta3')
C3 = K['third_ic']['C3']
Csym = sp.Symbol('Csym')
worst = 0
for seed in range(int(sys.argv[2]) if len(sys.argv) > 2 else 3):
    pt = rand_point(500 + seed); km, qm = pt['k'], pt['q']
    f = TT.three_mode_fields(pt, z1, z2, z3, amp=1)
    f = {nm: subs_abstract(v_, qm) for nm, v_ in f.items()}
    def S3(a, b):
        return 6 * coeff_at_origin(substitute(E[('E', a, b, 3)], f), (z1, z2, z3))
    lim = lambda e_: sp.nsimplify(sp.series(sp.expand(e_), t, 0, 1).removeO()) if False else sp.expand(sp.series(sp.expand(e_), t, 0, 1).removeO())
    S00_0 = lim(t**2 * S3(0, 0))
    Q3_0 = lim((S3(3, 3) - (S3(1, 1) + S3(2, 2)) / 2) / km**2)
    # zeta_nl to eps^3 with discrete fields at eta -> 0
    eps = sp.Symbol('eps')
    d1 = lim(f['delta1']); P1 = lim(f['Phi1']); d2 = lim(f['delta2']); P2 = lim(f['Phi2'])
    d3 = (S00_0 / 3 - 2 * (Csym - Q3_0)) * z1 * z2 * z3     # delta3 = eta^2 S00/3 - 2/3 eta^2 k^2 Phi3 - 2 eta Phi3' - 2 Psi3 at eta=0
    P3 = Csym * z1 * z2 * z3
    dl = eps * d1 + eps**2 * d2 / 2 + eps**3 * d3 / 6; Ph = eps * P1 + eps**2 * P2 / 2 + eps**3 * P3 / 6
    zeta = sp.log(1 + dl) / 4 + sp.log(1 - 2 * Ph) / 2
    z3c = 6 * sp.series(zeta, eps, 0, 4).removeO().coeff(eps, 3)
    z3c = coeff_at_origin(sp.expand(z3c), (z1, z2, z3))
    sol = sp.solve(z3c, Csym)[0]
    pm_m = sp.sqrt(km**2 + qm**2 - 2 * km * qm * pt['mu']); pm_p = sp.sqrt(km**2 + qm**2 + 2 * km * qm * pt['mu'])
    c1 = C3.subs({k: km, q: qm, mu: pt['mu']}) + C3.subs({k: km, q: qm, mu: -pt['mu']})
    d = abs(sp.N(sol - c1, 30)); worst = max(worst, d)
    print("point", pt, " C3sym mine =", sp.N(sol, 20), " C1 image sum =", sp.N(c1, 20), " |diff| = %.2e" % d)
print("zero-mode shifts (Phi2,Psi2) =", sPhi, sPsi, "; alpha =", alpha, " WORST |C3sym mine - C1|:", worst)
