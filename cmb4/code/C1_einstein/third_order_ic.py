"""third_order_ic.py -- tadpoles (angle-averaged k->0 limits of second-order scalar kernels) and the
third-order regular-mode amplitude C3 from zeta_nl^(3)(eta->0) = 0 (contracted form)."""
import sympy as sp, pickle, time
from tpoly import *
from fourier import *
t0 = time.time()
S2 = {kk: sp.parse_expr(vv) for kk, vv in pickle.load(open('second_order.pkl', 'rb')).items()}
R2 = {kk: sp.parse_expr(vv) for kk, vv in pickle.load(open('second_order_raw.pkl', 'rb')).items()}
TP = lambda kk: TF['Phi'](kk, eta)
X = lambda kk: kk * eta / sp.sqrt(3)
TPhi_cf = lambda kk: -2 * (sp.sin(X(kk)) - X(kk) * sp.cos(X(kk))) / X(kk)**3
def closed(expr, syms=(q, p, k)):
    expr = expr.doit()
    for kk in syms:
        expr = expr.subs(TP(kk), TPhi_cf(kk))
    return expr.doit()

# ---------- 1. tadpole of the nonlocal Q2: angle-averaged k->0 limit of K_Q(q, k-q)
KQ = closed(S2['K_Q'])
KQ = KQ.subs(p, sp.sqrt(k**2 + q**2 - 2 * k * q * mu))
ser = sp.series(KQ, k, 0, 1).removeO()          # keep k^-2 ... k^0
ser = sp.expand(ser)
avg = {n: sp.simplify(sp.integrate(ser.coeff(k, n), (mu, -1, 1)) / 2) for n in (-2, -1, 0)}
print("angle-averaged coefficients of K_Q at k^-2, k^-1, k^0:", avg)
assert avg[-2] == 0 and avg[-1] == 0
tadQ2 = sp.simplify(avg[0])
print("tadQ2(q,eta) = <Q2> integrand per P(q) =", tadQ2, "   eta->0:", sp.limit(tadQ2, eta, 0))
# ---------- 2. local tadpoles from the mu-form (raw) kernels: k->0, p->q
def tad_local(expr):
    e = closed(expr.subs(k, 0).subs(p, q), syms=(q,))
    e = sp.simplify(e)
    assert not e.has(mu), e
    return e
tadS00 = tad_local(R2['K_S00'])
tad_tr_m_00 = tad_local(sp.expand(R2['ker_tr'] - R2['K_S00']))
tadz2q = tad_local(S2['K_zeta2_quad'].subs(p, q).subs(k, 0)) if False else None
print("tadS00 =", tadS00, "  eta^2 tadS00 ->", sp.limit(eta**2 * tadS00, eta, 0))
tadSPhi = sp.simplify(-tad_tr_m_00 / 6 + sp.diff(tadQ2, eta) / eta)
print("tadSPhi =", tadSPhi, "   eta->0:", sp.limit(tadSPhi, eta, 0))
etap = sp.Symbol('etap', positive=True)
G0 = etap / 3 * (1 - etap**3 / eta**3)
integrand = sp.simplify(G0 * tadSPhi.subs(eta, etap))
tadPhi2_int = sp.integrate(integrand, (etap, 0, eta))
tadPhi2 = sp.simplify(-sp.Rational(8, 9) + tadPhi2_int)
print("tadPhi2(eta) =", tadPhi2)
print("   series:", sp.series(tadPhi2, eta, 0, 4))
# verify ODE: Phi'' + 4 Phi'/eta = tadSPhi (k=0)
assert sp.simplify(sp.diff(tadPhi2, eta, 2) + 4 / eta * sp.diff(tadPhi2, eta) - tadSPhi) == 0
taddelta2 = sp.simplify(eta**2 / 3 * tadS00 - 2 * eta * sp.diff(tadPhi2, eta) - 2 * (tadPhi2 - tadQ2))
print("taddelta2(eta) =", taddelta2, "  eta->0:", sp.limit(taddelta2, eta, 0))
tads = {'tadQ2': tadQ2, 'tadS00': tadS00, 'tadSPhi': tadSPhi, 'tadPhi2': tadPhi2, 'taddelta2': taddelta2}
print("tadpoles done (%.0fs)" % (time.time() - t0))

# ---------- 3. zeta3 contracted kernel and C3
import third_order as TO
Phi = field_series('Phi'); delta = field_series('delta')
zeta_nl = (const(1) + delta).log1p_series() / 4 + (const(1) - Phi * 2).log1p_series() / 2
z3 = zeta_nl.order_part(3) * 6
lin_z3 = linear_part(z3, 3)
ker_z3 = TO.on_shell(TO.contracted_kernel(z3))
print("zeta3 linear part:", lin_z3, " contracted kernel terms:", len(sp.Add.make_args(ker_z3)))
# small-eta models
a_, b_, c_ = sp.symbols('a b c', positive=True)
def at(expr, A, B, C):
    return expr.subs({q: a_, p: b_, k: c_}, simultaneous=True).subs({a_: A, b_: B, c_: C}, simultaneous=True)
KPhi2_0 = S2['K_Phi2_0']
# K_SPhi at eta=0 in invariant form
SPhi0 = sp.limit(closed(S2['K_SPhi']), eta, 0)
SPhi0 = sp.simplify(SPhi0)
print("K_SPhi(eta=0) invariant:", SPhi0)
A, B, C = q, k, p
KPhi2_small = at(KPhi2_0, A, B, C) + eta**2 * (sp.Rational(3, 10) * at(SPhi0, A, B, C) - C**2 * at(KPhi2_0, A, B, C) / 30)
Kh_small = sp.Rational(8, 9) * eta**2
def small_eta(expr, order=6):
    """replace abstract second-order functions and T_Phi by their small-eta expansions, p -> |k-q|."""
    expr = expr.subs(TO.KPhi2(A, B, C, eta), KPhi2_small).subs(TO.Kh(A, B, C, eta), Kh_small)
    for nm, fn in tads.items():
        expr = expr.subs(sp.Function(nm)(eta), sp.series(fn, eta, 0, order).removeO())
    expr = expr.doit()
    Xs = lambda kk: sp.series(TPhi_cf(kk), eta, 0, order).removeO()
    for kk in (q, p, k):
        expr = expr.subs(sp.Derivative(TP(kk), eta), sp.diff(Xs(kk), eta)).subs(TP(kk), Xs(kk))
    expr = expr.doit()
    return expr
Phi3, de3 = amp('Phi', 3), amp('delta', 3)
Q3 = sp.Function('Q3')(eta)
C3 = sp.Symbol('C3')
# limits needed: Q3c(0), lim eta^2 S00c
Q3_0 = sp.limit(sp.expand(small_eta(TO.K3['Q'])).subs(p, sp.sqrt(k**2 + q**2 - 2 * k * q * mu)), eta, 0)
print("Q3c(eta->0) =", sp.factor(Q3_0), "(%.0fs)" % (time.time() - t0))
S00_0 = sp.limit(sp.expand(eta**2 * small_eta(TO.K3['S00'])).subs(p, sp.sqrt(k**2 + q**2 - 2 * k * q * mu)), eta, 0)
print("lim eta^2 S00c^(3) =", sp.factor(S00_0))
SPhi3_0 = sp.limit(sp.expand(small_eta(TO.K3['SPhi'])).subs(p, sp.sqrt(k**2 + q**2 - 2 * k * q * mu)), eta, 0)
print("S_Phi^(3)c(eta->0) finite:", SPhi3_0.is_finite, SPhi3_0 != sp.zoo)
# delta3c(0) = S00_0/3 - 2 (C3 - Q3_0);  Phi3c(0) = C3;  zeta3 = lin_z3 + ker_z3
z3_0 = lin_z3.subs(de3, S00_0 / 3 - 2 * (C3 - Q3_0)).subs(Phi3, C3)
kz = sp.limit(sp.expand(small_eta(ker_z3)).subs(p, sp.sqrt(k**2 + q**2 - 2 * k * q * mu)), eta, 0)
z3_0 = sp.expand(z3_0 + kz)
C3_sol = sp.factor(sp.solve(z3_0, C3)[0])
print("==> C3(k; q, mu) = Phi3c(k, eta->0) kernel per P(q):", C3_sol)
C3_inv = sp.factor(sp.simplify(C3_sol.subs(mu, (k**2 + q**2 - p**2) / (2 * k * q))))
print("    invariant form:", C3_inv)
# tadpole part of C3 (terms proportional to T(k)-type local pieces) and the angle average
C3_avg = sp.simplify(sp.integrate(C3_sol, (mu, -1, 1)) / 2)
print("    angle average over mu:", C3_avg)
out = {'tads': {kk: sp.srepr(vv) for kk, vv in tads.items()}, 'C3': sp.srepr(C3_sol), 'C3_inv': sp.srepr(C3_inv),
       'Q3_0': sp.srepr(Q3_0), 'S00_0': sp.srepr(S00_0), 'SPhi0_2nd': sp.srepr(SPhi0), 'KPhi2_small': sp.srepr(KPhi2_small),
       'zeta3_lin': sp.srepr(lin_z3), 'zeta3_contracted': sp.srepr(ker_z3)}
pickle.dump(out, open('third_order_ic.pkl', 'wb'))
with open('third_order_ic.tex', 'w') as fh:
    for kk, vv in tads.items():
        fh.write("\\paragraph{%s}\n\\begin{dmath}\n%s\n\\end{dmath}\n" % (kk, sp.latex(vv)))
    fh.write("\\paragraph{C3}\n\\begin{dmath}\n%s\n\\end{dmath}\n" % sp.latex(C3_sol))
    fh.write("\\paragraph{C3 invariant}\n\\begin{dmath}\n%s\n\\end{dmath}\n" % sp.latex(C3_inv))
print("saved third_order_ic.pkl/.tex (%.0fs)" % (time.time() - t0))
