"""second_order_solve.py -- Green's functions, vector sector (vorticity-free), tensor (KT check),
superhorizon amplitude C_Phi from zeta_nl^(2)=0, comparison with PUB (3.12) and BMR (7.30/7.31)."""
import sympy as sp, pickle, random, time
from tpoly import *
from fourier import *
t0 = time.time()
R = pickle.load(open('second_order_raw.pkl', 'rb')); R = {kk: sp.sympify(vv) if isinstance(vv, str) else vv for kk, vv in R.items()}
R = {kk: eval(sp.srepr(vv)) if False else vv for kk, vv in R.items()}
# re-parse srepr strings into live objects
R = {kk: sp.parse_expr(vv) if isinstance(vv, str) else vv for kk, vv in pickle.load(open('second_order_raw.pkl', 'rb')).items()}
def g(name): return R[name]
TP = lambda kk: TF['Phi'](kk, eta)
X = lambda kk: kk * eta / sp.sqrt(3)
TPhi_cf = lambda kk: -2 * (sp.sin(X(kk)) - X(kk) * sp.cos(X(kk))) / X(kk)**3
def closed(expr):
    expr = expr.doit()
    for kk in (q, p, k):
        expr = expr.subs(TP(kk), TPhi_cf(kk))
    return expr.doit()
def numcheck(expr, nvals=4, tol=1e-25, extra=()):
    """evaluate expr at random numerical points (mu, q, k, eta) with p = |k-q|."""
    for _ in range(nvals):
        vals = {mu: sp.Rational(random.randint(-90, 90), 100), q: sp.Rational(random.randint(10, 300), 100),
                k: sp.Rational(random.randint(10, 300), 100), eta: sp.Rational(random.randint(10, 300), 100)}
        vals[p] = sp.sqrt(vals[k]**2 + vals[q]**2 - 2 * vals[k] * vals[q] * vals[mu])
        for s, v in extra: vals[s] = v
        val = sp.N(expr.subs(vals), 40)
        if abs(val) > tol:
            return False, val
    return True, 0

# ===================== 1. Green's functions
Phi2 = sp.Function('Phi2')(eta); etap = sp.Symbol('etap', positive=True); S = sp.Function('S')
x = X(k); xp = k * etap / sp.sqrt(3)
u1 = (sp.sin(x) - x * sp.cos(x)) / x**3; u2 = (sp.cos(x) + x * sp.sin(x)) / x**3
W = sp.simplify(u1 * sp.diff(u2, eta) - sp.diff(u1, eta) * u2)
print("Wronskian u1 u2' - u1' u2 =", W)
G_Phi = sp.simplify((u2 * u1.subs(eta, etap) - u1 * u2.subs(eta, etap)) / W.subs(eta, etap))
print("G_Phi(eta,eta') =", G_Phi)
# verify: L[int_0^eta G S] = S  <=> G(eta,eta)=0, d_eta G|_{eta'=eta} = 1, L_eta G = 0
assert sp.simplify(G_Phi.subs(etap, eta)) == 0
assert sp.simplify(sp.diff(G_Phi, eta).subs(etap, eta)) == 1
LG = sp.diff(G_Phi, eta, 2) + 4 / eta * sp.diff(G_Phi, eta) + k**2 / 3 * G_Phi
assert sp.simplify(LG) == 0
print("VERIFIED G_Phi: Phi2 = int_0^eta G_Phi(eta,eta') S_Phi(eta') deta' + C_Phi ghat(k eta),  ghat = 3 u1 -> 1")
G_h = etap / eta * sp.sin(k * (eta - etap)) / k
assert sp.simplify(G_h.subs(etap, eta)) == 0 and sp.simplify(sp.diff(G_h, eta).subs(etap, eta)) == 1
assert sp.simplify(sp.diff(G_h, eta, 2) + 2 / eta * sp.diff(G_h, eta) + k**2 * G_h) == 0
print("VERIFIED G_h = (eta'/eta) sin(k(eta-eta'))/k  for h'' + 2h'/eta + k^2 h = S")
G_om = etap**2 / eta**2
print("G_om = (eta'/eta)^2 for om' + 2 om/eta = S  (regular solution; homogeneous mode eta^-2 excluded)")

# ===================== 2. Vector sector: vorticity-free condition
d = load_polys('E.pkl'); u_dn = d['u_dn']
uhat_x = u_dn[1] * eta_pow(-1)                # uhat_i = u_i / a
delta_s = field_series('delta')
enth = (const(1) + delta_s).sqrt_series().sqrt_series()     # (1+delta)^{1/4} = 1 + delta_T  (enthalpy/ T-hat weight)
# first test: curl(uhat)=0  (fails, see report); correct statement: curl[(1+delta_T) uhat_i] = 0 (relativistic Kelvin)
P2_naive = uhat_x.order_part(2) * 2
P2 = (enth * uhat_x).order_part(2) * 2
lin_ux = linear_part(P2, 2); ker_ux = quadratic_kernel(P2)
print("uhat_x^(2) linear part:", lin_ux)
v2x, omx = amp('v1', 2), amp('om1', 2)
# transverse 1-form must vanish (k x uhat = 0):  v2_x + om2_x + K_ux = 0  and 0x constraint
S0x = sp.Symbol('S0x'); Kux = sp.Symbol('Kux')
sol = sp.solve([lin_ux + Kux, g('lin0x') + S0x], [v2x, omx], dict=True)[0]
print("vorticity-free + 0i constraint =>", sol)
om_alg = sol[omx].subs({S0x: g('K_S0x'), Kux: ker_ux}); om_alg = sp.expand(om_alg.doit())
# on-shell first-order substitution, then check the om ODE from E_xz
import second_order as SO   # reuse on_shell (re-runs quickly)
om_alg = SO.on_shell(om_alg)
resid = SO.on_shell(sp.diff(om_alg, eta) - g('sol_omp').subs(sp.Symbol('Sxz'), g('K_Sxz')).subs(omx, om_alg).doit())
ok, val = numcheck(closed(resid))
print("omega2 (algebraic, from vorticity=0 + 0i) satisfies the ij-vector ODE:", ok, val)
assert ok
sq = sp.sqrt(1 - mu**2)
K_om = sp.expand(sp.cancel(om_alg / (q * sq)))     # om2_i = int P_ij(k) q_j K_om(q,p) zeta zeta
assert not K_om.has(sq)
print("K_om (mu form) =", K_om)
v2T = SO.on_shell(sol[v2x].subs({S0x: g('K_S0x'), Kux: ker_ux}).doit())
K_vT = sp.expand(sp.cancel(v2T / (q * sq)))
print("K_vT (transverse velocity, mu form) =", K_vT)

# ===================== 3. Tensor: compare with Kohri-Terada / Baumann
# our eq: h_p'' + 2 h_p'/eta + k^2 h_p = -kerT  with (xx-yy) component;  h2_xx - h2_yy = 2 h_p
# KT target (our conventions, h2 = h^KT): h2'' + 2h2'/eta - lap h2 = 4 Lambda_ij,kl [2 d_k Phi1 d_l Phi1 + d_k(eta Phi1' + Phi1) d_l(eta Phi1' + Phi1)]
kv = sp.Matrix([0, 0, k]); qv = sp.Matrix(QV.vec); pv = sp.Matrix(PV.vec)
Pm = sp.eye(3) - kv * kv.T / k**2
def Lam(Sm):
    return sp.Matrix(3, 3, lambda i, j: sum(Pm[i, a] * Pm[j, b] * Sm[a, b] for a in range(3) for b in range(3)) - sp.Rational(1, 2) * Pm[i, j] * sum(Pm[a, b] * Sm[a, b] for a in range(3) for b in range(3)))
A = lambda kk: eta * sp.Derivative(TP(kk), eta) + TP(kk)
Sm = sp.Matrix(3, 3, lambda i, j: -(2 * (qv[i] * pv[j] + pv[i] * qv[j]) / 2 * TP(q) * TP(p) + (qv[i] * pv[j] + pv[i] * qv[j]) / 2 * A(q) * A(p)))
KT_rhs = 4 * Lam(Sm)
KT_xxyy = sp.expand(KT_rhs[0, 0] - KT_rhs[1, 1])          # = 2 * RHS for h_p
ours_xxyy = sp.expand(-2 * g('kerT'))                         # RHS of (h2_xx - h2_yy)'' + ... ; h2_xx-h2_yy = 2 h_p
diff = SO.on_shell(sp.expand(KT_xxyy - ours_xxyy))
ok, val = numcheck(closed(diff))
print("Kohri-Terada/Baumann tensor source reproduced:", ok, val)
assert ok
K_h = sp.expand(sp.cancel(-2 * g('kerT') / (q**2 * (1 - mu**2))))   # h2_ij = Lambda_ij,kl q_k q_l * int G_h K_h
print("K_h (mu form) =", K_h)

# ===================== 4. zeta_nl^(2) and the superhorizon amplitude
Phi = field_series('Phi'); delta = field_series('delta')
zeta_nl = (const(1) + delta).log1p_series() / 4 + (const(1) - Phi * 2).log1p_series() / 2
z2 = zeta_nl.order_part(2) * 2
print("zeta_nl^(2) =", z2.to_sympy())
z3 = zeta_nl.order_part(3) * 6
print("zeta_nl^(3) =", z3.to_sympy())
# Fourier: zeta2 kernel = delta2/4 - Phi2 + quadratic kernel of (first-order)^2 terms
K_z2quad = quadratic_kernel(z2)
K_z2quad = SO.on_shell(K_z2quad)
print("zeta2 quadratic kernel:", K_z2quad)
# superhorizon limits with closed forms
K_Q0 = sp.limit(closed(g('K_Q')), eta, 0)
print("K_Q(eta->0) =", sp.factor(K_Q0))
S00lim = sp.limit(closed(sp.expand(eta**2 * g('K_S00'))), eta, 0)
print("lim eta^2 K_S00 =", sp.factor(S00lim))
# S_Phi small-eta behaviour
SPhi_cf = closed(g('K_SPhi'))
ser = sp.series(SPhi_cf, eta, 0, 2).removeO()
print("K_SPhi series at eta->0:", sp.simplify(ser))
C = sp.Symbol('C_Phi')
# delta2(0) = lim[eta^2 S00/3 - (2/3)eta^2 k^2 Phi2 - 2 eta Phi2' - 2 Psi2] = S00lim/3 - 2 (C - K_Q0)
delta2_0 = S00lim / 3 - 2 * (C - K_Q0)
zeta2_0 = delta2_0 / 4 - C + sp.limit(closed(K_z2quad), eta, 0)
C_sol = sp.solve(zeta2_0, C)[0]
K_Phi2_0 = sp.factor(C_sol)
print("==> C_Phi(k;q,p) = Phi2(k, eta->0) kernel per zeta_q zeta_p =", K_Phi2_0)
# PUB (3.12b) in our conventions: Phi2(0) = -2 Phi1^2 - Delta^-1[dPhi dPhi] + 3 Delta^-2 dd[dPhi dPhi], Phi1 -> -2/3 zeta
qdp = k * q * mu - q**2; kdq = k * q * mu; kdp = k**2 - k * q * mu
PUB_Phi2 = sp.Rational(4, 9) * (-2 - qdp / k**2 + 3 * kdq * kdp / k**4)
PUB_QmP = sp.Rational(4, 9) * (4 + 3 * qdp / k**2 - 9 * kdq * kdp / k**4)    # Psi2 - Phi2 (ours) = PUB Phi2-Psi2
print("PUB (3.12b) Phi2(0):", sp.factor(PUB_Phi2), "   difference ours-PUB:", sp.simplify(K_Phi2_0 - PUB_Phi2))
print("PUB (3.12c-b) Psi2-Phi2 (ours) :", sp.factor(PUB_QmP), "  difference ours(-K_Q0) - PUB:", sp.simplify(-K_Q0 - PUB_QmP))
F = kdq * kdp / k**2 - qdp / 3
print("BMR F/k^2 form of our Q2(0) per Phi1Phi1:", sp.simplify(K_Q0 / sp.Rational(4, 9)), " vs  -4 + 9F/k^2:", sp.simplify(K_Q0 / sp.Rational(4, 9) - (-4 + 9 * F / k**2)))
print("BMR (7.30) in linear variables would be -4 + 33 F/k^2 ; difference  =", sp.simplify((-4 + 33 * F / k**2) - K_Q0 / sp.Rational(4, 9)))

# ===================== 5. invariant (k,q,p) forms and save
def inv(e):
    e = sp.expand(e).subs(mu, (k**2 + q**2 - p**2) / (2 * k * q))
    return sp.collect(sp.expand(e), [TP(q), TP(p), sp.Derivative(TP(q), eta), sp.Derivative(TP(p), eta)], sp.factor)
out = {}
for nm in ('K_Q', 'K_S00', 'K_S0z', 'K_SPhi'):
    out[nm] = inv(g(nm))
out['K_om'] = inv(K_om); out['K_vT'] = inv(K_vT); out['K_h'] = inv(K_h)
out['K_delta2_expr'] = g('sol_de2'); out['K_v2z_expr'] = g('sol_v2z')
out['K_Phi2_0'] = inv(K_Phi2_0); out['K_Q0'] = inv(K_Q0); out['S00lim'] = inv(S00lim)
out['K_zeta2_quad'] = inv(K_z2quad)
out['G_Phi'] = G_Phi; out['G_h'] = G_h; out['G_om'] = G_om; out['ghat'] = 3 * u1
out['zeta2_position'] = z2.to_sympy(); out['zeta3_position'] = z3.to_sympy()
out['u1'] = u1; out['u2'] = u2
pickle.dump({kk: sp.srepr(vv) for kk, vv in out.items()}, open('second_order.pkl', 'wb'))
with open('second_order.tex', 'w') as fh:
    fh.write("% second-order kernels, per zeta_q zeta_p, p = k - q, frame-independent form in magnitudes k,q,p\n")
    for kk, vv in out.items():
        fh.write("\\paragraph{%s}\n\\begin{dmath}\n%s\n\\end{dmath}\n" % (kk.replace('_', '\\_'), sp.latex(vv)))
print("saved second_order.pkl / .tex (%.1fs)" % (time.time() - t0))
for nm in ('K_Q', 'K_SPhi', 'K_om', 'K_h', 'K_Phi2_0'):
    print(nm, "=", out[nm])
