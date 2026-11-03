"""
V2, script 2: squeezed-limit predictions for K_2 and K_3 (SPEC conventions) from the adiabatic-mode
coordinate transformation, plus numerical evaluation.

Logic (see report_V2.md):
  * zeta_L = c + q.x (x measured from the observer's label point, the origin) is, to O(q) and to all
    orders in c for the c part, a diffeomorphism of FRW (script adiabatic_mode_rd.py, PART B/D).
  * T_obs(n) = T_o exp(-S), S = int sigma_ee dl, the emission surface T_phys=T_*, the observer worldline
    (tagged fluid element), its clock (proper time from the T_phys=T_i crossing), its FW/canonical
    tetrad and the monopole normalisation are all defined covariantly.  Hence
        Theta[zeta_L + zeta_S](n) = Theta[zeta~_S ; observer at x~_o](R n),
    zeta~_S(x~) = zeta_S(x(x~)) (scalar relabelling: proven from the uniform-density-slice metric),
    x~_o = -(eta_i^2/6) e^{-c} q (tagged worldline in tilde coords), R = 1 + O(q^2, q S).
  * Plane-wave rule: Theta^(1)[P(x) e^{ip.x}] = P(-i d_p) K_1(p).  This turns the relabelling into
    differential operators on K_1:
        K_2(k1,k2;n) = -D K_1(k2,n) - S_{k1} K_1(k2,n) + (eta_i^2/6)(k1.k2) K_1(k2,n) + O(k1^2)
        D = k2.d_{k2},   S_p = (k2.d_{k2})(p.d_{k2}) - (1/2)(p.k2) Lap_{k2}.
"""
import sympy as sp
import mpmath as mp
import random

print("=== 1. Operator identities ===")
# (a) plane-wave rule & relabelling -> operators.  Test: Theta^(1)[zeta] := int_k zeta_k K(k).  For
# zeta = P(x) e^{ip.x}, zeta_k = P(i d_k)(2pi)^3 delta(k-p) -> Theta = P(-i d_p) K(p) (integration by parts).
# relabelling  delta zeta = -xi(x).grad e^{ipx} = -i (p.xi(x)) e^{ipx}  with p held fixed inside P:
p1, p2, p3, q1, q2, q3, c = sp.symbols('p1 p2 p3 q1 q2 q3 c', real=True)
X = sp.symbols('X1:4', real=True)   # stands for x^a, to be replaced by -i d/dp_a
pv, qv = [p1, p2, p3], [q1, q2, q3]
Kf = sp.Function('K')(p1, p2, p3)
def apply_poly(P):
    """P polynomial in X (with p held fixed): replace X^a -> -i d/dp_a acting on K."""
    P = sp.Poly(sp.expand(P), *X)
    out = 0
    for monom, coeff in P.terms():
        e = Kf
        for a, m in enumerate(monom):
            for _ in range(m):
                e = -sp.I*sp.diff(e, pv[a])
        out += coeff*e
    return sp.expand(out)
xdot = lambda u, v: sum(u[i]*v[i] for i in range(3))
xi_dil = [c*X[a] for a in range(3)]
xi_sct = [xdot(qv, X)*X[a] - sp.Rational(1, 2)*xdot(X, X)*qv[a] for a in range(3)]
Ddil = apply_poly(-sp.I*xdot(pv, xi_dil))
Dsct = apply_poly(-sp.I*xdot(pv, xi_sct))
Dop = sum(pv[a]*sp.diff(Kf, pv[a]) for a in range(3))
Sop = lambda u: sum(pv[j]*sp.diff(sum(u[i]*sp.diff(Kf, pv[i]) for i in range(3)), pv[j]) for j in range(3)) \
    - sp.Rational(1, 2)*xdot(u, pv)*sum(sp.diff(Kf, pv[i], 2) for i in range(3))
print(" dilation:  delta K = -c D K ?", sp.simplify(Ddil + c*Dop) == 0)
print(" SCT:       delta K = i S_q K ?", sp.simplify(Dsct - sp.I*Sop(qv)) == 0)
# with zeta_L = zeta_{k1} e^{i k1.x}: c = zeta_{k1}, q = i k1 zeta_{k1}  =>  i S_q -> -S_{k1}
assert sp.simplify(Ddil + c*Dop) == 0 and sp.simplify(Dsct - sp.I*Sop(qv)) == 0

# (b) the operator S_p acting on F(k, mu), mu = khat.n  (a = p.khat, b = p.n):
#   S_p F = (a k/2) F_kk - a F_k + (b - mu a) F_kmu + ((2 a mu - b)/k) F_mu - (a (1-mu^2)/(2k)) F_mumu
kk, mm = sp.symbols('k mu', positive=True)
k1s, k2s, k3s = sp.symbols('k1 k2 k3', real=True)
kv = [k1s, k2s, k3s]
Ftest = sp.sin(2*kk)*mm**3 + sp.exp(kk*mm)*kk**2 + sp.cos(3*kk*mm/2)/kk
Kabs = sp.sqrt(k1s**2 + k2s**2 + k3s**2); Mu = k3s/Kabs        # n = zhat
f = Ftest.subs({kk: Kabs, mm: Mu})
dd = lambda e, i: sp.diff(e, kv[i])
pdf = sum(pv[i]*dd(f, i) for i in range(3))
lhs = sum(kv[j]*dd(pdf, j) for j in range(3)) - sp.Rational(1, 2)*xdot(pv, kv)*sum(dd(dd(f, i), i) for i in range(3))
a_, b_ = xdot(pv, kv)/Kabs, p3
def Dk(i, j):
    e = Ftest
    if i: e = sp.diff(e, kk, i)
    if j: e = sp.diff(e, mm, j)
    return e.subs({kk: Kabs, mm: Mu})
rhs = (a_*Kabs/2)*Dk(2, 0) - a_*Dk(1, 0) + (b_ - Mu*a_)*Dk(1, 1) + ((2*a_*Mu - b_)/Kabs)*Dk(0, 1) - (a_*(1 - Mu**2)/(2*Kabs))*Dk(0, 2)
random.seed(1)
errs = [abs(sp.N((lhs - rhs).subs({s: sp.Rational(random.randint(-9, 9), 7) for s in (k1s, k2s, k3s, p1, p2, p3)}), 30)) for _ in range(5)]
print(" S_p in (k,mu) variables: max residual at random points =", max(errs))
assert max(errs) < 1e-25

# (c) composite c then q (to O(cq)): kernel = (1 + i S_q) e^{-cD} [e^{ik.x~_o} K]; check e^{-cD} S_q = e^{c} S_q e^{-cD}
# (S_q lowers the k-degree by one) on a test function, to O(c)
g = sp.Function('g')
Ktest = sp.exp(-(p1**2 + 2*p2**2 + 3*p3**2)/5)*(p1 + p2*p3)
lam = sp.Symbol('lam', positive=True)
Sq_on = lambda e: sum(pv[j]*sp.diff(sum(qv[i]*sp.diff(e, pv[i]) for i in range(3)), pv[j]) for j in range(3)) \
    - sp.Rational(1, 2)*xdot(qv, pv)*sum(sp.diff(e, pv[i], 2) for i in range(3))
lhs_c = Sq_on(Ktest).subs({p1: lam*p1, p2: lam*p2, p3: lam*p3}, simultaneous=True)
rhs_c = Sq_on(Ktest.subs({p1: lam*p1, p2: lam*p2, p3: lam*p3}, simultaneous=True))/lam
print(" homogeneity S_q[K(lam k)] = lam (S_qK)(lam k):", sp.simplify(lhs_c - rhs_c) == 0)

print("\n=== 2. Observer worldline, tetrad, emission surface in the tilde frame ===")
eta, etai = sp.symbols('eta eta_i', positive=True)
x1, x2, x3 = sp.symbols('x1 x2 x3', real=True)
xs = [x1, x2, x3]
zL = c + xdot(qv, xs)
xi = [zL*xs[i] - sp.Rational(1, 2)*xdot(xs, xs)*qv[i] - eta**2/6*qv[i] for i in range(3)]
# tagged worldline, label 0 at eta_init (eta_init = eta_i(1 + c/3) at the observer, irrelevant at O(q)):
etainit = etai*(1 + c/3)
xo = [(eta**2 - etainit**2)/6*qv[i] for i in range(3)]          # v_L = eta q/3
xt = [sp.expand((xs[i] + xi[i]).subs(dict(zip(xs, xo)), simultaneous=True)) for i in range(3)]
# keep O(q) only (drop q^2): q-degree
qq = sp.Symbol('qq')
xt1 = [sp.expand(sp.series(e.subs({q1: qq*q1, q2: qq*q2, q3: qq*q3}), qq, 0, 2).removeO().subs(qq, 1)) for e in xt]
print(" x~_o(eta) =", [sp.factor(sp.series(e, c, 0, 1).removeO()) for e in xt1], "(time independent: fluid at rest in tilde frame)")
# antisymmetric part of the spatial Jacobian at the observer = rotation of the canonical triad
J = sp.Matrix(3, 3, lambda i, j: sp.diff(xs[i] + xi[i], xs[j]))
A = ((J - J.T)/2).subs(dict(zip(xs, xo)), simultaneous=True)
print(" antisymmetric Jacobian at observer (O(q) part):", sp.simplify(A.applyfunc(lambda e: sp.series(e.subs({q1: qq*q1, q2: qq*q2, q3: qq*q3}), qq, 0, 3).removeO().coeff(qq, 1))))
print(" -> no rotation of the canonical triad at O(q); the boost part is absorbed because the triad is built orthogonal")
print("    to the PHYSICAL fluid velocity (same worldline in both charts) -> no net aberration.")

print("\n=== 3. First-order kernel K_1(k,n) (SPEC conventions, observer at origin) ===")
k, mu, eo, es = sp.symbols('k mu eta_o eta_s', positive=True)
def TPhi(e):
    x = k*e/sp.sqrt(3)
    return -2*(sp.sin(x) - x*sp.cos(x))/x**3
H = 1/eta
TP = TPhi(eta)
TdT = -(1/(6*H**2))*(k**2*TP + 3*H*(sp.diff(TP, eta) + H*TP))   # delta_T = delta/4
Tv = -(sp.diff(TP, eta) + H*TP)/(2*H**2)                          # v^i = d_i vpot, vpot_k = Tv zeta_k
chi = eo - es
loc_e = sp.exp(sp.I*k*mu*chi)*((TdT + TP).subs(eta, es) - sp.I*k*mu*Tv.subs(eta, es))
loc_o = sp.I*k*mu*Tv.subs(eta, eo)
isw_integrand = 2*sp.diff(TP, eta)*sp.exp(sp.I*k*mu*(eo - eta))
j0 = lambda z: sp.sin(z)/z
j1 = lambda z: (sp.sin(z) - z*sp.cos(z))/z**2
mono_loc = j0(k*chi)*(TdT + TP).subs(eta, es) + k*Tv.subs(eta, es)*j1(k*chi)
mono_isw_integrand = 2*sp.diff(TP, eta)*j0(k*(eo - eta))
# small-k check: F - F_0 = O(k^2): symbolic series of the transfer functions + numerical scaling test below
Floc = loc_e + loc_o
print(" T_dT+T_Phi, T_v/eta at small k:", sp.series((TdT + TP), k, 0, 3).removeO(), ",", sp.series(Tv/eta, k, 0, 3).removeO())
print(" ISW integrand small-k:", sp.series(sp.diff(TP, eta), k, 0, 2).removeO(), " (O(k^2))")
print(" => local part: (T_dT+T_Phi)(1 + i k mu chi) - i k mu eta_s/3 + i k mu eta_o/3 = -1/3 + O(k^2) since chi = eta_o - eta_s;")
print("    the constant -1/3 is removed by the monopole (and is cancelled by the observer's own [dT+Psi]_o in S_1).")

# numerical machinery
mp.mp.dps = 25
def lamb(expr):
    return sp.lambdify((k, mu, eo, es, eta), expr, modules='mpmath')
derivs = {}
for (i, j) in [(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (0, 2)]:
    def dkm(e):
        if i: e = sp.diff(e, k, i)
        if j: e = sp.diff(e, mu, j)
        return e
    derivs[(i, j)] = (lamb(dkm(Floc)), lamb(dkm(isw_integrand)),
                      lamb(sp.diff(mono_loc, k, i)) if j == 0 else None,
                      lamb(sp.diff(mono_isw_integrand, k, i)) if j == 0 else None)
def Fd(i, j, kv_, muv, eov, esv):
    """d^i_k d^j_mu of the monopole-subtracted K_1."""
    fl, fi, ml, mi = derivs[(i, j)]
    val = fl(kv_, muv, eov, esv, 0) + mp.quad(lambda t: fi(kv_, muv, eov, esv, t), mp.linspace(esv, eov, 9))
    if j == 0:
        val -= ml(kv_, muv, eov, esv, 0) + mp.quad(lambda t: mi(kv_, muv, eov, esv, t), mp.linspace(esv, eov, 9))
    return val

def K1(kv_, muv, eov, esv):
    return Fd(0, 0, kv_, muv, eov, esv)
def DK1(kv_, muv, eov, esv):
    return kv_*Fd(1, 0, kv_, muv, eov, esv)
def SK1(a, b, kv_, muv, eov, esv):
    """S_p K_1 for p with p.khat2 = a, p.n = b."""
    return (a*kv_/2)*Fd(2, 0, kv_, muv, eov, esv) - a*Fd(1, 0, kv_, muv, eov, esv) + (b - muv*a)*Fd(1, 1, kv_, muv, eov, esv) \
        + ((2*a*muv - b)/kv_)*Fd(0, 1, kv_, muv, eov, esv) - (a*(1 - muv**2)/(2*kv_))*Fd(0, 2, kv_, muv, eov, esv)

def K2_squeezed(k1vec, k2vec, nvec, eov, esv, eiv=0):
    """Prediction for K_2(k1,k2;n) through O(k1)."""
    k2m = mp.sqrt(sum(x**2 for x in k2vec)); k2h = [x/k2m for x in k2vec]
    muv = sum(k2h[i]*nvec[i] for i in range(3))
    a = sum(k1vec[i]*k2h[i] for i in range(3)); b = sum(k1vec[i]*nvec[i] for i in range(3))
    K = K1(k2m, muv, eov, esv)
    order0 = -DK1(k2m, muv, eov, esv)
    order1 = -SK1(a, b, k2m, muv, eov, esv) + (mp.mpf(eiv)**2/6)*a*k2m*K
    return order0, order1, K

if __name__ == "__main__":
    print("\n=== 4. Numerical predictions ===")
    eov, esv = mp.mpf(1), mp.mpf('0.2')
    for kv_ in (mp.mpf('1e-2'), mp.mpf('2e-2')):
        print(" small-k scaling: |K_1(k, mu=0.7)| at k =", kv_, ":", mp.nstr(abs(K1(kv_, mp.mpf('0.7'), eov, esv)), 10))
    # finite-difference check of D vs eta-scaling: D K_1 = (eta_o d_eta_o + eta_s d_eta_s) K_1
    kv_, muv = mp.mpf(12), mp.mpf('0.3')
    lhs = DK1(kv_, muv, eov, esv)
    rhs = mp.diff(lambda s: K1(kv_, muv, eov*s, esv*s), 1)
    print(" D K_1 =", mp.nstr(lhs, 12), "  (eta_o d_eo + eta_s d_es) K_1 =", mp.nstr(rhs, 12))
    # finite-difference check of the (k,mu) form of S_p against direct differentiation in 3D
    k2vec = [mp.mpf(3), mp.mpf(-4), mp.mpf(11)]; nvec = [mp.mpf(0), mp.mpf(0), mp.mpf(1)]; pvec = [mp.mpf('0.7'), mp.mpf('-0.2'), mp.mpf('0.5')]
    def Kvec(v):
        m = mp.sqrt(sum(x**2 for x in v)); return K1(m, v[2]/m, eov, esv)
    h = mp.mpf('1e-4')
    def grad(fun, v):
        return [(fun([v[j] + (h if j == i else 0) for j in range(3)]) - fun([v[j] - (h if j == i else 0) for j in range(3)]))/(2*h) for i in range(3)]
    pdK = lambda v: sum(pvec[i]*g for i, g in enumerate(grad(Kvec, v)))
    kdpdK = sum(k2vec[i]*g for i, g in enumerate(grad(pdK, k2vec)))
    lap = sum((Kvec([k2vec[j] + (h if j == i else 0) for j in range(3)]) - 2*Kvec(k2vec) + Kvec([k2vec[j] - (h if j == i else 0) for j in range(3)]))/h**2 for i in range(3))
    direct = kdpdK - sum(pvec[i]*k2vec[i] for i in range(3))*lap/2
    k2m = mp.sqrt(sum(x**2 for x in k2vec))
    viaform = SK1(sum(pvec[i]*k2vec[i] for i in range(3))/k2m, pvec[2], k2m, k2vec[2]/k2m, eov, esv)
    print(" S_p K_1 direct 3D finite differences:", mp.nstr(direct, 8), "  (k,mu)-form:", mp.nstr(viaform, 8))
    print("\n Table: eta_o=1, eta_s=0.2; k1 = k1hat*|k1| with |k1|-coefficient shown (K_2 = O0 + O1 + O(k1^2))")
    print(" k2   mu2    K_1                      O0 = -D K_1                O1/|k1| (k1 || n)          O1/|k1| (k1 || k2)")
    for k2m in (mp.mpf(5), mp.mpf(20)):
        for muv in (mp.mpf('-0.6'), mp.mpf('0.1'), mp.mpf('0.8')):
            k2vec = [k2m*mp.sqrt(1 - muv**2), 0, k2m*muv]
            o0, o1n, K = K2_squeezed([0, 0, 1], k2vec, nvec, eov, esv)
            _, o1k, _ = K2_squeezed([x/k2m for x in k2vec], k2vec, nvec, eov, esv)
            print("", mp.nstr(k2m, 3), mp.nstr(muv, 2), mp.nstr(K, 8), mp.nstr(o0, 8), mp.nstr(o1n, 8), mp.nstr(o1k, 8))
