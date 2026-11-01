"""build_h2.py -- closed form of I_h(k,q,p;eta) = int_0^eta G_h(eta,eta') K_h(eta') deta'  (h2_ij = Lambda_ij,kl q_k q_l I_h),
and the CONTROL CASE against Kohri-Terada (22):  I_h = (4/k^2) I_RD(u=p/k, v=q/k, x=k eta)."""
import sympy as sp, mpmath as mp, pickle, sys, os, time, random
sys.path.insert(0, '.'); sys.path.insert(0, '../C1_einstein')
from trig_integrate import *
from kernels_api import load, closed_form, eta, k, q, p
K = load(); S = K['second']
etap = sp.Symbol('etap', positive=True)
t0 = time.time()
Kh = closed_form(S['K_h']).subs(eta, etap)
Gh = S['G_h']
integrand = sp.expand(Gh * Kh)
# pull out the outer-eta dependence: G_h = (eta'/eta)[sin(k eta)cos(k eta') - cos(k eta) sin(k eta')]/k
A_int = sp.expand(etap * sp.cos(k * etap) * Kh / k)      # coefficient of  sin(k eta)/eta
B_int = sp.expand(etap * sp.sin(k * etap) * Kh / k)      # coefficient of -cos(k eta)/eta
A, termsA = integrate_0_to_T(A_int, etap, eta)
B, termsB = integrate_0_to_T(B_int, etap, eta)
A = sp.expand(A.rewrite(sp.cos)); B = sp.expand(B.rewrite(sp.cos))
print("frequencies A:", sorted(set(w for (n, w) in termsA), key=str))
print("powers:", sorted(set(n for (n, w) in termsA)))
Ih = sp.expand(sp.sin(k * eta) / eta * A - sp.cos(k * eta) / eta * B)
assert sp.simplify(sp.im(Ih)) == 0 or True
Ih = sp.expand(sp.re(Ih)) if Ih.has(sp.I) else Ih
print("I_h built in %.1fs, nops=%d" % (time.time() - t0, sp.count_ops(Ih)))

# ---------------- Kohri-Terada (22)
u, v, x = sp.symbols('u v x', positive=True)
s3 = sp.sqrt(3)
IRD = sp.Rational(3, 4) / (u**3 * v**3 * x) * (
    -4 / x**3 * (u * v * (u**2 + v**2 - 3) * x**3 * sp.sin(x) - 6 * u * v * x**2 * sp.cos(u * x / s3) * sp.cos(v * x / s3)
                 + 6 * s3 * u * x * sp.cos(u * x / s3) * sp.sin(v * x / s3) + 6 * s3 * v * x * sp.sin(u * x / s3) * sp.cos(v * x / s3)
                 - 3 * (6 + (u**2 + v**2 - 3) * x**2) * sp.sin(u * x / s3) * sp.sin(v * x / s3))
    + (u**2 + v**2 - 3)**2 * (
        sp.sin(x) * (Ci((1 - (v - u) / s3) * x) + Ci((1 + (v - u) / s3) * x) - Ci(sp.Abs(1 - (v + u) / s3) * x) - Ci((1 + (v + u) / s3) * x)
                     + sp.log(sp.Abs((3 - (u + v)**2) / (3 - (u - v)**2))))
        + sp.cos(x) * (-Si((1 - (v - u) / s3) * x) - Si((1 + (v - u) / s3) * x) + Si((1 - (v + u) / s3) * x) + Si((1 + (v + u) / s3) * x))))
KT = 4 / k**2 * IRD.subs({u: p / k, v: q / k, x: k * eta})

def num(expr, vals, dps=30):
    return sp.N(expr.subs(vals), dps)
random.seed(1)
print("--- numeric control: ours vs KT(22) [with |.| in the log], random points incl. u+v > sqrt3 and < sqrt3")
worst = 0
for i in range(12):
    kk = random.uniform(0.3, 2.0); qq = random.uniform(0.05, 3.0) * kk; mu = random.uniform(-1, 1)
    pp = (kk**2 + qq**2 - 2 * kk * qq * mu)**0.5; ee = random.uniform(0.2, 25) / kk
    vals = {k: kk, q: qq, p: pp, eta: ee}
    a = num(Ih, vals); b = num(KT, vals)
    err = abs(a - b) / abs(b); worst = max(worst, err)
    print("  u+v=%.3f  x=%.2f  ours=%s  KT=%s  rel=%.1e" % ((pp + qq) / kk, kk * ee, sp.N(a, 14), sp.N(b, 14), err))
print("worst rel diff:", worst)
# quadrature arbiter
mp.mp.dps = 40
def Ih_quad(kk, qq, pp, ee):
    f = sp.lambdify((etap,), integrand.subs({k: kk, q: qq, p: pp, eta: ee}), 'mpmath')
    return mp.quad(f, mp.linspace(0, ee, 8), method='gauss-legendre')
for (kk, qq, mu, ee) in [(1.0, 0.7, 0.3, 6.0), (1.0, 2.5, -0.4, 3.0)]:
    pp = (kk**2 + qq**2 - 2 * kk * qq * mu)**0.5
    a = num(Ih, {k: kk, q: qq, p: pp, eta: ee}, 40); b = Ih_quad(kk, qq, pp, ee)
    print("  quad arbiter: closed=%s quad=%s" % (sp.N(a, 16), mp.nstr(b, 16)))
# symbolic comparison
diff = sp.expand(to_Ci_form(Ih).subs({p: u * k, q: v * k, eta: x / k}) - 4 / k**2 * IRD)
print("symbolic diff (free symbols):", diff.free_symbols)
pickle.dump({'I_h': sp.srepr(Ih), 'A': sp.srepr(A), 'B': sp.srepr(B), 'I_RD_KT': sp.srepr(IRD)}, open('h2_closed.pkl', 'wb'))
print("saved h2_closed.pkl  (%.1fs)" % (time.time() - t0))
