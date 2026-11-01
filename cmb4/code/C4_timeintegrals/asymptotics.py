"""asymptotics.py -- eta -> infinity forms of the particular integrals P (Phi2) and I_h (h2).
Cin(z) -> EulerGamma + ln|z| - sin z/z - cos z/z^2 + ...,   Si(z) -> (pi/2) sgn z - cos z/z - sin z/z^2 + ...
Signs of the frequencies are fixed by the triangle inequalities k <= q+p, p <= k+q, q <= k+p:
   Phi2: (k+q+p), (k+q-p), (k-q+p) >= 0,  (k-q-p) <= 0.     h2: k+(q+p)/s3, k+(q-p)/s3, k-(q-p)/s3 > 0;  k-(q+p)/s3 has either sign.
Result kept to relative order 1/eta^2 beyond leading (i.e. terms eta^-2, eta^-3, eta^-4 for Phi2; eta^-1, eta^-2 for h2)."""
import sympy as sp, pickle, sys, time, random
sys.path.insert(0, '.'); sys.path.insert(0, '../C1_einstein')
from trig_integrate import *
from kernels_api import load, eta, k, q, p
L = {'Cin': Cin}
P = sp.parse_expr(pickle.load(open('phi2_closed.pkl', 'rb'))['P'], local_dict=L)
Ih = sp.parse_expr(pickle.load(open('h2_closed.pkl', 'rb'))['I_h'], local_dict=L)
t0 = time.time()
s3 = sp.sqrt(3)
def sgn_phys(w):
    """sign of a frequency combination in the physical domain, or None if indefinite"""
    w = sp.expand(w * s3)
    ck, cq, cp = [w.coeff(s) for s in (k, q, p)]
    if (ck, cq, cp) in [(1, 1, 1), (1, 1, -1), (1, -1, 1), (-1, 1, 1)]: return 1
    if (ck, cq, cp) in [(-1, -1, -1), (-1, -1, 1), (-1, 1, -1), (1, -1, -1)]: return -1
    return None
def asym(e, nmax, sgn=sgn_phys):
    e = sp.expand(e.rewrite(sp.cos)) if e.has(sp.I) else sp.expand(e)
    e = sp.expand(sp.re(e)) if e.has(sp.I) else e
    rep = {}
    M = 2
    for f in e.atoms(Cin):
        z = f.args[0]; w = sp.cancel(z / eta); s = sgn(w)
        absw = s * w if s is not None else sp.Abs(w)
        rep[f] = sp.EulerGamma + sp.log(absw * eta) - sum((-1)**j * sp.factorial(2*j) * sp.sin(z) / z**(2*j+1) for j in range(M)) \
                 - sum((-1)**j * sp.factorial(2*j+1) * sp.cos(z) / z**(2*j+2) for j in range(M))
    for f in e.atoms(Si):
        z = f.args[0]; w = sp.cancel(z / eta); s = sgn(w)
        sg = s if s is not None else sp.sign(w)
        rep[f] = sp.pi / 2 * sg - sum((-1)**j * sp.factorial(2*j) * sp.cos(z) / z**(2*j+1) for j in range(M)) \
                 + sum((-1)**j * sp.factorial(2*j+1) * sp.sin(z) / z**(2*j+2) for j in range(M))
    e = sp.expand(e.xreplace(rep))
    # keep powers eta^n with n >= -nmax ; log(eta) terms must cancel
    out = 0
    for term in sp.Add.make_args(e):
        pw = 0
        for f in sp.Mul.make_args(term):
            if f == eta: pw += 1
            elif f.is_Pow and f.base == eta: pw += int(f.exp)
        if pw >= -nmax: out += term
    out = sp.expand(sp.expand_log(out, force=True))
    assert not out.has(sp.log(eta)), "log(eta) survives"
    return out
def sgn_h(w):
    w = sp.expand(w * s3); ck, cq, cp = [w.coeff(s) for s in (k, q, p)]
    ck = ck / s3
    if (ck, cq, cp) in [(1, 1, 1), (1, 1, -1), (1, -1, 1)]: return 1
    if (ck, cq, cp) in [(-1, -1, -1), (-1, -1, 1), (-1, 1, -1)]: return -1
    return None      # k - (q+p)/sqrt3 and its negative: indefinite -> sign()
Ih_as = asym(Ih, 2, sgn_h)
print("h2 asymptotics done (%.0fs), nops %d" % (time.time() - t0, sp.count_ops(Ih_as)))
P_as = asym(P, 3)
print("Phi2 asymptotics done (%.0fs), nops %d" % (time.time() - t0, sp.count_ops(P_as)))
# leading terms
def leading(e, n):
    out = 0
    for term in sp.Add.make_args(sp.expand(e)):
        pw = 0
        for f in sp.Mul.make_args(term):
            if f == eta: pw += 1
            elif f.is_Pow and f.base == eta: pw += int(f.exp)
        if pw == n: out += term
    return out
Ih_lead = sp.collect(leading(Ih_as, -1), [sp.sin(k*eta), sp.cos(k*eta)], sp.factor)
P_lead = sp.collect(leading(P_as, -2), P_as.atoms(sp.sin, sp.cos), sp.factor)
print("h2 leading (1/eta):", Ih_lead)
print("Phi2 leading (1/eta^2):", P_lead)
# KT (25) control
u, v, x = sp.symbols('u v x', positive=True)
KT25 = 3 * (u**2 + v**2 - 3) / (4 * u**3 * v**3 * x) * (sp.sin(x) * (-4 * u * v + (u**2 + v**2 - 3) * sp.log(sp.Abs((3 - (u + v)**2) / (3 - (u - v)**2))))
                                                   - sp.pi * (u**2 + v**2 - 3) * sp.Heaviside(v + u - s3) * sp.cos(x))
KT25 = 4 / k**2 * KT25.subs({u: p / k, v: q / k, x: k * eta})
random.seed(5)
print("control vs KT(25):")
for i in range(4):
    kk = sp.Rational(random.randint(3, 20), 10); qq = sp.Rational(random.randint(1, 60), 20) * kk
    mu = sp.Rational(random.randint(-19, 19), 20); pp = sp.sqrt(kk**2 + qq**2 - 2 * kk * qq * mu); ee = sp.Rational(1000 + random.randint(0, 500), 10) / kk
    vals = {k: kk, q: qq, p: pp, eta: ee}
    a = sp.N(Ih_lead.subs(vals), 30); b = sp.N(KT25.subs(vals), 30)
    full = sp.re(sp.N(Ih.subs(vals), 40)); fa = sp.N(Ih_as.subs(vals), 30)
    print("  u+v=%s x=%s: lead=%s KT25=%s | full=%s asym(3 orders)=%s" % (sp.N((qq+pp)/kk, 4), sp.N(kk*ee, 5), sp.N(a, 12), sp.N(b, 12), sp.N(full, 15), sp.N(fa, 15)))
print("Phi2: full vs asymptotic at large eta:")
for i in range(4):
    kk = sp.Rational(random.randint(3, 20), 10); qq = sp.Rational(random.randint(1, 60), 20) * kk
    mu = sp.Rational(random.randint(-19, 19), 20); pp = sp.sqrt(kk**2 + qq**2 - 2 * kk * qq * mu); ee = sp.Rational(1000 + random.randint(0, 500), 10) / kk
    vals = {k: kk, q: qq, p: pp, eta: ee}
    full = sp.re(sp.N(P.subs(vals), 50)); fa = sp.N(P_as.subs(vals), 30); fl = sp.N(P_lead.subs(vals), 30)
    print("  x=%s: full=%s  asym=%s  lead=%s  rel(asym)=%s rel(lead)=%s" % (sp.N(kk*ee, 5), sp.N(full, 15), sp.N(fa, 15), sp.N(fl, 15), sp.N(abs(fa/full-1), 3), sp.N(abs(fl/full-1), 3)))
pickle.dump({'Phi2_asym': sp.srepr(P_as), 'Phi2_lead': sp.srepr(P_lead), 'h2_asym': sp.srepr(Ih_as), 'h2_lead': sp.srepr(Ih_lead)}, open('asymptotics.pkl', 'wb'))
print("saved asymptotics.pkl (%.0fs)" % (time.time() - t0))
