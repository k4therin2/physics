"""check_series_ode.py -- (1) eta->0 series of the closed forms (regularity: no 1/eta^n, no log eta), compared with the
Frobenius series from the ODE;  (2) ODE residual by substitution of the closed form, evaluated at random points (30 digits)."""
import sympy as sp, mpmath as mp, pickle, sys, time, random
sys.path.insert(0, '.'); sys.path.insert(0, '../C1_einstein')
from trig_integrate import *
from kernels_api import load, closed_form, eta, k, q, p
K = load(); S = K['second']
L = {'Cin': Cin}
P = sp.parse_expr(pickle.load(open('phi2_closed.pkl', 'rb'))['P'], local_dict=L)
Ih = sp.parse_expr(pickle.load(open('h2_closed.pkl', 'rb'))['I_h'], local_dict=L)
SPhi = closed_form(S['K_SPhi']); Kh = closed_form(S['K_h'])
t0 = time.time()

def min_power(e):
    e = sp.expand(e)
    ps = set()
    for term in sp.Add.make_args(e):
        pw = 0
        for f in sp.Mul.make_args(term):
            if f == eta: pw += 1
            elif f.is_Pow and f.base == eta: pw += int(f.exp)
        ps.add(pw)
    return min(ps), max(ps)

def series_closed(e, order):
    """series of e in eta to O(eta^order) by replacing every f(w eta) by its Taylor polynomial."""
    e = sp.expand(e.rewrite(sp.cos)) if e.has(sp.I) else sp.expand(e)
    e = sp.expand(sp.re(e)) if e.has(sp.I) else e
    mn, mx = min_power(e)
    N = order - mn + 1
    rep = {}
    for f in e.atoms(Cin): rep[f] = cin_series(f.args[0], N)
    for f in e.atoms(Si): rep[f] = si_series(f.args[0], N)
    for f in e.atoms(sp.sin): rep[f] = sp.series(f, eta, 0, N).removeO()
    for f in e.atoms(sp.cos): rep[f] = sp.series(f, eta, 0, N).removeO()
    assert not e.atoms(sp.log)
    e = sp.expand(e.xreplace(rep))
    poly = {}
    for term in sp.Add.make_args(e):
        pw = 0; c = sp.S.One
        for f in sp.Mul.make_args(term):
            if f == eta: pw += 1
            elif f.is_Pow and f.base == eta: pw += int(f.exp)
            else: c *= f
        if pw < order:
            poly[pw] = poly.get(pw, 0) + c
    return {pw: sp.factor(sp.cancel(c)) for pw, c in sorted(poly.items())}, (mn, mx)

# ---- Frobenius series from the ODEs
def frob_phi(order):
    s = sp.series(SPhi, eta, 0, order - 1).removeO()
    sc = {n: sp.cancel(s.coeff(eta, n)) for n in range(0, order - 1)}
    a = {0: 0, 1: 0}
    for n in range(2, order):
        a[n] = sp.cancel((sc.get(n - 2, 0) - k**2 / 3 * a.get(n - 2, 0)) / (n * (n + 3)))
    return {n: sp.factor(v) for n, v in a.items()}
def frob_h(order):
    s = sp.series(Kh, eta, 0, order - 1).removeO()
    sc = {n: sp.cancel(s.coeff(eta, n)) for n in range(0, order - 1)}
    a = {0: 0, 1: 0}
    for n in range(2, order):
        a[n] = sp.cancel((sc.get(n - 2, 0) - k**2 * a.get(n - 2, 0)) / (n * (n + 1)))
    return {n: sp.factor(v) for n, v in a.items()}

res = {}
for nm, e, frob in (('h2', Ih, frob_h), ('Phi2', P, frob_phi)):
    ser, (mn, mx) = series_closed(e, 5)
    print("%s: eta powers in closed form %d..%d; series coefficients (%.0fs):" % (nm, mn, mx, time.time() - t0))
    fr = frob(5)
    ok = True
    for pw in sorted(set(ser) | set(fr)):
        c1 = ser.get(pw, 0); c2 = fr.get(pw, 0)
        d = sp.simplify(c1 - c2)
        ok &= (d == 0)
        print("   eta^%d: %s   | Frobenius: %s  | diff: %s" % (pw, c1, c2, d))
    print("%s regular (no negative powers, no log) and matches ODE Frobenius series: %s" % (nm, ok))
    res[nm] = {pw: sp.srepr(c) for pw, c in ser.items()}
pickle.dump(res, open('series_small_eta.pkl', 'wb'))

# ---- ODE residual by substitution (numerical, 30 digits)
print("ODE residuals:")
dP = sp.diff(P, eta); ddP = sp.diff(dP, eta)
resid_phi = ddP + 4 / eta * dP + k**2 / 3 * P - SPhi
dI = sp.diff(Ih, eta); ddI = sp.diff(dI, eta)
resid_h = ddI + 2 / eta * dI + k**2 * Ih - Kh
random.seed(3)
for i in range(4):
    kk = sp.Rational(random.randint(3, 20), 10); qq = sp.Rational(random.randint(1, 60), 20) * kk
    mu = sp.Rational(random.randint(-19, 19), 20); pp = sp.sqrt(kk**2 + qq**2 - 2 * kk * qq * mu)
    ee = sp.Rational(random.randint(2, 200), 10) / kk
    vals = {k: kk, q: qq, p: pp, eta: ee}
    rphi = sp.re(sp.N(resid_phi.subs(vals), 30)); rh = sp.re(sp.N(resid_h.subs(vals), 30))
    sphi = sp.N(SPhi.subs(vals), 30); sh = sp.N(Kh.subs(vals), 30)
    print("  k eta=%s: |resid_Phi|/|S_Phi| = %s ,  |resid_h|/|K_h| = %s" % (sp.N(kk * ee, 4), sp.N(abs(rphi / sphi), 3), sp.N(abs(rh / sh), 3)))
print("done %.0fs" % (time.time() - t0))
