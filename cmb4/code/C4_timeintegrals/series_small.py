"""series_small.py -- eta->0 series of the second-order fields to O(eta^4) (V2 to O(eta^5)), from the ODE Frobenius
recursion + C1's algebraic constraints, and check against the closed forms (mpmath, 120 digits) at small eta."""
import sympy as sp, mpmath as mp, pickle, sys, time
sys.path.insert(0, '.'); sys.path.insert(0, '../C1_einstein')
from trig_integrate import *
from kernels_api import load, closed_form, eta, k, q, p
t0 = time.time()
S = load()['second']; L = {'Cin': Cin}
SPhi = closed_form(S['K_SPhi']); Kh = closed_form(S['K_h'])
NORD = 7    # keep eta^0..eta^6 in Phi2, h2 (report to eta^4)
def ser(e, n=NORD):
    return sp.expand(sp.series(e, eta, 0, n).removeO())
s = ser(SPhi, NORD - 1)
a = {0: 0, 1: 0}
for n in range(2, NORD):
    a[n] = sp.factor(sp.cancel((s.coeff(eta, n - 2) - k**2 / 3 * a[n - 2]) / (n * (n + 3))))
sh = ser(Kh, NORD - 1); b = {0: 0, 1: 0}
for n in range(2, NORD):
    b[n] = sp.factor(sp.cancel((sh.coeff(eta, n - 2) - k**2 * b[n - 2]) / (n * (n + 1))))
K0 = S['K_Phi2_0']
gh = ser(S['ghat'])
Pser = sum(a[n] * eta**n for n in a)
Phi2s = sp.expand(Pser + K0 * gh)
dPhi2s = sp.diff(Phi2s, eta)
Qs = ser(closed_form(S['K_Q'])); S00s = ser(eta**2 * closed_form(S['K_S00'])); S0zs = ser(closed_form(S['K_S0z']) / (sp.I * k), NORD + 1)
Psi2s = sp.expand(Phi2s - Qs)
delta2s = sp.expand(S00s / 3 - sp.Rational(2, 3) * eta**2 * k**2 * Phi2s - 2 * eta * dPhi2s - 2 * Psi2s)
xx = lambda kk: kk * eta / sp.sqrt(3)
TPh = lambda kk: closed_form(sp.Function('T_Phi')(kk, eta)) if False else -2 * (sp.sin(xx(kk)) - xx(kk) * sp.cos(xx(kk))) / xx(kk)**3
Tdel = lambda kk: -2 * eta**2 * kk**2 * TPh(kk) / 3 - 2 * eta * sp.diff(TPh(kk), eta) - 2 * TPh(kk)
deltaT2s = sp.expand(delta2s / 4 - sp.Rational(3, 16) * ser(Tdel(q) * Tdel(p)))
V2s = sp.expand(-eta**2 / 4 * (2 * (dPhi2s + Psi2s / eta) + S0zs))
out = {}
def report(nm, e, nmax):
    e = sp.expand(e)
    cs = {n: sp.factor(sp.cancel(e.coeff(eta, n))) for n in range(0, nmax + 1)}
    out[nm] = cs
    print(nm + ":")
    for n, c in cs.items():
        if c != 0: print("   eta^%d: %s" % (n, c))
report('P', Pser, 6); report('Phi2', Phi2s, 4); report('Psi2', Psi2s, 4); report('Q2', Qs, 4); report('delta2', delta2s, 4)
report('deltaT2', deltaT2s, 4); report('V2', V2s, 5); report('I_h', sum(b[n] * eta**n for n in b), 6)
pickle.dump({nm: {n: sp.srepr(c) for n, c in cs.items()} for nm, cs in out.items()}, open('series_small_eta.pkl', 'wb'))
print("series built (%.0fs)" % (time.time() - t0))
# ---- check vs closed forms at small eta (120 digits): (closed - series_to_eta^6)/eta^7 must stay bounded (odd powers vanish -> eta^8)
P = sp.parse_expr(pickle.load(open('phi2_closed.pkl', 'rb'))['P'], local_dict=L)
Ih = sp.parse_expr(pickle.load(open('h2_closed.pkl', 'rb'))['I_h'], local_dict=L)
mp.mp.dps = 120
def cin_mp(z):
    z = abs(z)
    return mp.euler + mp.log(z) - mp.ci(z) if z != 0 else mp.mpf(0)
NS = {'Cin': cin_mp, 'Si': mp.si}
fP = sp.lambdify((k, q, p, eta), P, modules=[NS, 'mpmath']); fI = sp.lambdify((k, q, p, eta), Ih, modules=[NS, 'mpmath'])
fPs = sp.lambdify((k, q, p, eta), Pser, 'mpmath'); fIs = sp.lambdify((k, q, p, eta), sum(b[n] * eta**n for n in b), 'mpmath')
print("check closed vs series (P to eta^6, I_h to eta^6):  (closed-series)/eta^8")
for kk, qq, pp in [(mp.mpf(1), mp.mpf('0.8'), mp.mpf('0.75')), (mp.mpf('1.3'), mp.mpf('2.1'), mp.mpf('1.7'))]:
    for ee in [mp.mpf('1e-1'), mp.mpf('1e-2'), mp.mpf('1e-3')]:
        dP = mp.re(fP(kk, qq, pp, ee)) - fPs(kk, qq, pp, ee); dI = mp.re(fI(kk, qq, pp, ee)) - fIs(kk, qq, pp, ee)
        print("   k=%s q=%s p=%s eta=%s:  P: %s   I_h: %s   (P itself %s)" % (kk, qq, pp, mp.nstr(ee, 2), mp.nstr(dP / ee**8, 10), mp.nstr(dI / ee**8, 10), mp.nstr(mp.re(fP(kk, qq, pp, ee)), 8)))
print("done %.0fs" % (time.time() - t0))
