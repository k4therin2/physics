import sympy as sp, mpmath as mp, pickle, sys, random, time
sys.path.insert(0, '.'); sys.path.insert(0, '../C1_einstein')
from trig_integrate import *
from kernels_api import load, closed_form, eta, k, q, p
d = pickle.load(open('phi2_closed.pkl', 'rb'))
P = sp.parse_expr(d['P'], local_dict={'Cin': Cin})
K = load(); S = K['second']
etap = sp.Symbol('etap', positive=True)
SPhi = closed_form(S['K_SPhi']).subs(eta, etap)
integrand = sp.expand(S['G_Phi'] * SPhi)
mp.mp.dps = 60
random.seed(7)
worst = 0; t0 = time.time()
pts = []
for i in range(8):
    kk = sp.Rational(random.randint(3, 20), 10); qq = sp.Rational(random.randint(1, 60), 20) * kk
    mu = sp.Rational(random.randint(-19, 19), 20)
    pp = sp.sqrt(kk**2 + qq**2 - 2 * kk * qq * mu); ee = sp.Rational(random.randint(2, 200), 10) / kk
    pts.append((kk, qq, mu, pp, ee))
pts.append((sp.S(1), sp.Rational(1, 2), sp.S(1), sp.Rational(1, 2), sp.S(7)))          # collinear q+p = k (resonance w=0 of Phi2)
pts.append((sp.S(1), sp.S(3), sp.S(1), sp.S(2), sp.S(5)))                               # collinear k+p = q
pts.append((sp.S(1), sp.sqrt(3)/2, 1/sp.sqrt(3), sp.sqrt(3)/2, sp.S(9)))               # q+p = sqrt3 k (h2 resonance; nothing special for Phi2)
for (kk, qq, mu, pp, ee) in pts:
    vals = {k: kk, q: qq, p: pp, eta: ee}
    a = sp.re(sp.N(P.subs(vals), 50))
    f = sp.lambdify((etap,), integrand.subs(vals), 'mpmath')
    b = mp.quad(f, mp.linspace(0, mp.mpf(sp.N(ee, 60)), 40), method='gauss-legendre', maxdegree=10)
    err = abs(mp.mpf(str(sp.N(a, 40))) - b) / (abs(b) + mp.mpf('1e-30'))
    worst = max(worst, err)
    print("k=%s q/k=%s mu=%s keta=%s : closed=%s quad=%s rel=%s  (%.0fs)" % (kk, sp.N(qq / kk, 4), sp.N(mu, 4), sp.N(kk * ee, 4), sp.N(a, 18), mp.nstr(b, 18), mp.nstr(err, 3), time.time() - t0))
print("WORST rel err:", mp.nstr(worst, 3))
