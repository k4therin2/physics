import sympy as sp, pickle, sys, random
sys.path.insert(0, '.'); sys.path.insert(0, '../C1_einstein')
from trig_integrate import *
from kernels_api import eta, k, q, p
L = {'Cin': Cin}
d = pickle.load(open('asymptotics.pkl', 'rb'))
P = sp.parse_expr(pickle.load(open('phi2_closed.pkl', 'rb'))['P'], local_dict=L)
Ih = sp.parse_expr(pickle.load(open('h2_closed.pkl', 'rb'))['I_h'], local_dict=L)
P_as = sp.parse_expr(d['Phi2_asym'], local_dict=L); Ih_as = sp.parse_expr(d['h2_asym'], local_dict=L)
random.seed(5)
cfg = []
for i in range(4):
    kk = sp.Rational(random.randint(3, 20), 10); qq = sp.Rational(random.randint(1, 60), 20) * kk
    mu = sp.Rational(random.randint(-19, 19), 20); pp = sp.sqrt(kk**2 + qq**2 - 2 * kk * qq * mu); ee = sp.Rational(1000 + random.randint(0, 500), 10) / kk
    cfg.append((kk, qq, mu, pp))
kk, qq, mu, pp = cfg[1]
print("config:", kk, qq, mu, sp.N(pp), "u+v=", sp.N((qq + pp) / kk), " frequencies h2:", [sp.N(kk + s1 * qq / sp.sqrt(3) + s2 * pp / sp.sqrt(3), 4) for s1 in (1, -1) for s2 in (1, -1)])
for X in [50, 100, 200, 400, 800]:
    ee = sp.Rational(X) / kk
    vals = {k: kk, q: qq, p: pp, eta: ee}
    full = sp.re(sp.N(Ih.subs(vals), 40)); fa = sp.N(Ih_as.subs(vals), 30)
    fullP = sp.re(sp.N(P.subs(vals), 50)); faP = sp.N(P_as.subs(vals), 30)
    print("x=%d: h2 full=%s asym=%s diff*x^3=%s | Phi2 full=%s asym=%s diff*x^4=%s" % (X, sp.N(full, 12), sp.N(fa, 12), sp.N((full - fa) * X**3, 6), sp.N(fullP, 12), sp.N(faP, 12), sp.N((fullP - faP) * X**4, 6)))
