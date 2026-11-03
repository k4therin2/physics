"""asymptotics_compact.py -- compact leading eta->infinity forms of Phi2 (incl. homogeneous piece) and I_h,
from asymptotics.pkl (verified by probe_asym2.py: remainder scales as eta^-(nmax+1) up to x = 1e5)."""
import sympy as sp, pickle, sys
sys.path.insert(0, '.'); sys.path.insert(0, '../C1_einstein')
from trig_integrate import *
from kernels_api import eta, k, q, p, load
L = {'Cin': Cin}
S = load()['second']
d = pickle.load(open('asymptotics.pkl', 'rb'))
Pl = sp.parse_expr(d['Phi2_lead'], local_dict=L); Il = sp.parse_expr(d['h2_lead'], local_dict=L)
s3 = sp.sqrt(3); x = k * eta / s3
l = sp.symbols('l1:5')
logmap = {sp.log((k + q + p) / s3): l[0], sp.log((q + p - k) / s3): l[1], sp.log((k + q - p) / s3): l[2], sp.log((k + p - q) / s3): l[3]}
def packlogs(e):
    e = sp.expand(sp.expand_log(e, force=True))
    rep = {}
    for a in e.atoms(sp.log):
        arg = sp.expand(a.args[0] * s3)
        for key, sym in logmap.items():
            if sp.expand(key.args[0] * s3 - arg) == 0: rep[a] = sym
    return sp.expand(e.xreplace(rep))
Lam = sp.Symbol('Lambda')     # Lambda := ln[(k+q+p)(q+p-k)/((k+q-p)(k+p-q))]
e = packlogs(Pl)
# homogeneous piece: K_Phi2_0 * ghat, ghat = 3 u1 -> -3 cos x / x^2 + O(x^-3)
gl = sp.limit(S['ghat'].subs(eta, sp.sqrt(3) * sp.Symbol('X') / k) * sp.Symbol('X')**2 + 3 * sp.cos(sp.Symbol('X')), sp.Symbol('X'), sp.oo) if False else None
res = {}
for trig in (sp.sin(x), sp.cos(x)):
    c = sp.expand(e.coeff(trig) * eta**2)
    cl = [sp.factor(c.coeff(li)) for li in l]
    rest = sp.factor(sp.expand(c - sum(ci * li for ci, li in zip(cl, l))))
    print(trig, "log coeffs:", [sp.simplify(ci / cl[0]) for ci in cl])
    res[trig] = (cl[0], rest)
    print("   coeff of Lambda:", cl[0]); print("   rest:", rest)
# full Phi2 leading = P_lead + K_Phi2_0 * (-3 cos x / x^2)
K0 = S['K_Phi2_0']
cosc = sp.factor(sp.expand(res[sp.cos(x)][1] + K0 * (-9 / k**2)))
print("Phi2 total cos(x)/eta^2 non-log coefficient:", cosc)
Ie = packlogs(Il)
for trig in (sp.sin(k * eta), sp.cos(k * eta)):
    c = sp.expand(Ie.coeff(trig) * eta)
    print("h2", trig, sp.factor(c))
