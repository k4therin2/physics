import sympy as sp, mpmath as mp, pickle, sys, time
sys.path.insert(0, '.'); sys.path.insert(0, '../C1_einstein')
from trig_integrate import *
from kernels_api import load, closed_form, eta, k, q, p
L = {'Cin': Cin}
C = {kk: sp.parse_expr(vv, local_dict=L) for kk, vv in pickle.load(open('collinear_closed.pkl', 'rb')).items()}
P = sp.parse_expr(pickle.load(open('phi2_closed.pkl', 'rb'))['P'], local_dict=L)
Ih = sp.parse_expr(pickle.load(open('h2_closed.pkl', 'rb'))['I_h'], local_dict=L)
K = load(); S = K['second']
etap = sp.Symbol('etap', positive=True)
SPhi = closed_form(S['K_SPhi']).subs(eta, etap); Kh = closed_form(S['K_h']).subs(eta, etap)
intPhi = sp.expand(S['G_Phi'] * SPhi); inth = sp.expand(S['G_h'] * Kh)
mp.mp.dps = 50
def quad(integrand, vals):
    f = sp.lambdify((etap,), integrand.subs(vals), 'mpmath')
    return mp.quad(f, mp.linspace(0, mp.mpf(sp.N(vals[eta], 50)), 40), method='gauss-legendre', maxdegree=10)
cases = [('Phi2_q+p=k', intPhi, {k: 1, q: sp.Rational(1, 2), p: sp.Rational(1, 2), eta: 7}),
         ('Phi2_q+p=k', intPhi, {k: sp.Rational(3, 2), q: sp.Rational(1, 5), p: sp.Rational(13, 10), eta: sp.Rational(5, 2)}),
         ('Phi2_q=k+p', intPhi, {k: 1, q: 3, p: 2, eta: 5}),
         ('Phi2_p=k+q', intPhi, {k: 1, q: sp.Rational(3, 4), p: sp.Rational(7, 4), eta: 4}),
         ('h2_q+p=sqrt3k', inth, {k: 1, q: sp.sqrt(3) / 2, p: sp.sqrt(3) / 2, eta: 9}),
         ('h2_q+p=sqrt3k', inth, {k: 1, q: sp.sqrt(3) / 4, p: 3 * sp.sqrt(3) / 4, eta: 3})]
for nm, integrand, vals in cases:
    a = sp.re(sp.N(C[nm].subs(vals), 40)); b = quad(integrand, vals)
    # continuity of the generic Cin-form: evaluate at a point displaced by 1e-12 from the degenerate one (40-digit arithmetic)
    gen = P if nm.startswith('Phi2') else Ih
    vd = dict(vals); vd[p] = vals[p] + sp.Rational(1, 10**12)
    c = sp.re(sp.N(gen.subs(vd), 60))
    print("%s %s: limit-form=%s quad=%s generic(p+1e-12)=%s  rel(lim-quad)=%s" % (nm, {str(s): str(v) for s, v in vals.items()}, sp.N(a, 18), mp.nstr(b, 18), sp.N(c, 18), mp.nstr(abs(mp.mpf(str(a)) - b) / abs(b), 3)))
print("log terms in h2 resonance form:", C['h2_q+p=sqrt3k'].atoms(sp.log))
