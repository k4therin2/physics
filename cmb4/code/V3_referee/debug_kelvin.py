import sys, pickle, sympy as sp
from pw_tools import *
C1 = '/home/user/physics/cmb4/code/C1_einstein'; sys.path.insert(0, C1)
import kernels_api as KA; K = KA.load()
from fourier import eta, k, q, p, mu
raw = pickle.load(open('indep_E.pkl', 'rb'))
u_x2 = sp.parse_expr(raw[('u_dn', 1, 2)]); u_x1 = sp.parse_expr(raw[('u_dn', 1, 1)])
print("u_x1 =", u_x1); print("u_x2 =", u_x2)
zq, zp = sp.symbols('zeta_q zeta_p')
for seed in range(2):
    pt = rand_point(seed)
    kv, qv, pv = vec3(pt['k'], pt['q'], pt['mu']); pm = sp.sqrt(pt['k']**2 + pt['q']**2 - 2*pt['k']*pt['q']*pt['mu'])
    f = {}
    for nm, val in mode(zq, qv, pt['q']).items(): f[nm] = f.get(nm, 0) + val
    for nm, val in mode(zp, pv, pm).items(): f[nm] = f.get(nm, 0) + val
    qperp = pt['q']*sp.sqrt(1-pt['mu']**2)
    quad_engine = coeff_at_origin(substitute(2*u_x2 + sp.Rational(1,2)*f['delta1']*u_x1, f), (zq, zp))/t
    hand = coeff_at_origin(sp.expand(2*(f['delta1']/4 - 3*f['Phi1'])*f['v11']), (zq, zp))
    Kom_c = KA.closed_form(K['second']['K_om']).subs({k: pt['k'], q: pt['q'], p: pm, mu: pt['mu']}).subs(eta, t)
    KvT_c = KA.closed_form(K['second']['K_vT']).subs({k: pt['k'], q: pt['q'], p: pm, mu: pt['mu']}).subs(eta, t)
    N = lambda e: sp.N(e.subs(t, pt['eta']), 25)
    print("quad_engine - hand:", N(quad_engine - hand))
    print("Kelvin residual with hand quad:", N(qperp*(Kom_c + KvT_c) + hand))
    # alternative hypotheses: which weight of delta_T v and Phi v makes it vanish?
    a, b = sp.symbols('a b')
    alt = coeff_at_origin(sp.expand(2*(a*f['delta1']/4 + b*f['Phi1'])*f['v11']), (zq, zp))
    print("residual(a,b) =", sp.N(N(qperp*(Kom_c + KvT_c) + alt), 12))
