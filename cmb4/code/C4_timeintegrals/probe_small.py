import sys, numpy as np, sympy as sp, mpmath as mp
sys.path.insert(0, '.'); sys.path.insert(0, '../C1_einstein')
import second_order_closed as SO
from second_order_closed import *
from kernels_api import closed_form, load
S = load()['second']; etap = sp.Symbol('etap', positive=True)
integrand = sp.expand(S['G_Phi'] * closed_form(S['K_SPhi']).subs(eta, etap))
mp.mp.dps = 60
kk, qq, pp = sp.S(1), sp.Rational(4, 5), sp.Rational(3, 4)
for z in [sp.Rational(1, 2), 1, 2]:
    ee = z * sp.sqrt(3) / (kk + qq + pp)
    vals = {k: kk, q: qq, p: pp, eta: ee}
    f = sp.lambdify((etap,), integrand.subs(vals), 'mpmath')
    quad = mp.quad(f, mp.linspace(0, mp.mpf(sp.N(ee, 60)), 40), method='gauss-legendre', maxdegree=10)
    exact = sp.re(sp.N(P_expr.subs(vals), 60))
    ser = SO._series_phi(float(kk), float(qq), float(pp), float(sp.N(ee, 17)))[0]
    print("z=%s: quad=%s  exactN=%s  series(double)=%.16e" % (z, mp.nstr(quad, 20), sp.N(exact, 20), ser))
