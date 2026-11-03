"""probe_asym2.py -- does (full - asym) scale like eta^-(nmax+1)?  mpmath-lambdified, 60 digits, large x."""
import sympy as sp, pickle, sys, mpmath as mp, time
sys.path.insert(0, '.'); sys.path.insert(0, '../C1_einstein')
from trig_integrate import *
from kernels_api import eta, k, q, p
L = {'Cin': Cin}
t0=time.time()
d = pickle.load(open('asymptotics.pkl', 'rb'))
P = sp.parse_expr(pickle.load(open('phi2_closed.pkl', 'rb'))['P'], local_dict=L)
Ih = sp.parse_expr(pickle.load(open('h2_closed.pkl', 'rb'))['I_h'], local_dict=L)
P_as = sp.parse_expr(d['Phi2_asym'], local_dict=L); Ih_as = sp.parse_expr(d['h2_asym'], local_dict=L)
mp.mp.dps = 80
def cin_mp(z): return mp.euler + mp.log(abs(z)) - mp.ci(abs(z)) if z != 0 else mp.mpf(0)
NS = {'Cin': cin_mp, 'Si': mp.si}
f = {nm: sp.lambdify((k, q, p, eta), e, modules=[NS, 'mpmath']) for nm, e in (('P', P), ('Pa', P_as), ('I', Ih), ('Ia', Ih_as))}
print("lambdified %.0fs" % (time.time()-t0))
for (kk, qq, pp) in [(mp.mpf(1), mp.mpf('0.8'), mp.mpf('0.75')), (mp.mpf('1.9'), mp.mpf('0.19'), mp.mpf('1.81248448269220')), (mp.mpf(1), mp.mpf(2), mp.mpf('1.5'))]:
    print("config k,q,p =", kk, qq, pp)
    for X in [100, 1000, 10000, 100000]:
        ee = mp.mpf(X) / kk
        a, b = f['P'](kk, qq, pp, ee), f['Pa'](kk, qq, pp, ee)
        c, d_ = f['I'](kk, qq, pp, ee), f['Ia'](kk, qq, pp, ee)
        print("  x=%6d  Phi2: full*x^2=%s (full-as)*x^4=%s | h2: full*x=%s (full-as)*x^3=%s" % (X, mp.nstr(mp.re(a)*X**2, 10), mp.nstr(mp.re(a-b)*X**4, 8), mp.nstr(mp.re(c)*X, 10), mp.nstr(mp.re(c-d_)*X**3, 8)))
