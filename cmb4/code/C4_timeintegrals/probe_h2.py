import sympy as sp, mpmath as mp, pickle, sys
sys.path.insert(0, '.'); sys.path.insert(0, '../C1_einstein')
from trig_integrate import *
from kernels_api import load, closed_form, eta, k, q, p
d = pickle.load(open('h2_closed.pkl', 'rb'))
Ih = sp.parse_expr(d['I_h'], local_dict={'Cin': Cin}); IRD = sp.parse_expr(d['I_RD_KT'])
u, v, x = sp.symbols('u v x', positive=True)
K = load(); S = K['second']
etap = sp.Symbol('etap', positive=True)
Kh = closed_form(S['K_h']).subs(eta, etap); integrand = sp.expand(S['G_h'] * Kh)
# exact rational inputs, 50 digits
for (kk, qq, pp, ee) in [(1, sp.Rational(7, 10), sp.Rational(349, 1000), 6), (1, sp.Rational(5, 2), sp.Rational(3, 1), 3), (sp.Rational(3,2), sp.Rational(1,2), sp.Rational(1,1), 4)]:
    vals = {k: kk, q: qq, p: pp, eta: ee}
    a = sp.N(Ih.subs(vals), 50)
    b = sp.N((4 / k**2 * IRD.subs({u: p / k, v: q / k, x: k * eta})).subs(vals), 50)
    mp.mp.dps = 60
    f = sp.lambdify((etap,), integrand.subs(vals), 'mpmath')
    c = mp.quad(f, mp.linspace(0, ee, 40), method='gauss-legendre', maxdegree=10)
    c2 = mp.quad(f, mp.linspace(mp.mpf(1)/1000, ee, 40))   # tanh-sinh from a small cutoff
    # series for [0, 1e-3]: integrand ~ a1 eta' + ...; get from series of integrand
    print("u+v=%s  closed=%s\n           KT=%s\n           GL=%s\n       TS(cut)=%s" % (sp.N((qq+pp)/kk, 5), sp.N(a, 25), sp.N(b, 25), mp.nstr(c, 25), mp.nstr(c2, 25)))
