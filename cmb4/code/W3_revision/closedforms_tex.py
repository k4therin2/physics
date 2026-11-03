"""W3: extract printable closed forms for App. closedforms from C4's pickles and re-verify them by quadrature.
I_h = int_0^eta G_h K_h  (special-function skeleton + trigonometric part), Phi2 skeleton coefficient D."""
import pickle, sys, sympy as sp, mpmath as mp
sys.path.insert(0, '../C4_timeintegrals')
from trig_integrate import Cin
L = {'Cin': Cin}
k, p, q, eta = sp.symbols('k p q eta', positive=True)
d = pickle.load(open('../C4_timeintegrals/h2_compact.pkl', 'rb'))
s3 = sp.sqrt(3)
w = [k + (q + p)/s3, k + (q - p)/s3, k - (q - p)/s3, k - (q + p)/s3]
sf = 0
Fn = {'Cin': Cin, 'Si': sp.Si}
for kk, v in d['sf'].items():
    f, j = eval(kk)
    e = sp.parse_expr(v, local_dict=L).replace(sp.re, lambda a: a).replace(sp.im, lambda a: 0)
    e = sp.simplify(e)
    print('  sf', f, 'w%d' % (j + 1), ':', e)
    sf += e*Fn[f](w[j]*eta)
trig = sum(sp.parse_expr(a, local_dict=L)*sp.parse_expr(b, local_dict=L) for a, b in d['trig'].items())
Ih = sf + trig
# quadrature reference
def T(x): return -2*(mp.sin(x) - x*mp.cos(x))/x**3
def Tp(a, t): return mp.diff(lambda tt: T(a*tt/mp.sqrt(3)), t)
def Kh(qq, pp, t):
    Tq, Tpp = T(qq*t/mp.sqrt(3)), T(pp*t/mp.sqrt(3)); dq, dp = Tp(qq, t), Tp(pp, t)
    return 4*(3*Tq*Tpp + t*(dq*Tpp + Tq*dp) + t**2*dq*dp)
mp.mp.dps = 30
def Cin_n(z): return mp.euler + mp.log(abs(z)) - mp.ci(abs(z)) if z != 0 else mp.mpf(0)
fn = sp.lambdify((k, p, q, eta), Ih, modules=[{'Cin': Cin_n, 'Si': mp.si}, 'mpmath'])
for (kk, pp, qq, ee) in [(1.0, 0.7, 0.9, 3.3), (1.3, 0.4, 1.2, 7.1), (0.8, 1.5, 1.1, 2.0)]:
    ref = mp.quad(lambda t: (t/ee)*mp.sin(kk*(ee - t))/kk*Kh(qq, pp, t), [0, ee/2, ee])
    val = fn(mp.mpf(kk), mp.mpf(pp), mp.mpf(qq), mp.mpf(ee))
    print('I_h closed vs quadrature at', (kk, pp, qq, ee), ':', mp.nstr(val, 15), mp.nstr(ref, 15), 'rel', mp.nstr(abs(val/ref - 1), 3))
# LaTeX of the trigonometric part, grouped by trig monomial and power of eta
xk, xq, xp = sp.symbols('x_k x_q x_p')
out = []
for a, b in sorted(d['trig'].items()):
    m = sp.parse_expr(a, local_dict=L).subs({eta*k: sp.Symbol('kη')})
    c = sp.expand(sp.parse_expr(b, local_dict=L))
    m = sp.parse_expr(a, local_dict=L)
    mm = m.subs({sp.sqrt(3)*eta*p/3: sp.Symbol('x_p'), sp.sqrt(3)*eta*q/3: sp.Symbol('x_q'), eta*k: sp.Symbol('s_k')})
    terms = {}
    for t in sp.Add.make_args(c):
        n = sp.Poly(sp.expand(t*eta**30), eta).degree() - 30
        terms[n] = terms.get(n, 0) + t/eta**n
    body = r' \\ &\qquad '.join(r'+\Big[%s\Big]\eta^{%d}' % (sp.latex(sp.factor(sp.expand(v))), n) for n, v in sorted(terms.items()))
    out.append(r'&%s\times\Big\{%s\Big\}' % (sp.latex(mm), body))
open('Ih_trig.tex', 'w').write(r'\\' '\n'.join(out) + '\n')
print('trig monomials:', len(out))
