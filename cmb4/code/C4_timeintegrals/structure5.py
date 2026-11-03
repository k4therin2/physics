"""structure5.py -- compact real form of A(eta), B(eta) (P = u2 A - u1 B) built from the monomial tables termsA/B:
   A = D [Si(w1 eta) - Si(w2 eta) - Si(w3 eta) - Si(w4 eta)] + sum_j [C^A_j(eta) cos(w_j eta) + S^A_j(eta) sin(w_j eta)] + c^A_0,
   and likewise B (Cin skeleton); C, S Laurent polynomials in eta.  Rebuild check vs stored P.  -> phi2_compact.pkl"""
import sympy as sp, pickle, sys, time
sys.path.insert(0, '.'); sys.path.insert(0, '../C1_einstein')
from trig_integrate import *
from kernels_api import eta, k, q, p
t0 = time.time(); L = {'Cin': Cin}; s3 = sp.sqrt(3)
d = pickle.load(open('phi2_closed.pkl', 'rb'))
W = [sp.expand(w) for w in [(k + q + p) / s3, (k + q - p) / s3, (k - q + p) / s3, (q + p - k) / s3]]
def canon(w):
    w = sp.expand(w)
    for j, wj in enumerate(W):
        if sp.expand(w - wj) == 0: return j, 1
        if sp.expand(w + wj) == 0: return j, -1
    raise ValueError(w)
E = sp.symbols('E0:4'); Ec = sp.symbols('F0:4')       # e^{+i w_j eta}, e^{-i w_j eta}
out = {}
for nm in ('A', 'B'):
    terms = {eval(kk, {'sqrt': sp.sqrt, 'k': k, 'q': q, 'p': p}): sp.parse_expr(v) for kk, v in d['terms' + nm].items()}
    sf = {}; el = 0
    for (n, w), c in terms.items():
        fin, sing = integrate_term(n, w, eta)
        fin = sp.expand(c * fin)
        j, sg = canon(w)
        fin = fin.replace(lambda f: isinstance(f, sp.exp), lambda f: (E[j] if sg > 0 else Ec[j]) if sp.expand(f.args[0] - sp.I * w * eta) == 0 else f)
        for term in sp.Add.make_args(fin):
            fs = list(term.atoms(Cin) | term.atoms(Si))
            if fs:
                f = fs[0]; jj, s_ = canon(f.args[0] / eta); cc = term / f
                if f.func == Si: cc *= s_
                sf[(f.func.__name__, jj)] = sf.get((f.func.__name__, jj), 0) + cc
            elif term.has(sp.log):
                sf[('log', 0)] = sf.get(('log', 0), 0) + term
            else:
                el += term
    sf = {kk: sp.factor(sp.expand(v)) for kk, v in sf.items()}
    sf = {kk: v for kk, v in sf.items() if v != 0}
    el = sp.expand(el)
    C = {}; Sn = {}
    for j in range(4):
        cp = el.coeff(E[j]); cm = el.coeff(Ec[j])
        C[j] = sp.collect(sp.expand(cp + cm), eta, sp.factor); Sn[j] = sp.collect(sp.expand(sp.I * (cp - cm)), eta, sp.factor)
        assert not C[j].has(sp.I) and not Sn[j].has(sp.I), (nm, j)
    c0 = sp.factor(sp.expand(el.subs({**{e: 0 for e in E}, **{e: 0 for e in Ec}})))
    print(nm, "skeleton:", sf, " const:", c0, "(%.0fs)" % (time.time() - t0), flush=True)
    for j in range(4):
        print("   cos(w%d eta): %s\n   sin(w%d eta): %s" % (j + 1, str(C[j])[:300], j + 1, str(Sn[j])[:300]))
    out[nm] = (sf, C, Sn, c0)
x = k * eta / s3
U1 = (sp.sin(x) - x * sp.cos(x)) / x**3; U2 = (sp.cos(x) + x * sp.sin(x)) / x**3
Fn = {'Cin': Cin, 'Si': Si}
def rebuild(nm):
    sf, C, Sn, c0 = out[nm]
    return sum(v * Fn[f](W[j] * eta) for (f, j), v in sf.items() if f != 'log') + sum(C[j] * sp.cos(W[j] * eta) + Sn[j] * sp.sin(W[j] * eta) for j in range(4)) + c0
Pre = U2 * rebuild('A') - U1 * rebuild('B')
P = sp.parse_expr(d['P'], local_dict=L)
for vals in ({k: sp.Rational(11, 10), q: sp.Rational(7, 10), p: sp.Rational(9, 10), eta: sp.Rational(37, 10)},
             {k: sp.Rational(1), q: sp.Rational(5, 2), p: sp.Rational(2), eta: sp.Rational(23, 3)}):
    print("rebuild check:", sp.N(Pre.subs(vals), 30), sp.N(sp.re(P.subs(vals)), 30))
pickle.dump({nm: {'sf': {str(a): sp.srepr(b) for a, b in sf.items()}, 'C': {j: sp.srepr(v) for j, v in C.items()},
                  'S': {j: sp.srepr(v) for j, v in Sn.items()}, 'c0': sp.srepr(c0)} for nm, (sf, C, Sn, c0) in out.items()}, open('phi2_compact.pkl', 'wb'))
print("done %.0fs" % (time.time() - t0))
