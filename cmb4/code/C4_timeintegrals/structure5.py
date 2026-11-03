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
def prim(n, w, Z):
    """finite part of int_eps^eta t^n e^{iwt} dt (same convention as trig_integrate.integrate_term), Z := e^{i w eta}; n <= 0"""
    if n == 0: return (Z - 1) / (sp.I * w)
    if n == -1: return -Cin(w * eta) + sp.I * Si(w * eta) + sp.log(eta)
    return eta**(n + 1) * Z / (n + 1) - sp.I * w / (n + 1) * prim(n + 1, w, Z) - (sp.I * w)**(-n - 1) / (sp.factorial(-n - 1) * (n + 1))
for nm in ('A', 'B'):
    terms = {eval(kk, {'sqrt': sp.sqrt, 'k': k, 'q': q, 'p': p}): sp.parse_expr(v) for kk, v in d['terms' + nm].items()}
    sf = {}; buckets = {}      # (trig key, eta power) -> list of coefficient terms
    for (n, w), c in terms.items():
        j, sg = canon(w)
        Zs = E[j] if sg > 0 else Ec[j]
        fin = prim(n, w, Zs)
        for term in sp.Add.make_args(sp.expand(c * fin)):
            fs = list(term.atoms(Cin) | term.atoms(Si))
            if fs:
                f = fs[0]; jj, s_ = canon(f.args[0] / eta); cc = term / f
                if f.func == Si: cc *= s_
                sf.setdefault((f.func.__name__, jj), []).append(cc)
                continue
            if term.has(sp.log):
                sf.setdefault(('log', 0), []).append(term); continue
            tk = [e for e in E + Ec if term.has(e)]
            assert len(tk) <= 1, term
            tk = tk[0] if tk else 1
            rest = term / tk
            pw = rest.as_coeff_exponent(eta)[1]
            buckets.setdefault((tk, pw), []).append(rest / eta**pw)
    print("  integrated (%.0fs)" % (time.time() - t0), flush=True)
    sf = {kk: sp.factor(sp.Add(*v)) for kk, v in sf.items()}
    sf = {kk: v for kk, v in sf.items() if v != 0}
    bc = {kk: sp.Add(*v) for kk, v in buckets.items()}
    C = {}; Sn = {}
    for j in range(4):
        pws = sorted({pw for (tk, pw) in bc if tk in (E[j], Ec[j])})
        C[j] = sp.Add(*[sp.factor(bc.get((E[j], pw), 0) + bc.get((Ec[j], pw), 0)) * eta**pw for pw in pws])
        Sn[j] = sp.Add(*[sp.factor(sp.I * (bc.get((E[j], pw), 0) - bc.get((Ec[j], pw), 0))) * eta**pw for pw in pws])
        assert not C[j].has(sp.I) and not Sn[j].has(sp.I), (nm, j)
    c0 = sp.factor(sp.Add(*[v * eta**pw for (tk, pw), v in bc.items() if tk == 1]))
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
