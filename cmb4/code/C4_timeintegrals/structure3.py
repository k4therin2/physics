"""structure3.py -- canonical real form of the Phi2 particular integral:
   P = u2(eta) A(eta) - u1(eta) B(eta),   A = A_sf + A_el,  B = B_sf + B_el,
   special-function skeleton in terms of the four POSITIVE (physical-domain) frequencies
   w1 = (k+q+p)/s3, w2 = (k+q-p)/s3, w3 = (k-q+p)/s3, w4 = (q+p-k)/s3 ;
   elementary parts collected as  sum_j [ C_j(eta) cos(w_j eta) + S_j(eta) sin(w_j eta) ] + const, C_j,S_j Laurent polys in eta,
   and tested for the sign-flip (master) structure."""
import sympy as sp, pickle, sys, time
sys.path.insert(0, '.'); sys.path.insert(0, '../C1_einstein')
from trig_integrate import *
from kernels_api import eta, k, q, p
t0 = time.time(); L = {'Cin': Cin}; s3 = sp.sqrt(3)
d = pickle.load(open('phi2_closed.pkl', 'rb'))
W = [(k + q + p) / s3, (k + q - p) / s3, (k - q + p) / s3, (q + p - k) / s3]
Wn = [sp.expand(w) for w in W]
def canon(w):
    w = sp.expand(w)
    for j, wj in enumerate(Wn):
        if sp.expand(w - wj) == 0: return j, 1
        if sp.expand(w + wj) == 0: return j, -1
    raise ValueError(w)
def real_form(e):
    e = sp.expand(e)
    sf = {}; el = {}; const = 0
    for term in sp.Add.make_args(e):
        fs = list(term.atoms(Cin) | term.atoms(Si))
        trig = [f for f in term.atoms(sp.exp) if f.args[0].has(eta)] + [f for f in term.atoms(sp.sin, sp.cos)]
        if fs:
            f = fs[0]; j, sg = canon(f.args[0] / eta); c = term / f
            if f.func == Si: c = c * sg
            sf[(f.func.__name__, j)] = sf.get((f.func.__name__, j), 0) + c
        else:
            el_term = term
            const += el_term
    return sf, const
out = {}
for nm in ('A', 'B'):
    e = sp.parse_expr(d[nm], local_dict=L)
    e = sp.expand(e.replace(lambda f: isinstance(f, sp.exp), lambda f: sp.cos(sp.expand(f.args[0] / sp.I)) + sp.I * sp.sin(sp.expand(f.args[0] / sp.I))))
    e = sp.expand((e + e.subs(sp.I, -sp.I)) / 2)      # real part (Cin, Si real; symbols positive)
    sf, rest = real_form(e)
    sf = {kk: sp.factor(sp.expand(v)) for kk, v in sf.items()}
    print(nm, "skeleton:", {kk: v for kk, v in sf.items() if v != 0}, "(%.0fs)" % (time.time() - t0))
    # elementary rest: collect trig(w_j eta) with Laurent coefficients
    rest = sp.expand(sp.expand_trig(rest)) if False else sp.expand(rest)
    coll = {}
    for term in sp.Add.make_args(rest):
        tr = [f for f in term.atoms(sp.sin, sp.cos)]
        assert len(tr) <= 1, term
        if tr:
            f = tr[0]; j, sg = canon(f.args[0] / eta)
            key = (f.func.__name__, j); c = term / f * (sg if f.func == sp.sin else 1)
        else:
            key = ('1', -1); c = term
        coll[key] = coll.get(key, 0) + c
    coll = {kk: sp.collect(sp.expand(v), eta, sp.factor) for kk, v in coll.items()}
    for kk, v in sorted(coll.items(), key=str):
        print("   %s: %s" % (kk, str(v)[:400]))
    out[nm] = (sf, coll)
pickle.dump({nm: ({str(a): sp.srepr(b) for a, b in sf.items()}, {str(a): sp.srepr(b) for a, b in coll.items()}) for nm, (sf, coll) in out.items()}, open('phi2_realform.pkl', 'wb'))
# reconstruct and compare numerically with the stored P
u1 = sp.Symbol('u1'); u2 = sp.Symbol('u2')
def rebuild(nm):
    sf, coll = out[nm]
    e = 0
    for (fn, j), c in sf.items(): e += c * (Cin if fn == 'Cin' else Si)(Wn[j] * eta)
    for (fn, j), c in coll.items(): e += c * ({'sin': sp.sin, 'cos': sp.cos}[fn](Wn[j] * eta) if fn != '1' else 1)
    return e
x = k * eta / s3
U1 = (sp.sin(x) - x * sp.cos(x)) / x**3; U2 = (sp.cos(x) + x * sp.sin(x)) / x**3
Pre = U2 * rebuild('A') - U1 * rebuild('B')
P = sp.parse_expr(d['P'], local_dict=L)
vals = {k: sp.Rational(11, 10), q: sp.Rational(7, 10), p: sp.Rational(9, 10), eta: sp.Rational(37, 10)}
print("rebuild check:", sp.N(Pre.subs(vals), 30), sp.N(P.subs(vals), 30))
print("done %.0fs" % (time.time() - t0))
