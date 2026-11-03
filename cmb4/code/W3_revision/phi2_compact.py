"""structure4.py -- compact real form of the closed-form particular integrals.
Phi2:  P = P_sf + P_el,
   P_sf = D(k,q,p) { u2(x) [Si(w1 eta) - Si(w2 eta) - Si(w3 eta) - Si(w4 eta)] - u1(x) Sum_j s_j Cin(w_j eta) }   (signs from data)
   P_el = sum over the 8 products  t_k t_q t_p  (t in {sin, cos} of x_k, x_q, x_p; x_a = a eta/sqrt3)  of  Laurent polynomials in eta.
h2: same with x_k -> k eta (sin/cos(k eta)) and x_q, x_p as above.
Output: phi2_compact.pkl / h2_compact.pkl and LaTeX snippets C4_appendix_*.tex."""
import sympy as sp, pickle, sys, time
sys.path.insert(0, '../C4_timeintegrals'); sys.path.insert(0, '../C1_einstein')
from trig_integrate import *
from kernels_api import eta, k, q, p
t0 = time.time(); L = {'Cin': Cin}; s3 = sp.sqrt(3)
X = {'k': k * eta / s3, 'q': q * eta / s3, 'p': p * eta / s3}
W = [(k + q + p) / s3, (k + q - p) / s3, (k - q + p) / s3, (q + p - k) / s3]
def canon(w, Wl):
    w = sp.expand(w)
    for j, wj in enumerate(Wl):
        if sp.expand(w - wj) == 0: return j, 1
        if sp.expand(w + wj) == 0: return j, -1
    raise ValueError(w)
def realpart(e):
    e = sp.expand(e)
    e = sp.expand(e.replace(lambda f: isinstance(f, sp.exp), lambda f: sp.cos(sp.expand(f.args[0] / sp.I)) + sp.I * sp.sin(sp.expand(f.args[0] / sp.I))))
    return sp.expand((e + e.subs(sp.I, -sp.I)) / 2)
def compact(e, Wl, xk):
    e = realpart(e)
    sf = {}; el = 0
    for term in sp.Add.make_args(e):
        fs = list(term.atoms(Cin) | term.atoms(Si))
        if fs:
            f = fs[0]; j, sg = canon(f.args[0] / eta, Wl); c = term / f
            if f.func == Si: c *= sg
            key = (f.func.__name__, j); sf[key] = sf.get(key, 0) + c
        else:
            el += term
    print("  split done (%.0fs)" % (time.time() - t0), flush=True)
    sf = {kk: sp.collect(sp.expand(v), [sp.sin(xk), sp.cos(xk)], sp.factor) for kk, v in sf.items()}
    # elementary part: expand trig of sums into products of trig(x_k), trig(x_q), trig(x_p)
    yk, yq, yp = sp.symbols('y_k y_q y_p')
    el = el.subs({k: yk * s3 / eta, q: yq * s3 / eta, p: yp * s3 / eta}, simultaneous=True) if False else el
    el = sp.expand(sp.expand_trig(el))
    print("  expand_trig done (%.0fs)" % (time.time() - t0), flush=True)
    trig = {}
    for term in sp.Add.make_args(el):
        tr = tuple(sorted([str(f) for f in sp.Mul.make_args(term) if isinstance(f, (sp.sin, sp.cos)) or (f.is_Pow and isinstance(f.base, (sp.sin, sp.cos)))]))
        c = sp.Mul(*[f for f in sp.Mul.make_args(term) if not (isinstance(f, (sp.sin, sp.cos)) or (f.is_Pow and isinstance(f.base, (sp.sin, sp.cos))))])
        mono = sp.Mul(*[f for f in sp.Mul.make_args(term) if (isinstance(f, (sp.sin, sp.cos)) or (f.is_Pow and isinstance(f.base, (sp.sin, sp.cos))))])
        trig[mono] = trig.get(mono, 0) + c
    trig = {m: sp.collect(sp.expand(c), eta, sp.factor) for m, c in trig.items()}
    trig = {m: c for m, c in trig.items() if c != 0}
    return sf, trig
res = {}
for nm, fn, key, Wl, xk in (
                            ('Phi2', '../C4_timeintegrals/phi2_closed.pkl', 'P', W, X['k']),):
    e = sp.parse_expr(pickle.load(open(fn, 'rb'))[key], local_dict=L)
    print(nm, flush=True)
    sf, trig = compact(e, Wl, xk)
    for kk, v in sf.items(): print("   SF %s: %s" % (kk, v))
    print("   %d elementary trig monomials:" % len(trig))
    for m, c in trig.items(): print("     %s : %d terms, eta powers %s" % (m, len(sp.Add.make_args(sp.expand(c))), sorted({sp.degree(t, eta) if t.has(eta) and sp.Poly(t * eta**20, eta).is_monomial else None for t in sp.Add.make_args(sp.expand(c * eta**20))}, key=str)[:3]))
    # rebuild check
    Fn = {'Cin': Cin, 'Si': Si}
    reb = sum(v * Fn[f](Wl[j] * eta) for (f, j), v in sf.items()) + sum(m * c for m, c in trig.items())
    vals = {k: sp.Rational(11, 10), q: sp.Rational(7, 10), p: sp.Rational(9, 10), eta: sp.Rational(37, 10)}
    print("   rebuild check: %s vs %s" % (sp.N(reb.subs(vals), 25), sp.N(sp.re(e.subs(vals)), 25)), flush=True)
    res[nm] = (sf, trig)
    pickle.dump({kk: sp.srepr(v) for kk, v in [('sf', {str(a): sp.srepr(b) for a, b in sf.items()}), ('trig', {sp.srepr(a): sp.srepr(b) for a, b in trig.items()})]} if False else
                {'sf': {str(a): sp.srepr(b) for a, b in sf.items()}, 'trig': {sp.srepr(a): sp.srepr(b) for a, b in trig.items()}, 'W': [sp.srepr(w) for w in Wl]},
                open(nm.lower() + '_compact.pkl', 'wb'))
print("done %.0fs" % (time.time() - t0))
