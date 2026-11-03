"""s2_third.py (V6): 1/p^2 pole of C1's contracted third-order sources K3c_X(k, q=k, mu -> 1) (p = |k-q| -> 0).
KPhi2(q,k,p) -> a/p^2 + b + c p^2 (a = -tau/eta: the gauge pole from s1), Kh -> h0 + h1 p^2 (finite), mu = 1 - p^2/(2k^2), q = k.
Output: Laurent coefficients p^-4, p^-2 (and check that b, c, h0, h1 do not enter the singular part)."""
import sys, os, time, pickle
sys.path.insert(0, '../C1_einstein'); sys.path.insert(0, '../A3_third_order')
import sympy as sp
from kernels_api import load
from fourier import eta, k, q, p, mu, TF
K = load(); K3 = K['third']
# symbols for the abstract objects
Tk, dTk, ddTk = sp.symbols('Tk dTk ddTk')
P0, P1, P2, H0, H1, H2, H3 = sp.symbols('P0 P1 P2 H0 H1 H2 H3')
tP0, tP1, tP2, tQ0, tQ1, td0 = sp.symbols('tP0 tP1 tP2 tQ0 tQ1 td0')
a0, a1, a2, b0, b1, b2, c0, c1, c2 = sp.symbols('a0 a1 a2 b0 b1 b2 c0 c1 c2')
h00, h01, h02, h03, h10, h11, h12, h13 = sp.symbols('h00 h01 h02 h03 h10 h11 h12 h13')
KPhi2, Kh = sp.Function('KPhi2'), sp.Function('Kh')
tP, tQ, td = sp.Function('tadPhi2')(eta), sp.Function('tadQ2')(eta), sp.Function('taddelta2')(eta)
def prep(expr):
    expr = expr.doit()
    rep = {}
    for kk in (k, q):     # q = k: T_Phi(q) = T_Phi(k)
        rep[sp.Derivative(TF['Phi'](kk, eta), (eta, 2))] = ddTk
        rep[sp.Derivative(TF['Phi'](kk, eta), eta)] = dTk
        rep[TF['Phi'](kk, eta)] = Tk
    F = KPhi2(q, k, p, eta); H = Kh(q, k, p, eta)
    rep[sp.Derivative(F, (eta, 2))] = P2; rep[sp.Derivative(F, eta)] = P1; rep[F] = P0
    rep[sp.Derivative(H, (eta, 3))] = H3; rep[sp.Derivative(H, (eta, 2))] = H2; rep[sp.Derivative(H, eta)] = H1; rep[H] = H0
    rep[sp.Derivative(tP, (eta, 2))] = tP2; rep[sp.Derivative(tP, eta)] = tP1; rep[tP] = tP0
    rep[sp.Derivative(tQ, eta)] = tQ1; rep[tQ] = tQ0; rep[td] = td0
    keys = sorted(rep.keys(), key=lambda e: -(sum(c for _, c in e.variable_count) if isinstance(e, sp.Derivative) else 0))
    for kk in keys: expr = expr.xreplace({kk: rep[kk]})
    assert not expr.atoms(sp.core.function.AppliedUndef), expr.atoms(sp.core.function.AppliedUndef)
    return expr
sub = {P0: a0 / p**2 + b0 + c0 * p**2, P1: a1 / p**2 + b1 + c1 * p**2, P2: a2 / p**2 + b2 + c2 * p**2,
       H0: h00 + h10 * p**2, H1: h01 + h11 * p**2, H2: h02 + h12 * p**2, H3: h03 + h13 * p**2}
out = {}
for nm in ('Q', 'S00', 'S0z', 'SPhi'):
    t0 = time.time()
    e = prep(K3[nm]).subs(q, k)
    e = e.subs(mu, 1 - p**2 / (2 * k**2))
    e = e.xreplace(sub)
    coeffs = {}
    terms = sp.Add.make_args(sp.expand(e))
    # each term: rational in p times monomial in the symbols; collect by power of p
    acc = {}
    for t in terms:
        num, den = sp.fraction(sp.factor_terms(t))
        # den should be a monomial in p, k, eta (C1 kernels have no other denominators at q=k after mu-substitution?)
        acc.setdefault('n', []).append(t)
    e = sp.Add(*terms)
    e = sp.together(e)
    num, den = sp.fraction(e)
    num = sp.expand(num); den = sp.factor(den)
    print(nm, 'denominator:', den, ' (%.0fs)' % (time.time() - t0), flush=True)
    # den = monomial * polynomial in p possibly; series via polynomial division in p
    dpoly = sp.Poly(den, p); npoly = sp.Poly(num, p)
    lowd = min(m[0] for m in dpoly.monoms())
    d0 = sp.Poly(sp.expand(den / p**lowd), p)
    # n/d = p^-lowd * n / d0 with d0(0) != 0: series of 1/d0 to needed order
    N = 8
    d0c = [d0.coeff_monomial(p**i) for i in range(N)]
    assert d0c[0] != 0
    inv = [1 / d0c[0]]
    for i in range(1, N):
        inv.append(-sum(d0c[j] * inv[i - j] for j in range(1, i + 1)) / d0c[0])
    nc = [npoly.coeff_monomial(p**i) for i in range(N + lowd)]
    for n in range(-lowd, 1):
        idx = n + lowd
        c = sum(nc[j] * inv[idx - j] for j in range(idx + 1))
        c = sp.simplify(sp.expand(c))
        if c != 0: coeffs[n] = c
    out[nm] = coeffs
    print(nm, {n: (sorted(map(str, c.free_symbols))) for n, c in coeffs.items() if n < 0}, ' (%.0fs)' % (time.time() - t0), flush=True)
pickle.dump({nm: {n: sp.srepr(c) for n, c in d.items()} for nm, d in out.items()}, open('s2_poles.pkl', 'wb'))
