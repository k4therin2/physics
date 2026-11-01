"""check_third.py -- V3: third-order contracted kernels K3c (C1 third_order.pkl) vs the independent E^{(3)} realised on
three discrete modes zeta1 e^{iq.x}, zeta2 e^{-iq.x}, zeta3 e^{ik.x} with explicit second-order fields for the pairs.
Identity tested:  6 * [coefficient of zeta1 zeta2 zeta3 in [eps^3]E_ab]  ==  K3c(k,q,mu;p=|k-q|) + K3c(k,q,-mu;p=|k+q|)
(sum over the two images of q, with the tadpole pair (q,-q) represented by abstract tadPhi2/tadQ2/taddelta2(eta)).
Abstract KPhi2(a,b,c,eta), Kh(a,b,c,eta), tadpoles and their eta-derivatives are given random numerical values
(the same on both sides).  Verified second-order kernels (K_Q, K_S00, K_S0z, K_om, K_vT) are taken from C1."""
import sys, pickle, random, sympy as sp
from pw_tools import *
C1 = '/home/user/physics/cmb4/code/C1_einstein'; sys.path.insert(0, C1)
import kernels_api as KA; K = KA.load()
from fourier import eta, k, q, p, mu
E = E_mine()
S2 = K['second']; K3 = K['third']
KPhi2 = sp.Function('KPhi2'); Kh = sp.Function('Kh')
tadPhi2, tadQ2, taddelta2 = [sp.Function(n)(t) for n in ('tadPhi2', 'tadQ2', 'taddelta2')]
PREC = 30
z1, z2, z3 = sp.symbols('zeta1 zeta2 zeta3')
at = lambda expr, A, B, C: KA.closed_form(expr).subs({q: A, p: B, k: C}, simultaneous=True).subs(eta, t)

def pair_fields(zm, zn, qm, qn, am, an):
    """second-order fields (eps^2 coefficient = X2/2) sourced by the pair of modes (m,n); amplitudes 2*K per zeta_m zeta_n."""
    Kv = qm + qn; Km = sp.sqrt(sum(c**2 for c in Kv))
    ph = sp.exp(I * (Kv[0] * x + Kv[1] * y + Kv[2] * z))
    pre = 2 * zm * zn * ph / 2          # /2: eps^2 coefficient is X2/2
    out = {}
    if sp.simplify(Km) == 0:
        out['Phi2'] = pre * tadPhi2; out['Psi2'] = pre * (tadPhi2 - tadQ2); out['delta2'] = pre * taddelta2
        return out
    KQ = at(S2['K_Q'], am, an, Km); KS00 = at(S2['K_S00'], am, an, Km); KS0z = at(S2['K_S0z'], am, an, Km)
    KP = KPhi2(am, an, Km, t); KPs = KP - KQ
    Kde = t**2 / 3 * KS00 - sp.Rational(2, 3) * t**2 * Km**2 * KP - 2 * t * sp.diff(KP, t) - 2 * KPs
    KV = -(t**2 / 4) * (2 * (sp.diff(KP, t) + KPs / t) + KS0z / (I * Km))
    Kom = at(S2['K_om'], am, an, Km); KvT = at(S2['K_vT'], am, an, Km); Khh = Kh(am, an, Km, t)
    Pm = sp.eye(3) - sp.Matrix(Kv) * sp.Matrix(Kv).T / Km**2
    Lam = lambda i, j, a, b: Pm[i, a] * Pm[j, b] - sp.Rational(1, 2) * Pm[i, j] * Pm[a, b]
    qmv = sp.Matrix(qm)
    out['Phi2'] = pre * KP; out['Psi2'] = pre * KPs; out['delta2'] = pre * Kde
    for i in range(3):
        out['v%d2' % (i + 1)] = pre * (I * Kv[i] * KV + sum(Pm[i, j] * qmv[j] for j in range(3)) * KvT)
        out['om%d2' % (i + 1)] = pre * sum(Pm[i, j] * qmv[j] for j in range(3)) * Kom
    hn = {(0, 0): 'h112', (0, 1): 'h122', (0, 2): 'h132', (1, 1): 'h222', (1, 2): 'h232', (2, 2): 'h332'}
    for (i, j), nm in hn.items():
        out[nm] = pre * sum(Lam(i, j, a, b) * qmv[a] * qmv[b] for a in range(3) for b in range(3)) * Khh
    return out

def numeric_map(expr, rng):
    """map abstract function values / derivatives to random numbers, keyed by (name, numeric args, deriv order)."""
    reps = {}
    for d in expr.atoms(sp.Derivative):
        f = d.expr
        if isinstance(f, sp.core.function.AppliedUndef) and f.func.__name__ in ('KPhi2', 'Kh', 'tadPhi2', 'tadQ2', 'taddelta2'):
            key = (f.func.__name__, tuple(round(float(a), 10) for a in f.args if a not in (t, eta)), sum(n for v_, n in d.variable_count))
            reps[d] = key
    for f in expr.atoms(sp.core.function.AppliedUndef):
        if f.func.__name__ in ('KPhi2', 'Kh', 'tadPhi2', 'tadQ2', 'taddelta2'):
            key = (f.func.__name__, tuple(round(float(a), 10) for a in f.args if a not in (t, eta)), 0)
            reps[f] = key
    return reps
VALS = {}
def val(key, rng):
    if key not in VALS:
        VALS[key] = sp.Rational(rng.randint(-50, 50), 7)
    return VALS[key]

comps = {'S00': (0, 0), 'S0z': (0, 3), 'S0x': (0, 1), 'Sxz': (1, 3)}
worst = 0
for seed in range(3):
    rng = random.Random(100 + seed)
    pt = rand_point(seed)
    kv, qv, pv = vec3(pt['k'], pt['q'], pt['mu'])
    qm = pt['q']; km = pt['k']
    pm_minus = sp.sqrt(km**2 + qm**2 - 2 * km * qm * pt['mu']); pm_plus = sp.sqrt(km**2 + qm**2 + 2 * km * qm * pt['mu'])
    f = {}
    for nm, v_ in mode(z1, qv, qm).items(): f[nm] = f.get(nm, 0) + v_
    for nm, v_ in mode(z2, -qv, qm).items(): f[nm] = f.get(nm, 0) + v_
    for nm, v_ in mode(z3, kv, km).items(): f[nm] = f.get(nm, 0) + v_
    for (zm, zn, qa, qb, am, an) in [(z1, z2, qv, -qv, qm, qm), (z1, z3, qv, kv, qm, km), (z2, z3, -qv, kv, qm, km)]:
        for nm, v_ in pair_fields(zm, zn, qa, qb, am, an).items(): f[nm] = f.get(nm, 0) + v_
    # my side
    mine = {}
    for ab in [(0, 0), (0, 3), (0, 1), (1, 1), (2, 2), (3, 3), (1, 3)]:
        mine[ab] = 6 * coeff_at_origin(substitute(E[('E',) + ab + (3,)], f), (z1, z2, z3))
    def ev(expr):
        reps = numeric_map(expr, rng)
        e = expr.subs({d_: val(kk_, rng) for d_, kk_ in reps.items()})
        return sp.N(e.subs(t, pt['eta']), PREC)
    # C1 side: K3c at (mu, p-) + K3c at (-mu, p+); odd-in-q_x components (0x, xz) flip sign of sqrt(1-mu^2) in the image
    s_sym = sp.Symbol('s_perp')
    def c1(expr, sign_s):
        e = expr.replace(lambda w: isinstance(w, sp.Pow) and w.base == 1 - mu**2, lambda w: s_sym**(2 * w.exp))
        e1 = KA.closed_form(e).subs({k: km, q: qm}).subs({p: pm_minus, mu: pt['mu']}).subs(s_sym, sp.sqrt(1 - pt['mu']**2)).subs(eta, t)
        e2 = KA.closed_form(e).subs({k: km, q: qm}).subs({p: pm_plus, mu: -pt['mu']}).subs(s_sym, sign_s * sp.sqrt(1 - pt['mu']**2)).subs(eta, t)
        return e1 + e2
    print("point", seed, pt)
    checks = {'S00': (mine[(0, 0)], c1(K3['S00'], 1)), 'S0z': (mine[(0, 3)], c1(K3['S0z'], 1)),
              'Q (=S_TL/k^2)': ((mine[(3, 3)] - (mine[(1, 1)] + mine[(2, 2)]) / 2) / km**2, c1(K3['Q'], 1)),
              'Sxx-Syy': (mine[(1, 1)] - mine[(2, 2)], c1(K3['Sxxyy'], 1)),
              'S0x (image s->-s)': (mine[(0, 1)], c1(K3['S0x'], -1)), 'Sxz (image s->-s)': (mine[(1, 3)], c1(K3['Sxz'], -1))}
    for nm, (a, b) in checks.items():
        va, vb = ev(a), ev(b); d = abs(va - vb); worst = max(worst, d)
        print("   %-20s mine=%s  C1=%s  |diff|=%.2e" % (nm, sp.N(va, 12), sp.N(vb, 12), d))
print("WORST third-order contracted-kernel discrepancy:", worst)
