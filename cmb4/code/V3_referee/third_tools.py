"""third_tools.py -- V3: discrete three-mode realisation helpers for the third-order checks."""
import sys, sympy as sp
from pw_tools import *
C1 = '/home/user/physics/cmb4/code/C1_einstein'
if C1 not in sys.path: sys.path.insert(0, C1)
import kernels_api as KA
K = KA.load()
from fourier import eta, k, q, p, mu
S2 = K['second']
KPhi2 = sp.Function('KPhi2'); Kh = sp.Function('Kh')
tadPhi2, tadQ2, taddelta2 = [sp.Function(n)(t) for n in ('tadPhi2', 'tadQ2', 'taddelta2')]
at = lambda expr, A, B, C: KA.closed_form(expr).subs({q: A, p: B, k: C}, simultaneous=True).subs(eta, t)

def pair_fields(zm, zn, qm, qn, am, an, amp=1):
    """second-order fields sourced by the pair (m,n): X2 discrete amplitude = 2 K per zeta_m zeta_n; returned scaled by amp
    (amp = 1/2 gives the eps^2 coefficient X2/2 used inside E; amp = 1 gives X2 itself)."""
    Kv = qm + qn; Km = sp.sqrt(sum(c**2 for c in Kv))
    ph = sp.exp(I * (Kv[0] * x + Kv[1] * y + Kv[2] * z))
    pre = 2 * zm * zn * ph * amp
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

VALS = {}
def numeric_map(expr, rng):
    reps = {}
    names = ('KPhi2', 'Kh', 'tadPhi2', 'tadQ2', 'taddelta2')
    for d in expr.atoms(sp.Derivative):
        f = d.expr
        if isinstance(f, sp.core.function.AppliedUndef) and f.func.__name__ in names:
            key = (f.func.__name__, tuple(round(float(a), 10) for a in f.args if a not in (t, eta)), sum(n for v_, n in d.variable_count))
            reps[d] = key
    for f in expr.atoms(sp.core.function.AppliedUndef):
        if f.func.__name__ in names:
            key = (f.func.__name__, tuple(round(float(a), 10) for a in f.args if a not in (t, eta)), 0)
            reps[f] = key
    return reps
def val(key, rng):
    if key not in VALS:
        VALS[key] = sp.Rational(rng.randint(-50, 50), 7)
    return VALS[key]
def ev(expr, rng, eta0, prec=30):
    reps = numeric_map(expr, rng)
    e = expr.subs({d_: val(kk_, rng) for d_, kk_ in reps.items()}, simultaneous=True)
    return sp.N(e.subs(t, eta0), prec)

s_sym = sp.Symbol('s_perp')
def c1_image_sum(expr, pt, sign_s):
    """C1 kernel at (mu, p=|k-q|) plus its q -> -q image (mu -> -mu, p -> |k+q|, sqrt(1-mu^2) -> sign_s * sqrt)."""
    km, qm = pt['k'], pt['q']
    pm_minus = sp.sqrt(km**2 + qm**2 - 2 * km * qm * pt['mu']); pm_plus = sp.sqrt(km**2 + qm**2 + 2 * km * qm * pt['mu'])
    e = expr.replace(lambda w: isinstance(w, sp.Pow) and w.base == 1 - mu**2, lambda w: s_sym**(2 * w.exp))
    e1 = KA.closed_form(e).subs({k: km, q: qm}).subs({p: pm_minus, mu: pt['mu']}).subs(s_sym, sp.sqrt(1 - pt['mu']**2)).subs(eta, t)
    e2 = KA.closed_form(e).subs({k: km, q: qm}).subs({p: pm_plus, mu: -pt['mu']}).subs(s_sym, sign_s * sp.sqrt(1 - pt['mu']**2)).subs(eta, t)
    return e1 + e2

def three_mode_fields(pt, z1, z2, z3, amp=sp.Rational(1, 2)):
    kv, qv, pv = vec3(pt['k'], pt['q'], pt['mu'])
    qm = pt['q']; km = pt['k']
    f = {}
    for nm, v_ in mode(z1, qv, qm).items(): f[nm] = f.get(nm, 0) + v_
    for nm, v_ in mode(z2, -qv, qm).items(): f[nm] = f.get(nm, 0) + v_
    for nm, v_ in mode(z3, kv, km).items(): f[nm] = f.get(nm, 0) + v_
    for (zm, zn, qa, qb, am, an) in [(z1, z2, qv, -qv, qm, qm), (z1, z3, qv, kv, qm, km), (z2, z3, -qv, kv, qm, km)]:
        for nm, v_ in pair_fields(zm, zn, qa, qb, am, an, amp).items(): f[nm] = f.get(nm, 0) + v_
    return f
