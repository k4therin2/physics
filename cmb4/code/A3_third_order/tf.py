"""tf.py (A3) -- numerical transfer functions and C1 kernels as numpy callables.
Conventions: SPEC + report_C1.  a = eta, x = k eta/sqrt3.
First order per zeta_k:  Phi1 = Psi1 = T_Phi,  delta_T1 = T_dT (= T_delta/4),  v1^i = i k_i T_V,  T_V = eta gV(x).
Second order (C1, per zeta_q zeta_p with K = q + p; invariant symbols (k,q,p) = (|K|,|q1|,|q2|)): KER2[name](K,q,p,eta).
Third order contracted (C1 third_order.pkl): K3c_X(k,q,mu;eta) with abstract KPhi2(q,k,p,eta), Kh(q,k,p,eta) and tadpoles
tadPhi2(eta), tadQ2(eta), taddelta2(eta) -> lambdified with those objects (and their eta-derivatives) as explicit arguments.
"""
import sys, os, pickle
import numpy as np
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
C1 = os.path.join(HERE, '..', 'C1_einstein')
sys.path.insert(0, C1)
from kernels_api import load
from fourier import eta as ETA, k as KS, q as QS, p as PS, mu as MUS, TF

SQ3 = np.sqrt(3.0)
_x = sp.Symbol('x', positive=True)
EXPR = {'Phi': -2 * (sp.sin(_x) - _x * sp.cos(_x)) / _x**3,
        'dT': (-sp.cos(_x) + 2 * sp.sin(_x) / _x + 2 * sp.cos(_x) / _x**2 - 2 * sp.sin(_x) / _x**3),
        'gV': (_x**2 * sp.sin(_x) + 2 * _x * sp.cos(_x) - 2 * sp.sin(_x)) / _x**3}
NSER, NDER = 40, 6
_ser, _fun = {}, {}
for nm, f in EXPR.items():
    s = sp.series(f, _x, 0, NSER).removeO()
    _ser[nm] = [float(s.coeff(_x, n)) for n in range(NSER)]
    _fun[nm] = [sp.lambdify(_x, sp.diff(f, _x, d), 'numpy') for d in range(NDER + 1)]

def _poly_d(c, x, d):
    cc = list(c)
    for _ in range(d):
        cc = [n * cc[n] for n in range(1, len(cc))]
    out = np.zeros_like(x)
    for a in reversed(cc):
        out = out * x + a
    return out

def fx(nm, x, d=0):
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    small = x < 1.5
    if np.any(small): out[small] = _poly_d(_ser[nm], x[small], d)
    if np.any(~small): out[~small] = _fun[nm][d](x[~small])
    return out

def tphi(k, eta, d=0):
    k = np.asarray(k, float); eta = np.asarray(eta, float)
    return fx('Phi', k * eta / SQ3 + 0 * (k + eta), d) * (k / SQ3)**d

def tdt(k, eta, d=0):
    k = np.asarray(k, float); eta = np.asarray(eta, float)
    return fx('dT', k * eta / SQ3 + 0 * (k + eta), d) * (k / SQ3)**d

def tdelta(k, eta, d=0):
    return 4 * tdt(k, eta, d)

def tv(k, eta, d=0):
    """T_V = eta gV(x) and eta-derivatives: d^d (eta g) = eta g^(d) + d g^(d-1)."""
    k = np.asarray(k, float); eta = np.asarray(eta, float)
    x = k * eta / SQ3 + 0 * (k + eta)
    g = lambda dd: fx('gV', x, dd) * (k / SQ3)**dd
    return eta * g(d) + (d * g(d - 1) if d > 0 else 0)

def first(name, k, eta, d=0):
    if name in ('Phi', 'Psi'): return tphi(k, eta, d)
    if name == 'dT': return tdt(k, eta, d)
    if name == 'delta': return tdelta(k, eta, d)
    if name == 'V': return tv(k, eta, d)
    raise ValueError(name)

# ---------------------------------------------------------------- second order (C1)
KC1 = load()
_K2 = KC1['second']
_T, _dT, _ddT = sp.Function('TT'), sp.Function('dTT'), sp.Function('ddTT')
def _prep(expr, moms=(QS, PS, KS)):
    expr = expr.doit()
    for kk in moms:
        expr = expr.subs(sp.Derivative(TF['Phi'](kk, ETA), (ETA, 2)), _ddT(kk, ETA))
        expr = expr.subs(sp.Derivative(TF['Phi'](kk, ETA), ETA), _dT(kk, ETA))
        expr = expr.subs(TF['Phi'](kk, ETA), _T(kk, ETA))
    return expr
_mods = [{'TT': lambda kk, e: tphi(kk, e, 0), 'dTT': lambda kk, e: tphi(kk, e, 1), 'ddTT': lambda kk, e: tphi(kk, e, 2)}, 'numpy']
N2 = ['K_Q', 'K_S00', 'K_S0z', 'K_SPhi', 'K_om', 'K_vT', 'K_h', 'K_Phi2_0']
K2SYM = {nm: _prep(_K2[nm]) for nm in N2}
for nm in N2: assert not K2SYM[nm].has(sp.Derivative), nm
KER2 = {nm: sp.lambdify((KS, QS, PS, ETA), K2SYM[nm], modules=_mods) for nm in N2}
# eta-derivative of K_h (for h2''')
KER2['dK_h'] = sp.lambdify((KS, QS, PS, ETA), _prep(sp.diff(_K2['K_h'], ETA)), modules=_mods)
KER2['dK_SPhi'] = sp.lambdify((KS, QS, PS, ETA), _prep(sp.diff(_K2['K_SPhi'], ETA)), modules=_mods)
_etap = sp.Symbol('etap', positive=True)
G_Phi = sp.lambdify((KS, ETA, _etap), _K2['G_Phi'], 'numpy')
dG_Phi = sp.lambdify((KS, ETA, _etap), sp.diff(_K2['G_Phi'], ETA), 'numpy')
G_h = sp.lambdify((KS, ETA, _etap), _K2['G_h'], 'numpy')
dG_h = sp.lambdify((KS, ETA, _etap), sp.diff(_K2['G_h'], ETA), 'numpy')
def ghat(K, eta, d=0):
    """regular homogeneous solution of the Phi operator normalised to 1 at eta=0: ghat = -(3/2) T_Phi."""
    return -1.5 * tphi(K, eta, d)

# ---------------------------------------------------------------- tadpoles (C1 third_order_ic, angle-averaged k->0)
_tads = KC1['third_ic']['tads']
_qs, _es = sp.Symbol('q', positive=True), sp.Symbol('eta', positive=True)
TADSYM = {nm: _tads[nm] for nm in ('tadPhi2', 'tadQ2', 'taddelta2', 'tadS00', 'tadSPhi')}
def _series_or_closed(expr, nder=3, nser=24):
    """numpy callable f(q, eta, d) with small-x Taylor branch (x = q eta) to avoid cancellations."""
    zs = sp.Symbol('zz', positive=True)
    e = expr.subs(_es, zs / _qs)            # function of z = q eta only (dimensionless) -- check
    e = sp.simplify(e)
    assert _qs not in e.free_symbols, expr
    fl = [sp.lambdify(zs, sp.diff(e, zs, d), 'numpy') for d in range(nder + 1)]
    s = sp.series(e, zs, 0, nser).removeO()
    cs = [float(s.coeff(zs, n)) for n in range(nser)]
    def f(q, eta, d=0):
        q = np.asarray(q, float); eta = np.asarray(eta, float)
        z = q * eta + 0 * (q + eta)
        out = np.empty_like(z)
        sm = z < 2.0
        if np.any(sm): out[sm] = _poly_d(cs, z[sm], d)
        if np.any(~sm): out[~sm] = fl[d](z[~sm])
        return out * q**d
    return f
TAD = {nm: _series_or_closed(TADSYM[nm]) for nm in ('tadPhi2', 'tadQ2', 'taddelta2')}

# ---------------------------------------------------------------- third order contracted (C1)
_K3 = KC1['third']
_KPhi2, _Kh = sp.Function('KPhi2'), sp.Function('Kh')
_tP, _tQ, _td = sp.Function('tadPhi2')(ETA), sp.Function('tadQ2')(ETA), sp.Function('taddelta2')(ETA)
ARGS3 = sp.symbols('Tk dTk ddTk Tq dTq ddTq P0 P1 P2 H0 H1 H2 H3 tP0 tP1 tP2 tQ0 tQ1 td0')
def _prep3(expr):
    expr = expr.doit()
    Tk, dTk, ddTk, Tq, dTq, ddTq, P0, P1, P2, H0, H1, H2, H3, tP0, tP1, tP2, tQ0, tQ1, td0 = ARGS3
    rep = {}
    for kk, (a0, a1, a2) in ((KS, (Tk, dTk, ddTk)), (QS, (Tq, dTq, ddTq))):
        rep[sp.Derivative(TF['Phi'](kk, ETA), (ETA, 2))] = a2
        rep[sp.Derivative(TF['Phi'](kk, ETA), ETA)] = a1
        rep[TF['Phi'](kk, ETA)] = a0
    F = _KPhi2(QS, KS, PS, ETA); H = _Kh(QS, KS, PS, ETA)
    rep[sp.Derivative(F, (ETA, 2))] = P2; rep[sp.Derivative(F, ETA)] = P1; rep[F] = P0
    rep[sp.Derivative(H, (ETA, 3))] = H3; rep[sp.Derivative(H, (ETA, 2))] = H2; rep[sp.Derivative(H, ETA)] = H1; rep[H] = H0
    rep[sp.Derivative(_tP, (ETA, 2))] = tP2; rep[sp.Derivative(_tP, ETA)] = tP1; rep[_tP] = tP0
    rep[sp.Derivative(_tQ, ETA)] = tQ1; rep[_tQ] = tQ0; rep[_td] = td0
    # order matters: replace highest derivatives first
    keys = sorted(rep.keys(), key=lambda e: -(sum(c for _, c in e.variable_count) if isinstance(e, sp.Derivative) else 0))
    for kk in keys:
        expr = expr.xreplace({kk: rep[kk]})
    assert not expr.has(sp.Derivative) and not expr.atoms(sp.core.function.AppliedUndef), expr.atoms(sp.core.function.AppliedUndef)
    return expr
N3 = ['Q', 'S00', 'S0z', 'SPhi']
K3SYM = {nm: _prep3(_K3[nm]) for nm in N3}
KER3 = {nm: sp.lambdify((KS, QS, PS, MUS, ETA) + tuple(ARGS3), K3SYM[nm], 'numpy') for nm in N3}
_ic = KC1['third_ic']
C3SYM = _ic['C3']
C3 = sp.lambdify((KS, QS, MUS), C3SYM, 'numpy')

if __name__ == '__main__':
    import mpmath as mp
    for kk, ee in ((0.3, 0.01), (1.0, 2.0), (2.0, 7.5)):
        x = mp.mpf(kk) * ee / mp.sqrt(3)
        print(kk, ee, float(tphi(kk, ee)), float(-2 * (mp.sin(x) - x * mp.cos(x)) / x**3), float(tv(kk, ee, 2)),
              float(mp.diff(lambda e: e * ((e*kk/mp.sqrt(3))**2*mp.sin(e*kk/mp.sqrt(3)) + 2*(e*kk/mp.sqrt(3))*mp.cos(e*kk/mp.sqrt(3)) - 2*mp.sin(e*kk/mp.sqrt(3)))/(e*kk/mp.sqrt(3))**3, ee, 2)))
    print({nm: float(TAD[nm](1.0, 1e-4)) for nm in TAD}, ' expect tadPhi2(0)=-8/9, tadQ2(0)=-16/9, taddelta2(0)=16/9')
    print({nm: float(TAD[nm](1.0, 3.0)) for nm in TAD})
    print({nm: complex(KER3[nm](1.0, 1.0, 0.7, 0.755, 2.0, *np.linspace(0.1, 1.9, 19))) for nm in N3})
    print('C3(q->0)=', C3(1.0, 1e-6, 0.3), 'expect -32/9 =', -32 / 9)
