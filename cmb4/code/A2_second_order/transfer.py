"""transfer.py -- numerical first-order transfer functions (C1, SPEC §5) and the C1 second-order kernels
as vectorised numpy callables.  Units: a = eta, x = k eta/sqrt3.

First order (per zeta_k):  Phi1 = Psi1 = T_Phi,  delta_T1 = T_dT = T_delta/4,  v1^i = i k_i T_V.
Small-x branch: Taylor series (generated with sympy, 16 terms) for x < 1, closed form otherwise.
Second order (C1 conventions, per zeta_q zeta_p, K = q + p; q,p,K magnitudes):
  Q2 = Phi2 - Psi2 = K_Q ;  Phi2'' + 4Phi2'/eta + K^2 Phi2/3 = K_SPhi ;  Phi2(eta->0) = K_Phi2_0 ;
  delta2 = (eta^2/3) K_S00 - (2/3) eta^2 K^2 Phi2 - 2 eta Phi2' - 2 Psi2 ;
  v2_z   = (eta/4)(-eta K_S0z - 2 i eta K Phi2' - 2 i K Psi2)   [component along K-hat]  => V2 = v2_z/(iK);
  om2_i  = P_ij(K) q_j K_om ;  v2T_i = P_ij(K) q_j K_vT ;  h2_ij = Lambda_ij,kl(K) q_k q_l int_0^eta G_h K_h.
"""
import sys, os, pickle
import numpy as np
import sympy as sp
C1 = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'C1_einstein')
sys.path.insert(0, C1)
from kernels_api import load, closed_form
from fourier import eta as ETA, k as KS, q as QS, p as PS, TF

SQ3 = np.sqrt(3.0)
_xs = sp.Symbol('x', positive=True)
_TPhi_x = -2 * (sp.sin(_xs) - _xs * sp.cos(_xs)) / _xs**3
_TdT_x = (-sp.cos(_xs) + 2 * sp.sin(_xs) / _xs + 2 * sp.cos(_xs) / _xs**2 - 2 * sp.sin(_xs) / _xs**3)
_gV_x = (_xs**2 * sp.sin(_xs) + 2 * _xs * sp.cos(_xs) - 2 * sp.sin(_xs)) / _xs**3     # T_V = eta * gV(x)
NSER = 34

def _series_coeffs(f):
    s = sp.series(f, _xs, 0, NSER).removeO()
    return [float(s.coeff(_xs, n)) for n in range(NSER)]

_cPhi, _cdT, _cgV = _series_coeffs(_TPhi_x), _series_coeffs(_TdT_x), _series_coeffs(_gV_x)
_fPhi = [sp.lambdify(_xs, sp.diff(_TPhi_x, _xs, d), 'numpy') for d in range(3)]
_fdT = [sp.lambdify(_xs, sp.diff(_TdT_x, _xs, d), 'numpy') for d in range(3)]
_fgV = [sp.lambdify(_xs, sp.diff(_gV_x, _xs, d), 'numpy') for d in range(3)]

def _poly_d(c, x, d):
    """d-th derivative of sum c_n x^n at x (Horner)."""
    cc = list(c)
    for _ in range(d):
        cc = [n * cc[n] for n in range(1, len(cc))]
    out = np.zeros_like(x)
    for a in reversed(cc):
        out = out * x + a
    return out

def _fx(c, fl, x, d):
    x = np.asarray(x, dtype=float)
    small = x < 1.0
    out = np.empty_like(x)
    if np.any(small):
        out[small] = _poly_d(c, x[small], d)
    if np.any(~small):
        out[~small] = fl[d](x[~small])
    return out

def tphi(k, eta, d=0):
    """d^d/deta^d T_Phi(k,eta)."""
    k = np.asarray(k, dtype=float); eta = np.asarray(eta, dtype=float)
    x = k * eta / SQ3 + 0 * (k + eta)
    return _fx(_cPhi, _fPhi, x, d) * (k / SQ3)**d

def tdt(k, eta, d=0):
    k = np.asarray(k, dtype=float); eta = np.asarray(eta, dtype=float)
    x = k * eta / SQ3 + 0 * (k + eta)
    return _fx(_cdT, _fdT, x, d) * (k / SQ3)**d

def tdelta(k, eta, d=0):
    return 4 * tdt(k, eta, d)

def tv(k, eta, d=0):
    """T_V = eta gV(x) and eta-derivatives (d <= 2)."""
    k = np.asarray(k, dtype=float); eta = np.asarray(eta, dtype=float)
    x = k * eta / SQ3 + 0 * (k + eta)
    g = lambda dd: _fx(_cgV, _fgV, x, dd) * (k / SQ3)**dd
    if d == 0: return eta * g(0)
    if d == 1: return g(0) + eta * g(1)
    if d == 2: return 2 * g(1) + eta * g(2)
    raise ValueError

# ---------------- second-order kernels (C1) ----------------
_K = load()['second']
_DTPhi = sp.Function('DT_Phi')
_D2TPhi = sp.Function('D2T_Phi')
def _prep(expr):
    expr = expr.doit()
    for kk in (QS, PS, KS):
        expr = expr.subs(sp.Derivative(TF['Phi'](kk, ETA), (ETA, 2)), _D2TPhi(kk, ETA))
        expr = expr.subs(sp.Derivative(TF['Phi'](kk, ETA), ETA), _DTPhi(kk, ETA))
    assert not expr.has(sp.Derivative), expr
    return expr
_mods = [{'T_Phi': lambda kk, e: tphi(kk, e, 0), 'DT_Phi': lambda kk, e: tphi(kk, e, 1), 'D2T_Phi': lambda kk, e: tphi(kk, e, 2)}, 'numpy']
_names = ['K_Q', 'K_S00', 'K_S0z', 'K_SPhi', 'K_om', 'K_vT', 'K_h', 'K_Phi2_0']
KER = {nm: sp.lambdify((KS, QS, PS, ETA), _prep(_K[nm]), modules=_mods) for nm in _names}
KER['K_Q_d'] = sp.lambdify((KS, QS, PS, ETA), _prep(sp.diff(_K['K_Q'].doit(), ETA)), modules=_mods)   # d/deta K_Q (exact)
KER_sym = {nm: _K[nm] for nm in _names}
# Green's functions and homogeneous mode
_etap = sp.Symbol('etap', positive=True)
GPhi_sym = _K['G_Phi']; Gh_sym = _K['G_h']; ghat_sym = _K['ghat']
G_Phi = sp.lambdify((KS, ETA, _etap), GPhi_sym, 'numpy')
dG_Phi = sp.lambdify((KS, ETA, _etap), sp.diff(GPhi_sym, ETA), 'numpy')
G_h = sp.lambdify((KS, ETA, _etap), Gh_sym, 'numpy')
dG_h = sp.lambdify((KS, ETA, _etap), sp.diff(Gh_sym, ETA), 'numpy')
def ghat(K, eta, d=0):
    """ghat = 3 u1 = -(3/2) T_Phi (regular homogeneous solution normalised to 1 at eta=0)."""
    return -1.5 * tphi(K, eta, d)
def G_h0(eta, etap):
    """K -> 0 limit of G_h: (eta'/eta)(eta - eta')."""
    return etap / eta * (eta - etap)
def dG_h0(eta, etap):
    return etap**2 / eta**2

if __name__ == '__main__':
    import mpmath as mp
    # spot checks against mpmath high precision
    for kk, ee in ((0.3, 0.01), (1.0, 2.0), (2.0, 7.5), (0.05, 3.0)):
        x = mp.mpf(kk) * ee / mp.sqrt(3)
        TP = -2 * (mp.sin(x) - x * mp.cos(x)) / x**3
        print(kk, ee, 'T_Phi', float(tphi(kk, ee)), float(TP),
              ' T_dT', float(tdt(kk, ee)), float(-mp.cos(x) + 2 * mp.sin(x) / x + 2 * mp.cos(x) / x**2 - 2 * mp.sin(x) / x**3),
              ' T_V', float(tv(kk, ee)), float(ee * (x**2 * mp.sin(x) + 2 * x * mp.cos(x) - 2 * mp.sin(x)) / x**3))
    # derivative check by finite differences
    h = 1e-5
    for f in (tphi, tdt, tv):
        for kk, ee in ((0.7, 0.5), (1.3, 4.0)):
            fd = (f(kk, ee + h) - f(kk, ee - h)) / (2 * h)
            print(f.__name__, 'd/deta FD', fd, 'analytic', f(kk, ee, 1))
    # kernels
    print({nm: complex(KER[nm](1.3, 0.7, 1.1, 2.3)) for nm in _names})
    # compare with sympy evaluation
    vals = {KS: 1.3, QS: 0.7, PS: 1.1, ETA: 2.3}
    print({nm: complex(closed_form(_K[nm]).subs(vals).evalf()) for nm in _names})
    print('G_Phi(eta,eta)=', G_Phi(1.3, 2.0, 2.0), ' dG_Phi(eta,eta)=', dG_Phi(1.3, 2.0, 2.0), ' ghat(0)=', ghat(1.3, 1e-6))
