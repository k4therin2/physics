"""second_order_closed.py -- C4: closed-form second-order regular solutions in RD (Poisson gauge, a = eta),
per zeta_q zeta_p (p = k - q), built from C1's kernels (../C1_einstein) and the C4 term-by-term integration engine.

Conventions (SPEC.md, report_C1.md): X = eps X1 + eps^2 X2/2 + ...;  X2_k = int d^3q/(2pi)^3 K_X(q, k-q; eta) zeta_q zeta_{k-q};
k = k zhat, q = (q sqrt(1-mu^2), 0, q mu), p = |k - q| = sqrt(k^2 + q^2 - 2 k q mu).   x_k = k eta/sqrt3.
    Phi2   = int_0^eta G_Phi S_Phi + K_Phi2_0 ghat(k eta)          (closed form: Si, Cin + elementary)
    Psi2   = Phi2 - Q2                                             (Q2 = K_Q algebraic)
    delta2 = (eta^2/3) S00 - (2/3) eta^2 k^2 Phi2 - 2 eta Phi2' - 2 Psi2
    V2     = -(eta^2/4) [ 2 (Phi2' + Psi2/eta) + S0z/(i k) ]        (v2_i = d_i V2 + v2^T_i ;  Fourier d_i -> i k_i)
    omega2 : omega2_i = P_ij(k) q_j K_om  ->  k-frame x-component  omega2_x = q sqrt(1-mu^2) K_om   (K_om = i * real)
    vT2    : v2^T_i   = P_ij(k) q_j K_vT  ->  v2T_x = q sqrt(1-mu^2) K_vT                          (K_vT = i * real)
    deltaT2 = delta2/4 - (3/16) T_delta(q) T_delta(p)              (1 + delta_T = (1+delta)^{1/4})
    h2     : h2_ij = Lambda_ij,kl q_k q_l I_h,  I_h = int_0^eta G_h K_h = (4/k^2) I_RD^{KT}(u=p/k, v=q/k, x=k eta);
             in the k-frame only the + polarisation is sourced:  h2_xx = -h2_yy = q^2(1-mu^2)/2 I_h,  h2_xy = 0;
             helicity/plus amplitude (e^+_ij = (xx - yy)/sqrt2, e e = 1):  h2_plus = e^+_ij h2_ij = q^2 (1-mu^2)/sqrt2 I_h.
Special functions:  Si(z);  Cin(z) = EulerGamma + ln|z| - Ci(|z|) = int_0^z (1 - cos t)/t dt  (entire, even).

API:  every function F in (Phi2, dPhi2, Psi2, delta2, deltaT2, V2, omega2, vT2, h2_I, h2_plus) is called as
      F(k, q, mu, eta):  sympy expression if any argument is a sympy object, else a fast float/complex (numpy-vectorised)
      evaluation: ODE Taylor series for (k+q+p)eta/sqrt3 <= 4; closed form where it is well conditioned (all frequencies x eta >= 6,
      q/k and p/k in [1/5, 5]); otherwise composite Gauss-Legendre quadrature of the Green's-function integral (collinear edges,
      h2 resonance q+p ~ sqrt3 k, hierarchical triangles).  Accuracy vs 50-digit reference: see verify_module.log.
      F_expr(k, q, p, eta) sympy expressions in the invariants are also exported as module attributes (e.g. Phi2_expr).
"""
import os, sys, pickle, numpy as np, sympy as sp
from scipy.special import sici
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'C1_einstein'))
from trig_integrate import Cin
from kernels_api import load, closed_form, eta, k, q, p, mu

# ------------------------------------------------------------------ symbolic
_K = load(); _S = _K['second']
_L = {'Cin': Cin}
_P = sp.parse_expr(pickle.load(open(os.path.join(HERE, 'phi2_closed.pkl'), 'rb'))['P'], local_dict=_L)
_Ih = sp.parse_expr(pickle.load(open(os.path.join(HERE, 'h2_closed.pkl'), 'rb'))['I_h'], local_dict=_L)
_COLL = {kk: sp.parse_expr(vv, local_dict=_L) for kk, vv in pickle.load(open(os.path.join(HERE, 'collinear_closed.pkl'), 'rb')).items()}
def _real(e):
    e = sp.expand(e)
    return sp.expand(sp.re(e)) if e.has(sp.I) else e
_P = _real(_P); _Ih = _real(_Ih)
s3 = sp.sqrt(3)
x_of = lambda kk: kk * eta / s3
T_Phi = lambda kk: -2 * (sp.sin(x_of(kk)) - x_of(kk) * sp.cos(x_of(kk))) / x_of(kk)**3
T_delta = lambda kk: -2 * eta**2 * kk**2 * T_Phi(kk) / 3 - 2 * eta * sp.diff(T_Phi(kk), eta) - 2 * T_Phi(kk)
T_V = lambda kk: -eta * (eta * sp.diff(T_Phi(kk), eta) + T_Phi(kk)) / 2
ghat = _S['ghat']                                     # 3 u1(k eta) -> 1
K_Phi2_0 = _S['K_Phi2_0']                             # Phi2(k, eta->0) kernel (PUB 3.12b)
K_Q = closed_form(_S['K_Q']); K_S00 = closed_form(_S['K_S00']); K_S0z = closed_form(_S['K_S0z'])
K_om = closed_form(_S['K_om']); K_vT = closed_form(_S['K_vT']); K_h = closed_form(_S['K_h']); K_SPhi = closed_form(_S['K_SPhi'])

P_expr = _P                                           # particular integral int_0^eta G_Phi S_Phi
Phi2_expr = _P + K_Phi2_0 * ghat
dPhi2_expr = sp.diff(Phi2_expr, eta)
Psi2_expr = Phi2_expr - K_Q
delta2_expr = eta**2 / 3 * K_S00 - sp.Rational(2, 3) * eta**2 * k**2 * Phi2_expr - 2 * eta * dPhi2_expr - 2 * Psi2_expr
V2_expr = -eta**2 / 4 * (2 * (dPhi2_expr + Psi2_expr / eta) + K_S0z / (sp.I * k))
deltaT2_expr = delta2_expr / 4 - sp.Rational(3, 16) * T_delta(q) * T_delta(p)
omega2_expr = K_om          # times P_ij q_j
vT2_expr = K_vT             # times P_ij q_j
h2_I_expr = _Ih
Phi2_collinear = {kk: vv + K_Phi2_0 * ghat for kk, vv in _COLL.items() if kk.startswith('Phi2')}   # keys 'Phi2_q+p=k', 'Phi2_q=k+p', 'Phi2_p=k+q'
h2_I_resonant = _COLL['h2_q+p=sqrt3k']                 # q + p = sqrt3 k (contains log(eta))

def p_of(kk, qq, mm):
    return sp.sqrt(kk**2 + qq**2 - 2 * kk * qq * mm)

# ------------------------------------------------------------------ numerics
_EG = 0.5772156649015329
def cin_np(z):
    z = np.abs(np.asarray(z, dtype=float))
    out = np.empty_like(z)
    small = z < 0.5
    zs = z[small]
    acc = np.zeros_like(zs); term = np.ones_like(zs); z2 = zs * zs
    # Cin(z) = sum_{j>=1} (-1)^{j+1} z^{2j} / (2j (2j)!)
    fact = 1.0
    for j in range(1, 14):
        fact *= (2 * j - 1) * (2 * j)
        term = term * z2
        acc += (-1)**(j + 1) * term / (2 * j * fact)
    out[small] = acc
    zb = z[~small]
    out[~small] = _EG + np.log(zb) - sici(zb)[1]
    return out
def si_np(z):
    return sici(np.asarray(z, dtype=float))[0]
_NS = {'Cin': cin_np, 'Si': si_np, 'Ci': lambda z: sici(np.abs(z))[1]}
_lam_cache = {}
def _lam(name, expr):
    if name not in _lam_cache:
        _lam_cache[name] = sp.lambdify((k, q, p, eta), expr, modules=[_NS, 'numpy'], cse=True)
    return _lam_cache[name]

# --- numeric first-order transfer functions with series protection
def _tphi_np(kk, ee):
    x = kk * ee / np.sqrt(3.0); x = np.asarray(x, dtype=float)
    out = np.empty_like(x); small = np.abs(x) < 0.3
    xs = x[small]; x2 = xs * xs
    # T = sum_m tau_m x^{2m}, tau_m = -4 (m+1) (-1)^m / (2m+3)!
    acc = np.zeros_like(xs); f = 6.0; pw = np.ones_like(xs)
    for m in range(0, 12):
        if m > 0:
            f *= (2 * m + 2) * (2 * m + 3); pw = pw * x2
        acc += -4 * (m + 1) * (-1)**m / f * pw
    out[small] = acc
    xb = x[~small]
    out[~small] = -2 * (np.sin(xb) - xb * np.cos(xb)) / xb**3
    return out
def _dtphi_np(kk, ee):
    """d T_Phi / d eta"""
    x = np.asarray(kk * ee / np.sqrt(3.0), dtype=float)
    out = np.empty_like(x); small = np.abs(x) < 0.3
    xs = x[small]; x2 = xs * xs
    acc = np.zeros_like(xs); f = 6.0; pw = np.ones_like(xs)
    for m in range(0, 12):
        if m > 0:
            f *= (2 * m + 2) * (2 * m + 3); pw = pw * x2
            acc += -4 * (m + 1) * (-1)**m / f * 2 * m * pw / xs
    out[small] = acc
    xb = x[~small]
    out[~small] = (6 * np.sin(xb) - 6 * xb * np.cos(xb) - 2 * xb**2 * np.sin(xb)) / xb**4
    return out * kk / np.sqrt(3.0)

# --- algebraic kernels as functions of (T, T') to avoid trig cancellation at small eta
_Tq, _Tp, _dTq, _dTp = sp.symbols('Tq Tp dTq dTp')
def _abstract(expr):
    e = expr.doit()
    e = e.subs({sp.Derivative(sp.Function('T_Phi')(q, eta), eta): _dTq, sp.Derivative(sp.Function('T_Phi')(p, eta), eta): _dTp})
    e = e.subs({sp.Function('T_Phi')(q, eta): _Tq, sp.Function('T_Phi')(p, eta): _Tp})
    return sp.expand(e)
_alg = {nm: sp.lambdify((k, q, p, eta, _Tq, _Tp, _dTq, _dTp), _abstract(_S[nm]), 'numpy')
        for nm in ('K_Q', 'K_S00', 'K_S0z', 'K_om', 'K_vT', 'K_h', 'K_SPhi')}
def _alg_eval(nm, kk, qq, pp, ee):
    return _alg[nm](kk, qq, pp, ee, _tphi_np(qq, ee), _tphi_np(pp, ee), _dtphi_np(qq, ee), _dtphi_np(pp, ee))

# --- ODE Taylor series (small eta) for the particular integrals.  Coefficients from the sources' Taylor series.
_poly_SPhi = sp.Poly(_abstract(_S['K_SPhi']), _Tq, _Tp, _dTq, _dTp)
_mono_SPhi = [(m, [sp.lambdify((k, q, p), sp.expand(c).coeff(eta, j), 'numpy') for j in range(3)]) for m, c in _poly_SPhi.terms()]
_NSER = 44
def _tseries(qq, N):
    """Taylor coefficients (in eta) of T_Phi(q, eta) and T_Phi'(q, eta), arrays of length N"""
    qq = np.asarray(qq, dtype=float)
    T = np.zeros((N,) + qq.shape); dT = np.zeros((N,) + qq.shape)
    f = 6.0; pw = np.ones_like(qq); q2 = qq * qq / 3.0
    for m in range(0, N // 2):
        if m > 0:
            f *= (2 * m + 2) * (2 * m + 3); pw = pw * q2
        T[2 * m] = -4 * (m + 1) * (-1)**m / f * pw
        if 2 * m - 1 >= 0:
            dT[2 * m - 1] = 2 * m * T[2 * m]
    return T, dT
def _mul(a, b):
    N = a.shape[0]; out = np.zeros_like(a)
    for i in range(N):
        out[i] = np.sum(a[:i + 1] * b[i::-1], axis=0)
    return out
def _series_source_phi(kk, qq, pp, N):
    Tq, dTq = _tseries(qq, N); Tp, dTp = _tseries(pp, N)
    s = np.zeros_like(Tq)
    base = {_Tq: Tq, _Tp: Tp, _dTq: dTq, _dTp: dTp}
    for mono, cfs in _mono_SPhi:
        prod = None
        for sym, e in zip((_Tq, _Tp, _dTq, _dTp), mono):
            for _ in range(e):
                prod = base[sym] if prod is None else _mul(prod, base[sym])
        for j in range(3):
            c = np.asarray(cfs[j](kk, qq, pp), dtype=float) * np.ones_like(kk)
            if np.any(c != 0):
                s[j:] += c * prod[:N - j]
    return s
def _series_phi(kk, qq, pp, ee, N=_NSER):
    """particular integral P and dP/deta by the Frobenius recursion a_n n(n+3) = s_{n-2} - (k^2/3) a_{n-2}"""
    kk, qq, pp, ee = np.broadcast_arrays(*[np.asarray(a, dtype=float) for a in (kk, qq, pp, ee)])
    s = _series_source_phi(kk, qq, pp, N)
    a = np.zeros_like(s)
    for n in range(2, N):
        a[n] = (s[n - 2] - kk**2 / 3 * a[n - 2]) / (n * (n + 3))
    P = np.zeros_like(ee); dP = np.zeros_like(ee)
    for n in range(N - 1, 1, -1):
        P = P * ee + a[n]; dP = dP * ee + n * a[n]
    return P * ee**2, dP * ee
def _series_h(kk, qq, pp, ee, N=_NSER):
    kk, qq, pp, ee = np.broadcast_arrays(*[np.asarray(a, dtype=float) for a in (kk, qq, pp, ee)])
    Tq, dTq = _tseries(qq, N); Tp, dTp = _tseries(pp, N)
    Aq = Tq.copy(); Aq[1:] += dTq[:-1]            # eta T' + T
    Ap = Tp.copy(); Ap[1:] += dTp[:-1]
    s = 4 * (2 * _mul(Tq, Tp) + _mul(Aq, Ap))
    a = np.zeros_like(s)
    for n in range(2, N):
        a[n] = (s[n - 2] - kk**2 * a[n - 2]) / (n * (n + 1))
    I = np.zeros_like(ee)
    for n in range(N - 1, 1, -1):
        I = I * ee + a[n]
    return I * ee**2

Z_SWITCH = 4.0      # ODE series used for (k+q+p) eta/sqrt3 <= Z_SWITCH (calibrated in calibrate_numerics.py)
W_SAFE = 6.0        # closed form used only if every frequency x eta >= W_SAFE ...
RATIO_SAFE = 5.0    # ... and 1/RATIO_SAFE <= q/k, p/k <= RATIO_SAFE; otherwise double-precision Gauss-Legendre quadrature of the
                    # Green's-function integral (verify_module.py: the closed form loses digits by cancellation of the 1/w and
                    # 1/(p^3 q^3) terms near the collinear edges / h2 resonance and for hierarchical triangles).
_GL_X, _GL_W = np.polynomial.legendre.leggauss(40)
def _quad_panels(z):
    return np.maximum(2, np.ceil((2 * z + 1) / 20.0)).astype(int)
def _gw_nodes(ee, M):
    """composite Gauss-Legendre nodes/weights on [0, ee] with M panels (ee: array, same M for all)"""
    edges = np.linspace(0, 1, M + 1)
    t = ((edges[:-1, None] + edges[1:, None]) / 2 + (edges[1:, None] - edges[:-1, None]) / 2 * _GL_X[None, :]).ravel()
    w = ((edges[1:, None] - edges[:-1, None]) / 2 * _GL_W[None, :]).ravel()
    return ee[:, None] * t[None, :], ee[:, None] * w[None, :]
def _u1(kk, ee): return -0.5 * _tphi_np(kk, ee)
def _du1(kk, ee): return -0.5 * _dtphi_np(kk, ee)
def _u2(kk, ee):
    x = kk * ee / np.sqrt(3.0); return (np.cos(x) + x * np.sin(x)) / x**3
def _du2(kk, ee):
    x = kk * ee / np.sqrt(3.0); return (x**2 * np.cos(x) - 3 * np.cos(x) - 3 * x * np.sin(x)) / x**4 * kk / np.sqrt(3.0)
def _quad_phi(kk, qq, pp, ee):
    """P and dP/deta by quadrature:  P = u2(eta) int (u1/W) S - u1(eta) int (u2/W) S"""
    P = np.empty_like(ee); dP = np.empty_like(ee)
    z = (kk + qq + pp) * ee / np.sqrt(3.0); Ms = _quad_panels(z)
    for M in np.unique(Ms):
        sel = Ms == M
        k_, q_, p_, e_ = kk[sel], qq[sel], pp[sel], ee[sel]
        t, w = _gw_nodes(e_, M)
        K_, Q_, P_ = [np.repeat(a[:, None], t.shape[1], axis=1) for a in (k_, q_, p_)]
        S = _alg_eval('K_SPhi', K_, Q_, P_, t)
        xp = K_ * t / np.sqrt(3.0)
        u1W = -t * np.sin(xp) + K_ / np.sqrt(3.0) * t**2 * np.cos(xp)
        u2W = -t * np.cos(xp) - K_ / np.sqrt(3.0) * t**2 * np.sin(xp)
        A = np.sum(w * u1W * S, axis=1); B = np.sum(w * u2W * S, axis=1)
        P[sel] = _u2(k_, e_) * A - _u1(k_, e_) * B
        dP[sel] = _du2(k_, e_) * A - _du1(k_, e_) * B
    return P, dP
def _quad_h(kk, qq, pp, ee):
    I = np.empty_like(ee)
    z = (kk + qq + pp) * ee / np.sqrt(3.0) + kk * ee; Ms = _quad_panels(z)
    for M in np.unique(Ms):
        sel = Ms == M
        k_, q_, p_, e_ = kk[sel], qq[sel], pp[sel], ee[sel]
        t, w = _gw_nodes(e_, M)
        K_, Q_, P_ = [np.repeat(a[:, None], t.shape[1], axis=1) for a in (k_, q_, p_)]
        Kh = _alg_eval('K_h', K_, Q_, P_, t)
        E_ = e_[:, None]
        I[sel] = np.sum(w * (t / E_) * np.sin(K_ * (E_ - t)) / K_ * Kh, axis=1)
    return I
def _safe(kk, qq, pp, ee, freqs):
    ok = (qq / kk >= 1 / RATIO_SAFE) & (qq / kk <= RATIO_SAFE) & (pp / kk >= 1 / RATIO_SAFE) & (pp / kk <= RATIO_SAFE)
    for f in freqs:
        ok &= np.abs(f) * ee >= W_SAFE
    return ok
def _P_num(kk, qq, pp, ee):
    kk, qq, pp, ee = [np.array(a, dtype=float) for a in np.broadcast_arrays(*[np.asarray(a, dtype=float) for a in (kk, qq, pp, ee)])]
    z = (kk + qq + pp) * ee / np.sqrt(3.0)
    P = np.empty_like(ee); dP = np.empty_like(ee)
    s3_ = np.sqrt(3.0)
    ser = z <= Z_SWITCH
    cf = ~ser & _safe(kk, qq, pp, ee, [(kk + qq - pp) / s3_, (kk - qq + pp) / s3_, (qq + pp - kk) / s3_])
    qd = ~ser & ~cf
    if np.any(cf):
        P[cf] = _lam('P', P_expr)(kk[cf], qq[cf], pp[cf], ee[cf])
        dP[cf] = _lam('dP', sp.diff(P_expr, eta))(kk[cf], qq[cf], pp[cf], ee[cf])
    if np.any(ser):
        P[ser], dP[ser] = _series_phi(kk[ser], qq[ser], pp[ser], ee[ser])
    if np.any(qd):
        P[qd], dP[qd] = _quad_phi(kk[qd], qq[qd], pp[qd], ee[qd])
    return P, dP
def _Ih_num(kk, qq, pp, ee):
    kk, qq, pp, ee = [np.array(a, dtype=float) for a in np.broadcast_arrays(*[np.asarray(a, dtype=float) for a in (kk, qq, pp, ee)])]
    z = (kk + qq + pp) * ee / np.sqrt(3.0)
    I = np.empty_like(ee); s3_ = np.sqrt(3.0)
    ser = z <= Z_SWITCH
    cf = ~ser & _safe(kk, qq, pp, ee, [kk - (qq + pp) / s3_, kk - (qq - pp) / s3_])
    qd = ~ser & ~cf
    if np.any(cf):
        I[cf] = _lam('Ih', h2_I_expr)(kk[cf], qq[cf], pp[cf], ee[cf])
    if np.any(ser):
        I[ser] = _series_h(kk[ser], qq[ser], pp[ser], ee[ser])
    if np.any(qd):
        I[qd] = _quad_h(kk[qd], qq[qd], pp[qd], ee[qd])
    return I
_ghat = sp.lambdify((k, eta), ghat, 'numpy'); _dghat = sp.lambdify((k, eta), sp.diff(ghat, eta), 'numpy')
_KPhi20 = sp.lambdify((k, q, p), K_Phi2_0, 'numpy')
def _ghat_np(kk, ee):
    x = kk * ee / np.sqrt(3.0); out = np.empty_like(np.asarray(x, dtype=float) * np.ones_like(ee))
    return 3 * _tphi_np(kk, ee) / (-2)     # ghat = 3 u1 = -(3/2) T_Phi
def _dghat_np(kk, ee):
    return -1.5 * _dtphi_np(kk, ee)

# ------------------------------------------------------------------ public API
def _is_sym(*args):
    return any(isinstance(a, sp.Basic) for a in args)
def _sub(expr, kk, qq, mm, ee):
    return expr.subs({k: kk, q: qq, p: p_of(kk, qq, mm), eta: ee}, simultaneous=True)
def _num_args(kk, qq, mm, ee):
    kk, qq, mm, ee = np.broadcast_arrays(*[np.asarray(a, dtype=float) for a in (kk, qq, mm, ee)])
    pp = np.sqrt(kk**2 + qq**2 - 2 * kk * qq * mm)
    return kk, qq, pp, mm, ee

def Phi2(kk, qq, mm, ee):
    if _is_sym(kk, qq, mm, ee): return _sub(Phi2_expr, kk, qq, mm, ee)
    kk, qq, pp, mm, ee = _num_args(kk, qq, mm, ee)
    P, dP = _P_num(kk, qq, pp, ee)
    return P + _KPhi20(kk, qq, pp) * _ghat_np(kk, ee)
def dPhi2(kk, qq, mm, ee):
    if _is_sym(kk, qq, mm, ee): return _sub(dPhi2_expr, kk, qq, mm, ee)
    kk, qq, pp, mm, ee = _num_args(kk, qq, mm, ee)
    P, dP = _P_num(kk, qq, pp, ee)
    return dP + _KPhi20(kk, qq, pp) * _dghat_np(kk, ee)
def Q2(kk, qq, mm, ee):
    if _is_sym(kk, qq, mm, ee): return _sub(K_Q, kk, qq, mm, ee)
    kk, qq, pp, mm, ee = _num_args(kk, qq, mm, ee)
    return _alg_eval('K_Q', kk, qq, pp, ee)
def Psi2(kk, qq, mm, ee):
    if _is_sym(kk, qq, mm, ee): return _sub(Psi2_expr, kk, qq, mm, ee)
    return Phi2(kk, qq, mm, ee) - Q2(kk, qq, mm, ee)
def delta2(kk, qq, mm, ee):
    if _is_sym(kk, qq, mm, ee): return _sub(delta2_expr, kk, qq, mm, ee)
    kk, qq, pp, mm, ee = _num_args(kk, qq, mm, ee)
    P, dP = _P_num(kk, qq, pp, ee); c = _KPhi20(kk, qq, pp)
    Ph = P + c * _ghat_np(kk, ee); dPh = dP + c * _dghat_np(kk, ee)
    Ps = Ph - _alg_eval('K_Q', kk, qq, pp, ee)
    return ee**2 / 3 * _alg_eval('K_S00', kk, qq, pp, ee) - 2.0 / 3 * ee**2 * kk**2 * Ph - 2 * ee * dPh - 2 * Ps
def _tdelta_np(kk, ee):
    return -2 * ee**2 * kk**2 * _tphi_np(kk, ee) / 3 - 2 * ee * _dtphi_np(kk, ee) - 2 * _tphi_np(kk, ee)
def deltaT2(kk, qq, mm, ee):
    if _is_sym(kk, qq, mm, ee): return _sub(deltaT2_expr, kk, qq, mm, ee)
    kk, qq, pp, mm, ee = _num_args(kk, qq, mm, ee)
    return delta2(kk, qq, mm, ee) / 4 - 3.0 / 16 * _tdelta_np(qq, ee) * _tdelta_np(pp, ee)
def V2(kk, qq, mm, ee):
    if _is_sym(kk, qq, mm, ee): return _sub(V2_expr, kk, qq, mm, ee)
    kk, qq, pp, mm, ee = _num_args(kk, qq, mm, ee)
    P, dP = _P_num(kk, qq, pp, ee); c = _KPhi20(kk, qq, pp)
    Ph = P + c * _ghat_np(kk, ee); dPh = dP + c * _dghat_np(kk, ee)
    Ps = Ph - _alg_eval('K_Q', kk, qq, pp, ee)
    S0z = _alg_eval('K_S0z', kk, qq, pp, ee)
    return np.real(-ee**2 / 4 * (2 * (dPh + Ps / ee) + S0z / (1j * kk)))
def omega2(kk, qq, mm, ee):
    """k-frame x-component  omega2_x = q sqrt(1-mu^2) K_om   (imaginary, as C1's K_om)"""
    if _is_sym(kk, qq, mm, ee): return qq * sp.sqrt(1 - mm**2) * _sub(K_om, kk, qq, mm, ee)
    kk, qq, pp, mm, ee = _num_args(kk, qq, mm, ee)
    return qq * np.sqrt(1 - mm**2) * _alg_eval('K_om', kk, qq, pp, ee)
def vT2(kk, qq, mm, ee):
    """k-frame x-component  v2T_x = q sqrt(1-mu^2) K_vT   (imaginary)"""
    if _is_sym(kk, qq, mm, ee): return qq * sp.sqrt(1 - mm**2) * _sub(K_vT, kk, qq, mm, ee)
    kk, qq, pp, mm, ee = _num_args(kk, qq, mm, ee)
    return qq * np.sqrt(1 - mm**2) * _alg_eval('K_vT', kk, qq, pp, ee)
def h2_I(kk, qq, mm, ee):
    """I_h = int_0^eta G_h K_h = (4/k^2) I_RD(u=p/k, v=q/k, k eta)"""
    if _is_sym(kk, qq, mm, ee): return _sub(h2_I_expr, kk, qq, mm, ee)
    kk, qq, pp, mm, ee = _num_args(kk, qq, mm, ee)
    return _Ih_num(kk, qq, pp, ee)
def h2_plus(kk, qq, mm, ee):
    """e^+_ij h2_ij = q^2 (1-mu^2)/sqrt2 * I_h  (h2_xx = -h2_yy = q^2(1-mu^2)/2 I_h, h2_xy = 0)"""
    if _is_sym(kk, qq, mm, ee): return qq**2 * (1 - mm**2) / sp.sqrt(2) * _sub(h2_I_expr, kk, qq, mm, ee)
    kk, qq, pp, mm, ee = _num_args(kk, qq, mm, ee)
    return qq**2 * (1 - mm**2) / np.sqrt(2.0) * _Ih_num(kk, qq, pp, ee)
h2 = h2_plus

__all__ = ['Phi2', 'dPhi2', 'Q2', 'Psi2', 'delta2', 'deltaT2', 'V2', 'omega2', 'vT2', 'h2_I', 'h2_plus', 'h2',
           'Phi2_expr', 'P_expr', 'dPhi2_expr', 'Psi2_expr', 'delta2_expr', 'deltaT2_expr', 'V2_expr', 'omega2_expr', 'vT2_expr',
           'h2_I_expr', 'Phi2_collinear', 'h2_I_resonant', 'K_Phi2_0', 'ghat', 'K_Q', 'K_S00', 'K_S0z', 'K_om', 'K_vT', 'K_h', 'K_SPhi',
           'T_Phi', 'T_delta', 'T_V', 'Cin', 'k', 'q', 'p', 'mu', 'eta', 'p_of', 'Z_SWITCH']

# ------------------------------------------------------------------ limits (lazy loaders)
def small_eta_series():
    """dict field -> {power n: coefficient(k,q,p)} for Phi2, Psi2, Q2, delta2, deltaT2 (to eta^4), V2 (eta^5), P, I_h (eta^6);
    from series_small.py (Frobenius + constraints, checked against the closed forms at 120 digits)."""
    d = pickle.load(open(os.path.join(HERE, 'series_small_eta.pkl'), 'rb'))
    return {nm: {n: sp.parse_expr(c) for n, c in cs.items()} for nm, cs in d.items()}
def large_eta_asymptotics():
    """dict: 'Phi2_asym' (particular integral P through eta^-3), 'Phi2_lead' (P at eta^-2), 'h2_asym' (I_h through eta^-2),
    'h2_lead' (I_h at eta^-1 = (4/k^2) KT(25)).  Add K_Phi2_0*ghat for the full Phi2 (ghat is elementary)."""
    d = pickle.load(open(os.path.join(HERE, 'asymptotics.pkl'), 'rb'))
    return {kk: sp.parse_expr(v, local_dict=_L) for kk, v in d.items()}
__all__ += ['small_eta_series', 'large_eta_asymptotics']
