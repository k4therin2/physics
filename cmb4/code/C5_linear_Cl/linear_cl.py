"""
linear_cl.py -- agent C5: first-order temperature transfer function Delta_l(k), C_l^(11), closed-form ISW
integral, monochromatic collapse and the sky monopole <m1>, in SPEC conventions (radiation-only GR, a = eta).

CONVENTIONS (SPEC §5, §6, §9a; C2 report):
  Fourier  f(x) = int d^3k/(2pi)^3 e^{ik.x} f_k,  <zeta_k zeta_k'> = (2pi)^3 delta^3(k+k') P(k).
  n = direction the observer LOOKS at (SPEC §6: photon momentum phat_o ∝ uhat_o - n^a e_a, propagation direction -n),
  observer at x_o = 0, background ray x(eta) = + n chi(eta), chi(eta) = eta_o - eta, emission at eta_*, chi_* = eta_o - eta_*.
  (The C5 brief wrote x_e = x_o - n(eta_o - eta_*); that has dx/deta = +n and is inconsistent with the Doppler sign
   -n.v of SPEC §9a. The sign of the Doppler term relative to the geometry is PHYSICAL -- see rayleigh_check() and the
   report -- so we follow SPEC §6 as C2 does. With the other sign of n one gets a_lm -> (-1)^l a_lm, Delta_l unchanged,
   PROVIDED the Doppler term is simultaneously written +n.v.)
  First-order map (SPEC §9a):
     Theta1(n) = [dT1 + Psi1 - n.v1]_e - [dT1 + Psi1 - n.v1]_o + int_{eta_*}^{eta_o} (Phi1' + Psi1') d eta .
  First-order solution (C1, verified here):  x = k eta/sqrt3,
     Phi1 = Psi1 = zeta T_Phi(x),  T_Phi = -2 j_1(x)/x = -2(sin x - x cos x)/x^3
     dT1 = delta1/4 = zeta T_dT(x),  T_dT = -cos x + 2 j_0(x) + T_Phi(x)
     v1_i = d_i V1, V1 = zeta T_V,  k T_V = sqrt3 (sin x - 2 j_1(x))        (T_V = eta (x^2 sin x + 2x cos x - 2 sin x)/x^3)
     2 d_eta T_Phi = (4k/sqrt3) j_2(x)/x .
  RESULT:  a_lm = 4 pi i^l int d^3k/(2pi)^3 zeta_k Y*_lm(khat) Delta_l(k),
     Delta_l(k) = G_l(eta_*) - G_l(eta_o) + I_l(k),
     G_l(eta)   = A(x) j_l(k chi) - kT_V(x) j_l'(k chi),   A := T_dT + T_Phi,  x = k eta/sqrt3,  chi = eta_o - eta,
     I_l(k)     = int_{eta_*}^{eta_o} 2 d_eta T_Phi  j_l(k(eta_o - eta)) d eta
                = (4/sqrt3) int_{s_*}^{S} [j_2(s/sqrt3)/(s/sqrt3)] j_l(S - s) ds,   s_* = k eta_*, S = k eta_o.
     G_l(eta_o) = A(x_o) delta_{l0} - kT_V(x_o) delta_{l1}/3   (observer monopole: l = 0 only; observer Doppler: l = 1 only).
     C_l^(11) = (2/pi) int k^2 dk P(k) Delta_l(k)^2 ;  monochromatic P = 2 pi^2 sigma^2 delta(k-k0)/k0^2:  C_l = 4 pi sigma^2 Delta_l(k0)^2.
     Sky monopole <m1> = int d^3k/(2pi)^3 zeta_k M(k),  M(k) = Delta_0(k) (observer monopole INCLUDED).

Run `python3 linear_cl.py` for the full test suite (log: test_log.txt).
"""
import sys, os, functools, time
import numpy as np
import sympy as sp
import mpmath as mp
from scipy.special import spherical_jn, sph_harm_y, sici, eval_legendre
from scipy.integrate import quad
from numpy.polynomial.legendre import leggauss

HERE = os.path.dirname(os.path.abspath(__file__))
SQ3 = sp.sqrt(3)

# ----------------------------------------------------------------------------------------------------------------
# 1. symbols and first-order transfer functions
# ----------------------------------------------------------------------------------------------------------------
k, eta, eta_s, eta_o = sp.symbols('k eta eta_s eta_o', positive=True)     # eta_s = eta_*, eta_o observer time
x = sp.Symbol('x', positive=True)                                         # x = k eta / sqrt3
s, S, sst, u = sp.symbols('s S s_s u', positive=True)                       # s = k eta, S = k eta_o, s_s = k eta_*, u = S - s
chi = sp.Symbol('chi', positive=True)                                      # chi = eta_o - eta_*
sigma2, k0 = sp.symbols('sigma2 k0', positive=True)

j0 = lambda z: sp.sin(z) / z
j1 = lambda z: sp.sin(z) / z**2 - sp.cos(z) / z
j2 = lambda z: (3 / z**3 - 1 / z) * sp.sin(z) - 3 * sp.cos(z) / z**2

T_Phi = -2 * j1(x) / x                                  # = -2 (sin x - x cos x)/x^3
T_dT = -sp.cos(x) + 2 * j0(x) + T_Phi                   # = delta1/(4 zeta)
A_SW = T_dT + T_Phi                                     # Sachs-Wolfe-type emission source (dT + Psi)_e / zeta
kT_V = SQ3 * (sp.sin(x) - 2 * j1(x))                    # k T_V  (T_V = eta (x^2 sin x + 2x cos x - 2 sin x)/x^3)
T_V_eta = eta * (x**2 * sp.sin(x) + 2 * x * sp.cos(x) - 2 * sp.sin(x)) / x**3
dTPhi_deta_over_k = (4 / SQ3) * j2(x) / x               # (1/k) * 2 d_eta T_Phi  = (4/sqrt3) j_2(x)/x   [ISW source]

def X(kk, ee):
    return kk * ee / SQ3


def verify_transfer_functions(log=print):
    """Check the closed forms above against C1's first_order.pkl and the first-order equations."""
    sys.path.insert(0, os.path.join(HERE, '..', 'C1_einstein'))
    import kernels_api
    K1 = kernels_api.load()['first']
    xs = sp.Symbol('x', positive=True)
    ok = True
    for name, mine, theirs in (('T_Phi', T_Phi, K1['T_Phi']), ('T_deltaT', T_dT, K1['T_deltaT']),
                               ('T_V', T_V_eta.subs(eta, SQ3 * x / k), K1['T_V'])):
        d = sp.simplify(mine - theirs.subs(xs, x))
        log("  %-9s matches C1: %s" % (name, d == 0)); ok &= (d == 0)
    # equations (C1 M2): Phi'' + 4 Phi'/eta + k^2 Phi/3 = 0 ; delta = -(2/3)eta^2k^2 Phi - 2 eta Phi' - 2 Phi ; V = -(eta^2/2)(Phi' + Phi/eta)
    Phi = T_Phi.subs(x, X(k, eta))
    e1 = sp.simplify(sp.diff(Phi, eta, 2) + 4 * sp.diff(Phi, eta) / eta + k**2 * Phi / 3)
    e2 = sp.simplify(4 * T_dT.subs(x, X(k, eta)) - (-sp.Rational(2, 3) * eta**2 * k**2 * Phi - 2 * eta * sp.diff(Phi, eta) - 2 * Phi))
    e3 = sp.simplify(T_V_eta.subs(x, X(k, eta)) - (-(eta**2 / 2) * (sp.diff(Phi, eta) + Phi / eta)))
    # continuity for radiation: delta' = -(4/3) d_i v^i + 4 Phi' with d_i v^i -> -k^2 V
    e4 = sp.simplify(sp.diff(4 * T_dT.subs(x, X(k, eta)), eta) - (sp.Rational(4, 3) * k**2 * T_V_eta.subs(x, X(k, eta)) + 4 * sp.diff(Phi, eta)))
    e5 = sp.simplify(2 * sp.diff(Phi, eta) - k * dTPhi_deta_over_k.subs(x, X(k, eta)))
    for nm, e in (('ODE', e1), ('00 constraint', e2), ('0i constraint', e3), ('continuity', e4), ('ISW source', e5)):
        log("  %-14s residual = %s" % (nm, e)); ok &= (e == 0)
    lims = (sp.limit(A_SW, x, 0), sp.limit(kT_V / x, x, 0), sp.limit(T_dT, x, 0), sp.limit(T_Phi, x, 0))
    log("  limits x->0: A = T_dT + T_Phi -> %s, kT_V/x -> %s, T_dT -> %s, T_Phi -> %s" % lims)
    ok &= lims == (-sp.Rational(1, 3), SQ3 / 3, sp.Rational(1, 3), -sp.Rational(2, 3))
    return ok


# ----------------------------------------------------------------------------------------------------------------
# 2. Rayleigh rules and their numerical verification
#    e^{i k chi n.khat} = 4 pi sum_lm i^l j_l(k chi) Y_lm(n) Y*_lm(khat)
#    (n.khat) e^{i k chi n.khat} = (1/(i k)) d/dchi e^{...}  ->  n.khat -> -i j_l'(k chi)
# ----------------------------------------------------------------------------------------------------------------
def _sphere_grid(nth=96, nph=192):
    mu, w = leggauss(nth)
    th = np.arccos(mu)
    ph = 2 * np.pi * np.arange(nph) / nph
    TH, PH = np.meshgrid(th, ph, indexing='ij')
    W = np.outer(w, np.full(nph, 2 * np.pi / nph))
    return TH, PH, W


def alm_numeric(f, l, m, grid=None):
    """a_lm = int dOmega Y*_lm(n) f(theta, phi) by Gauss-Legendre x trapezoid quadrature."""
    TH, PH, W = grid or _sphere_grid()
    return np.sum(W * np.conj(sph_harm_y(l, m, TH, PH)) * f(TH, PH))


def rayleigh_check(kchi=5.3, khat=(0.3, -0.5, 0.81), lmax=6, log=print):
    """Numerically verify  e^{i kchi n.khat} -> 4pi i^l j_l Y*_lm(khat)  and  (n.khat) e^{...} -> 4pi i^l (-i j_l') Y*_lm(khat)."""
    khat = np.array(khat, float); khat /= np.linalg.norm(khat)
    thk, phk = np.arccos(khat[2]), np.arctan2(khat[1], khat[0])
    grid = _sphere_grid()
    TH, PH, W = grid
    nvec = np.stack([np.sin(TH) * np.cos(PH), np.sin(TH) * np.sin(PH), np.cos(TH)], -1)
    mu = nvec @ khat
    worst = 0.0
    for l in range(lmax + 1):
        for m in range(-l, l + 1):
            Yk = np.conj(sph_harm_y(l, m, thk, phk))
            a0 = alm_numeric(lambda T, P: np.exp(1j * kchi * mu), l, m, grid)
            a1 = alm_numeric(lambda T, P: mu * np.exp(1j * kchi * mu), l, m, grid)
            p0 = 4 * np.pi * 1j**l * Yk * spherical_jn(l, kchi)
            p1 = 4 * np.pi * 1j**l * Yk * (-1j) * spherical_jn(l, kchi, derivative=True)
            worst = max(worst, abs(a0 - p0), abs(a1 - p1))
    # observer Doppler: n.khat (no phase) is pure l=1 with a_1m = 4pi i (1/3) Y*_1m(khat)  [since j_1'(0) = 1/3]
    for l in range(0, 4):
        for m in range(-l, l + 1):
            a = alm_numeric(lambda T, P: mu + 0j, l, m, grid)
            pred = 4 * np.pi * 1j * np.conj(sph_harm_y(1, m, thk, phk)) / 3 if l == 1 else 0.0
            worst = max(worst, abs(a - pred))
    log("  Rayleigh rules (1 -> i^l j_l, n.khat -> -i i^l j_l', observer n.khat -> l=1 only, 1/3): max |err| = %.2e" % worst)
    return worst < 1e-10


# ----------------------------------------------------------------------------------------------------------------
# 3. Closed-form ISW integral  I_l(k) = (4/sqrt3) int_{s_*}^{S} [j_2(s/sqrt3)/(s/sqrt3)] j_l(S - s) ds
#    general-l algorithm:
#      (i)  j_l(u) = C_l(u) sin u + D_l(u) cos u  (Laurent polynomials, A&S 10.1.8-9), j_2(x)/x likewise in x = s/sqrt3
#      (ii) product -> sum:  phases  cos/sin(a_± s - S),  a_± = 1 ± 1/sqrt3   (no resonance: a_- != 0)
#      (iii) partial fractions 1/(s^p u^q), u = S - s   (binomial formula)
#      (iv) antiderivatives  C_m[a](t) = int cos(a t)/t^m dt,  S_m[a](t) = int sin(a t)/t^m dt  by the recursion
#             C_1 = Ci(a t), S_1 = Si(a t),
#             C_m = -cos(a t)/((m-1) t^{m-1}) - a/(m-1) S_{m-1},   S_m = -sin(a t)/((m-1) t^{m-1}) + a/(m-1) C_{m-1}
#      (v)  I_l = [F_s]_{s_*}^{S} + [F_u]_{u = S - s_*}^{u -> 0^+}, the u -> 0 endpoint (observer) by its finite part:
#             poles and log(u) cancel because the integrand is regular at s = S (j_l entire).
# ----------------------------------------------------------------------------------------------------------------
A_PLUS, A_MINUS = 1 + 1 / SQ3, 1 - 1 / SQ3
PHASES = {'+': A_PLUS, '-': A_MINUS}


def jl_trig(l, var):
    """Laurent polynomials (C, D) with j_l(var) = C sin(var) + D cos(var)."""
    e = sp.expand(sp.expand_func(sp.jn(l, var)))
    C, D = e.coeff(sp.sin(var)), e.coeff(sp.cos(var))
    assert sp.simplify(e - C * sp.sin(var) - D * sp.cos(var)) == 0
    return C, D


def laurent_terms(expr, var):
    """{exponent: coefficient} of a Laurent polynomial in var."""
    out = {}
    for t in sp.Add.make_args(sp.expand(expr)):
        c, e = t.as_coeff_exponent(var)
        out[int(e)] = out.get(int(e), 0) + c
    return out


def partial_fractions(p, q):
    """1/(s^p (S-s)^q) = sum_i c_i/s^i + sum_j d_j/(S-s)^j ; returns ({i: c_i}, {j: d_j}) with S-dependence."""
    cs = {i: sp.binomial(p + q - i - 1, q - 1) * S**(-(p + q - i)) for i in range(1, p + 1)}
    ds = {j: sp.binomial(p + q - j - 1, p - 1) * S**(-(p + q - j)) for j in range(1, q + 1)}
    return cs, ds


@functools.lru_cache(maxsize=None)
def CS_antiderivative(kind, m, a_key, var):
    """C_m[a](var) (kind 'c') or S_m[a](var) (kind 's'), recursion above; a = PHASES[a_key]."""
    a = PHASES[a_key]
    if m == 1:
        return sp.Ci(a * var) if kind == 'c' else sp.Si(a * var)
    if kind == 'c':
        return -sp.cos(a * var) / ((m - 1) * var**(m - 1)) - a / (m - 1) * CS_antiderivative('s', m - 1, a_key, var)
    return -sp.sin(a * var) / ((m - 1) * var**(m - 1)) + a / (m - 1) * CS_antiderivative('c', m - 1, a_key, var)


def isw_integrand_phases(l):
    """Coefficients P_a, Q_a (rational in s, u) with
       (4/sqrt3)[j_2(x)/x] j_l(u) = sum_a [P_a cos(a s - S) + Q_a sin(a s - S)],  x = s/sqrt3, u = S - s."""
    C2, D2 = jl_trig(2, x)
    A2 = sp.expand((C2 / x).subs(x, s / SQ3)); B2 = sp.expand((D2 / x).subs(x, s / SQ3))
    Cl, Dl = jl_trig(l, u)
    pref = 4 / SQ3
    P = {'+': pref * (A2 * Cl + B2 * Dl) / 2, '-': pref * (-A2 * Cl + B2 * Dl) / 2}
    Q = {'+': pref * (A2 * Dl - B2 * Cl) / 2, '-': pref * (-(A2 * Dl + B2 * Cl)) / 2}
    return P, Q


def _check_phase_decomposition(l, log=print):
    P, Q = isw_integrand_phases(l)
    direct = (4 / SQ3) * (j2(x) / x).subs(x, s / SQ3) * sp.expand_func(sp.jn(l, u))
    rec = sum(P[a] * sp.cos(PHASES[a] * s - S) + Q[a] * sp.sin(PHASES[a] * s - S) for a in PHASES)
    vals = {s: sp.Rational(7, 3), S: sp.Rational(19, 4)}
    d = (direct - rec).subs(u, S - s).subs(vals).evalf(30)
    log("  product-to-sum identity l=%d: residual %s" % (l, d))
    return abs(d) < 1e-25


@functools.lru_cache(maxsize=None)
def isw_antiderivatives(l):
    """Return (F_s(s), F_u(u)) sympy expressions: antiderivative (w.r.t. s) of the ISW integrand, split into the
       part with poles at s = 0 (function of s) and the part with poles at u = 0 (function of u = S - s)."""
    P, Q = isw_integrand_phases(l)
    coef = {}   # (var, kind, m, a_key) -> coefficient (function of S)

    def add(key, c):
        coef[key] = coef.get(key, 0) + c

    for a_key, a in PHASES.items():
        b = (a - 1) * S                                     # a s - S = b - a u
        for trig, R in (('cos', P[a_key]), ('sin', Q[a_key])):
            terms_s = laurent_terms(R, s)
            for es, cs_ in terms_s.items():
                for eu, c in laurent_terms(cs_, u).items():
                    p, q = -es, -eu
                    assert p >= 1 and q >= 1, (p, q)
                    cs, ds = partial_fractions(p, q)
                    for i, ci in cs.items():                 # c ci trig(a s - S)/s^i
                        cc = c * ci
                        if trig == 'cos':                    # cos(as-S) = cos S cos(as) + sin S sin(as)
                            add(('s', 'c', i, a_key), cc * sp.cos(S)); add(('s', 's', i, a_key), cc * sp.sin(S))
                        else:                                # sin(as-S) = cos S sin(as) - sin S cos(as)
                            add(('s', 's', i, a_key), cc * sp.cos(S)); add(('s', 'c', i, a_key), -cc * sp.sin(S))
                    for j, dj in ds.items():                 # c dj trig(b - a u)/u^j ; d s = -d u
                        cc = -c * dj
                        if trig == 'cos':                    # cos(b-au) = cos b cos(au) + sin b sin(au)
                            add(('u', 'c', j, a_key), cc * sp.cos(b)); add(('u', 's', j, a_key), cc * sp.sin(b))
                        else:                                # sin(b-au) = sin b cos(au) - cos b sin(au)
                            add(('u', 'c', j, a_key), cc * sp.sin(b)); add(('u', 's', j, a_key), -cc * sp.cos(b))
    Fs = sum((c * CS_antiderivative(kind, m, a_key, s) for (var, kind, m, a_key), c in coef.items() if var == 's'), sp.Integer(0))
    Fu = sum((c * CS_antiderivative(kind, m, a_key, u) for (var, kind, m, a_key), c in coef.items() if var == 'u'), sp.Integer(0))
    return sp.expand(Fs), sp.expand(Fu)


def finite_part(F, var, check=True):
    """Constant term of the Laurent/log expansion of F (sum of terms  c var^n trig(a var),  c var^n,  c Ci(a var), c Si(a var))
       as var -> 0^+.  Returns (finite_part, singular_dict) where singular_dict holds the coefficients of var^{-p} and log(var)
       (all of which must vanish for a regular antiderivative)."""
    fp = sp.Integer(0); sing = {}
    for t in sp.Add.make_args(sp.expand(F)):
        c, n = sp.Integer(1), 0; kind = None; a = None
        for f in sp.Mul.make_args(t):
            if f.has(var):
                if f.is_Pow and f.base == var:
                    n += int(f.exp)
                elif f == var:
                    n += 1
                elif isinstance(f, (sp.sin, sp.cos, sp.Ci, sp.Si)):
                    kind = type(f); a = sp.simplify(f.args[0] / var); assert not a.has(var)
                else:
                    raise ValueError("unexpected factor %s" % f)
            else:
                c *= f
        # Taylor/Laurent data of the var-dependent factor
        if kind is None:
            series = {0: sp.Integer(1)}
        elif kind is sp.cos:
            series = {2 * j: (-1)**j * a**(2 * j) / sp.factorial(2 * j) for j in range(0, (-n) // 2 + 1)}
        elif kind is sp.sin:
            series = {2 * j + 1: (-1)**j * a**(2 * j + 1) / sp.factorial(2 * j + 1) for j in range(0, (-n) // 2 + 1)}
        elif kind is sp.Si:
            series = {2 * j + 1: (-1)**j * a**(2 * j + 1) / ((2 * j + 1) * sp.factorial(2 * j + 1)) for j in range(0, (-n) // 2 + 1)}
        elif kind is sp.Ci:        # Ci(a v) = EulerGamma + log(a) + log(v) + sum (-1)^j (a v)^{2j}/(2j (2j)!)
            assert n == 0
            series = {0: sp.EulerGamma + sp.log(a)}
            sing['log'] = sing.get('log', 0) + c
        for j, cj in series.items():
            e = n + j
            if e == 0:
                fp += c * cj
            elif e < 0:
                sing[e] = sing.get(e, 0) + c * cj
    if check:
        for key, v in sing.items():
            vv = sp.simplify(v)
            if vv != 0:
                # fall back to a numerical test at random points
                num = max(abs(sp.N(vv.subs(S, val), 30)) for val in (sp.Rational(11, 7), sp.Rational(23, 3), sp.pi))
                assert num < 1e-25, ("singular term does not cancel", key, vv)
    return sp.expand(fp), sing


@functools.lru_cache(maxsize=None)
def I_l_closed(l):
    """Closed form of I_l(k) as a sympy expression in s_s = k eta_* and S = k eta_o (elementary + Ci/Si)."""
    Fs, Fu = isw_antiderivatives(l)
    Fu0, _ = finite_part(Fu, u)
    expr = (Fs.subs(s, S) - Fs.subs(s, sst)) + (Fu0 - Fu.subs(u, S - sst))
    return sp.expand(expr)


@functools.lru_cache(maxsize=None)
def I_l_closed_eta0(l):
    """Closed form of I_l with eta_* -> 0 (lower endpoint regular, finite part at s -> 0)."""
    Fs, Fu = isw_antiderivatives(l)
    Fs0, _ = finite_part(Fs, s)
    Fu0, _ = finite_part(Fu, u)
    return sp.expand(Fs.subs(s, S) - Fs0 + Fu0 - Fu.subs(u, S))


def I_l_closed_collected(l):
    """I_l organised as  R + sum_{±} [ coefficients x Ci/Si of the six arguments a_± s_*, a_± S, a_± (S - s_*) ]."""
    e = I_l_closed(l)
    trans = [f for f in e.atoms(sp.Ci, sp.Si)]
    return sp.collect(e, trans)


# --- numerical evaluators of I_l -----------------------------------------------------------------------------------
def _mp_jl(l, z):
    z = mp.mpf(z)
    if z == 0:
        return mp.mpf(1) if l == 0 else mp.mpf(0)
    return mp.sqrt(mp.pi / (2 * z)) * mp.besselj(l + mp.mpf(1) / 2, z)


def I_l_quad_mp(l, kk, es, eo, dps=30):
    """I_l(k) by mpmath quadrature of the eta-integrand 2 d_eta T_Phi j_l(k(eta_o - eta))."""
    with mp.workdps(dps):
        kk, es, eo = mp.mpf(kk), mp.mpf(es), mp.mpf(eo)
        f = lambda e: kk * (4 / mp.sqrt(3)) * _mp_j2_over_x(kk * e / mp.sqrt(3)) * _mp_jl(l, kk * (eo - e))
        npts = int(kk * (eo - es) / 2) + 2
        pts = [es + (eo - es) * i / npts for i in range(npts + 1)]
        return mp.quad(f, pts)


def _mp_j2_over_x(xx):
    if xx == 0:
        return mp.mpf(0)
    return ((3 / xx**3 - 1 / xx) * mp.sin(xx) - 3 * mp.cos(xx) / xx**2) / xx


@functools.lru_cache(maxsize=None)
def _I_l_closed_mp_fn(l):
    return sp.lambdify((sst, S), I_l_closed(l), 'mpmath')


def I_l_closed_mp(l, kk, es, eo, dps=30):
    with mp.workdps(dps):
        return _I_l_closed_mp_fn(l)(mp.mpf(kk) * mp.mpf(es), mp.mpf(kk) * mp.mpf(eo))


def I_l_quad_np(l, kk, es, eo, nodes=None):
    """Vectorised (over l array) Gauss-Legendre evaluation of I_l(k); robust default numerical route."""
    l = np.atleast_1d(l)
    n = nodes or int(40 + 6 * kk * (eo - es))
    t, w = leggauss(n)
    e = 0.5 * (eo - es) * t + 0.5 * (eo + es); w = 0.5 * (eo - es) * w
    xx = kk * e / np.sqrt(3)
    src = kk * (4 / np.sqrt(3)) * spherical_jn(2, xx) / xx
    return np.array([np.sum(w * src * spherical_jn(ll, kk * (eo - e))) for ll in l])


# ----------------------------------------------------------------------------------------------------------------
# 4. Delta_l(k), C_l^(11), monochromatic, monopole
# ----------------------------------------------------------------------------------------------------------------
def G_l_sym(l, ee):
    """G_l(eta) = A(x) j_l(k chi) - kT_V(x) j_l'(k chi),  chi = eta_o - eta (sympy)."""
    z = sp.Symbol('z', positive=True)
    jl = sp.expand_func(sp.jn(l, z)); jlp = sp.diff(jl, z)
    xx = X(k, ee)
    return (A_SW.subs(x, xx) * jl.subs(z, k * (eta_o - ee)) - kT_V.subs(x, xx) * jlp.subs(z, k * (eta_o - ee)))


def Delta_l_sym(l, isw='symbolic'):
    """Delta_l(k) for general l >= 0 as a sympy expression (observer terms included: l=0 monopole, l=1 Doppler)."""
    em = G_l_sym(l, eta_s)
    obs = sp.Integer(0)
    if l == 0:
        obs = A_SW.subs(x, X(k, eta_o))
    if l == 1:
        obs = -kT_V.subs(x, X(k, eta_o)) / 3
    if isw == 'symbolic':
        I = sp.Integral(2 * sp.diff(T_Phi.subs(x, X(k, eta)), eta) * sp.expand_func(sp.jn(l, k * (eta_o - eta))), (eta, eta_s, eta_o))
    else:
        I = I_l_closed(l).subs({sst: k * eta_s, S: k * eta_o})
    return em - obs + I


def Delta_l_num(l, kk, es, eo, isw='quad'):
    """Numerical Delta_l(k) (float, vectorised over l) incl. observer terms for l = 0, 1.  isw: 'quad' | 'closed' | 'none'."""
    l = np.atleast_1d(l).astype(int)
    xs_ = kk * es / np.sqrt(3); xo = kk * eo / np.sqrt(3); c = kk * (eo - es)
    A = lambda xx: -np.cos(xx) + 2 * spherical_jn(0, xx) - 4 * spherical_jn(1, xx) / xx
    kTV = lambda xx: np.sqrt(3) * (np.sin(xx) - 2 * spherical_jn(1, xx))
    em = A(xs_) * spherical_jn(l, c) - kTV(xs_) * spherical_jn(l, c, derivative=True)
    obs = np.where(l == 0, A(xo), 0.0) + np.where(l == 1, -kTV(xo) / 3, 0.0)
    if isw == 'quad':
        I = I_l_quad_np(l, kk, es, eo)
    elif isw == 'closed':
        I = np.array([float(I_l_closed_mp(int(ll), kk, es, eo)) for ll in l])
    else:
        I = 0.0
    return em - obs + I


def C11_general(l, P, kmin, kmax, es, eo, **kw):
    """C_l^(11) = (2/pi) int k^2 dk P(k) Delta_l(k)^2  (scipy quad; P a python callable)."""
    f = lambda kk: kk**2 * P(kk) * float(Delta_l_num(l, kk, es, eo)[0])**2
    val, err = quad(f, kmin, kmax, limit=2000, **kw)
    return 2 / np.pi * val


def C11_mono(l, sig2, kk0, es, eo, **kw):
    """Monochromatic P(k) = 2 pi^2 sigma^2 delta(k-k0)/k0^2  =>  C_l^(11) = 4 pi sigma^2 Delta_l(k0)^2."""
    return 4 * np.pi * sig2 * Delta_l_num(l, kk0, es, eo, **kw)**2


def M_sym(isw='closed'):
    """Sky-monopole kernel M(k) = Delta_0(k)  (<m1> = int d^3k/(2pi)^3 zeta_k M(k))."""
    return Delta_l_sym(0, isw)


def M_num(kk, es, eo, **kw):
    return Delta_l_num(0, kk, es, eo, **kw)[0]


def var_m1_mono(sig2, kk0, es, eo):
    """<<m1>^2> = int d^3k/(2pi)^3 P |M|^2 -> sigma^2 M(k0)^2 (monochromatic)."""
    return sig2 * M_num(kk0, es, eo)**2


# ----------------------------------------------------------------------------------------------------------------
# 5. real-space check of C_l^(11) (monochromatic), including observer terms and ISW
# ----------------------------------------------------------------------------------------------------------------
def kernel_K(kk, khat, nvec, es, eo, neta=160):
    """First-order kernel K(k, n) (Theta1 = int d^3k/(2pi)^3 zeta_k K): emission + observer + ISW, x_e = +n chi_*."""
    mu = nvec @ khat
    xs_ = kk * es / np.sqrt(3); xo = kk * eo / np.sqrt(3); c = kk * (eo - es)
    A = lambda xx: -np.cos(xx) + 2 * spherical_jn(0, xx) - 4 * spherical_jn(1, xx) / xx
    kTV = lambda xx: np.sqrt(3) * (np.sin(xx) - 2 * spherical_jn(1, xx))
    Ke = (A(xs_) - 1j * kTV(xs_) * mu) * np.exp(1j * c * mu)          # [dT + Psi - n.v]_e,  n.v = i k mu T_V zeta
    Ko = A(xo) - 1j * kTV(xo) * mu                                     # [dT + Psi - n.v]_o
    t, w = leggauss(neta); e = 0.5 * (eo - es) * t + 0.5 * (eo + es); w = 0.5 * (eo - es) * w
    xx = kk * e / np.sqrt(3); src = kk * (4 / np.sqrt(3)) * spherical_jn(2, xx) / xx
    Ki = np.sum(w[:, None] * src[:, None] * np.exp(1j * kk * np.outer(eo - e, mu)), 0)
    return Ke - Ko + Ki


def realspace_check(kk0=4.0, es=0.5, eo=2.0, lmax=30, log=print):
    """xi(n1.n2) = sigma^2 <K(k0 khat, n1) K*(k0 khat, n2)>_khat  vs  sum_l (2l+1)/(4pi) C_l P_l(n1.n2), C_l = 4pi sigma^2 Delta_l^2."""
    TH, PH, W = _sphere_grid(64, 128)
    khat = np.stack([np.sin(TH) * np.cos(PH), np.sin(TH) * np.sin(PH), np.cos(TH)], -1).reshape(-1, 3)
    W = W.reshape(-1) / (4 * np.pi)
    ls = np.arange(lmax + 1)
    Cl = 4 * np.pi * Delta_l_num(ls, kk0, es, eo)**2            # sigma^2 = 1
    n1 = np.array([0, 0, 1.0]); worst = 0
    for theta in (0.0, 0.4, 1.1, 2.0, np.pi):
        n2 = np.array([np.sin(theta), 0, np.cos(theta)])
        K1 = np.array([kernel_K(kk0, kh, n1, es, eo) for kh in khat])
        K2 = np.array([kernel_K(kk0, kh, n2, es, eo) for kh in khat])
        xi = np.sum(W * K1 * np.conj(K2))
        pred = np.sum((2 * ls + 1) / (4 * np.pi) * Cl * eval_legendre(ls, np.cos(theta)))
        worst = max(worst, abs(xi - pred))
        log("  theta=%.2f  xi_direct=% .10f%+.1ej   harmonic sum=% .10f" % (theta, xi.real, xi.imag, pred))
    log("  real-space vs harmonic C_l^(11) (incl. observer terms, l=0,1, and ISW): max diff %.2e" % worst)
    return worst < 1e-8


# ----------------------------------------------------------------------------------------------------------------
# 6. limits
# ----------------------------------------------------------------------------------------------------------------
def sw_limit_sym(l):
    """k eta_* -> 0 at fixed k eta_o:  Delta_l -> -(1/3) j_l(k eta_o) + I_l(k; eta_* = 0)  (l >= 2)."""
    z = sp.Symbol('z', positive=True)
    return -sp.Rational(1, 3) * sp.expand_func(sp.jn(l, z)).subs(z, S) + I_l_closed_eta0(l)


def plateau_check(log=print):
    """Scale-invariant P = 2 pi^2 A_s/k^3 with Delta_l = -(1/3) j_l(k chi):  C_l = (4 pi A_s/9) int dk/k j_l^2 = (4 pi A_s/9)/(2 l (l+1))
       => l(l+1) C_l/(2 pi) = A_s/9   (RD analogue of the matter-era A_s/25)."""
    ok = True
    for l in (2, 3, 7):
        val = quad(lambda y: spherical_jn(l, y)**2 / y, 1e-6, 400, limit=4000)[0]
        val += quad(lambda y: spherical_jn(l, y)**2 / y, 400, 20000, limit=20000)[0]
        log("  int dy/y j_%d(y)^2 = %.8f  vs 1/(2l(l+1)) = %.8f" % (l, val, 1 / (2 * l * (l + 1))))
        ok &= abs(val - 1 / (2 * l * (l + 1))) < 1e-5
    return ok


# ----------------------------------------------------------------------------------------------------------------
# 7. test driver
# ----------------------------------------------------------------------------------------------------------------
def run_tests(lmax_closed=6):
    lines = []

    def log(msg=''):
        print(msg); lines.append(msg)
    t0 = time.time()
    log("== C5 linear_cl test suite ==")
    log("[1] transfer functions vs C1 and first-order equations")
    ok1 = verify_transfer_functions(log); log("  -> %s" % ok1)
    log("[2] Rayleigh rules (numerical angular quadrature)")
    ok2 = rayleigh_check(log=log); log("  -> %s" % ok2)
    log("[3] ISW closed form: product-to-sum, partial fractions, finite part, vs mpmath quadrature (12 digits)")
    ok3 = all(_check_phase_decomposition(l, log) for l in (0, 1, 2, 5))
    for p, q in ((2, 1), (4, 3), (3, 7)):
        cs, ds = partial_fractions(p, q)
        rec = sum(c / s**i for i, c in cs.items()) + sum(d / (S - s)**j for j, d in ds.items())
        r = sp.simplify(rec - 1 / (s**p * (S - s)**q)); log("  partial fractions (p,q)=(%d,%d): residual %s" % (p, q, r)); ok3 &= (r == 0)
    cases = [(1.0, 0.3, 2.0), (3.7, 0.2, 2.5), (0.4, 1.0, 9.0), (6.0, 0.5, 4.0), (2.2, 0.05, 1.0)]
    worst = 0
    for l in range(0, lmax_closed + 1):
        t1 = time.time(); e = I_l_closed(l); nterms = len(sp.Add.make_args(e))
        # finite part sanity: F_u at tiny u (high precision) must equal the finite part
        Fs, Fu = isw_antiderivatives(l); Fu0, sing = finite_part(Fu, u)
        with mp.workdps(60):
            fU = sp.lambdify((u, S), Fu, 'mpmath'); f0 = sp.lambdify((S,), Fu0, 'mpmath')
            dfp = abs(fU(mp.mpf('1e-20'), mp.mpf('3.1')) - f0(mp.mpf('3.1')))
        for (kk, es, eo) in cases:
            a = I_l_closed_mp(l, kk, es, eo); b = I_l_quad_mp(l, kk, es, eo)
            c = I_l_quad_np(l, kk, es, eo)[0]
            rel = abs(a - b) / max(abs(b), mp.mpf('1e-30')); worst = max(worst, rel)
            worst = max(worst, abs(float(a) - c) / max(abs(c), 1e-30) * 1e-4)   # numpy route ~1e-12 anyway
        log("  l=%d: %d terms, singular pieces cancelled (%d checked), finite-part check %.1e, max rel err closed vs quad %.1e  [%.1fs]"
            % (l, nterms, len(sing), dfp, worst, time.time() - t1))
    ok3 &= worst < 1e-12
    # eta_* -> 0 closed form vs quadrature
    for l in (2, 3):
        a = sp.lambdify((S,), I_l_closed_eta0(l), 'mpmath')(mp.mpf('5.5'))
        b = I_l_quad_mp(l, 1.0, 0.0, 5.5)
        log("  eta_*=0 closed form l=%d: %s vs quad %s" % (l, mp.nstr(a, 15), mp.nstr(b, 15))); ok3 &= abs(a - b) < 1e-12
    log("  -> %s" % ok3)
    log("[4] C_l^(11) real-space check (monochromatic; verifies prefactor 4 pi sigma^2 |Delta_l|^2, observer l=0,1 terms, ISW)")
    ok4 = realspace_check(log=log); log("  -> %s" % ok4)
    log("[5] limits")
    # k->0 of Delta_l (dilation): M(k) = Delta_0 -> 0 as O(k^2); Delta_1 -> O(k^3)
    for kk in (1e-2, 1e-3):
        log("  k=%.0e: M(k)=Delta_0 = %.3e (expect O(k^2)),  Delta_1 = %.3e (expect O(k^3)),  Delta_2 = %.3e" %
            (kk, M_num(kk, 0.5, 2.0), Delta_l_num(1, kk, 0.5, 2.0)[0], Delta_l_num(2, kk, 0.5, 2.0)[0]))
    ok5 = abs(M_num(1e-3, 0.5, 2.0)) < 1e-5 and abs(Delta_l_num(1, 1e-3, 0.5, 2.0)[0]) < 1e-8
    Msmall = sp.series(Delta_l_sym(0, 'closed').subs({eta_s: sp.Rational(1, 2), eta_o: 2}), k, 0, 3).removeO()
    log("  M(k) small-k series (eta_*=1/2, eta_o=2): %s" % sp.nsimplify(sp.simplify(Msmall)))
    ok5 &= plateau_check(log)
    # k -> infinity: Delta_l ~ -cos x_* j_l - sqrt3 sin x_* j_l'  (ISW O(1/k))
    for kk in (200.0, 2000.0):
        es, eo = 0.5, 2.0; xs_ = kk * es / np.sqrt(3); c = kk * (eo - es)
        full = Delta_l_num(2, kk, es, eo)[0]
        asym = -np.cos(xs_) * spherical_jn(2, c) - np.sqrt(3) * np.sin(xs_) * spherical_jn(2, c, derivative=True)
        log("  k=%.0f: Delta_2 = %.6e, asymptotic = %.6e, diff*k^2 = %.3f" % (kk, full, asym, (full - asym) * kk**2))
    log("  -> %s" % ok5)
    log("[6] example numbers (sigma^2 = 1, eta_* = 0.5, eta_o = 2, k0 = 4): C_l^(11) = 4 pi Delta_l(k0)^2")
    ls = np.arange(0, 9); D = Delta_l_num(ls, 4.0, 0.5, 2.0)
    for l, d in zip(ls, D):
        log("   l=%d  Delta_l = % .8f   C_l = %.8f" % (l, d, 4 * np.pi * d**2))
    log("  <m1> variance (mono): sigma^2 M(k0)^2 = %.8f" % var_m1_mono(1.0, 4.0, 0.5, 2.0))
    allok = ok1 and ok2 and ok3 and ok4 and ok5
    log("== ALL PASS: %s   (%.0fs) ==" % (allok, time.time() - t0))
    with open(os.path.join(HERE, 'test_log.txt'), 'w') as fh:
        fh.write("\n".join(lines) + "\n")
    return allok


if __name__ == '__main__':
    run_tests()
