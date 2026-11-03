"""kernel2.py -- Fourier-space map kernels K_1(k;n), K_2(k1,k2;n) assembled from C2's line-of-sight term lists
(L_terms.pkl, n = z-hat 'split' lists, shear form AND direct form) and C1's transfer functions / kernels.

Conventions (SPEC + C2): ray xbar(eta) = x_o + n chi, chi = eta_o - eta, x_o = 0, n = observation direction.
Fourier: f(x) = int d^3k/(2pi)^3 e^{ik.x} f_k, so d_j -> i k_j and a field on the ray carries e^{+i (k.n) chi}.
Theta1(n) = int d^3k/(2pi)^3 zeta_k K_1(k;n);  Theta2(n) = int int d^3k1 d^3k2/(2pi)^6 zeta_k1 zeta_k2 K_2(k1,k2;n),
Theta = eps Th1 + eps^2 Th2/2 (SPEC).  With Lambda = eps L1 + eps^2 L2/2, m2 = L2 + L1^2, Theta2 = m2 - 2<m1> m1 (l>=1):
   K_2 = 2 K_{L2/2} + K_1(k1)K_1(k2) - Delta_0(k1) K_1(k2) - Delta_0(k2) K_1(k1),   K_{L2/2} = symmetrised kernel of C2's list
   (C2's L[2] list IS the coefficient of eps^2, i.e. L2/2).  Delta_0 = sky-monopole transfer function of m1 = L1.
The engine works in the frame n = z-hat: for each direction n an orthonormal transverse pair (e1,e2) is built and every
momentum is represented by (e1.k, e2.k, n.k).  Time integrals: Chebyshev-Lobatto grid on [eta_*, eta_o] (spectral cumulative
integrals for the nested ray functionals), Gauss-Legendre for the Green's-function integrals int_0^eta and the observer-worldline
integrals int_0^eta_o.
"""
import os, re, pickle
import numpy as np
import sympy as sp
from numpy.polynomial import chebyshev as Ch, legendre as Lg
from scipy.special import spherical_jn
from transfer import tphi, tdt, tdelta, tv, KER, G_Phi, dG_Phi, G_h, dG_h, ghat, G_h0, dG_h0

HERE = os.path.dirname(os.path.abspath(__file__))
C2DIR = os.path.join(HERE, '..', 'C2_observable')
LT = pickle.load(open(os.path.join(C2DIR, 'L_terms.pkl'), 'rb'))
ETA, XS, YS, ZS = sp.symbols('eta x y z')
ETA_STAR = sp.Symbol('eta_star')

# ------------------------------------------------------------------ geometry / quadrature
class Geometry:
    def __init__(self, eta_s=2.0, eta_o=8.0, M=128, Ng=160, Nw=240):
        self.eta_s, self.eta_o, self.M = float(eta_s), float(eta_o), M
        xg = -np.cos(np.pi * np.arange(M) / (M - 1))                    # ascending Chebyshev-Lobatto nodes in [-1,1]
        self.t = eta_s + (xg + 1) * (eta_o - eta_s) / 2
        self.chi = self.eta_o - self.t
        V = Ch.chebvander(xg, M - 1)
        Cinv = np.linalg.inv(V)
        B = np.zeros((M, M))
        for j in range(M):
            e = np.zeros(M); e[j] = 1.0
            ci = Ch.chebint(e)
            B[:, j] = Ch.chebval(1.0, ci) - Ch.chebval(xg, ci)
        self.Cup = (B @ Cinv) * (eta_o - eta_s) / 2                   # (Cup f)_i = int_{t_i}^{eta_o} f
        self.w = self.Cup[0].copy()                                     # int_{eta_s}^{eta_o} f = w . f
        # observer worldline int_0^eta_o
        u, wu = Lg.leggauss(Nw)
        self.tw = (u + 1) * eta_o / 2; self.ww = wu * eta_o / 2
        # Green's integrals int_0^{t_i}: GL nodes per i
        ug, wg = Lg.leggauss(Ng)
        self.gs = (ug[None, :] + 1) * self.t[:, None] / 2             # (M, Ng)
        self.gw = wg[None, :] * self.t[:, None] / 2
    def cum(self, f):
        """int_{t}^{eta_o} along the last axis of f (..., M)."""
        return f @ self.Cup.T
    def integ(self, f):
        return f @ self.w

def transverse_basis(nv):
    """orthonormal (e1, e2, n) for an array of unit vectors nv (N,3)."""
    nv = np.asarray(nv, dtype=float)
    a = np.zeros_like(nv); a[:, 0] = 1.0
    bad = np.abs(nv[:, 0]) > 0.9
    a[bad] = 0; a[bad, 1] = 1.0
    e1 = a - (np.sum(a * nv, axis=1))[:, None] * nv
    e1 /= np.linalg.norm(e1, axis=1)[:, None]
    e2 = np.cross(nv, e1)
    return e1, e2

# ------------------------------------------------------------------ parsing of C2 monomials
_re_scalar = re.compile(r'^(Psi|Phi|dT)(\d)$')
_re_vec = re.compile(r'^(v|om)_(\d)_(\d)$')
_re_h = re.compile(r'^h_(\d)(\d)_(\d)$')
_re_fun = re.compile(r'^(k1|Dx1)_(\d)$')
def parse_factor(f):
    """-> descriptor tuple or None (numeric)."""
    der = [0, 0, 0, 0]
    if isinstance(f, sp.Derivative):
        for v_, c_ in f.variable_count:
            der[(ETA, XS, YS, ZS).index(v_)] += int(c_)
        f = f.expr
    if isinstance(f, sp.core.function.AppliedUndef):
        nm = f.func.__name__
        m = _re_scalar.match(nm)
        if m: return ('field', m.group(1), None, int(m.group(2)), tuple(der))
        m = _re_vec.match(nm)
        if m: return ('field', m.group(1), int(m.group(2)) - 1, int(m.group(3)), tuple(der))
        m = _re_h.match(nm)
        if m: return ('field', 'h', (int(m.group(1)) - 1, int(m.group(2)) - 1), int(m.group(3)), tuple(der))
        m = _re_fun.match(nm)
        if m:
            assert der == [0, 0, 0, 0]
            when = 'eta' if f.args[0] == ETA else 'eta_star'
            return ('fun', m.group(1), int(m.group(2)) - 1, when)
        raise ValueError(nm)
    if isinstance(f, sp.Symbol):
        s = f.name
        if s in ('deta_o1', 'deta_e1'): return ('sym', s, None)
        if s.startswith('dxo1_'): return ('sym', 'dxo1', int(s[-1]) - 1)
        if s == 'eta': return ('eta',)
        if s == 'eta_star': return ('const', None)   # handled as numeric below
        raise ValueError(s)
    return None

def order_of(desc):
    if desc[0] == 'field': return desc[3]
    if desc[0] in ('fun', 'sym'): return 1
    return 0

def parse_term(expr):
    """-> list of (coeff(float, may contain eta_star numerically later), [descriptors])"""
    out = []
    for mono in sp.Add.make_args(sp.expand(expr)):
        coeff = 1.0; descs = []
        for f in sp.Mul.make_args(mono):
            base, ex = (f.base, int(f.exp)) if isinstance(f, sp.Pow) else (f, 1)
            d = parse_factor(base)
            if d is None:
                coeff *= complex(base**ex)
            elif d[0] == 'const':
                descs += [('eta_star',)] * ex
            else:
                descs += [d] * ex
        out.append((coeff, descs))
    return out

# ------------------------------------------------------------------ momentum slot (first-order objects)
class Slot:
    """one momentum: magnitude kmag and z-frame components kv (N,3) for each direction n."""
    def __init__(self, geom, kvec, nv):
        self.g = geom
        kvec = np.asarray(kvec, dtype=float)
        self.kmag = float(np.linalg.norm(kvec))
        e1, e2 = transverse_basis(nv)
        self.kv = np.stack([e1 @ kvec, e2 @ kvec, nv @ kvec], axis=1)     # (N,3)
        self.N = self.kv.shape[0]
        kz = self.kv[:, 2]
        self.E_ray = np.exp(1j * kz[:, None] * geom.chi[None, :])          # (N,M)
        self.E_emi = np.exp(1j * kz * (geom.eta_o - geom.eta_s))
        self.cache = {}
    # scalar amplitudes (depend on magnitude only)
    def amp(self, name, T, d):
        k = self.kmag
        if name in ('Psi', 'Phi'): return tphi(k, T, d)
        if name == 'dT': return tdt(k, T, d)
        raise ValueError(name)
    def field(self, name, comp, der, where):
        key = ('field', name, comp, der, where)
        if key in self.cache: return self.cache[key]
        g = self.g
        me, mx, my, mz = der
        if where == 'ray': T, phase = g.t, self.E_ray
        elif where == 'obs': T, phase = g.eta_o, 1.0
        elif where == 'emi': T, phase = g.eta_s, self.E_emi
        if name == 'v':
            base = 1j * self.kv[:, comp][:, None] * np.atleast_1d(tv(self.kmag, T, me))[None, :] if where == "ray" else 1j * self.kv[:, comp] * tv(self.kmag, T, me)
        else:
            a = self.amp(name, T, me)
            base = np.broadcast_to(a[None, :], (self.N, g.M)).copy() if where == 'ray' else np.full(self.N, a)
        fac = (1j * self.kv[:, 0])**mx * (1j * self.kv[:, 1])**my * (1j * self.kv[:, 2])**mz
        res = base * (fac[:, None] if where == 'ray' else fac) * phase
        self.cache[key] = res
        return res
    # observer-event shifts (direction independent, x_o = 0)
    def deta_o1(self):
        g = self.g
        return -np.sum(g.ww * g.tw * tphi(self.kmag, g.tw)) / g.eta_o
    def dxo1(self, c):
        g = self.g
        return 1j * self.kv[:, c] * np.sum(g.ww * tv(self.kmag, g.tw))
    def deta_e1(self):
        return self.g.eta_s * tdt(self.kmag, self.g.eta_s) * self.E_emi
    def sym(self, name, comp):
        if name == 'deta_o1': return np.full(self.N, self.deta_o1())
        if name == 'dxo1': return self.dxo1(comp)
        if name == 'deta_e1': return self.deta_e1()
    # ray functionals (n = z-hat frame): k1_perp = i k_perp [T_V(eta_o) - 2 int_eta^eta_o T_Phi e^{i k_z chi'}], k1_z = -2 T_Phi e^{..}
    def k1(self, c):
        key = ('k1', c)
        if key in self.cache: return self.cache[key]
        g = self.g
        if c < 2:
            res = 1j * self.kv[:, c][:, None] * (tv(self.kmag, g.eta_o) - 2 * g.cum(tphi(self.kmag, g.t)[None, :] * self.E_ray))
        else:
            res = -2 * tphi(self.kmag, g.t)[None, :] * self.E_ray
        self.cache[key] = res
        return res
    def Dx1(self, c):
        key = ('Dx1', c)
        if key in self.cache: return self.cache[key]
        g = self.g
        res = self.dxo1(c)[:, None] - g.cum(self.k1(c))
        if c == 2: res = res + self.deta_o1()
        self.cache[key] = res
        return res
    def fun(self, name, c, when):
        arr = self.k1(c) if name == 'k1' else self.Dx1(c)
        return arr[:, 0] if when == 'eta_star' else arr

# ------------------------------------------------------------------ second-order fields for a momentum pair
class Pair:
    """second-order fields X2(K; k1,k2) on the grid, K = k1 + k2; zero_mode=True: exact K = 0 (plane wave) prescription."""
    def __init__(self, geom, k1vec, k2vec, nv, zero_mode=False):
        self.g = g = geom
        k1vec = np.asarray(k1vec, float); k2vec = np.asarray(k2vec, float)
        self.q, self.p = np.linalg.norm(k1vec), np.linalg.norm(k2vec)
        Kvec = k1vec + k2vec
        self.K = np.linalg.norm(Kvec)
        self.zero_mode = zero_mode or self.K < 1e-13
        e1, e2 = transverse_basis(nv)
        self.Kv = np.stack([e1 @ Kvec, e2 @ Kvec, nv @ Kvec], axis=1)
        self.k1v = np.stack([e1 @ k1vec, e2 @ k1vec, nv @ k1vec], axis=1)
        self.k2v = np.stack([e1 @ k2vec, e2 @ k2vec, nv @ k2vec], axis=1)
        self.N = self.Kv.shape[0]
        self.E_ray = np.exp(1j * self.Kv[:, 2][:, None] * g.chi[None, :])
        self.E_emi = np.exp(1j * self.Kv[:, 2] * (g.eta_o - g.eta_s))
        self.cache = {}
        self._scalars()
    def _scalars(self):
        g, K, q, p = self.g, self.K, self.q, self.p
        t = g.t
        if self.zero_mode:
            # K = 0: scalar and vector zero modes are monopole/dipole (dropped); tensor zero mode from the ODE h''+2h'/eta = -2 S^TL
            Kh = KER['K_h'](1e-30 if K == 0 else K, q, p, g.gs)
            self.H = np.sum(g.gw * G_h0(t[:, None], g.gs) * Kh, axis=1)
            self.dH = np.sum(g.gw * dG_h0(t[:, None], g.gs) * Kh, axis=1)
            z = np.zeros(g.M)
            self.Phi2 = self.dPhi2 = self.Q2 = self.Psi2 = self.delta2 = self.dT2 = self.V2 = z
            self.KvT = self.Kom = z
            return
        SPhi = KER['K_SPhi'](K, q, p, g.gs)                       # (M,Ng)
        self.Phi2 = np.sum(g.gw * G_Phi(K, t[:, None], g.gs) * SPhi, axis=1) + KER['K_Phi2_0'](K, q, p, t) * ghat(K, t)
        self.dPhi2 = np.sum(g.gw * dG_Phi(K, t[:, None], g.gs) * SPhi, axis=1) + KER['K_Phi2_0'](K, q, p, t) * ghat(K, t, 1)
        self.Q2 = KER['K_Q'](K, q, p, t)
        self.Psi2 = self.Phi2 - self.Q2
        self.delta2 = (t**2 / 3) * KER['K_S00'](K, q, p, t) - (2 / 3) * t**2 * K**2 * self.Phi2 - 2 * t * self.dPhi2 - 2 * self.Psi2
        self.dT2 = self.delta2 / 4 - (3 / 16) * tdelta(q, t) * tdelta(p, t)
        v2z = (t / 4) * (-t * KER['K_S0z'](K, q, p, t) - 2j * t * K * self.dPhi2 - 2j * K * self.Psi2)
        self.V2 = v2z / (1j * K)
        self.KvT = KER['K_vT'](K, q, p, t); self.Kom = KER['K_om'](K, q, p, t)
        Kh = KER['K_h'](K, q, p, g.gs)
        self.H = np.sum(g.gw * G_h(K, t[:, None], g.gs) * Kh, axis=1)
        self.dH = np.sum(g.gw * dG_h(K, t[:, None], g.gs) * Kh, axis=1)
    # angular structures (N,...)
    def Pk1(self):
        if self.zero_mode: return self.k1v
        Kk1 = (self.K**2 + self.q**2 - self.p**2) / 2
        return self.k1v - self.Kv * Kk1 / self.K**2
    def Pk2(self):
        if self.zero_mode: return self.k2v
        Kk2 = (self.K**2 + self.p**2 - self.q**2) / 2
        return self.k2v - self.Kv * Kk2 / self.K**2
    def LamQQ(self, i, j):
        """Lambda_ij,kl(K) k1_k k1_l (N,) ; zero mode: traceless projector (k1_i k1_j - delta_ij q^2/3)."""
        if self.zero_mode:
            return self.k1v[:, i] * self.k1v[:, j] - (self.q**2 / 3 if i == j else 0.0)
        Pk1 = self.Pk1()
        Pij = (1.0 if i == j else 0.0) - self.Kv[:, i] * self.Kv[:, j] / self.K**2
        qPq = self.q**2 - ((self.K**2 + self.q**2 - self.p**2) / 2)**2 / self.K**2
        return Pk1[:, i] * Pk1[:, j] - 0.5 * Pij * qPq
    def profile(self, name, comp, me):
        """time profile (N,M) of the second-order field component before spatial derivatives/phase."""
        M = self.g.M
        bc = lambda a: np.broadcast_to(a[None, :], (self.N, M))
        if name == 'Phi': return bc(self.Phi2 if me == 0 else self.dPhi2)
        if name == 'Psi': return bc(self.Psi2 if me == 0 else self.dPhi2 - np.gradient(self.Q2, self.g.t))  # (Q2' not needed in L2)
        if name == 'dT': assert me == 0; return bc(self.dT2)
        if name == 'v':
            assert me == 0
            vL = 1j * self.Kv[:, comp][:, None] * self.V2[None, :]
            vT = 0.5 * (self.Pk1()[:, comp][:, None] * self.KvT[None, :] + self.Pk2()[:, comp][:, None] * KER['K_vT'](self.K, self.p, self.q, self.g.t)[None, :]) if not self.zero_mode else 0.0
            return vL + vT
        if name == 'om':
            assert me == 0
            if self.zero_mode: return np.zeros((self.N, M))
            return 0.5 * (self.Pk1()[:, comp][:, None] * self.Kom[None, :] + self.Pk2()[:, comp][:, None] * KER['K_om'](self.K, self.p, self.q, self.g.t)[None, :])
        if name == 'h':
            i, j = comp
            return self.LamQQ(i, j)[:, None] * (self.H if me == 0 else self.dH)[None, :]
        raise ValueError(name)
    def split_profiles(self, name, comp, me):
        """sub-channel decomposition of v (longitudinal / transverse) for the channel bookkeeping."""
        if name != 'v': return {name: self.profile(name, comp, me)}
        vL = 1j * self.Kv[:, comp][:, None] * self.V2[None, :]
        if self.zero_mode: return {'V2': vL}
        vT = 0.5 * (self.Pk1()[:, comp][:, None] * self.KvT[None, :] + self.Pk2()[:, comp][:, None] * KER['K_vT'](self.K, self.p, self.q, self.g.t)[None, :])
        return {'V2': vL, 'v2T': vT}
    def field(self, name, comp, der, where, sub=False):
        me, mx, my, mz = der
        fac = (1j * self.Kv[:, 0])**mx * (1j * self.Kv[:, 1])**my * (1j * self.Kv[:, 2])**mz
        profs = self.split_profiles(name, comp, me) if sub else {name: self.profile(name, comp, me)}
        out = {}
        for key, prof in profs.items():
            if where == 'ray': out[key] = prof * fac[:, None] * self.E_ray
            elif where == 'obs': out[key] = prof[:, -1] * fac
            elif where == 'emi': out[key] = prof[:, 0] * fac * self.E_emi
        return out

# ------------------------------------------------------------------ term evaluation
def eval_first(slot, desc, where):
    kind = desc[0]
    if kind == 'field': return slot.field(desc[1], desc[2], desc[4], where)
    if kind == 'fun': return slot.fun(desc[1], desc[2], desc[3])
    if kind == 'sym': return slot.sym(desc[1], desc[2])
    raise ValueError(desc)

def eval_termlist(terms, geom, slots, pair, where_of):
    """terms: list of C2 term dicts; slots: (Slot_a, Slot_b); pair: Pair; returns dict channel -> (N,) kernel (coefficient of eps^2 i.e. L2/2)."""
    sa, sb = slots
    N = sa.N
    out = {}
    for term in terms:
        where = term['where']
        for coeff, descs in parse_term(term['expr']):
            num = [d for d in descs if d[0] in ('eta', 'eta_star')]
            objs = [d for d in descs if d[0] not in ('eta', 'eta_star')]
            orders = [order_of(d) for d in objs]
            cfac = coeff
            for d in num:
                if d[0] == 'eta_star': cfac = cfac * geom.eta_s
            etafac = sum(1 for d in num if d[0] == 'eta')
            def finish(arr):
                if etafac:
                    T = geom.t if where == 'ray' else (geom.eta_o if where == 'obs' else geom.eta_s)
                    arr = arr * (T**etafac if np.ndim(arr) == 1 else (T**etafac)[None, :])
                arr = arr * cfac
                return geom.integ(arr) if where == 'ray' else arr
            if sorted(orders) == [1, 1]:
                A, B = objs
                val = 0.5 * (eval_first(sa, A, where) * eval_first(sb, B, where) + eval_first(sb, A, where) * eval_first(sa, B, where))
                name = term['name']
                out[name] = out.get(name, 0) + finish(val)
            elif orders == [2]:
                d = objs[0]
                for sub, arr in pair.field(d[1], d[2], d[4], where, sub=True).items():
                    name = '[%s] 2nd:%s' % (where, sub)
                    out[name] = out.get(name, 0) + finish(arr)
            else:
                raise ValueError((orders, objs))
    return out

def Delta_l(l, kmag, geom):
    """type-A transfer function of m1 = L1 (standard form), with observer terms (l = 0, 1 only).  Real."""
    g = geom; k = kmag; chi_s = g.eta_o - g.eta_s
    val = (tdt(k, g.eta_s) + tphi(k, g.eta_s)) * spherical_jn(l, k * chi_s) - k * tv(k, g.eta_s) * spherical_jn(l, k * chi_s, derivative=True) \
        + g.integ(2 * tphi(k, g.t, 1) * spherical_jn(l, k * g.chi))
    if l == 0: val -= (tdt(k, g.eta_o) + tphi(k, g.eta_o))
    if l == 1: val += k * tv(k, g.eta_o) / 3        # +(n.v)_o: (i n.k) -> j_1prime(0) = 1/3
    return val

class MapKernels:
    """K_1 and K_2 (both C2 forms) for one momentum pair and an array of directions."""
    def __init__(self, geom, k1vec, k2vec, nv, form='direct', zero_mode=False):
        self.g = geom; self.form = form
        nv = np.asarray(nv, float)
        self.nv = nv
        self.sa = Slot(geom, k1vec, nv); self.sb = Slot(geom, k2vec, nv)
        self.pair = Pair(geom, k1vec, k2vec, nv, zero_mode=zero_mode)
        self.L1list = LT['L_direct_split'][1] if form == 'direct' else LT['L_split'][1]
        self.L2list = LT['L_direct_split'][2] if form == 'direct' else LT['L_split'][2]
    def K1(self, slot):
        """K_1(k;n) from the L1 term list (first-order terms, momentum of the slot)."""
        N = slot.N; tot = np.zeros(N, complex)
        for term in self.L1list:
            where = term['where']
            for coeff, descs in parse_term(term['expr']):
                objs = [d for d in descs if d[0] not in ('eta', 'eta_star')]
                assert len(objs) == 1 and order_of(objs[0]) == 1, objs
                arr = eval_first(slot, objs[0], where) * coeff
                tot += self.g.integ(arr) if where == 'ray' else arr
        return tot
    def channels(self):
        """dict: channel -> kernel contribution to K_2 (already including the factor 2 of K_2 = 2 K_{L2/2} + ...)."""
        ch = eval_termlist(self.L2list, self.g, (self.sa, self.sb), self.pair, None)
        ch = {k_: 2 * v_ for k_, v_ in ch.items()}
        K1a, K1b = self.K1(self.sa), self.K1(self.sb)
        D0a, D0b = Delta_l(0, self.sa.kmag, self.g), Delta_l(0, self.sb.kmag, self.g)
        ch['L1*L1'] = K1a * K1b
        ch['-2<m1>m1'] = -(D0a * K1b + D0b * K1a)
        self._K1 = (K1a, K1b)
        return ch
    def K2(self):
        ch = self.channels()
        return sum(ch.values()), ch

if __name__ == '__main__':
    g = Geometry()
    rng = np.random.default_rng(1)
    nv = rng.normal(size=(5, 3)); nv /= np.linalg.norm(nv, axis=1)[:, None]
    k1 = np.array([0.7, -0.2, 0.4]); k2 = np.array([-0.3, 0.9, 0.5])
    # K_1 check: engine (direct-form L1 list) vs type-A Delta_l expansion
    mk = MapKernels(g, k1, k2, nv)
    K1 = mk.K1(mk.sa)
    from scipy.special import sph_harm_y
    kk = np.linalg.norm(k1); kh = k1 / kk
    th_n = np.arccos(nv[:, 2]); ph_n = np.arctan2(nv[:, 1], nv[:, 0]); th_k = np.arccos(kh[2]); ph_k = np.arctan2(kh[1], kh[0])
    rec = np.zeros(len(nv), complex)
    for l in range(0, 40):
        Dl = Delta_l(l, kk, g)
        for m in range(-l, l + 1):
            rec += 4 * np.pi * 1j**l * Dl * sph_harm_y(l, m, th_n, ph_n) * np.conj(sph_harm_y(l, m, th_k, ph_k))
    print('K_1 engine:', K1); print('K_1 type-A :', rec); print('max diff', np.max(np.abs(K1 - rec)))
    for form in ('direct', 'shear'):
        mk = MapKernels(g, k1, k2, nv, form=form)
        K2, ch = mk.K2()
        print(form, 'K2 =', K2)
        for k_, v_ in ch.items(): print('   %-60s %s' % (k_, v_[:2]))
