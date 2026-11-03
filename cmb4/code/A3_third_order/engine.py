"""engine.py (A3) -- polarisation engine for the map kernels K_1, K_2, K_3 evaluated from C2's line-of-sight term lists.

ALGEBRA.  zeta = a e^{i k_0.x} + b e^{i k_1.x} + c e^{i k_2.x} with nilpotent a,b,c (a^2=b^2=c^2=0).  Every quantity is an element
of R[a,b,c]/(a^2,b^2,c^2): a dict {mask: array}, mask = bit set of slots (1=a, 2=b, 4=c).  The abc (mask 7) component of the
cubic map Th3 is the polarisation coefficient 6 K3sym(k0,k1,k2;n); for (k0,k1,k2) = (k,q,-q):
     kappa(k;n) := 3 int d^3q/(2pi)^3 P(q) K3sym(k,q,-q;n) = (1/2) int d^3q/(2pi)^3 P(q) coef_abc(k,q,-q;n)
     = <d Th3/d zeta_k> (Gaussian response), and C^(13)_l = (2/pi) int k^2 dk P(k) Delta_l(k) Delta^kappa_l(k)*.
Second-order objects of a slot pair carry AB-COEFFICIENTS (K_X(q1,q2)+K_X(q2,q1)); third-order fields carry the abc coefficient.

GEOMETRY (C2 convention): background ray xbar(eta) = n chi, chi = eta_o - eta, observer at x = 0; fields on the ray carry
e^{i (k.n) chi}.  For every direction n an orthonormal frame (e1, e2, n) is used and the C2 'split' term lists (n = z-hat) are
evaluated with frame components; the second-order direction map k2_init (canonical Gram-Schmidt triad = global-coordinate
dependent) is evaluated in GLOBAL components with general n and rotated into the frame.
RAY FUNCTIONALS (re-derived; NOTE C2's `definitions` dict has the wrong sign of the integrals, los.py itself uses D k = G):
   k1(eta) = k1_init - int_eta^eta_o G1 ;  Dx1 = dxo1 + n deta_o1 - int_eta^eta_o k1
   k2/2 = k2o/2 - int_eta^eta_o [G2/2 + Dx1.grad G1] ;  Dx2/2 = dxo2/2 + n deta_o2/2 - int_eta^eta_o k2/2 - deta_o1 k1(eta_o)
OBSERVER (re-derived from SPEC 6, eta_i -> 0):  X1(eta) = int_0^eta v1(eta',0),  X2/2 = int_0^eta [v2/2 + X1.grad v1];
   deta_o1 = -(1/eta_o) int_0^eta_o eta Psi1 ;  deta_o2/2 = -(1/eta_o)[deta_o1^2/2 + deta_o1 eta_o Psi1(eta_o) + int_0^eta_o eta f2/2],
   f2/2 = Psi2/2 - Psi1^2/2 - v1^2/2 + X1.grad Psi1 ;  dxo1 = X1(eta_o),  dxo2/2 = X2(eta_o)/2 + deta_o1 v1(eta_o).
"""
import os, re, pickle
import numpy as np
import sympy as sp
from numpy.polynomial import chebyshev as Ch
import tf
from grid import TGrid, pair_profiles, zero_profiles

HERE = os.path.dirname(os.path.abspath(__file__))
C2DIR = os.path.join(HERE, '..', 'C2_observable')
LT = pickle.load(open(os.path.join(C2DIR, 'L_terms.pkl'), 'rb'))
import sys; sys.path.append(C2DIR)
RC = pickle.load(open(os.path.join(C2DIR, 'ray_cache.pkl'), 'rb'))
ETA, XS, YS, ZS = sp.symbols('eta x y z')
NSYM = RC['n']
NZ = {NSYM[0]: 0, NSYM[1]: 0, NSYM[2]: 1}

# ------------------------------------------------------------------ algebra
def amul(A, B):
    out = {}
    for ma, va in A.items():
        for mb, vb in B.items():
            if ma & mb: continue
            m = ma | mb
            out[m] = out[m] + va * vb if m in out else va * vb
    return out
def aadd(A, B, s=1.0):
    out = dict(A)
    for m, v in B.items():
        out[m] = out[m] + s * v if m in out else s * v
    return out
def ascale(A, s):
    return {m: s * v for m, v in A.items()}
def amap(A, f):
    return {m: f(v) for m, v in A.items()}
def apart(A, maxorder):
    return {m: v for m, v in A.items() if bin(m).count('1') <= maxorder}

# ------------------------------------------------------------------ parsing of C2 monomials
_re_scalar = re.compile(r'^(Psi|Phi|dT)(\d)$')
_re_vec = re.compile(r'^(v|om)_(\d)_(\d)$')
_re_h = re.compile(r'^h_(\d)(\d)_(\d)$')
_re_fun = re.compile(r'^(k|Dx)(\d)_(\d)$')
def parse_factor(f):
    der = [0, 0, 0, 0]
    if isinstance(f, sp.Derivative):
        for v_, c_ in f.variable_count:
            der[(ETA, XS, YS, ZS).index(v_)] += int(c_)
        f = f.expr
    if isinstance(f, sp.core.function.AppliedUndef):
        nm = f.func.__name__
        m = _re_scalar.match(nm)
        if m: return ('f', m.group(1), None, int(m.group(2)), tuple(der))
        m = _re_vec.match(nm)
        if m: return ('f', m.group(1), int(m.group(2)) - 1, int(m.group(3)), tuple(der))
        m = _re_h.match(nm)
        if m:
            i, j = int(m.group(1)) - 1, int(m.group(2)) - 1
            return ('f', 'h', (min(i, j), max(i, j)), int(m.group(3)), tuple(der))
        m = _re_fun.match(nm)
        if m:
            assert der == [0, 0, 0, 0], f
            when = 'eta' if f.args[0] == ETA else 'eta_star'
            return ('fun', m.group(1), int(m.group(2)), int(m.group(3)) - 1, when)
        raise ValueError(nm)
    if isinstance(f, sp.Symbol):
        s = f.name
        mm = re.match(r'^(deta_o|deta_e)(\d)$', s)
        if mm: return ('sym', mm.group(1), int(mm.group(2)), None)
        mm = re.match(r'^dxo(\d)_(\d)$', s)
        if mm: return ('sym', 'dxo', int(mm.group(1)), int(mm.group(2)) - 1)
        if s == 'eta': return ('eta',)
        if s == 'eta_star': return ('eta_star',)
        if s in ('n1', 'n2', 'n3'): return ('n', int(s[1]) - 1)
        raise ValueError(s)
    return None

_parse_cache = {}
def parse_expr(expr):
    key = sp.srepr(expr)
    if key in _parse_cache: return _parse_cache[key]
    out = []
    for mono in sp.Add.make_args(sp.expand(expr)):
        coeff = 1.0 + 0j; descs = []
        for f in sp.Mul.make_args(mono):
            base, ex = (f.base, int(f.exp)) if isinstance(f, sp.Pow) else (f, 1)
            d = parse_factor(base)
            if d is None: coeff *= complex(base**ex)
            else: descs += [d] * ex
        out.append((coeff, descs))
    _parse_cache[key] = out
    return out

# ------------------------------------------------------------------ ray grid (Lobatto on [eta_*, eta_o])
class RayGrid:
    def __init__(self, eta_s=2.0, eta_o=8.0, M=64):
        self.eta_s, self.eta_o, self.M = float(eta_s), float(eta_o), M
        xg = -np.cos(np.pi * np.arange(M) / (M - 1))
        self.t = eta_s + (xg + 1) * (eta_o - eta_s) / 2
        self.chi = eta_o - self.t
        V = Ch.chebvander(xg, M - 1); Cinv = np.linalg.inv(V)
        B = np.zeros((M, M))
        for j in range(M):
            e = np.zeros(M); e[j] = 1.0
            ci = Ch.chebint(e)
            B[:, j] = Ch.chebval(1.0, ci) - Ch.chebval(xg, ci)
        self.Cup = (B @ Cinv) * (eta_o - eta_s) / 2        # (Cup f)_i = int_{t_i}^{eta_o} f
        self.w = self.Cup[0].copy()
    def cum(self, f):
        return f @ self.Cup.T
    def integ(self, f):
        return f @ self.w

def frames(nv):
    """orthonormal (e1, e2, n) with e1 x e2 = n; returns E (Nn,3,3), rows e1,e2,n."""
    nv = np.asarray(nv, float)
    a = np.zeros_like(nv); a[:, 1] = 1.0
    bad = np.abs(nv[:, 1]) > 0.9
    a[bad] = 0; a[bad, 0] = 1.0
    e2 = a - np.sum(a * nv, 1)[:, None] * nv; e2 /= np.linalg.norm(e2, axis=1)[:, None]
    e1 = np.cross(e2, nv)
    return np.stack([e1, e2, nv], axis=1)

# ------------------------------------------------------------------ the evaluation context
class Ctx:
    """slots: list of 3 arrays (B,3) of global momenta (slot 0 may be (1,3)); nv: (Nn,3) directions.
    third: optional dict of contracted third-order field profiles on the global grid (attached to slot 0, mask 7)."""
    def __init__(self, slots, nv, g=None, ray=None, third=None, tadshift=None):
        self.g = g or TGrid()
        self.ray = ray or RayGrid(eta_o=self.g.E)
        self.nv = np.asarray(nv, float); self.Nn = len(self.nv)
        self.E = frames(self.nv)
        self.slots = [np.atleast_2d(np.asarray(s, float)) for s in slots]
        self.B = max(s.shape[0] for s in self.slots)
        self.slots = [np.broadcast_to(s, (self.B, 3)) for s in self.slots]
        self.nslot = len(self.slots)
        self.third = third
        self.tadshift = tadshift
        r, g = self.ray, self.g
        self.I_ray = g.interp_matrix(r.t); self.I_obs = g.interp_matrix([r.eta_o]); self.I_emi = g.interp_matrix([r.eta_s])
        self.cache = {}
        self.sl = []
        for s in range(self.nslot):
            kg = self.slots[s]
            d = dict(kg=kg, kmag=np.linalg.norm(kg, axis=1))
            d['kv'] = np.einsum('nab,Bb->Bna', self.E, kg)         # frame components (B,Nn,3)
            self.sl.append(d)
        self.pairs = {}
        for s1 in range(self.nslot):
            for s2 in range(s1 + 1, self.nslot):
                self._make_pair(s1, s2)
    # ---------------- helpers
    def phase(self, kz, where):
        if where == 'ray': return np.exp(1j * kz[..., None] * self.ray.chi)
        if where == 'emi': return np.exp(1j * kz * (self.ray.eta_o - self.ray.eta_s))
        return np.ones_like(kz) + 0j
    def times(self, where):
        return {'ray': self.ray.t, 'obs': self.ray.eta_o, 'emi': self.ray.eta_s, 'wl': self.g.t}[where]
    def from_global(self, prof, where):
        """global-grid node values (B, N) -> (B, M) / (B,) / (B, N)."""
        if where == 'ray': return prof @ self.I_ray.T
        if where == 'obs': return (prof @ self.I_obs.T)[:, 0]
        if where == 'emi': return (prof @ self.I_emi.T)[:, 0]
        return prof
    def bshape(self, arr_bn, prof, where):
        """combine angular factor (B,Nn) with time profile (B,M)/(B,)/(B,N)."""
        if where in ('ray', 'wl'): return arr_bn[..., None] * prof[:, None, :]
        return arr_bn * prof[:, None]
    # ---------------- pairs
    def _make_pair(self, s1, s2):
        mask = (1 << s1) | (1 << s2)
        k1, k2 = self.slots[s1], self.slots[s2]
        Kg = k1 + k2
        Km = np.linalg.norm(Kg, axis=1)
        q1m, q2m = np.linalg.norm(k1, axis=1), np.linalg.norm(k2, axis=1)
        zero = Km < 1e-12
        P = dict(s=(s1, s2), Kg=Kg, Km=Km, zero=zero)
        prof = None
        if np.any(~zero):
            trip = np.round(np.stack([np.where(zero, 1.0, Km), q1m, q2m], 1), 13)
            u, inv = np.unique(trip, axis=0, return_inverse=True)
            pp = pair_profiles(self.g, u[:, 0], u[:, 1], u[:, 2])
            prof = {nm_: v[inv.ravel()] for nm_, v in pp.items()}
        if np.any(zero):
            zp = zero_profiles(self.g, q1m)
            if self.tadshift is not None:
                for nm_, val in self.tadshift(self.g.t, q1m).items(): zp[nm_] = zp[nm_] + val
            if prof is None: prof = zp
            else:
                for nm_ in prof: prof[nm_] = np.where(zero[:, None], zp[nm_], prof[nm_])
        P['prof'] = prof
        P['Kv'] = np.einsum('nab,Bb->Bna', self.E, Kg)
        P['k1v'] = self.sl[s1]['kv']; P['k2v'] = self.sl[s2]['kv']
        Ks2 = np.where(zero, 1.0, Km**2)[:, None, None]
        P['Pq1'] = np.where(zero[:, None, None], 0.0, P['k1v'] - P['Kv'] * np.sum(P['Kv'] * P['k1v'], -1, keepdims=True) / Ks2)
        P['Pq2'] = np.where(zero[:, None, None], 0.0, P['k2v'] - P['Kv'] * np.sum(P['Kv'] * P['k2v'], -1, keepdims=True) / Ks2)
        # global-basis versions (for k2_init)
        kg1 = np.broadcast_to(k1[:, None, :], (self.B, self.Nn, 3)); kg2 = np.broadcast_to(k2[:, None, :], (self.B, self.Nn, 3))
        KgB = kg1 + kg2
        P['Kv_g'], P['k1v_g'] = KgB, kg1
        P['Pq1_g'] = np.where(zero[:, None, None], 0.0, kg1 - KgB * np.sum(KgB * kg1, -1, keepdims=True) / Ks2)
        P['Pq2_g'] = np.where(zero[:, None, None], 0.0, kg2 - KgB * np.sum(KgB * kg2, -1, keepdims=True) / Ks2)
        self.pairs[mask] = P
    # ---------------- field evaluation
    def field(self, desc, where, basis='frame'):
        key = ('f', desc, where, basis)
        if key in self.cache: return self.cache[key]
        _, name, comp, order, der = desc
        me, mx, my, mz = der
        T = self.times(where)
        out = {}
        if order == 1:
            for s in range(self.nslot):
                d = self.sl[s]
                kv = d['kv'] if basis == 'frame' else np.broadcast_to(d['kg'][:, None, :], (self.B, self.Nn, 3))
                kmag = d['kmag'][0]
                assert np.allclose(d['kmag'], kmag)
                if name == 'v':
                    ang = 1j * kv[..., comp]; prof = tf.tv(kmag, T, me)
                elif name in ('Phi', 'Psi', 'dT', 'delta'):
                    ang = np.ones(kv.shape[:2], complex); prof = tf.first(name, kmag, T, me)
                else:
                    continue   # no first-order om, h
                ang = ang * (1j * kv[..., 0])**mx * (1j * kv[..., 1])**my * (1j * kv[..., 2])**mz
                ph = self.phase(d['kv'][..., 2], where)
                if where in ('ray', 'wl'):
                    out[1 << s] = ang[..., None] * np.asarray(prof)[None, None, :] * (ph if where == 'ray' else 1.0)
                else:
                    out[1 << s] = ang * prof * ph
        elif order == 2:
            for mask, P in self.pairs.items():
                arr = self.pair_field(P, name, comp, der, where, basis)
                if arr is not None: out[mask] = arr
        elif order == 3:
            if self.third is not None and self.nslot == 3:
                arr = self.third_field(name, comp, der, where)
                if arr is not None: out[7] = arr
            if name == 'dT' and der == (0, 0, 0, 0):
                # delta_T3 = delta3/4 - (9/16) delta1 delta2 + (21/64) delta1^3  (from (1+delta)^(1/4)); delta3/4 is in 'dT3lin'
                d1 = self.field(('f', 'delta', None, 1, der), where); d2 = self.field(('f', 'delta', None, 2, der), where)
                out = aadd(out, aadd(ascale(amul(d1, d2), -9 / 16), ascale(amul(amul(d1, d1), d1), 21 / 64)))
        self.cache[key] = out
        return out
    def pair_field(self, P, name, comp, der, where, basis):
        me, mx, my, mz = der
        pr = P['prof']
        g = lambda nm: self.from_global(pr[nm], where)
        Kv = P['Kv'] if basis == 'frame' else P['Kv_g']
        Pq1 = P['Pq1'] if basis == 'frame' else P['Pq1_g']
        Pq2 = P['Pq2'] if basis == 'frame' else P['Pq2_g']
        k1v = P['k1v'] if basis == 'frame' else P['k1v_g']
        dd = {0: '', 1: 'd', 2: 'dd'}[me]
        if name in ('Phi', 'Psi', 'dT', 'delta'):
            nm = {'Phi': {0: 'Phi', 1: 'dPhi', 2: 'ddPhi'}, 'Psi': {0: 'Psi', 1: 'dPsi'}, 'dT': {0: 'dT', 1: 'ddT'}, 'delta': {0: 'delta'}}[name][me]
            arr = self.bshape(np.ones(Kv.shape[:2], complex), g(nm), where)
        elif name == 'v':
            arr = self.bshape(1j * Kv[..., comp], g(dd + 'V'), where) + self.bshape(Pq1[..., comp] + 0j, g(dd + 'vT12'), where) \
                + self.bshape(Pq2[..., comp] + 0j, g(dd + 'vT21'), where)
        elif name == 'om':
            arr = self.bshape(Pq1[..., comp] + 0j, g(dd + 'om12'), where) + self.bshape(Pq2[..., comp] + 0j, g(dd + 'om21'), where)
        elif name == 'h':
            i, j = comp
            Ks2 = np.where(P['zero'], 1.0, P['Km']**2)[:, None]
            Pij = (1.0 if i == j else 0.0) - Kv[..., i] * Kv[..., j] / Ks2
            q1sq = np.sum(k1v * k1v, -1); Kq1 = np.sum(Kv * k1v, -1)
            qPq = q1sq - Kq1**2 / Ks2
            Lam = Pq1[..., i] * Pq1[..., j] - 0.5 * Pij * qPq
            Lam = np.where(P['zero'][:, None], 0.0, Lam)
            arr = self.bshape(Lam + 0j, g(dd + 'H'), where)
        else:
            raise ValueError(name)
        if mx or my or mz:
            if np.all(P['zero']): return None
            arr = arr * self._expand((1j * Kv[..., 0])**mx * (1j * Kv[..., 1])**my * (1j * Kv[..., 2])**mz, where)
        if where in ('ray', 'emi'):
            arr = arr * self.phase(P['Kv'][..., 2], where)
        return arr
    def _expand(self, a, where):
        return a[..., None] if where in ('ray', 'wl') else a
    def third_field(self, name, comp, der, where):
        """contracted third-order fields (q-averaged, x2 for the abc coefficient) attached to slot 0."""
        me, mx, my, mz = der
        th = self.third
        d = self.sl[0]; kv = d['kv']; kmag = d['kmag'][0]
        dd = {0: '', 1: 'd'}[me]
        if name in ('Phi', 'Psi', 'dT'):
            key = dd + {'Phi': 'Phi', 'Psi': 'Psi', 'dT': 'dT3lin'}[name]
            ang = np.ones(kv.shape[:2], complex)
        elif name == 'v':
            key = dd + 'V'; ang = 1j * kv[..., comp]
        else:
            return None       # om3, h3: contracted transverse/TT parts vanish by isotropy
        prof = self.from_global(np.atleast_2d(th[key]), where)
        ang = ang * (1j * kv[..., 0])**mx * (1j * kv[..., 1])**my * (1j * kv[..., 2])**mz
        arr = self.bshape(ang, prof, where)
        if where in ('ray', 'emi'): arr = arr * self.phase(kv[..., 2], where)
        return arr
    # ---------------- generic expression evaluation
    def evaluate(self, expr, where, basis='frame', maxorder=3):
        tot = {}
        T = self.times(where)
        for coeff, descs in parse_expr(expr):
            el = {0: coeff}
            for d in descs:
                if d[0] == 'eta': el = {m: v * T for m, v in el.items()}; continue
                if d[0] == 'eta_star': el = ascale(el, self.ray.eta_s); continue
                if d[0] == 'n': el = ascale(el, self.nv[:, d[1]]) if basis == 'global' else ascale(el, 1.0 if d[1] == 2 else 0.0); continue
                el = apart(amul(el, self.obj(d, where, basis)), maxorder)
                if not el: break
            tot = aadd(tot, el)
        return tot
    def obj(self, d, where, basis='frame'):
        if d[0] == 'f': return self.field(d, where, basis)
        if d[0] == 'fun':
            _, kind, order, comp, when = d
            arr = self.functional(kind, order, comp)
            if where == 'emi' or when == 'eta_star': return amap(arr, lambda v: v[..., 0])
            assert where == 'ray', (d, where)
            return arr
        if d[0] == 'sym':
            _, name, order, comp = d
            return self.symbol(name, order, comp)
        raise ValueError(d)
    # ---------------- observer / emission symbols
    def symbol(self, name, order, comp):
        key = ('sym', name, order, comp)
        if key in self.cache: return self.cache[key]
        g = self.g
        if name == 'deta_o' and order == 1:
            Psi = self.field(('f', 'Psi', None, 1, (0, 0, 0, 0)), 'wl')
            res = amap(Psi, lambda v: -g.integ(v * g.t) / g.E)
        elif name == 'dxo' and order == 1:
            res = amap(self.field(('f', 'v', comp, 1, (0, 0, 0, 0)), 'wl'), lambda v: g.integ(v))
        elif name == 'deta_e' and order == 1:
            res = ascale(self.field(('f', 'dT', None, 1, (0, 0, 0, 0)), 'emi'), self.ray.eta_s)
        elif name == 'deta_o' and order == 2:
            d1 = self.symbol('deta_o', 1, None)
            Psi1o = self.field(('f', 'Psi', None, 1, (0, 0, 0, 0)), 'obs')
            f2 = self.worldline_f2half()
            inner = amap(f2, lambda v: g.integ(v * g.t))
            br = aadd(aadd(ascale(amul(d1, d1), 0.5), ascale(amul(d1, Psi1o), g.E)), inner)
            res = ascale(apart(br, 2), -2.0 / g.E)        # deta_o2 itself (= 2 x deta_o2/2)
        elif name == 'dxo' and order == 2:
            d1 = self.symbol('deta_o', 1, None)
            v1o = self.field(('f', 'v', comp, 1, (0, 0, 0, 0)), 'obs')
            integrand = aadd(ascale(self.field(('f', 'v', comp, 2, (0, 0, 0, 0)), 'wl'), 0.5), self.X1grad(('f', 'v', comp, 1)))
            res = ascale(apart(aadd(amap(integrand, lambda v: g.integ(v)), amul(d1, v1o)), 2), 2.0)
        elif name == 'deta_e' and order == 2:
            ex = LT['emission'][sp.Symbol('deta_e2')].subs(NZ)
            res = apart(self.evaluate(ex, 'emi'), 2)
        else:
            raise ValueError(key)
        self.cache[key] = res
        return res
    def X1(self, comp):
        """observer worldline displacement X1_c(eta) = int_0^eta v1_c(eta',0) on global nodes (B,Nn,N)."""
        key = ('X1', comp)
        if key in self.cache: return self.cache[key]
        g = self.g
        v = self.field(('f', 'v', comp, 1, (0, 0, 0, 0)), 'wl')
        # cumulative integral from 0 via Chebyshev coefficients
        VN = Ch.chebvander(g.xg, g.N)
        def cum0(f):
            c = f @ g.Vinv.T
            ci = Ch.chebint(c, lbnd=-1, axis=-1) * g.E / 2
            return ci @ VN.T
        res = amap(v, cum0)
        self.cache[key] = res
        return res
    def X1grad(self, base):
        """sum_j X1_j d_j F1 on the worldline for first-order field descriptor base = ('f', name, comp, 1)."""
        tot = {}
        for j in range(3):
            der = [0, 0, 0, 0]; der[j + 1] = 1
            F = self.field(base + (tuple(der),), 'wl')
            tot = aadd(tot, amul(self.X1(j), F))
        return tot
    def worldline_f2half(self):
        f = ascale(self.field(('f', 'Psi', None, 2, (0, 0, 0, 0)), 'wl'), 0.5)
        P1 = self.field(('f', 'Psi', None, 1, (0, 0, 0, 0)), 'wl')
        f = aadd(f, ascale(amul(P1, P1), -0.5))
        for c in range(3):
            v = self.field(('f', 'v', c, 1, (0, 0, 0, 0)), 'wl')
            f = aadd(f, ascale(amul(v, v), -0.5))
        f = aadd(f, self.X1grad(('f', 'Psi', None, 1)))
        return apart(f, 2)
    # ---------------- ray functionals
    def functional(self, kind, order, comp):
        key = ('fun', kind, order, comp)
        if key in self.cache: return self.cache[key]
        r = self.ray
        if kind == 'k' and order == 1:
            k1i = self.evaluate(RC['k_init'][comp].c[1].subs(NZ), 'obs')
            G1 = self.evaluate(RC['G1'][comp].subs(NZ), 'ray')
            res = aadd(amap(k1i, lambda v: v[..., None] * np.ones(r.M)), amap(G1, r.cum), -1.0)
        elif kind == 'Dx' and order == 1:
            res = amap(self.symbol('dxo', 1, comp), lambda v: v[..., None] * np.ones(r.M))
            if comp == 2: res = aadd(res, amap(self.symbol('deta_o', 1, None), lambda v: v[..., None] * np.ones(r.M)))
            res = aadd(res, amap(self.functional('k', 1, comp), r.cum), -1.0)
        elif kind == 'k' and order == 2:
            res = ascale(self.k2half(comp), 2.0)
        elif kind == 'Dx' and order == 2:
            h = ascale(self.symbol('dxo', 2, comp), 0.5)
            if comp == 2: h = aadd(h, ascale(self.symbol('deta_o', 2, None), 0.5))
            h = amap(h, lambda v: v[..., None] * np.ones(r.M))
            h = aadd(h, amap(self.k2half(comp), r.cum), -1.0)
            k1o = amap(self.functional('k', 1, comp), lambda v: v[..., -1])
            h = aadd(h, amap(apart(amul(self.symbol('deta_o', 1, None), k1o), 2), lambda v: v[..., None] * np.ones(r.M)), -1.0)
            res = ascale(apart(h, 2), 2.0)
        else:
            raise ValueError(key)
        self.cache[key] = res
        return res
    def k2half(self, comp):
        key = ('k2half', comp)
        if key in self.cache: return self.cache[key]
        r = self.ray
        # integrand G2/2 + Dx1.grad G1  (G2half contains abstract k1_i(eta))
        integ = self.evaluate(RC['G2half'][comp].subs(NZ), 'ray', maxorder=2)
        G1e = RC['G1'][comp].subs(NZ)
        for j, X in enumerate((XS, YS, ZS)):
            integ = aadd(integ, apart(amul(self.functional('Dx', 1, j), self.evaluate(sp.diff(G1e, X), 'ray')), 2))
        # initial value at the background observer event
        k2o = self.k2init_frame(comp)
        k1i = RC['k_init'][comp].c[1].subs(NZ)
        d1 = self.symbol('deta_o', 1, None)
        k2o = aadd(k2o, amul(d1, self.evaluate(sp.diff(k1i, ETA), 'obs')))
        for j, X in enumerate((XS, YS, ZS)):
            k2o = aadd(k2o, amul(self.symbol('dxo', 1, j), self.evaluate(sp.diff(k1i, X), 'obs')))
        k2o = aadd(k2o, amul(d1, self.evaluate(RC['G1'][comp].subs(NZ), 'obs')), -1.0)
        k2o = apart(k2o, 2)
        res = aadd(amap(k2o, lambda v: v[..., None] * np.ones(r.M)), amap(apart(integ, 2), r.cum), -1.0)
        self.cache[key] = res
        return res
    def k2init_frame(self, comp):
        key = ('k2init',)
        if key not in self.cache:
            glob = [apart(self.evaluate(RC['k_init'][i].c[2], 'obs', basis='global'), 2) for i in range(3)]
            fr = []
            for a in range(3):
                tot = {}
                for i in range(3):
                    tot = aadd(tot, ascale(glob[i], self.E[:, a, i]))
                fr.append(tot)
            self.cache[key] = fr
        return self.cache[key][comp]
    # ---------------- term lists
    def termlist(self, terms, maxorder=3):
        """sum over C2 term dicts -> algebra element of point arrays (B,Nn); ray terms integrated."""
        tot = {}; parts = {}
        for term in terms:
            where = term['where']
            el = self.evaluate(term['expr'], where, maxorder=maxorder)
            if where == 'ray': el = amap(el, self.ray.integ)
            parts[term['name']] = el
            tot = aadd(tot, el)
        return tot, parts
