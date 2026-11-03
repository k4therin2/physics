"""pw_map.py -- Stage C of A1: the full map Theta_1,2,3 for zeta = eps cos x from C2's line-of-sight term lists.

Inputs: pw_fields.pkl (Poisson-gauge mode functions from C1's equations, pw_symbolic.py + pw_numeric.py),
        C2_observable/L_terms.pkl (general-n shear-form term lists L[1..3], emission shifts),
        C2_observable/ray_cache.pkl (G1, G2half, k_init: the ray functionals' building blocks).
Geometry: fields depend on x only; sky direction n = (mu, sqrt(1-mu^2), 0) (mu = n.xhat, N1's polar axis +x);
background ray xbar(eta) = x_o + n (etabar_o - eta) (C2 convention, e = -n).
Ray functionals k1, Dx1, k2, Dx2 (C2 definitions) by Chebyshev spectral integration on Lobatto nodes in
[eta_*, etabar_o] (vectorised over mu, float64); observer clock/displacement in mpmath (exact series coefficients
of the perturbative expansion, quadrature of the exact mode functions).
"""
import sympy as sp, numpy as np, pickle, os, sys, time, math
from numpy.polynomial import chebyshev as C
from mpmath import mp, mpf
HERE = os.path.dirname(os.path.abspath(__file__))
C2 = os.path.join(HERE, '..', 'C2_observable')
sys.path.insert(0, HERE); sys.path.insert(0, C2)
from pw_eval import PlaneWave, MODES

X4 = sp.symbols('eta x y z')
ZERO_FIELDS = ('v_2_', 'v_3_', 'om_', 'h_12_', 'h_13_', 'h_23_')

# ------------------------------------------------------------------ C2 inputs
def load_c2():
    d = pickle.load(open(os.path.join(C2, 'L_terms.pkl'), 'rb'))
    r = pickle.load(open(os.path.join(C2, 'ray_cache.pkl'), 'rb'))
    return d, r

def planar_reduce(e):
    """impose x-only dependence, n3 = 0, vanishing fields; returns expanded expression."""
    e = sp.sympify(e)
    rep = {}
    for dd in e.atoms(sp.Derivative):
        vc = dict((str(v), c) for v, c in dd.variable_count)
        nm = str(dd.expr.func)
        if vc.get('y', 0) or vc.get('z', 0) or nm.startswith(ZERO_FIELDS):
            rep[dd] = 0
    e = e.xreplace(rep)
    rep = {f: 0 for f in e.atoms(sp.core.function.AppliedUndef) if str(f.func).startswith(ZERO_FIELDS)}
    e = e.xreplace(rep).xreplace({sp.Symbol('n3'): 0})
    return sp.expand(e)

FUNCTIONALS = ('k1_', 'k2_', 'Dx1_', 'Dx2_')
def jet_key(a):
    if isinstance(a, sp.Derivative):
        vc = dict((str(v), int(c)) for v, c in a.variable_count)
        return (str(a.expr.func), vc.get('eta', 0), vc.get('x', 0))
    return (str(a.func), 0, 0)

class NpExpr:
    """compile a (planar-reduced) C2 expression for numpy evaluation."""
    def __init__(self, expr):
        expr = planar_reduce(expr)
        self.expr = expr
        args = []; rep = {}
        for a in sorted(expr.atoms(sp.Derivative), key=str):
            k = jet_key(a); s = sp.Symbol('J_%s_%d_%d' % k); rep[a] = s; args.append(('jet', k, s))
        e2 = expr.xreplace(rep)
        for f in sorted(e2.atoms(sp.core.function.AppliedUndef), key=str):
            nm = str(f.func)
            if nm.startswith(FUNCTIONALS):
                s = sp.Symbol('F_' + nm); rep[f] = s; args.append(('fun', nm, s))
            else:
                k = (nm, 0, 0); s = sp.Symbol('J_%s_0_0' % nm); rep[f] = s; args.append(('jet', k, s))
        e2 = expr.xreplace(rep)
        used = {s for _, _, s in args}
        for s in sorted(e2.free_symbols - used, key=str):
            args.append(('sym', str(s), s))
        self.args = args
        self.nterms = len(sp.Add.make_args(e2))
        self.f = sp.lambdify([s for _, _, s in args], e2, modules='numpy', cse=True)
    def __call__(self, ctx):
        vals = []
        for kind, k, s in self.args:
            if kind == 'jet': vals.append(ctx.jet(*k))
            elif kind == 'fun': vals.append(ctx.fun(k))
            else: vals.append(ctx.sym(k))
        return self.f(*vals) if vals else self.f()

# ------------------------------------------------------------------ numeric field tables
def parse_name(name):
    if name.startswith('v_'):
        _, i, n = name.split('_'); return ('v', int(i), int(n))
    if name.startswith('h_'):
        _, ij, n = name.split('_'); return ('h', ij, int(n))
    return (name[:-1], None, int(name[-1]))

class Tables:
    """d^k A(eta) for every mode amplitude at a list of eta values (computed from the mpmath Chebyshev Funs)."""
    def __init__(self, pw, etas, KT=4):
        self.etas = etas
        self.A = {}
        for nm, F in pw.reg.items():
            self.A[nm] = [np.array([float(F.val(k, mpf(e))) for e in etas]) for k in range(KT + 1)]
        self.KT = KT

class Ctx:
    """evaluation context: eta indices (array, broadcast with X), X positions, functionals, symbols."""
    def __init__(self, tab, ieta, X, funs, syms):
        self.tab, self.ieta, self.X, self.funs, self.syms = tab, ieta, X, funs, syms
    def amp(self, nm, kt):
        assert kt <= self.tab.KT, (nm, kt)
        return self.tab.A[nm][kt][self.ieta]
    def jet(self, name, kt, kx):
        kind, idx, n = parse_name(name)
        X = self.X
        if kind == 'h':
            if n != 2 or idx not in ('11', '22', '33') or kx: return np.zeros_like(X)
            H = self.amp('H2_0', kt) + 0 * X
            return H if idx == '11' else -H / 2
        if kind == 'v':
            assert idx == 1
            tot = 0 * X
            for m in MODES[n]:
                if m == 0: continue
                tot = tot - m * self.amp('V%d_%d' % (n, m), kt) * m ** kx * np.sin(m * X + kx * np.pi / 2)
            return tot
        if kind in ('Phi', 'Psi', 'dT', 'delta'):
            tot = 0 * X
            for m in MODES[n]:
                A = self.amp('%s%d_%d' % (kind, n, m), kt)
                if m == 0:
                    if kx == 0: tot = tot + A
                else:
                    tot = tot + A * m ** kx * np.cos(m * X + kx * np.pi / 2)
            return tot
        raise ValueError(name)
    def fun(self, nm): return self.funs[nm]
    def sym(self, nm): return self.syms[nm]

# ------------------------------------------------------------------ Chebyshev-Lobatto machinery on [a, b]
def lobatto(a, b, N):
    j = np.arange(N + 1)
    t = np.cos(np.pi * j / N)           # t_0 = 1 -> b (observer), t_N = -1 -> a (emission)
    return (a + b) / 2 + (b - a) / 2 * t, t

def cheb_coef(vals, N):
    """vals[..., j] at Lobatto nodes t_j = cos(pi j/N) -> Chebyshev coefficients along last axis."""
    j = np.arange(N + 1); k = j[:, None]
    w = np.ones(N + 1); w[0] = w[-1] = 0.5
    T = np.cos(np.pi * k * j[None, :] / N)          # T[k, j]
    c = (2.0 / N) * np.tensordot(vals * w, T.T, axes=([-1], [0]))
    c[..., 0] /= 2; c[..., -1] /= 2
    return c

def integral_to_top(vals, a, b, N):
    """F(eta_j) = int_{eta_j}^{b} f deta for f given at Lobatto nodes (last axis)."""
    c = cheb_coef(vals, N)
    ci = C.chebint(c, axis=-1) * (b - a) / 2
    j = np.arange(N + 1); t = np.cos(np.pi * j / N)
    P = C.chebval(t, np.moveaxis(ci, -1, 0)) if ci.ndim > 1 else C.chebval(t, ci)
    Ptop = C.chebval(1.0, np.moveaxis(ci, -1, 0)) if ci.ndim > 1 else C.chebval(1.0, ci)
    return Ptop[..., None] - P

def total_integral(vals, a, b, N):
    return integral_to_top(vals, a, b, N)[..., -1]

# ------------------------------------------------------------------ observer (mpmath, exact perturbative coefficients)
def observer(pw, xL, etao, eta_i, lower0=True):
    """Clock and displacement of the fluid element.  X(eta) = xL + eps a1 + eps^2 a2 (coefficient form);
    a1(eta) = int_{l}^{eta} v1(eta', xL), a2(eta) = int_l^eta [v2/2 + a1 d_x v1];  l = 0 (label = position at eta->0,
    which equals N1's comoving label since the gauge map is O(eta^2)) or eta_i (SPEC §6 literally).
    Clock: Delta tau = int_{eta_init}^{eta_o} eta/uhat^0 deta = (etao^2 - eta_i^2)/2, initial event (1+dT)/eta_init = 1/eta_i.
    Returns dict with deta_o1, deta_o2 (C2 normalisation: deta_o = eps deta_o1 + eps^2 deta_o2/2), dxo1, dxo2 (same), T_o series."""
    mp.dps = 30
    J = lambda nm, kt, kx, e, x: pw.jet(nm, kt, kx, 0, 0, e, x)
    l = mpf(0) if lower0 else mpf(eta_i)
    xL = mpf(xL); etao = mpf(etao); eta_i = mpf(eta_i)
    a1 = lambda e: mp.quad(lambda s: J('v_1_1', 0, 0, s, xL), [l, e])
    def f1(e): return J('Psi1', 0, 0, e, xL)                                     # eps coefficient of 1/uhat0
    def f2(e):                                                                    # eps^2 coefficient
        A1 = a1(e)
        return (J('Psi2', 0, 0, e, xL) / 2 + A1 * J('Psi1', 0, 1, e, xL) - J('v_1_1', 0, 0, e, xL) ** 2 / 2
                - J('Psi1', 0, 0, e, xL) ** 2 / 2)
    I1 = mp.quad(lambda e: e * f1(e), [eta_i, 2, 4, 6, etao] if etao > 6 else [eta_i, 2, etao])
    I2 = mp.quad(lambda e: e * f2(e), [eta_i, 2, 4, 6, etao] if etao > 6 else [eta_i, 2, etao])
    # initial event: eta_init = eta_i (1 + dT(eta_init, X(eta_init)))
    e1 = eta_i * J('dT1', 0, 0, eta_i, xL)
    e2 = eta_i * (J('dT2', 0, 0, eta_i, xL) / 2 + e1 * J('dT1', 1, 0, eta_i, xL) + a1(eta_i) * J('dT1', 0, 1, eta_i, xL))
    D1 = (eta_i * e1 - I1) / etao
    D2 = -(D1 ** 2 / 2 + f1(etao) * etao * D1 + I2 - eta_i * e2 - e1 ** 2 / 2 - f1(eta_i) * eta_i * e1) / etao
    A1o = a1(etao)
    A2o = mp.quad(lambda e: J('v_1_2', 0, 0, e, xL) / 2 + a1(e) * J('v_1_1', 0, 1, e, xL), [l, 2, etao])
    dxo1 = A1o
    dxo2half = A2o + D1 * J('v_1_1', 0, 0, etao, xL)
    # T_o = (1 + dT(eta_o, X(eta_o)))/eta_o as eps series (coefficients of eps^0..2)
    dT1o = J('dT1', 0, 0, etao, xL); dT2o = J('dT2', 0, 0, etao, xL)
    T0 = 1 / etao
    T1 = dT1o / etao - D1 / etao ** 2
    T2 = (dT2o / 2 + D1 * J('dT1', 1, 0, etao, xL) + A1o * J('dT1', 0, 1, etao, xL)) / etao - dT1o * D1 / etao ** 2 \
         - D2 / etao ** 2 + D1 ** 2 / etao ** 3
    return dict(deta_o1=D1, deta_o2=2 * D2, dxo1=dxo1, dxo2=2 * dxo2half, eta_init=(e1, e2), T_o=(T0, T1, T2),
                eta_o=(etao, D1, D2))

# ------------------------------------------------------------------ the map
class MapBuilder:
    def __init__(self, pw, d, r, verbose=True):
        self.pw, self.d, self.r = pw, d, r
        t0 = time.time()
        Ls = d['L']
        self.terms = {m: [(t['name'], t['where'], NpExpr(t['expr'])) for t in Ls[m]] for m in (1, 2, 3)}
        self.G1 = [NpExpr(g) for g in r['G1']]
        k1i = [r['k_init'][i].c[1] for i in range(3)]
        k2ih = [r['k_init'][i].c[2] for i in range(3)]
        self.k1_init = [NpExpr(e) for e in k1i]
        dxo1 = sp.symbols('dxo1_1:4'); deta_o1 = sp.Symbol('deta_o1')
        X = X4
        k2oh = [k2ih[i] + deta_o1 * sp.diff(k1i[i], X[0]) + sum(dxo1[j] * sp.diff(k1i[i], X[j + 1]) for j in range(3))
                - deta_o1 * r['G1'][i] for i in range(3)]
        self.k2o_half = [NpExpr(e) for e in k2oh]
        Dx1 = [sp.Function('Dx1_%d' % i)(X[0]) for i in (1, 2, 3)]
        G2t = [r['G2half'][i] + sum(Dx1[j] * sp.diff(r['G1'][i], X[j + 1]) for j in range(3)) for i in range(3)]
        self.G2tot_half = [NpExpr(e) for e in G2t]
        em = d['emission']; de1, de2 = d['symbols']['deta_e']
        self.deta_e1 = NpExpr(em[de1]); self.deta_e2 = NpExpr(em[de2])
        if verbose:
            print('compiled: L1 %d, L2/2 %d, L3/6 %d terms (planar-reduced) (%.1fs)' % tuple(
                [sum(e.nterms for _, _, e in self.terms[m]) for m in (1, 2, 3)] + [time.time() - t0]))

    def run(self, xL, etao, etas, eta_i=0.01, nmu=96, N=64, lower0=True, verbose=True, obs=None):
        t0 = time.time()
        pw = self.pw
        obs = obs or observer(pw, xL, etao, eta_i, lower0)
        mu, wmu = np.polynomial.legendre.leggauss(nmu)
        etan, tn = lobatto(etas, etao, N)
        tab = Tables(pw, list(etan))
        n1 = mu[:, None]; n2 = np.sqrt(1 - mu ** 2)[:, None]
        Xray = xL + n1 * (etao - etan[None, :])
        syms = dict(n1=n1, n2=n2, eta_star=etas, eta_o_bar=etao,
                    deta_o1=float(obs['deta_o1']), deta_o2=float(obs['deta_o2']),
                    dxo1_1=float(obs['dxo1']), dxo1_2=0.0, dxo1_3=0.0, dxo2_1=float(obs['dxo2']), dxo2_2=0.0, dxo2_3=0.0)
        idx = np.arange(N + 1)[None, :]
        ray = Ctx(tab, idx, Xray, {}, syms)
        syms_o = dict(syms); syms_o['n1'] = mu; syms_o['n2'] = np.sqrt(1 - mu ** 2)
        octx = Ctx(tab, 0, xL + 0 * mu, {}, syms_o)
        Xe = xL + mu * (etao - etas)
        ectx = Ctx(tab, N, Xe, {}, dict(syms_o))
        sh = (nmu, N + 1)
        bc = lambda v: np.broadcast_to(v, sh).astype(float) if np.ndim(v) else np.full(sh, float(v))
        bco = lambda v: np.broadcast_to(v, (nmu,)).astype(float) if np.ndim(v) else np.full(nmu, float(v))
        F = {}
        # --- k1
        for i in range(3):
            G = bc(self.G1[i](ray)); k1o = bco(self.k1_init[i](octx))
            F['k1_%d' % (i + 1)] = k1o[:, None] - integral_to_top(G, etas, etao, N)     # dk/deta = G (C2's ray.definitions print '+': sign slip)
        n = [n1, n2, 0 * n1]
        dxo1 = [syms['dxo1_1'], 0.0, 0.0]; dxo2 = [syms['dxo2_1'], 0.0, 0.0]
        for i in range(3):
            F['Dx1_%d' % (i + 1)] = dxo1[i] + n[i] * syms['deta_o1'] - integral_to_top(F['k1_%d' % (i + 1)], etas, etao, N)
        ray.funs = F
        octx.funs = {k: v[:, 0] for k, v in F.items()}
        for i in range(3):
            G2 = bc(self.G2tot_half[i](ray)); k2o = bco(self.k2o_half[i](octx))
            F['k2_%d' % (i + 1)] = 2 * (k2o[:, None] - integral_to_top(G2, etas, etao, N))
        for i in range(3):
            F['Dx2_%d' % (i + 1)] = 2 * (dxo2[i] / 2 + n[i] * syms['deta_o2'] / 2 - syms['deta_o1'] * F['k1_%d' % (i + 1)][:, :1]
                                         - integral_to_top(F['k2_%d' % (i + 1)] / 2, etas, etao, N))
        ectx.funs = {k: v[:, N] for k, v in F.items()}
        de1 = bco(self.deta_e1(ectx)); de2 = bco(self.deta_e2(ectx))
        ectx.syms['deta_e1'] = de1; ectx.syms['deta_e2'] = de2
        # --- L terms
        out = {}
        groups = {}
        for m in (1, 2, 3):
            tot = np.zeros(nmu)
            for name, where, ex in self.terms[m]:
                if where == 'ray':
                    v = total_integral(bc(ex(ray)), etas, etao, N)
                elif where == 'obs':
                    v = bco(ex(octx))
                else:
                    v = bco(ex(ectx))
                groups[(m, name)] = v * math.factorial(m)
                tot = tot + v
            out['L%d' % m] = tot * math.factorial(m)
        L1, L2, L3 = out['L1'], out['L2'], out['L3']
        m1 = L1; m2 = L2 + L1 ** 2; m3 = L3 + 3 * L1 * L2 + L1 ** 3
        avg = lambda f: 0.5 * np.sum(wmu * f)
        a1, a2, a3 = avg(m1), avg(m2), avg(m3)
        Th1 = m1 - a1
        Th2 = m2 - a2 - 2 * a1 * (m1 - a1)
        Th3 = m3 - a3 - 3 * a1 * (m2 - a2) - 3 * a2 * (m1 - a1) + 6 * a1 ** 2 * (m1 - a1)
        res = dict(mu=mu, w=wmu, L=(L1, L2, L3), m=(m1, m2, m3), mono=(a1, a2, a3), Th=(Th1, Th2, Th3), groups=groups,
                   funs=F, obs=obs, etan=etan, de=(de1, de2))
        if verbose: print('map built: xL=%s etao=%s nmu=%d N=%d (%.1fs)' % (xL, etao, nmu, N, time.time() - t0))
        return res

def alm(mu, w, f, lmax=12):
    from scipy.special import eval_legendre
    return np.array([math.sqrt(4 * math.pi * (2 * l + 1)) * 0.5 * np.sum(w * f * eval_legendre(l, mu)) for l in range(lmax + 1)])
