"""pw_eval.py -- evaluation layer for the plane-wave solution (A1).
PlaneWave: jets  d_eta^kt d_x^kx X_n (eta, x)  of every Poisson-gauge field in C2's naming
  (Phi_n, Psi_n, dT_n, v_i_n, om_i_n, h_ij_n), built from the mode functions of pw_numeric.py;
  also uhat^0 order parts and the temperature contrast dT_n = n! [eps^n]((1+delta)^{1/4} - 1).
NSer: numeric truncated eps-series (coefficients of eps^0..eps^3, mpmath numbers).
expr_evaluator: compiles a sympy expression of C2 (fields as Functions of (eta,x,y,z), symbols) into a callable
  that takes a point (eta, x) and a dict of symbol values; the jets are supplied by a PlaneWave.
"""
import sympy as sp, os, sys, pickle, math
from mpmath import mp, mpf, mpc, sin, cos, exp, pi, fac
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from cheb import Cheb, Fun, nodes
from pw_numeric import load_fields, make_evaluator, MODES, amp, eta, gauge_c

eps = sp.Symbol('eps')
NORD = 3

# ------------------------------------------------------------------ numeric eps-series
class NSer:
    __slots__ = ('c',)
    def __init__(self, c):
        c = [mpf(v) if not isinstance(v, (mpf, mpc)) else v for v in c]
        self.c = c[:NORD + 1] + [mpf(0)] * (NORD + 1 - len(c))
    @staticmethod
    def const(a): return NSer([a])
    def __add__(self, o):
        o = o if isinstance(o, NSer) else NSer([o]); return NSer([a + b for a, b in zip(self.c, o.c)])
    __radd__ = __add__
    def __neg__(self): return NSer([-a for a in self.c])
    def __sub__(self, o): return self + (-(o if isinstance(o, NSer) else NSer([o])))
    def __rsub__(self, o): return NSer([o]) + (-self)
    def __mul__(self, o):
        if not isinstance(o, NSer):
            return NSer([a * o for a in self.c])
        r = [mpf(0)] * (NORD + 1)
        for i, a in enumerate(self.c):
            for j in range(NORD + 1 - i):
                r[i + j] += a * o.c[j]
        return NSer(r)
    __rmul__ = __mul__
    def apply(self, derivs):
        """f(self) where derivs[k] = f^{(k)}(c0), k = 0..NORD."""
        y = NSer([0] + self.c[1:])
        r = NSer([derivs[0]]); p = NSer([1])
        for k in range(1, NORD + 1):
            p = p * y
            r = r + p * (derivs[k] / fac(k))
        return r
    def inv(self):
        c0 = self.c[0]
        return self.apply([(-1) ** k * fac(k) / c0 ** (k + 1) for k in range(NORD + 1)])
    def __truediv__(self, o):
        return self * (o.inv() if isinstance(o, NSer) else 1 / mpf(o))
    def exp(self):
        e0 = exp(self.c[0]); return self.apply([e0] * (NORD + 1))
    def log(self):
        c0 = self.c[0]
        return self.apply([mp.log(c0)] + [(-1) ** (k + 1) * fac(k - 1) / c0 ** k for k in range(1, NORD + 1)])
    def pow(self, a):
        c0 = self.c[0]
        d = []
        for k in range(NORD + 1):
            coef = mpf(1)
            for j in range(k): coef *= (a - j)
            d.append(coef * c0 ** (a - k))
        return self.apply(d)
    def __repr__(self): return 'NSer(%s)' % ', '.join(mp.nstr(v, 15) for v in self.c)

# ------------------------------------------------------------------ fields
def parse_name(name):
    """C2 field function name -> (kind, index, order).  kind in Phi,Psi,dT,v,om,h."""
    if name.startswith('v_'):
        _, i, n = name.split('_'); return ('v', int(i), int(n))
    if name.startswith('om_'):
        _, i, n = name.split('_'); return ('om', int(i), int(n))
    if name.startswith('h_'):
        _, ij, n = name.split('_'); return ('h', ij, int(n))
    return (name[:-1], None, int(name[-1]))

class PlaneWave:
    """zeta = eps cos x; mode functions from pw_fields.pkl."""
    def __init__(self, fn=None):
        d = load_fields(fn)
        self.reg = d['reg']; self.consts = d['consts']; self.series = d['series']; self.meta = d['meta']
        self.uhat0 = d['uhat0']
        self.L = mpf(d['meta']['L'])
        self.gauge_c = sp.parse_expr(d['meta']['GAUGE_C'])
        self._build_dT()
        self._build_uhat0()
    # modes: dict kind -> {order: {m: (Fun, trig)}} with trig 'cos' or 'sin'; v^x = d_x V
    def _amp(self, kind, n, m):
        return self.reg.get('%s%d_%d' % (kind, n, m))
    def _build_dT(self):
        """dT_n mode functions from delta: dT = (1+delta)^{1/4} - 1 (series in eps)."""
        z = sp.Symbol('z')
        delta_ser = 0
        for n in (1, 2, 3):
            for m in MODES[n]:
                A = amp('delta', n, m)
                delta_ser += eps ** n / sp.factorial(n) * A * ((z ** m + z ** (-m)) / 2 if m > 0 else 1)   # (bug fixed: m=0 gave 2A)
        dT = sp.series((1 + delta_ser) ** sp.Rational(1, 4) - 1, eps, 0, NORD + 1).removeO()
        self.dT_expr = {}
        xs = nodes(0, self.L, len(self.reg['Phi1_1'].ch.c) - 1)
        for n in (1, 2, 3):
            cn = sp.expand(dT.coeff(eps, n) * sp.factorial(n))
            P = sp.Poly(sp.expand(cn * z ** 4), z)
            for (e,), c in P.terms():
                m = e - 4
                if m < 0: continue
                A = sp.expand(c * (2 if m > 0 else 1))
                self.dT_expr[(n, m)] = A
                ev = make_evaluator(A, self.reg)
                self.reg['dT%d_%d' % (n, m)] = Fun(Cheb.fit(0, self.L, [ev(x) for x in xs]))
    def _build_uhat0(self):
        """uhat^0 = 1 + sum_n eps^n u0_n(eta,x): evaluators for the z-coefficients."""
        self.u0 = {}
        for n in (1, 2, 3):
            for m, ex in self.uhat0[n].items():
                self.u0[(n, m)] = make_evaluator(ex, self.reg)
    # ---- jets
    def jet(self, name, kt, kx, ky, kz, et, x):
        """d_eta^kt d_x^kx d_y^ky d_z^kz of C2's field `name` at (eta, x)  (fields independent of y,z)."""
        if ky or kz:
            return mpf(0)
        kind, idx, n = parse_name(name)
        if kind in ('om',):
            return mpf(0)
        if kind == 'h':
            if n != 2 or idx not in ('11', '22', '33'):
                return mpf(0)
            if kx: return mpf(0)
            H = self.reg['H2_0'].val(kt, et)
            return H if idx == '11' else -H / 2
        if kind == 'v':
            if idx != 1:
                return mpf(0)
            tot = mpf(0)
            for m in MODES[n]:
                if m == 0: continue
                F = self.reg['V%d_%d' % (n, m)]
                # v^x = d_x [V cos(mx)] = -m V sin(mx);  d_x^kx sin(mx) = m^kx sin(mx + kx pi/2)
                tot += -m * F.val(kt, et) * mpf(m) ** kx * sin(m * x + kx * pi / 2)
            return tot
        if kind in ('Phi', 'Psi', 'dT'):
            tot = mpf(0)
            for m in MODES[n]:
                F = self.reg['%s%d_%d' % (kind, n, m)]
                if m == 0:
                    if kx: continue
                    tot += F.val(kt, et)
                else:
                    tot += F.val(kt, et) * mpf(m) ** kx * cos(m * x + kx * pi / 2)
            return tot
        raise ValueError(name)
    def uhat0_jet(self, n, kt, kx, et, x):
        """d_eta^kt d_x^kx of the eps^n coefficient of uhat^0 at (eta, x)  (n>=1)."""
        tot = mpf(0)
        for m in sorted({abs(m) for m in self.uhat0[n]}):
            if m == 0:
                if kx: continue
                tot += self._u0_val(n, 0, kt, et)
            else:
                A = 2 * self._u0_val(n, m, kt, et)      # real amplitude of cos(mx)
                tot += A * mpf(m) ** kx * cos(m * x + kx * pi / 2)
        return tot
    def _u0_val(self, n, m, kt, et):
        key = (n, m)
        if not hasattr(self, '_u0fun'):
            self._u0fun = {}
        if key not in self._u0fun:
            xs = nodes(0, self.L, len(self.reg['Phi1_1'].ch.c) - 1)
            ev = self.u0[key]
            self._u0fun[key] = Fun(Cheb.fit(0, self.L, [ev(x) for x in xs]))
        return self._u0fun[key].val(kt, et)

# ------------------------------------------------------------------ C2 expression compiler
X4 = sp.symbols('eta x y z')
def jet_key(d):
    """Derivative or AppliedUndef -> (name, kt, kx, ky, kz)."""
    if isinstance(d, sp.Derivative):
        f = d.expr; cnt = {}
        for v, c in d.variable_count: cnt[str(v)] = int(c)
        return (str(f.func), cnt.get('eta', 0), cnt.get('x', 0), cnt.get('y', 0), cnt.get('z', 0))
    return (str(d.func), 0, 0, 0, 0)

class Compiled:
    """callable(point_jets, symvals) for a C2 expression. point_jets: dict jet_key -> value; symvals: dict symbol-name -> value."""
    def __init__(self, expr, extra_funcs=()):
        """extra_funcs: names of AppliedUndef that are NOT fields (ray functionals k1_i(eta), Dx1_i(eta_star), ...);
        they are mapped to symbols named by str(f) (e.g. 'k1_1(eta)')."""
        expr = sp.expand(expr) if not isinstance(expr, (int, float)) else sp.sympify(expr)
        self.jets = []
        rep = {}
        for d in expr.atoms(sp.Derivative):
            key = jet_key(d)
            if key[0] in extra_funcs:
                raise ValueError('derivative of functional %s' % key[0])
            s = sp.Symbol('J_%s_%d%d%d%d' % key); rep[d] = s
            if key not in self.jets: self.jets.append(key)
        e2 = expr.xreplace(rep)
        rep2 = {}
        self.funcs = []
        for f in e2.atoms(sp.core.function.AppliedUndef):
            nm = str(f.func)
            if nm in extra_funcs:
                s = sp.Symbol('F_' + str(f).replace('(', '_').replace(')', '').replace(',', '_').replace(' ', ''))
                rep2[f] = s; self.funcs.append((str(f), s))
            else:
                key = jet_key(f); s = sp.Symbol('J_%s_%d%d%d%d' % key); rep2[f] = s
                if key not in self.jets: self.jets.append(key)
        e2 = e2.xreplace(rep2)
        self.syms = sorted(e2.free_symbols - {sp.Symbol('J_%s_%d%d%d%d' % k) for k in self.jets} - {s for _, s in self.funcs}, key=str)
        self.jetsyms = [sp.Symbol('J_%s_%d%d%d%d' % k) for k in self.jets]
        self.funcsyms = [s for _, s in self.funcs]
        self.funcnames = [nm for nm, _ in self.funcs]
        self.symnames = [str(s) for s in self.syms]
        self.nterms = len(sp.Add.make_args(e2))
        self.f = sp.lambdify(self.jetsyms + self.funcsyms + self.syms, e2, modules='mpmath', cse=True)
    def __call__(self, jets, funcs, symvals):
        return self.f(*[jets[k] for k in self.jets], *[funcs[nm] for nm in self.funcnames], *[symvals[s] for s in self.symnames])

def point_jets(pw, keys, et, x):
    return {k: pw.jet(k[0], k[1], k[2], k[3], k[4], et, x) for k in keys}
