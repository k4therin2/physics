"""
tpoly.py -- truncated polynomial arithmetic in eps for perturbative GR (C1 engine).

A `Poly` is a sparse polynomial in "jet symbols" J(field, order, m) = d^m X_order,
where m=(nt,nx,ny,nz) is a derivative multi-index, with coefficients that are
rational numbers times integer powers of eta.  Each jet carries an eps-order
(the `order` label of the field), so the eps-order of a monomial is the sum of
the orders of its jets; products are truncated at MAXORD=3.

Convention (SPEC §4): X = eps X1 + eps^2 X2/2 + eps^3 X3/6.
A Poly at eps-order n therefore holds the coefficient of eps^n (including the 1/n!).

Monomial key: (order, eta_exp, mono) with mono a sorted tuple of jet ids (with repetition).
"""
from fractions import Fraction
import itertools, pickle
import sympy as sp

MAXORD = 3
COORDS = ('t', 'x', 'y', 'z')   # t = conformal time eta

# ---------------------------------------------------------------- jets --------
class JetRegistry:
    def __init__(self):
        self.ids = {}        # (field, order, m) -> id
        self.info = []       # id -> (field, order, m)
    def get(self, field, order, m):
        key = (field, order, tuple(m))
        if key not in self.ids:
            self.ids[key] = len(self.info)
            self.info.append(key)
        return self.ids[key]
    def order(self, jid):
        return self.info[jid][1]
    def field(self, jid):
        return self.info[jid][0]
    def multi(self, jid):
        return self.info[jid][2]

REG = JetRegistry()

def jet(field, order, m=(0, 0, 0, 0)):
    return Poly({(order, 0, (REG.get(field, order, m),)): Fraction(1)})

def eta_pow(n, c=1):
    return Poly({(0, n, ()): Fraction(c)})

def const(c):
    return Poly({(0, 0, ()): Fraction(c)}) if c != 0 else Poly({})

# ---------------------------------------------------------------- Poly --------
class Poly:
    def __init__(self, d=None):
        self.d = d if d is not None else {}

    # --- basic
    def copy(self):
        return Poly(dict(self.d))
    def is_zero(self):
        return not self.d
    def __bool__(self):
        return bool(self.d)
    def nterms(self):
        return len(self.d)

    def __add__(self, other):
        if not isinstance(other, Poly):
            other = const(other)
        d = dict(self.d)
        for k, v in other.d.items():
            nv = d.get(k, 0) + v
            if nv == 0:
                d.pop(k, None)
            else:
                d[k] = nv
        return Poly(d)
    __radd__ = __add__

    def __neg__(self):
        return Poly({k: -v for k, v in self.d.items()})

    def __sub__(self, other):
        if not isinstance(other, Poly):
            other = const(other)
        return self + (-other)
    def __rsub__(self, other):
        return (-self) + other

    def scale(self, c):
        c = Fraction(c)
        if c == 0:
            return Poly({})
        return Poly({k: v * c for k, v in self.d.items()})

    def __mul__(self, other):
        if not isinstance(other, Poly):
            return self.scale(other)
        d = {}
        for (o1, e1, m1), c1 in self.d.items():
            for (o2, e2, m2), c2 in other.d.items():
                o = o1 + o2
                if o > MAXORD:
                    continue
                key = (o, e1 + e2, tuple(sorted(m1 + m2)))
                nv = d.get(key, 0) + c1 * c2
                if nv == 0:
                    d.pop(key, None)
                else:
                    d[key] = nv
        return Poly(d)
    __rmul__ = __mul__

    def __truediv__(self, c):
        return self.scale(Fraction(1) / Fraction(c))

    def __pow__(self, n):
        r = const(1)
        for _ in range(n):
            r = r * self
        return r

    # --- order bookkeeping
    def order_part(self, n):
        return Poly({k: v for k, v in self.d.items() if k[0] == n})
    def truncate(self, n):
        return Poly({k: v for k, v in self.d.items() if k[0] <= n})
    def max_order(self):
        return max((k[0] for k in self.d), default=-1)
    def min_order(self):
        return min((k[0] for k in self.d), default=99)

    # --- derivatives
    def diff(self, var):
        """partial derivative w.r.t. coordinate var in COORDS (t=eta)."""
        i = COORDS.index(var)
        d = {}
        def add(key, val):
            nv = d.get(key, 0) + val
            if nv == 0:
                d.pop(key, None)
            else:
                d[key] = nv
        for (o, e, mono), c in self.d.items():
            if i == 0 and e != 0:
                add((o, e - 1, mono), c * e)
            # product rule over jets; group identical jets
            seen = {}
            for j in mono:
                seen[j] = seen.get(j, 0) + 1
            for j, mult in seen.items():
                f, od, m = REG.info[j]
                m2 = list(m); m2[i] += 1
                j2 = REG.get(f, od, tuple(m2))
                lst = list(mono)
                lst.remove(j)
                lst.append(j2)
                add((o, e, tuple(sorted(lst))), c * mult)
        return Poly(d)

    def dt(self): return self.diff('t')
    def dx(self): return self.diff('x')
    def dy(self): return self.diff('y')
    def dz(self): return self.diff('z')
    def sd(self, i):
        """spatial derivative, i in 0,1,2 -> x,y,z;  i == -1 -> t"""
        return self.diff(COORDS[i + 1])
    def lap(self):
        return self.dx().dx() + self.dy().dy() + self.dz().dz()

    # --- series functions (background must be a single monomial eta^n * c)
    def background(self):
        b = self.order_part(0)
        if len(b.d) > 1:
            raise ValueError("background not a single monomial: %s" % b.d)
        return b
    def inverse(self):
        b = self.background()
        if not b.d:
            raise ZeroDivisionError("no background")
        ((o, e, m), c), = b.d.items()
        binv = Poly({(0, -e, ()): 1 / c})
        pert = (self - b) * binv          # order >= 1
        r = const(1)
        term = const(1)
        for n in range(1, MAXORD + 1):
            term = term * (-pert)
            r = r + term
        return r * binv
    def sqrt_series(self):
        b = self.background()
        ((o, e, m), c), = b.d.items()
        if e % 2 != 0:
            raise ValueError("odd eta power in sqrt")
        cs = Fraction(c).limit_denominator()
        import math
        num, den = cs.numerator, cs.denominator
        rn, rd = math.isqrt(num), math.isqrt(den)
        if rn * rn != num or rd * rd != den:
            raise ValueError("non-square background %s" % c)
        sb = Poly({(0, e // 2, ()): Fraction(rn, rd)})
        pert = (self - b) * b.inverse()   # x, order >= 1
        # sqrt(1+x) = 1 + x/2 - x^2/8 + x^3/16
        r = const(1) + pert / 2 - (pert * pert) / 8 + (pert * pert * pert) / 16
        return r * sb
    def log1p_series(self):
        """log(self) for self = 1 + x, x of order>=1 (background must be exactly 1)."""
        b = self.background()
        ((o, e, m), c), = b.d.items()
        assert e == 0 and c == 1, "log1p needs background 1"
        x = self - b
        return x - (x * x) / 2 + (x * x * x) / 3

    # --- conversion to sympy
    def to_sympy(self, eta=None, funcs=None):
        """Return a sympy expression: jets -> Derivative(Function(field_order)(eta,x,y,z)).
        If funcs is given, it is a dict field -> callable(order, m) -> sympy expr (custom substitution)."""
        eta = eta if eta is not None else sp.Symbol('eta')
        X = sp.symbols('x y z')
        cache = {}
        def jet_expr(j):
            if j in cache:
                return cache[j]
            f, od, m = REG.info[j]
            if funcs is not None:
                ex = funcs(f, od, m)
            else:
                F = sp.Function('%s%d' % (f, od))(eta, *X)
                ex = F
                for v, n in zip((eta,) + X, m):
                    if n:
                        ex = sp.Derivative(ex, (v, n))
            cache[j] = ex
            return ex
        terms = []
        for (o, e, mono), c in self.d.items():
            t = sp.Rational(c.numerator, c.denominator) * eta ** e
            for j in mono:
                t = t * jet_expr(j)
            terms.append(t)
        return sp.Add(*terms)

    def __repr__(self):
        return "Poly(%d terms, orders %s)" % (len(self.d), sorted({k[0] for k in self.d}))

# ----------------------------------------------------------- persistence ------
def save_polys(obj, fname):
    """Save a (nested) structure of Polys together with the jet registry."""
    def conv(o):
        if isinstance(o, Poly):
            return ('POLY', {k: (v.numerator, v.denominator) for k, v in o.d.items()})
        if isinstance(o, dict):
            return {k: conv(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return type(o)(conv(v) for v in o)
        return o
    with open(fname, 'wb') as fh:
        pickle.dump({'registry': REG.info, 'data': conv(obj)}, fh)

def load_polys(fname):
    with open(fname, 'rb') as fh:
        blob = pickle.load(fh)
    # merge registry: map saved ids to current ids
    idmap = {}
    for i, (f, od, m) in enumerate(blob['registry']):
        idmap[i] = REG.get(f, od, m)
    def conv(o):
        if isinstance(o, tuple) and len(o) == 2 and o[0] == 'POLY':
            return Poly({(k[0], k[1], tuple(sorted(idmap[j] for j in k[2]))): Fraction(*v)
                         for k, v in o[1].items()})
        if isinstance(o, dict):
            return {k: conv(v) for k, v in o.items()}
        if isinstance(o, list):
            return [conv(v) for v in o]
        if isinstance(o, tuple):
            return tuple(conv(v) for v in o)
        return o
    return conv(blob['data'])

# ------------------------------------------------------------- helpers --------
def field_series(name, orders=(1, 2, 3)):
    """X = eps X1 + eps^2 X2/2 + eps^3 X3/6 as a Poly (eps absorbed in the order label)."""
    import math
    r = Poly({})
    for n in orders:
        r = r + jet(name, n) / math.factorial(n)
    return r
