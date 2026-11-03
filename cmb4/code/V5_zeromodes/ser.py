# V5 copy of C2_observable/ser.py with N = 2
"""Truncated power series in eps (order <= N) with sympy coefficients.
Convention (SPEC §4): X = eps X1 + eps^2 X2/2 + eps^3 X3/6, so the stored
coefficient list is c = [X0, X1, X2/2, X3/6]."""
import sympy as sp
N = 2
eps = sp.Symbol('eps')

def _S(o):
    return o if isinstance(o, Ser) else Ser([sp.sympify(o)])

class Ser:
    __slots__ = ('c',)
    def __init__(self, c, n=N):
        c = list(c)
        self.c = [sp.sympify(t) for t in c[:n+1]] + [sp.Integer(0)]*(n+1-len(c))
    @staticmethod
    def const(a): return Ser([sp.sympify(a)])
    def __add__(self, o): o = _S(o); return Ser([a+b for a, b in zip(self.c, o.c)])
    __radd__ = __add__
    def __neg__(self): return Ser([-a for a in self.c])
    def __sub__(self, o): return self + (-_S(o))
    def __rsub__(self, o): return _S(o) + (-self)
    def __mul__(self, o):
        o = _S(o); r = [sp.Integer(0)]*(N+1)
        for i, a in enumerate(self.c):
            if a == 0: continue
            for j in range(N+1-i):
                b = o.c[j]
                if b == 0: continue
                r[i+j] += a*b
        return Ser([sp.expand(t) for t in r])
    __rmul__ = __mul__
    def __truediv__(self, o): return self * _S(o).inv()
    def __rtruediv__(self, o): return _S(o) * self.inv()
    def diff(self, *args): return Ser([sp.diff(a, *args) for a in self.c])
    def expand(self): return Ser([sp.expand(a) for a in self.c])
    def subs(self, *a, **k): return Ser([t.subs(*a, **k) for t in self.c])
    def apply(self, f):
        """f: sympy expression in symbol t; returns f(self) as truncated series (c0 must be eps-free)."""
        t = sp.Symbol('t_apply')
        c0 = self.c[0]
        y = Ser([0] + self.c[1:])
        r = Ser.const(f.subs(t, c0))
        p = Ser.const(1)
        for k in range(1, N+1):
            p = p*y
            r = r + p*(sp.diff(f, t, k).subs(t, c0)/sp.factorial(k))
        return r
    def inv(self):
        t = sp.Symbol('t_apply'); return self.apply(1/t)
    def sqrt(self):
        t = sp.Symbol('t_apply'); return self.apply(sp.sqrt(t))
    def log(self):
        t = sp.Symbol('t_apply'); return self.apply(sp.log(t))
    def exp(self):
        t = sp.Symbol('t_apply'); return self.apply(sp.exp(t))
    def tot(self): return sum(eps**i*a for i, a in enumerate(self.c))
    def order(self, i): return self.c[i]
    def iszero(self): return all(sp.expand(a) == 0 for a in self.c)
    def __repr__(self): return 'Ser(%s)' % self.c

def dot4(a, b, g):
    """g_{mu nu} a^mu b^nu for lists of Ser."""
    r = Ser.const(0)
    for m in range(4):
        for n in range(4):
            if g[m][n].iszero(): continue
            r = r + g[m][n]*a[m]*b[n]
    return r
