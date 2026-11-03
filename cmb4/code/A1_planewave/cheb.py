"""cheb.py -- Chebyshev series in mpmath (A1).  Used as the high-precision numerical EVALUATION layer for the
exact expressions: entire functions of eta are represented by Chebyshev-Gauss interpolants on [a,b]; derivatives
and antiderivatives are exact operations on the coefficient vector (spectral integration).  All arithmetic in mpmath
at the current working precision.

Cheb(a, b, coef): f(x) = sum_j c_j T_j(t),  t = (2x - a - b)/(b - a).
nodes(a, b, N): the N+1 Chebyshev-Gauss nodes (interior points, ascending).
Cheb.fit(a, b, vals): interpolant through the values at nodes(a, b, len(vals)-1).
"""
from mpmath import mp, mpf, cos, pi, fsum

_node_cache = {}
_T_cache = {}

def nodes(a, b, N):
    """N+1 Chebyshev-Gauss nodes on [a,b] (ascending)."""
    key = (mpf(a), mpf(b), N, mp.prec)
    if key not in _node_cache:
        a, b = mpf(a), mpf(b)
        t = [cos(pi * (2 * k + 1) / (2 * (N + 1))) for k in range(N + 1)]      # descending in t
        xs = [(b + a) / 2 + (b - a) / 2 * tk for tk in t]
        _node_cache[key] = list(reversed(xs))
    return _node_cache[key]

def _Tmatrix(N):
    """T_j(t_k) for the Gauss nodes t_k (ascending order matching nodes())."""
    key = (N, mp.prec)
    if key not in _T_cache:
        t = [cos(pi * (2 * k + 1) / (2 * (N + 1))) for k in range(N + 1)]
        t = list(reversed(t))
        T = [[mpf(1)] * (N + 1), list(t)]
        for j in range(2, N + 1):
            T.append([2 * t[k] * T[j - 1][k] - T[j - 2][k] for k in range(N + 1)])
        _T_cache[key] = T
    return _T_cache[key]

class Cheb:
    __slots__ = ('a', 'b', 'c')
    def __init__(self, a, b, coef):
        self.a, self.b, self.c = mpf(a), mpf(b), [mpf(v) for v in coef]
    @classmethod
    def fit(cls, a, b, vals):
        N = len(vals) - 1
        T = _Tmatrix(N)
        c = []
        for j in range(N + 1):
            s = fsum(vals[k] * T[j][k] for k in range(N + 1)) * 2 / (N + 1)
            c.append(s)
        c[0] /= 2
        return cls(a, b, c)
    @classmethod
    def from_function(cls, a, b, f, N):
        xs = nodes(a, b, N)
        return cls.fit(a, b, [f(x) for x in xs])
    def __call__(self, x):
        t = (2 * mpf(x) - self.a - self.b) / (self.b - self.a)
        c = self.c
        b1 = b2 = mpf(0)
        t2 = 2 * t
        for j in range(len(c) - 1, 0, -1):
            b1, b2 = t2 * b1 - b2 + c[j], b1
        return t * b1 - b2 + c[0]
    def deriv(self):
        c = self.c; N = len(c) - 1
        d = [mpf(0)] * (N + 2)
        for j in range(N, 0, -1):
            d[j - 1] = d[j + 1] + 2 * j * c[j]
        d[0] /= 2
        scale = 2 / (self.b - self.a)
        return Cheb(self.a, self.b, [v * scale for v in d[:N + 1]])
    def antideriv(self, x0=None):
        """F with F' = f and F(x0) = 0 (x0 defaults to a)."""
        c = self.c; N = len(c) - 1
        cc = list(c) + [mpf(0), mpf(0)]
        C = [mpf(0)] * (N + 2)
        for j in range(1, N + 2):
            C[j] = (cc[j - 1] - cc[j + 1]) / (2 * j)
        C[1] = (2 * cc[0] - cc[2]) / 2     # j=1: (2 c0 - c2)/2  (c0 is the full constant term)
        scale = (self.b - self.a) / 2
        F = Cheb(self.a, self.b, [v * scale for v in C])
        x0 = self.a if x0 is None else x0
        F.c[0] -= F(x0)
        return F
    def integral(self, x1=None, x2=None):
        F = self.antideriv()
        x1 = self.a if x1 is None else x1; x2 = self.b if x2 is None else x2
        return F(x2) - F(x1)
    def __add__(self, o):
        if isinstance(o, Cheb):
            n = max(len(self.c), len(o.c))
            return Cheb(self.a, self.b, [(self.c[j] if j < len(self.c) else 0) + (o.c[j] if j < len(o.c) else 0) for j in range(n)])
        c = list(self.c); c[0] += o; return Cheb(self.a, self.b, c)
    __radd__ = __add__
    def __neg__(self): return Cheb(self.a, self.b, [-v for v in self.c])
    def __sub__(self, o): return self + (-o)
    def __rsub__(self, o): return (-self) + o
    def scale(self, s): return Cheb(self.a, self.b, [v * s for v in self.c])
    def tail(self, m=5):
        """size of the last m coefficients relative to the largest (convergence diagnostic)."""
        mx = max(abs(v) for v in self.c)
        return max(abs(v) for v in self.c[-m:]) / mx if mx else mpf(0)

class Fun:
    """A function of eta on [a,b] with cached derivatives: Fun.d(k)(x)."""
    def __init__(self, cheb):
        self.ch = cheb; self._d = {0: cheb}
    def d(self, k):
        if k not in self._d:
            self._d[k] = self.d(k - 1).deriv()
        return self._d[k]
    def __call__(self, x): return self.ch(x)
    def val(self, k, x): return self.d(k)(x)
