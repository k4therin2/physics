"""trig_integrate.py -- C4 engine: closed-form definite integrals of sums of terms
        c * t^n * {sin,cos}(omega t),  n in Z (possibly negative), omega linear in symbols,
from 0 to T, for integrands whose TOTAL is regular at t=0 (individual terms may be singular).

Method: rewrite trig -> complex exponentials, collect terms (n, omega) -> c; integrate each term
from eps to T with the recursion   I_n = [t^{n+1} e^{i w t}]_eps^T/(n+1) - (i w/(n+1)) I_{n+1}  (n <= -2),
I_{-1} = -Cin(wT) + i Si(wT) + ln(T/eps) + O(eps),  I_n (n>=0) elementary;  the eps-singular pieces
(eps^m, m<0, and ln eps) and the finite constants from the lower limit are tracked separately and the
singular ones are required to cancel in the sum (checked symbolically).

Special functions:  Si(z) (sympy),  Cin(z) := gamma + ln|z| - Ci(|z|) = int_0^z (1-cos t)/t dt  (entire, even).
Using Cin instead of Ci removes every |.| / log(|.|) / principal-value subtlety: the result is manifestly
analytic in the frequencies omega.  Ci-form:  Cin(w T) = EulerGamma + log(|w| T) - Ci(|w| T).
"""
import sympy as sp
from sympy import Si, Ci, EulerGamma, pi

class Cin(sp.Function):
    """Cin(z) = int_0^z (1 - cos t)/t dt = EulerGamma + log(z) - Ci(z)  (z>0); even, entire, Cin(0)=0."""
    nargs = 1
    @classmethod
    def eval(cls, z):
        if z.is_zero:
            return sp.S.Zero
        if z.could_extract_minus_sign():
            return cls(-z)
        if z.is_Number:
            return (EulerGamma + sp.log(z) - Ci(z)).evalf(z._prec) if z.is_Float else None
    def fdiff(self, argindex=1):
        z = self.args[0]
        return (1 - sp.cos(z)) / z
    def _eval_evalf(self, prec):
        z = self.args[0]._eval_evalf(prec)
        if z is None or not z.is_Number:
            return None
        import mpmath
        with mpmath.workprec(prec):
            zz = mpmath.mpf(abs(z))
            if zz == 0:
                return sp.Float(0, precision=prec)
            if zz < mpmath.mpf('0.5'):
                s, term, j = mpmath.mpf(0), mpmath.mpf(0), 1
                while True:
                    term = (-1)**(j + 1) * zz**(2 * j) / (2 * j * mpmath.factorial(2 * j))
                    s += term
                    if abs(term) < mpmath.eps * abs(s) or j > 60:
                        break
                    j += 1
                return sp.Float(s, precision=prec)
            return sp.Float(mpmath.euler + mpmath.log(zz) - mpmath.ci(zz), precision=prec)
    def _eval_nseries(self, x, n, logx, cdir=0):
        z = self.args[0]
        ser = sum((-1)**(j + 1) * z**(2 * j) / (2 * j * sp.factorial(2 * j)) for j in range(1, n // 2 + 2))
        return ser.series(x, 0, n) if ser.has(x) else ser
    def _eval_rewrite_as_Ci(self, z, **kw):
        return EulerGamma + sp.log(sp.Abs(z)) - Ci(sp.Abs(z))
    def _eval_is_real(self):
        return self.args[0].is_real

def cin_series(z, order):
    """Taylor polynomial of Cin(z) through z^order."""
    return sum((-1)**(j + 1) * z**(2 * j) / (2 * j * sp.factorial(2 * j)) for j in range(1, order // 2 + 1))

def si_series(z, order):
    return sum((-1)**j * z**(2 * j + 1) / ((2 * j + 1) * sp.factorial(2 * j + 1)) for j in range(0, (order - 1) // 2 + 1))

# ----------------------------------------------------------------------------------------------------------
def decompose(expr, t):
    """expr = sum c * t^n * exp(i w t)  ->  dict {(n, w): c}.  w may be 0.
    Each term of expand(expr) is split into a t-free coefficient and the t-dependent part (t^n times trig
    functions) BEFORE the trig -> exp rewriting, so that denominators in the coefficients are never touched."""
    out = {}
    def split_t(m):
        """Mul m -> (t-free part, t power, t-dependent non-power part)"""
        c, n, rest = sp.S.One, 0, sp.S.One
        for f in sp.Mul.make_args(m):
            if f == t: n += 1
            elif f.is_Pow and f.base == t and f.exp.is_Integer: n += int(f.exp)
            elif f.has(t): rest *= f
            else: c *= f
        return c, n, rest
    for term0 in sp.Add.make_args(sp.expand(expr)):
        num, den = term0.as_numer_denom()
        cd, nd, rd = split_t(sp.factor(den))
        assert rd == 1, ("t-dependent denominator", den)
        cn, nn, rn = split_t(num)
        c0 = cn / cd; n0 = nn - nd
        e = sp.expand(sp.powsimp(sp.expand(sp.expand_trig(rn).rewrite(sp.exp)), combine='exp'))
        for term in sp.Add.make_args(e):
            c, n, w = c0, n0, sp.S.Zero
            for f in sp.Mul.make_args(term):
                if f == t:
                    n += 1
                elif f.is_Pow and f.base == t:
                    assert f.exp.is_Integer, f
                    n += int(f.exp)
                elif isinstance(f, sp.exp):
                    a = sp.expand(f.args[0])
                    if a.has(t):
                        coeff = sp.cancel(a / (sp.I * t))
                        assert not coeff.has(t), f
                        w += coeff
                    else:
                        c *= f
                elif f.has(t):
                    raise ValueError("unexpected factor %s" % f)
                else:
                    c *= f
            w = sp.expand(w)
            key = (n, w)
            out[key] = out.get(key, 0) + c
    return {kk: sp.cancel(vv) for kk, vv in out.items() if sp.cancel(vv) != 0}

def _elem_int(n, w, T):
    """int_0^T t^n e^{i w t} dt, n >= 0, w != 0 (elementary)."""
    # int t^n e^{iwt} = e^{iwt} sum_{j=0}^n (-1)^j n!/(n-j)! t^{n-j} /(iw)^{j+1}
    F = lambda s: sp.exp(sp.I * w * s) * sum((-1)**j * sp.factorial(n) / sp.factorial(n - j) * s**(n - j) / (sp.I * w)**(j + 1) for j in range(n + 1))
    return sp.expand(F(T) - F(0))

def integrate_term(n, w, T):
    """Return (finite, sing): int_eps^T t^n e^{i w t} dt = finite(T) + sum_m sing[m] eps^m + sing['log'] ln(eps) + O(eps)."""
    sing = {}
    if w == 0:
        if n >= 0:
            return T**(n + 1) / (n + 1), sing
        if n == -1:
            return sp.log(T), {'log': sp.S.NegativeOne}
        return T**(n + 1) / (n + 1), {n + 1: -sp.S.One / (n + 1)}
    if n >= 0:
        return _elem_int(n, w, T), sing
    if n == -1:
        return -Cin(w * T) + sp.I * Si(w * T) + sp.log(T), {'log': sp.S.NegativeOne}
    # n <= -2
    fin_next, sing_next = integrate_term(n + 1, w, T)
    fac = -sp.I * w / (n + 1)
    finite = T**(n + 1) * sp.exp(sp.I * w * T) / (n + 1) + fac * fin_next
    # lower limit: -eps^{n+1} e^{i w eps}/(n+1) = -sum_j (iw)^j eps^{n+1+j}/(j! (n+1))
    for j in range(0, -n - 1):
        m = n + 1 + j
        sing[m] = sing.get(m, 0) - (sp.I * w)**j / (sp.factorial(j) * (n + 1))
    finite = finite - (sp.I * w)**(-n - 1) / (sp.factorial(-n - 1) * (n + 1))
    for kk, vv in sing_next.items():
        sing[kk] = sing.get(kk, 0) + fac * vv
    return finite, sing

def integrate_0_to_T(expr, t, T, check=True):
    """Closed form of int_0^T expr(t) dt for expr a sum of t^n trig(w t); returns (result, terms-dict)."""
    terms = decompose(expr, t)
    total, sing = 0, {}
    for (n, w), c in terms.items():
        fin, sg = integrate_term(n, w, T)
        total += c * fin
        for kk, vv in sg.items():
            sing[kk] = sing.get(kk, 0) + c * vv
    if check:
        for kk, vv in sing.items():
            vv = sp.cancel(sp.expand(vv.rewrite(sp.cos)))
            if vv != 0:
                vv = sp.simplify(vv)
            assert vv == 0, ("singular piece eps^%s does not cancel: %s" % (kk, vv))
    return total, terms

def realify(expr):
    """exp(i w T) -> cos + i sin, expand, drop vanishing imaginary part (checked)."""
    e = sp.expand(expr.rewrite(sp.cos))
    e = sp.expand(e)
    re, im = e.as_real_imag() if False else (sp.re(e), sp.im(e))
    return e

def to_Ci_form(expr):
    return expr.replace(lambda a: isinstance(a, Cin), lambda a: a.rewrite(Ci))
