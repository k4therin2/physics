"""
fourier.py -- Fourier-space substitution for the C1 engine.
Convention: f(x) = int d^3k/(2pi)^3 e^{ik.x} f_k, so d_j -> i k_j.
Frame: k = (0,0,k);  q = (q sqrt(1-mu^2), 0, q mu);  p = k - q, |p| = p = sqrt(k^2+q^2-2 k q mu).
First-order fields: X1_k = zeta_k T_X(k, eta) for X in Phi, Psi, delta;  v1^i = d_i V1, V1_k = zeta_k T_V(k, eta).
Order-n unknown amplitudes: Function('Xn')(eta) (scalars), vn_x, vn_y, vn_z, omn_x, omn_y, hn_p, hn_x.
"""
import sympy as sp
import itertools, math
from tpoly import REG, Poly
from fractions import Fraction

eta = sp.Symbol('eta', positive=True)
k, q, p = sp.symbols('k q p', positive=True)
mu = sp.Symbol('mu', real=True)
I = sp.I

class Mom:
    def __init__(self, name, vec, mag):
        self.name, self.vec, self.mag = name, list(vec), mag
    def neg(self):
        return Mom('-' + self.name, [-c for c in self.vec], self.mag)

KV = Mom('k', [0, 0, k], k)
QV = Mom('q', [q * sp.sqrt(1 - mu**2), 0, q * mu], q)
PV = Mom('p', [-q * sp.sqrt(1 - mu**2), 0, k - q * mu], p)      # p = k - q
ZERO = Mom('0', [0, 0, 0], 0)

TF = {n: sp.Function('T_' + n) for n in ('Phi', 'Psi', 'delta', 'V')}

def jet_first(field, m, mom):
    """Fourier image of d^m X1 with momentum mom (per unit zeta)."""
    if field.startswith('v'):
        comp = int(field[1]) - 1
        base = I * mom.vec[comp] * TF['V'](mom.mag, eta)
    elif field in ('Phi', 'Psi', 'delta'):
        base = TF[field](mom.mag, eta)
    else:
        raise ValueError("no first-order %s" % field)
    if m[0]:
        base = sp.diff(base, eta, m[0])
    for i in range(3):
        if m[i + 1]:
            base = base * (I * mom.vec[i]) ** m[i + 1]
    return base

def amp(field, n):
    """Order-n unknown amplitude in the k-frame (k along z)."""
    if field in ('Phi', 'Psi', 'delta'):
        return sp.Function('%s%d' % (field, n))(eta)
    if field.startswith('v'):
        return sp.Function('v%d_%s' % (n, 'xyz'[int(field[1]) - 1]))(eta)
    if field.startswith('om'):
        c = 'xyz'[int(field[2]) - 1]
        return sp.Function('om%d_%s' % (n, c))(eta) if c != 'z' else sp.Integer(0)
    if field.startswith('h'):
        i, j = field[1], field[2]
        if (i, j) == ('1', '1'): return sp.Function('h%d_p' % n)(eta)
        if (i, j) == ('2', '2'): return -sp.Function('h%d_p' % n)(eta)
        if (i, j) == ('1', '2'): return sp.Function('h%d_x' % n)(eta)
        return sp.Integer(0)      # h13, h23, h33 vanish for TT with k along z
    raise ValueError(field)

def jet_amp(field, n, m, mom=KV):
    base = amp(field, n)
    if m[0]:
        base = sp.diff(base, eta, m[0])
    for i in range(3):
        if m[i + 1]:
            base = base * (I * mom.vec[i]) ** m[i + 1]
    return base

def monomial_kernel(jets, moms, lower):
    """Symmetrised kernel for product of jets with momenta moms (same length). lower(jet_id, mom)->expr."""
    r = len(jets)
    tot = 0
    for perm in itertools.permutations(range(r)):
        t = 1
        for i, j in enumerate(jets):
            t = t * lower(j, moms[perm[i]])
        tot += t
    return tot / math.factorial(r)

def linear_part(poly, n):
    """Fourier image of the terms linear in order-n jets (amplitudes in k-frame)."""
    tot = 0
    for (o, e, mono), c in poly.d.items():
        if o != n or len(mono) != 1:
            continue
        f, od, m = REG.info[mono[0]]
        assert od == n
        tot += sp.Rational(c.numerator, c.denominator) * eta**e * jet_amp(f, n, m)
    return sp.expand(tot)

def quadratic_kernel(poly, moms=(QV, PV)):
    """Kernel K(q,p) of the terms quadratic in first-order jets (per zeta_q zeta_p)."""
    tot = 0
    for (o, e, mono), c in poly.d.items():
        if len(mono) != 2 or any(REG.order(j) != 1 for j in mono):
            continue
        def lower(j, mom):
            f, od, m = REG.info[j]
            return jet_first(f, m, mom)
        tot += sp.Rational(c.numerator, c.denominator) * eta**e * monomial_kernel(mono, moms, lower)
    return tot

def simplify_kernel(expr):
    """Replace sqrt(1-mu^2)^2 and express via k,q,p: mu = (k^2+q^2-p^2)/(2kq)."""
    expr = sp.expand(expr)
    return expr
