"""pw_tools.py -- V3: plane-wave substitution into the independent E_{mu nu} (indep_E.pkl).
Convention f(x) = int d^3k/(2pi)^3 e^{ik.x} f_k  (SPEC).  Discrete modes zeta_m e^{i q_m.x}: the coefficient of
zeta_q zeta_p e^{i(q+p).x} in [eps^2]E (for q != p) equals the symmetric kernel K_S of E^{(2)} = 2[eps^2]E
(two orderings x 1/2 from the 2 = 2!).  Similarly at third order (see check_third.py).
"""
import sympy as sp, pickle, os, random
HERE = os.path.dirname(os.path.abspath(__file__))
t, x, y, z = sp.symbols('t x y z')
X = (t, x, y, z)
I = sp.I
_E = None
def E_mine():
    global _E
    if _E is None:
        raw = pickle.load(open(os.path.join(HERE, 'indep_E.pkl'), 'rb'))
        _E = {kk: sp.parse_expr(vv) for kk, vv in raw.items()}
    return _E

# closed-form linear transfer functions (per zeta), a = eta = t
def xx(kk):
    return kk * t / sp.sqrt(3)
def T_Phi(kk):
    xv = xx(kk)
    return -2 * (sp.sin(xv) - xv * sp.cos(xv)) / xv**3
def T_delta(kk):
    TP = T_Phi(kk)
    return -sp.Rational(2, 3) * t**2 * kk**2 * TP - 2 * t * sp.diff(TP, t) - 2 * TP
def T_V(kk):
    TP = T_Phi(kk)
    return -(t**2 / 2) * (sp.diff(TP, t) + TP / t)

def vec3(kmag, qmag, mu):
    """frame: k along z, q in xz-plane; returns (kvec, qvec, pvec=k-q) as sympy Matrices."""
    s = sp.sqrt(1 - mu**2)
    kv = sp.Matrix([0, 0, kmag]); qv = sp.Matrix([qmag * s, 0, qmag * mu]); pv = kv - qv
    return kv, qv, pv

def mode(zeta, qv, qmag):
    """dict field-name -> contribution of one first-order mode (per zeta, times e^{i q.x})."""
    ph = sp.exp(I * (qv[0] * x + qv[1] * y + qv[2] * z))
    d = {}
    d['Phi1'] = zeta * T_Phi(qmag) * ph
    d['Psi1'] = zeta * T_Phi(qmag) * ph
    d['delta1'] = zeta * T_delta(qmag) * ph
    for i in range(3):
        d['v%d1' % (i + 1)] = zeta * I * qv[i] * T_V(qmag) * ph
    return d

def substitute(expr, fields):
    """replace generic Function('name')(t,x,y,z) by explicit expressions (derivatives evaluated), others -> 0."""
    reps = {}
    for f in expr.atoms(sp.Function):
        if isinstance(f, sp.core.function.AppliedUndef):
            nm = f.func.__name__
            reps[f] = fields.get(nm, sp.Integer(0))
    return expr.subs(reps).doit()

def coeff_at_origin(expr, monomial_syms):
    """coefficient of prod(monomial_syms) (each to first power) at x=y=z=0."""
    e = sp.expand(expr.subs({x: 0, y: 0, z: 0}))
    for s in monomial_syms:
        e = e.coeff(s, 1)
    for s in monomial_syms:
        e = e.subs(s, 0)
    return e

def rand_point(seed=None):
    rng = random.Random(seed)
    return {'k': sp.Rational(rng.randint(20, 300), 100), 'q': sp.Rational(rng.randint(20, 300), 100),
            'mu': sp.Rational(rng.randint(-90, 90), 100), 'eta': sp.Rational(rng.randint(20, 400), 100)}
