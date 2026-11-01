"""kernels_api.py -- convenience loader for the C1 results.
Usage:  from kernels_api import *;  K = load();  K['second']['K_h'], K['third']['SPhi'], ...
closed_form(expr) inserts the closed-form T_Phi; invariant kernels use symbols k,q,p = |q1+q2|,|q1|,|q2|.
"""
import sympy as sp, pickle, os
from fourier import eta, k, q, p, mu, TF, I
HERE = os.path.dirname(os.path.abspath(__file__))
def _load(fn):
    return {kk: (sp.parse_expr(vv) if isinstance(vv, str) else vv) for kk, vv in pickle.load(open(os.path.join(HERE, fn), 'rb')).items()}
def load():
    K = {'first': _load('first_order.pkl'), 'second': _load('second_order.pkl'), 'third': _load('third_order.pkl')}
    ic = pickle.load(open(os.path.join(HERE, 'third_order_ic.pkl'), 'rb'))
    K['third_ic'] = {kk: (sp.parse_expr(vv) if isinstance(vv, str) else {k2: sp.parse_expr(v2) for k2, v2 in vv.items()}) for kk, vv in ic.items()}
    return K
X = lambda kk: kk * eta / sp.sqrt(3)
T_Phi_closed = lambda kk: -2 * (sp.sin(X(kk)) - X(kk) * sp.cos(X(kk))) / X(kk)**3
def closed_form(expr, syms=(q, p, k)):
    expr = expr.doit()
    for kk in syms:
        expr = expr.subs(TF['Phi'](kk, eta), T_Phi_closed(kk))
    return expr.doit()
def to_invariants(expr):
    """mu-form -> magnitudes (k,q,p)"""
    return sp.expand(expr).subs(mu, (k**2 + q**2 - p**2) / (2 * k * q))
