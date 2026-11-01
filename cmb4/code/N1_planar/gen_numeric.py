"""Turn the sympy expressions (eqs_generated.pkl, ics_generated_N.pkl) into numpy-callable functions.
Also rewrites the background-cancelling combinations in terms of dalh = alh - 1 (exact identities, verified).
"""
import sympy as sp, pickle, numpy as np

d = pickle.load(open('eqs_generated.pkl', 'rb'))
S = d['symbols']
t, Ah, Bh, alh, Ax, Bx, Axx, Bxx, pA, pB, pAx, pBx, ze, zx, zxx = [S[k] for k in
    ('t', 'Ah', 'Bh', 'alh', 'Ax', 'Bx', 'Axx', 'Bxx', 'pA', 'pB', 'pAx', 'pBx', 'zeta', 'zx', 'zxx')]
dalh = sp.Symbol('dalh', real=True)   # = alh - 1, computed with expm1 to avoid cancellation at small t

# exact rewrites (verified below):  1/(2 alh t^2) - 1/(2 alh^3 t^2) = dalh (2 + dalh) / (2 alh^3 t^2)
bg = 1/(2*alh*t**2) - 1/(2*alh**3*t**2)
bg_new = dalh*(2 + dalh)/(2*alh**3*t**2)
assert sp.simplify((bg - bg_new).subs(dalh, alh - 1)) == 0
pA_t = sp.expand(d['pA_t']) - bg + bg_new
pB_t = sp.expand(d['pB_t']) - bg + bg_new
assert sp.simplify((pA_t - d['pA_t']).subs(dalh, alh - 1)) == 0
assert sp.simplify((pB_t - d['pB_t']).subs(dalh, alh - 1)) == 0
Ham_new = (3*dalh*(2 + dalh) + 2*alh**3*pB*t + alh**4*pB**2*t**2)/(alh**4*t**4) + 2*(pA + pB)/(alh*t**3) + 2*pA*pB/t**2 \
    - 2*Bxx/(Ah**2*Bh*t**2) - Bx**2/(Ah**2*Bh**2*t**2) + 2*Ax*Bx/(Ah**3*Bh*t**2)
assert sp.simplify((Ham_new - d['Ham']).subs(dalh, alh - 1)) == 0
Mom = d['Mom']
args = (t, Ah, Bh, alh, dalh, Ax, Bx, Axx, Bxx, pA, pB, pAx, pBx, ze, zx, zxx)
f_pA_t = sp.lambdify(args, pA_t, 'numpy')
f_pB_t = sp.lambdify(args, pB_t, 'numpy')
f_Ham = sp.lambdify(args, Ham_new, 'numpy')
f_Mom = sp.lambdify(args, Mom, 'numpy')
# size of the largest individual term of Ham (for relative constraint diagnostics)
f_Ham_terms = sp.lambdify(args, [sp.Abs(term) for term in sp.Add.make_args(sp.expand(Ham_new))], 'numpy')
f_Mom_terms = sp.lambdify(args, [sp.Abs(term) for term in sp.Add.make_args(sp.expand(Mom))], 'numpy')

def load_ics(order):
    """returns (order, list of (a_n, b_n) numpy functions of (z0..z_order)) from ics_generated_{order}.pkl"""
    p = pickle.load(open(f'ics_generated_{order}.pkl', 'rb'))
    x = sp.Symbol('x', real=True)
    zeta = sp.Function('zeta')(x)
    zs = sp.symbols(f'z0:{order + 1}')
    rep = [(sp.Derivative(zeta, (x, n)), zs[n]) for n in range(order, 0, -1)] + [(zeta, zs[0])]
    def conv(e):
        e = sp.sympify(e)
        for a, b in rep:
            e = e.subs(a, b)
        return sp.lambdify(zs, e, 'numpy')
    funcs = {}
    for n in range(1, order + 1):
        funcs[n] = (conv(p['sol'][f'a{n}(x)']), conv(p['sol'][f'b{n}(x)']))
    return funcs
