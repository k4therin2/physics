"""Control: the V5 rebuild (first-order omega,h allowed) reduces to C2's L1, L2/2 (shear form) when omega1 = h1 = 0."""
import sympy as sp, pickle, sys
sys.path.insert(0, '.')
from fields import F, X
import ray
V5 = pickle.load(open('L_v5.pkl', 'rb'))
C2 = pickle.load(open('../C2_observable/L_terms.pkl', 'rb'))
zero = {}
for i in (1, 2, 3):
    zero[F('om_%d_1' % i)] = 0
    for j in range(i, 4):
        if (i, j) != (3, 3): zero[F('h_%d%d_1' % (i, j))] = 0
for m in (1, 2):
    for where in ('ray', 'obs', 'emi'):
        a = sp.expand(sp.sympify(sum(e for w, e in V5["L"][m] if w == where)).xreplace(zero).doit())
        b = sp.expand(sum(t['expr'] for t in C2['L'][m] if t['where'] == where))
        print('order', m, where, 'V5 - C2 =', ray.unit_rules(sp.expand(a - b)))
