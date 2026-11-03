import sympy as sp
from zm_invariance import *
a = sp.Symbol('a')
out = run('time rescaling eta_old = (1 + eps s a) eta', (1 + eps*s*a)*eta, [x, y, z], verbose=False)
L1ray = sp.sympify(sum(e for w, e in LG[1] if w == 'ray'))
verify(out, [eta*L1ray], 'time rescaling')
