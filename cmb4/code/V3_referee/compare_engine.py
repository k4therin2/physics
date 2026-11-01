"""compare_engine.py -- V3: term-by-term comparison of C1's E_{mu nu} (E.pkl) with the independent sympy computation."""
import sys, os, pickle, sympy as sp
C1 = '/home/user/physics/cmb4/code/C1_einstein'
sys.path.insert(0, C1)
from tpoly import load_polys
d = load_polys(os.path.join(C1, 'E.pkl'))
mine = pickle.load(open('indep_E.pkl', 'rb'))
t, x, y, z = sp.symbols('t x y z')
def funcs(f, od, m):
    ex = sp.Function('%s%d' % (f, od))(t, x, y, z)
    for var, n in zip((t, x, y, z), m):
        if n: ex = sp.Derivative(ex, (var, n))
    return ex
ok = True
for a in range(4):
    for b in range(a, 4):
        for n in (1, 2, 3):
            c1 = sp.expand(d['E'][a][b].order_part(n).to_sympy(eta=t, funcs=funcs))
            me = sp.expand(sp.parse_expr(mine[('E', a, b, n)]))
            diff = sp.expand(c1 - me)
            print("E_%d%d order %d: C1 terms %d, V3 terms %d, difference terms: %d" % (a, b, n, len(sp.Add.make_args(c1)), len(sp.Add.make_args(me)), 0 if diff == 0 else len(sp.Add.make_args(diff))))
            if diff != 0:
                ok = False; print("   DIFF:", diff)
for a in range(4):
    for n in (1, 2, 3):
        c1 = sp.expand(d['u_dn'][a].order_part(n).to_sympy(eta=t, funcs=funcs))
        me = sp.expand(sp.parse_expr(mine[('u_dn', a, n)]))
        if sp.expand(c1 - me) != 0:
            ok = False; print("u_dn", a, n, "DIFF", sp.expand(c1 - me))
print("ENGINE AGREEMENT (all components, orders 1-3, and u_mu):", ok)
