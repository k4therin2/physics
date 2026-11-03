import random, sympy as sp
from pw_tools import *
from third_tools import *
E = E_mine()
z1, z3 = sp.symbols('zeta1 zeta3')
pt = rand_point(0); rng = random.Random(7)
kv, qv, pv = vec3(pt['k'], pt['q'], pt['mu']); qm = pt['q']; km = pt['k']
base = {}
for nm, v_ in mode(z1, qv, qm).items(): base[nm] = base.get(nm, 0) + v_
for nm, v_ in mode(z3, kv, km).items(): base[nm] = base.get(nm, 0) + v_
pf = pair_fields(z1, z3, qv, kv, qm, km, amp=sp.Rational(1, 2))
def resid(names, comp):
    f = dict(base)
    for nm in names: f[nm] = pf[nm]
    r = coeff_at_origin(substitute(E[('E',) + comp + (2,)], f), (z1, z3))
    return ev(r, rng, pt['eta'], 15)
scal = ['Phi2', 'Psi2', 'delta2']; vel = ['v12', 'v22', 'v32']; om = ['om12', 'om22', 'om32']; hh = ['h112', 'h122', 'h132', 'h222', 'h232', 'h332']
for label, names in [('scalars only', scal), ('scalars+v', scal + vel), ('scalars+v+om', scal + vel + om), ('all', scal + vel + om + hh), ('h only', hh), ('om only', om)]:
    print("%-14s 00: %s   0x: %s   0z: %s   TL: %s" % (label, sp.N(resid(names, (0, 0)), 8), sp.N(resid(names, (0, 1)), 8), sp.N(resid(names, (0, 3)), 8),
          sp.N(resid(names, (3, 3)) - (resid(names, (1, 1)) + resid(names, (2, 2))) / 2, 8)))
# TT / transversality checks of the pair fields
Kv = qv + kv
h = sp.Matrix(3, 3, lambda i, j: pf[{(0,0):'h112',(0,1):'h122',(0,2):'h132',(1,1):'h222',(1,2):'h232',(2,2):'h332'}[(min(i,j),max(i,j))]])
print("trace h2:", sp.simplify(h.trace()), " K.h:", [sp.simplify(sum(Kv[i]*h[i,j] for i in range(3))) for j in range(3)])
print("K.om2:", sp.simplify(sum(Kv[i]*pf['om%d2'%(i+1)] for i in range(3))))
