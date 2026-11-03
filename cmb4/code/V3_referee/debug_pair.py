import random, sympy as sp
from pw_tools import *
from third_tools import *
E = E_mine()
z1, z2, z3 = sp.symbols('zeta1 zeta2 zeta3')
pt = rand_point(0); rng = random.Random(7)
kv, qv, pv = vec3(pt['k'], pt['q'], pt['mu']); qm = pt['q']; km = pt['k']
Kv = qv + kv; Km = sp.sqrt(sum(c**2 for c in Kv))
# first-order modes 1 (q) and 3 (k) only, no second-order fields: source for the pair (q,k)
f = {}
for nm, v_ in mode(z1, qv, qm).items(): f[nm] = f.get(nm, 0) + v_
for nm, v_ in mode(z3, kv, km).items(): f[nm] = f.get(nm, 0) + v_
S00_disc = coeff_at_origin(substitute(E[('E', 0, 0, 2)], f), (z1, z3))
S0_disc = [coeff_at_origin(substitute(E[('E', 0, i, 2)], f), (z1, z3)) for i in (1, 2, 3)]
N = lambda e: sp.N(e.subs(t, pt['eta']), 20)
print("discrete S00 for pair (q,k):", N(S00_disc), "   C1 K_S00 at (q,k,|q+k|):", N(at(S2['K_S00'], qm, km, Km)))
Khat = Kv / Km
SL = sum(Khat[i] * S0_disc[i] for i in range(3))
print("longitudinal S0 (Khat.S0):", N(SL), "   C1 K_S0z at (q,k,|q+k|):", N(at(S2['K_S0z'], qm, km, Km)))
# also in the standard frame (pair (q, k-q) with total k along z) for reference
f2 = {}
pm = sp.sqrt(km**2 + qm**2 - 2*km*qm*pt['mu'])
for nm, v_ in mode(z1, qv, qm).items(): f2[nm] = f2.get(nm, 0) + v_
for nm, v_ in mode(z3, pv, pm).items(): f2[nm] = f2.get(nm, 0) + v_
print("frame check S00 (q,k-q):", N(coeff_at_origin(substitute(E[('E',0,0,2)], f2), (z1, z3))), " C1:", N(at(S2['K_S00'], qm, pm, km)))
print("at() sanity: K_S00 symbols ->", at(S2['K_S00'], sp.Symbol('A'), sp.Symbol('B'), sp.Symbol('C')).free_symbols)
