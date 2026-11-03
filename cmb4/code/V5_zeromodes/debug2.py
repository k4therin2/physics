import sympy as sp
from zm_invariance import *
lam = sp.symbols('lam1:4')
fs = {}
for i in range(3):
    fs[kf[1][i]] = kf[1][i] + lam[i]*(eta_o - eta); fs[kf[1][i].subs(eta, eta_s)] = kf[1][i].subs(eta, eta_s) + lam[i]*(eta_o - eta_s)
    fs[Dxf[1][i]] = Dxf[1][i] - lam[i]*(eta_o - eta)**2/2; fs[Dxf[1][i].subs(eta, eta_s)] = Dxf[1][i].subs(eta, eta_s) - lam[i]*(eta_o - eta_s)**2/2
for where in ('ray', 'obs', 'emi'):
    old = sp.sympify(sum(e for w, e in LG[2] if w == where))
    print(where, 'functionals present:', sorted(str(a) for a in old.atoms(sp.core.function.AppliedUndef) if str(a.func).startswith(('k1', 'Dx1'))))
    d = unit_rules(sp.expand(old.xreplace(fs) - old))
    print(where, 'delta terms:', len(sp.Add.make_args(d)))
    if where == 'ray': print(sp.collect(d, lam)) if len(str(d)) < 3000 else print(str(d)[:3000])
