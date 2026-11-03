"""check_ops.py -- V3: claim (2) linear operators at order n (n=2 and 3 identical), from the independent E,
with transversality / TT imposed for k along z: om3 = 0, h13 = h23 = h33 = 0, h22 = -h11, h12 = hx."""
import sympy as sp
from pw_tools import *
E = E_mine()
kk = sp.Symbol('k', positive=True); ph = sp.exp(I * kk * z)
A = lambda nm: sp.Function('A_' + nm)(t)
for n in (2, 3):
    nm = lambda f: '%s%d' % (f, n)
    amps = {nm(f): A(f) * ph for f in ['Phi', 'Psi', 'delta', 'v1', 'v2', 'v3', 'om1', 'om2']}
    amps[nm('h11')] = A('hp') * ph; amps[nm('h22')] = -A('hp') * ph; amps[nm('h12')] = A('hx') * ph
    L = lambda a, b: sp.expand((sp.factorial(n) * substitute(E[('E', a, b, n)], amps)).subs({x: 0, y: 0, z: 0}))
    L00, L0x, L0z, Lxx, Lyy, Lzz, Lxz, Lxy = [L(*ab) for ab in [(0, 0), (0, 1), (0, 3), (1, 1), (2, 2), (3, 3), (1, 3), (1, 2)]]
    S = {s: sp.Symbol(s) for s in ['S00', 'S0x', 'S0z', 'Sxx', 'Syy', 'Szz', 'Sxz', 'Sxy']}
    ok = {}
    ok['TL: -k^2(Phi-Psi) + S_TL = 0'] = sp.simplify(Lzz - (Lxx + Lyy) / 2 + kk**2 * (A('Phi') - A('Psi'))) == 0
    de = sp.solve(L00 + S['S00'], A('delta'))[0]
    ok['00: delta = eta^2 S00/3 - ...'] = sp.simplify(de - (t**2 / 3 * S['S00'] - sp.Rational(2, 3) * t**2 * kk**2 * A('Phi') - 2 * t * A('Phi').diff(t) - 2 * A('Psi'))) == 0
    vz = sp.solve(L0z + S['S0z'], A('v3'))[0]
    ok['0z: V = -(eta^2/4)[2(Phi\'+Psi/eta) + S0z/(ik)]'] = sp.simplify(vz - I * kk * (-(t**2 / 4) * (2 * (A('Phi').diff(t) + A('Psi') / t) + S['S0z'] / (I * kk)))) == 0
    ok['0x: k^2 om/2 + 4(om+v)/eta^2 (=-S0x)'] = sp.simplify(L0x - (kk**2 * A('om1') / 2 + 4 * (A('om1') + A('v1')) / t**2)) == 0
    ok['xz: -(ik/2)(om\'+2om/eta)'] = sp.simplify(Lxz - (-(I * kk / 2) * (A('om1').diff(t) + 2 * A('om1') / t))) == 0
    ok['xx-yy: hp\'\'+2hp\'/eta+k^2hp'] = sp.simplify(Lxx - Lyy - (A('hp').diff(t, 2) + 2 * A('hp').diff(t) / t + kk**2 * A('hp'))) == 0
    # trace -> Phi ODE: (L_tr + S_tr) = 0 with delta from 00 and Psi = Phi - Q, Q = S_TL/k^2
    Q = sp.Function('Q')(t)
    Ltr = Lxx + Lyy + Lzz
    tr = (Ltr + sp.Symbol('Str')).subs(A('delta'), de).doit().subs(A('Psi'), A('Phi') - Q).doit()
    c = sp.simplify(tr.coeff(A('Phi').diff(t, 2)))
    lhs = sp.expand(tr / c)
    target = A('Phi').diff(t, 2) + 4 / t * A('Phi').diff(t) + kk**2 / 3 * A('Phi') - (-(sp.Symbol('Str') - S['S00']) / 6 + Q.diff(t) / t - kk**2 / 3 * Q)
    ok['trace: Phi\'\'+4Phi\'/eta+k^2Phi/3 = -(Str-S00)/6 + Q\'/eta - k^2Q/3  (norm=%s)' % c] = sp.simplify(lhs - target) == 0
    print("order", n)
    for kx, v_ in ok.items(): print("   %-75s %s" % (kx, v_))
