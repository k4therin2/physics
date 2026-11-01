"""check_second.py -- V3: second-order linear operators and source kernels from the independent E_{mu nu},
via explicit plane waves, compared numerically with C1's kernels (K_S00, K_S0z, K_Q, K_SPhi, K_h, K_om, K_vT)."""
import sys, os, sympy as sp, mpmath
from pw_tools import *
C1 = '/home/user/physics/cmb4/code/C1_einstein'
sys.path.insert(0, C1)
import kernels_api as KA
K = KA.load()
from fourier import eta, k, q, p, mu, TF
E = E_mine()
eps = sp.Symbol('eps')
def E2(a, b):
    return E[('E', min(a, b), max(a, b), 2)]
PREC = 30

# ---------------- 1. linear operators (second-order amplitudes x e^{ikz}, first order = 0)
kk = sp.Symbol('k', positive=True)
amps = {}
names = ['Phi2', 'Psi2', 'delta2', 'v12', 'v22', 'v32', 'om12', 'om22', 'om32', 'h112', 'h122', 'h132', 'h222', 'h232', 'h332']
ph = sp.exp(I * kk * z)
for nm in names:
    amps[nm] = sp.Function('A_' + nm)(t) * ph
def lin(a, b):
    return sp.expand((2 * substitute(E2(a, b), amps)).subs({x: 0, y: 0, z: 0}))
L00, L03, L01, L11, L22, L33, L13, L12 = [lin(*ab) for ab in [(0, 0), (0, 3), (0, 1), (1, 1), (2, 2), (3, 3), (1, 3), (1, 2)]]
A = lambda nm: sp.Function('A_' + nm)(t)
print("L00 =", L00)
print("L0z =", L03)
print("L0x =", L01)
print("L_TL = Lzz-(Lxx+Lyy)/2 =", sp.simplify(L33 - (L11 + L22) / 2))
print("L_xz =", L13)
print("L_xx - L_yy =", sp.expand(L11 - L22))
Ltr = sp.expand(L11 + L22 + L33)
# claimed operator forms (C1 report)
Q2 = sp.Function('Q2')(t)
claim_TL = -kk**2 * (A('Phi2') - A('Psi2'))
claim_00 = -3 / t**2 * (A('delta2') + sp.Rational(2, 3) * t**2 * kk**2 * A('Phi2') + 2 * t * A('Phi2').diff(t) + 2 * A('Psi2'))
# from report: delta_n = (eta^2/3) S00 - ... means L00 + S00 = 0 with L00 = -(3/eta^2) delta - ... ; we check the solve below instead
print("TL claim ok:", sp.simplify(L33 - (L11 + L22) / 2 - claim_TL) == 0)
sol_de = sp.solve(L00 + sp.Symbol('S00'), A('delta2'))[0]
print("delta2 from 00:", sol_de, "  matches C1:", sp.simplify(sol_de - (t**2 / 3 * sp.Symbol('S00') - sp.Rational(2, 3) * t**2 * kk**2 * A('Phi2') - 2 * t * A('Phi2').diff(t) - 2 * A('Psi2'))) == 0)
sol_vz = sp.solve(L03 + sp.Symbol('S0z'), A('v32'))[0]
claim_vz = I * kk * (-(t**2 / 4) * (2 * (A('Phi2').diff(t) + A('Psi2') / t) + sp.Symbol('S0z') / (I * kk)))
print("v2z from 0z:", sol_vz, "  matches C1 V_n formula:", sp.simplify(sol_vz - claim_vz) == 0)
# 0x transverse: (1/2) lap om - (4/eta^2)(om + vT) = S0x^T  ->  -k^2 om/2 - 4(om+v)/eta^2 + S0x = 0 ?
claim_0x = -kk**2 * A('om12') / 2 - 4 * (A('om12') + A('v12')) / t**2
print("0x operator matches C1 (-k^2 om/2 - 4(om+v)/eta^2):", sp.simplify(L01 - claim_0x) == 0)
# ij vector: om' + 2 om/eta = (2/(ik)) Sxz  <=>  L13 = -(ik/2)(om' + 2om/eta)
claim_xz = -(I * kk / 2) * (A('om12').diff(t) + 2 * A('om12') / t)
print("xz operator matches C1 (-(ik/2)(om'+2om/eta)):", sp.simplify(L13 - claim_xz) == 0)
# tensor: with h11 = -h22 = hp: (hxx - hyy)'' + 2(..)'/eta + k^2(..) = -2(Sxx - Syy)  <=> Lxx-Lyy = hp'' + 2hp'/eta + k^2 hp
hp = sp.Function('hp')(t)
LT = sp.expand((L11 - L22).subs({A('h112'): hp, A('h222'): -hp}).doit())
print("tensor operator (h11=-h22=hp):", LT, "  matches hp''+2hp'/eta+k^2 hp:", sp.simplify(LT - (hp.diff(t, 2) + 2 * hp.diff(t) / t + kk**2 * hp)) == 0)
# trace -> Phi ODE
odeL = sp.expand(((Ltr - L00) / 6).subs(A('Psi2'), A('Phi2') - Q2).doit())
odeL = sp.expand(odeL.subs(A('delta2'), sol_de.subs(A('Psi2'), A('Phi2') - Q2)).doit())
target = A('Phi2').diff(t, 2) + 4 / t * A('Phi2').diff(t) + kk**2 / 3 * A('Phi2') - Q2.diff(t) / t + kk**2 / 3 * Q2
print("(trace-00)/6 after eliminating delta2 (with S00 -> 0):", sp.simplify(odeL.subs(sp.Symbol('S00'), 0) - target) == 0,
      "; S00 coefficient:", sp.simplify(odeL.coeff(sp.Symbol('S00'))))
# NOTE: eliminating delta2 via the 00 equation brings S00 into the trace equation: S_Phi = -(S_tr - S_00)/6 ... check the sign/weight
# full: (Ltr + Str - L00 - S00)/6 = 0 with L00 = -S00 - (3/t^2)(delta2 + ...) => contribution of S00
print("   => with delta2 eliminated, the Phi ODE source is -(S_tr - S00)/6 + Q'/eta - k^2 Q/3 iff coefficient above is", -sp.Rational(1, 6) * 0, "(S00 enters only through delta2 elimination)")

# ---------------- 2. sources from two plane waves q, p (= k - q)
zq, zp = sp.symbols('zeta_q zeta_p')
def sources(pt):
    kv, qv, pv = vec3(pt['k'], pt['q'], pt['mu'])
    qm = pt['q']; pm = sp.sqrt(pt['k']**2 + pt['q']**2 - 2 * pt['k'] * pt['q'] * pt['mu'])
    f = {}
    for nm, val in mode(zq, qv, qm).items():
        f[nm] = f.get(nm, 0) + val
    for nm, val in mode(zp, pv, pm).items():
        f[nm] = f.get(nm, 0) + val
    S = {}
    for ab in [(0, 0), (0, 1), (0, 2), (0, 3), (1, 1), (2, 2), (3, 3), (1, 2), (1, 3), (2, 3)]:
        ex = substitute(E2(*ab), f)
        S[ab] = coeff_at_origin(ex, (zq, zp))       # still a function of t
    return S, pm
def num(e, pt):
    return sp.N(e.subs(t, pt['eta']), PREC)
def c1(expr, pt, pm):
    e = KA.closed_form(expr)
    return sp.N(e.subs({k: pt['k'], q: pt['q'], p: pm, mu: pt['mu'], eta: pt['eta']}), PREC)

worst = 0
for seed in range(4):
    pt = rand_point(seed)
    S, pm = sources(pt)
    kval = pt['k']; qperp = pt['q'] * sp.sqrt(1 - pt['mu']**2)
    res = {}
    S00 = S[(0, 0)]; S0z = S[(0, 3)]; S0x = S[(0, 1)]; Sxz = S[(1, 3)]
    STL = S[(3, 3)] - (S[(1, 1)] + S[(2, 2)]) / 2
    Str = S[(1, 1)] + S[(2, 2)] + S[(3, 3)]
    Qmine = STL / kval**2
    SPhi = -(Str - S00) / 6 + sp.diff(Qmine, t) / t - kval**2 / 3 * Qmine
    Kh = -2 * (S[(1, 1)] - S[(2, 2)]) / qperp**2
    res['K_S00'] = (num(S00, pt), c1(K['second']['K_S00'], pt, pm))
    res['K_S0z'] = (num(S0z, pt), c1(K['second']['K_S0z'], pt, pm))
    res['K_Q'] = (num(Qmine, pt), c1(K['second']['K_Q'], pt, pm))
    res['K_SPhi'] = (num(SPhi, pt), c1(K['second']['K_SPhi'], pt, pm))
    res['K_h'] = (num(Kh, pt), c1(K['second']['K_h'], pt, pm))
    # vanishing components
    res['S0y,Syz,Sxy=0'] = (num(S[(0, 2)], pt) + num(S[(2, 3)], pt) + num(S[(1, 2)], pt), 0)
    # omega: algebraic K_om must satisfy my xz equation: -(ik/2)(om' + 2 om/eta) + Sxz = 0 with om_x = qperp K_om
    Kom_c = KA.closed_form(K['second']['K_om']).subs({k: pt['k'], q: pt['q'], p: pm, mu: pt['mu']}).subs(eta, t)
    omx = qperp * Kom_c
    resid_xz = -(I * kval / 2) * (sp.diff(omx, t) + 2 * omx / t) + Sxz
    res['om ODE residual'] = (num(resid_xz, pt), 0)
    # transverse velocity from 0x: -k^2 om/2 - 4(om + v)/eta^2 + S0x = 0
    KvT_c = KA.closed_form(K['second']['K_vT']).subs({k: pt['k'], q: pt['q'], p: pm, mu: pt['mu']}).subs(eta, t)
    vx = qperp * KvT_c
    resid_0x = kval**2 * omx / 2 + 4 * (omx + vx) / t**2 + S0x   # Lin + S = 0, Lin from my E
    res['vT 0x residual'] = (num(resid_0x, pt), 0)
    # vorticity-free check: curl[(1+delta_T) uhat_i]^(2) = 0  with uhat_i = ghat_{i mu} uhat^mu:  transverse part of
    # (uhat_i)^(2) + delta_T1 v1_i  (discrete coefficient) must vanish: x-component (k along z): om2_x + v2_x + [delta_T1 v1_x]_{zq zp}
    dT_v = None
    kv, qv, pv = vec3(pt['k'], pt['q'], pt['mu'])
    f = {}
    for nm, val in mode(zq, qv, pt['q']).items(): f[nm] = f.get(nm, 0) + val
    for nm, val in mode(zp, pv, pm).items(): f[nm] = f.get(nm, 0) + val
    # (1+delta)^{1/4} uhat_x to second order: uhat_x = u_x / t ; u_dn from my engine order 2 includes om + v2 + quadratic
    u_x2 = sp.parse_expr(__import__('pickle').load(open('indep_E.pkl', 'rb'))[('u_dn', 1, 2)])   # coefficient of eps^2 of u_x
    u_x1 = sp.parse_expr(__import__('pickle').load(open('indep_E.pkl', 'rb'))[('u_dn', 1, 1)])
    # (1+delta)^{1/4} = 1 + eps delta1/4 + ...; product at eps^2: u_x2 + (delta1/4) u_x1
    quad = coeff_at_origin(substitute(2 * u_x2 + 2 * sp.Rational(1, 4) * f['delta1'] * u_x1, f), (zq, zp)) / t   # /t: uhat = u/a, x2 convention
    # linear part of 2*[eps^2](1+dT)uhat_x = (om2_x + v2_x) (t^2 * t^{-1}... ) : check normalisation from my engine
    linu = sp.expand(substitute(2 * u_x2, amps).subs({x: 0, y: 0, z: 0}) / t)
    res['vorticity (om2+v2 + quad) residual'] = (num(sp.expand(linu.subs({A('om12'): omx, A('v12'): vx}).doit() + quad / 2), pt), 0)
    print("point", seed, {kk2: str(v[0]) for kk2, v in pt.items()} if False else pt)
    for nm, (a, b) in res.items():
        d = abs(sp.N(a - b, PREC)); worst = max(worst, d)
        print("   %-40s mine=%s  C1=%s  |diff|=%.2e" % (nm, sp.N(a, 12), sp.N(b, 12), d))
print("linear part of 2[eps^2] uhat_x:", linu)
print("WORST DISCREPANCY over all second-order kernels:", worst)
