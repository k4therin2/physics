"""check_vector.py -- V3: claim (3).  From the independent E^{(2)}: 0x constraint + ij-vector ODE + Kelvin
curl[(1+delta_T) uhat_i] = 0, solved by me; compared with C1's K_om, K_vT.  Also PUB (A.17-18) and Lu et al. (20)."""
import sys, pickle, sympy as sp
from pw_tools import *
C1 = '/home/user/physics/cmb4/code/C1_einstein'
sys.path.insert(0, C1)
import kernels_api as KA
K = KA.load()
from fourier import eta, k, q, p, mu
E = E_mine()
raw = pickle.load(open('indep_E.pkl', 'rb'))
u_x2 = sp.parse_expr(raw[('u_dn', 1, 2)]); u_x1 = sp.parse_expr(raw[('u_dn', 1, 1)])
PREC = 30
zq, zp = sp.symbols('zeta_q zeta_p')
Aom, Av = sp.Function('A_om12')(t), sp.Function('A_v12')(t)
kk = sp.Symbol('k', positive=True)
worst = 0
for seed in range(4):
    pt = rand_point(seed)
    kv, qv, pv = vec3(pt['k'], pt['q'], pt['mu'])
    pm = sp.sqrt(pt['k']**2 + pt['q']**2 - 2 * pt['k'] * pt['q'] * pt['mu'])
    f = {}
    for nm, val in mode(zq, qv, pt['q']).items(): f[nm] = f.get(nm, 0) + val
    for nm, val in mode(zp, pv, pm).items(): f[nm] = f.get(nm, 0) + val
    kval = pt['k']; qperp = pt['q'] * sp.sqrt(1 - pt['mu']**2)
    S0x = coeff_at_origin(substitute(E[('E', 0, 1, 2)], f), (zq, zp))
    Sxz = coeff_at_origin(substitute(E[('E', 1, 3, 2)], f), (zq, zp))
    # quadratic part of 2[eps^2] of (1+delta)^{1/4} uhat_x, discrete coefficient of zeta_q zeta_p
    quad = coeff_at_origin(substitute(2 * u_x2 + sp.Rational(1, 2) * f['delta1'] * u_x1, f), (zq, zp)) / t
    # linear operators (my E): 0x: k^2 om/2 + 4(om+v)/t^2 + S0x = 0 ; Kelvin: om + v + quad = 0
    om_s, v_s = sp.symbols('om_s v_s')
    # continuous-kernel bookkeeping: discrete coefficient of zeta_q zeta_p in [eps^2]W = qperp(K_om+K_vT) + quad/2
    sol = sp.solve([kval**2 * om_s / 2 + 4 * (om_s + v_s) / t**2 + S0x, om_s + v_s + quad / 2], [om_s, v_s], dict=True)[0]
    om_mine = sp.simplify(sol[om_s]); v_mine = sp.simplify(sol[v_s])
    # xz ODE check on my algebraic omega:  -(ik/2)(om' + 2 om/t) + Sxz = 0
    resid = -(I * kval / 2) * (sp.diff(om_mine, t) + 2 * om_mine / t) + Sxz
    Kom_c = KA.closed_form(K['second']['K_om']).subs({k: pt['k'], q: pt['q'], p: pm, mu: pt['mu']}).subs(eta, t)
    KvT_c = KA.closed_form(K['second']['K_vT']).subs({k: pt['k'], q: pt['q'], p: pm, mu: pt['mu']}).subs(eta, t)
    omx = qperp * Kom_c; vx = qperp * KvT_c
    r0x = kval**2 * omx / 2 + 4 * (omx + vx) / t**2 + S0x
    vort = omx + vx + quad / 2
    N = lambda e: sp.N(e.subs(t, pt['eta']), PREC)
    out = {'my om vs C1 K_om*qperp': N(om_mine - omx), 'my vT vs C1 K_vT*qperp': N(v_mine - vx),
           'my om satisfies xz ODE': N(resid), 'C1 (om,vT) in my 0x eq': N(r0x), 'C1 (om,vT) in Kelvin eq': N(vort)}
    print("point", seed, pt)
    for nm, val in out.items():
        worst = max(worst, abs(val)); print("   %-32s %s" % (nm, sp.N(val, 12)))
print("WORST:", worst)
# ---- literature: Lu-Ananda-Clarkson-Maartens (0812.1349) eq (20), RD, Omega_m->1, w=1/3, H=1/eta, g_0i = -a^2 S_i:
#  lap S_i = 4 eta^2 { lap Phi d_i(Phi' + Phi/eta) }^V  ->  Fourier, per zeta_q zeta_p (symmetrised, P_ij q_j form)
TP = lambda kk_: sp.Function('T_Phi')(kk_, eta)
LAC = 4 * eta**2 * (-q**2) * TP(q) * (sp.diff(TP(p), eta) + TP(p) / eta) * (-I)      # lap Phi(q) x d_i Phi(p) -> (-q^2) T(q) (i p_i) (...) ; P_ij p_j = -P_ij q_j
LAC_sym = (LAC - LAC.subs({q: p, p: q}, simultaneous=True)) / 2                           # antisymmetric kernel multiplying P_ij q_j
S_LAC = sp.expand(LAC_sym / (-k**2))        # S_i = P_ij q_j S_LAC ; omega_i = -S_i (g_0i sign) if S is the eps^2-coefficient... see report
print("Lu et al (20) -> omega kernel (if omega = -S, S = full 2nd-order coefficient without 1/2):", sp.simplify(-S_LAC))
print("C1 K_om:", sp.expand(K['second']['K_om']))
print("ratio C1/(-S_LAC):", sp.simplify(sp.expand(K['second']['K_om']) / sp.expand(-S_LAC)))
