"""Check 4(b): dilation / constant-zeta mode (SPEC §9d).
Exact FRW (a = eta) in rescaled coordinates eta = lam*etat, x = mu*xt:
   ds^2 = lam^2 etat^2 [ -lam^2 detat^2 + mu^2 dxt^2 ] = a(etat)^2 ghat,  a(etat) = etat
   => ghat_00 = -lam^4, ghat_ij = lam^2 mu^2 delta_ij, omega = h = 0, v = 0, 1 + delta_T = 1/lam, uhat^0 = lam^-2.
   zeta_nl = ln(1+dT) + (1/2) ln(1-2Phi) = ln mu  (exact).
lam = 1 + eps a1 + eps^2 a2/2,  mu = exp(eps b1 + eps^2 b2/2).  First order fixed by Psi1 = Phi1 = -2c/3, dT1 = c/3.
"""
import sympy as sp, pickle, os
from ser import Ser, N, eps
import fields
from fields import build, X, eta, x, y, z, F, fname

c, a1, a2, b1, b2 = sp.symbols('c alpha1 alpha2 beta1 beta2')
lam = Ser([1, a1, a2/2, 0]); mu = Ser([0, b1, b2/2, 0]).exp()
Psi = (lam*lam*lam*lam - 1)*sp.Rational(1, 2)
Phi = (1 - lam*lam*mu*mu)*sp.Rational(1, 2)
dTc = lam.inv() - 1
zeta = (1 + dTc).log() + (1 - 2*Phi).log()*sp.Rational(1, 2)
print('Psi  =', Psi.c); print('Phi  =', Phi.c); print('dT   =', dTc.c); print('zeta_nl =', [sp.simplify(t) for t in zeta.c])
sol1 = sp.solve([Psi.c[1] + 2*c/3, Phi.c[1] + 2*c/3], [a1, b1], dict=True)[0]
print('first order: ', sol1, ' -> dT1 =', dTc.c[1].subs(sol1), '(SPEC: c/3)')
# zeta_nl^(2) = 0  => b2 = 0 ;  a2 stays FREE (residual constant time rescaling)
zeta2 = sp.expand(2*zeta.c[2].subs(sol1)); print('zeta_nl^(2) =', zeta2, ' => beta2 = 0')
sol = dict(sol1); sol[b2] = 0
consts = {'Psi1': Psi.c[1], 'Psi2': 2*Psi.c[2], 'Phi1': Phi.c[1], 'Phi2': 2*Phi.c[2], 'dT1': dTc.c[1], 'dT2': 2*dTc.c[2]}
consts = {k_: sp.expand(v_.subs(sol)) for k_, v_ in consts.items()}
print('second-order constants (alpha2 free):', consts)

# ---- (i) exact formula: sigma == 0 identically for constant fields -> S == 0, Lambda == 0 ----
d = build(verbose=False)
sub = {}
for nm, val in consts.items(): sub[F(nm)] = val
for nmf in ['Psi3', 'Phi3', 'dT3']: sub[F(nmf)] = 0
for i in (1, 2, 3):
    for n_ in (1, 2, 3): sub[fields.fname('v', n_, i)] = 0
    for n_ in (2, 3):
        sub[fields.fname('om', n_, i)] = 0
        for j in (1, 2, 3):
            if i <= j: sub[sp.Function('h_%d%d_%d' % (i, j, n_))(*X)] = 0
def const_sub(e): return sp.expand(sp.expand(e).subs(sub).doit())
sig = d['fl']['sig']
print('sigma_mu nu == 0 for constant fields through order 3:', all(const_sub(sig[m][q].c[o]) == 0 for m in range(4) for q in range(4) for o in range(N+1)))
print('uhat^0 for constant fields:', [const_sub(t) for t in d['fl']['u0'].c], ' (expect lam^-2 =', (lam*lam).inv().subs(sol).c, ')')

# ---- (ii) observer clock: fixed proper time from the crossing of T_phys = T_i ----
eta_i, eta_ob = sp.symbols('eta_i eta_o_bar', positive=True)
etat = sp.Symbol('etat')
u0c = (lam*lam).inv()                       # uhat^0 = lam^-2 (constant)
# d tau = a/uhat^0 detat = etat lam^2 detat ; eta_init from 1 + dT(eta_init) = eta_init/eta_i  -> eta_init = eta_i/lam
eta_init = lam.inv()*eta_i
# Delta tau = lam^2 (eta_o^2 - eta_init^2)/2  must equal  (eta_ob^2 - eta_i^2)/2  -> solve for eta_o as series
e_o = sp.Symbol('e_o'); lam2 = lam*lam
# unknown eta_o = eta_ob + eps d1 + eps^2 d2/2
d1, d2 = sp.symbols('d1 d2')
eta_o_ser = Ser([eta_ob, d1, d2/2, 0])
dtau = lam2*(eta_o_ser*eta_o_ser - eta_init*eta_init)*sp.Rational(1, 2) - Ser.const((eta_ob**2 - eta_i**2)/2)
s1 = sp.solve(dtau.c[1], d1)[0]; s2 = sp.solve(dtau.c[2].subs(d1, s1), d2)[0]
eta_o_sol = Ser([eta_ob, s1, s2/2, 0])
print('eta_o = ', [sp.simplify(t) for t in eta_o_sol.c], ' ; eta_ob/lam =', [sp.simplify(t) for t in (lam.inv()*eta_ob).c])
# local temperature at the observer relative to the background at etabar_o:
# T_o/Tbar(etabar_o) = (1+dT)(eta_ob/eta_o)
ratio = (1 + dTc)*(eta_o_sol.inv()*eta_ob)
print('T_o/Tbar(etabar_o) =', [sp.simplify(t) for t in ratio.c], '  (monopole unperturbed through 2nd order)')
# emission: 1+dT = eta_e/eta_*  -> eta_e = eta_*/lam ; a_e/a_o = eta_e/eta_o = eta_*/eta_ob exactly
# T_obs = T_* (Ehat_o/Ehat_e) (a_e/a_o), Ehat const (Gamma = 0, k const, u const) -> T_obs = T_* eta_*/eta_ob = Tbar(eta_ob): Theta == 0.

# ---- (iii) the assembled term lists with constant fields (needs L_terms.pkl) ----
if os.path.exists('L_terms.pkl'):
    with open('L_terms.pkl', 'rb') as f: LT = pickle.load(f)
    syms = LT['symbols']; n = syms['n']; eta_s = syms['eta_star']
    em = LT['emission']; de1, de2 = syms['deta_e']
    # emission shifts for constant fields: expect eta_*(1/lam - 1) expanded
    de1v = const_sub(em[de1]); de2v = const_sub(em[de2].subs(de1, em[de1]))
    print('deta_e1 =', de1v, ' deta_e2 =', de2v, ' ; eta_*(1/lam-1) series:', [sp.expand(t) for t in (eta_s*(lam.inv() - 1)).subs(sol).c])
    # all k_n, Dx_n functionals: substitute explicit constant-field values. k1 = k1_init (G1=0): radial -n(Phi1+Psi1), etc.
    import ray
    kf, Dxf = syms['kfun'], syms['Dxfun']
    kdefs = {}
    for i in range(3):
        kdefs[kf[1][i]] = const_sub(ray.k1_init[i]); kdefs[kf[1][i].subs(eta, eta_s)] = kdefs[kf[1][i]]
    for form in ('L', 'L_direct'):
        for m in (1, 2, 3):
            tot = 0
            for t_ in LT[form][m]:
                e = const_sub(t_['expr']).xreplace(kdefs)
                e = sp.expand(e.subs(de1, de1v).subs(de2, de2v))
                e = ray.unit_rules(e)
                tot += e
            print('%s order %d with constant (dilation) fields: %s' % (form, m, 'ZERO' if sp.expand(tot) == 0 else sp.expand(tot)))
else:
    print('L_terms.pkl not present yet; skipping (iii)')
