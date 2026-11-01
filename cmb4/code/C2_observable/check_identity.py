"""Milestone 1: independent verification of the exact formula (SPEC §2, §3).
(A) kinematic identity for generic p (exact, any p), hence dE/dlambda = -pp.nabla u = -(...):
      p^a p^b nabla_a u_b = sigma_pp + theta (p.p + E^2)/3 + E a.p ,  E = -u.p
(B) radiation identities (mod fluid conservation):  theta = -3 u.d ln T ,  a_mu = -h_mu^nu d_nu ln T
(C) eta-parametrised line-of-sight identity, k^mu = (1,k^i), k^i = -n^i + eps k1 + eps^2 k2/2 + eps^3 k3/6
    (abstract functions of eta), evolving by the (non-affine) geodesic eq dk^i/deta = -Gam^i_kk + Gam^0_kk k^i:
      -Gam^0_kk + d/deta ln(-u.k) - d/deta ln That + sigma_kk/(-u.k) + theta g_kk /(3(-u.k)) = 0  (mod conservation)
    This is exactly "(i) direct geodesic integration of ln(Ehat_o/Ehat_e)" minus "(ii) ln(That_o/That_e) - S"
    at the level of the integrand along the (perturbed) ray; the g_kk term vanishes on the null cone.
    No null condition and no |n|=1 is needed for (C), so n_i and k_n are completely free.
(D) S_1 reproduces the standard linear formula (SPEC §9a) via a total derivative along the background ray.
"""
import sympy as sp, time, pickle
from ser import Ser, N, eps
from fields import build, X, eta, x, y, z, F, fname, dT, v, Psi, Phi

d = build(verbose=False)
g, ginv, Gam, fl, rules = d['g'], d['ginv'], d['Gam'], d['fl'], d['rules']
uu, ul, Du, theta, sig, acc, hmix = fl['uu'], fl['ul'], fl['Du'], fl['theta'], fl['sig'], fl['acc'], fl['hmix']
lnT = (1 + dT).log()
results = {}

def sub_rules(e):
    return sp.expand(sp.expand(e).subs(rules))

def check(name, ser, use_rules=True, orders=range(1, N+1)):
    ok = True
    for n in orders:
        r = sp.expand(ser.c[n])
        if use_rules: r = sub_rules(r)
        ok &= (r == 0)
        print('   %s order %d: %s (%d terms before rules)' % (name, n, 'ZERO' if r == 0 else 'NONZERO', len(sp.Add.make_args(sp.expand(ser.c[n])))))
        if r != 0: print('      residual:', r)
    results[name] = ok
    return ok

t0 = time.time()
# ---------- (A) kinematic identity with a generic (non-null) p^mu, exact in p ----------
p = [Ser.const(sp.Symbol('p%d' % m)) for m in range(4)]
E = -sum((ul[m]*p[m] for m in range(4)), Ser.const(0))
pp = sum((g[m][n]*p[m]*p[n] for m in range(4) for n in range(4) if not g[m][n].iszero()), Ser.const(0))
lhs = sum((p[a]*p[b]*Du[a][b] for a in range(4) for b in range(4)), Ser.const(0))
rhs = sum((sig[m][n]*p[m]*p[n] for m in range(4) for n in range(4)), Ser.const(0)) \
      + theta*(pp + E*E)*sp.Rational(1, 3) + E*sum((acc[m]*p[m] for m in range(4)), Ser.const(0))
print('(A) kinematic identity pp.nabla u = sigma_pp + theta(p.p+E^2)/3 + E a.p   [exact, no rules]')
check('A_kinematic', lhs - rhs, use_rules=False)
print('   t=%.1fs' % (time.time()-t0))

# ---------- (B) radiation identities (need conservation) ----------
udlnT = sum((uu[m]*lnT.diff(X[m]) for m in range(4)), Ser.const(0))
print('(B1) theta + 3 u.d ln That = 0  [mod conservation]')
check('B_theta', theta + 3*udlnT)
print('(B2) a_mu + h_mu^nu d_nu ln That = 0  [mod conservation]')
okB2 = True
for m in range(4):
    e = acc[m] + sum((hmix[m][n]*lnT.diff(X[n]) for n in range(4)), Ser.const(0))
    okB2 &= check('B_acc_%d' % m, e)
results['B_acc'] = okB2
print('   t=%.1fs' % (time.time()-t0))

# ---------- (C) eta-parametrised identity with abstract perturbed k ----------
n = sp.symbols('n1:4')
kfun = {1: [sp.Function('k1_%d' % i)(eta) for i in (1, 2, 3)],
        2: [sp.Function('k2_%d' % i)(eta) for i in (1, 2, 3)],
        3: [sp.Function('k3_%d' % i)(eta) for i in (1, 2, 3)]}
kser = [Ser.const(1)] + [Ser([-n[i], kfun[1][i], kfun[2][i]/2, kfun[3][i]/6]) for i in range(3)]
def contract2(Tm):  # T_{mu nu} k^mu k^nu
    return sum((Tm[m][nn]*kser[m]*kser[nn] for m in range(4) for nn in range(4) if not Tm[m][nn].iszero()), Ser.const(0))
Gam0kk = contract2(Gam[0])
Gamikk = [contract2(Gam[i]) for i in (1, 2, 3)]
Gvec = [-Gamikk[i] + Gam0kk*kser[i+1] for i in range(3)]       # dk^i/deta
muk = -sum((ul[m]*kser[m] for m in range(4)), Ser.const(0))      # -u.k
gkk = contract2(g)
sigkk = contract2(sig)
def ddeta(S):
    """total eta-derivative along the ray of a Ser: d_eta + k^i d_i on field arguments, plus dk/deta -> G."""
    r = S.diff(eta)  # differentiates both explicit eta in fields and k_n(eta)
    for i in range(3):
        r = r + kser[i+1]*S.diff(X[i+1])
    # replace Derivative(k_n_i(eta), eta) by the geodesic RHS, order by order
    rep = {}
    for i in range(3):
        rep[sp.Derivative(kfun[1][i], eta)] = Gvec[i].c[1]
        rep[sp.Derivative(kfun[2][i], eta)] = 2*Gvec[i].c[2]
        rep[sp.Derivative(kfun[3][i], eta)] = 6*Gvec[i].c[3]
    return Ser([sp.expand(t.subs(rep)) for t in r.c])
lnmuk = muk.log()
ident = -Gam0kk + ddeta(lnmuk) - ddeta(lnT) + sigkk/muk + theta*gkk/(3*muk)
print('(C) LOS identity: -Gam^0_kk + d/deta[ln(-u.k) - ln That] + sigma_kk/(-u.k) + theta g_kk/(3(-u.k)) = 0  [mod conservation]')
check('C_los', ident)
print('   t=%.1fs' % (time.time()-t0))

# ---------- (D) first-order S_1 -> standard formula ----------
# I1 = sigma1_kk with background k=(1,-n). Claim: -I1 = -D[dT1 + Psi1 - n.v1] + Phi1' + Psi1' mod first-order rules,
# D = d_eta - n.d (total derivative along the background ray), so that
# L1 = -S1 = [dT+Psi-n.v]_e - [dT+Psi-n.v]_o + int (Phi'+Psi') deta.
I1 = sigkk.c[1].subs({kfun[1][i]: 0 for i in range(3)})
Q = fname('dT', 1) + fname('Psi', 1) - sum(n[i]*fname('v', 1, i+1) for i in range(3))
DQ = sp.diff(Q, eta) - sum(n[i]*sp.diff(Q, X[i+1]) for i in range(3))
resD = sub_rules(-I1 + DQ - sp.diff(fname('Phi', 1), eta) - sp.diff(fname('Psi', 1), eta))
# first-order rules only use |n|^2 = 1 ? (check: residual may be prop. to (n.n-1))
resD_unit = sp.expand(resD.subs(n[2]**2, 1 - n[0]**2 - n[1]**2))
print('(D) -I1 = -D[dT1+Psi1-n.v1] + Phi1\' + Psi1\'  [mod first-order conservation, |n|=1]:',
      'ZERO' if resD_unit == 0 else 'NONZERO: %s' % resD_unit)
results['D_linear'] = (resD_unit == 0)
print('   I1 =', I1)
print('   (residual before imposing |n|=1: %s)' % resD)

with open('check_identity_results.pkl', 'wb') as f:
    pickle.dump(dict(results=results, I1=I1), f)
print('SUMMARY:', results)
