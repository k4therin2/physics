"""Consistency checks of the line-of-sight machinery (los.py) and the expanded direct-vs-shear comparison.
T1 (machinery): for any boundary function B,  [B]_e^o  (shifted endpoints, shift_eta)  ==  int_{eta_e}^{eta_o} dB/deta
    (ray integral of Dtot(B) + endpoint corrections of the integral), after int_{eta_*}^{etabar_o} Dtot(Bq) -> Bq(o) - Bq(e).
    Must hold identically (no equations of motion) through eps^3.
T2 (direct vs shear, EXPANDED): Lambda^direct - Lambda^shear = assemble( R ) with
    R = -Gam0kk + Dtot[pos_expand(ln(-u.k) - ln That)] + pos_expand(sigma_kk/(-u.k) + theta g_kk/(3(-u.k)))
    (the theta g_kk term vanishes for the actual null k; it is needed because k1,k2 are kept as free functions).
    Every piece (ray / obs / emi, orders 1,2,3) must vanish modulo the conservation rules AND their derivatives.
T3 (exact P_0 identity used by Boubekeur et al. 2009 eq. 2.15-2.16):
    -Gam^0_kk  ==  -(1/2) d_eta g_{ab} k^a k^b /(1 + 2Psi - omega.k)  -  Dtot ln(1 + 2Psi - omega.k)   identically (any fields)
    [ln p^0 = ln(-p_0) - ln(1+2Psi-omega.k),  d ln(-p_0)/deta = -(1/2) d_eta g_kk/(1+2Psi-omega.k)].
"""
import sympy as sp, time, pickle
from ser import Ser, N, eps
from fields import build, X, eta, x, y, z, fname
import los
from los import (I_sigma, I_direct, B_direct, pos_expand, Dtot, assemble, at_obs, at_emi, muk, Gam0kk, theta, g, kser,
                 contract2, rules, S_sigma, S_direct_neg, n)
t0 = time.time()
res = {}

# ---- T1 ----
Bq = pos_expand(B_direct)
A = assemble(I_direct, B_direct)
Bp = assemble(I_direct + Dtot(Bq))
ok = True
for m in (1, 2, 3):
    ray_A = sum(e for w, e in A[m] if w == 'ray'); ray_B = sum(e for w, e in Bp[m] if w == 'ray')
    # ray difference must be exactly -Dtot(Bq).c[m]; replace its integral by endpoint values
    dray = sp.expand(ray_A - ray_B + Dtot(Bq).c[m])
    ok &= (dray == 0)
    obs = sum(e for w, e in A[m] if w == 'obs') - sum(e for w, e in Bp[m] if w == 'obs') - at_obs(Bq.c[m])
    emi = sum(e for w, e in A[m] if w == 'emi') - sum(e for w, e in Bp[m] if w == 'emi') + at_emi(Bq.c[m])
    r_obs, r_emi = sp.expand(obs), sp.expand(emi)
    print('T1 order %d: ray-diff %s, obs residual %s, emi residual %s' % (m, 'ZERO' if dray == 0 else 'NONZERO', 'ZERO' if r_obs == 0 else 'NONZERO(%d)' % len(sp.Add.make_args(r_obs)), 'ZERO' if r_emi == 0 else 'NONZERO(%d)' % len(sp.Add.make_args(r_emi))))
    ok &= (r_obs == 0) and (r_emi == 0)
res['T1_endpoint_machinery'] = ok
print('t=%.1f' % (time.time()-t0), flush=True)

# ---- extended conservation rules (derivatives of the rules) ----
dyn_funcs = set()
for k_ in rules: dyn_funcs.add(k_.expr.func)
def reduce_rules(e, maxit=12):
    """replace every eta-derivative (possibly mixed) of dT_n, v_n^i using the conservation rules, iteratively."""
    e = sp.expand(e)
    for it in range(maxit):
        rep = {}
        for D in e.atoms(sp.Derivative):
            if D.expr.func in dyn_funcs and eta in D.variables:
                vars_ = list(D.variables); vars_.remove(eta)
                base = sp.Derivative(D.expr, eta)
                if base in rules:
                    rep[D] = sp.diff(rules[base], *vars_) if vars_ else rules[base]
        if not rep: return e
        e = sp.expand(e.xreplace(rep))
    raise RuntimeError('rules did not converge')

# ---- T2 ----
gkk = contract2(g)
R = Ser.const(0) - Gam0kk + Dtot(Bq) + pos_expand(I_sigma + theta*gkk/(3*muk))
P = assemble(R)
ok = True
for m in (1, 2, 3):
    for w in ('ray', 'obs', 'emi'):
        tot = sum((e for w_, e in P[m] if w_ == w), sp.Integer(0))
        r = reduce_rules(tot)
        print('T2 order %d [%s]: %s  (%d terms before rules)' % (m, w, 'ZERO' if r == 0 else 'NONZERO(%d)' % len(sp.Add.make_args(r)), len(sp.Add.make_args(sp.expand(tot)))), flush=True)
        ok &= (r == 0)
res['T2_direct_vs_shear_expanded'] = ok
print('t=%.1f' % (time.time()-t0), flush=True)

# ---- T3 ----
om = [g[0][i+1] for i in range(3)]
Psi = (Ser.const(0) - g[0][0] - 1)*sp.Rational(1, 2)
den = 1 + 2*Psi - sum((om[i]*kser[i+1] for i in range(3)), Ser.const(0))
dgkk = sum((g[a][b].diff(eta)*kser[a]*kser[b] for a in range(4) for b in range(4) if not g[a][b].iszero()), Ser.const(0))
T3 = Ser.const(0) - Gam0kk + dgkk/den*sp.Rational(1, 2) + Dtot(den.log())
ok = all(sp.expand(T3.c[m]) == 0 for m in range(N+1))
print('T3 (P_0 form identity, exact, no rules):', 'ZERO' if ok else 'NONZERO', [len(sp.Add.make_args(sp.expand(T3.c[m]))) for m in range(N+1)])
res['T3_P0_identity'] = ok
print('SUMMARY', res, 't=%.1f' % (time.time()-t0))
with open('check_los_machinery_results.pkl', 'wb') as f: pickle.dump(res, f)
