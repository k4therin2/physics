"""Check 4(c): structural comparison of our second-order map with Boubekeur, Creminelli, D'Amico, Norena, Vernizzi 2009
(arXiv:0906.0980, eqs. 2.14-2.26). Their conventions -> ours: their Phi = our Psi, their Psi = our Phi, their gamma_ij = our h_ij,
their n-hat = our n (P^i = -P^0 n^i), observer and emitter at rest (v = 0), emission at FIXED coordinate time tau_e with
T_e(x_e) = Tbar(tau_e) (1 + delta_T(tau_e, x_e)), observer at the origin (dx_o = deta_o = 0), ln(1+int) -> int (MD: potentials
constant at first order).  Their eq. (2.25) at the level of ln T_obs:
   ln[T_obs/(Tbar_e a_e/a_o)] = (1/2) ln(1+2Psi)|_{e, perturbed position} - (1/2) ln(1+2Psi)|_o + ln(1+dT)|_{e, perturbed position}
                                + int_{eta_e}^{eta_o} deta (Phi' + Psi' + omega'.n - h'_nn/2)|_{background path}.
Ours (exact, P_0 form; identity T3 of check_los_machinery.py), with the same prescriptions:
   ln[T_obs/(Tbar(eta_*) eta_*/etabar_o)] = int (Psi' + Phi'|k|^2 + omega'.n - h'_kk/2)/(1+2Psi-omega.k) |_{perturbed path, perturbed k}
                                + [ln(-u.k) - ln(1+2Psi-omega.k)]_e^o + ln(1+dT)|_{e, perturbed position}.
We compute D = ours - theirs at first and second order (n = z-hat display) and then set d_eta Phi1 = d_eta Psi1 = 0.
"""
import sympy as sp, pickle, time
from ser import Ser, N, eps
from fields import build, X, eta, x, y, z, fname, dT
import los
from los import (pos_expand, Dtot, assemble, at_obs, at_emi, muk, Gam0kk, g, kser, contract2, n, deta_e1, deta_e2,
                 deta_o1, deta_o2, dxo1, dxo2, split_z, classify_split, SPLIT_PHYS, kfun, Dxfun, eta_s)
t0 = time.time()
om = [g[0][i+1] for i in range(3)]
Psi = (Ser.const(0) - g[0][0] - 1)*sp.Rational(1, 2)
den = 1 + 2*Psi - sum((om[i]*kser[i+1] for i in range(3)), Ser.const(0))
dgkk = sum((g[a][b].diff(eta)*kser[a]*kser[b] for a in range(4) for b in range(4) if not g[a][b].iszero()), Ser.const(0))
I_P0 = Ser.const(0) - dgkk/den*sp.Rational(1, 2)
B_P0 = muk.log() - den.log()
ours = assemble(I_P0, B_P0)
lnT_e = pos_expand((1 + dT).log())
# --- their expression ---
Phi = (1 - g[1][1] + g[0][0]*0)*0  # placeholder (not used)
Phi_s = Ser([0, fname('Phi', 1), fname('Phi', 2)/2, fname('Phi', 3)/6])
h_nn = sum((g[i+1][j+1]*n[i]*n[j] for i in range(3) for j in range(3)), Ser.const(0)) - (1 - 2*Phi_s)*sum(ni**2 for ni in n)
theirs_ray = (Phi_s.diff(eta) + Psi.diff(eta) + sum((om[i].diff(eta)*n[i] for i in range(3)), Ser.const(0)) - h_nn.diff(eta)*sp.Rational(1, 2))
theirs_emi = pos_expand((1 + 2*Psi).log()*sp.Rational(1, 2))
theirs_obs = Ser.const(0) - (1 + 2*Psi).log()*sp.Rational(1, 2)

vzero = {}
for i in (1, 2, 3):
    for m in (1, 2, 3): vzero[fname('v', m, i)] = 0
noobs = {deta_o1: 0, deta_o2: 0, deta_e1: 0, deta_e2: 0}
for s_ in dxo1 + dxo2: noobs[s_] = 0
def prep(e):
    e = sp.expand(sp.expand(e).xreplace(noobs).subs(vzero).doit())
    return split_z(e)
nodot = {}
def kill_dot(e):
    """set d_eta Phi1 = d_eta Psi1 = 0 (incl. mixed derivatives)"""
    rep = {}
    for D in sp.expand(e).atoms(sp.Derivative):
        if D.expr.func in (fname('Phi', 1).func, fname('Psi', 1).func) and eta in D.variables: rep[D] = 0
    return sp.expand(e.xreplace(rep))

out = {}
for m in (1, 2):
    for where in ('ray', 'obs', 'emi'):
        o_ = sum((e for w, e in ours[m] if w == where), sp.Integer(0))
        if where == 'emi': o_ += at_emi(lnT_e.c[m])
        if where == 'ray': t_ = theirs_ray.c[m]
        elif where == 'emi': t_ = at_emi(theirs_emi.c[m] + lnT_e.c[m])
        else: t_ = at_obs(theirs_obs.c[m])
        D = prep(o_ - t_)
        out[(m, where)] = D
        print('order %d [%s]: ours - theirs = %s' % (m, where, 'ZERO' if D == 0 else '%d terms' % len(sp.Add.make_args(D))))
        if D != 0:
            for key, val in classify_split(D).items(): print('      ', key, ':', val)
            print('       after setting Phi1\' = Psi1\' = 0:', kill_dot(D))
with open('lit_compare_results.pkl', 'wb') as f: pickle.dump(out, f)
print('t=%.1f' % (time.time()-t0))
