"""check_pub.py -- V3: compare MY second-order sources (discrete plane waves) with Pitrou-Uzan-Bernardeau (A.1)-(A.8)
(literature_extracts.md 2.2; PUB lapse Phi = our Psi, PUB curvature Psi = our Phi; X = X1 + X2/2; single fluid w=cs^2=1/3,
kappa rhobar a^2 = 3/eta^2, v_I = d_I V with V_PUB = our V1 (calibrated by (A.4) at first order)).
Expected (continuous kernels, see report_V3): discrete coeff of zeta_q zeta_p in  S1 = -K_S00 ; S3 = 2 K_Q ; S4 = -K_S0z/(ik);
S2: calibrated from the linear operators."""
import sympy as sp
from pw_tools import *
E = E_mine()
zq, zp = sp.symbols('zeta_q zeta_p')
PREC = 25
H = 1 / t; Hp = -1 / t**2
kap_rho = 3 / t**2            # kappa rhobar (a^2 absorbed)
w = sp.Rational(1, 3)
def grad(f_):
    return [sp.diff(f_, v_) for v_ in (x, y, z)]
def lap(f_):
    return sum(sp.diff(f_, v_, 2) for v_ in (x, y, z))
def dot(a, b):
    return sum(a[i] * b[i] for i in range(3))
for seed in range(3):
    pt = rand_point(seed)
    kv, qv, pv = vec3(pt['k'], pt['q'], pt['mu']); km = pt['k']; qm = pt['q']
    pm = sp.sqrt(km**2 + qm**2 - 2 * km * qm * pt['mu'])
    f = {}
    for nm, v_ in mode(zq, qv, qm).items(): f[nm] = f.get(nm, 0) + v_
    for nm, v_ in mode(zp, pv, pm).items(): f[nm] = f.get(nm, 0) + v_
    # our fields -> PUB names:  Phi_PUB (lapse) = our Psi1, Psi_PUB = our Phi1, V = our V1 (v_i = d_i V), delta = our delta1
    PhiP = f['Psi1']; PsiP = f['Phi1']
    V = zq * T_V(qm) * sp.exp(I * dot(qv, (x, y, z))) + zp * T_V(pm) * sp.exp(I * dot(pv, (x, y, z)))
    de = f['delta1']
    dV = grad(V); dPsi = grad(PsiP); dPhi = grad(PhiP)
    Psip = sp.diff(PsiP, t); Phip = sp.diff(PhiP, t)
    S1 = -8 * PsiP * lap(PsiP) - 3 * dot(dPsi, dPsi) - 3 * Psip**2 + kap_rho * (1 + w) * dot(dV, dV) - 12 * H * (PhiP - PsiP) * Psip - 12 * H**2 * PhiP**2
    S2 = (4 * H**2 * PhiP**2 + sp.Rational(8, 3) * PsiP * lap(PsiP) + 8 * H * (PhiP - PsiP) * Psip + 8 * Hp * PhiP**2 - Psip**2 + 2 * Phip * Psip
          + sp.Rational(2, 3) * dot(dPhi, dPhi) + sp.Rational(2, 3) * dot(dPhi, dPsi) + dot(dPsi, dPsi) + kap_rho / 3 * (1 + w) * dot(dV, dV)
          + sp.Rational(4, 3) * (PhiP - PsiP) * lap(PhiP) + 8 * H * (PhiP - PsiP) * Psip + 4 * (PhiP - PsiP) * sp.diff(PsiP, t, 2))
    # S3 with Delta^-1 -> -1/k^2 and (Delta Delta)^-1 d_I d_J [X_IJ] -> (k_I k_J / k^4) X_IJ acting on the e^{ik.x} mode (total momentum k)
    inv1 = -1 / km**2; invIJ = lambda Xij: sum(kv[i] * kv[j] * Xij(i, j) for i in range(3) for j in range(3)) / km**4
    dSum = [dPsi[i] + dPhi[i] for i in range(3)]
    ddPhi = lambda i, j: sp.diff(PhiP, (x, y, z)[i], (x, y, z)[j])
    S3 = (-4 * PsiP**2 - 2 * PhiP**2 - inv1 * (dot(dSum, dSum) + 2 * PsiP * lap(PhiP) + kap_rho * (1 + w) * dot(dV, dV))
          + 3 * invIJ(lambda i, j: dSum[i] * dSum[j] + 2 * PsiP * ddPhi(i, j) + kap_rho * (1 + w) * dV[i] * dV[j]))
    # S4: d_I^{-1} X_I -> -i k_I X_I / k^2 (longitudinal potential), c_s^2 = 1/3
    invI = lambda Xi: -I * sum(kv[i] * Xi(i) for i in range(3)) / km**2
    S4 = (2 * H * PhiP**2 - 4 * PsiP * Psip + 2 * invI(lambda i: Psip * dPsi[i])
          + kap_rho * invI(lambda i: (1 + w) * PsiP * dV[i] - (1 + sp.Rational(1, 3)) * de * dV[i]))
    co = lambda e_: coeff_at_origin(sp.expand(e_), (zq, zp))
    S1c, S2c, S3c, S4c = [co(e_) for e_ in (S1, S2, S3, S4)]
    # my discrete sources (= continuous kernels of E^(2))
    src = lambda a, b: coeff_at_origin(substitute(E[('E', a, b, 2)], f), (zq, zp))
    S00 = src(0, 0); S0z = src(0, 3); Sxx, Syy, Szz = src(1, 1), src(2, 2), src(3, 3)
    STL = Szz - (Sxx + Syy) / 2; Str = Sxx + Syy + Szz
    N = lambda e_: sp.N(e_.subs(t, pt['eta']), PREC)
    print("point", seed, pt)
    print("   S1 / (-K_S00)         =", sp.N(N(S1c) / N(-S00), 15))
    print("   S3 / (2 K_Q)          =", sp.N(N(S3c) / N(2 * STL / km**2), 15))
    print("   S4 / (-K_S0z/(ik))    =", sp.N(N(S4c) / N(-S0z / (I * km)), 15))
    # S2: PUB (A.2) is the trace-type equation. Our trace: E_ii^(2) = L_tr + S_tr with L_tr = 2[eps^2] linear (see check_second: L_tr).
    # Calibrate with the linear operators: L_tr(ours) vs PUB LHS2 = Psi'' + H^2 Phi + (1/3)lap(Phi - Psi) + H Phi' + 2H Psi' + 2H' Phi - (kappa rho/6) delta
    # (PUB names) -> in our names: Phi2'' + Psi2/t^2 + (1/3)lap(Psi2 - Phi2) + Psi2'/t + 2 Phi2'/t - 2 Psi2/t^2 - delta2/(2 t^2)
    kk = sp.Symbol('k', positive=True); ph = sp.exp(I * kk * z)
    A = lambda nm: sp.Function('A_' + nm)(t)
    amps = {nm: A(nm) * ph for nm in ['Phi2', 'Psi2', 'delta2']}
    Ltr = sp.expand(sum((2 * substitute(E[('E', i, i, 2)], amps)).subs({x: 0, y: 0, z: 0}) for i in (1, 2, 3)))
    LPUB = (A('Phi2').diff(t, 2) + A('Psi2') / t**2 - kk**2 / 3 * (A('Psi2') - A('Phi2')) + A('Psi2').diff(t) / t + 2 * A('Phi2').diff(t) / t - 2 * A('Psi2') / t**2 - A('delta2') / (2 * t**2))
    ratio = sp.simplify(Ltr / LPUB)
    print("   L_tr(ours) / LPUB2 =", ratio, "  (constant => (A.2) is our trace eq;  S_tr = -ratio * S2  => S2 / (-S_tr/ratio) below)")
    if ratio.is_number:
        print("   S2 / (-S_tr/ratio)    =", sp.N(N(S2c) / N(-Str / ratio), 15))
    else:
        print("   S2: no constant ratio; skipped")
