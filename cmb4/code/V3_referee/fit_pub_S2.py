"""fit_pub_S2.py -- V3: PUB (A.6) trace source S2 vs ours.  D := S2_PUB - (-S_tr/6) (discrete coefficients) fitted on
a basis of symmetric quadratic monomials (Phi1 = Psi1 on shell).  A sparse rational fit = transcription/typo in (A.6)."""
import sympy as sp, random
from pw_tools import *
E = E_mine()
zq, zp = sp.symbols('zeta_q zeta_p')
H = 1 / t
def grad(f_): return [sp.diff(f_, v_) for v_ in (x, y, z)]
def lap(f_): return sum(sp.diff(f_, v_, 2) for v_ in (x, y, z))
def dot(a, b): return sum(a[i] * b[i] for i in range(3))
rows = []; rhs = []
names = ['Phi^2/t^2', 'Phi Phi\'/t', 'Phi\'^2', 'Phi lap Phi', 'dPhi.dPhi', 't dPhi\'.dPhi', 't^2 dPhi\'.dPhi\'', 't Phi\' lapPhi']
for seed in range(12):
    pt = rand_point(1000 + seed); km, qm = pt['k'], pt['q']
    kv, qv, pv = vec3(km, qm, pt['mu']); pm = sp.sqrt(km**2 + qm**2 - 2 * km * qm * pt['mu'])
    f = {}
    for nm, v_ in mode(zq, qv, qm).items(): f[nm] = f.get(nm, 0) + v_
    for nm, v_ in mode(zp, pv, pm).items(): f[nm] = f.get(nm, 0) + v_
    PhiP = f['Psi1']; PsiP = f['Phi1']
    V = zq * T_V(qm) * sp.exp(I * dot(qv, (x, y, z))) + zp * T_V(pm) * sp.exp(I * dot(pv, (x, y, z)))
    dV = grad(V); dPsi = grad(PsiP); dPhi = grad(PhiP); Psip = sp.diff(PsiP, t); Phip = sp.diff(PhiP, t)
    S2 = (4 * H**2 * PhiP**2 + sp.Rational(8, 3) * PsiP * lap(PsiP) + 8 * H * (PhiP - PsiP) * Psip + 8 * (-1 / t**2) * PhiP**2 - Psip**2 + 2 * Phip * Psip
          + sp.Rational(2, 3) * dot(dPhi, dPhi) + sp.Rational(2, 3) * dot(dPhi, dPsi) + dot(dPsi, dPsi) + (3 / t**2) / 3 * sp.Rational(4, 3) * dot(dV, dV)
          + sp.Rational(4, 3) * (PhiP - PsiP) * lap(PhiP) + 8 * H * (PhiP - PsiP) * Psip + 4 * (PhiP - PsiP) * sp.diff(PsiP, t, 2))
    co = lambda e_: coeff_at_origin(sp.expand(e_), (zq, zp))
    src = lambda a, b: coeff_at_origin(substitute(E[('E', a, b, 2)], f), (zq, zp))
    Str = src(1, 1) + src(2, 2) + src(3, 3)
    D = co(S2) - (-Str / 6)
    P = PsiP; Pp = Psip; dP = dPsi; dPp = grad(Pp)
    basis = [P**2 / t**2, P * Pp / t, Pp**2, P * lap(P), dot(dP, dP), t * dot(dPp, dP), t**2 * dot(dPp, dPp), t * Pp * lap(P)]
    N = lambda e_: sp.N(e_.subs(t, pt['eta']), 40)
    rows.append([N(co(b)) for b in basis]); rhs.append(N(D))
import mpmath as mp
mp.mp.dps = 40
A = mp.matrix([[mp.mpf(sp.re(c)) for c in r] for r in rows]); b = mp.matrix([mp.mpf(sp.re(c)) for c in rhs])
sol = mp.lu_solve(A[:8, :], b[:8]) if False else None
# least squares
At = A.T; sol = mp.lu_solve(At * A, At * b)
res = A * sol - b
print("fit coefficients of D = S2_PUB + S_tr/6:")
for nm, c in zip(names, sol): print("   %-20s %s" % (nm, mp.nstr(c, 20)))
print("max residual:", mp.nstr(max(abs(r_) for r_ in res), 5), " (vs |D| ~", mp.nstr(max(abs(r_) for r_ in b), 5), ")")
