"""Principal-symbol / hyperbolicity analysis of the planar comoving-gauge ADM system.
Linearise about FRW (Ah=Bh=1, zeta=0, pA=pB=0): variables (a, b, pA, pB), Fourier e^{i kappa x}.
Second-order-in-space system: scaled variables (a, b, pA/kappa, pB/kappa); principal symbol M0 = coefficient of kappa.
Plain ADM: Jordan block at eigenvalue 0 (weakly hyperbolic).  Modified system:
   d_t pB  +=  cH * t^2 * Ham_n     (Ham_n = (G_tt - T_tt)/alpha^2, principal part -2 Bxx/(Ah^2 Bh t^2))
is strongly hyperbolic for cH < 0, cH != -1/6: eigenvalues +-i/sqrt(3) (sound), +-i sqrt(2|cH|) (constraint mode).
"""
import sympy as sp, pickle, numpy as np
d = pickle.load(open('eqs_generated.pkl', 'rb'))
S = d['symbols']
t, Ah, Bh, alh, Ax, Bx, Axx, Bxx, pA, pB, pAx, pBx, ze, zx, zxx = [S[k] for k in
    ('t', 'Ah', 'Bh', 'alh', 'Ax', 'Bx', 'Axx', 'Bxx', 'pA', 'pB', 'pAx', 'pBx', 'zeta', 'zx', 'zxx')]
eps = sp.Symbol('epsilon')
a, b, qA, qB = sp.symbols('a b qA qB')
kap = sp.Symbol('kappa', positive=True)
cH = sp.Symbol('c_H', real=True)
I = sp.I
sub = {Ah: 1 + eps*a, Bh: 1 + eps*b, Ax: eps*I*kap*a, Bx: eps*I*kap*b, Axx: -eps*kap**2*a, Bxx: -eps*kap**2*b,
       pA: eps*qA, pB: eps*qB, pAx: eps*I*kap*qA, pBx: eps*I*kap*qB, ze: 0, zx: 0, zxx: 0}
alh_e = ((1 + eps*a)*(1 + eps*b)**2)**sp.Rational(1, 3)
def lin(e):
    e = e.subs(alh, alh_e).subs(sub)
    return sp.simplify(sp.diff(e, eps).subs(eps, 0))
Y = sp.Matrix([a, b, qA, qB])
def analyse(cH_val, label):
    rhs = [lin(alh_e*(1 + eps*a)*eps*qA), lin(alh_e*(1 + eps*b)*eps*qB), lin(d['pA_t']),
           lin(d['pB_t'] + cH_val*t**2*d['Ham'])]
    M = sp.Matrix([[sp.diff(r, v) for v in Y] for r in rhs])
    Sc = sp.diag(1, 1, kap, kap)
    Ms = sp.expand(Sc.inv()*M*Sc)
    M0 = Ms.applyfunc(lambda e: sp.limit(e/kap, kap, sp.oo))
    print(f"\n=== {label} ===\nLinearised M(kappa,t):"); sp.pprint(M)
    print("Principal symbol M0 (coefficient of kappa):"); sp.pprint(M0)
    for val, mult, vecs in M0.eigenvects():
        print("  eigenvalue", val, " alg mult", mult, " geom mult", len(vecs))
    # numerical check of full M incl. lower-order terms: max real part of eigenvalues
    Mf = sp.lambdify((kap, t), M, 'numpy')
    worst = 0
    for kk in [1, 2, 4, 8, 16, 32, 64]:
        for tt in [0.01, 0.1, 1, 3, 10]:
            ev = np.linalg.eigvals(np.array(Mf(kk, tt), dtype=complex))
            worst = max(worst, np.max(np.abs(ev.real)))
    print("  max |Re eigenvalue| of full M over kappa in 1..64, t in 0.01..10:", worst)
    return M
analyse(0, "plain ADM (cH = 0)")
analyse(sp.Rational(-1, 2), "modified, cH = -1/2")
print("\nHam lin:", lin(d['Ham'])); print("Mom lin:", lin(d['Mom']))

# kappa-dependence of Re(eigenvalues) at fixed t (bounded in kappa <=> no exponential ill-posedness)
for cval, lab in [(0, 'plain'), (sp.Rational(-1, 2), 'cH=-1/2')]:
    rhs = [lin(alh_e*(1 + eps*a)*eps*qA), lin(alh_e*(1 + eps*b)*eps*qB), lin(d['pA_t']), lin(d['pB_t'] + cval*t**2*d['Ham'])]
    M = sp.Matrix([[sp.diff(r, v) for v in Y] for r in rhs])
    Mf = sp.lambdify((kap, t), M, 'numpy')
    for tt in [0.1, 1.0, 10.0]:
        row = []
        for kk in [1, 4, 16, 64, 256]:
            ev = np.linalg.eigvals(np.array(Mf(kk, tt), dtype=complex))
            row.append(np.max(ev.real))
        print(f"{lab:8s} t={tt:5.1f}  max Re(eig) for kappa=1,4,16,64,256:", np.array2string(np.array(row), precision=3))
