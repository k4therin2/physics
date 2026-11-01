"""Linear limit check against SPEC sec.5: T_Phi(x) = -2 (sin x - x cos x)/x^3, x = k eta/sqrt3.
Gauge-invariant comoving density contrast: delta_c = -(2/3) k^2 eta^2 Psi_1 = 4 zeta (sin x - x cos x)/x.
Poisson-gauge potentials from the linear gauge transformation (gauge_transform.py):
  Phi_P = -hB - (hA' - hB')/(k^2 t),  Psi_P = Pc + [t (hA'' - hB'') + (hA' - hB')]/(k^2 t),
with hA = ln Ah, hB = ln Bh, Pc = ln alh (linear order), hA' = alh pA, hB' = alh pB, hA'' = d_t(alh pA).
"""
import numpy as np, sys
from planar import PlanarRun
import gen_numeric as gn
ic = int(sys.argv[1]) if len(sys.argv) > 1 else 4
TPhi = lambda x: -2*(np.sin(x) - x*np.cos(x))/x**3
out = []
for eps in [1e-3, 1e-4]:
    r = PlanarRun(eps, N=32, ic_order=ic, nmu=8, verbose=False, t_obar=12.0)
    r.evolve_forward()
    sol = r.sol_forward
    worst = dict(dc=0, Phi=0, Psi=0)
    rows = []
    for i in range(len(sol.t)):
        t = sol.t[i]; y = sol.y[:-1, i]
        q = r.field_quantities(t, y)
        Ah, Bh, alh, pA, pB = q[1], q[2], q[3], q[9], q[10]
        dpA = gn.f_pA_t(*q); dpB = gn.f_pB_t(*q)
        alh_t = alh**2*(pA + 2*pB)/3
        hA1 = alh*pA; hB1 = alh*pB
        hA2 = alh_t*pA + alh*dpA; hB2 = alh_t*pB + alh*dpB
        k = 1.0
        Phi = -np.log(Bh) - (hA1 - hB1)/(k**2*t)
        Psi = np.log(alh) + (t*(hA2 - hB2) + (hA1 - hB1))/(k**2*t)
        dc = alh**-4 - 1
        xx = t/np.sqrt(3)
        # project on cos(x) (first harmonic) and divide by eps
        c1 = lambda f: 2*np.mean(f*np.cos(r.x))/eps
        e_dc = c1(dc) - 4*(np.sin(xx) - xx*np.cos(xx))/xx
        e_Phi = c1(Phi) - TPhi(xx); e_Psi = c1(Psi) - TPhi(xx)
        worst['dc'] = max(worst['dc'], abs(e_dc)); worst['Phi'] = max(worst['Phi'], abs(e_Phi)); worst['Psi'] = max(worst['Psi'], abs(e_Psi))
        if i % max(1, len(sol.t)//10) == 0 or i == len(sol.t) - 1:
            rows.append((t, c1(Phi), TPhi(xx), e_Phi, e_Psi, e_dc))
    print(f"eps={eps:.0e}: max |Phi_P/eps - T_Phi| = {worst['Phi']:.2e}, max |Psi_P/eps - T_Phi| = {worst['Psi']:.2e}, max |delta_c/eps - analytic| = {worst['dc']:.2e}   (expect O(eps) from 2nd order)")
    for row in rows:
        print("   t=%7.3f  Phi_P/eps=%+.10f  T_Phi=%+.10f  dPhi=%+.1e dPsi=%+.1e ddc=%+.1e" % row)
