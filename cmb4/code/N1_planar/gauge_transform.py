"""Linear gauge transformation comoving-orthogonal (this code) -> Poisson gauge, single mode cos(kx).
Metric here (linear): ds^2 = t^2[-(1+2 Psi_c)dt^2 + (1+2 hA)dx^2 + (1+2 hB)(dy^2+dz^2)],
hA = zeta+u, hB = zeta+v, Psi_c = (u+2v)/3  (alh = ((1+u)(1+v)^2)^{1/3}, Ah = e^zeta(1+u), Bh = e^zeta(1+v)).
xi^mu = (T, L, 0, 0): g~ = g - Lie_xi gbar, gbar = t^2 eta.  Poisson gauge: g~_tx = 0, E = 0 (g~_xx = g~_yy).
"""
import sympy as sp
t, x, k = sp.symbols('t x k', positive=True)
hA, hB, Pc = [sp.Function(n)(t) for n in ('hA', 'hB', 'Pc')]
T, L = [sp.Function(n)(t) for n in ('T', 'L')]
c, s = sp.cos(k*x), sp.sin(k*x)
coords = [t, x]
gbar = sp.diag(-t**2, t**2, t**2, t**2)
co4 = [t, x, sp.Symbol('y'), sp.Symbol('z')]
xi = [T*c, L*s, 0, 0]
dg = sp.diag(-2*t**2*Pc*c, 2*t**2*hA*c, 2*t**2*hB*c, 2*t**2*hB*c)
Lie = sp.zeros(4, 4)
for m in range(4):
    for n in range(4):
        Lie[m, n] = sum(xi[r]*sp.diff(gbar[m, n], co4[r]) + gbar[r, n]*sp.diff(xi[r], co4[m]) + gbar[m, r]*sp.diff(xi[r], co4[n]) for r in range(4))
gt = sp.simplify(dg - Lie)
cond_tx = sp.simplify(gt[0, 1]/s)
cond_E = sp.simplify((gt[1, 1] - gt[2, 2])/c)
solL = sp.solve(cond_E, L)[0]
solT = sp.solve(cond_tx.subs(L, solL).doit(), T)[0]
print("L =", solL, "   T =", solT)
Psi = sp.simplify(-gt[0, 0].subs({L: solL, T: solT}).doit()/(2*t**2*c))
Phi = sp.simplify(-gt[2, 2].subs({L: solL, T: solT}).doit()/(2*t**2*c))
print("Psi (Poisson) =", Psi)
print("Phi (Poisson) =", Phi)
# density contrast transforms: delta_P = delta_c - (rhobar'/rhobar) T = delta_c + 4 T/t  (rhobar ~ t^-4)
print("delta_Poisson = delta_c + 4*T/t with T =", solT)
