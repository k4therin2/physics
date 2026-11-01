"""check_dilation.py -- V3: the exact dilation family eta = A eta~, x = B x~ (Psi=(A^4-1)/2, Phi=(1-A^2B^2)/2,
delta = A^-4 - 1, v = 0) must solve the independent E_{mu nu} to O(eps^3) for arbitrary A(eps), B(eps).
Also: direct nonperturbative check that it is an exact solution (constant metric rescaling of FRW)."""
import sympy as sp
from pw_tools import *
E = E_mine()
eps = sp.Symbol('eps')
a1, a2, a3, b1, b2, b3 = sp.symbols('a1 a2 a3 b1 b2 b3')
A = 1 + eps * a1 + eps**2 * a2 / 2 + eps**3 * a3 / 6
B = 1 + eps * b1 + eps**2 * b2 / 2 + eps**3 * b3 / 6
exact = {'Psi': (A**4 - 1) / 2, 'Phi': (1 - A**2 * B**2) / 2, 'delta': A**(-4) - 1}
def coeff(expr, n):
    return sp.expand(sp.factorial(n) * sp.series(expr, eps, 0, 4).removeO().coeff(eps, n))
fields = {'%s%d' % (f, n): coeff(ex, n) for f, ex in exact.items() for n in (1, 2, 3)}
tot = 0
for n in (1, 2, 3):
    for a in range(4):
        for b in range(a, 4):
            # full eps^n coefficient of E_ab: sum over contributions; since E is stored per order with fields X_n generic,
            # the eps^n part = E[(E,a,b,n)] with X_m -> fields
            ex = sp.expand(substitute(E[('E', a, b, n)], fields))
            tot += abs(ex) if ex != 0 else 0
            if ex != 0:
                print("NONZERO", a, b, n, ex)
print("dilation family solves independent E to O(eps^3):", tot == 0)
# nonperturbative check: metric g = eta^2[-A^4 deta^2 + A^2 B^2 dx^2] with rho = rhobar(eta) A^-4 : exact FRW in rescaled coords
Asym, Bsym = sp.symbols('A B', positive=True)
g = sp.diag(-t**2 * Asym**4, t**2 * Asym**2 * Bsym**2, t**2 * Asym**2 * Bsym**2, t**2 * Asym**2 * Bsym**2)
ginv = g.inv()
Xc = (t, x, y, z)
Gam = [[[sum(ginv[a, d] * (sp.diff(g[d, c], Xc[b]) + sp.diff(g[d, b], Xc[c]) - sp.diff(g[b, c], Xc[d])) for d in range(4)) / 2 for c in range(4)] for b in range(4)] for a in range(4)]
Ric = sp.zeros(4, 4)
for b in range(4):
    for d in range(4):
        Ric[b, d] = sum(sp.diff(Gam[a][b][d], Xc[a]) - sp.diff(Gam[a][b][a], Xc[d]) + sum(Gam[a][a][e] * Gam[e][b][d] - Gam[a][d][e] * Gam[e][b][a] for e in range(4)) for a in range(4))
Rs = sum(ginv[a, b] * Ric[a, b] for a in range(4) for b in range(4))
G = Ric - g * Rs / 2
u_up = [1 / (t * Asym**2), 0, 0, 0]
u_dn = [sum(g[a, b] * u_up[b] for b in range(4)) for a in range(4)]
rho8 = 3 / t**4 * Asym**(-4)
T8 = sp.Matrix(4, 4, lambda a, b: rho8 * (sp.Rational(4, 3) * u_dn[a] * u_dn[b] + g[a, b] / 3))
print("exact (nonperturbative) check G - 8piG T = 0:", sp.simplify(G - T8) == sp.zeros(4, 4))
print("zeta_nl = ln(1+delta_T) + (1/2) ln(1-2Phi) =", sp.simplify(sp.log(Asym**-1) + sp.log(Asym * Bsym)))
