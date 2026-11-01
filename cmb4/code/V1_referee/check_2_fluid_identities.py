"""
V1 referee, attack item 2.

Derive, at a generic point of a generic metric (random rational coefficients, no symmetry),
the exact consequences of nabla_mu T^{mu nu} = 0 for T = (rho+P) u u + P g with
GENERAL rho(x), P(x):
   u_nu nabla_mu T^{mu nu}  = -[ u.d rho + (rho+P) theta ]            (energy)
   h^nu_lam nabla_mu T^{mu lam} = (rho+P) a^nu + h^{nu mu} d_mu P     (Euler)
Then ask: for a barotropic fluid rho(T), P(T), when do BOTH
   theta = -3 u.d ln T      and     a_mu = -h_mu^nu d_nu ln T
follow?  Answer derived symbolically below: iff  rho+P = T P'(T)  and  rho - 3P = const.
With const = 0 this is exactly P = rho/3, rho = alpha T^4 (then P = alpha T^4/3). No chemical
potential / particle-number assumption enters: T is *defined* by rho = alpha T^4 and the two
identities are pure consequences of P = rho/3 plus nabla.T = 0.
"""
import sympy as sp
import random
random.seed(7)
def rnd():
    return sp.Rational(random.randint(-9, 9), random.randint(3, 11))

x = sp.symbols('x0:4', real=True)
eta = sp.diag(-1, 1, 1, 1)
G0 = sp.zeros(4, 4)
for i in range(4):
    for j in range(i, 4):
        G0[i, j] = G0[j, i] = eta[i, j] + rnd()/4
G1 = [sp.zeros(4, 4) for _ in range(4)]
for a in range(4):
    for i in range(4):
        for j in range(i, 4):
            G1[a][i, j] = G1[a][j, i] = rnd()
g = G0 + sum((G1[a]*x[a] for a in range(4)), sp.zeros(4, 4))
ginv = g.inv()
Gam = [[[sp.Rational(1, 2)*sum(ginv[l, s]*(sp.diff(g[s, m], x[n]) + sp.diff(g[s, n], x[m]) - sp.diff(g[m, n], x[s]))
                                 for s in range(4)) for n in range(4)] for m in range(4)] for l in range(4)]

w = [None] + [rnd()/3 + sum(rnd()/3*x[a] for a in range(4)) for _ in range(3)]
u0 = sp.symbols('u0')
U = sp.Matrix([u0, w[1], w[2], w[3]])
sols = sp.solve((U.T*g*U)[0, 0] + 1, u0)
u0sol = [s for s in sols if s.subs({xx: 0 for xx in x}) > 0][0]
u_up = sp.Matrix([u0sol, w[1], w[2], w[3]]); u_dn = g*u_up

# generic scalar fields rho, P (linear, independent)
rho = sp.Rational(5, 2) + sum(rnd()*x[a] for a in range(4))
P = sp.Rational(3, 4) + sum(rnd()*x[a] for a in range(4))

Tup = (rho + P)*u_up*u_up.T + P*ginv          # T^{mu nu}
# nabla_mu T^{mu nu} = d_mu T^{mu nu} + Gam^mu_{mu l} T^{l nu} + Gam^nu_{mu l} T^{mu l}
divT = sp.Matrix([sum(sp.diff(Tup[m, n], x[m]) for m in range(4))
                  + sum(Gam[m][m][l]*Tup[l, n] for m in range(4) for l in range(4))
                  + sum(Gam[n][m][l]*Tup[m, l] for m in range(4) for l in range(4)) for n in range(4)])
at0 = {xx: 0 for xx in x}
ev0 = lambda e: sp.simplify(sp.sympify(e).subs(at0))
divT0 = divT.applyfunc(ev0)
g0 = G0; gi0 = ev0(ginv); u_up0 = ev0(u_up); u_dn0 = ev0(u_dn)
Du = sp.Matrix(4, 4, lambda m, n: ev0(sp.diff(u_dn[m], x[n]) - sum(Gam[l][n][m]*u_dn[l] for l in range(4))))
theta = sum(gi0[m, n]*Du[m, n] for m in range(4) for n in range(4))
acc_dn = sp.Matrix([sum(u_up0[n]*Du[m, n] for n in range(4)) for m in range(4)])
acc_up = gi0*acc_dn
h_mixed = sp.eye(4) + u_up0*u_dn0.T   # h^m_n
drho = sp.Matrix([ev0(sp.diff(rho, x[m])) for m in range(4)])
dP = sp.Matrix([ev0(sp.diff(P, x[m])) for m in range(4)])
rho0, P0 = ev0(rho), ev0(P)

energy = (u_dn0.T*divT0)[0, 0] + ((u_up0.T*drho)[0, 0] + (rho0 + P0)*theta)
print("energy:  u_nu nabla_mu T^{mu nu} + [u.drho + (rho+P) theta] == 0 :", sp.simplify(energy) == 0)
hdivT = sp.Matrix([sum(h_mixed[n, l]*divT0[l] for l in range(4)) for n in range(4)])
hdP_up = sp.Matrix([sum(gi0[n, a]*sum(h_mixed[b, a]*dP[b] for b in range(4)) for a in range(4)) for n in range(4)])
euler = hdivT - ((rho0 + P0)*acc_up + hdP_up)
print("Euler:   h^nu_l nabla_mu T^{mu l} - [(rho+P) a^nu + h^{nu mu} d_mu P] == 0 :", sp.simplify(euler) == sp.zeros(4, 1))

# ---- barotropic specialisation, purely thermodynamic algebra -------------------------
T = sp.symbols('T', positive=True)
rhoT = sp.Function('rho')(T); PT = sp.Function('P')(T)
# Euler: a = -h dP/(rho+P) = -h dlnT * [T P'/(rho+P)]   -> need T P' = rho+P
# Energy: theta = -u.drho/(rho+P) = -u.dlnT * [T rho'/(rho+P)] -> need T rho' = 3(rho+P)
print("Euler  identity a=-h dlnT   <=>  T P'(T) = rho+P")
print("Energy identity theta=-3u.dlnT <=> T rho'(T) = 3 (rho+P)")
# Solve the pair: rho + P = T P', T rho' = 3(rho+P) = 3 T P'  => rho' = 3P' => rho = 3P + c
c = sp.symbols('c')
sol = sp.dsolve(sp.Eq(T*sp.diff(PT, T), 3*PT + c + PT), PT)
print("General solution of the pair:", sol, " with rho = 3P + c")
# Check P = alpha T^4/3, rho = alpha T^4 satisfies both
alpha = sp.symbols('alpha', positive=True)
Pc = alpha*T**4/3; rhoc = alpha*T**4
print("P=rho/3, rho=alpha T^4: T P' - (rho+P) =", sp.simplify(T*sp.diff(Pc, T) - (rhoc + Pc)),
      ";  T rho' - 3(rho+P) =", sp.simplify(T*sp.diff(rhoc, T) - 3*(rhoc + Pc)))
# Does P = rho/3 ALONE suffice?  rho = alpha T^4 alone?
# P = rho/3 alone: Euler gives a = -(1/4) h dln rho, which equals -h dlnT only if rho ~ T^4. Need both.
Pgen = sp.Function('P')(T)
print("P=rho/3 alone: Euler needs T P' = 4P -> P ~ T^4; so rho ~ T^4 is forced (not independent).")
print("rho=alpha T^4 alone with P=w rho: Euler gives a = -4w/(1+w) h dlnT, energy theta = -3 u.dlnT * 4/(3(1+w));")
wsym = sp.symbols('w')
print("   both reduce to the claimed identities iff", sp.solve([4*wsym/(1 + wsym) - 1], wsym), "(w=1/3).")
