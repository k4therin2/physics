"""
V1 referee, attack items 1-3.

Pointwise, nonperturbative, index-explicit check of the chain
  (i)   dE/dlambda = -p^mu p^nu nabla_nu u_mu  = -E^2 (theta/3 + sigma_ee + a.e)   [omega drops]
  (ii)  energy/Euler for general (rho,P), then P=rho/3, rho=alpha T^4
  (iii) e^mu h_mu^nu = e^nu ;  theta/3 + a.e = -(1/E) dlnT/dlambda on shell
  =>    d/dlambda [ln E - ln T] = -E sigma_ee   (on shell)

Setting: a GENERIC metric g_{mu nu}(x) = G0 + G1_alpha x^alpha with random rational
coefficients (no symmetry), a generic unit timelike field u(x), a generic scalar T(x),
and a generic null p at the point x=0.  Everything is evaluated exactly (sympy rationals
and square roots).  The on-shell identity is checked as an *algebraic* identity:
    residual := dlnE/dl - dlnT/dl + sigma_ee
             == -[ (theta + 3 u.dlnT)/3 + (a + h.dlnT).e ]
i.e. residual is a specific linear combination of the energy-equation residual and the
Euler-equation residual, and vanishes iff both hold.  This is stronger than checking a
particular solution.
"""
import sympy as sp
import random

random.seed(20261101)
def rnd():
    return sp.Rational(random.randint(-9, 9), random.randint(3, 11))

x = sp.symbols('x0:4', real=True)
eta = sp.diag(-1, 1, 1, 1)

# ---- generic metric: G0 generic symmetric Lorentzian, linear x-dependence -----------
G0 = sp.Matrix(4, 4, lambda i, j: 0)
for i in range(4):
    for j in range(i, 4):
        G0[i, j] = G0[j, i] = eta[i, j] + rnd() / 4
G1 = [sp.Matrix(4, 4, lambda i, j: 0) for _ in range(4)]
for a in range(4):
    for i in range(4):
        for j in range(i, 4):
            G1[a][i, j] = G1[a][j, i] = rnd()
g = G0 + sum((G1[a] * x[a] for a in range(4)), sp.zeros(4, 4))
# check signature at the point
import numpy as np
ev = np.linalg.eigvalsh(np.array(G0.tolist(), dtype=float))
assert sum(1 for e in ev if e < 0) == 1, "G0 not Lorentzian"
ginv = g.inv()

# Christoffels Gamma^l_{mn}
def christoffel(g, ginv):
    Gam = [[[0]*4 for _ in range(4)] for _ in range(4)]
    for l in range(4):
        for m in range(4):
            for n in range(4):
                Gam[l][m][n] = sp.Rational(1, 2) * sum(
                    ginv[l, s] * (sp.diff(g[s, m], x[n]) + sp.diff(g[s, n], x[m]) - sp.diff(g[m, n], x[s]))
                    for s in range(4))
    return Gam
Gam = christoffel(g, ginv)

# ---- generic unit timelike u(x): spatial comps generic linear, u^0 from normalisation ----
w = [sp.Integer(0)] + [rnd()/3 + sum(rnd()/3 * x[a] for a in range(4)) for _ in range(3)]
# solve g_{mn} u^m u^n = -1 for u^0 (future directed root)
u0 = sp.symbols('u0')
U = sp.Matrix([u0, w[1], w[2], w[3]])
norm = (U.T * g * U)[0, 0] + 1
sols = sp.solve(norm, u0)
# pick the root that is positive at x=0
u0sol = [s for s in sols if s.subs({xx: 0 for xx in x}) > 0][0]
u_up = sp.Matrix([u0sol, w[1], w[2], w[3]])
u_dn = g * u_up

# ---- generic scalar T(x) ----
T = sp.Rational(7, 3) + sum(rnd() * x[a] for a in range(4))

# ---- generic null p at the point: spatial comps generic, p^0 from null condition ----
p0 = sp.symbols('p0')
ps = [p0, rnd(), rnd(), rnd()]
P = sp.Matrix(ps)
nullc = (P.T * G0 * P)[0, 0]
p0sol = [s for s in sp.solve(nullc, p0) if s > 0][0]
p_up = sp.Matrix([p0sol, ps[1], ps[2], ps[3]])

at0 = {xx: 0 for xx in x}
def ev0(e):
    return sp.simplify(sp.sympify(e).subs(at0))

# covariant derivative of u_mu: nabla_n u_m = d_n u_m - Gamma^l_{nm} u_l   (all at x=0)
Du = sp.zeros(4, 4)  # Du[m, n] = nabla_n u_m
for m in range(4):
    for n in range(4):
        Du[m, n] = ev0(sp.diff(u_dn[m], x[n]) - sum(Gam[l][n][m] * u_dn[l] for l in range(4)))

g0 = G0; gi0 = ev0(ginv)
u_up0 = ev0(u_up); u_dn0 = ev0(u_dn); p0v = p_up

# ---- (i): dE/dlambda from the geodesic equation, with E = -u.p --------------------
# dE/dlambda = p^n d_n(-u_m p^m) = -p^n p^m d_n u_m - u_m dp^m/dlambda,
# dp^m/dlambda = -Gamma^m_{ab} p^a p^b
dp = sp.Matrix([-sum(ev0(Gam[m][a][b]) * p0v[a] * p0v[b] for a in range(4) for b in range(4)) for m in range(4)])
dE_direct = ev0(-sum(p0v[n] * p0v[m] * sp.diff(u_dn[m], x[n]) for m in range(4) for n in range(4))
                - sum(u_dn0[m] * dp[m] for m in range(4)))
dE_claim = -sum(p0v[m] * p0v[n] * Du[m, n] for m in range(4) for n in range(4))
print("(i) dE/dlambda + p^m p^n nabla_n u_m == 0 :", sp.simplify(dE_direct - dE_claim) == 0)

# ---- kinematic decomposition at the point ----
E = -(u_dn0.T * p0v)[0, 0]
h_dn = g0 + u_dn0 * u_dn0.T              # h_{mn}
h_mixed = sp.eye(4) + u_up0 * u_dn0.T    # h^m_n  (row m, col n)
theta = sum(gi0[m, n] * Du[m, n] for m in range(4) for n in range(4))
acc = sp.Matrix([sum(u_up0[n] * Du[m, n] for n in range(4)) for m in range(4)])   # a_m = u^n nabla_n u_m
# projected gradient B_{mn} = h_m^a h_n^b nabla_b u_a   (note Du[a,b] = nabla_b u_a)
B = sp.zeros(4, 4)
for m in range(4):
    for n in range(4):
        B[m, n] = sum(h_mixed[a, m] * h_mixed[b, n] * Du[a, b] for a in range(4) for b in range(4))
sigma = (B + B.T) / 2 - theta / 3 * h_dn
omega = (B - B.T) / 2
# e^m = p^m/E - u^m
e_up = p0v / E - u_up0
e_dn = g0 * e_up
print("    u.u = -1, p.p = 0, u.e = 0, e.e = 1 :",
      [sp.simplify(v) for v in [(u_up0.T*g0*u_up0)[0,0], (p0v.T*g0*p0v)[0,0], (u_dn0.T*e_up)[0,0], (e_dn.T*e_up)[0,0]]])
# check decomposition reproduces Du
Du_rec = sigma + omega + theta/3*h_dn - acc * u_dn0.T    # nabla_n u_m = ... - a_m u_n
print("    decomposition nabla_n u_m = sigma+omega+theta/3 h - a_m u_n :", sp.simplify(Du - Du_rec) == sp.zeros(4,4))
print("    sigma traceless, sigma.u=0, omega.u=0, a.u=0 :",
      [sp.simplify(v) == 0 for v in [sum(gi0[m,n]*sigma[m,n] for m in range(4) for n in range(4)),
                                    *(sigma*u_up0), *(omega*u_up0), (acc.T*u_up0)[0,0]]])

sig_ee = (e_up.T * sigma * e_up)[0, 0]
om_ee = (e_up.T * omega * e_up)[0, 0]
a_e = (acc.T * e_up)[0, 0]
print("    omega_{mn} e^m e^n == 0 :", sp.simplify(om_ee) == 0)
print("    omega_{mn} p^m p^n == 0 :", sp.simplify((p0v.T*omega*p0v)[0,0]) == 0)
lhs = dE_claim
rhs = -E**2 * (theta/3 + sig_ee + a_e)
print("(i) dE/dlambda == -E^2 (theta/3 + sigma_ee + a.e) :", sp.simplify(lhs - rhs) == 0)

# ---- (iii): e^m h_m^n = e^n ----
e_proj = sp.Matrix([sum(e_up[m] * h_mixed[n, m] for m in range(4)) for n in range(4)])  # e^m h^n_m
print("(iii) e^m h_m^n == e^n :", sp.simplify(e_proj - e_up) == sp.zeros(4, 1))

# ---- (iii): off-shell identity. dlnT/dlambda = p^m d_m lnT ----
dlnT = sp.Matrix([ev0(sp.diff(T, x[m])) / ev0(T) for m in range(4)])      # d_m ln T (covector)
dlnT_dlam = (p0v.T * dlnT)[0, 0]
# fluid-equation residuals (zero on shell):
R_energy = theta + 3 * (u_up0.T * dlnT)[0, 0]                 # theta - (-3 u.dlnT)
hdlnT_up = sp.Matrix([sum(gi0[m, a] * sum(h_mixed[b, a] * dlnT[b] for b in range(4)) for a in range(4)) for m in range(4)])
R_euler = acc + g0 * hdlnT_up                                  # a_m + h_m^n d_n lnT
dlnE_dl = dE_claim / E**2
dlnT_dl = dlnT_dlam / E
residual = dlnE_dl - dlnT_dl + sig_ee
predicted = -(R_energy / 3 + (R_euler.T * e_up)[0, 0])
print("(iii) [dlnE/dl - dlnT/dl + sigma_ee] + [(theta+3u.dlnT)/3 + (a + h.dlnT).e] == 0 :",
      sp.simplify(residual - predicted) == 0)
print("      i.e. d/dl [ln E - ln T] = -sigma_ee exactly whenever the energy and Euler identities hold.")
print("      (sanity: residual is generically nonzero off shell:", sp.N(residual, 6), ")")
