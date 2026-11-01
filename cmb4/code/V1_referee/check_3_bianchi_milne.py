"""
V1 referee, attack item 4(a),(b).

(a) FRW radiation: shear-free, T_obs = T_o for all n (explicit redshift).
(b) Bianchi I, ds^2 = -dt^2 + sum_i a_i(t)^2 dx_i^2, with a radiation TEST fluid u = d_t
    (the claim uses only nabla.T=0 + geodesics, not Einstein's equations, so a test fluid in a
    fixed metric is a legitimate check; with a_i = t^{q_i}, sum q_i = 1 = sum q_i^2 this is also
    the vacuum Kasner solution so the fluid is a valid test fluid there).  Radiation fluid:
    rho ~ V^{-4/3}, T ~ V^{-1/3}, V = a1 a2 a3; a = 0 (Euler satisfied since T=T(t), h d_t = 0).
    Photon: p_i = const (x^i cyclic), E = sqrt(sum p_i^2/a_i^2).
    Compare EXACT T_obs = T_* E_o/E_e with the claimed T_o exp(-int sigma_ee dl),
    integrating along the actual ray with e^i(t) changing along it.
(c) Milne flow in Minkowski, u = x/tau: shear-free, T_obs isotropic = T(o). Exact algebra.
"""
import sympy as sp

t, te, to = sp.symbols('t t_e t_o', positive=True)
q1, q2, q3 = sp.symbols('q1 q2 q3', real=True)
p1, p2, p3 = sp.symbols('p1 p2 p3', real=True)
Tstar = sp.symbols('T_*', positive=True)

# ---------- (b) Bianchi I ----------
a = [t**q1, t**q2, t**q3]
V = a[0]*a[1]*a[2]
H = [sp.diff(ai, t)/ai for ai in a]
theta = sum(H)
# fluid temperature T(t) ~ V^{-1/3}, normalised so T(t_e) = T_*
Tf = Tstar*(V.subs(t, te)/V)**sp.Rational(1, 3)
# check energy conservation and Euler for this test fluid: theta = -3 u.dlnT (u = d_t), a = 0
print("Bianchi I: theta + 3 dlnT/dt == 0 :", sp.simplify(theta + 3*sp.diff(sp.log(Tf), t)) == 0)
# shear sigma_ij = a_i^2 (H_i - theta/3) delta_ij ; u = d_t has zero acceleration (geodesic in Bianchi I)
sigma = sp.diag(*[a[i]**2*(H[i] - theta/3) for i in range(3)])
pcov = [p1, p2, p3]
E = sp.sqrt(sum(pcov[i]**2/a[i]**2 for i in range(3)))
e_up = [pcov[i]/(a[i]**2*E) for i in range(3)]       # e^i = p^i/E
print("           e.e == 1 :", sp.simplify(sum(a[i]**2*e_up[i]**2 for i in range(3)) - 1) == 0)
sig_ee = sum(sigma[i, i]*e_up[i]**2 for i in range(3))
# dl = E dlambda = dt  (u = d_t):  claim says  d/dt [ln E - ln T] = -sigma_ee
lhs = sp.diff(sp.log(E), t) - sp.diff(sp.log(Tf), t)
print("           d/dt[lnE - lnT] + sigma_ee == 0 :", sp.simplify(lhs + sig_ee) == 0)
# Integrated statement with explicit numbers (Kasner exponents q=(2/3,2/3,-1/3) and a non-Kasner set)
import mpmath as mp
mp.mp.dps = 30
for qs in [(sp.Rational(2, 3), sp.Rational(2, 3), -sp.Rational(1, 3)), (sp.Rational(1, 2), sp.Rational(1, 3), sp.Rational(1, 5))]:
    for ps in [(1, 0, 0), (1, 2, 3), (0, 1, -1)]:
        sub = {q1: qs[0], q2: qs[1], q3: qs[2], p1: ps[0], p2: ps[1], p3: ps[2], te: 1, Tstar: 1}
        Eo_over_Ee = (E.subs(t, to)/E.subs(t, te)).subs(sub)
        Tobs_exact = (Tstar*Eo_over_Ee).subs(sub).subs(to, 5)
        To = Tf.subs(sub).subs(t, 5)
        S = mp.quad(sp.lambdify(t, sig_ee.subs(sub), 'mpmath'), [1, 5])
        Tobs_claim = mp.mpf(sp.N(To, 30))*mp.e**(-S)
        print(f"           q={tuple(map(str,qs))} p={ps}: T_obs exact = {mp.mpf(sp.N(Tobs_exact,25))}, claim T_o e^-S = {mp.mpf(sp.N(To,25))*mp.e**(-S)}, diff = {mp.nstr(mp.mpf(sp.N(Tobs_exact,30)) - Tobs_claim, 5)}")

# ---------- (a) FRW radiation (q1=q2=q3=q): sigma = 0, T_obs = T_o for every direction ----------
q = sp.symbols('q', positive=True)
subF = {q1: q, q2: q, q3: q}
print("FRW: sigma == 0 :", sp.simplify(sigma.subs(subF)) == sp.zeros(3, 3))
ratio = sp.simplify((Tstar*E.subs(t, to)/E.subs(t, te)).subs(subF) / Tf.subs(subF).subs(t, to))
print("FRW: T_* E_o/E_e / T_o == 1 for all p_i :", ratio == 1)

# ---------- (c) Milne flow u = x^mu/tau in Minkowski ----------
X = sp.Matrix(sp.symbols('X0:4', real=True))
eta = sp.diag(-1, 1, 1, 1)
tau = sp.sqrt(-(X.T*eta*X)[0, 0])
u_up = X/tau
u_dn = eta*u_up
Du = sp.Matrix(4, 4, lambda m, n: sp.diff(u_dn[m], X[n]))     # flat space, nabla = d
thetaM = sp.simplify(sum((eta.inv())[m, n]*Du[m, n] for m in range(4) for n in range(4)))
h = eta + u_dn*u_dn.T
hm = sp.eye(4) + u_up*u_dn.T
B = sp.Matrix(4, 4, lambda m, n: sum(hm[a, m]*hm[b, n]*Du[a, b] for a in range(4) for b in range(4)))
sigM = sp.simplify((B + B.T)/2 - thetaM/3*h)
accM = sp.simplify(sp.Matrix([sum(u_up[n]*Du[m, n] for n in range(4)) for m in range(4)]))
print("Milne: theta =", sp.simplify(thetaM*tau), "/tau ; sigma == 0 :", sigM == sp.zeros(4, 4), "; a == 0 :", accM == sp.zeros(4, 1))
# radiation test fluid: T ~ 1/tau. Photon p const, x_o = x_e + lam p: E_o/E_e = tau_e/tau_o exactly:
lam = sp.symbols('lam', positive=True)
Pn = sp.Matrix([sp.sqrt(sp.Integer(14)), 1, 2, 3])          # null
Xo = X + lam*Pn
Ee = -(u_dn.T*Pn)[0, 0]
Eo = -((eta*Xo/sp.sqrt(-(Xo.T*eta*Xo)[0, 0])).T*Pn)[0, 0]
tau_o = sp.sqrt(-(Xo.T*eta*Xo)[0, 0])
print("Milne: E_o/E_e * tau_o/tau_e == 1 :", sp.simplify(Eo/Ee*tau_o/tau - 1) == 0, " => T_obs = T_* tau_e/tau_o = T(o), isotropic.")
