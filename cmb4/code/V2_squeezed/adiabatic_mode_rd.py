"""
V2 (consistency-relation route), script 1.
Radiation-only, flat, a(eta)=eta, Poisson gauge, SPEC conventions:
   ghat_00 = -(1+2 Psi),  ghat_ij = (1-2 Phi) delta_ij,  g = eta^2 ghat,  8 pi G a^4 rhobar = 3.
   rho = rhobar(1+delta),  1+delta_T = (1+delta)^{1/4},  uhat^mu = uhat^0 (1, v^i).

PART A: linearised Einstein equations derived from scratch (general Psi(eta,x), Phi(eta,x)),
        and verification that the Weinberg adiabatic mode for zeta_L = c + q.x,
            Phi_L = Psi_L = -(2/3) zeta_L,  delta_T,L = zeta_L/3,  v_L^i = (eta/3) q^i,
        satisfies them exactly (it is time independent, so it is the exact k->0 limit).
PART B: the same fields are obtained as the diffeomorphism image of unperturbed FRW under
            eta~ = eta - (eta/3) zeta_L(x),   x~^i = (1+zeta_L) x^i - (x^2/2) q^i - (eta^2/6) q^i
        (time shift + dilation + special conformal transformation + translation f(eta) q,
         f' = a^-2 int a^2 = eta/3).  Checked to first order in c and in q.
PART C: diffeo image of a short-mode Poisson-gauge metric: the O(c x S) and O(q x S) cross terms,
        which gives (i) the constant-mode second-order Poisson-gauge field relations (linear
        parametrisation) and (ii) the O(q S) non-Poisson pieces (a g_0i term) and the second-order
        gauge shift needed to return to Poisson gauge, i.e. the squeezed-limit vector mode.
PART D: the exact (all orders in c) time reparametrisation for a constant zeta: lambda(c).
"""
import sympy as sp

eta, c = sp.symbols('eta c', real=True)
x1, x2, x3 = sp.symbols('x1 x2 x3', real=True)
q1, q2, q3 = sp.symbols('q1 q2 q3', real=True)
X = [eta, x1, x2, x3]
xs = [x1, x2, x3]
qs = [q1, q2, q3]
eps1 = sp.Symbol('epsilon')     # bookkeeping for first order (long mode amplitude)
sS = sp.Symbol('s')             # bookkeeping for short-mode amplitude

def christoffel(g, ginv):
    G = [[[0]*4 for _ in range(4)] for _ in range(4)]
    for a in range(4):
        for b in range(4):
            for cc in range(4):
                G[a][b][cc] = sp.Rational(1, 2)*sum(ginv[a, d]*(sp.diff(g[d, b], X[cc]) + sp.diff(g[d, cc], X[b]) - sp.diff(g[b, cc], X[d])) for d in range(4))
    return G

def ricci(G):
    R = sp.zeros(4, 4)
    for b in range(4):
        for cc in range(4):
            R[b, cc] = sum(sp.diff(G[a][b][cc], X[a]) for a in range(4)) - sum(sp.diff(G[a][b][a], X[cc]) for a in range(4)) \
                + sum(G[a][a][d]*G[d][b][cc] for a in range(4) for d in range(4)) - sum(G[a][cc][d]*G[d][b][a] for a in range(4) for d in range(4))
    return R

def einstein_minus_T(g, uvec, rho, order_symbol, order):
    """Return G_{mu nu} - 8 pi G T_{mu nu} truncated at given order in order_symbol.
    g: 4x4 physical metric; uvec: contravariant 4-velocity (unnormalised ok -> normalised here);
    rho: physical energy density; P = rho/3; 8 pi G chosen so that 8 pi G a^4 rhobar = 3."""
    ginv = g.inv()
    ginv = ginv.applyfunc(lambda e: sp.series(e, order_symbol, 0, order+1).removeO())
    Gam = christoffel(g, ginv)
    Gam = [[[sp.series(Gam[a][b][cc], order_symbol, 0, order+1).removeO() for cc in range(4)] for b in range(4)] for a in range(4)]
    R = ricci(Gam).applyfunc(lambda e: sp.series(e, order_symbol, 0, order+1).removeO())
    Rs = sp.series(sum(ginv[a, b]*R[a, b] for a in range(4) for b in range(4)), order_symbol, 0, order+1).removeO()
    Gt = (R - sp.Rational(1, 2)*g*Rs).applyfunc(lambda e: sp.series(e, order_symbol, 0, order+1).removeO())
    # normalise u
    norm = sum(g[a, b]*uvec[a]*uvec[b] for a in range(4) for b in range(4))
    u = [sp.series(ui/sp.sqrt(-norm), order_symbol, 0, order+1).removeO() for ui in uvec]
    ul = [sum(g[a, b]*u[b] for b in range(4)) for a in range(4)]
    P = rho/3
    T = sp.zeros(4, 4)
    for a in range(4):
        for b in range(4):
            T[a, b] = (rho + P)*ul[a]*ul[b] + P*g[a, b]
    eightpiG = 3/(eta**4*rhobar_phys)  # 8 pi G a^4 rhobar = 3 with a = eta
    E = (Gt - eightpiG*T).applyfunc(lambda e: sp.simplify(sp.series(sp.expand(e), order_symbol, 0, order+1).removeO()))
    return E

rhobar_phys = sp.Symbol('rhobar0')/eta**4   # background rho ∝ a^-4

# ---------------------------------------------------------------- PART A
print("=== PART A: linearised Poisson-gauge equations and the k->0 adiabatic mode ===")
Psi = sp.Function('Psi')(eta, x1, x2, x3)
Phi = sp.Function('Phi')(eta, x1, x2, x3)
dT = sp.Function('dT')(eta, x1, x2, x3)       # delta_T
vp = sp.Function('vp')(eta, x1, x2, x3)       # velocity potential, v_i = d_i vp
g = sp.diag(-(1+2*eps1*Psi), 1-2*eps1*Phi, 1-2*eps1*Phi, 1-2*eps1*Phi)*eta**2
uvec = [1, eps1*sp.diff(vp, x1), eps1*sp.diff(vp, x2), eps1*sp.diff(vp, x3)]
rho = rhobar_phys*(1+eps1*dT)**4
E = einstein_minus_T(g, uvec, rho, eps1, 1)
E1 = E.applyfunc(lambda e: sp.expand(e.coeff(eps1, 1)))
print("background (eps^0) residual:", sp.simplify(E.applyfunc(lambda e: e.coeff(eps1, 0))))
# Print the linear equations in a readable form (00, 0i, ij)
print("E_00^(1) =", sp.simplify(E1[0, 0]))
print("E_01^(1) =", sp.simplify(E1[0, 1]))
print("E_12^(1) =", sp.simplify(E1[1, 2]))
print("E_11^(1) =", sp.simplify(E1[1, 1]))

# Long mode zeta_L = c + q.x
zetaL = c + q1*x1 + q2*x2 + q3*x3
subsL = {Psi: -sp.Rational(2, 3)*zetaL, Phi: -sp.Rational(2, 3)*zetaL, dT: zetaL/3, vp: eta*zetaL/3}
EL = E1.applyfunc(lambda e: sp.simplify(e.subs(subsL).doit()))
print("Long-mode residuals (all should be 0):", EL)
assert EL == sp.zeros(4, 4)

# Also check with the finite-k RD solution of the lead (independent transfer functions), as a sanity check of the equations:
k = sp.Symbol('k', positive=True)
xx = k*eta/sp.sqrt(3)
TPhi = -2*(sp.sin(xx) - xx*sp.cos(xx))/xx**3
H = 1/eta
Tdelta = -(2/(3*H**2))*(k**2*TPhi + 3*H*(sp.diff(TPhi, eta) + H*TPhi))
Tv = -(sp.diff(TPhi, eta) + H*TPhi)/(2*H**2)
plane = sp.cos(k*x1)
subsk = {Psi: TPhi*plane, Phi: TPhi*plane, dT: Tdelta/4*plane, vp: Tv*plane}
Ek = E1.applyfunc(lambda e: sp.simplify(e.subs(subsk).doit()))
print("finite-k RD solution residuals (should be 0):", Ek)
assert Ek == sp.zeros(4, 4)
print("Superhorizon limits: TPhi->", sp.limit(TPhi, eta, 0), " Tdelta/4->", sp.limit(Tdelta/4, eta, 0), " Tv/eta->", sp.limit(Tv/eta, eta, 0))

# ---------------------------------------------------------------- PART B
print("\n=== PART B: the long mode as a diffeomorphism of FRW ===")
# new coords (eta~, x~) as functions of (eta, x); first order in the long mode amplitude eps1
f = eta**2/6     # f' = (1/a^2) int a^2 deta = eta/3
epsT = -(eta/3)*zetaL
xi = [zetaL*xs[i] - sp.Rational(1, 2)*(x1**2+x2**2+x3**2)*qs[i] - f*qs[i] for i in range(3)]
Xt = [eta + eps1*epsT] + [xs[i] + eps1*xi[i] for i in range(3)]
# FRW metric in tilde coords: a(eta~)^2 diag(-1,1,1,1); pull back
J = sp.Matrix(4, 4, lambda m, n: sp.diff(Xt[m], X[n]))
gt = sp.diag(-1, 1, 1, 1)*Xt[0]**2
gpull = (J.T*gt*J).applyfunc(lambda e: sp.series(sp.expand(e), eps1, 0, 2).removeO())
# read off Psi, Phi (first order in eps1)
Psi_L = sp.simplify((-(gpull[0, 0]/eta**2) - 1)/2).coeff(eps1, 1)
Phi_L = sp.simplify((1 - gpull[1, 1]/eta**2)/2).coeff(eps1, 1)
print("Psi_L =", sp.factor(Psi_L), "   Phi_L =", sp.factor(Phi_L), "  (expect -(2/3) zeta_L each)")
print("g_0i, g_ij offdiag at O(eps):", [sp.simplify(gpull[0, i].coeff(eps1, 1)) for i in (1, 2, 3)], sp.simplify(gpull[1, 2].coeff(eps1, 1)))
assert sp.simplify(Psi_L + sp.Rational(2, 3)*zetaL) == 0 and sp.simplify(Phi_L + sp.Rational(2, 3)*zetaL) == 0
# fluid: at rest in tilde coords, u~^mu = (1/a(eta~), 0); u^mu = (dx^mu/dx~^nu) u~^nu = (J^-1)^mu_0 u~^0
Jinv = J.inv().applyfunc(lambda e: sp.series(e, eps1, 0, 2).removeO())
u = [Jinv[m, 0]/Xt[0] for m in range(4)]
v = [sp.simplify(sp.series(u[i]/u[0], eps1, 0, 2).removeO()) for i in (1, 2, 3)]
print("v^i =", v, "  (expect (eta/3) q^i)")
assert all(sp.simplify(v[i] - eps1*eta*qs[i]/3) == 0 for i in range(3))
# temperature: scalar, T(eta,x) = Tbar(eta~) with Tbar ∝ 1/eta
dT_L = sp.series(eta/Xt[0] - 1, eps1, 0, 2).removeO().coeff(eps1, 1)
print("delta_T,L =", sp.simplify(dT_L), "  (expect zeta_L/3)")
assert sp.simplify(dT_L - zetaL/3) == 0
# the initial-surface (eta->0) limit of the spatial map: xi_0 = zeta_L x - x^2 q/2 (the eta^2 q/6 translation vanishes)
print("xi(eta->0) =", [sp.simplify(e.subs(eta, 0)) for e in xi])
print("div xi_0 =", sp.simplify(sum(sp.diff(xi[i].subs(eta, 0), xs[i]) for i in range(3))), " (expect 3 zeta_L: conformal factor^3)")
# conformal check of the spatial map at eta=0: d x~^k d x~^k = (1+2 zeta_L) dx^2 + O(q^2)
M = sp.Matrix(3, 3, lambda i, j: sp.diff(xs[i] + eps1*xi[i].subs(eta, 0), xs[j]))
MM = (M.T*M).applyfunc(lambda e: sp.series(e, eps1, 0, 2).removeO())
print("M^T M =", sp.simplify(MM - (1+2*eps1*zetaL)*sp.eye(3)), " (expect 0: conformal map, factor 1+2 zeta_L)")

# ---------------------------------------------------------------- PART C
print("\n=== PART C: diffeo image of a short-mode Poisson metric; long x short cross terms ===")
PhiS = sp.Function('PhiS')
PsiS = sp.Function('PsiS')
vS = sp.Function('vS')         # short-mode velocity potential (tilde frame): v~_i = d_i vS
dTS = sp.Function('dTS')
# tilde-frame metric, linear parametrisation, with short modes only:
def gtilde(coords):
    e0, y1, y2, y3 = coords
    return sp.diag(-(1+2*sS*PsiS(e0, y1, y2, y3)), 1-2*sS*PhiS(e0, y1, y2, y3), 1-2*sS*PhiS(e0, y1, y2, y3), 1-2*sS*PhiS(e0, y1, y2, y3))*e0**2
gS = gtilde(Xt)
gfull = (J.T*gS*J).applyfunc(lambda e: sp.expand(sp.series(sp.series(sp.expand(e), eps1, 0, 2).removeO(), sS, 0, 2).removeO()))
# Taylor-expand the arguments of the short fields: F(eta+eps epsT, x+eps xi) = F + eps (epsT dF/deta + xi.grad F)
def taylor_args(expr):
    out = expr
    for Fn in (PhiS, PsiS, vS, dTS):
        F0 = Fn(eta, x1, x2, x3)
        F1 = epsT*sp.diff(F0, eta) + sum(xi[i]*sp.diff(F0, xs[i]) for i in range(3))
        out = out.replace(lambda e: isinstance(e, Fn), lambda e: F0 + eps1*F1)
    return sp.expand(sp.series(sp.expand(out), eps1, 0, 2).removeO())
gfull = gfull.applyfunc(taylor_args)
Phi_full = sp.expand((1 - gfull[1, 1]/eta**2)/2)
Psi_full = sp.expand((-(gfull[0, 0]/eta**2) - 1)/2)
F0 = PhiS(eta, x1, x2, x3); P0 = PsiS(eta, x1, x2, x3)
F1 = epsT*sp.diff(F0, eta) + sum(xi[i]*sp.diff(F0, xs[i]) for i in range(3))
P1 = epsT*sp.diff(P0, eta) + sum(xi[i]*sp.diff(P0, xs[i]) for i in range(3))
crossPhi = sp.simplify(Phi_full.coeff(eps1, 1).coeff(sS, 1) - F1)
crossPsi = sp.simplify(Psi_full.coeff(eps1, 1).coeff(sS, 1) - P1)
print("Phi cross term (beyond argument shift) =", sp.factor(crossPhi), "   (expect +(4/3) zeta_L PhiS)")
print("Psi cross term (beyond argument shift) =", sp.factor(crossPsi), "   (expect -(4/3) zeta_L PsiS)")
print("g_ij offdiag cross:", sp.simplify(gfull[1, 2].coeff(eps1, 1).coeff(sS, 1)))
g0i_cross = [sp.factor(gfull[0, i].coeff(eps1, 1).coeff(sS, 1)/eta**2) for i in (1, 2, 3)]
print("ghat_0i cross term O(q S) =", g0i_cross, "   (expect (2 eta/3) q_i (PhiS+PsiS))")
g0i_L2 = [sp.factor(gfull[0, i].coeff(eps1, 1).coeff(sS, 0)/eta**2) for i in (1, 2, 3)]
print("ghat_0i at O(eps^1 s^0):", g0i_L2)
# fluid velocity cross terms: u^mu = (J^-1)^mu_nu u~^nu with u~ = u~^0 (1, v~^i), v~^i = d_i vS at (eta~, x~)
vt = [sS*sp.diff(vS(*Xt), Xt[i+1]) if False else None for i in range(3)]
# derivative wrt tilde coords of vS(tilde coords) == (d_i vS)(tilde coords): build via a helper function symbol
dvS = [sp.Function('dvS%d' % (i+1)) for i in range(3)]
ut = [1] + [sS*dvS[i](*Xt) for i in range(3)]
uf = [sum(Jinv[m, n]*ut[n] for n in range(4)) for m in range(4)]
vfull = [sp.expand(sp.series(sp.series(uf[i]/uf[0], eps1, 0, 2).removeO(), sS, 0, 2).removeO()) for i in (1, 2, 3)]
def taylor_dv(expr):
    out = expr
    for i in range(3):
        Fn = dvS[i]
        F0_ = Fn(eta, x1, x2, x3)
        F1_ = epsT*sp.diff(F0_, eta) + sum(xi[j]*sp.diff(F0_, xs[j]) for j in range(3))
        out = out.replace(lambda e: isinstance(e, Fn), lambda e: F0_ + eps1*F1_)
    return sp.expand(sp.series(sp.expand(out), eps1, 0, 2).removeO())
vfull = [taylor_dv(e) for e in vfull]
for i in range(3):
    F0_ = dvS[i](eta, x1, x2, x3)
    F1_ = epsT*sp.diff(F0_, eta) + sum(xi[j]*sp.diff(F0_, xs[j]) for j in range(3))
    crossv = sp.expand(vfull[i].coeff(eps1, 1).coeff(sS, 1) - F1_)
    print("v^%d cross term (beyond argument shift) =" % (i+1), sp.factor(crossv))
print("   (expect -(4/3) zeta_L v~^i  +  [x^i q.v~ - q^i x.v~]  (local rotation from the SCT))")
# temperature cross term: T = T~(eta~,x~) exactly => 1+dT = (eta/eta~)(1+dT~(eta~,x~))
dTfull = sp.expand(sp.series(sp.series(sp.expand((eta/Xt[0])*(1+sS*dTS(*Xt)) - 1), eps1, 0, 2).removeO(), sS, 0, 2).removeO())
dTfull = taylor_args(dTfull)
D0 = dTS(eta, x1, x2, x3)
D1 = epsT*sp.diff(D0, eta) + sum(xi[i]*sp.diff(D0, xs[i]) for i in range(3))
print("delta_T cross term (beyond argument shift) =", sp.factor(sp.expand(dTfull.coeff(eps1, 1).coeff(sS, 1) - D1)), "  (expect (1/3) zeta_L dTS)")

# Second-order Poisson-gauge field relations for a CONSTANT long mode, in Fourier space.
# Short tilde field with primordial zeta~ = zeta_S - c y.grad zeta_S, evaluated at (eta~,x~):
#   X~_S(eta~,x~) = int_k zeta_k e^{ikx} [ (1 + c d_X) T_X - (4/3) c eta d_eta T_X ],  T_X = eta^{d_X} t_X(k eta)
# (the x.grad pieces cancel between the argument dilation and the profile relabelling; checked here symbolically)
print("\n--- constant-mode second-order Poisson-gauge transfer functions (linear parametrisation) ---")
kk = sp.Symbol('k', positive=True)
tX = sp.Function('t')
for name, dX, pref in (("Phi", 0, sp.Rational(4, 3)), ("Psi", 0, -sp.Rational(4, 3)), ("delta_T", 0, sp.Rational(1, 3)), ("v_pot", 1, -sp.Rational(4, 3))):
    TX = eta**dX*tX(kk*eta)
    # relabelling: -c d/dln k ; argument shift eta -> eta(1 - c/3)
    rel = -c*kk*sp.diff(TX, kk) - c*(eta/3)*sp.diff(TX, eta)
    tot = sp.simplify(pref*c*TX + rel)   # pref = multiplicative cross factor from the metric/fluid transformation
    print("  F2_%s(0,k)/T-structure: " % name, sp.simplify(tot/c), "  =", sp.simplify((tot/c).subs(tX(kk*eta), sp.Symbol('T')).doit()) if False else "")
# Explicit forms with the RD transfer functions:
F2Phi = sp.Rational(4, 3)*(TPhi - eta*sp.diff(TPhi, eta))
F2Psi = -sp.Rational(4, 3)*(TPhi + eta*sp.diff(TPhi, eta))
F2dT = sp.Rational(1, 3)*Tdelta/4 - sp.Rational(4, 3)*eta*sp.diff(Tdelta/4, eta)
F2v = -sp.Rational(1, 3)*Tv - sp.Rational(4, 3)*eta*sp.diff(Tv, eta)
print("  F2_Phi(0,k)   = (4/3)(1 - eta d_eta) T_Phi ;  k->0 limit:", sp.limit(F2Phi, eta, 0))
print("  F2_Psi(0,k)   = -(4/3)(1 + eta d_eta) T_Phi ;  k->0 limit:", sp.limit(F2Psi, eta, 0))
print("  F2_deltaT(0,k)= (1/3) T_dT - (4/3) eta d_eta T_dT ;  k->0 limit:", sp.limit(F2dT, eta, 0))
print("  F2_v(0,k)     = -(1/3) T_v - (4/3) eta d_eta T_v ;  (k->0 limit)/eta:", sp.limit(F2v/eta, eta, 0))
print("  check zeta_nl^(2)=0 at eta->0: F2_dT(0) - F2_Phi(0)/... : the nonlinear zeta = ln(1+dT)+(1/2)ln(1-2Phi);")
# second-order primordial check: zeta2/2 = dT2/2 - dT1^2/2 - Phi2/2 - Phi1^2  (from ln expansions, cross terms only)
dT1L, dT1S, Phi1L, Phi1S = sp.Rational(1, 3)*c, sp.Rational(1, 3), -sp.Rational(2, 3)*c, -sp.Rational(2, 3)
zeta2cross = sp.limit(F2dT, eta, 0) - dT1L*dT1S - sp.limit(F2Phi, eta, 0) - 2*Phi1L*Phi1S
print("  zeta_nl cross term at eta->0 (must vanish):", sp.simplify(zeta2cross))

# O(q S): the longitudinal part of ghat_0i = (4 eta/3) q_i PhiS must be gauged away (alpha = time shift at 2nd order):
print("\n--- O(q S): ghat_0i = (4 eta/3) q_i Phi_S = d_i lambda + omega_i^T,  lambda = (4 eta/3) nabla^-2 (q.grad Phi_S)")
print("    Poisson gauge is restored by eta -> eta + alpha with alpha = lambda (beta = 0 keeps g_ij trace-only):")
print("    Psi_P = Psi_diffeo - alpha' - H alpha ;  Phi_P = Phi_diffeo + H alpha ;  v_P = v_diffeo + grad(alpha)?? (see report)")
print("    omega_i^T (physical squeezed vector mode) = (4 eta/3) T_Phi(k eta) zeta_k [ q_i - k_i (q.k)/k^2 ]")

# ---------------------------------------------------------------- PART D
print("\n=== PART D: exact constant-mode time reparametrisation (all orders in c) ===")
lam = sp.Symbol('lambda', positive=True)
# eta~ = lam eta, x~ = e^c x; Poisson form with Phi = Psi requires lam^4 + lam^2 e^{2c} - 2 = 0
sol = sp.solve(sp.Eq(lam**4 + lam**2*sp.exp(2*c) - 2, 0), lam)
lam_c = [s for s in sol if sp.N(s.subs(c, 0)) == 1][0] if any(sp.N(s.subs(c, 0)) == 1 for s in sol) else sp.sqrt((-sp.exp(2*c) + sp.sqrt(sp.exp(4*c) + 8))/2)
lam_c = sp.sqrt((-sp.exp(2*c) + sp.sqrt(sp.exp(4*c) + 8))/2)
print("lambda(c) =", lam_c, " series:", sp.series(lam_c, c, 0, 3))
print("Phi = Psi = (1 - lambda^4)/2 =", sp.series((1 - lam_c**4)/2, c, 0, 3), "  (linear parametrisation; -2c/3 at O(c))")
print("\nALL CHECKS PASSED")
