"""
N1: sympy derivation of the exact Einstein + radiation-fluid equations in PLANAR symmetry,
comoving-orthogonal (fluid-rest) gauge.   (see README.md / report for the physics summary)

Metric (physical):  ds^2 = -alpha^2 dt^2 + A^2 dx^2 + B^2 (dy^2 + dz^2),  fields (t,x) only.
Fluid: u^mu = (1/alpha, 0, 0, 0) (comoving, hypersurface orthogonal: OK since irrotational).
Radiation: P = rho/3, rho = alpha_rad T^4.  Units: 8 pi G = 1, conformal background temperature
That_bar = 1, Friedmann 8 pi G a^4 rhobar = 3  =>  8 pi G rho = 3 That^4 / t^4, That = t*T.

Derived below (exact):
  Euler  : d_x ln alpha = -d_x ln T  => alpha T = F(t); time relabelling fixes F=1: alpha = 1/T.
  Energy : d_t ln T = -(1/3) d_t ln(A B^2)  => T (A B^2)^{1/3} = C(x) = e^{zeta(x)} (regular sol.)
  Hence  alpha = e^{-zeta(x)} (A B^2)^{1/3} EXACTLY, and 8 pi G rho = 3/alpha^4.
Hatted variables: A = t Ah, B = t Bh, alpha = t alh, alh = e^{-zeta}(Ah Bh^2)^{1/3}, That = 1/alh.
Momenta: pA = (1/alh) d_t ln Ah, pB = (1/alh) d_t ln Bh  (= -K^x_x, -K^y_y of conformal metric).
"""
import sympy as sp
import pickle

t, x, y, z = sp.symbols('t x y z', positive=True)
A = sp.Function('A')(t, x)
B = sp.Function('B')(t, x)
al = sp.Function('alpha')(t, x)
rho = sp.Function('rho')(t, x)
coords = [t, x, y, z]
g = sp.diag(-al**2, A**2, B**2, B**2)
ginv = g.inv()

def christoffel(g, ginv, coords):
    n = len(coords)
    Gam = [[[0]*n for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            for c in range(n):
                Gam[a][b][c] = sp.simplify(sum(ginv[a, d]*(sp.diff(g[d, b], coords[c]) + sp.diff(g[d, c], coords[b]) - sp.diff(g[b, c], coords[d])) for d in range(n))/2)
    return Gam

def ricci(Gam, coords):
    n = len(coords)
    Rm = sp.zeros(n, n)
    for b in range(n):
        for c in range(n):
            expr = 0
            for a in range(n):
                expr += sp.diff(Gam[a][b][c], coords[a]) - sp.diff(Gam[a][b][a], coords[c])
                for d in range(n):
                    expr += Gam[a][a][d]*Gam[d][b][c] - Gam[a][c][d]*Gam[d][b][a]
            Rm[b, c] = sp.simplify(expr)
    return Rm

print("Christoffels / Ricci ...", flush=True)
Gam = christoffel(g, ginv, coords)
Ric = ricci(Gam, coords)
Rs = sp.simplify(sum(ginv[a, b]*Ric[a, b] for a in range(4) for b in range(4)))
G = sp.simplify(Ric - Rs*g/2)

u_up = sp.Matrix([1/al, 0, 0, 0])
u_dn = g*u_up
p = rho/3
T_dn = (rho + p)*u_dn*u_dn.T + p*g
EFE = sp.simplify(G - T_dn)                 # G_{mu nu} - 8 pi G T_{mu nu}, 8 pi G = 1
T_up = (rho + p)*u_up*u_up.T + p*ginv
divT = []
for nu in range(4):
    expr = 0
    for mu in range(4):
        expr += sp.diff(T_up[mu, nu], coords[mu])
        for lam in range(4):
            expr += Gam[mu][mu][lam]*T_up[lam, nu] + Gam[nu][mu][lam]*T_up[mu, lam]
    divT.append(sp.simplify(expr))
chk_energy = sp.simplify(divT[0]*al**2/rho - (sp.diff(rho, t)/rho + sp.Rational(4, 3)*(sp.diff(A, t)/A + 2*sp.diff(B, t)/B)))
chk_euler = sp.simplify(divT[1]*A**2/rho - (sp.Rational(4, 3)*sp.diff(al, x)/al + sp.Rational(1, 3)*sp.diff(rho, x)/rho))
assert chk_energy == 0 and chk_euler == 0
print("energy eq: d_t ln rho = -(4/3) d_t ln(AB^2);  Euler: d_x ln alpha = -(1/4) d_x ln rho  [verified]")

# ---- shear (projected) ----
nab_u = sp.zeros(4, 4)
for m in range(4):
    for n in range(4):
        nab_u[m, n] = sp.diff(u_dn[n], coords[m]) - sum(Gam[l][m][n]*u_dn[l] for l in range(4))
h_mix = sp.eye(4) + u_up*u_dn.T          # h^mu_nu
h_dn = g + u_dn*u_dn.T
theta = sp.simplify(sum(ginv[a, b]*nab_u[a, b] for a in range(4) for b in range(4)))
sym = (nab_u + nab_u.T)/2
proj = sp.simplify(h_mix.T*sym*h_mix)    # h^a_mu h^b_nu sym_ab
sigma = sp.simplify(proj - theta*h_dn/3)
s_sym = (sp.diff(A, t)/A - sp.diff(B, t)/B)/al
assert sp.simplify(sigma[1, 1] - sp.Rational(2, 3)*A**2*s_sym) == 0
assert sp.simplify(sigma[2, 2] + sp.Rational(1, 3)*B**2*s_sym) == 0
assert all(sp.simplify(sigma[i, j]) == 0 for i in range(4) for j in range(4) if i != j)
assert sp.simplify(sigma[0, 0]) == 0
acc = sp.simplify(u_up.T*nab_u)   # a_nu = u^mu nabla_mu u_nu
print("shear verified: sigma_xx = 2/3 A^2 s, sigma_yy = sigma_zz = -1/3 B^2 s, s = (d_t lnA - d_t lnB)/alpha;",
      " theta =", sp.simplify(theta), "; a_x =", sp.simplify(acc[1]))

# ---- replace by symbols (hatted variables) ----
Ah, Bh, alh, Ax, Bx, Axx, Bxx, pA, pB, pAx, pBx = sp.symbols('Ah Bh alh Ax Bx Axx Bxx pA pB pAx pBx', real=True)
Att, Btt = sp.symbols('Att Btt')
ze, zx, zxx = sp.symbols('zeta zx zxx', real=True)
# alh = e^{-zeta}(Ah Bh^2)^{1/3}: exact derivative identities
At = alh*Ah*pA
Bt = alh*Bh*pB
alh_x = alh*(-zx + (Ax/Ah + 2*Bx/Bh)/3)
alh_xx = alh_x*(-zx + (Ax/Ah + 2*Bx/Bh)/3) + alh*(-zxx + (Axx/Ah - Ax**2/Ah**2 + 2*Bxx/Bh - 2*Bx**2/Bh**2)/3)
alh_t = alh*(At/Ah + 2*Bt/Bh)/3
Atx = alh_x*Ah*pA + alh*Ax*pA + alh*Ah*pAx
Btx = alh_x*Bh*pB + alh*Bx*pB + alh*Bh*pBx
alh_tx = alh_x*(At/Ah + 2*Bt/Bh)/3 + alh*(Atx/Ah - At*Ax/Ah**2 + 2*Btx/Bh - 2*Bt*Bx/Bh**2)/3
alh_tt = sp.Symbol('alh_tt')  # placeholder, substituted after Att,Btt known (appears? check)
rep = [
    (sp.Derivative(A, t, t), 2*At + t*Att), (sp.Derivative(B, t, t), 2*Bt + t*Btt),
    (sp.Derivative(A, t, x), Ax + t*Atx), (sp.Derivative(B, t, x), Bx + t*Btx),
    (sp.Derivative(A, x, x), t*Axx), (sp.Derivative(B, x, x), t*Bxx),
    (sp.Derivative(A, t), Ah + t*At), (sp.Derivative(B, t), Bh + t*Bt),
    (sp.Derivative(A, x), t*Ax), (sp.Derivative(B, x), t*Bx),
    (sp.Derivative(al, t, t), 2*alh_t + t*alh_tt),
    (sp.Derivative(al, t, x), alh_x + t*alh_tx),
    (sp.Derivative(al, x, x), t*alh_xx),
    (sp.Derivative(al, t), alh + t*alh_t), (sp.Derivative(al, x), t*alh_x),
    (A, t*Ah), (B, t*Bh), (al, t*alh), (rho, 3/(t*alh)**4)]
def R(e):
    for a, b in rep:
        e = e.subs(a, b)
    return e
Ham = sp.simplify(R(EFE[0, 0]))
Mom = sp.simplify(R(EFE[0, 1]))
Exx = sp.simplify(R(EFE[1, 1]))
Eyy = sp.simplify(R(EFE[2, 2]))
assert sp.simplify(R(EFE[3, 3]) - Eyy) == 0
for nm, e in [('Ham', Ham), ('Mom', Mom), ('Exx', Exx), ('Eyy', Eyy)]:
    print(nm, 'has Att,Btt,alh_tt:', e.has(Att), e.has(Btt), e.has(alh_tt))
assert not (Ham.has(Att) or Ham.has(Btt) or Mom.has(Att) or Mom.has(Btt))
# alh_tt = d_t[ alh (At/Ah + 2 Bt/Bh)/3 ]; At = alh Ah pA -> At/Ah = alh pA ; d_t(alh pA) = alh_t pA + alh pA_t
# so alh_tt = alh_t*(alh pA + 2 alh pB)/3 + alh*(alh_t (pA+2pB) + alh (pAt + 2 pBt))/3 -> couples to pA_t, pB_t.
pAt, pBt = sp.symbols('pAt pBt')
alh_tt_e = alh_t*alh*(pA + 2*pB)/3 + alh*(alh_t*(pA + 2*pB) + alh*(pAt + 2*pBt))/3
# Att = d_t(alh Ah pA) = alh_t Ah pA + alh At pA + alh Ah pAt
Att_e = alh_t*Ah*pA + alh*At*pA + alh*Ah*pAt
Btt_e = alh_t*Bh*pB + alh*Bt*pB + alh*Bh*pBt
Exx2 = Exx.subs({Att: Att_e, Btt: Btt_e, alh_tt: alh_tt_e})
Eyy2 = Eyy.subs({Att: Att_e, Btt: Btt_e, alh_tt: alh_tt_e})
sol = sp.solve([Exx2, Eyy2], [pAt, pBt], dict=True)
assert len(sol) == 1
pA_t = sp.simplify(sol[0][pAt])
pB_t = sp.simplify(sol[0][pBt])
# scale constraints nicely: Ham ~ G_tt - T_tt ; divide by alpha^2 to get energy-density-like units
Ham_n = sp.simplify(Ham/(t*alh)**2)
Mom_n = sp.simplify(Mom)
print("\n d_t pA =", pA_t)
print("\n d_t pB =", pB_t)
print("\n Ham/alpha^2 =", Ham_n)
print("\n Mom =", Mom_n)

syms = dict(t=t, Ah=Ah, Bh=Bh, alh=alh, Ax=Ax, Bx=Bx, Axx=Axx, Bxx=Bxx, pA=pA, pB=pB, pAx=pAx, pBx=pBx, zeta=ze, zx=zx, zxx=zxx)
with open('eqs_generated.pkl', 'wb') as f:
    pickle.dump({'pA_t': pA_t, 'pB_t': pB_t, 'Ham': Ham_n, 'Mom': Mom_n, 'symbols': syms}, f)
print("saved eqs_generated.pkl")
