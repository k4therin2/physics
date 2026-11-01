"""
einstein.py -- Einstein tensor and fluid stress tensor to O(eps^3) in Poisson gauge, RD, a(eta)=eta.
g_{mu nu} = eta^2 ghat_{mu nu};  ghat_00=-(1+2Psi), ghat_0i=omega_i, ghat_ij=(1-2Phi)delta_ij+h_ij.
Fluid: rho = rhobar(1+delta), P=rho/3, 8 pi G rhobar = 3/eta^4, u^mu = eta^{-1} uhat^0 (1, v^i).
E_{mu nu} = G_{mu nu} - 8 pi G T_{mu nu}.
Output: E.pkl (Poly), E_sympy.pkl (sympy srepr strings), E.tex
"""
import time, pickle, sys, os
import sympy as sp
from tpoly import *

t0 = time.time()
def log(*a):
    print("[%7.1fs]" % (time.time() - t0), *a); sys.stdout.flush()

# ---------------- fields
Psi = field_series('Psi'); Phi = field_series('Phi'); delta = field_series('delta')
om = [field_series('om%d' % i, orders=(2, 3)) for i in (1, 2, 3)]
v = [field_series('v%d' % i) for i in (1, 2, 3)]
hnames = {(0, 0): 'h11', (0, 1): 'h12', (0, 2): 'h13', (1, 1): 'h22', (1, 2): 'h23', (2, 2): 'h33'}
def hcomp(i, j):
    return field_series(hnames[(min(i, j), max(i, j))], orders=(2, 3))

# conformal metric ghat
gh = [[None] * 4 for _ in range(4)]
gh[0][0] = -(const(1) + Psi * 2)
for i in range(3):
    gh[0][i + 1] = gh[i + 1][0] = om[i]
    for j in range(3):
        gh[i + 1][j + 1] = ((const(1) - Phi * 2) if i == j else const(0)) + hcomp(i, j)
# impose tracelessness h33 = -h11-h22 ? NO: keep h33 independent here; tracelessness/transversality imposed later
# (we keep all six h components so the equations are valid for any symmetric h; constraints imposed at decomposition stage)

# physical metric g = eta^2 ghat
g = [[gh[a][b] * eta_pow(2) for b in range(4)] for a in range(4)]

# inverse of ghat as series: ghat = eta0 + H, ghat^{-1} = eta0^{-1} sum_n (-H eta0^{-1})^n
eta0 = [[const(-1 if a == 0 else 1) if a == b else const(0) for b in range(4)] for a in range(4)]
eta0inv = eta0  # Minkowski inverse = itself
H = [[gh[a][b] - eta0[a][b] for b in range(4)] for a in range(4)]
def matmul(A, B):
    return [[sum((A[a][c] * B[c][b] for c in range(4)), Poly({})) for b in range(4)] for a in range(4)]
M = matmul(H, eta0inv)            # H eta0^{-1}
term = [[const(1) if a == b else const(0) for b in range(4)] for a in range(4)]
ghinv = [[const(1) if a == b else const(0) for b in range(4)] for a in range(4)]
negM = [[-M[a][b] for b in range(4)] for a in range(4)]
for n in range(1, MAXORD + 1):
    term = matmul(term, negM)
    ghinv = [[ghinv[a][b] + term[a][b] for b in range(4)] for a in range(4)]
ghinv = matmul(eta0inv, ghinv)    # eta0^{-1} (1 - M + M^2 - ...)
# check: ghat * ghinv = identity to O(eps^3)
chk = matmul(gh, ghinv)
for a in range(4):
    for b in range(4):
        diff = chk[a][b] - (const(1) if a == b else const(0))
        assert diff.is_zero(), (a, b, diff.to_sympy())
log("inverse metric OK")
ginv = [[ghinv[a][b] * eta_pow(-2) for b in range(4)] for a in range(4)]

# ---------------- Christoffels of g (direct)
dg = [[[g[a][b].diff(COORDS[c]) for c in range(4)] for b in range(4)] for a in range(4)]  # dg[a][b][c] = d_c g_ab
Gam = [[[None] * 4 for _ in range(4)] for _ in range(4)]
for a in range(4):
    for b in range(4):
        for c in range(b, 4):
            s = Poly({})
            for d in range(4):
                if ginv[a][d].is_zero():
                    continue
                s = s + ginv[a][d] * (dg[d][c][b] + dg[d][b][c] - dg[b][c][d])
            Gam[a][b][c] = Gam[a][c][b] = s / 2
log("Christoffels done; max terms", max(Gam[a][b][c].nterms() for a in range(4) for b in range(4) for c in range(4)))

# ---------------- Ricci
Ric = [[None] * 4 for _ in range(4)]
for b in range(4):
    for d in range(b, 4):
        s = Poly({})
        for a in range(4):
            s = s + Gam[a][b][d].diff(COORDS[a]) - Gam[a][b][a].diff(COORDS[d])
            for e in range(4):
                s = s + Gam[a][a][e] * Gam[e][b][d] - Gam[a][d][e] * Gam[e][b][a]
        Ric[b][d] = Ric[d][b] = s
        log("Ricci", b, d, s.nterms())
Rs = sum((ginv[a][b] * Ric[a][b] for a in range(4) for b in range(4)), Poly({}))
G = [[Ric[a][b] - g[a][b] * Rs / 2 for b in range(4)] for a in range(4)]
log("Einstein tensor done")

# ---------------- fluid
uh_up = [const(1)] + v  # uhat^mu / uhat^0
norm = sum((gh[a][b] * uh_up[a] * uh_up[b] for a in range(4) for b in range(4)), Poly({}))  # = -(uhat^0)^{-2}
uh0 = (-norm).inverse().sqrt_series()       # uhat^0
uhat = [uh0 * uh_up[a] for a in range(4)]   # uhat^mu
u_up = [uhat[a] * eta_pow(-1) for a in range(4)]     # u^mu = uhat^mu / a
u_dn = [sum((g[a][b] * u_up[b] for b in range(4)), Poly({})) for a in range(4)]
# check normalisation g u u = -1
nn = sum((g[a][b] * u_up[a] * u_up[b] for a in range(4) for b in range(4)), Poly({}))
assert (nn + const(1)).is_zero(), nn.to_sympy()
log("u normalisation OK")
# 8 pi G T_mu nu = (3/eta^4)(1+delta) [ (4/3) u_mu u_nu + (1/3) g_mu nu ]
pref = (const(1) + delta) * eta_pow(-4, 3)
T8 = [[pref * (u_dn[a] * u_dn[b] * Fraction(4, 3) + g[a][b] * Fraction(1, 3)) for b in range(4)] for a in range(4)]
E = [[G[a][b] - T8[a][b] for b in range(4)] for a in range(4)]
log("E done")

# background check
for a in range(4):
    for b in range(4):
        bg = E[a][b].order_part(0)
        assert bg.is_zero(), ("background not solved", a, b, bg.to_sympy())
log("BACKGROUND SOLVED EXACTLY (order-0 part of all E_{mu nu} vanishes)")

# also conservation: nabla_mu T^{mu nu}, computed with T^{mu nu} = ginv T ginv
Tup = [[sum((ginv[a][c] * T8[c][d] * ginv[d][b] for c in range(4) for d in range(4)), Poly({})) for b in range(4)] for a in range(4)]
cons = []
for nu in range(4):
    s = Poly({})
    for mu in range(4):
        s = s + Tup[mu][nu].diff(COORDS[mu])
        for lam in range(4):
            s = s + Gam[mu][mu][lam] * Tup[lam][nu] + Gam[nu][mu][lam] * Tup[mu][lam]
    cons.append(s)
    assert s.order_part(0).is_zero()
log("conservation equations computed; background OK")

out = {'E': E, 'cons': cons, 'uhat0': uh0, 'u_dn': u_dn, 'Gam': Gam, 'G': G, 'T8': T8}
save_polys(out, 'E.pkl')
log("saved E.pkl")

# sympy / latex dumps per order
with open('E_sympy.pkl', 'wb') as fh:
    blob = {}
    for a in range(4):
        for b in range(a, 4):
            for n in (1, 2, 3):
                blob[('E', a, b, n)] = sp.srepr(E[a][b].order_part(n).to_sympy())
    for nu in range(4):
        for n in (1, 2, 3):
            blob[('cons', nu, n)] = sp.srepr(cons[nu].order_part(n).to_sympy())
    pickle.dump(blob, fh)
with open('E.tex', 'w') as fh:
    fh.write("% Einstein + fluid equations E_{mu nu} = G_{mu nu} - 8 pi G T_{mu nu}, coefficient of eps^n (includes 1/n!)\n")
    for a in range(4):
        for b in range(a, 4):
            for n in (1, 2, 3):
                ex = E[a][b].order_part(n).to_sympy()
                fh.write("\n\\subsection*{$E_{%d%d}$ at order $\\epsilon^%d$ (%d terms)}\n" % (a, b, n, len(ex.args)))
                fh.write("\\begin{dmath}\n" + sp.latex(ex) + "\n\\end{dmath}\n")
log("wrote E.tex;  term counts:")
for a in range(4):
    for b in range(a, 4):
        print("E_%d%d:" % (a, b), [E[a][b].order_part(n).nterms() for n in (1, 2, 3)])
