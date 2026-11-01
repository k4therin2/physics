"""Milestone 2: the null geodesic of ghat through second order (SPEC §3, §6).
Parametrisation by eta: k^mu = dx^mu/deta = (1, k^i), background k^i = -n^i (propagation e = -n),
background ray  xbar(eta) = xbar_o + n (etabar_o - eta).
Non-affine geodesic equation:   dk^i/deta = G^i := -Gam^i_{ab} k^a k^b + Gam^0_{ab} k^a k^b k^i.
Perturbed ray through the perturbed observer event (eta_o, x_o) = (etabar_o + deta_o, xbar_o + dx_o):
   x(eta) = x_o + n (eta_o - eta) + int_{eta_o}^{eta} dk   =>   Dx(eta) := x(eta) - xbar(eta)
   Dx(eta) = dx_o + n deta_o - int_eta^{eta_o} dk(eta') deta'.
Everything is expressed through functions ON THE BACKGROUND RAY (fields at (eta, xbar(eta))):
   k1(eta)   = k1_init|_obar + int_eta^{etabar_o} G1 deta'
   Dx1(eta)  = dx_o1 + n deta_o1 - int_eta^{etabar_o} k1 deta'
   k2(eta)/2 = k2_init/2 + deta_o1 d_eta k1_init + dx_o1^i d_i k1_init - deta_o1 G1(etabar_o)
               + int_eta^{etabar_o} [ G2/2 + Dx1^i d_i G1 ] deta'      (G2/2 contains the k1-dependent pieces)
   Dx2(eta)/2 = dx_o2/2 + n deta_o2/2 - int_eta^{etabar_o} k2/2 deta' - deta_o1 k1(etabar_o)
Initial data (SPEC §6): phat_o ∝ uhat_o - n^a e_a, canonical triad e_a = Gram-Schmidt of (d_1,d_2,d_3) wrt ghat,
orthogonal to uhat:  k_init^i = (uhat^i - n^a e_a^i)/(uhat^0 - n^a e_a^0).
"""
import sympy as sp, pickle, time
from ser import Ser, N, eps, dot4
from fields import build, X, eta, x, y, z, F, fname, v

d = build(verbose=False)
g, ginv, Gam, fl = d['g'], d['ginv'], d['Gam'], d['fl']
uu, ul = fl['uu'], fl['ul']
n = sp.symbols('n1:4')
nn = sum(ni**2 for ni in n)

def unit_rules(e):
    """impose |n|=1: reduce e (polynomial in n3) modulo n3^2 - (1 - n1^2 - n2^2)."""
    e = sp.expand(e)
    if not e.has(n[2]): return e
    return sp.expand(sp.rem(sp.Poly(e, n[2]), sp.Poly(n[2]**2 - (1 - n[0]**2 - n[1]**2), n[2])).as_expr())

# ---------------- canonical triad (Gram-Schmidt, to order N) ----------------
def vec_const(i):  # coordinate basis vector d_i, i=1..3
    return [Ser.const(1 if m == i else 0) for m in range(4)]
def gdot(a, b): return dot4(a, b, g)
def vadd(a, b, s=1): return [a[m] + s*b[m] for m in range(4)]
def vscale(a, s): return [a[m]*s for m in range(4)]
tri = []
for i in (1, 2, 3):
    w = vec_const(i)
    w = vadd(w, vscale(uu, gdot(w, uu)))          # project orthogonal to u: w + (g(w,u)) u
    for e_ in tri:
        w = vadd(w, vscale(e_, gdot(w, e_)), -1)   # subtract projections on previous e's
    nrm = gdot(w, w).sqrt()
    tri.append(vscale(w, nrm.inv()))
triad = tri  # triad[a][mu] = e_a^mu as Ser

def check_triad():
    ok = True
    for a in range(3):
        ok &= gdot(triad[a], uu).iszero()
        for b in range(3):
            ok &= (gdot(triad[a], triad[b]) - (1 if a == b else 0)).iszero()
    return ok

# ---------------- direction map ----------------
num = [uu[m] - sum((triad[a][m]*n[a] for a in range(3)), Ser.const(0)) for m in range(4)]
k_init = [num[i+1]/num[0] for i in range(3)]   # k^i(eta_o) as Ser in fields at the observer event
# null check (needs |n|=1)
def check_null():
    kk = [Ser.const(1)] + k_init
    gk = gdot(kk, kk)
    return all(unit_rules(gk.c[m]) == 0 for m in range(N+1))

# ---------------- geodesic RHS with abstract k_n(eta) ----------------
kfun = {1: [sp.Function('k1_%d' % i)(eta) for i in (1, 2, 3)],
        2: [sp.Function('k2_%d' % i)(eta) for i in (1, 2, 3)]}
kser = [Ser.const(1)] + [Ser([-n[i], kfun[1][i], kfun[2][i]/2, 0]) for i in range(3)]
def contract2(Tm):
    return sum((Tm[m][q]*kser[m]*kser[q] for m in range(4) for q in range(4) if not Tm[m][q].iszero()), Ser.const(0))
Gam0kk = contract2(Gam[0])
Gvec = [-contract2(Gam[i+1]) + Gam0kk*kser[i+1] for i in range(3)]   # dk^i/deta, Ser (orders 1..3 used)
G1 = [Gvec[i].c[1] for i in range(3)]              # pure first-order fields
G2half = [Gvec[i].c[2] for i in range(3)]          # (1/2)G2: second-order fields at bg point + k1 x first-order
G3sixth = [Gvec[i].c[3] for i in range(3)]         # (1/6)G3 (k2 and k1k1 pieces), not needed for L3 (enters Dx3 only)

# ---------------- symbolic functionals (documentation / explicit integral forms) ----------------
etap = sp.Symbol('etap'); eta_o = sp.Symbol('eta_o_bar'); eta_s = sp.Symbol('eta_star')
xo = sp.symbols('xo1:4')
deta_o1, deta_o2 = sp.symbols('deta_o1 deta_o2')
dxo1 = sp.symbols('dxo1_1:4'); dxo2 = sp.symbols('dxo2_1:4')
Dxfun = {1: [sp.Function('Dx1_%d' % i)(eta) for i in (1, 2, 3)],
         2: [sp.Function('Dx2_%d' % i)(eta) for i in (1, 2, 3)]}
def on_ray(e, tau=None):
    """evaluate an expression in (eta,x,y,z) on the background ray at time tau (default: eta)."""
    tau = eta if tau is None else tau
    return e.subs({x: xo[0] + n[0]*(eta_o - tau), y: xo[1] + n[1]*(eta_o - tau), z: xo[2] + n[2]*(eta_o - tau)}, simultaneous=True).subs(eta, tau)
at_obs = {x: xo[0], y: xo[1], z: xo[2], eta: eta_o}
k1_init = [k_init[i].c[1] for i in range(3)]
k2_init_half = [k_init[i].c[2] for i in range(3)]
definitions = {}
definitions['k1'] = [sp.Eq(kfun[1][i], k1_init[i].subs(at_obs) + sp.Integral(on_ray(G1[i], etap), (etap, eta, eta_o))) for i in range(3)]
definitions['Dx1'] = [sp.Eq(Dxfun[1][i], dxo1[i] + n[i]*deta_o1 - sp.Integral(kfun[1][i].subs(eta, etap), (etap, eta, eta_o))) for i in range(3)]
G2_total_half = [G2half[i] + sum(Dxfun[1][j]*sp.diff(G1[i], X[j+1]) for j in range(3)) for i in range(3)]
definitions['k2half'] = [sp.Eq(kfun[2][i]/2,
        k2_init_half[i].subs(at_obs) + deta_o1*sp.diff(k1_init[i], eta).subs(at_obs)
        + sum(dxo1[j]*sp.diff(k1_init[i], X[j+1]).subs(at_obs) for j in range(3)) - deta_o1*G1[i].subs(at_obs)
        + sp.Integral(on_ray(G2_total_half[i], etap).subs({kfun[1][j]: kfun[1][j].subs(eta, etap) for j in range(3)}).subs({Dxfun[1][j]: Dxfun[1][j].subs(eta, etap) for j in range(3)}), (etap, eta, eta_o))) for i in range(3)]
definitions['Dx2half'] = [sp.Eq(Dxfun[2][i]/2, dxo2[i]/2 + n[i]*deta_o2/2 - sp.Integral(kfun[2][i].subs(eta, etap)/2, (etap, eta, eta_o)) - deta_o1*kfun[1][i].subs(eta, eta_o)) for i in range(3)]

if __name__ == '__main__':
    t0 = time.time()
    print('triad orthonormal & orthogonal to u through order %d:' % N, check_triad())
    print('k_init null (|n|=1) through order %d:' % N, check_null())
    print('e_1^mu to 2nd order:')
    for m in range(4): print('   e_1^%d = ' % m, triad[0][m].c[0], '+ eps*(', triad[0][m].c[1], ') + eps^2*(', triad[0][m].c[2], ')')
    print('e_2^mu first order:', [triad[1][m].c[1] for m in range(4)])
    print('e_2^mu second order (x2):', [sp.factor(2*triad[1][m].c[2]) for m in range(4)])
    print('k1_init^i:'); [print('   ', sp.factor(unit_rules(k1_init[i]))) for i in range(3)]
    print('k2_init^i/2 (term counts):', [len(sp.Add.make_args(k2_init_half[i])) for i in range(3)])
    print('k2_init^1/2 =', unit_rules(k2_init_half[0]))
    print('G1^i (dk1/deta):'); [print('   ', sp.collect(G1[i], n)) for i in range(3)]
    print('G2half term counts:', [len(sp.Add.make_args(G2half[i])) for i in range(3)])
    with open('ray_cache.pkl', 'wb') as f:
        pickle.dump(dict(triad=triad, k_init=k_init, G1=G1, G2half=G2half, G3sixth=G3sixth, definitions=definitions, n=n), f)
    print('t=%.1fs' % (time.time()-t0))
