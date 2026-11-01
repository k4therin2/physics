"""Abstract Poisson-gauge fields (SPEC §4), conformal metric ghat, fluid uhat,
Christoffels, shear sigmahat, radiation stress tensor and the fluid conservation
equations solved as substitution rules for d_eta delta_T and d_eta v^i, order by order.
All quantities are truncated eps-series (ser.Ser) with coefficients c_n = X_n/n!.
Index 0 = eta (conformal time), 1..3 = x,y,z. Signature (-+++). Background ghat = Minkowski.
"""
import sympy as sp, pickle, os, time
from ser import Ser, N, eps, dot4

eta, x, y, z = sp.symbols('eta x y z')
X = (eta, x, y, z)
def F(name): return sp.Function(name)(*X)

fac = {1: sp.Integer(1), 2: sp.Rational(1, 2), 3: sp.Rational(1, 6)}
def field(name, start=1):
    c = [sp.Integer(0)]*(N+1)
    for n in range(start, N+1):
        c[n] = fac[n]*F('%s%d' % (name, n))
    return Ser(c)

# scalar potentials, temperature contrast
Psi = field('Psi'); Phi = field('Phi'); dT = field('dT')
# vector omega_i (starts at eps^2), tensor h_ij (starts at eps^2, traceless imposed via h33=-h11-h22)
om = [field('om_%d_' % i, start=2) for i in (1, 2, 3)]
def _h(i, j):
    i, j = sorted((i, j)); return field('h_%d%d_' % (i, j), start=2)
h = [[None]*3 for _ in range(3)]
for i in range(3):
    for j in range(3):
        h[i][j] = _h(i+1, j+1)
h[2][2] = -(h[0][0] + h[1][1])
# velocity v^i
v = [field('v_%d_' % i) for i in (1, 2, 3)]
# convenient names of the order-n field functions
def fname(base, n, *idx):
    return F(base + ('%d' % n if not idx else '_' + ''.join(str(i) for i in idx) + '_%d' % n))

# ---------- metric ----------
delta = sp.eye(3)
g = [[None]*4 for _ in range(4)]
g[0][0] = -(1 + 2*Psi)
for i in range(3):
    g[0][i+1] = g[i+1][0] = om[i]
    for j in range(3):
        g[i+1][j+1] = (1 - 2*Phi)*delta[i, j] + h[i][j]
eta_mink = sp.diag(-1, 1, 1, 1)

def matmul(A, B):
    return [[sum((A[i][k]*B[k][j] for k in range(4)), Ser.const(0)) for j in range(4)] for i in range(4)]

def inverse_metric(g):
    # g = eta + dg ; g^{-1} = sum_k (-eta^{-1} dg)^k eta^{-1}
    dg = [[g[i][j] - Ser.const(eta_mink[i, j]) for j in range(4)] for i in range(4)]
    etainv = [[Ser.const(eta_mink[i, j]) for j in range(4)] for i in range(4)]
    M = [[-sum((etainv[i][k]*dg[k][j] for k in range(4)), Ser.const(0)) for j in range(4)] for i in range(4)]
    res = [[Ser.const(eta_mink[i, j]) for j in range(4)] for i in range(4)]
    P = [[Ser.const(1 if i == j else 0) for j in range(4)] for i in range(4)]
    for k in range(1, N+1):
        P = matmul(P, M)
        PE = matmul(P, etainv)
        res = [[res[i][j] + PE[i][j] for j in range(4)] for i in range(4)]
    return res

def christoffel(g, ginv):
    dg = [[[g[b][c].diff(X[a]) for c in range(4)] for b in range(4)] for a in range(4)]  # dg[a][b][c] = d_a g_bc
    Gam = [[[None]*4 for _ in range(4)] for _ in range(4)]
    for a in range(4):
        for b in range(4):
            for c in range(b, 4):
                s = Ser.const(0)
                for d in range(4):
                    if ginv[a][d].iszero(): continue
                    s = s + ginv[a][d]*(dg[b][d][c] + dg[c][d][b] - dg[d][b][c])
                Gam[a][b][c] = Gam[a][c][b] = s*sp.Rational(1, 2)
    return Gam

def fluid(g, ginv, Gam):
    # uhat^mu = u0 (1, v^i); u0 from normalisation
    vv = [Ser.const(1)] + v
    nrm = -dot4(vv, vv, g)              # = 1 + O(eps)
    u0 = nrm.apply(1/sp.sqrt(sp.Symbol('t_apply')))
    uu = [u0*vv[m] for m in range(4)]
    ul = [sum((g[m][n]*uu[n] for n in range(4)), Ser.const(0)) for m in range(4)]
    # nabla_a u_b
    Du = [[ul[b].diff(X[a]) - sum((Gam[c][a][b]*ul[c] for c in range(4)), Ser.const(0)) for b in range(4)] for a in range(4)]
    theta = sum((ginv[a][b]*Du[a][b] for a in range(4) for b in range(4) if not ginv[a][b].iszero()), Ser.const(0))
    # projector h_mu^a = delta + u_mu u^a ; h_mu nu = g + u u
    hmix = [[Ser.const(1 if m == a else 0) + ul[m]*uu[a] for a in range(4)] for m in range(4)]
    hlow = [[g[m][n] + ul[m]*ul[n] for n in range(4)] for m in range(4)]
    # B_mu nu = h_mu^a h_nu^b Du_ab
    tmp = [[sum((hmix[m][a]*Du[a][b] for a in range(4)), Ser.const(0)) for b in range(4)] for m in range(4)]
    B = [[sum((tmp[m][b]*hmix[n][b] for b in range(4)), Ser.const(0)) for n in range(4)] for m in range(4)]
    sig = [[(B[m][n] + B[n][m])*sp.Rational(1, 2) - theta*hlow[m][n]*sp.Rational(1, 3) for n in range(4)] for m in range(4)]
    # acceleration a_mu = u^a nabla_a u_mu ; vorticity (lower) omega_mu nu = h h D_[a u_b]
    acc = [sum((uu[a]*Du[a][m] for a in range(4)), Ser.const(0)) for m in range(4)]
    vort = [[(B[m][n] - B[n][m])*sp.Rational(1, 2) for n in range(4)] for m in range(4)]
    return dict(u0=u0, uu=uu, ul=ul, Du=Du, theta=theta, hmix=hmix, hlow=hlow, sig=sig, acc=acc, vort=vort)

def stress(g, ginv, fl):
    # That^{mu nu} = rho (4/3 u u + g^{-1}/3), rho = (1+dT)^4 (rhohat_bar = 1)
    rho = (1 + dT)*(1 + dT); rho = rho*rho
    uu = fl['uu']
    T = [[rho*(uu[m]*uu[n]*sp.Rational(4, 3) + ginv[m][n]*sp.Rational(1, 3)) for n in range(4)] for m in range(4)]
    return T

def conservation(T, Gam):
    # nabla_mu T^{mu nu}
    out = []
    for nu in range(4):
        s = Ser.const(0)
        for mu in range(4):
            s = s + T[mu][nu].diff(X[mu])
            for lam in range(4):
                s = s + Gam[mu][mu][lam]*T[lam][nu] + Gam[nu][mu][lam]*T[mu][lam]
        out.append(s)
    return out

def solve_rules(cons):
    """Solve nabla T = 0 order by order for d_eta dT_n, d_eta v_n^i. Returns dict of rules (all orders),
    with lower-order rules already substituted into higher-order right-hand sides."""
    rules = {}
    for n in range(1, N+1):
        unknowns = [sp.Derivative(fname('dT', n), eta)] + [sp.Derivative(fname('v', n, i), eta) for i in (1, 2, 3)]
        dummies = sp.symbols('U0:4')
        eqs = []
        for nu in range(4):
            e = sp.expand(cons[nu].c[n].subs(rules))
            e = e.subs(dict(zip(unknowns, dummies)))
            eqs.append(e)
        sol = sp.solve(eqs, dummies, dict=True)
        assert len(sol) == 1, sol
        sol = sol[0]
        for u_, d_ in zip(unknowns, dummies):
            rules[u_] = sp.expand(sol[d_])
        print('  conservation rules at order', n, 'solved;', [len(sp.Add.make_args(rules[u_])) for u_ in unknowns], 'terms')
    return rules

CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fields_cache.pkl')

def build(verbose=True):
    if os.path.exists(CACHE):
        with open(CACHE, 'rb') as f: return pickle.load(f)
    t0 = time.time()
    ginv = inverse_metric(g); print('ginv', time.time()-t0) if verbose else None
    Gam = christoffel(g, ginv); print('Gamma', time.time()-t0) if verbose else None
    fl = fluid(g, ginv, Gam); print('fluid', time.time()-t0) if verbose else None
    T = stress(g, ginv, fl); cons = conservation(T, Gam); print('conservation', time.time()-t0) if verbose else None
    rules = solve_rules(cons); print('rules', time.time()-t0) if verbose else None
    data = dict(g=g, ginv=ginv, Gam=Gam, fl=fl, T=T, cons=cons, rules=rules)
    with open(CACHE, 'wb') as f: pickle.dump(data, f)
    return data

if __name__ == '__main__':
    d = build()
    print('u0:', d['fl']['u0'].c[1])
    print('sigma_11 order 1:', d['fl']['sig'][1][1].c[1])
    print('theta order1:', d['fl']['theta'].c[1])
    for k_, v_ in d['rules'].items():
        if 'dT1' in str(k_) or 'v_1_1' in str(k_): print(k_, '=', v_)
