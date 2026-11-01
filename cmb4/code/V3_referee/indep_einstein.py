"""indep_einstein.py -- V3 referee: independent sympy computation of E_{mu nu} = G_{mu nu} - 8 pi G T_{mu nu}
for g = eta^2 ghat, ghat_00 = -(1+2Psi), ghat_0i = om_i, ghat_ij = (1-2Phi) delta_ij + h_ij, radiation fluid,
a = eta, 8 pi G a^4 rhobar = 3.  Fields X = eps X1 + eps^2 X2/2 + eps^3 X3/6 with generic sympy Functions of (t,x,y,z).
Everything is a plain sympy expression, truncated in eps after every product (no custom polynomial engine).
"""
import sympy as sp, time, pickle, sys, os
from functools import reduce

eps = sp.Symbol('eps')
t, x, y, z = sp.symbols('t x y z')
X = (t, x, y, z)
ORDER = int(os.environ.get('V3_ORDER', '3'))

def trunc(e, n=None):
    n = ORDER if n is None else n
    e = sp.expand(e)
    return sum((e.coeff(eps, m) * eps**m for m in range(n + 1)), sp.Integer(0))

def F(name):
    return sp.Function(name)(*X)

def series_field(name, orders=(1, 2, 3)):
    return sum((eps**n / sp.factorial(n) * F('%s%d' % (name, n)) for n in orders), sp.Integer(0))

Psi = series_field('Psi'); Phi = series_field('Phi'); delta = series_field('delta')
om = [series_field('om%d' % i, (2, 3)) for i in (1, 2, 3)]
v = [series_field('v%d' % i) for i in (1, 2, 3)]
hn = {(0, 0): 'h11', (0, 1): 'h12', (0, 2): 'h13', (1, 1): 'h22', (1, 2): 'h23', (2, 2): 'h33'}
def h(i, j):
    return series_field(hn[(min(i, j), max(i, j))], (2, 3))

gh = sp.zeros(4, 4)
gh[0, 0] = -(1 + 2 * Psi)
for i in range(3):
    gh[0, i + 1] = gh[i + 1, 0] = om[i]
    for j in range(3):
        gh[i + 1, j + 1] = (1 - 2 * Phi if i == j else 0) + h(i, j)
g = t**2 * gh

def build(order=ORDER):
    global ORDER
    ORDER = order
    t0 = time.time()
    eta0 = sp.diag(-1, 1, 1, 1)
    H = gh - eta0
    M = H * eta0                       # H eta0^{-1}
    term = sp.eye(4); ghinv = sp.eye(4)
    for n in range(1, order + 1):
        term = (term * (-M)).applyfunc(trunc)
        ghinv = ghinv + term
    ghinv = (eta0 * ghinv).applyfunc(trunc)
    chk = (gh * ghinv).applyfunc(trunc)
    assert chk == sp.eye(4), "inverse metric wrong"
    ginv = ghinv / t**2
    print("[%.1fs] inverse ok" % (time.time() - t0)); sys.stdout.flush()
    dg = [[[sp.diff(g[a, b], X[c]) for c in range(4)] for b in range(4)] for a in range(4)]
    Gam = [[[None] * 4 for _ in range(4)] for _ in range(4)]
    for a in range(4):
        for b in range(4):
            for c in range(b, 4):
                s = sum((ginv[a, d] * (dg[d][c][b] + dg[d][b][c] - dg[b][c][d]) for d in range(4) if ginv[a, d] != 0), sp.Integer(0))
                Gam[a][b][c] = Gam[a][c][b] = trunc(s / 2)
    print("[%.1fs] Christoffels" % (time.time() - t0)); sys.stdout.flush()
    Ric = sp.zeros(4, 4)
    for b in range(4):
        for d in range(b, 4):
            s = sp.Integer(0)
            for a in range(4):
                s += sp.diff(Gam[a][b][d], X[a]) - sp.diff(Gam[a][b][a], X[d])
                for e in range(4):
                    s += Gam[a][a][e] * Gam[e][b][d] - Gam[a][d][e] * Gam[e][b][a]
            Ric[b, d] = Ric[d, b] = trunc(s)
            print("[%.1fs] Ric %d%d" % (time.time() - t0, b, d)); sys.stdout.flush()
    Rs = trunc(sum((ginv[a, b] * Ric[a, b] for a in range(4) for b in range(4)), sp.Integer(0)))
    G = (Ric - g * Rs / 2).applyfunc(trunc)
    print("[%.1fs] Einstein tensor" % (time.time() - t0)); sys.stdout.flush()
    # fluid
    uu = [sp.Integer(1)] + v                      # uhat^mu / uhat^0
    norm = trunc(-sum((gh[a, b] * uu[a] * uu[b] for a in range(4) for b in range(4)), sp.Integer(0)))  # (uhat^0)^-2
    n_ = norm - 1
    uh0 = trunc(1 - n_ / 2 + sp.Rational(3, 8) * n_**2 - sp.Rational(5, 16) * n_**3)         # (1+n)^(-1/2)
    assert trunc(uh0**2 * norm) == 1
    u_up = [trunc(uh0 * uu[a] / t) for a in range(4)]
    u_dn = [trunc(sum((g[a, b] * u_up[b] for b in range(4)), sp.Integer(0))) for a in range(4)]
    assert trunc(sum((u_dn[a] * u_up[a] for a in range(4)), sp.Integer(0))) == -1
    T8 = sp.zeros(4, 4)
    for a in range(4):
        for b in range(a, 4):
            T8[a, b] = T8[b, a] = trunc(3 / t**4 * (1 + delta) * (sp.Rational(4, 3) * u_dn[a] * u_dn[b] + g[a, b] / 3))
    E = (G - T8).applyfunc(trunc)
    for a in range(4):
        for b in range(4):
            assert E[a, b].coeff(eps, 0) == 0, "background not solved"
    print("[%.1fs] E done, background OK" % (time.time() - t0)); sys.stdout.flush()
    return dict(E=E, G=G, T8=T8, Gam=Gam, ginv=ginv, u_dn=u_dn, u_up=u_up, uh0=uh0)

if __name__ == '__main__':
    R = build(ORDER)
    out = {('E', a, b, n): sp.srepr(R['E'][a, b].coeff(eps, n)) for a in range(4) for b in range(a, 4) for n in range(1, ORDER + 1)}
    out.update({('u_dn', a, n): sp.srepr(R['u_dn'][a].coeff(eps, n)) for a in range(4) for n in range(1, ORDER + 1)})
    out[('uh0',)] = sp.srepr(R['uh0'])
    pickle.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'indep_E.pkl'), 'wb'))
    for a in range(4):
        for b in range(a, 4):
            print("E_%d%d term counts:" % (a, b), [len(sp.Add.make_args(R['E'][a, b].coeff(eps, n))) for n in range(1, ORDER + 1)])
