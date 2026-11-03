"""V5 gauge engine: exact coordinate transformation of the perturbed conformal-frame solution (C2/SPEC conventions).
Old coordinates as functions of the new ones:  eta_old = E(eta,x),  x_old^i = Xi(eta,x)   (sympy expressions, polynomial
in eps; the zero-mode amplitude is eps*s).  Physical metric g = a^2 ghat with a = eta (in EACH coordinate system):
    g_new_{mu nu}(x) = J^a_mu J^b_nu E^2 ghat_old_ab(x_old),  ghat_new = g_new/eta^2,  J^a_mu = d x_old^a / d x_new^mu
Fluid: T_phys = That/a is a scalar  =>  1 + dT_new = (eta/E)(1 + dT_old(x_old)) ;  u is a vector => v_new^i = (J^-1 u)^i/(J^-1 u)^0.
Returns the new Poisson-type fields as eps-series (coefficients c_n = X_n/n!) in terms of OLD field functions evaluated at
the new coordinate point, explicit coordinates and the generator parameters."""
import sympy as sp
from ser import Ser, N, eps
import fields
from fields import X, eta, x, y, z, F, build
_d = build(verbose=False)

def ser_of(expr):
    """eps-polynomial expression -> Ser (truncated at N)."""
    e = sp.expand(expr)
    return Ser([e.coeff(eps, k) if k > 0 else e.subs(eps, 0) for k in range(N+1)])

def shift(Q, Delta):
    """Q(x + Delta), Delta = list of 4 Ser starting at O(eps)."""
    out = Q; term = Q
    for m in range(1, N+1):
        term = sum((Delta[a]*term.diff(X[a]) for a in range(4)), Ser.const(0))*sp.Rational(1, m)
        if term.iszero(): break
        out = out + term
    return out

def transform(E, Xs):
    old = [ser_of(E)] + [ser_of(Xi) for Xi in Xs]
    Delta = [old[a] - Ser.const(X[a]) for a in range(4)]
    assert all(D.c[0] == 0 for D in Delta), 'zero-mode map must be identity at eps^0'
    J = [[old[a].diff(X[m]) for m in range(4)] for a in range(4)]          # J[a][mu]
    fac = (old[0]*(1/Ser.const(eta)))*(old[0]*(1/Ser.const(eta)))           # (E/eta)^2
    gold = [[shift(fields.g[a][b], Delta) for b in range(4)] for a in range(4)]
    gnew = [[None]*4 for _ in range(4)]
    for m in range(4):
        for q in range(m, 4):
            s_ = Ser.const(0)
            for a in range(4):
                for b in range(4):
                    if gold[a][b].iszero(): continue
                    s_ = s_ + J[a][m]*J[b][q]*gold[a][b]
            gnew[m][q] = gnew[q][m] = s_*fac
    # inverse Jacobian (J = 1 + O(eps)) by Neumann series
    M = [[(J[a][b] - Ser.const(1 if a == b else 0)) for b in range(4)] for a in range(4)]
    Jinv = [[Ser.const(1 if a == b else 0) for b in range(4)] for a in range(4)]
    P = [[Ser.const(1 if a == b else 0) for b in range(4)] for a in range(4)]
    for k in range(1, N+1):
        P = [[sum((P[a][c]*M[c][b] for c in range(4)), Ser.const(0)) for b in range(4)] for a in range(4)]
        Jinv = [[Jinv[a][b] + P[a][b]*((-1)**k) for b in range(4)] for a in range(4)]
    uold = [shift(_d['fl']['uu'][a], Delta) for a in range(4)]
    unew = [sum((Jinv[m][a]*uold[a] for a in range(4)), Ser.const(0)) for m in range(4)]
    vnew = [unew[i+1]/unew[0] for i in range(3)]
    dTnew = (Ser.const(eta)/old[0])*(1 + shift(fields.dT, Delta)) - 1
    Psi = (gnew[0][0] + 1)*sp.Rational(-1, 2)
    om = [gnew[0][i+1] for i in range(3)]
    Phi = (gnew[1][1] + gnew[2][2] + gnew[3][3] - 3)*sp.Rational(-1, 6)
    h = [[gnew[i+1][j+1] - (1 - 2*Phi)*(1 if i == j else 0) for j in range(3)] for i in range(3)]
    return dict(Psi=Psi, Phi=Phi, dT=dTnew, v=vnew, om=om, h=h, J=J, Jinv=Jinv, old=old, Delta=Delta)

def substitution(new, s=None):
    """dict: old field Function application -> new field expression (coefficient X_n = n! c_n), truncated O(s^2) if s given."""
    fac = {1: 1, 2: 2}
    def tr(e):
        e = sp.expand(e)
        if s is None: return e
        return sp.Add(*[t for t in sp.Add.make_args(e) if sp.degree(t, s) <= 1])
    sub = {}
    for n_ in (1, 2):
        sub[F('Psi%d' % n_)] = tr(fac[n_]*new['Psi'].c[n_])
        sub[F('Phi%d' % n_)] = tr(fac[n_]*new['Phi'].c[n_])
        sub[F('dT%d' % n_)] = tr(fac[n_]*new['dT'].c[n_])
        for i in range(3):
            sub[F('v_%d_%d' % (i+1, n_))] = tr(fac[n_]*new['v'][i].c[n_])
            sub[F('om_%d_%d' % (i+1, n_))] = tr(fac[n_]*new['om'][i].c[n_])
            for j in range(i, 3):
                if (i, j) == (2, 2): continue
                sub[F('h_%d%d_%d' % (i+1, j+1, n_))] = tr(fac[n_]*new['h'][i][j].c[n_])
    return sub
