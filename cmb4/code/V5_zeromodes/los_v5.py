"""V5: C2's line-of-sight assembly (C2_observable/los.py, shear form only) rebuilt at N = 2 with FIRST-ORDER
omega_i and h_ij allowed (needed to superpose a constant shear h1 = -S or a vector zero mode on first-order fields).
Physics/algorithm identical to C2 (positional expansion, endpoint shifts, observer/emission functionals); the
direct form and display splits are dropped.  Output: L_v5.pkl with L[m] = list of (where, expr) for the coefficient of
eps^m in Lambda = -S (m = 1, 2), plus the emission shifts and symbols."""
import sympy as sp, pickle, time
from ser import Ser, N, eps
from fields import build, X, eta, x, y, z, F, fname, dT
import ray
from ray import n, kfun, Dxfun, G1, G2half, k1_init, k2_init_half, deta_o1, deta_o2, dxo1, dxo2
t0 = time.time()
d = build(verbose=False)
sig, ul = d['fl']['sig'], d['fl']['ul']
kser = ray.kser
eta_s = sp.Symbol('eta_star'); eta_o = sp.Symbol('eta_o_bar')
deta_e1, deta_e2 = sp.symbols('deta_e1 deta_e2')
def contract2(Tm):
    return sum((Tm[m][q]*kser[m]*kser[q] for m in range(4) for q in range(4) if not Tm[m][q].iszero()), Ser.const(0))
muk = -sum((ul[m]*kser[m] for m in range(4)), Ser.const(0))
I_sigma = contract2(sig)/muk
Dx = [Ser([0, Dxfun[1][i], Dxfun[2][i]/2]) for i in range(3)]
def pos_expand(Q):
    out = Q; term = Q
    for m in range(1, N+1):
        term = sum((Dx[i]*term.diff(X[i+1]) for i in range(3)), Ser.const(0))*sp.Rational(1, m)
        if term.iszero(): break
        out = out + term
    return out
G2tot_half = [G2half[i] + sum(Dxfun[1][j]*sp.diff(G1[i], X[j+1]) for j in range(3)) for i in range(3)]
Drep = {}
for i in range(3):
    Drep[sp.Derivative(kfun[1][i], eta)] = G1[i]
    Drep[sp.Derivative(kfun[2][i], eta)] = 2*G2tot_half[i]
    Drep[sp.Derivative(Dxfun[1][i], eta)] = kfun[1][i]
    Drep[sp.Derivative(Dxfun[2][i], eta)] = kfun[2][i]
def Dtot_expr(e):
    r = sp.diff(e, eta) - sum(n[i]*sp.diff(e, X[i+1]) for i in range(3))
    return sp.expand(r.xreplace(Drep))
def Dtot(Q, upto=N): return Ser([Dtot_expr(t) if i <= upto else sp.Integer(0) for i, t in enumerate(Q.c)])
def trunc(Q, upto): return Ser([t if i <= upto else sp.Integer(0) for i, t in enumerate(Q.c)])
def shift_eta(Q, dser):
    out = Q; term = Q; p = Ser.const(1)
    for m in range(1, N+1):
        term = Dtot(term, upto=N-m); p = p*dser
        add = term*p*sp.Rational(1, sp.factorial(m))
        if add.iszero(): break
        out = out + add
    return out
sub_obs_fun = {}
for i in range(3):
    sub_obs_fun[kfun[1][i]] = k1_init[i]
    sub_obs_fun[Dxfun[1][i]] = dxo1[i] + n[i]*deta_o1
sub_emi_fun = {}
for i in range(3):
    for f_ in (kfun[1][i], kfun[2][i], Dxfun[1][i], Dxfun[2][i]):
        sub_emi_fun[f_] = f_.subs(eta, eta_s)
def at_obs(e): return sp.expand(sp.expand(e).xreplace(sub_obs_fun))
def at_emi(e): return sp.expand(sp.expand(e).xreplace(sub_emi_fun))
RT = pos_expand(dT)
dser_e = Ser([0, deta_e1, deta_e2/2])
lhs = shift_eta(RT, dser_e); rhs = dser_e*(1/eta_s)
emission = {}
emission[deta_e1] = sp.solve(sp.expand(at_emi(lhs.c[1] - rhs.c[1])), deta_e1)[0]
emission[deta_e2] = sp.expand(sp.solve(sp.expand(at_emi(lhs.c[2] - rhs.c[2]).subs(deta_e1, emission[deta_e1])), deta_e2)[0])
dser_o = Ser([0, deta_o1, deta_o2/2])
def assemble(I):
    Fq = pos_expand(I)
    pieces = {1: [], 2: []}
    for m in (1, 2): pieces[m].append(('ray', Fq.c[m]))
    for sign, dser, where, ev in ((1, dser_o, 'obs', at_obs), (-1, dser_e, 'emi', at_emi)):
        term = Fq; p = Ser.const(1)
        for j in range(1, N+1):
            p = p*dser
            term = Dtot(term, upto=N-j) if j > 1 else trunc(term, N-j)
            add = term*p*sp.Rational(1, sp.factorial(j))
            for m in (1, 2):
                if add.c[m] != 0: pieces[m].append((where, sign*ev(add.c[m])))
    return pieces
S = assemble(I_sigma)
L = {m: [(w, sp.expand(-e)) for w, e in S[m]] for m in (1, 2)}
for m in (1, 2): print('L coefficient eps^%d:' % m, [(w, len(sp.Add.make_args(e))) for w, e in L[m]])
with open('L_v5.pkl', 'wb') as f:
    pickle.dump(dict(L=L, emission=emission, I_sigma=I_sigma.c,
                     symbols=dict(n=n, eta_star=eta_s, eta_o_bar=eta_o, deta_o=(deta_o1, deta_o2), dxo=(dxo1, dxo2),
                                  deta_e=(deta_e1, deta_e2), kfun=kfun, Dxfun=Dxfun)), f)
print('saved L_v5.pkl t=%.1f' % (time.time()-t0))
