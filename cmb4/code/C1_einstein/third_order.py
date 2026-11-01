"""third_order.py -- contracted (13-type) third-order sources and equations.
<X3_k zeta_k'> = (2pi)^3 delta(k+k') P(k) X3c(k,eta),  X3c = int d^3q/(2pi)^3 P(q) K3c_X(k,q;eta).
For every component of E^{(3)} = 6 [eps^3] E we build K3c (per P(q), NOT symmetrised under q->-q; only |k-q| =: p appears):
  cubic  A B C      ->  A(q)B(-q)C(k) + A(q)C(-q)B(k) + B(q)C(-q)A(k)
  bilinear A1 B2    ->  2 A1(q) B2[(-q,k); total k-q]  +  A1(k) <B2>(eta)   (tadpole = angle-averaged k->0 limit)
Second-order fields: Q2, delta2, V2 (longitudinal), v2T, om2 explicit; Phi2 = KPhi2(a,b,c,eta) and h2 = Kh(a,b,c,eta) abstract,
(a,b,c) = (|q1|,|q2|,|q1+q2|); tadpoles tadPhi2(eta), tadQ2(eta), taddelta2(eta) abstract (explicit forms in third_order_ic.py).
"""
import sympy as sp, pickle, time, itertools
from tpoly import *
from fourier import *
t0 = time.time()
S2 = {kk: sp.parse_expr(vv) for kk, vv in pickle.load(open('second_order.pkl', 'rb')).items()}
d = load_polys('E.pkl'); E = d['E']
TP = lambda kk: TF['Phi'](kk, eta)
a_, b_, c_ = sp.symbols('a b c', positive=True)
KPhi2 = sp.Function('KPhi2'); Kh = sp.Function('Kh')
tadPhi2, tadQ2, taddelta2, tadS00 = [sp.Function(n)(eta) for n in ('tadPhi2', 'tadQ2', 'taddelta2', 'tadS00')]

def at(expr, A, B, C):
    """evaluate an invariant second-order kernel (symbols q,p,k = |q1|,|q2|,|q1+q2|) at (A,B,C)."""
    return expr.subs({q: a_, p: b_, k: c_}, simultaneous=True).subs({a_: A, b_: B, c_: C}, simultaneous=True)

# ---- second-order fields for the pair (q1,q2) = (-q, k), total momentum k2 = k - q (vector PV, magnitude p)
A, B, C = q, k, p
KQ_pair = at(S2['K_Q'], A, B, C)
KS00_pair = at(S2['K_S00'], A, B, C)
KS0z_pair = at(S2['K_S0z'], A, B, C)          # longitudinal (k2-hat . S_0) component
KPhi2_pair = KPhi2(A, B, C, eta)
KPsi2_pair = KPhi2_pair - KQ_pair
Kdelta2_pair = eta**2 / 3 * KS00_pair - sp.Rational(2, 3) * eta**2 * C**2 * KPhi2_pair - 2 * eta * sp.diff(KPhi2_pair, eta) - 2 * KPsi2_pair
KV2_pair = -(eta**2 / 4) * (2 * (sp.diff(KPhi2_pair, eta) + KPsi2_pair / eta) + KS0z_pair / (I * C))   # v2_i = i k2_i V2 + transverse
Kom_pair = at(S2['K_om'], A, B, C); KvT_pair = at(S2['K_vT'], A, B, C); Kh_pair = Kh(A, B, C, eta)
k2v = sp.Matrix(PV.vec); q1v = -sp.Matrix(QV.vec)        # q1 = -q
Pm = sp.eye(3) - k2v * k2v.T / p**2
Lam = lambda i, j, kk, l: Pm[i, kk] * Pm[j, l] - sp.Rational(1, 2) * Pm[i, j] * Pm[kk, l]
om2_vec = [sum(Pm[i, j] * q1v[j] for j in range(3)) * Kom_pair for i in range(3)]
vT_vec = [sum(Pm[i, j] * q1v[j] for j in range(3)) * KvT_pair for i in range(3)]
v2_vec = [I * k2v[i] * KV2_pair + vT_vec[i] for i in range(3)]
h2_mat = [[sum(Lam(i, j, kk, l) * q1v[kk] * q1v[l] for kk in range(3) for l in range(3)) * Kh_pair for j in range(3)] for i in range(3)]
hidx = {'h11': (0, 0), 'h12': (0, 1), 'h13': (0, 2), 'h22': (1, 1), 'h23': (1, 2), 'h33': (2, 2)}
second_pair = {'Phi': KPhi2_pair, 'Psi': KPsi2_pair, 'delta': Kdelta2_pair}
for i in range(3):
    second_pair['v%d' % (i + 1)] = v2_vec[i]; second_pair['om%d' % (i + 1)] = om2_vec[i]
for nm, (i, j) in hidx.items():
    second_pair[nm] = h2_mat[i][j]
second_tad = {'Phi': tadPhi2, 'Psi': tadPhi2 - tadQ2, 'delta': taddelta2}

def jet_second(field, m, mom=PV):
    base = second_pair[field]
    if m[0]:
        base = sp.diff(base, eta, m[0])
    for i in range(3):
        if m[i + 1]:
            base = base * (I * mom.vec[i]) ** m[i + 1]
    return base
def jet_second_tad(field, m):
    if any(m[1:]) or field not in second_tad:
        return sp.Integer(0)
    return sp.diff(second_tad[field], eta, m[0]) if m[0] else second_tad[field]

def contracted_kernel(P):
    tot = 0
    for (o, e, mono), cf in P.d.items():
        if len(mono) < 2:
            continue
        cfs = sp.Rational(cf.numerator, cf.denominator) * eta**e
        orders = [REG.order(j) for j in mono]
        if len(mono) == 2:
            j1, j2 = mono if orders[0] == 1 else mono[::-1]
            f1, _, m1 = REG.info[j1]; f2, _, m2 = REG.info[j2]
            assert REG.order(j1) == 1 and REG.order(j2) == 2
            tot += cfs * (2 * jet_first(f1, m1, QV) * jet_second(f2, m2) + jet_first(f1, m1, KV) * jet_second_tad(f2, m2))
        else:
            assert orders == [1, 1, 1]
            for i in range(3):
                others = [mono[j] for j in range(3) if j != i]
                fi, _, mi = REG.info[mono[i]]; fa, _, ma = REG.info[others[0]]; fb, _, mb = REG.info[others[1]]
                tot += cfs * jet_first(fi, mi, KV) * jet_first(fa, ma, QV) * jet_first(fb, mb, QV.neg())
    return tot

# first-order on-shell relations
def on_shell(expr):
    expr = expr.doit()
    for kk in (q, p, k):
        TPk = TP(kk)
        expr = expr.subs(TF['Psi'](kk, eta), TPk)
        expr = expr.subs(TF['delta'](kk, eta), -sp.Rational(2, 3) * eta**2 * kk**2 * TPk - 2 * eta * sp.Derivative(TPk, eta) - 2 * TPk)
        expr = expr.subs(TF['V'](kk, eta), -(eta**2 / 2) * (sp.Derivative(TPk, eta) + TPk / eta)).doit()
        for nn in (4, 3, 2):
            d2 = -4 / eta * sp.Derivative(TPk, eta) - kk**2 / 3 * TPk
            expr = expr.subs(sp.Derivative(TPk, (eta, nn)), sp.diff(d2, eta, nn - 2)).doit()
        expr = expr.subs(sp.Derivative(TPk, (eta, 2)), -4 / eta * sp.Derivative(TPk, eta) - kk**2 / 3 * TPk).doit()
    return sp.expand(expr)

N = 3
lin, ker = {}, {}
comps = [(0, 0), (0, 1), (0, 2), (0, 3), (1, 1), (2, 2), (3, 3), (1, 2), (1, 3), (2, 3)]
for (a, b) in comps:
    P = E[a][b].order_part(N) * 6
    lin[(a, b)] = linear_part(P, N)
    ker[(a, b)] = on_shell(contracted_kernel(P))
    print("component", a, b, "linear:", lin[(a, b)], "| contracted source terms:", len(sp.Add.make_args(ker[(a, b)])), "(%.0fs)" % (time.time() - t0))

Phi3, Psi3, de3 = amp('Phi', 3), amp('Psi', 3), amp('delta', 3)
Q3 = sp.Function('Q3')(eta)
# consistency: linear operators identical to lower orders
assert sp.simplify(sp.expand(lin[(3, 3)] - (lin[(1, 1)] + lin[(2, 2)]) / 2) + k**2 * (Phi3 - Psi3)) == 0
odeL = sp.expand(sp.expand((lin[(1, 1)] + lin[(2, 2)] + lin[(3, 3)] - lin[(0, 0)]) / 6).subs(Psi3, Phi3 - Q3).doit())
assert sp.simplify(odeL - (sp.Derivative(Phi3, (eta, 2)) + 4 / eta * sp.Derivative(Phi3, eta) + k**2 / 3 * Phi3 - sp.Derivative(Q3, eta) / eta + k**2 / 3 * Q3)) == 0
print("linear operators at third order identical to second order: OK")
K3 = {}
K3['Q'] = sp.expand((ker[(3, 3)] - (ker[(1, 1)] + ker[(2, 2)]) / 2) / k**2)                      # Q3c = K3['Q']
K3['S00'] = ker[(0, 0)]; K3['S0z'] = ker[(0, 3)]; K3['S0x'] = ker[(0, 1)]; K3['S0y'] = ker[(0, 2)]
K3['Sxz'] = ker[(1, 3)]; K3['Syz'] = ker[(2, 3)]; K3['Sxy'] = ker[(1, 2)]; K3['Sxxyy'] = sp.expand(ker[(1, 1)] - ker[(2, 2)])
K3['SPhi'] = sp.expand(-(ker[(1, 1)] + ker[(2, 2)] + ker[(3, 3)] - ker[(0, 0)]) / 6 + sp.diff(K3['Q'], eta) / eta - k**2 / 3 * K3['Q'])
print("S0y, Syz, Sxy contracted sources vanish identically:", K3['S0y'] == 0, K3['Syz'] == 0, K3['Sxy'] == 0)
# delta3c, V3c in terms of Phi3c, Q3c:
K3['delta3_expr'] = eta**2 / 3 * sp.Symbol('S00c') - sp.Rational(2, 3) * eta**2 * k**2 * Phi3 - 2 * eta * sp.Derivative(Phi3, eta) - 2 * (Phi3 - Q3)
K3['V3_expr'] = -(eta**2 / 4) * (2 * (sp.Derivative(Phi3, eta) + (Phi3 - Q3) / eta) + sp.Symbol('S0zc') / (I * k))
pickle.dump({kk: sp.srepr(vv) for kk, vv in K3.items()}, open('third_order.pkl', 'wb'))
with open('third_order.tex', 'w') as fh:
    fh.write("% contracted third-order kernels K3c(k,q,mu;eta) per P(q); p=|k-q|; KPhi2(a,b,c,eta), Kh(a,b,c,eta): second-order Phi2/h2 kernels at (|q1|,|q2|,|q1+q2|)\n")
    fh.write("% tadPhi2, tadQ2, taddelta2: angle-averaged k->0 tadpoles (third_order_ic.py)\n")
    for kk, vv in K3.items():
        fh.write("\\paragraph{%s} (%d terms)\n\\begin{dmath}\n%s\n\\end{dmath}\n" % (kk.replace('_', '\\_'), len(sp.Add.make_args(vv)), sp.latex(vv)))
print("saved third_order.pkl/.tex; term counts:", {kk: len(sp.Add.make_args(vv)) for kk, vv in K3.items()}, "(%.0fs)" % (time.time() - t0))
