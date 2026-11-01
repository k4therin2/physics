"""check_wick.py -- V3: brute-force Wick contraction on a toy cubic/bilinear expression to test C1's contraction rules
  cubic A1 B1 C1  -> A(q)B(-q)C(k) + A(q)C(-q)B(k) + B(q)C(-q)A(k)
  bilinear A1 B2  -> 2 A(q) B2[(-q,k)] + A(k) <B2>
Discrete Gaussian field on a momentum lattice; kernels are random functions of momenta (A(q) depends on vector q)."""
import itertools, random, sympy as sp
random.seed(1)
# momentum lattice: 1D integers -N..N, zeta_k with zeta_{-k} = conj (real field); Gaussian with <zeta_k zeta_k'> = delta_{k+k'} P(|k|)
N = 3
moms = list(range(-N, N + 1))
P = {kk: sp.Rational(random.randint(1, 9), random.randint(1, 9)) for kk in range(0, N + 1)}
P = {kk: P[abs(kk)] for kk in moms}
def Av(name, kk):   # random "transfer function" depending on the signed momentum
    return sp.Symbol('%s_%s' % (name, str(kk).replace('-', 'm')))
# Wick: <zeta_{k1} zeta_{k2} ... zeta_{kn}> = sum over pairings prod delta(k_i + k_j) P(k_i)
def wick(ks):
    if not ks: return sp.Integer(1)
    k0 = ks[0]; tot = 0
    for j in range(1, len(ks)):
        if ks[0] + ks[j] == 0:
            rest = ks[1:j] + ks[j + 1:]
            tot += P[k0] * wick(rest)
    return tot
# fields: A1(x) = sum_q A(q) zeta_q e^{iqx};  product A1 B1 C1 at momentum k: sum_{q1+q2+q3=k} A(q1)B(q2)C(q3) zeta zeta zeta
k = 2
def cubic_corr(k):
    tot = 0
    for q1, q2 in itertools.product(moms, moms):
        q3 = k - q1 - q2
        if q3 not in moms: continue
        tot += Av('A', q1) * Av('B', q2) * Av('C', q3) * wick([q1, q2, q3, -k])
    return sp.expand(tot)
# C1 rule, per P(q): sum_q P(q) [A(q)B(-q)C(k) + A(q)C(-q)B(k) + B(q)C(-q)A(k)] times P(k) (the zeta_k' pairing)
rule = sp.expand(P[k] * sum(P[q] * (Av('A', q) * Av('B', -q) * Av('C', k) + Av('A', q) * Av('C', -q) * Av('B', k) + Av('B', q) * Av('C', -q) * Av('A', k)) for q in moms))
print("cubic: brute-force Wick == C1 rule:", sp.expand(cubic_corr(k) - rule) == 0)
# bilinear: B2(x) = sum_{q1,q2} KB(q1,q2) zeta_q1 zeta_q2 e^{i(q1+q2)x}, KB symmetric
def KB(q1, q2):
    a, b = sorted([q1, q2]); return sp.Symbol('KB_%s_%s' % (str(a).replace('-', 'm'), str(b).replace('-', 'm')))
def bilinear_corr(k):
    tot = 0
    for q1, q2 in itertools.product(moms, moms):
        q3 = k - q1 - q2
        if q3 not in moms: continue
        tot += Av('A', q1) * KB(q2, q3) * wick([q1, q2, q3, -k])
    return sp.expand(tot)
rule_b = sp.expand(P[k] * sum(P[q] * (2 * Av('A', q) * KB(-q, k) + Av('A', k) * KB(q, -q)) for q in moms))
print("bilinear: brute-force Wick == C1 rule (2 A(q) B2(-q,k) + A(k) <B2>):", sp.expand(bilinear_corr(k) - rule_b) == 0)
# note: the Wick pairing with k1 = -k uses q3 = k - q1 - q2 = k paired with zeta_{-k}; disconnected pieces with k=0 excluded since k != 0
