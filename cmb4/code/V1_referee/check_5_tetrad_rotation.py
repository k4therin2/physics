"""
V1 referee, attack item 6: observer tetrad orientation and C_l.

SPEC sec 3 says a rigid rotation of the sky drops out of C_l exactly, realisation by
realisation.  The rotation is realisation-DEPENDENT (FW transport along a perturbed worldline),
R(eps) = exp(eps A), A = i (r1 Jx + r2 Jy + r3 Jz) in the spin-l representation (anti-Hermitian).
Test for l = 2 with exact arithmetic (sqrt(6) entries of J+-):
  a_lm(eps) = sum_m' D_{m m'}(eps) a'_lm'(eps),  a' = eps a1 + eps^2/2 a2 + eps^3/6 a3.
  1. sum_m |a_lm|^2 equals sum_m |a'_lm|^2 order by order in eps (D unitary).
  2. The eps^4 coefficient 1/4 sum|a2|^2 + 1/3 Re sum a1 a3* is invariant,
     but sum|a2|^2 and Re sum a1 a3* are NOT separately invariant.
Consequence: C^(22)_l and C^(13)_l must be computed in the SAME frame; only the combination
1/4 C^(22) + 1/3 Re C^(13) is tetrad-independent.  Likewise the monopole normalisation (sky
average) is rotation invariant, so T(o) cancels in Theta regardless of frame.
"""
import sympy as sp, random
random.seed(3)
eps = sp.symbols('epsilon', real=True)
l = 2; ms = list(range(-l, l+1)); N = len(ms)
Jz = sp.diag(*ms)
Jp = sp.zeros(N, N); Jm = sp.zeros(N, N)
for i, m in enumerate(ms):
    if m + 1 in ms: Jp[ms.index(m+1), i] = sp.sqrt(l*(l+1) - m*(m+1))
    if m - 1 in ms: Jm[ms.index(m-1), i] = sp.sqrt(l*(l+1) - m*(m-1))
Jx = (Jp + Jm)/2; Jy = (Jp - Jm)/(2*sp.I)
print("[Jx,Jy] = i Jz :", sp.simplify(Jx*Jy - Jy*Jx - sp.I*Jz) == sp.zeros(N))
r = [sp.Rational(random.randint(1, 9), 7) for _ in range(3)]
A = sp.I*(r[0]*Jx + r[1]*Jy + r[2]*Jz)
print("A anti-Hermitian:", sp.simplify(A + A.H) == sp.zeros(N))
D = sum((eps**k*A**k/sp.factorial(k) for k in range(5)), sp.zeros(N))   # exp(eps A) to O(eps^4)

def rc(): return sp.Rational(random.randint(-9, 9), 5) + sp.I*sp.Rational(random.randint(-9, 9), 5)
a1 = sp.Matrix([rc() for _ in ms]); a2 = sp.Matrix([rc() for _ in ms]); a3 = sp.Matrix([rc() for _ in ms])
a_can = eps*a1 + eps**2/2*a2 + eps**3/6*a3
a_rot = (D*a_can).applyfunc(sp.expand)
coef = lambda vec, k: sp.Matrix([sp.expand(c).coeff(eps, k) for c in vec])
b1, b2, b3 = coef(a_rot, 1), 2*coef(a_rot, 2), 6*coef(a_rot, 3)
print("rotated a^(1) == canonical a^(1):", sp.simplify(b1 - a1) == sp.zeros(N, 1))
print("rotated a^(2) == canonical a^(2):", sp.simplify(b2 - a2) == sp.zeros(N, 1), " (expected False: a^(2) picks up 2 A a^(1))")
P22 = lambda v: sp.nsimplify(sp.simplify(sum(sp.expand(c*sp.conjugate(c)) for c in v)))
P13 = lambda v1, v3: sp.nsimplify(sp.simplify(sp.re(sp.expand(sum(v1[i]*sp.conjugate(v3[i]) for i in range(N))))))
print("sum_m |a2|^2           canonical:", P22(a2), "  rotated:", P22(b2))
print("Re sum_m a1 a3*        canonical:", P13(a1, a3), "  rotated:", P13(b1, b3))
tc = P22(a2)/4 + P13(a1, a3)/3; tr = P22(b2)/4 + P13(b1, b3)/3
print("1/4 sum|a2|^2 + 1/3 Re sum a1 a3*  canonical:", tc, "  rotated:", tr, "  invariant:", sp.simplify(tc - tr) == 0)
full = sp.expand(sum(sp.expand(c*sp.conjugate(c)) for c in a_rot) - sum(sp.expand(c*sp.conjugate(c)) for c in a_can))
print("sum_m |a_lm|^2 (rotated - canonical) through eps^4:", sp.simplify(sum(full.coeff(eps, k)*eps**k for k in range(5))))
