"""squeezed_check3.py -- q->0 limit of the third-order amplitude C3 vs the exact dilation family.
With zeta = c + zeta_s, the exact solution is the short-mode solution in rescaled coordinates eta = A eta~, x = B x~:
1 - 2 Phi_new = A^2 B^2 (1 - 2 Phi_s)  =>  the c^2 zeta_s coefficient of Phi(eta->0) is [c^2](A^2B^2) * (-2/3).
Using the angle-averaged k->0 member (a1=-1/3, b1=1, b2=1, a2=1/9 from tadPhi2(0)=-8/9):  A^2B^2 = e^{4c/3}.
Prediction: C3(q->0) = 3! * [c^2 zeta_s] Phi = 6 * (8/9) * (-2/3) = -32/9."""
import sympy as sp
from kernels_api import *
K = load(); C3 = K['third_ic']['C3']
c = sp.Symbol('c')
A = 1 - c / 3 + sp.Rational(1, 9) * c**2 / 2; B = 1 + c + c**2 / 2
coef = sp.series(A**2 * B**2, c, 0, 3).removeO().coeff(c, 2)
pred = 6 * coef * sp.Rational(-2, 3)
lim = sp.limit(C3, q, 0)
print("C3(q->0) =", lim, "   dilation prediction =", pred)
assert sp.simplify(lim - pred) == 0
print("SQUEEZED-LIMIT CHECK OF C3: PASS")
