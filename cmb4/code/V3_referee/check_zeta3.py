"""check_zeta3.py -- V3: expansion of zeta_nl = 1/4 ln(1+delta) + 1/2 ln(1-2 Phi) to third order (X = eps X1 + eps^2 X2/2 + eps^3 X3/6),
and equivalence with SPEC's ln(1+delta_T) + 1/2 ln(1-2Phi), 1+delta_T = (1+delta)^{1/4}."""
import sympy as sp
eps = sp.Symbol('eps'); d1, d2, d3, P1, P2, P3 = sp.symbols('delta1 delta2 delta3 Phi1 Phi2 Phi3')
dl = eps*d1 + eps**2*d2/2 + eps**3*d3/6; Ph = eps*P1 + eps**2*P2/2 + eps**3*P3/6
z = sp.series(sp.log((1 + dl)**sp.Rational(1, 4)) + sp.log(1 - 2*Ph)/2, eps, 0, 4).removeO()
z2 = sp.expand(2*z.coeff(eps, 2)); z3 = sp.expand(6*z.coeff(eps, 3))
print("zeta^(2) =", z2, " claim ok:", sp.expand(z2 - (d2/4 - d1**2/4 - P2 - 2*P1**2)) == 0)
print("zeta^(3) =", z3, " claim ok:", sp.expand(z3 - (d3/4 - sp.Rational(3, 4)*d1*d2 + d1**3/2 - P3 - 6*P1*P2 - 8*P1**3)) == 0)
