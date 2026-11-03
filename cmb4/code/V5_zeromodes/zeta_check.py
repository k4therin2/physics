"""V5: which residual k = 0 generators preserve the initial-data statement zeta_nl = eps zeta_1 (SPEC §5), h(eta->0) = 0 and the
observer label?  zeta_nl = lim_{eta->0} [ln(1+dT) + (1/2) ln(1-2Phi)] evaluated with the transformed fields (gauge.py), linear in the
generator amplitude s, through eps^2 (generator counted as first order; cross terms s x first-order fields)."""
import sympy as sp
from ser import Ser, eps
from fields import X, eta, x, y, z, F
import gauge
s = sp.Symbol('s')
xv = sp.Matrix([x, y, z])
def zeta_change(E, Xs):
    new = gauge.transform(E, Xs)
    zeta = (1 + new['dT']).log() + (1 - 2*new['Phi']).log()*sp.Rational(1, 2)
    out = []
    for o in (1, 2):
        e = sp.expand(zeta.c[o]).coeff(s, 1)
        e = sp.expand(e.doit()).subs(eta, 0)
        out.append(sp.simplify(e))
    h0 = [sp.simplify(sp.expand(new['h'][i][j].c[1]).coeff(s, 1).subs(eta, 0)) for i in range(3) for j in range(i, 3)]
    lab = [sp.simplify(sp.expand(Xs[i] - X[i+1]).coeff(eps, 1).coeff(s, 1).subs({eta: 0, x: 0, y: 0, z: 0})) for i in range(3)]
    return out, h0, lab
a1, a2, a3, b = sp.symbols('a1 a2 a3 b')
S11, S12, S13, S22, S23 = sp.symbols('S11 S12 S13 S22 S23'); Sm = sp.Matrix([[S11, S12, S13], [S12, S22, S23], [S13, S23, -S11 - S22]])
w0 = sp.symbols('w0_1:4'); w1 = sp.symbols('w1_1:4'); c1 = sp.symbols('c1_1:4'); c0 = sp.symbols('c0_1:4')
C11, C12, C13, C22, C23 = sp.symbols('C11 C12 C13 C22 C23'); C = sp.Matrix([[C11, C12, C13], [C12, C22, C23], [C13, C23, -C11 - C22]])
f1, g0, g1 = sp.symbols('f1 g0 g1')
gens = {
 'time reparametrisation F = a1 eta + a2 eta^2 + a3 eta^3': (eta + eps*s*(a1*eta + a2*eta**2 + a3*eta**3), [x, y, z]),
 'dilation b': (eta, [(1 + eps*s*b)*xi for xi in (x, y, z)]),
 'constant shear S': (eta, list(xv - eps*s*Sm*xv/2)),
 'rotation w(eta) = w0 + w1 eta': (eta, list(xv + eps*s*sp.Matrix([w0[i] + w1[i]*eta for i in range(3)]).cross(xv))),
 'translation c(eta) = c1 eta': (eta, [X[i+1] + eps*s*c1[i]*eta for i in range(3)]),
 'translation c0 (constant)': (eta, [X[i+1] + eps*s*c0[i] for i in range(3)]),
 'tidal f = f1 eta, g = g0 + g1 eta': (eta + eps*s*f1*eta*(xv.T*C*xv)[0], list(xv + eps*s*(g0 + g1*eta)*C*xv)),
}
for name, (E, Xs) in gens.items():
    (z1, z2), h0, lab = zeta_change(E, Xs)
    print('%-58s d zeta^(1) = %s ;  d zeta^(2)/2 = %s' % (name, z1, z2))
    print('%-58s h1(eta->0) zm = %s ; label shift = %s' % ('', h0, lab))
