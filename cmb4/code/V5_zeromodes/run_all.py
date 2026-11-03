"""V5 step 2 driver: every residual k = 0 generator, superposed on arbitrary first-order fields (all of Psi1, Phi1, dT1,
v1, omega1, h1 generic functions), map compared through second order.  Log: run_all.log"""
import sympy as sp, sys, time
from zm_invariance import *
L1ray = sp.sympify(sum(e for w, e in LG[1] if w == 'ray'))
results = {}
def gen(name, E, Xs, naive=False):
    t1 = time.time()
    out = run(name, E, Xs, verbose=True, naive=naive)
    alpha = sp.expand(sp.expand(E - eta).coeff(eps, 1).coeff(s, 1))
    W_list = [alpha*L1ray] if alpha != 0 else []
    ok, sol = verify(out, W_list, name)
    results[name] = ok
    print('  (%.0fs)' % (time.time() - t1), flush=True)
    return out, sol
a, b = sp.symbols('a b')
C11, C12, C13, C22, C23 = sp.symbols('C11 C12 C13 C22 C23'); C = sp.Matrix([[C11, C12, C13], [C12, C22, C23], [C13, C23, -C11 - C22]])
S11, S12, S13, S22, S23 = sp.symbols('S11 S12 S13 S22 S23'); Sm = sp.Matrix([[S11, S12, S13], [S12, S22, S23], [S13, S23, -S11 - S22]])
xv = sp.Matrix([x, y, z])
which = sys.argv[1:] or ['time', 'dil', 'shear', 'rot', 'transl', 'tidal', 'transl0', 'naive_h']
if 'time' in which: gen('time rescaling', (1 + eps*s*a)*eta, [x, y, z])
if 'dil' in which: gen('dilation', eta, [(1 + eps*s*b)*xi for xi in (x, y, z)])
if 'shear' in which:
    out, sol = gen('constant anisotropic shear x_old = (1 - S/2) x', eta, list(xv - eps*s*Sm*xv/2))
    print('  expected rotation from Gram-Schmidt: A_ab = S_ab/2 (a<b) -> omega = (A_23, -A_13, A_12)-type; fitted above')
if 'rot' in which:
    w0 = sp.symbols('w0_1:4'); w1 = sp.symbols('w1_1:4')
    w = sp.Matrix([w0[i] + w1[i]*eta for i in range(3)])
    gen('time-dependent rotation x_old = x + w(eta) x x', eta, list(xv + eps*s*w.cross(xv)))
if 'transl' in which:
    c1 = sp.symbols('c1_1:4'); c2 = sp.symbols('c2_1:4')
    gen('time-dependent translation c(eta) = c1 eta + c2 eta^2 (c(0)=0: label kept)', eta, [X[i+1] + eps*s*(c1[i]*eta + c2[i]*eta**2) for i in range(3)])
if 'tidal' in which:
    f1, f2, g0, g1, g2, g3 = sp.symbols('f1 f2 g0 g1 g2 g3')
    fq = f1*eta + f2*eta**2; gq = g0 + g1*eta + g2*eta**2 + g3*eta**3
    gen('quadratic (tidal) mode: eta_old = eta + f(eta) x.C.x, x_old = x + g(eta) C x', eta + eps*s*fq*(xv.T*C*xv)[0], list(xv + eps*s*gq*C*xv))
if 'transl0' in which:
    c0 = sp.symbols('c0_1:4')
    gen('NEGATIVE CONTROL: constant translation (changes the observer label)', eta, [X[i+1] + eps*s*c0[i] for i in range(3)])
if 'treparam' in which:
    a1, a2, a3 = sp.symbols('a1 a2 a3')
    gen('general time reparametrization eta_old = eta + F(eta), F = a1 eta + a2 eta^2 + a3 eta^3', eta + eps*s*(a1*eta + a2*eta**2 + a3*eta**3), [x, y, z])
if 'grad' in which:
    b1 = sp.symbols('b1:4'); f1, f2, g1, g2, g3 = sp.symbols('f1 f2 g1 g2 g3')
    bv = sp.Matrix(b1); fq = f1*eta + f2*eta**2; gq = g1*eta + g2*eta**2 + g3*eta**3
    gen('dipolar (gradient) mode: eta_old = eta + f(eta) b.x, x_old = x + g(eta) b', eta + eps*s*fq*(bv.T*xv)[0], list(xv + eps*s*gq*bv))
if 'naive_h' in which:
    gen('NEGATIVE CONTROL: constant h1 = -S superposed WITHOUT transforming the first-order fields', eta, list(xv - eps*s*Sm*xv/2), naive=True)
print('\nSUMMARY:', results)
