"""Sensitivity test of the invariance harness: add to the zero-mode ray functionals a spurious but kinematically consistent
displacement (K += beta'(eta) with beta(eta_o) = 0, D += beta) -- e.g. the effect of the documented C2 sign slip
k1 = k1_init + int G1 -- and check that the test then FAILS for a generator with G1 != 0 (translation c2) and for a
constant offset."""
import sympy as sp
from zm_invariance import *
c1 = sp.symbols('c1_1:4'); c2 = sp.symbols('c2_1:4'); lam = sp.symbols('lam1:4')
E, Xs = eta, [X[i+1] + eps*s*(c1[i]*eta + c2[i]*eta**2) for i in range(3)]
out = run('translation + spurious K = lam (eta_o - eta)', E, Xs, verbose=False, Kextra=[lam[i]*(eta_o - eta) for i in range(3)])
verify(out, [], 'spurious-linear')
out = run('translation + spurious K = lam (const)', E, Xs, verbose=False, Kextra=[lam[i] for i in range(3)])
verify(out, [], 'spurious-const')
