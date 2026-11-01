"""build_phi2.py -- closed form of the particular integral  P(k,q,p;eta) = int_0^eta G_Phi(eta,eta') K_SPhi(eta') deta'
   written as  P = u2(eta) A(eta) - u1(eta) B(eta),  A = int_0^eta u1(eta') S(eta')/W(eta'),  B = int_0^eta u2 S/W,
   W(eta') = -3 sqrt3/(eta'^4 k^3),  u1 = (sin x - x cos x)/x^3, u2 = (cos x + x sin x)/x^3, x = k eta/sqrt3.
   Full regular solution:  K_Phi2 = P + K_Phi2_0(k,q,p) ghat(k eta),  ghat = 3 u1."""
import sympy as sp, pickle, sys, time
sys.path.insert(0, '.'); sys.path.insert(0, '../C1_einstein')
from trig_integrate import *
from kernels_api import load, closed_form, eta, k, q, p
K = load(); S = K['second']
etap = sp.Symbol('etap', positive=True)
t0 = time.time()
SPhi = closed_form(S['K_SPhi']).subs(eta, etap)
u1, u2 = S['u1'], S['u2']
W = -3 * sp.sqrt(3) / (eta**4 * k**3)
G_check = sp.simplify((u2 * u1.subs(eta, etap) - u1 * u2.subs(eta, etap)) / W.subs(eta, etap) - S['G_Phi'])
assert G_check == 0
A_int = sp.expand((u1 / W).subs(eta, etap) * SPhi)
B_int = sp.expand((u2 / W).subs(eta, etap) * SPhi)
print("integrands built (%.1fs); nterms A_int=%d" % (time.time() - t0, len(sp.Add.make_args(A_int))))
A, termsA = integrate_0_to_T(A_int, etap, eta)
print("A done (%.1fs): %d (n,w) terms; powers %s" % (time.time() - t0, len(termsA), sorted(set(n for n, w in termsA))))
print("frequencies:", sorted(set(w for (n, w) in termsA), key=str))
B, termsB = integrate_0_to_T(B_int, etap, eta)
print("B done (%.1fs): %d terms" % (time.time() - t0, len(termsB)))
A = sp.expand(A.rewrite(sp.cos)); B = sp.expand(B.rewrite(sp.cos))
for nm, e in (('A', A), ('B', B)):
    im = sp.expand(sp.im(e)) if e.has(sp.I) else 0
    print(nm, "has I:", e.has(sp.I), " nops:", sp.count_ops(e))
P = sp.expand(u2 * A - u1 * B)
print("P built (%.1fs), nops=%d" % (time.time() - t0, sp.count_ops(P)))
pickle.dump({'A': sp.srepr(A), 'B': sp.srepr(B), 'P': sp.srepr(P), 'termsA': {str(kk): sp.srepr(vv) for kk, vv in termsA.items()},
             'termsB': {str(kk): sp.srepr(vv) for kk, vv in termsB.items()}}, open('phi2_closed.pkl', 'wb'))
print("saved phi2_closed.pkl (%.1fs)" % (time.time() - t0))
