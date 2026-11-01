"""bianchi.py -- check  g^{mu a} nabla_a E_{mu nu} + g_{nu b} cons^b == 0  identically to O(eps^3).
(nabla^mu G_{mu nu}=0 is an identity; cons^b = nabla_mu (8piG T)^{mu b}.)  This verifies the engine's G and
shows that the time derivatives of the constraints are implied by the evolution + conservation equations
at every order (SPEC 9c / task control check (b))."""
import time
from tpoly import *
d = load_polys('E.pkl'); E, cons, Gam = d['E'], d['cons'], d['Gam']
import einstein_metric as em   # metric helpers (g, ginv) rebuilt cheaply
g, ginv = em.g, em.ginv
t0 = time.time()
for nu in range(4):
    tot = Poly({})
    for a in range(4):
        for m in range(4):
            if ginv[m][a].is_zero():
                continue
            cov = E[m][nu].diff(COORDS[a])
            for c in range(4):
                cov = cov - Gam[c][a][m] * E[c][nu] - Gam[c][a][nu] * E[m][c]
            tot = tot + ginv[m][a] * cov
    for b in range(4):
        tot = tot + g[nu][b] * cons[b]
    print("nu=%d: residual terms by order:" % nu, [tot.order_part(n).nterms() for n in range(4)], "(%.1fs)" % (time.time() - t0))
    assert tot.is_zero(), "Bianchi identity FAILED"
print("BIANCHI IDENTITY VERIFIED to O(eps^3): nabla^mu E_{mu nu} = -(8 pi G) nabla^mu T_{mu nu} identically.")
