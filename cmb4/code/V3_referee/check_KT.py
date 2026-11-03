"""check_KT.py -- V3: K_h (C1) vs Kohri-Terada (20) f_RD transcribed from literature_extracts/raw text.
Mapping: KT g_ij > (1/2) h^KT, ours g_ij > h = eps^2 h2/2 => h^KT = h2;  KT h'' + 2H h' + k^2 h = 4 S_k,
S_k = e^{ij} q_i q_j [2 Phi Phi + (Phi + eta Phi')(...)] (w = 1/3), Phi = T zeta (T -> -2/3).  KT's f is written with
their Phi(x) -> 1 normalisation and the (2/3)^2 absorbed:  f_KT(v,u,x) = (4/9)[2 PP + (P+xP')(P+xP')]|_{P->1} = [2TT + ...]_{T->-2/3}.
So expect K_h(q,p,k;eta) = 4 f_RD(v=q/k, u=p/k, x=k eta) exactly."""
import sys, mpmath as mp, sympy as sp
C1 = '/home/user/physics/cmb4/code/C1_einstein'; sys.path.insert(0, C1)
import kernels_api as KA; K = KA.load()
from fourier import eta, k, q, p, mu
Kh = sp.lambdify((k, q, p, eta), KA.closed_form(K['second']['K_h']), 'mpmath')
mp.mp.dps = 30
def fRD(v, u, x):
    s3 = mp.sqrt(3); cu, su, cv, sv = mp.cos(u*x/s3), mp.sin(u*x/s3), mp.cos(v*x/s3), mp.sin(v*x/s3)
    return 12/(u**3*v**3*x**6)*(18*u*v*x**2*cu*cv + (54 - 6*(u**2+v**2)*x**2 + u**2*v**2*x**4)*su*sv
                                + 2*s3*u*x*(v**2*x**2 - 9)*cu*sv + 2*s3*v*x*(u**2*x**2 - 9)*su*cv)
import random
rng = random.Random(5); worst = 0
for i in range(6):
    kk, qq, mm, ee = [mp.mpf(rng.uniform(a, b)) for a, b in ((0.2, 3), (0.2, 3), (-0.95, 0.95), (0.1, 5))]
    pp = mp.sqrt(kk**2 + qq**2 - 2*kk*qq*mm)
    a = Kh(kk, qq, pp, ee); b = 4*fRD(qq/kk, pp/kk, kk*ee)
    worst = max(worst, abs(a - b)/abs(b)); print("K_h = %s   4 f_RD = %s" % (mp.nstr(a, 18), mp.nstr(b, 18)))
print("worst relative difference:", mp.nstr(worst, 5))
