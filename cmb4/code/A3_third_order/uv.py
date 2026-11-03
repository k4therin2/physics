"""uv.py (A3): large-q behaviour of the q-integrand kappa~(q) = (1/2)<coef_abc(k,q,-q)>_qhat at fixed k = 1 (per P(q)).
Usage: python3 uv.py Nmu Nphi Nt q1 q2 ..."""
import numpy as np, sys, time
from numpy.polynomial import legendre as Lg
from thirteen import *
Nmu, Nphi, Nt = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
qs = [float(a) for a in sys.argv[4:]]
k0 = 1.0; lmax = 6
nu, wnu, nv = nu_grid(16)
xm, wm = Lg.leggauss(Nmu)
sky = (Nt, 2 * Nt)
t0 = time.time()
for qm in qs:
    Mt = M2_tad(qm, sky)
    acc = 0; chacc = {}
    for mu, w in zip(xm, wm):
        ch = kappa_mu(k0, qm, mu, nv, Nphi=Nphi, M2tad=Mt, sky=sky, chan=True)
        for k_, v in ch.items(): chacc[k_] = chacc.get(k_, 0) + 0.5 * w * v
    acc = sum(chacc.values())
    Dk = project(acc, nu, wnu, lmax)
    print('q=%g  Delta^kappa~_l (l=0..%d):' % (qm, lmax), np.round(Dk.real, 5), '(%.0fs)' % (time.time() - t0), flush=True)
    big = sorted(chacc.items(), key=lambda kv: -abs(project(kv[1], nu, wnu, 2)[2]))[:6]
    for k_, v in big: print('      %-60s l=2: %+.4e' % (k_, project(v, nu, wnu, 2)[2].real), flush=True)
