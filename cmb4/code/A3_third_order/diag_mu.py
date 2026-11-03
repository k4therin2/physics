"""diag_mu.py: behaviour of the per-mu integrand of kappa near mu -> +1 (|k-q| -> 0) and mu -> -1 (|k+q| -> 0), k = q = 1."""
import numpy as np, sys, time
from thirteen import *
nu, wnu, nv = nu_grid(16)
t0 = time.time()
Mt = M2_tad(1.0)
print('M2tad', Mt, flush=True)
side = float(sys.argv[1])
for e in (1e-1, 1e-2, 1e-3, 1e-4):
    mu = side * (1 - e)
    ch = kappa_mu(1.0, 1.0, mu, nv, Nphi=8, M2tad=Mt, chan=True)
    tot = sum(ch.values())
    D = project(tot, nu, wnu, 6)
    big = sorted(ch.items(), key=lambda kv: -np.max(np.abs(project(kv[1], nu, wnu, 4)[2:])))[:5]
    print('mu=%+.5f  Dkappa_l(l=0..6) =' % mu, np.round(D.real if True else D, 4), np.round(D.imag, 4), ' x e:', np.round(np.abs(D[2:5]) * e, 4), '(%.0fs)' % (time.time() - t0), flush=True)
    for k_, v in big: print('      ', k_, np.round(project(v, nu, wnu, 4)[2:], 3))
