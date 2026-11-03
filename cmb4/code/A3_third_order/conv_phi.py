import numpy as np, sys, time
from thirteen import *
nu, wnu, nv = nu_grid(24)
Mt = M2_tad(1.0)
mu = float(sys.argv[1])
for Nphi in (8, 16, 32, 48):
    t0 = time.time()
    v = kappa_mu(1.0, 1.0, mu, nv, Nphi=Nphi, M2tad=Mt)
    print('mu=%g Nphi=%d  Dk_l(l=2..6):' % (mu, Nphi), np.round(project(v, nu, wnu, 6)[2:].real, 6), '(%.0fs)' % (time.time() - t0), flush=True)
