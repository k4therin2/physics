"""n2_point (V6): per-channel Legendre coefficients (l=0..8) of A3's monochromatic mu-integrand at mu = side*(1-e)."""
import sys, pickle, time
sys.path.insert(0, '../A3_third_order')
import numpy as np
from thirteen import *
e = float(sys.argv[1]); side = int(sys.argv[2]) if len(sys.argv) > 2 else 1
nu, wnu, nv = nu_grid(32)
t0 = time.time()
Mt = M2_tad(1.0)
ch = kappa_mu(1.0, 1.0, side * (1 - e), nv, Nphi=16, M2tad=Mt, chan=True)
out = {k_: project(v, nu, wnu, 8) for k_, v in ch.items()}
pickle.dump(out, open('n2_e%.6g_s%+d.pkl' % (e, side), 'wb'))
print('e=%g side=%d done (%.0fs)' % (e, side, time.time() - t0))
