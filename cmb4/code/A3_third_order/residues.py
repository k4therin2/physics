"""residues.py (A3): mu -> 1 behaviour of the monochromatic integrand, channel by channel.
For each channel c: Delta^kappa_mu,c,l(mu) at 1-mu = e; residue R_c,l = lim e Delta_c,l (coefficient of 1/(1-mu)); finite part.
(mu -> -1 is identical: integrand even in mu.)"""
import numpy as np, time, pickle
from thirteen import *
nu, wnu, nv = nu_grid(24)
Mt = M2_tad(1.0)
res = {}
for e in (1e-4, 1e-5, 1e-6):
    t0 = time.time()
    ch = kappa_mu(1.0, 1.0, 1 - e, nv, Nphi=16, M2tad=Mt, chan=True)
    res[e] = {k_: project(v, nu, wnu, 6).real for k_, v in ch.items()}
    print('done e=%g (%.0fs)' % (e, time.time() - t0), flush=True)
pickle.dump(res, open('residues.pkl', 'wb'))
names = list(res[1e-4].keys())
print('%-70s %s' % ('channel', 'residue e*Delta_l (l=0..6) from e=1e-6'))
for k_ in names:
    r = 1e-6 * res[1e-6][k_]
    if np.max(np.abs(r)) > 1e-6: print('%-70s %s' % (k_, np.round(r, 5)))
print('%-70s %s' % ('SUM of residues', np.round(sum(1e-6 * res[1e-6][k_] for k_ in names), 7)))
for e in (1e-4, 1e-5, 1e-6):
    print('total Delta^kappa_mu,l at 1-mu=%g:' % e, np.round(sum(res[e][k_] for k_ in names), 5))
