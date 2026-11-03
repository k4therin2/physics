"""X1: superhorizon runs at a different emission time (eta_o = 8 fixed, eta_* = ES), to locate the O(k0^6) difference in 1/3 Re C13.
Usage: python3 sh_ratio.py ES x   (runs both 22 and 13 at k0 = x/8; (13) with Nmu, Nphi, Nn = 8, 8, 16 as in W4)"""
import sys, os, pickle, numpy as np
from numpy.polynomial import legendre as Lg
HERE = os.path.dirname(os.path.abspath(__file__))
ES, x = float(sys.argv[1]), float(sys.argv[2]); EO = 8.0; k0 = x / EO
sys.path.insert(0, os.path.join(HERE, '..', 'A2_second_order'))
import cl22
r22 = cl22.run(k0=k0, eta_s=ES, eta_o=EO, nK=16, nth=16, lmax=3, tag=None, verbose=False)
d = os.path.join(HERE, '..', 'A3_third_order'); sys.path.insert(0, d); os.chdir(d)
import engine
engine.RayGrid.__init__.__defaults__ = (ES, EO, 64)
from thirteen import kappa_mu, M2_tad, Delta_l, nu_grid, project, RAY
assert abs(RAY.eta_s - ES) < 1e-14
nu, wnu, nv = nu_grid(16); xg, w = Lg.leggauss(8); mus = (xg + 1) / 2; wmu = w / 2
Mt = M2_tad(k0)
tot = sum(wmu[i] * sum(kappa_mu(k0, k0, mu, nv, Nphi=8, M2tad=Mt, chan=True).values()) for i, mu in enumerate(mus))
Dk = project(tot, nu, wnu, 3); Dl = np.array([Delta_l(l, k0) for l in range(4)]); C13 = 4 * np.pi * Dl * np.conj(Dk)
out = dict(ES=ES, x=x, C11=r22['C11'], C22=r22['C22'], C13=C13, Delta=Dl)
pickle.dump(out, open(os.path.join(HERE, 'shr_%g_%g.pkl' % (ES, x)), 'wb'))
l = 2
print('ES=%g x=%g  C11=%.10e  (C22/4)/C11=%.10f  (ReC13/3)/C11=%.10f' % (ES, x, r22['C11'][l], r22['C22'][l] / 4 / r22['C11'][l], C13[l].real / 3 / r22['C11'][l]))
