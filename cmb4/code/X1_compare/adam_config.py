"""X1: evaluate our C11, 1/4 C22, 1/3 Re C13 (l = 2..4) on Adam's fixed-ratio trajectory x_* = x_o/3 (his App. A tables),
x = k t/sqrt3 with t = conformal time.  Scaling: results depend on (k0 eta_*, k0 eta_o) only, so we keep eta_o = 8, eta_* = 8/3 and set
k0 = sqrt3 x_o / 8.  CAVEAT: our evaluators use eta_i -> 0 (observer initial event at the big bang), Adam uses x_I = x_o/6; the clock
enters at fourth order, so 1/4 C22 and 1/3 Re C13 differ by an O(eta_i^2)-suppressed amount (C11 is eta_i independent).
Usage: python3 adam_config.py 22 x_o [nK nth]  |  python3 adam_config.py 13 x_o [Nmu Nphi Nn]"""
import sys, os, time, pickle
import numpy as np
from numpy.polynomial import legendre as Lg
HERE = os.path.dirname(os.path.abspath(__file__))
mode, xo = sys.argv[1], float(sys.argv[2]); extra = [int(v) for v in sys.argv[3:]]
EO, ES = 8.0, 8.0 / 3.0
k0 = np.sqrt(3) * xo / EO
t0 = time.time()
if mode == '22':
    nK, nth = (extra + [24, 24][len(extra):])[:2]
    sys.path.insert(0, os.path.join(HERE, '..', 'A2_second_order'))
    import cl22
    r = cl22.run(k0=k0, eta_s=ES, eta_o=EO, nK=nK, nth=nth, lmax=4, tag=None, verbose=False)
    r['settings'] = (nK, nth)
    pickle.dump(r, open(os.path.join(HERE, 'ad22_%g_%d_%d.pkl' % (xo, nK, nth)), 'wb'))
    for l in range(2, 5):
        print('x_o=%g l=%d  C11=%.6e  C22/4=%.6e  (%.0fs)' % (xo, l, r['C11'][l], r['C22'][l] / 4, time.time() - t0))
else:
    Nmu, Nphi, Nn = (extra + [16, 8, 24][len(extra):])[:3]
    d = os.path.join(HERE, '..', 'A3_third_order'); sys.path.insert(0, d); os.chdir(d)
    import engine
    engine.RayGrid.__init__.__defaults__ = (ES, EO, 64)     # emission time eta_* = 8/3 instead of 2 (read-only patch in memory)
    from thirteen import kappa_mu, M2_tad, Delta_l, nu_grid, project, RAY
    assert abs(RAY.eta_s - ES) < 1e-14 and abs(RAY.eta_o - EO) < 1e-14
    nu, wnu, nv = nu_grid(Nn)
    xg, w = Lg.leggauss(Nmu); mus = (xg + 1) / 2; wmu = w / 2
    Mt = M2_tad(k0)
    chans = [kappa_mu(k0, k0, mu, nv, Nphi=Nphi, M2tad=Mt, chan=True) for mu in mus]
    names = list(chans[0].keys())
    tot = sum(sum(wmu[i] * chans[i][nm] for i in range(Nmu)) for nm in names)
    Dk = project(tot, nu, wnu, 4)
    Dl = np.array([Delta_l(l, k0) for l in range(5)])
    C13 = 4 * np.pi * Dl * np.conj(Dk)
    pickle.dump(dict(xo=xo, k0=k0, settings=(Nmu, Nphi, Nn), Delta=Dl, Dkappa=Dk, C13=C13),
                open(os.path.join(HERE, 'ad13_%g_%d_%d_%d.pkl' % (xo, Nmu, Nphi, Nn)), 'wb'))
    for l in range(2, 5):
        print('x_o=%g l=%d  Delta=%.8e  C11=%.6e  C13/3=%.6e  Im=%.1e  (%.0fs)' % (xo, l, Dl[l], 4 * np.pi * Dl[l]**2, C13[l].real / 3, C13[l].imag, time.time() - t0))
