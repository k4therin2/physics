"""superhorizon.py (W4): superhorizon limit of the monochromatic O(eps^4) coefficient.
Uses A2's cl22.run and A3's thirteen.kappa_mu read-only (no A2/A3 file modified).
eta_* = 2, eta_o = 8 fixed; k0 = x/eta_o with x = k0 eta_o (the results depend on k0 eta_*, k0 eta_o only).
Usage: python3 superhorizon.py 22 x [nK nth]      -> sh22_<x>.pkl  (C11_l, C22_l)
       python3 superhorizon.py 13 x [Nmu Nphi Nn] -> sh13_<x>.pkl  (Delta_l, Delta^kappa_l, C13_l)"""
import sys, os, time, pickle
import numpy as np
from numpy.polynomial import legendre as Lg
HERE = os.path.dirname(os.path.abspath(__file__))
mode, x = sys.argv[1], float(sys.argv[2])
extra = [int(v) for v in sys.argv[3:]]
k0 = x / 8.0
t0 = time.time()
if mode == '22':
    nK, nth = (extra + [16, 16][len(extra):])[:2]
    d = os.path.join(HERE, '..', 'A2_second_order'); sys.path.insert(0, d)
    import cl22
    r = cl22.run(k0=k0, eta_s=2.0, eta_o=8.0, nK=nK, nth=nth, lmax=6, tag=None, verbose=True)
    r['settings'] = (nK, nth)
    pickle.dump(r, open(os.path.join(HERE, 'sh22_%g_%d_%d.pkl' % (x, nK, nth)), 'wb'))
    for l in range(2, 7):
        print('l=%d  C11=%.6e  C22=%.6e  C22/4=%.6e' % (l, r['C11'][l], r['C22'][l], r['C22'][l] / 4))
else:
    Nmu, Nphi, Nn = (extra + [8, 8, 16][len(extra):])[:3]
    d = os.path.join(HERE, '..', 'A3_third_order'); sys.path.insert(0, d); os.chdir(d)
    from thirteen import kappa_mu, M2_tad, Delta_l, nu_grid, project
    nu, wnu, nv = nu_grid(Nn)
    xg, w = Lg.leggauss(Nmu); mus = (xg + 1) / 2; wmu = w / 2
    Mt = M2_tad(k0)
    chans = []
    for i, mu in enumerate(mus):
        chans.append(kappa_mu(k0, k0, mu, nv, Nphi=Nphi, M2tad=Mt, chan=True))
        print('mu %d/%d = %.5f (%.0fs)' % (i + 1, Nmu, mu, time.time() - t0), flush=True)
    names = list(chans[0].keys())
    kap = {nm: sum(wmu[i] * chans[i][nm] for i in range(Nmu)) for nm in names}
    tot = sum(kap.values())
    lmax = 6
    Dk = project(tot, nu, wnu, lmax)
    Dl = np.array([Delta_l(l, k0) for l in range(lmax + 1)])
    C13 = 4 * np.pi * Dl * np.conj(Dk)
    chanD = {nm: project(v, nu, wnu, lmax) for nm, v in kap.items()}
    r = dict(x=x, k0=k0, settings=(Nmu, Nphi, Nn), M2tad=Mt, Dkappa=Dk, Delta=Dl, C13=C13, chanD=chanD, mus=mus)
    pickle.dump(r, open(os.path.join(HERE, 'sh13_%g_%d_%d.pkl' % (x, Nmu, Nphi)), 'wb'))
    for l in range(2, 7):
        print('l=%d  Delta=%.6e  Dkappa=%.6e  C13/3=%.6e  Im=%.1e  maxchan=%.3e' % (l, Dl[l], Dk[l].real, C13[l].real / 3, C13[l].imag,
              max(abs(v[l]) for v in chanD.values())))
