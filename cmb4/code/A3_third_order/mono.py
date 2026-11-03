"""mono.py (A3): monochromatic C^(13).  P = 2 pi^2 sigma^2 delta(k-k0)/k0^2, sigma = 1, k0 = 1, eta_* = 2, eta_o = 8.
kappa(k0;nu) = (1/2) int_{-1}^{1} dmu kappa_mu(mu;nu) = int_0^1 dmu kappa_mu  (integrand even in mu: (k,q,-q) symmetric),
Delta^kappa_l = (1/(2 i^l)) int dnu P_l kappa,   C13_l = 4 pi Delta_l(k0) (Delta^kappa_l)*   [units sigma^4].
Usage: python3 mono.py Nmu Nphi Nn tag"""
import numpy as np, sys, time, pickle
from numpy.polynomial import legendre as Lg
from thirteen import *
Nmu, Nphi, Nn, tag = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
k0 = 1.0
nu, wnu, nv = nu_grid(Nn)
x, w = Lg.leggauss(Nmu); mus = (x + 1) / 2; wmu = w / 2
t0 = time.time()
Mt = M2_tad(k0)
out = dict(mus=mus, wmu=wmu, nu=nu, wnu=wnu, M2tad=Mt, chans=[], Nphi=Nphi)
for i, mu in enumerate(mus):
    ch = kappa_mu(k0, k0, mu, nv, Nphi=Nphi, M2tad=Mt, chan=True)
    out['chans'].append(ch)
    print('mu %d/%d = %.5f done (%.0fs)' % (i + 1, Nmu, mu, time.time() - t0), flush=True)
    pickle.dump(out, open('mono_%s.pkl' % tag, 'wb'))
names = list(out['chans'][0].keys())
kap = {nm: sum(wmu[i] * out['chans'][i][nm] for i in range(Nmu)) for nm in names}
tot = sum(kap.values())
lmax = 20
Dk = project(tot, nu, wnu, lmax)
Dl = np.array([Delta_l(l, k0) for l in range(lmax + 1)])
C13 = 4 * np.pi * Dl * np.conj(Dk)
out.update(kappa=kap, Dkappa=Dk, Delta=Dl, C13=C13)
pickle.dump(out, open('mono_%s.pkl' % tag, 'wb'))
print(' l   Delta_l            Delta^kappa_l                 C13_l/sigma^4 (Re, Im)')
for l in range(lmax + 1):
    print('%2d  %+.9e  %+.9e %+.2e i  %+.9e %+.2e' % (l, Dl[l], Dk[l].real, Dk[l].imag, C13[l].real, C13[l].imag))
