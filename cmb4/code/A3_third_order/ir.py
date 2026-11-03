"""ir.py (A3): two-long-mode (q -> 0) limit of the contracted kernel.
Prediction (dilation, cf. A2's K2(k1->0) = -k.d_k K1): for zeta = c + zeta_s the map is Th[zeta_s](x -> e^c x), so
Th3 contains 3 c^2 (k.d_k)^2 K1 zeta_s, i.e. kappa~(q) := (1/2)<coef_abc(k,q,-q)>_qhat -> 3 (k d_k)^2 K1(k;n),
Delta^kappa~_l -> 3 (k d/dk)^2 Delta_l(k).  The O(q) (gradient/SCT) part is odd in qhat and drops out of the qhat average."""
import numpy as np, sys, time
from numpy.polynomial import legendre as Lg
from thirteen import *
k0 = 1.0
nu, wnu, nv = nu_grid(16)
qs = [float(a) for a in sys.argv[1:]] or [0.05]
h = 1e-3
lmax = 6
D = lambda k: np.array([Delta_l(l, k) for l in range(lmax + 1)])
pred = 3 * ((D(k0 * np.exp(h)) - 2 * D(k0) + D(k0 * np.exp(-h))) / h**2)     # (d/dlnk)^2
print('prediction 3 (k d/dk)^2 Delta_l, l=0..%d:' % lmax, np.round(pred, 5), flush=True)
xm, wm = Lg.leggauss(8)
t0 = time.time()
for qm in qs:
    Mt = M2_tad(qm)
    acc = 0
    for mu, w in zip(xm, wm):
        acc = acc + 0.5 * w * kappa_mu(k0, qm, mu, nv, Nphi=8, M2tad=Mt)
    Dk = project(acc, nu, wnu, lmax)
    print('q=%g  Delta^kappa~_l:' % qm, np.round(Dk.real, 5), np.round(Dk.imag, 5), ' ratio to prediction:', np.round((Dk / pred).real, 4), '(%.0fs)' % (time.time() - t0), flush=True)
