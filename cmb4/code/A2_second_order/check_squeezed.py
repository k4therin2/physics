"""check_squeezed.py -- CHECK (i): constant mode.  zeta_nl = c + zeta_s is the coordinate rescaling x~ = e^{c} x of the
configuration zeta~_s(x~) = zeta_s(e^{-c} x~) (spatial metric e^{2zeta} dx^2 at eta -> 0; observer at x = 0, proper time,
emission temperature are invariant).  Hence Theta[c + zeta_s] = Theta_1[zeta~_s] + ...  and, with Theta = eps Th1 + eps^2 Th2/2,
     lim_{k1 -> 0} K_2(k1, k; n) = - k . grad_k K_1(k; n) = - k d/dk K_1(k khat; n)      (for l >= 1; consistency relation).
So the squeezed limit does NOT vanish channel-by-channel nor in total; it equals the dilation of the linear kernel.
(The brief's 'must vanish' holds for the zeta-independent part only after subtracting this rescaling: the constant mode is
unobservable in the sense that it can be absorbed in a relabelling of the short-mode wavenumbers.)  Verified numerically
here for l >= 1 (l = 0 also, up to the normalisation convention)."""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kernel2 import Geometry, MapKernels
from sphere import Sky
g = Geometry(2.0, 8.0)
sky = Sky(40); nv = sky.nv; L = 30
e = np.array([0.3, 0.5, -0.8]); e /= np.linalg.norm(e)
for k2 in (np.array([-0.3, 0.9, 0.5]), np.array([1.2, -2.0, 0.7])):
    h = 1e-4
    mk = MapKernels(g, k2 * (1 + h), k2 * (1 - h), nv)
    dil = -(mk.K1(mk.sa) - mk.K1(mk.sb)) / (2 * h)          # -k d/dk K_1
    Pd = sky.power(dil, L)
    print('|k2| = %.3f  dilation -k dK1/dk: l>=2 rms %.4e, l=1 %.4e, l=0 %.4e' % (np.linalg.norm(k2), np.sqrt(np.sum(Pd[2:])), np.sqrt(Pd[1]), np.sqrt(Pd[0])))
    for s in (1e-2, 1e-3, 1e-4):
        K, ch = MapKernels(g, s * e, k2, nv).K2()
        R = sky.power(K - dil, L)
        print('   s=%.0e  |K_2(s e,k2) - (-k dK1/dk)|:  l>=2 %.3e  l=1 %.3e  l=0 %.3e' % (s, np.sqrt(np.sum(R[2:])), np.sqrt(R[1]), np.sqrt(R[0])))
