"""Residual test against N1's raw scan (no N1 fitting): R(eps) = Theta_N1(eps) - [eps Th1 + eps^2 Th2/2 + eps^3 Th3/6]_A1 must be O(eps^4)."""
import numpy as np, sys
from n1ref import load
tag = sys.argv[1] if len(sys.argv) > 1 else 'to8_xo0'
m = np.load('map_%s_N64.npz' % tag); s = load(tag)
for i in np.argsort(s['eps']):
    e = s['eps'][i]
    R = s['Theta'][i] - (e * m['Th1'] + e**2 / 2 * m['Th2'] + e**3 / 6 * m['Th3'])
    print('eps=% .2f  max|R| = %.3e   max|R|/eps^4 = %.4f' % (e, np.max(np.abs(R)), np.max(np.abs(R)) / e**4))
