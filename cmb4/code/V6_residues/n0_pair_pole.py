"""n0 (V6): numerical check of the second-order pole with A3's pair profiles (C1 Green integrals):
P^2 x (AB/2) -> (-tau/eta, tau' + tau/eta, -4 tau/eta, -tau/eta [delta_T], -tau [V]) as P -> 0, |q1| = |q2| = 1."""
import sys; sys.path.insert(0, '../A3_third_order')
import numpy as np
from grid import TGrid, pair_profiles
from taufun import tau
g = TGrid(); t = g.t
ta, dta = tau(1.0, t), tau(1.0, t, 1)
exp = {'Phi': -ta / t, 'Psi': dta + ta / t, 'delta': -4 * ta / t, 'dT': -ta / t, 'V': -ta}
for P in (1e-2, 1e-3, 1e-4):
    pr = pair_profiles(g, np.array([P]), np.array([1.0]), np.array([1.0]))
    print('P=%g ' % P + '  '.join('%s: %.2e' % (nm, np.max(np.abs(P**2 * pr[nm][0].real / 2 - e))) for nm, e in exp.items()), '  (scale %.2f)' % np.max(np.abs(exp['Psi'])))
