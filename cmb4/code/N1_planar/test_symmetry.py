"""Symmetry checks: (i) zeta -> -zeta == shift x -> x + pi  =>  Theta(-eps; x_o=0) = Theta(+eps; x_o=pi).
(ii) parity about x_o = 0: Theta(mu) = Theta(-mu)."""
import numpy as np
from planar import PlanarRun
P = dict(N=48, ic_order=6, nmu=64, verbose=False)
a = PlanarRun(-0.05, x_o=0.0, **P).run()
b = PlanarRun(+0.05, x_o=np.pi, **P).run()
print("shift symmetry: max|Theta(-eps,0) - Theta(+eps,pi)| =", np.max(np.abs(a.Theta - b.Theta)), " dt_o =", a.t_o - b.t_o, " dT_o/T_o =", a.T_o/b.T_o - 1)
print("parity (x_o=0): max|Theta(mu) - Theta(-mu)| =", np.max(np.abs(a.Theta - a.Theta[::-1])))
c = PlanarRun(+0.05, x_o=np.pi/2, **P).run()
print("generic observer x_o=pi/2: parity breaking max|Theta(mu)-Theta(-mu)| =", np.max(np.abs(c.Theta - c.Theta[::-1])), " a_10 =", c.alm(c.Theta, 2)[1])
