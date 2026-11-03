"""
plateau_isw.py -- agent C5: the RD 'Sachs-Wolfe plateau' with and without ISW.
Limit k eta_* -> 0 at fixed S = k eta_o (emission long before horizon entry of every observed mode):
   Delta_l -> -(1/3) j_l(S) + I_l^0(S),   I_l^0(S) = (4/sqrt3) int_0^S [j_2(s/sqrt3)/(s/sqrt3)] j_l(S-s) ds   (l >= 2),
a function of S alone.  For P = 2 pi^2 A_s / k^3:  C_l = 4 pi A_s int dS/S Delta_l(S)^2  -> exact plateau (eta_o-independent)
   l(l+1) C_l / (2 pi A_s) = 2 l(l+1) int_0^inf dS/S Delta_l(S)^2 .
SW alone gives 1/9; ISW (decay of Phi = Psi from -2zeta/3 to 0 after horizon entry, int 2 T_Phi' = +4/3) adds to it.
Large-l: the ISW source is concentrated at s = O(1) << S ~ l, so Delta_l ~ -(1/3) j_l(S) + (4/3) j_l(S) = j_l(S) (+ phase
smearing corrections O(1/l^{1/3}) ...) => l(l+1) C_l/(2 pi A_s) -> 1 ; checked numerically.
"""
import numpy as np
from scipy.special import spherical_jn
from linear_cl import composite_gl


def I0_grid(l, Sgrid):
    out = np.empty_like(Sgrid)
    for i, Sv in enumerate(Sgrid):
        e, w = composite_gl(0.0, Sv, 1.0, per_panel=40)
        xx = e / np.sqrt(3)
        src = (4 / np.sqrt(3)) * np.where(xx > 1e-3, spherical_jn(2, np.maximum(xx, 1e-3)) / np.maximum(xx, 1e-3), xx / 15)
        out[i] = np.sum(w * src * spherical_jn(l, Sv - e))
    return out


def plateau(l, Smax=None, n=6000):
    Smax = Smax or (60 + 40 * l)
    # log-spaced at small S (integrand ~ S^{2l-1}), linear later
    S1 = np.geomspace(1e-3, 2.0, 400); S2 = np.linspace(2.0, Smax, n)[1:]
    Sg = np.concatenate([S1, S2])
    I0 = I0_grid(l, Sg)
    D_sw = -spherical_jn(l, Sg) / 3
    D = D_sw + I0
    f = lambda y: np.trapezoid(y / Sg, Sg)
    tail = lambda y: 0.0
    full = 2 * l * (l + 1) * f(D**2)
    sw = 2 * l * (l + 1) * f(D_sw**2)
    isw = 2 * l * (l + 1) * f(I0**2)
    return sw, isw, full


if __name__ == '__main__':
    print("  l   SW only (->1/9=0.1111)   ISW only    SW+ISW   [l(l+1)C_l/(2 pi A_s), eta_* -> 0, RD]")
    for l in (2, 3, 4, 6, 10, 20, 40):
        sw, isw, full = plateau(l)
        print("  %2d   %.5f                  %.5f     %.5f" % (l, sw, isw, full))
