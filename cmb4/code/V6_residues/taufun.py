"""taufun.py (V6): the pole gauge function tau(eta) of the long second-order mode (pair |q1|=|q2|=s, total P -> 0):
   KPhi2 = -tau/(eta P^2) + O(P^0),  (eta^2 tau)' = -eta^2 A_Q,  A_Q = s^2 [2T^2 + (T + eta T')^2],
   tau = -sqrt3 s [x^3(2x - sin 2x) - 8x^2 cos^2 x + 8x sin 2x + 4 cos 2x - 4]/x^5,  x = s eta/sqrt3.
numpy callable tau(s, eta, d) (d = 0..3 eta-derivatives) with a small-x Taylor branch."""
import pickle, os
import numpy as np, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
_x, _s = sp.symbols('x s', positive=True)
_tau = eval(pickle.load(open(os.path.join(HERE, 's1_tau.pkl'), 'rb'))['tau_x'], sp.__dict__.copy()).subs(_s, 1) / sp.sqrt(3)  # tau/(sqrt3 s)
# tau(s,eta) = sqrt3 s * g(x), g = _tau ; d/deta = (s/sqrt3) d/dx
_ser = sp.series(_tau, _x, 0, 40).removeO()
_cs = [float(_ser.coeff(_x, n)) for n in range(40)]
_f = [sp.lambdify(_x, sp.diff(_tau, _x, d), 'numpy') for d in range(4)]
def _g(x, d):
    x = np.asarray(x, float); out = np.empty_like(x)
    sm = x < 1.5
    if np.any(sm):
        cc = list(_cs)
        for _ in range(d): cc = [n * cc[n] for n in range(1, len(cc))]
        v = np.zeros_like(x[sm])
        for a in reversed(cc): v = v * x[sm] + a
        out[sm] = v
    if np.any(~sm): out[~sm] = _f[d](x[~sm])
    return out
def tau(s, eta, d=0):
    eta = np.asarray(eta, float)
    return np.sqrt(3) * s * _g(s * eta / np.sqrt(3), d) * (s / np.sqrt(3))**d
if __name__ == '__main__':
    import mpmath as mp
    for e in (0.01, 0.5, 2.0, 8.0):
        xx = mp.mpf(e) / mp.sqrt(3)
        ex = -mp.sqrt(3) * (xx**3 * (2 * xx - mp.sin(2 * xx)) - 8 * xx**2 * mp.cos(xx)**2 + 8 * xx * mp.sin(2 * xx) + 4 * mp.cos(2 * xx) - 4) / xx**5
        print(e, float(tau(1.0, e)), float(ex), float(tau(1.0, e, 1)), float(mp.diff(lambda z: -mp.sqrt(3) * ((z/mp.sqrt(3))**3 * (2 * z/mp.sqrt(3) - mp.sin(2 * z/mp.sqrt(3))) - 8 * (z/mp.sqrt(3))**2 * mp.cos(z/mp.sqrt(3))**2 + 8 * z/mp.sqrt(3) * mp.sin(2 * z/mp.sqrt(3)) + 4 * mp.cos(2 * z/mp.sqrt(3)) - 4) / (z/mp.sqrt(3))**5, e)))
