"""n1ref.py -- Taylor coefficients of N1 scan quantities (T_o, t_o, <T_obs>, S(mu)) by even/odd least squares in eps."""
import numpy as np, os
N1 = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'N1_planar', 'output')
def load(tag):
    return np.load(os.path.join(N1, 'scan_%s.npz' % tag), allow_pickle=True)
def taylor(eps, Q, emax=0.08, nodd=3, neven=3):
    """Q[i,...] at eps[i] (+-) -> coefficients c0..c_{2max} (Q = sum c_n eps^n)."""
    eps = np.asarray(eps); Q = np.asarray(Q, dtype=float)
    pos = np.where((eps > 0) & (eps <= emax + 1e-12))[0]
    out = {}
    e = eps[pos]
    neg = np.array([np.where(np.isclose(eps, -x))[0][0] for x in e])
    O = (Q[pos] - Q[neg]) / 2; E = (Q[pos] + Q[neg]) / 2
    Ao = np.stack([e ** (2 * j + 1) for j in range(nodd)], 1)
    Ae = np.stack([e ** (2 * j) for j in range(neven)], 1)
    sh = Q.shape[1:]
    co = np.linalg.lstsq(Ao, O.reshape(len(e), -1), rcond=None)[0]
    ce = np.linalg.lstsq(Ae, E.reshape(len(e), -1), rcond=None)[0]
    for j in range(nodd): out[2 * j + 1] = co[j].reshape(sh)
    for j in range(neven): out[2 * j] = ce[j].reshape(sh)
    return out
