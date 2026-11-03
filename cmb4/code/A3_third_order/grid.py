"""grid.py (A3) -- time grids, Chebyshev interpolation, Green's-function integrals, second-order pair profiles.
Global grid: Chebyshev first-kind nodes on (0, E], E = eta_o (no node at eta = 0).  Every time function of a
Fourier mode is stored as its node values; it is evaluated elsewhere by spectral (Chebyshev) interpolation.
Green's integrals  F(eta_i) = int_0^{eta_i} G(K; eta_i, eta') S(eta') deta'  by Gauss-Legendre with Ng nodes on [0, eta_i];
G_Phi uses a small-K*eta Taylor branch (the closed form u1/u2 cancels catastrophically for K eta << 1).
"""
import numpy as np
import sympy as sp
from numpy.polynomial import chebyshev as Ch, legendre as Lg
import tf

SQ3 = np.sqrt(3.0)

# ---------------- G_Phi with series branch
_K, _e, _ep = sp.symbols('K eta etap', positive=True)
_GPhi = tf._K2['G_Phi'].subs({tf.KS: _K, tf.ETA: _e, sp.Symbol('etap', positive=True): _ep})
NGS = 28
_ser = sp.series(_GPhi, _K, 0, NGS).removeO()
_serd = sp.diff(_ser, _e)
_Gser = sp.lambdify((_K, _e, _ep), _ser, 'numpy')
_dGser = sp.lambdify((_K, _e, _ep), _serd, 'numpy')
_Gcl = sp.lambdify((_K, _e, _ep), _GPhi, 'numpy')
_dGcl = sp.lambdify((_K, _e, _ep), sp.diff(_GPhi, _e), 'numpy')
def GPhi(K, e, ep, d=0):
    K, e, ep = np.broadcast_arrays(np.asarray(K, float), np.asarray(e, float), np.asarray(ep, float))
    out = np.empty(K.shape)
    sm = K * e < 2.0
    fs, fc = (_Gser, _Gcl) if d == 0 else (_dGser, _dGcl)
    if np.any(sm): out[sm] = fs(K[sm], e[sm], ep[sm])
    if np.any(~sm): out[~sm] = fc(K[~sm], e[~sm], ep[~sm])
    return out
def Gh(K, e, ep, d=0):
    K, e, ep = np.broadcast_arrays(np.asarray(K, float), np.asarray(e, float), np.asarray(ep, float))
    tau = e - ep
    if d == 0:
        return ep / e * tau * np.sinc(K * tau / np.pi)
    # d/deta [ (ep/eta) sin(K tau)/K ] = -(ep/eta^2) sin(K tau)/K + (ep/eta) cos(K tau)
    return -ep / e**2 * tau * np.sinc(K * tau / np.pi) + ep / e * np.cos(K * tau)

class TGrid:
    def __init__(self, E=8.0, N=96, Ng=72):
        self.E, self.N, self.Ng = float(E), N, Ng
        j = np.arange(N)
        xg = -np.cos(np.pi * (j + 0.5) / N)                      # ascending first-kind nodes in (-1,1)
        self.xg = xg
        self.t = (xg + 1) * E / 2
        self.V = Ch.chebvander(xg, N - 1)
        self.Vinv = np.linalg.inv(self.V)
        # quadrature weights int_0^E f = w . f
        e = Ch.chebint(np.eye(N), lbnd=-1)                         # columns: integrals of T_n
        self.w = (Ch.chebval(1.0, e) @ self.Vinv) * E / 2
        u, wu = Lg.leggauss(Ng)
        self.s = (u[None, :] + 1) * self.t[:, None] / 2            # (N,Ng)
        self.gw = wu[None, :] * self.t[:, None] / 2
        self.Isub = self.interp_matrix(self.s.ravel())             # (N*Ng, N)
    def interp_matrix(self, pts):
        x = 2 * np.asarray(pts, float) / self.E - 1
        return Ch.chebvander(x, self.N - 1) @ self.Vinv
    def sub(self, f):
        """node values (..., N) -> values at Green sub-nodes (..., N, Ng)."""
        return (f @ self.Isub.T).reshape(f.shape[:-1] + (self.N, self.Ng))
    def green(self, kind, K, Ssub, d=0):
        """kind 'Phi' or 'h'; K (...,) ; Ssub (..., N, Ng) -> (..., N)."""
        K = np.asarray(K, float)[..., None, None]
        G = GPhi(K, self.t[:, None], self.s, d) if kind == 'Phi' else Gh(K, self.t[:, None], self.s, d)
        return np.sum(self.gw * G * Ssub, axis=-1)
    def integ(self, f):
        return f @ self.w

# ---------------- extra lambdified kernels (eta derivatives of algebraic second-order kernels)
_d = {}
for nm in ('K_Q', 'K_S00', 'K_S0z', 'K_om', 'K_vT'):
    _d['d' + nm] = sp.lambdify((tf.KS, tf.QS, tf.PS, tf.ETA), tf._prep(sp.diff(tf._K2[nm], tf.ETA)), modules=tf._mods)
KER2 = dict(tf.KER2); KER2.update(_d)

def pair_profiles(g, K, q1, q2):
    """Second-order fields of a momentum pair (q1,q2), K = |q1+q2| > 0, as AB-COEFFICIENTS (= K_X(q1,q2) + K_X(q2,q1)),
    node values on the global grid g.  K, q1, q2: arrays (B,).  Returns dict of (B, N) arrays:
      Phi, dPhi, Psi, dPsi, dT, ddT (eta-derivative of delta_T2), V, dV (longitudinal: v2_i = i K_i V + transverse),
      vT12, dvT12 (coefficient of P_ij(K) q1_j), vT21, dvT21 (coefficient of P_ij q2_j), om12, dom12, om21, dom21,
      H, dH (h2_ij = Lambda_ij,kl(K) q1_k q1_l H)."""
    K, q1, q2 = [np.atleast_1d(np.asarray(a, float)) for a in (K, q1, q2)]
    t = g.t[None, :]
    Kb, ab, bb = K[:, None], q1[:, None], q2[:, None]
    Ks, as_, bs = K[:, None, None], q1[:, None, None], q2[:, None, None]
    out = {}
    sym = lambda nm, A, B, C, T: KER2[nm](A, B, C, T) + KER2[nm](A, C, B, T)
    S = sym('K_SPhi', Ks, as_, bs, g.s)
    Phi = g.green('Phi', K, S) + (KER2['K_Phi2_0'](Kb, ab, bb, t) + KER2['K_Phi2_0'](Kb, bb, ab, t)) * tf.ghat(Kb, t)
    dPhi = g.green('Phi', K, S, 1) + (KER2['K_Phi2_0'](Kb, ab, bb, t) + KER2['K_Phi2_0'](Kb, bb, ab, t)) * tf.ghat(Kb, t, 1)
    St = sym('K_SPhi', Kb, ab, bb, t)
    ddPhi = St - 4 * dPhi / t - Kb**2 / 3 * Phi
    Q, dQ = sym('K_Q', Kb, ab, bb, t), sym('dK_Q', Kb, ab, bb, t)
    S00, dS00 = sym('K_S00', Kb, ab, bb, t), sym('dK_S00', Kb, ab, bb, t)
    S0z, dS0z = sym('K_S0z', Kb, ab, bb, t), sym('dK_S0z', Kb, ab, bb, t)
    Psi, dPsi = Phi - Q, dPhi - dQ
    delta = t**2 / 3 * S00 - 2 / 3 * t**2 * Kb**2 * Phi - 2 * t * dPhi - 2 * Psi
    ddelta = 2 * t / 3 * S00 + t**2 / 3 * dS00 - 4 / 3 * t * Kb**2 * Phi - 2 / 3 * t**2 * Kb**2 * dPhi - 2 * dPhi - 2 * t * ddPhi - 2 * dPsi
    d1a, d1b = tf.tdelta(ab, t), tf.tdelta(bb, t)
    dT = delta / 4 - 3 / 16 * 2 * d1a * d1b
    ddT = ddelta / 4 - 3 / 16 * 2 * (tf.tdelta(ab, t, 1) * d1b + d1a * tf.tdelta(bb, t, 1))
    # V = v_parallel/(iK):  V = -(eta^2/4)[2(Phi' + Psi/eta) + S0z/(iK)]
    V = -(t**2 / 4) * (2 * (dPhi + Psi / t) + S0z / (1j * Kb))
    dV = -(t / 2) * (2 * (dPhi + Psi / t) + S0z / (1j * Kb)) - (t**2 / 4) * (2 * (ddPhi + dPsi / t - Psi / t**2) + dS0z / (1j * Kb))
    out.update(Phi=Phi, dPhi=dPhi, ddPhi=ddPhi, Psi=Psi, dPsi=dPsi, delta=delta, dT=dT, ddT=ddT, V=V, dV=dV, Q=Q)
    out['vT12'], out['dvT12'] = KER2['K_vT'](Kb, ab, bb, t), KER2['dK_vT'](Kb, ab, bb, t)
    out['vT21'], out['dvT21'] = KER2['K_vT'](Kb, bb, ab, t), KER2['dK_vT'](Kb, bb, ab, t)
    out['om12'], out['dom12'] = KER2['K_om'](Kb, ab, bb, t), KER2['dK_om'](Kb, ab, bb, t)
    out['om21'], out['dom21'] = KER2['K_om'](Kb, bb, ab, t), KER2['dK_om'](Kb, bb, ab, t)
    Sh = sym('K_h', Ks, as_, bs, g.s)
    out['H'], out['dH'] = g.green('h', K, Sh), g.green('h', K, Sh, 1)
    Sht = sym('K_h', Kb, ab, bb, t)
    out['ddH'] = Sht - 2 / t * out['dH'] - Kb**2 * out['H']
    return out

def zero_profiles(g, q):
    """(q,-q) pair (K = 0): C1 angle-averaged tadpoles, AB-coefficient = 2 x tadpole density."""
    q = np.atleast_1d(np.asarray(q, float))[:, None]; t = g.t[None, :]
    tP, dtP = tf.TAD['tadPhi2'](q, t), tf.TAD['tadPhi2'](q, t, 1)
    ddtP = tf.TAD['tadPhi2'](q, t, 2)
    tQ, dtQ = tf.TAD['tadQ2'](q, t), tf.TAD['tadQ2'](q, t, 1)
    td, dtd = tf.TAD['taddelta2'](q, t), tf.TAD['taddelta2'](q, t, 1)
    z = np.zeros_like(tP + 0j)
    out = dict(Phi=2 * tP, dPhi=2 * dtP, ddPhi=2 * ddtP, Psi=2 * (tP - tQ), dPsi=2 * (dtP - dtQ), delta=2 * td, Q=2 * tQ,
               dT=2 * td / 4 - 3 / 16 * 2 * tf.tdelta(q, t)**2, ddT=2 * dtd / 4 - 3 / 16 * 4 * tf.tdelta(q, t) * tf.tdelta(q, t, 1),
               V=z, dV=z, vT12=z, dvT12=z, vT21=z, dvT21=z, om12=z, dom12=z, om21=z, dom21=z, H=z, dH=z, ddH=z)
    return out

if __name__ == '__main__':
    g = TGrid()
    # Green's function checks: Phi operator on f = Green integral of S=1 (constant source): F'' + 4F'/eta + K^2 F/3 = 1
    K = np.array([1e-3, 0.5, 1.7])
    S = np.ones((3, g.N, g.Ng))
    F = g.green('Phi', K, S); dF = g.green('Phi', K, S, 1)
    # finite-difference check of the ODE at a few nodes via spectral differentiation
    c = (F @ g.Vinv.T); dc = Ch.chebder(c.T).T * 2 / g.E; ddc = Ch.chebder(c.T, 2).T * (2 / g.E)**2
    Fp = Ch.chebvander(g.xg, g.N - 1)[:, :g.N - 1] @ dc.T; Fpp = Ch.chebvander(g.xg, g.N - 1)[:, :g.N - 2] @ ddc.T
    res = Fpp.T + 4 * Fp.T / g.t - 1 + K[:, None]**2 / 3 * F
    print('Phi Green ODE residual (interior nodes):', np.max(np.abs(res[:, 10:-10])), ' dF vs spectral:', np.max(np.abs(dF - Fp.T)[:, 5:-5]))
    # tadpole consistency: taddelta2 vs (eta^2/3) tadS00 - 2 eta tadPhi2' - 2 (tadPhi2 - tadQ2)
    q = 1.0; t = g.t
    lhs = tf.TAD['taddelta2'](q, t)
    tS00 = sp.lambdify((tf._qs, tf._es), tf.TADSYM['tadS00'], 'numpy')(q, t)
    rhs = t**2 / 3 * tS00 - 2 * t * tf.TAD['tadPhi2'](q, t, 1) - 2 * (tf.TAD['tadPhi2'](q, t) - tf.TAD['tadQ2'](q, t))
    print('tadpole delta2 constraint residual:', np.max(np.abs(lhs - rhs)[t > 0.5]))
    # angle-averaged K->0 limit of pair profiles vs tadpoles: q1 = q zhat, q2 = -q1 + Kvec, average over Kvec directions
    for Kmag in (0.02, 0.005):
        mus, wm = Lg.leggauss(40)
        q2 = np.sqrt(q**2 + Kmag**2 - 2 * q * Kmag * mus)          # Kvec.q1hat = Kmag mu
        pr = pair_profiles(g, np.full(40, Kmag), np.full(40, q), q2)
        z = zero_profiles(g, q)
        for nm in ('Phi', 'Psi', 'delta', 'dT'):
            avg = 0.5 * np.sum(wm[:, None] * pr[nm], axis=0)
            print('K=%g %s: angle-avg pair / zero-mode at eta=2,8:' % (Kmag, nm), np.interp([2, 8], t, avg.real), np.interp([2, 8], t, z[nm][0].real))
