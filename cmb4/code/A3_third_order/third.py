"""third.py (A3) -- contracted third-order fields X3c(k, q, mu; eta) from C1's K3c sources (channel iii).
Per P(q) (C1 normalisation: <X3_k zeta_k'> = (2pi)^3 delta P(k) X3c,  X3c = int d^3q/(2pi)^3 P(q) K3c_X):
   Q3c = K3c_Q ;  Phi3c = int_0^eta G_Phi(k) K3c_SPhi + C3(k;q,mu) ghat(k eta) ;  Psi3c = Phi3c - Q3c ;
   delta3c = (eta^2/3) K3c_S00 - (2/3) eta^2 k^2 Phi3c - 2 eta Phi3c' - 2 Psi3c ;
   V3c = -(eta^2/4)[2 (Phi3c' + Psi3c/eta) + K3c_S0z/(i k)]   (v3_i = i k_i V3c; transverse/TT contracted parts vanish by isotropy).
Inputs: KPhi2(q,k,p,eta) = C1's per-pair Phi2 kernel = (AB-coefficient)/2 of the pair (q1,q2)=(q,k) [symmetric kernels, checked],
Kh likewise; tadpoles tadPhi2, tadQ2, taddelta2 (C1, angle-averaged) -- optionally shifted (gauge test)."""
import numpy as np
from numpy.polynomial import chebyshev as Ch
import tf
from grid import TGrid, pair_profiles

def spectral_d(g, f):
    c = f @ g.Vinv.T
    dc = Ch.chebder(c, axis=-1) * 2 / g.E
    return dc @ Ch.chebvander(g.xg, g.N - 2).T

def X3c(g, k, q, mu, tads=None, C3shift=None):
    """k, q scalars; mu array (B,) -> dict of (B, N) global-node profiles per P(q)."""
    mu = np.atleast_1d(np.asarray(mu, float)); B = len(mu)
    p = np.sqrt(k**2 + q**2 - 2 * k * q * mu)
    pr = pair_profiles(g, p, np.full(B, q), np.full(B, k))
    sub = lambda f: g.sub(f)
    t = g.t
    if tads is None:
        tads = {nm: (lambda d, nm=nm: tf.TAD[nm](q, t, d)) for nm in ('tadPhi2', 'tadQ2', 'taddelta2')}
    def args(where):
        """argument arrays for KER3 at global nodes ('n') or Green sub-nodes ('s')."""
        T = t if where == 'n' else g.s
        conv = (lambda f: f) if where == 'n' else sub
        bc = lambda a: np.broadcast_to(a, (B,) + np.shape(T))
        P0 = conv(pr['Phi'] / 2); P1 = conv(pr['dPhi'] / 2); P2 = conv(pr['ddPhi'] / 2)
        H0 = conv(pr['H'] / 2); H1 = conv(pr['dH'] / 2); H2 = conv(pr['ddH'] / 2)
        # H3 from differentiating h'' + (2/eta) h' + p^2 h = K_h:  h''' = K_h' - (2/eta) h'' + (2/eta^2) h' - p^2 h'
        pb = p.reshape((B,) + (1,) * np.ndim(T))
        dKh = 0.5 * (tf.KER2['dK_h'](pb, q, k, T) + tf.KER2['dK_h'](pb, k, q, T))
        H3 = dKh - 2 / T * H2 + 2 / T**2 * H1 - pb**2 * H1
        tP = [conv(np.broadcast_to(tads['tadPhi2'](d), (B, g.N))) for d in (0, 1, 2)]
        tQ = [conv(np.broadcast_to(tads['tadQ2'](d), (B, g.N))) for d in (0, 1)]
        td0 = conv(np.broadcast_to(tads['taddelta2'](0), (B, g.N)))
        Tk = [bc(tf.tphi(k, T, d)) for d in (0, 1, 2)]; Tq = [bc(tf.tphi(q, T, d)) for d in (0, 1, 2)]
        return (pb, mu.reshape(pb.shape), T, Tk[0], Tk[1], Tk[2], Tq[0], Tq[1], Tq[2], P0, P1, P2, H0, H1, H2, H3,
                tP[0], tP[1], tP[2], tQ[0], tQ[1], td0)
    def ev(nm, where):
        a = args(where)
        pb, mub, T = a[0], a[1], a[2]
        return np.broadcast_to(tf.KER3[nm](k, q, pb, mub, T, *a[3:]), (B,) + np.shape(T))
    S_sub = ev('SPhi', 's')
    C3v = tf.C3(k, q, mu)
    if C3shift is not None: C3v = C3v + C3shift
    Phi = g.green('Phi', np.full(B, k), S_sub) + C3v[:, None] * tf.ghat(k, t)[None, :]
    dPhi = g.green('Phi', np.full(B, k), S_sub, 1) + C3v[:, None] * tf.ghat(k, t, 1)[None, :]
    Q = ev('Q', 'n'); S00 = ev('S00', 'n'); S0z = ev('S0z', 'n')
    Psi = Phi - Q
    dPsi = dPhi - spectral_d(g, Q)
    delta = t**2 / 3 * S00 - 2 / 3 * t**2 * k**2 * Phi - 2 * t * dPhi - 2 * Psi
    V = -(t**2 / 4) * (2 * (dPhi + Psi / t) + S0z / (1j * k))
    out = dict(Phi=Phi, dPhi=dPhi, Psi=Psi, dPsi=dPsi, delta=delta, V=V, dV=spectral_d(g, V), Q=Q, S0z=S0z, S00=S00, C3=C3v)
    return out

if __name__ == '__main__':
    import time
    g = TGrid()
    # symmetry of the second-order scalar kernels used as KPhi2/Kh (AB/2 = K_X requires K_X(q1,q2) = K_X(q2,q1))
    for nm in ('K_SPhi', 'K_Phi2_0', 'K_h', 'K_Q', 'K_S00', 'K_S0z'):
        a = tf.KER2[nm](1.3, 0.7, 1.1, 2.3); b = tf.KER2[nm](1.3, 1.1, 0.7, 2.3)
        print('symmetry %-9s' % nm, abs(a - b))
    t0 = time.time()
    k = q = 1.0
    eps = np.array([1e-1, 1e-2, 1e-3, 1e-4])
    for side in (+1, -1):
        mu = side * (1 - eps)
        X = X3c(g, k, q, mu)
        i8 = -1
        print('side mu -> %+d: (1-|mu|) =' % side, eps)
        for nm in ('Phi', 'Psi', 'V', 'delta', 'C3'):
            v = X[nm][:, i8] if X[nm].ndim == 2 else X[nm]
            print('   %-6s' % nm, np.round(v, 6), '  x(1-|mu|):', np.round(v * eps, 6))
    print('time', time.time() - t0)
