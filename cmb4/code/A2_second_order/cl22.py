"""cl22.py -- C_l^(22) for the MONOCHROMATIC spectrum P(k) = 2 pi^2 sigma^2 delta(k - k0)/k0^2 (results in units sigma^4).

Exact reduction (SPEC §8, C3 (C1) + mono rules):  with K_{2,lm}(k1,k2) = int dOmega_n Y*_lm(n) K_2(k1,k2;n),
   C^(22)_l = (2/(2l+1)) sum_m int d^3k1 d^3k2/(2pi)^6 P P |K_{2,lm}|^2          [Wick factor 2: two pairings, K_2 symmetric]
            = (sigma^4/(2l+1)) int_{-1}^{1} dmu  sum_m |K_{2,lm}(k0 khat1, k0 khat2)|^2 ,  mu = khat1.khat2,
   and C_l contains (1/4) C^(22) at O(eps^4).  The m-sum is rotation invariant, so for each mu we use the frame Khat = z,
   k1,2 = k0 (+-sin a, 0, cos a), cos a = K/(2 k0), K = k0 sqrt(2(1+mu)).  Integration variable K in (0, 2k0) (dmu = K dK/k0^2),
   Gauss-Legendre in K: the integrand is analytic in K (not in mu) near mu = -1, and bounded there for l >= 2 (k0_limits.py).
   K_{2,lm} = sky quadrature (GL nth x 2nth) of the engine kernel; the full map kernel is split into channels and the
   channel Gram matrix G_cc'(l) = (1/(2l+1)) int dmu sum_m K^c_lm K^c'*_lm is accumulated, so that C^(22)_l = sum_cc' G_cc'.
Usage: python3 cl22.py k0 eta_s eta_o nK nth lmax tag"""
import sys, os, pickle, time
import numpy as np
from numpy.polynomial import legendre as Lg
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kernel2 import Geometry, MapKernels, Delta_l
from sphere import Sky

def group_of(name):
    if name in ('L1*L1', '-2<m1>m1'): return 'B: L1^2 - 2<m1>m1'
    if '2nd:' in name:
        f = name.split('2nd:')[1]
        return {'Phi': 'C: Phi2,Psi2 (ISW+SW)', 'Psi': 'C: Phi2,Psi2 (ISW+SW)', 'dT': 'C: dT2 (intrinsic)', 'V2': 'C: V2 (Doppler, long.)',
                'v2T': 'C: om2+v2T (vector)', 'om': 'C: om2+v2T (vector)', 'h': 'C: h2 (tensor)'}[f]
    if 'local | fields:1+1' in name:
        return 'C: ray local 1x1' if name.startswith('[ray]') else 'B: endpoint local 1x1'
    if 'lensing' in name: return 'D: lensing (transverse Dx1)'
    if 'time delay' in name: return 'D: time delay (radial Dx1)'
    if 'transverse direction' in name: return 'D: aberration+deflection (k1_perp)'
    if 'emission-surface' in name: return 'D: emission-surface shift'
    if 'clock' in name: return 'D: observer clock'
    if 'displacement' in name: return 'D: observer displacement'
    raise ValueError(name)

def run(k0=1.0, eta_s=2.0, eta_o=8.0, nK=40, nth=48, lmax=24, tag=None, verbose=True):
    t0 = time.time()
    g = Geometry(eta_s, eta_o)
    sky = Sky(nth)
    LM = [(l, m) for l in range(lmax + 1) for m in range(-l, l + 1)]
    Ymat = np.array([np.conj(sky.Y(l, m)) * sky.w for (l, m) in LM])          # (nlm, Npix)
    lidx = np.array([l for l, m in LM])
    xK, wK = Lg.leggauss(nK)
    Ks = k0 * (xK + 1); wKs = k0 * wK                                          # K in (0, 2k0)
    groups = None; G = None; Ktot = None
    lim = {}
    for j, (K, wk) in enumerate(zip(Ks, wKs)):
        wmu = wk * K / k0**2                                                   # dmu = K dK / k0^2
        ca = K / (2 * k0); sa = np.sqrt(1 - ca**2)
        k1 = k0 * np.array([sa, 0, ca]); k2 = k0 * np.array([-sa, 0, ca])
        _, ch = MapKernels(g, k1, k2, sky.nv, form='direct').K2()
        if groups is None:
            groups = sorted(set(group_of(c) for c in ch))
            G = np.zeros((len(groups), len(groups), lmax + 1))
        A = np.zeros((len(groups), len(LM)), complex)
        for c, v in ch.items():
            A[groups.index(group_of(c))] += Ymat @ v
        # Gram: sum_m A_g A_g'^*
        P = np.einsum('gi,hi->ghi', A, np.conj(A)).real                        # (ng, ng, nlm)
        for l in range(lmax + 1):
            G[:, :, l] += wmu * P[:, :, lidx == l].sum(axis=2)
        tot = A.sum(axis=0)
        lim[K] = np.array([np.sum(np.abs(tot[lidx == l])**2) for l in range(lmax + 1)])
        if verbose and (j % 8 == 0):
            print('  node %d/%d K=%.4f  (%.0fs)' % (j + 1, nK, K, time.time() - t0), flush=True)
    ls = np.arange(lmax + 1)
    G /= (2 * ls + 1)[None, None, :]
    C22 = G.sum(axis=(0, 1))
    C11 = np.array([4 * np.pi * Delta_l(l, k0, g)**2 for l in range(lmax + 1)])
    res = dict(k0=k0, eta_s=eta_s, eta_o=eta_o, nK=nK, nth=nth, lmax=lmax, groups=groups, G=G, C22=C22, C11=C11,
               integrand=lim, Delta=[Delta_l(l, k0, g) for l in range(lmax + 1)])
    if tag:
        pickle.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cl22_%s.pkl' % tag), 'wb'))
    if verbose: print('done %.0fs' % (time.time() - t0))
    return res

if __name__ == '__main__':
    a = sys.argv[1:]
    k0, es, eo = float(a[0]), float(a[1]), float(a[2])
    nK, nth, lmax = int(a[3]), int(a[4]), int(a[5]); tag = a[6]
    r = run(k0, es, eo, nK, nth, lmax, tag)
    print(' l    C11/sigma^2        C22/sigma^4        (1/4)C22/sigma^4')
    for l in range(1, lmax + 1):
        print('%3d  %.10e  %.10e  %.10e' % (l, r['C11'][l], r['C22'][l], r['C22'][l] / 4))
