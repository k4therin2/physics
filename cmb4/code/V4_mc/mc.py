"""mc.py -- V4: Monte-Carlo check of A2's monochromatic C_l^(22) (and C_l^(11)) via Gaussian realisations.

Field: zeta(x) = sum_{j=1}^N (sigma/sqrt N)[a_j cos(k_j.x) + b_j sin(k_j.x)] = sum_j [c_j e^{i k_j.x} + c.c.],
       c_j = (sigma/sqrt N)(a_j - i b_j)/2,  |k_j| = k0, isotropic random directions, a_j,b_j ~ N(0,1).
With SPEC Fourier convention zeta_k = (2pi)^3 sum_q c_q delta^3(k-q), q in {+-k_j}, c_{-k_j} = c_j^*, so
   Th1(n) = sum_q c_q K_1(q;n),  Th2(n) = sum_{q,q'} c_q c_q' K_2(q,q';n)    (all ordered pairs, incl. q'=q and q'=-q)
          = sum_{i,j} 2 Re[ c_i c_j K_2(k_i,k_j) + c_i c_j^* K_2(k_i,-k_j) ]     (K_2 symmetric, K_2(-a,-b) = K_2(a,b)^*).
A2's K_2 already contains the normalisation term -Delta_0(k1)K_1(k2) - Delta_0(k2)K_1(k1), whose pair sum is exactly
-2<m1>m1 with <m1> = sum_q c_q Delta_0(q) the realisation's sky monopole (checked numerically below).
Per direction set: the real-harmonic coefficients of every pair kernel are computed once (sky quadrature GL nth x 2nth),
then R Gaussian amplitude draws are made.  Estimators:
   full : Th2 as above (the physical map of the finite-mode field); E = C22 (1 + O(1/N)) [j=j' terms: k2 = +-k1 not isotropic]
   off  : Th2 without the i=j terms; E = (N-1)/N * C22 exactly (pairs of distinct modes have independent isotropic
          directions; the q1=-q2 pairing needs i=j and is absent) -> reported x N/(N-1).
Usage: python3 mc.py N nsets R nth seed nproc tag"""
import sys, os, time, pickle
import numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'A2_second_order'))
from kernel2 import Geometry, MapKernels, Delta_l
from sphere import Sky

K0, ETA_S, ETA_O, LMAX = 1.0, 2.0, 8.0, 12
_G = _SKY = _YR = None

def real_Y(sky, lmax):
    """real orthonormal harmonics (rows lm, l<=lmax) times quadrature weights -> projection matrix."""
    rows, ls = [], []
    for l in range(lmax + 1):
        for m in range(-l, l + 1):
            if m == 0: y = sky.Y(l, 0).real
            elif m > 0: y = np.sqrt(2) * (-1)**m * sky.Y(l, m).real
            else: y = np.sqrt(2) * (-1)**m * sky.Y(l, -m).imag
            rows.append(y * sky.w); ls.append(l)
    return np.array(rows), np.array(ls)

def init(nth):
    global _G, _SKY, _YR
    _G = Geometry(ETA_S, ETA_O); _SKY = Sky(nth); _YR, _ = real_Y(_SKY, LMAX)

def job(args):
    kind, k1, k2 = args
    mk = MapKernels(_G, k1, k2, _SKY.nv, form='direct')
    if kind == 'K1':
        return _YR @ mk.K1(mk.sa)
    K2, ch = mk.K2()
    return _YR @ K2, _YR @ ch['-2<m1>m1']

def direction_set(rng, N):
    v = rng.normal(size=(N, 3)); return K0 * v / np.linalg.norm(v, axis=1)[:, None]

def run(N, nsets, R, nth, seed, nproc, tag):
    rng = np.random.default_rng(seed)
    _, ls = real_Y(Sky(nth), LMAX)
    D0 = Delta_l(0, K0, Geometry(ETA_S, ETA_O))
    out = dict(N=N, R=R, nth=nth, seed=seed, sets=[])
    pool = Pool(nproc, initializer=init, initargs=(nth,))
    t0 = time.time()
    iu = [(i, j) for i in range(N) for j in range(i, N)]
    for s in range(nsets):
        kv = direction_set(rng, N)
        jobs = [('K1', kv[i], kv[i]) for i in range(N)]
        jobs += [('A', kv[i], kv[j]) for i, j in iu] + [('B', kv[i], -kv[j]) for i, j in iu]
        res = pool.map(job, jobs, chunksize=4)
        K1 = np.array(res[:N])                                   # (N, nlm) complex
        A = np.zeros((N, N, len(ls)), complex); B = np.zeros_like(A); An = np.zeros_like(A); Bn = np.zeros_like(A)
        for n_, (i, j) in enumerate(iu):
            A[i, j], An[i, j] = res[N + n_]; B[i, j], Bn[i, j] = res[N + len(iu) + n_]
            if i != j:
                A[j, i], An[j, i] = A[i, j], An[i, j]; B[j, i], Bn[j, i] = np.conj(B[i, j]), np.conj(Bn[i, j])
        # amplitude draws
        a = rng.normal(size=(R, N)); b = rng.normal(size=(R, N))
        c = (a - 1j * b) / (2 * np.sqrt(N))                     # sigma = 1
        T1 = 2 * np.einsum('ri,il->rl', c, K1).real
        T2 = 2 * (np.einsum('ri,rj,ijl->rl', c, c, A) + np.einsum('ri,rj,ijl->rl', c, np.conj(c), B)).real
        dg = np.arange(N)
        T2d = 2 * (np.einsum('ri,il->rl', c * c, A[dg, dg]) + np.einsum('ri,il->rl', np.abs(c)**2, B[dg, dg])).real
        T2off = T2 - T2d
        # normalisation-term check: pair sum of the -2<m1>m1 channel vs -2 <m1> Th1 with <m1> = sum_q c_q Delta_0
        Tn = 2 * (np.einsum('ri,rj,ijl->rl', c, c, An) + np.einsum('ri,rj,ijl->rl', c, np.conj(c), Bn)).real
        m1mono = 2 * (c.sum(axis=1)).real * D0
        normchk = np.max(np.abs(Tn - (-2 * m1mono[:, None] * T1))) / np.max(np.abs(Tn))
        zeta0 = 2 * c.sum(axis=1).real                           # zeta(0): <zeta^2> = sigma^2 check
        pw = lambda T: np.array([[np.sum(T[r, ls == l]**2) / (2 * l + 1) for l in range(LMAX + 1)] for r in range(len(T))])
        # exact amplitude average (Wick, conditional on the directions) of the 'off' estimator, x N/(N-1) applied in analyze
        W = (np.abs(A)**2 + np.abs(B)**2).sum(axis=(0, 1)) - (np.abs(A[dg, dg])**2 + np.abs(B[dg, dg])**2).sum(axis=0)
        W = 4 * W / (2 * N)**2                                   # ordered q,q' over +-: 2|A|^2+2|B|^2 per (i,j); Wick factor 2
        wick = np.array([W[ls == l].sum() / (2 * l + 1) for l in range(LMAX + 1)])
        st = dict(kv=kv, wickoff=wick, C11=pw(T1), C22=pw(T2), C22off=pw(T2off), mean2=T2.mean(0), mean2off=T2off.mean(0),
                  zeta0sq=np.mean(zeta0**2), normchk=normchk)
        out['sets'].append(st)
        pickle.dump(out, open(os.path.join(HERE, 'mc_%s.pkl' % tag), 'wb'))
        print('set %d/%d  %.0fs  normchk %.1e  <zeta^2> %.3f  C22[2..6] %s' % (s + 1, nsets, time.time() - t0, normchk,
              st['zeta0sq'], np.round(st['C22'].mean(0)[2:7] / 4, 3)), flush=True)
    pool.close()

if __name__ == '__main__':
    a = sys.argv[1:]
    run(int(a[0]), int(a[1]), int(a[2]), int(a[3]), int(a[4]), int(a[5]), a[6])
