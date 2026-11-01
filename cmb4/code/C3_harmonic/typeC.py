"""
typeC.py -- second-order line-of-sight maps with a convolution kernel (type C), their C^(22), and the C x B cross term.

Theta_C(n) = int d^3k1 d^3k2/(2pi)^6 zeta_{k1} zeta_{k2} int d eta K(k1,k2,K,eta) Q(n; k1hat,k2hat,Khat) e^{-i K.n chi},  K = k1+k2,
   Q = sum_{abc} q_{abc}(k1,k2,mu) (n.k1hat)^a (n.k2hat)^b (n.Khat)^c ,   mu = k1hat.k2hat ,  degree d = max(a+b+c).
   (Transverse / TT projectors contracted with n are polynomials of this form, e.g. n^i P_ij(Khat) k1hat^j = (n.k1hat) - (n.Khat)(Khat.k1hat).)

Kernel in harmonic space:  a_lm = int int zeta zeta K_{2,lm}(k1,k2),   K_{2,lm}(k1,k2) := int dOmega_n Y*_lm(n) int d eta K Q e^{-iK.n chi}.
General frame (Rayleigh + repeated angular coupling):
   K_{2,lm}(k1,k2) = int d eta K(k1,k2,K,eta) sum_{L=|l-d|}^{l+d} sum_M 4 pi (-i)^L j_L(K chi) Y*_{LM}(Khat) < lm | Q | LM >,
   < lm | Q | LM > := int dOmega Y*_lm Q Y_LM  -- a FINITE sum of 3j symbols obtained by applying the coupling rule
   (n.vhat) Y_{lm}(n) = (4pi/3) sum_{mu} Y*_{1mu}(vhat) sum_{l'=l+-1, m'} G^{l' l 1}_{m' m mu} Y_{l'm'}(n)   (a+b+c) times.
   Only |L - l| <= d contributes (triangle rule), so the L-sum is finite for polynomial Q.
Frame Khat = z, k1,k2 in the x-z plane (cos a1 = (k1+k2 mu)/K, sin a1 = k2 sqrt(1-mu^2)/K; cos a2 = (k2+k1 mu)/K, sin a2 = -k1 sqrt(1-mu^2)/K):
   Ktilde_lm(k1,k2,mu) = int d eta K sum_L (-i)^L sqrt(4pi(2L+1)) j_L(K chi) < lm | Q | L0 >,   |m| <= d.
Isotropy:  sum_m |K_{2,lm}(k1,k2)|^2 is invariant under joint rotations => a function of (k1,k2,mu) only =>
   int dOmega_1 dOmega_2 (.) = 8 pi^2 int_{-1}^{1} dmu (.),  and with SPEC §8 (K_2 symmetrised in k1<->k2):

   C^(22)_l = (2/(2l+1)) sum_m int d^3k1 d^3k2/(2pi)^6 P P |K_{2,lm}|^2
            = 1/(4 pi^4 (2l+1)) int k1^2 dk1 k2^2 dk2 P(k1) P(k2) int_{-1}^{1} dmu sum_m |Ktilde_lm(k1,k2,mu)|^2 .            (C1)

   If the kernel handed over is NOT symmetric, use |K|^2 -> K_lm(k1,k2) [K_lm(k1,k2) + K_lm(k2,k1)]^* / 2 ... i.e. symmetrise first.
Cross term with a type-B map B = A A' (symmetrised kernel K^B(k1,k2;n) = (1/2)[K^A(k1;n)K^{A'}(k2;n) + (1<->2)]):
   (1/(2l+1)) sum_m <C_lm B*_lm> = 1/(4 pi^4 (2l+1)) int k1^2 dk1 k2^2 dk2 P P int dmu sum_m Ktilde^C_lm Ktilde^{B*}_lm,      (C2)
   Ktilde^B_lm(k1,k2,mu) = (1/2) sum_{l1 m1 l2 m2} G^{l l1 l2}_{m m1 m2} [ Atil_{l1m1}(k1) A'til_{l2m2}(k2) + A'til_{l1m1}(k1) Atil_{l2m2}(k2) ],
   Atil_{l1m1}(k) := 4 pi (-i)^{l1} Y*_{l1m1}(khat) Delta^A_{l1}(k)  (khat at the frame angles; the l1,l2 sums converge by Bessel decay).
Monochromatic P = 2 pi^2 sigma^2 delta(k-k0)/k0^2  (so <zeta^2> = sigma^2):
   C^(22)_l = sigma^4/(2l+1) int_{-1}^{1} dmu sum_m |Ktilde_lm(k0,k0,mu)|^2 ,   K = k0 sqrt(2(1+mu)).                           (C1-mono)
"""
import numpy as np
from harm_core import gaunt, Y, jl, w3j000, SHT, cl_from_alm
from typeB import Cl_product, w_plain


# ------------------------------------------------------------------ harmonic polynomials: multiply by (n.vhat)
def mult_dot(coeffs, vhat):
    """coeffs: dict {(l,m): c} for f(n) = sum c Y_lm(n); returns dict for (n.vhat) f(n)."""
    th, ph = np.arccos(np.clip(vhat[2], -1, 1)), np.arctan2(vhat[1], vhat[0])
    y1 = {mu: np.conj(Y(1, mu, th, ph)) for mu in (-1, 0, 1)}
    out = {}
    for (l, m), c in coeffs.items():
        if c == 0:
            continue
        for mu in (-1, 0, 1):
            mp = m + mu
            for lp in (l - 1, l + 1):
                if lp < 0 or abs(mp) > lp:
                    continue
                g = gaunt(lp, mp, l, m, 1, mu)
                if g == 0:
                    continue
                out[(lp, mp)] = out.get((lp, mp), 0.0) + (4 * np.pi / 3) * y1[mu] * g * c
    return out


def Q_times_YLM(Qterms, L, M, k1hat, k2hat, Khat):
    """harmonic coefficients of Q(n) Y_LM(n);  Qterms = list of (coef, a, b, c)."""
    total = {}
    for coef, a, b, c in Qterms:
        if coef == 0:
            continue
        f = {(L, M): 1.0}
        for _ in range(a):
            f = mult_dot(f, k1hat)
        for _ in range(b):
            f = mult_dot(f, k2hat)
        for _ in range(c):
            f = mult_dot(f, Khat)
        for key, v in f.items():
            total[key] = total.get(key, 0.0) + coef * v
    return total


def K2_lm_general(l, k1, k2, Kfun, Qfun, chi_nodes, chi_w, mlist=None):
    """K_{2,lm}(k1,k2) (general frame).  k1,k2: 3-vectors.  Kfun(k1,k2,K,chi) -> kernel value at that chi
    (chi_nodes/chi_w: quadrature of int d eta -> sum_w; pass chi_nodes=[chi], chi_w=[1] for a single-time map or a delta term).
    Qfun(k1,k2,mu) -> list of (coef,a,b,c).  Returns dict m -> K_lm."""
    Kv = k1 + k2
    Kn = np.linalg.norm(Kv)
    k1n, k2n = np.linalg.norm(k1), np.linalg.norm(k2)
    mu = k1 @ k2 / (k1n * k2n)
    Qterms = Qfun(k1n, k2n, mu)
    d = max(a + b + c for _, a, b, c in Qterms)
    mlist = range(-l, l + 1) if mlist is None else mlist
    out = {m: 0.0 for m in mlist}
    if Kn < 1e-14:
        Khat = np.array([0, 0, 1.0]); thK, phK = 0.0, 0.0
    else:
        Khat = Kv / Kn; thK, phK = np.arccos(np.clip(Khat[2], -1, 1)), np.arctan2(Khat[1], Khat[0])
    k1hat, k2hat = k1 / k1n, k2 / k2n
    for L in range(max(0, l - d), l + d + 1):
        radial = sum(w * Kfun(k1n, k2n, Kn, chi) * jl(L, Kn * chi) for chi, w in zip(chi_nodes, chi_w))
        if radial == 0:
            continue
        for M in range(-L, L + 1):
            ylm = np.conj(Y(L, M, thK, phK))
            if abs(ylm) < 1e-300:
                continue
            f = Q_times_YLM(Qterms, L, M, k1hat, k2hat, Khat)
            for m in mlist:
                out[m] += 4 * np.pi * (-1j) ** L * radial * ylm * f.get((l, m), 0.0)
    return out


def frame_vectors(k1, k2, mu):
    """Khat = z frame, k1,k2 in the x-z plane"""
    K = np.sqrt(k1 ** 2 + k2 ** 2 + 2 * k1 * k2 * mu)
    s = np.sqrt(max(1 - mu ** 2, 0.0))
    v1 = np.array([k2 * s, 0.0, k1 + k2 * mu]) * (k1 / K)
    v2 = np.array([-k1 * s, 0.0, k2 + k1 * mu]) * (k2 / K)
    return v1, v2


def K2_lm_frame(l, k1, k2, mu, Kfun, Qfun, chi_nodes, chi_w):
    v1, v2 = frame_vectors(k1, k2, mu)
    return K2_lm_general(l, v1, v2, Kfun, Qfun, chi_nodes, chi_w)


def C22_typeC_mono(l, k0, sigma2, Kfun, Qfun, chi_nodes, chi_w, nmu=64):
    """(C1-mono): sigma^4/(2l+1) int dmu sum_m |Ktilde_lm(k0,k0,mu)|^2 ; Gauss-Legendre in mu."""
    mus, ws = np.polynomial.legendre.leggauss(nmu)
    tot = 0.0
    for mu, w in zip(mus, ws):
        K = K2_lm_frame(l, k0, k0, mu, Kfun, Qfun, chi_nodes, chi_w)
        tot += w * sum(abs(v) ** 2 for v in K.values())
    return sigma2 ** 2 / (2 * l + 1) * tot


def KB_lm_frame(l, k1, k2, mu, DeltaA, DeltaAp, lmax_in):
    """type-B kernel in the Khat=z frame (symmetrised), DeltaA(l,k) transfer functions"""
    v1, v2 = frame_vectors(k1, k2, mu)
    def atil(Delta, v, l1, m1):
        vn = np.linalg.norm(v); vh = v / vn
        return 4 * np.pi * (-1j) ** l1 * np.conj(Y(l1, m1, np.arccos(np.clip(vh[2], -1, 1)), np.arctan2(vh[1], vh[0]))) * Delta(l1, vn)
    out = {m: 0.0 for m in range(-l, l + 1)}
    for l1 in range(lmax_in + 1):
        for l2 in range(abs(l - l1), min(l + l1, lmax_in) + 1):
            if (l + l1 + l2) % 2:
                continue
            for m1 in range(-l1, l1 + 1):
                for m in out:
                    m2 = m - m1
                    if abs(m2) > l2:
                        continue
                    g = gaunt(l, m, l1, m1, l2, m2)
                    if g == 0:
                        continue
                    out[m] += 0.5 * g * (atil(DeltaA, v1, l1, m1) * atil(DeltaAp, v2, l2, m2)
                                         + atil(DeltaAp, v1, l1, m1) * atil(DeltaA, v2, l2, m2))
    return out


# ------------------------------------------------------------------ tests
def _toyK(k1, k2, K, chi):
    return np.exp(-0.1 * K ** 2) * (1 + 0.3 * k1 * k2) / (1 + chi)


def _toyQ(k1, k2, mu):
    # 1 + 0.7 (n.k1)(n.k2) - 0.4 (n.K)^2 + 0.2 mu (n.k1)^2 (n.K)
    return [(1.0, 0, 0, 0), (0.7, 1, 1, 0), (-0.4, 0, 0, 2), (0.2 * mu, 2, 0, 1)]


def test_typeC_matrix_elements(lmax=4, verbose=True):
    """finite-sum K_{2,lm} vs direct 2D quadrature over n (random vectors, two chi values)."""
    rng = np.random.default_rng(4)
    sht = SHT(lmax + 8, ntheta=48, nphi=64)
    n = sht.n_vectors()
    out = []
    chi_nodes, chi_w = [0.9, 2.4], [0.6, 0.4]
    for trial in range(2):
        k1 = rng.normal(size=3); k2 = rng.normal(size=3)
        K = k1 + k2; Kn = np.linalg.norm(K)
        k1n, k2n = np.linalg.norm(k1), np.linalg.norm(k2); mu = k1 @ k2 / (k1n * k2n)
        Qt = _toyQ(k1n, k2n, mu)
        Qn = sum(c * (n @ (k1 / k1n)) ** a * (n @ (k2 / k2n)) ** b * (n @ (K / Kn)) ** cc for c, a, b, cc in Qt)
        f = sum(w * _toyK(k1n, k2n, Kn, chi) * np.exp(-1j * (n @ K) * chi) for chi, w in zip(chi_nodes, chi_w)) * Qn
        alm = sht.map2alm(f, lmax)
        err = 0.0
        for l in range(lmax + 1):
            Kl = K2_lm_general(l, k1, k2, _toyK, _toyQ, chi_nodes, chi_w)
            for m in range(-l, l + 1):
                err = max(err, abs(Kl[m] - alm[(l, m)]))
        out.append(("type-C K_{2,lm} finite 3j-sum vs direct n-quadrature (trial %d)" % trial, err))
        # isotropy: sum_m |K_lm|^2 equals the Khat=z frame value
        err2 = 0.0
        for l in range(1, lmax + 1):
            Kl = K2_lm_general(l, k1, k2, _toyK, _toyQ, chi_nodes, chi_w)
            Kf = K2_lm_frame(l, k1n, k2n, mu, _toyK, _toyQ, chi_nodes, chi_w)
            s1 = sum(abs(v) ** 2 for v in Kl.values()); s2 = sum(abs(v) ** 2 for v in Kf.values())
            err2 = max(err2, abs(s1 - s2) / s2)
        out.append(("sum_m |K_lm|^2 : random frame vs Khat=z frame (isotropy reduction)", err2))
    if verbose:
        for nme, e in out:
            print("  %-80s max err = %.2e" % (nme, e))
    return out


def test_typeC_vs_typeB(lmax=4, verbose=True):
    """K = f(k1) g(k2), Q = 1, single chi: the type-C map is the product A A' of two type-A maps.
    (C1-mono) must reproduce (B1) with monochromatic spectra C^{XY}_l = 4 pi sigma^2 Delta^X_l Delta^Y_l.
    Also tests the C x B cross formula (C2) [C == B here, so cross == auto]."""
    k0, chi, sigma2 = 1.0, 3.0, 1.0
    f = lambda k: 1.0 + 0.5 * k
    g = lambda k: 2.0 - 0.3 * k
    Kfun = lambda k1, k2, K, c: 0.5 * (f(k1) * g(k2) + g(k1) * f(k2))   # symmetrised
    Qfun = lambda k1, k2, mu: [(1.0, 0, 0, 0)]
    lmax_in = 24
    DA = lambda l, k: f(k) * jl(l, k * chi)
    DAp = lambda l, k: g(k) * jl(l, k * chi)
    ls = np.arange(lmax_in + 1)
    CAA = 4 * np.pi * sigma2 * np.array([DA(l, k0) ** 2 for l in ls])
    CBB = 4 * np.pi * sigma2 * np.array([DAp(l, k0) ** 2 for l in ls])
    CAB = 4 * np.pi * sigma2 * np.array([DA(l, k0) * DAp(l, k0) for l in ls])
    typeB = Cl_product(CAA, CBB, CAB, CAB, lmax, lmax_in, w=w_plain)
    out = []
    errC, errX = 0.0, 0.0
    for l in range(1, lmax + 1):
        c22 = C22_typeC_mono(l, k0, sigma2, Kfun, Qfun, [chi], [1.0], nmu=80)
        errC = max(errC, abs(c22 - typeB[l]) / typeB[l])
        # cross (C2) with the type-B frame kernel
        mus, ws = np.polynomial.legendre.leggauss(40)
        tot = 0.0
        for mu, w in zip(mus, ws):
            KC = K2_lm_frame(l, k0, k0, mu, Kfun, Qfun, [chi], [1.0])
            KB = KB_lm_frame(l, k0, k0, mu, DA, DAp, lmax_in)
            tot += w * sum(KC[m] * np.conj(KB[m]) for m in KC)
        cross = sigma2 ** 2 / (2 * l + 1) * tot
        errX = max(errX, abs(cross - typeB[l]) / typeB[l])
    out.append(("type-C (C1-mono) vs type-B (B1) for K=f g, Q=1 (Bessel addition identity), l=1..%d" % lmax, errC))
    out.append(("C x B cross formula (C2) vs (B1)  (same toy)", abs(errX)))
    if verbose:
        for nme, e in out:
            print("  %-80s max rel err = %.2e" % (nme, e))
    return out


if __name__ == "__main__":
    print("type-C tests")
    test_typeC_matrix_elements()
    test_typeC_vs_typeB()
