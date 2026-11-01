"""
thirteen.py -- the (13) term: <a^(1)_lm a^(3)*_lm> for Theta^(3) built from (i) three type-A factors, (ii) type-A x type-C,
(iii) a genuine third-order kernel K3 contracted with P(q); plus the SPEC §7 normalisation pieces.

C^(13)_l := (1/(2l+1)) sum_m <Th1_lm Th3*_lm>  with  Th3 = m3 - 3<m1> m2 - 3<m2> m1 + 6<m1>^2 m1   (SPEC §7)   =>
   C^(13)_l = T3_l - 3 T_{<1>2,l} - 3 T_{<2>1,l} + 6 <<m1>^2> C^(11)_l                                                   (E0)
   <<m1>^2> = (1/2pi^2) int k^2 dk P Delta_0(k)^2  (Delta_0 = l=0 transfer function of m1).

RULE A-L (generalised type-A rule; verified in test_rule_AL):  a kernel  h(k,eta) P_L(n.khat) e^{-i k.n chi}  gives
   Delta_l(k) = int d eta h(k,eta) sum_{Lam=|l-L|}^{l+L} (2Lam+1) (l L Lam;000)^2 (-i)^{Lam-l} j_Lam(k chi)         (A-L)
   (L=0: j_l; L=1: i j_l'; reproduces rule A-s after converting monomials to Legendre polynomials).
   Any kernel that is a function of (n.khat) only can be written this way, which is what every q-contraction produces.

KEY q-ANGULAR REDUCTION (used by (ii),(iii) and by the <m1>-type terms with a type-C m2):
   int dOmega_q/(4pi) g(qhat.khat) (n.qhat)^b e^{-i q.n Dchi}
     = sum_{L'} (-i)^{L'} (2L'+1) j_{L'}(q Dchi) sum_{L''} beta_{L' b L''} g_{L''} P_{L''}(n.khat),                        (E1)
   g_{L''} := (1/2) int_{-1}^{1} dmu P_{L''}(mu) g(mu),   x^b P_{L'}(x) = sum_{L''} beta_{L' b L''} P_{L''}(x)   (finite),
   because int dOmega_q/(4pi) P_{L''}(qhat.n) g(qhat.khat) = g_{L''} P_{L''}(n.khat)  (addition theorem).
   For Dchi = 0 only L' = 0 survives: a FINITE Legendre polynomial of degree b in n.khat with coefficients = Legendre moments of g.
   A remaining factor (n.khat)^a is absorbed by re-expanding nu^a P_{L''}(nu) in Legendre polynomials.  Then apply (A-L).

(iii) genuine third-order field   Theta^(3) ⊃ int zeta zeta zeta K3(k1,k2,k3,eta) Q3(n;...) e^{-i(k1+k2+k3).n chi}  (K3 symmetric):
   contracted kernel  K3c(k,eta;n) = int d^3q/(2pi)^3 P(q) K3(k,q,-q,eta) Q3(n;khat,qhat,-qhat) e^{-i k.n chi}
   = [ int q^2 dq/(2pi^2) P(q) sum_{L} h_L(k,q,eta) P_L(n.khat) ] e^{-i k.n chi}      (E1 with Dchi = 0; depends on |k| only)
   => Delta^{(3c)}_l(k) = int q^2dq/(2pi^2) P(q) int d eta sum_L h_L sum_Lam (2Lam+1)(l L Lam;000)^2 (-i)^{Lam-l} j_Lam(k chi)
   T3^{(iii)}_l = 3 (2/pi) int k^2 dk P(k) Delta^{m1}_l(k) Delta^{(3c)}_l(k)^*        (factor 3 = three pairings; if K3 is not
   symmetrised replace 3 K3(k,q,-q) -> K3(k,q,-q) + K3(q,k,-q) + K3(q,-q,k)).
   Monochromatic: Delta^{(3c)}_l(k0) = sigma^2 sum_L h_L(k0,k0) [...],  h_L from (1/2) int dmu P_L(mu) K3(k0,k0,mu),
   |k+-q| = k0 sqrt(2(1+-mu)) inside K3;  T3^{(iii)} = 12 pi sigma^2 Delta_l(k0) Delta^{(3c)}_l(k0)^*.
(i) three type-A factors  m3 ⊃ A B C  (with possible gradient contractions):
   T3^{(i)}_l = C^{m1 A}_l xi^{BC}(0) + C^{m1 B}_l xi^{AC}(0) + C^{m1 C}_l xi^{AB}(0),
   xi^{XY}(0) := <X(n)Y(n)> = sum_{l'} (2l'+1)/(4pi) C^{XY}_{l'} = (1/2pi^2) int k^2 dk P(k) int d eta d eta' x(k,eta) y(k,eta')^* j_0(k|chi-chi'|)
   (for angular-gradient contractions grad B.grad C -> xi^{grad B grad C}(0) = sum l'(l'+1)(2l'+1)/(4pi) C^{BC}_{l'};
    for grad B grad C : grad grad A -> -l(l+1) R_{BC} C^{m1 A}_l with R_{BC} = (1/2) xi^{grad B grad C}(0), as in the lensing control).
(ii) type-A x type-C   m3 ⊃ A(n) C(n):
   T3^{(ii)}_l = C^{m1 A}_l <<C>> + 2 (2/pi) int k^2 dk P(k) Delta^{m1}_l(k) Delta^{AC}_l(k)^* ,
   <<C>> = int d^3k/(2pi)^3 P(k) K^C(k,-k;n)  (the K=0 tadpole mean of the type-C map; direction independent),
   Delta^{AC}_l(k): type-A transfer function of the effective kernel int d^3q/(2pi)^3 P(q) K^A(q;n) K^C(k,-q;n)
     = int q^2dq/(2pi^2) P(q) int d eta_A d eta_C g_A(q,eta_A) K(k,q,|k-q|,eta_C)
       x [ qhat-average of Q_A(n,qhat) Q_C(n;khat,-qhat) e^{-i q.n (chi_A - chi_C)} ]  e^{-i k.n chi_C}   -> (E1) then (A-L).
Normalisation pieces (type-C or general m2 with sky-average kernel M2(q1,q2) := K_{2,00}(q1,q2)/sqrt(4pi)):
   T_{<2>1,l} = <<m2>> C^(11)_l + (1/pi^3) int k^2dk k'^2dk' P P Delta_l(k) Delta_l(k')^* int dmu P_l(mu) M2(k khat, -k' khat')^*   (E2)
       <<m2>> = (1/2pi^2) int k^2 dk P(k) M2(k,-k)      [K=0 tadpole of the second-order kernel]
   T_{<1>2,l} = 2 (1/2pi^2) int q^2 dq P(q) Delta_0(q) (2/pi) int k^2 dk P(k) Delta_l(k) Deltabar^{m2}_l(k;q)^*                    (E3)
       Deltabar^{m2}_l(k;q): type-A transfer function (in k) of the qhat-averaged kernel int dOmega_q/(4pi) K2(k,q;n),
       = (1/2) int dmu K_{2,l0}(k zhat, q)/((-i)^l sqrt(4pi(2l+1)))  (frame khat = z, mu = cos theta_q)
       (for a type-B m2 = AA': Deltabar_l(k;q) = (1/2)[Delta^A_l(k) Delta^{A'}_0(q) + Delta^{A'}_l(k) Delta^A_0(q)], recovering typeB (ii)).
   For the (22) term the same objects appear:  (1/4)C^(22) = (1/4)[C^{m2m2} - 4 Re T^{(22)}_{2<1>1} + 4<<m1>^2> C^(11)],
       T^{(22)}_{2<1>1,l} = (1/(2l+1)) sum_m <m2_lm (<m1> m1_lm)^*> = conj of (E3) with m2 in place of the (13)-m2.
   Monochromatic: <<m1>^2> = sigma^2 Delta_0(k0)^2; <<m2>> = sigma^2 M2(k0 khat,-k0 khat);
       (E2) cross piece = 4 pi sigma^4 |Delta_l(k0)|^2 int dmu P_l(mu) M2(k0 khat,-k0 khat')^*  (|k-k'| = k0 sqrt(2(1-mu)) inside M2);
       (E3) = 8 pi sigma^4 Delta_0(k0) Delta_l(k0) Deltabar_l(k0;k0)^*.
"""
import numpy as np
from numpy.polynomial import legendre as L
from scipy.special import eval_legendre
from harm_core import w3j000, jl, SHT, Y, jl_d, monomial_legendre_coeffs


def rule_AL(l, Lleg, x):
    """radial function for P_L(n.khat) e^{-i k.n chi}: sum_Lam (2Lam+1)(l L Lam;000)^2 (-i)^(Lam-l) j_Lam(x)"""
    s = 0.0
    for Lam in range(abs(l - Lleg), l + Lleg + 1):
        c = w3j000(l, Lleg, Lam) ** 2
        if c:
            s = s + (2 * Lam + 1) * c * (-1j) ** (Lam - l) * jl(Lam, x)
    return s


def legendre_moments(g, Lmax, nmu=200):
    """g_L = (1/2) int dmu P_L(mu) g(mu), L = 0..Lmax (Gauss-Legendre)"""
    mus, ws = np.polynomial.legendre.leggauss(nmu)
    gv = np.array([g(m) for m in mus])
    return np.array([0.5 * np.sum(ws * gv * eval_legendre(Lg, mus)) for Lg in range(Lmax + 1)])


def leg_times_monomial(Lp, b):
    """x^b P_Lp(x) = sum_Lpp beta P_Lpp(x);  returns array beta[Lpp] (exact via GL quadrature of a polynomial)"""
    n = Lp + b + 2
    mus, ws = np.polynomial.legendre.leggauss(n)
    f = mus ** b * eval_legendre(Lp, mus)
    return np.array([(2 * Lpp + 1) / 2 * np.sum(ws * f * eval_legendre(Lpp, mus)) for Lpp in range(Lp + b + 1)])


def qavg_legendre(g, b, qDchi, a=0, Lpmax=None, nmu=200):
    """(E1): int dOmega_q/(4pi) g(qhat.khat) (n.qhat)^b e^{-i q.n Dchi}  times (n.khat)^a
    returned as Legendre coefficients h[L] of P_L(n.khat).  Lpmax: truncation of the Rayleigh sum (only L'=0 if qDchi=0)."""
    if qDchi == 0:
        Lpmax = 0
    elif Lpmax is None:
        Lpmax = int(qDchi + 12 + 3 * qDchi ** (1 / 3))
    gL = legendre_moments(g, Lpmax + b, nmu)
    h = np.zeros(Lpmax + b + a + 1, dtype=complex)
    for Lp in range(Lpmax + 1):
        pref = (-1j) ** Lp * (2 * Lp + 1) * jl(Lp, qDchi)
        beta = leg_times_monomial(Lp, b)
        for Lpp, bb in enumerate(beta):
            if bb == 0:
                continue
            coef = pref * bb * gL[Lpp]
            if a == 0:
                h[Lpp] += coef
            else:
                gam = leg_times_monomial(Lpp, a)
                h[:len(gam)] += coef * gam
    return h


def Delta_from_legendre(l, x, h):
    """Delta_l for a kernel sum_L h[L] P_L(n.khat) e^{-i k.n chi}  (x = k chi), by rule (A-L)"""
    return sum(hL * rule_AL(l, Lg, x) for Lg, hL in enumerate(h) if hL != 0)


# ------------------------------------------------------------------ tests
def test_rule_AL(lmax=5, verbose=True):
    rng = np.random.default_rng(9)
    sht = SHT(lmax + 8, ntheta=48, nphi=64)
    n = sht.n_vectors()
    khat = rng.normal(size=3); khat /= np.linalg.norm(khat)
    thk, phk = np.arccos(khat[2]), np.arctan2(khat[1], khat[0])
    x = 4.4
    nu = n @ khat
    err = 0.0
    for Lg in range(0, 4):
        alm = sht.map2alm(eval_legendre(Lg, nu) * np.exp(-1j * x * nu), lmax)
        for l in range(lmax + 1):
            for m in range(-l, l + 1):
                pred = 4 * np.pi * (-1j) ** l * np.conj(Y(l, m, thk, phk)) * rule_AL(l, Lg, x)
                err = max(err, abs(alm[(l, m)] - pred))
    if verbose:
        print("  %-80s max abs err = %.2e" % ("rule (A-L): P_L(n.khat) e^{-ik.n chi} -> sum (2Lam+1)(l L Lam;000)^2 (-i)^(Lam-l) j_Lam", err))
    return err


def test_qavg_reduction(verbose=True):
    """(E1) vs brute-force quadrature over qhat, for Dchi = 0 and Dchi != 0, including an extra (n.khat)^a factor."""
    rng = np.random.default_rng(3)
    g = lambda mu: np.exp(0.8 * mu) / (1.3 + mu)          # toy kernel of qhat.khat (e.g. depends on |k-q|)
    khat = np.array([0.3, -0.5, 0.81]); khat /= np.linalg.norm(khat)
    # quadrature over qhat
    xg, wg = np.polynomial.legendre.leggauss(80)
    phis = 2 * np.pi * np.arange(128) / 128
    TH, PH = np.meshgrid(np.arccos(xg), phis, indexing='ij')
    W = (wg[:, None] * (2 * np.pi / 128)) / (4 * np.pi)
    qh = np.stack([np.sin(TH) * np.cos(PH), np.sin(TH) * np.sin(PH), np.cos(TH)], -1)
    gq = g(qh @ khat)
    out = []
    for (b, a, qD) in ((2, 0, 0.0), (3, 1, 0.0), (1, 0, 2.7), (2, 2, 3.9)):
        h = qavg_legendre(g, b, qD, a=a)
        err = 0.0
        for t in range(5):
            nv = rng.normal(size=3); nv /= np.linalg.norm(nv)
            direct = np.sum(W * gq * (qh @ nv) ** b * np.exp(-1j * (qh @ nv) * qD)) * (nv @ khat) ** a
            pred = sum(hL * eval_legendre(Lg, nv @ khat) for Lg, hL in enumerate(h))
            err = max(err, abs(direct - pred))
        out.append(("q-average reduction (E1): b=%d a=%d qDchi=%.1f" % (b, a, qD), err))
    if verbose:
        for nme, e in out:
            print("  %-80s max abs err = %.2e" % (nme, e))
    return out


if __name__ == "__main__":
    print("(13) reductions")
    test_rule_AL()
    test_qavg_reduction()
