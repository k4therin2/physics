"""
lensing.py -- CONTROL CASE: full-sky lensed CMB power spectrum to O(C^phi) (Hu 2000 eqs 62-63; Lewis-Challinor 2006)
reproduced from the SPEC (22)/(13) bookkeeping, plus the general type-D (nested / unequal-time) reduction.

Remapped map  Theta~(n) = Theta(n + grad phi) = Theta + grad phi . grad Theta + (1/2) grad^i phi grad^j phi grad_i grad_j Theta + ...
Identify with SPEC:  Theta = eps Th1 + eps^2 Th2/2 + eps^3 Th3/6  with Th1 = Theta, phi both O(eps):
      Th2 = 2 grad phi . grad Theta ,      Th3 = 6 * (1/2) grad phi grad phi : grad grad Theta = 3 grad phi grad phi : grad grad Theta.
SPEC:  C_l = eps^2 C^(11) + eps^4 [ 1/4 C^(22) + 1/3 Re C^(13) ],  C^(22) = (1/(2l+1)) sum_m <|Th2_lm|^2>,
       C^(13) = (1/(2l+1)) sum_m <Th1_lm Th3*_lm>.   So  1/4 C^(22) = <|delta1_lm|^2>,  1/3 Re C^(13) = 2 Re <Theta_lm delta2*_lm>,
       which is exactly the (first-order)^2 + 2 x (zeroth x second-order) decomposition of Lewis-Challinor (4.11).  [bookkeeping OK]

(22): a^{delta1}_lm = sum I^{m m1 m2}_{l l1 l2} phi_{l1m1} Theta_{l2m2},  I = (1/2)[l1(l1+1)+l2(l2+1)-l(l+1)] G  (Hu 47,50)
      => (B1) with w = w_grad:
      1/4 C^(22)_l = sum_{l1 l2} (2l1+1)(2l2+1)/(4pi) (l l1 l2;000)^2 w^2 [ C^{phi phi}_{l1} C^{TT}_{l2} + C^{phi T}_{l1} C^{T phi}_{l2} ]
                   = sum_{l1 l2} (F_{l l1 l2})^2/(2l+1) [ C^{phi phi}_{l1} C^{TT}_{l2} + C^{phi T}_{l1} C^{T phi}_{l2} ],
      F_{l l1 l2} = (1/2)[l1(l1+1)+l2(l2+1)-l(l+1)] sqrt((2l+1)(2l1+1)(2l2+1)/(4pi)) (l l1 l2;000)   (Hu 52).  EXACT MATCH with Hu (62)
      for the first term; the second term is the extra piece present when <Theta phi> != 0 (Hu/LC assume it vanishes).
(13): <Theta_lm delta2*_lm>:  Wick on <Theta (1/2) grad^i phi grad^j phi grad_i grad_j Theta>:
      (a) <grad^i phi grad^j phi>(same point) <Theta grad_i grad_j Theta> : 2D isotropy gives <grad_i phi grad_j phi> = g_ij R,
          R := (1/2) <|grad phi|^2> = (1/2) sum_l l(l+1)(2l+1)/(4pi) C^{phi phi}_l ;  grad^2 Y_lm = -l(l+1) Y_lm
          =>  <Theta_lm delta2*_lm> = -(1/2) l(l+1) R C^{TT}_l ,   1/3 Re C^(13) = - l(l+1) R C^{TT}_l .   EXACT MATCH with Hu (63).
      (b) the two pairings <Theta grad phi><grad phi grad grad Theta> involve the coincident-point tensor
          <grad_i grad_j phi grad_k Theta> = 0 (no parity-even isotropic rank-3 tensor on the tangent plane) and
          grad_i <grad^j phi grad_j Theta> = 0 (constant) => they vanish identically even when <Theta phi> != 0.
Total:  C~_l = C_l [1 - l(l+1) R] + sum_{l1l2} F^2/(2l+1) [C^phi_{l1} C_{l2} + C^{phi T}_{l1} C^{T phi}_{l2}]   = Hu (62)-(63) for C^{phiT}=0.

SPEC §8 prefactors re-derived (Theta^(r) = int prod d^3k_i/(2pi)^3 zeta...zeta K_r, K_r symmetric):
   C^(22): <zeta_1 zeta_2 zeta_3* zeta_4*> connected = 2 pairings, each giving int P P |K_{2,lm}|^2 (K symmetric) => factor 2. OK.
   C^(13): <zeta_k zeta*_{q1} zeta*_{q2} zeta*_{q3}> = 3 pairings k<->q_i, remaining pair q_j+q_k = 0 => factor 3,
           C^(13) = (3/(2l+1)) sum_m int d^3k/(2pi)^3 P(k) K_{1,lm}(k) [int d^3q/(2pi)^3 P(q) K_{3,lm}(k,q,-q)]^*.
           SPEC writes K*_1 K_3, i.e. the complex conjugate (= <a^(3) a^(1)*>); harmless since only Re C^(13) enters. OK.
           Note the (13) "tadpole" pairing q_j = -q_k is NOT removed: it is the whole (13) term.

TYPE D (nested, unequal-time):  Theta_D(n) = int d eta1 int_{eta1}^{eta_o} d eta2 W(eta1,eta2) (delta_ij - n_i n_j) d_i X(eta1,-n chi1) d_j Y(eta2,-n chi2)
  Angular-gradient identity:  (delta - nn)_{ij} d_i X d_j Y |_{x=-n chi} = (1/(chi1 chi2)) grad_n X~ . grad_n Y~ ,  X~(n) := X(eta1, -n chi1)
  (shell maps; a_lm of X~ = 4 pi (-i)^l int zeta_k Y*_lm(khat) x(k,eta1) j_l(k chi1)).  Hence
     a^D_lm = int d eta1 d eta2 W/(chi1 chi2) sum w_grad(l,l1,l2) G^{l l1 l2}_{m m1 m2} x~_{l1m1}(eta1) y~_{l2m2}(eta2)
  and (B1) with unequal-time shell spectra  C^{xx}_l(eta,eta') = (2/pi) int k^2 dk P x(k,eta) x(k,eta')^* j_l(k chi) j_l(k chi'):
     C^{DD}_l = sum_{l1 l2} (2l1+1)(2l2+1)/(4pi) (000)^2 w_grad^2 int d eta1 d eta2 d eta1' d eta2' W W'/(chi1 chi2 chi1' chi2')
                [ C^{xx}_{l1}(eta1,eta1') C^{yy}_{l2}(eta2,eta2') + C^{xy}_{l1}(eta1,eta2') C^{yx}_{l2}(eta2,eta1') ].
  Radial contraction n.dX n.dY: n.grad X(-n chi) = -d/dchi X~  => shell a_lm with j_l(k chi) -> -k j_l'(k chi), plain Gaunt (w = 1).
  Fourier form (same object as a two-time type-C map):  Q = -k1 k2 [mu - (n.k1hat)(n.k2hat)],  e^{-i k1.n chi1 - i k2.n chi2};
  the identity between the two routes is checked numerically in test_typeD_identity().
"""
import numpy as np
from harm_core import w3j000, SHT, gaunt, cl_from_alm, cross_cl_from_alm, Y, jl
from typeB import Cl_product, w_grad, w_plain, correlated_realisation, _toy_spectra


def F_hu(l, l1, l2):
    return (0.5 * (l1 * (l1 + 1) + l2 * (l2 + 1) - l * (l + 1))
            * np.sqrt((2 * l + 1) * (2 * l1 + 1) * (2 * l2 + 1) / (4 * np.pi)) * w3j000(l, l1, l2))


def R_lens(Cphi):
    l = np.arange(len(Cphi))
    return 0.5 * np.sum(l * (l + 1) * (2 * l + 1) / (4 * np.pi) * Cphi)


def lensed_cl_hu(CT, Cphi, lmax, lmax_in, CphiT=None):
    """Hu (62)-(63) [+ correlated-term extension]: returns (C~_l, 22-part, 13-part)"""
    CphiT = np.zeros(lmax_in + 1) if CphiT is None else CphiT
    c22 = np.zeros(lmax + 1)
    for l in range(lmax + 1):
        for l1 in range(lmax_in + 1):
            for l2 in range(abs(l - l1), min(l + l1, lmax_in) + 1):
                if (l + l1 + l2) % 2:
                    continue
                c22[l] += F_hu(l, l1, l2) ** 2 / (2 * l + 1) * (Cphi[l1] * CT[l2] + CphiT[l1] * CphiT[l2])
    l = np.arange(lmax + 1)
    c13 = -l * (l + 1) * R_lens(Cphi) * CT[:lmax + 1]
    return CT[:lmax + 1] + c22 + c13, c22, c13


def test_lensing_control(nreal=400, lmax=4, lmax_in=6, verbose=True):
    rng = np.random.default_rng(5)
    CT, Cphi, CphiT = _toy_spectra(lmax_in, rng)
    Cphi = Cphi * 0.05  # weak lensing
    CphiT = CphiT * np.sqrt(0.05)
    out = []
    # (1) algebraic identity: my (B1) with w_grad == Hu's F^2/(2l+1) convolution (correlated case included)
    mine = Cl_product(Cphi, CT, CphiT, CphiT, lmax, lmax_in, w=w_grad)
    hu, c22, c13 = lensed_cl_hu(CT, Cphi, lmax, lmax_in, CphiT)
    out.append(("1/4 C^(22) from (B1,w_grad)  ==  sum F^2/(2l+1) C^phi C^T (+corr.)   [exact]", np.max(np.abs(mine[1:] - c22[1:]) / c22[1:])))
    # (2) Monte Carlo of the Taylor-expanded remapping on the sphere (band-limited exact quadrature)
    sht = SHT(3 * lmax_in + 3)
    acc22 = np.zeros(lmax + 1); acc13 = np.zeros(lmax + 1); acc0 = np.zeros(lmax + 1)
    for r in range(nreal):
        t, p = correlated_realisation(sht, CT, Cphi, CphiT, rng, lmax_in)
        gT, gP = sht.gradient(t).real, sht.gradient(p).real
        d1 = np.sum(gP * gT, axis=-1)                      # grad phi . grad Theta
        d1_alm = sht.map2alm(d1, 2 * lmax_in)
        g_d1 = sht.gradient(d1_alm).real
        gp2_alm = sht.map2alm(np.sum(gP * gP, axis=-1), 2 * lmax_in)
        g_gp2 = sht.gradient(gp2_alm).real
        # delta2 = 1/2 [ grad phi . grad(grad phi . grad Theta) - 1/2 grad|grad phi|^2 . grad Theta ]
        d2 = 0.5 * (np.sum(gP * g_d1, axis=-1) - 0.5 * np.sum(g_gp2 * gT, axis=-1))
        d2_alm = sht.map2alm(d2, lmax)
        acc22 += cl_from_alm(d1_alm, lmax)
        acc13 += 2 * cross_cl_from_alm(t, d2_alm, lmax).real
        acc0 += cl_from_alm(t, lmax)
    acc22 /= nreal; acc13 /= nreal; acc0 /= nreal
    out.append(("MC: <|delta1_lm|^2>  vs  1/4 C^(22) (Hu 22-term + corr.)  (l=1..%d, nreal=%d)" % (lmax, nreal),
                np.max(np.abs(acc22[1:] - c22[1:]) / c22[1:])))
    out.append(("MC: 2Re<Theta delta2*>  vs  -l(l+1) R C_l  (Hu 13-term)   (l=1..%d, nreal=%d)" % (lmax, nreal),
                np.max(np.abs(acc13[1:] - c13[1:]) / np.abs(c13[1:]))))
    out.append(("MC: sanity, unlensed C_l recovered", np.max(np.abs(acc0[1:] - CT[1:lmax + 1]) / CT[1:lmax + 1])))
    if verbose:
        for n, e in out:
            print("  %-85s max rel dev = %.2e" % (n, e))
    return out


def test_typeD_identity(lmax=4, lmax_in=22, verbose=True):
    """Fourier transverse-gradient contraction (type-C-like polynomial Q, two radii) == angular-gradient Gaunt route."""
    rng = np.random.default_rng(2)
    sht = SHT(lmax + 8, ntheta=48, nphi=64)
    n = sht.n_vectors()
    k1 = rng.normal(size=3); k2 = rng.normal(size=3)
    chi1, chi2 = 1.7, 2.9
    f = (np.sum(-(k1 * k2)[None, None, :], axis=-1) + (n @ k1) * (n @ k2)) * np.exp(-1j * (n @ k1) * chi1 - 1j * (n @ k2) * chi2)
    alm = sht.map2alm(f, lmax)
    def shell(k, chi, l, m):
        kk = np.linalg.norm(k); kh = k / kk
        return 4 * np.pi * (-1j) ** l * np.conj(Y(l, m, np.arccos(kh[2]), np.arctan2(kh[1], kh[0]))) * jl(l, kk * chi)
    err = 0.0
    for l in range(lmax + 1):
        for m in range(-l, l + 1):
            s = 0.0
            for l1 in range(lmax_in + 1):
                for m1 in range(-l1, l1 + 1):
                    m2 = m - m1
                    for l2 in range(abs(l - l1), min(l + l1, lmax_in) + 1):
                        if abs(m2) > l2 or (l + l1 + l2) % 2:
                            continue
                        g = gaunt(l, m, l1, m1, l2, m2)
                        if g == 0:
                            continue
                        s += w_grad(l, l1, l2) * g * shell(k1, chi1, l1, m1) * shell(k2, chi2, l2, m2)
            s /= chi1 * chi2
            err = max(err, abs(s - alm[(l, m)]) / max(1e-300, abs(alm[(l, m)])))
    if verbose:
        print("  %-85s max rel err = %.2e" % ("type-D: Fourier (delta-nn) dX dY  ==  (1/chi1chi2) w_grad Gaunt shell-map route", err))
    return err


if __name__ == "__main__":
    print("lensing control case")
    test_lensing_control()
    test_typeD_identity()
