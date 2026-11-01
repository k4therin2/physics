"""
typeB.py -- products of two maps at the same n (type B) and the "monopole times map" terms.

Theta_B(n) = A(n) A'(n)   =>   a^B_lm = sum_{l1 m1 l2 m2} a_{l1 m1} a'_{l2 m2} G^{l l1 l2}_{m m1 m2}
(G = Gaunt, see harm_core).  More generally a "weighted product" a^B_lm = sum w(l,l1,l2) G a a' covers
the angular-gradient contraction  grad A . grad A'  with  w = (1/2)[l1(l1+1) + l2(l2+1) - l(l+1)]  (Hu 2000 eq. 47/50).

Wick (Gaussian a_lm, cross spectra C^{XY}_l := (1/(2l+1)) sum_m <a^X a^{Y*}>, <a^X_{lm} a^Y_{l'm'}> = (-1)^m C^{XY}_l delta_{ll'} delta_{m,-m'}):
 <a^B_lm a^{B'*}_lm> with B = A A', B' = C D:
   pairing 1: <a_{l1m1} a^{C*}_{l3m3}> <a'_{l2m2} a^{D*}_{l4m4}>  -> C^{AC}_{l1} C^{A'D}_{l2}
   pairing 2: <a_{l1m1} a^{D*}_{l4m4}> <a'_{l2m2} a^{C*}_{l3m3}>  -> C^{AD}_{l1} C^{A'C}_{l2}
   pairing 3: <a_{l1m1} a'_{l2m2}> <a^{C*} a^{D*}>                 -> only l = 0 (sum_{m1} (-1)^{m1} G^{l l1 l1}_{m m1 -m1} ∝ delta_{l0})
 and  sum_{m m1 m2} |G^{l l1 l2}_{m m1 m2}|^2 = (2l+1)(2l1+1)(2l2+1)/(4pi) (l l1 l2;000)^2 ,  so for l >= 1

   C^{BB'}_l = sum_{l1 l2} (2l1+1)(2l2+1)/(4pi) (l l1 l2;0 0 0)^2 w(l,l1,l2) w'(l,l1,l2)
               x [ C^{AC}_{l1} C^{A'D}_{l2} + C^{AD}_{l1} C^{A'C}_{l2} ]                                   (B1)

 (parity: only l+l1+l2 even contributes; w symmetric in l1<->l2 is assumed for pairing 2 [true for 1 and for the
  gradient weight]; for a non-symmetric weight replace w w' -> w(l,l1,l2) w'(l,l2,l1) in the second term.)
 Special case B = B' = A A':  C^{BB}_l = sum (2l1+1)(2l2+1)/(4pi) (000)^2 [ C^{AA}_{l1} C^{A'A'}_{l2} + C^{AA'}_{l1} C^{A'A}_{l2} ].

Monopole-times-map terms (SPEC §7:  Th2_lm = m2_lm - 2 <m1> m1_lm, etc.).  <X> := sky average = x_00/sqrt(4pi).
 Treat <X>.1 as a map with only l=0:  C^{<X>.1, Y}_{l} = 4 pi <<X><Y>> delta_{l0}  (<<X><Y>> = C^{XY}_0/(4pi)).
 Inserting into (B1) (only l1 = 0 or l2 = 0 survive, (l 0 l2;000)^2 = delta_{l l2}/(2l+1)):
   (i)   (1/(2l+1)) sum_m <|<X> y_lm|^2>              = <<X>^2> C^{YY}_l                      (l >= 1)
   (ii)  (1/(2l+1)) sum_m <(AA')_lm (<X> y_lm)^*>     = <<A><X>> C^{A'Y}_l + <<A'><X>> C^{AY}_l
   (iii) (1/(2l+1)) sum_m <(<X> y_lm)(<Z> w_lm)^*>    = <<X><Z>> C^{YW}_l        (l >= 1; the other pairings vanish by isotropy)
 In k-space (type-A maps): <<X><Y>> = (1/2pi^2) int k^2 dk P(k) Delta^X_0(k) Delta^Y_0(k)^*.
 All of these are exact consequences of (B1); they are also checked by Monte Carlo below.
"""
import numpy as np
from numpy.polynomial import legendre as L
from harm_core import w3j000, SHT, gaunt, cl_from_alm, cross_cl_from_alm


def w_plain(l, l1, l2):
    return 1.0


def w_grad(l, l1, l2):
    """angular-gradient contraction weight: int Y*_lm grad Y_{l1m1}.grad Y_{l2m2} = w * Gaunt"""
    return 0.5 * (l1 * (l1 + 1) + l2 * (l2 + 1) - l * (l + 1))


def Cl_product(CAC, CApD, CAD, CApC, lmax, lmax_in, w=w_plain, wp=None):
    """Formula (B1).  Inputs are arrays indexed by l (length >= lmax_in+1).  Returns C^{BB'}_l for l=0..lmax
    (the l=0 entry omits the disconnected pairing 3, i.e. it is the connected part only)."""
    wp = w if wp is None else wp
    out = np.zeros(lmax + 1)
    for l in range(lmax + 1):
        s = 0.0
        for l1 in range(lmax_in + 1):
            for l2 in range(abs(l - l1), min(l + l1, lmax_in) + 1):
                if (l + l1 + l2) % 2:
                    continue
                g = (2 * l1 + 1) * (2 * l2 + 1) / (4 * np.pi) * w3j000(l, l1, l2) ** 2
                s += g * (w(l, l1, l2) * wp(l, l1, l2) * CAC[l1] * CApD[l2]
                          + w(l, l1, l2) * wp(l, l2, l1) * CAD[l1] * CApC[l2])
        out[l] = s
    return out


# ------------------------------------------------------------------ real-space (correlation function) route
def xi_from_cl(Cl):
    """Legendre coefficients c_l = (2l+1)/(4pi) C_l so that xi(c) = sum c_l P_l(c)"""
    l = np.arange(len(Cl))
    return (2 * l + 1) / (4 * np.pi) * np.asarray(Cl)


def cl_from_xi_samples(xi_vals, c, wgl, lmax):
    """C_l = 2 pi int_{-1}^{1} xi(c) P_l(c) dc  by Gauss-Legendre"""
    return np.array([2 * np.pi * np.sum(wgl * xi_vals * L.legval(c, np.eye(lmax + 1)[l])) for l in range(lmax + 1)])


def Cl_product_realspace(CAC, CApD, CAD, CApC, lmax, w=w_plain, ngl=400):
    """Exact real-space evaluation of the connected product spectrum.
    w_plain : xi_B(c) = xi_AC xi_A'D + xi_AD xi_A'C
    w_grad  : xi_B(c) = xi_1'' xi_2'' (1-c^2)^2 - c(1-c^2)(xi_1'' xi_2' + xi_1' xi_2'') + (1+c^2) xi_1' xi_2' for each pairing,
              where xi_1 = xi_AC, xi_2 = xi_A'D (pairing 1) and xi_1 = xi_AD, xi_2 = xi_A'C (pairing 2)
              [from  grad^(1)_i grad^(2)_j xi(n1.n2) = xi'' t_i t'_j + xi' (P1 P2)_ij ].
    """
    c, wgl = np.polynomial.legendre.leggauss(ngl)
    def pair(C1, C2):
        a1, a2 = xi_from_cl(C1), xi_from_cl(C2)
        if w is w_plain:
            return L.legval(c, a1) * L.legval(c, a2)
        d1, d2 = L.legder(a1), L.legder(a2)
        dd1, dd2 = L.legder(a1, 2), L.legder(a2, 2)
        x1, x2 = L.legval(c, d1), L.legval(c, d2)
        y1, y2 = L.legval(c, dd1), L.legval(c, dd2)
        return y1 * y2 * (1 - c ** 2) ** 2 - c * (1 - c ** 2) * (y1 * x2 + x1 * y2) + (1 + c ** 2) * x1 * x2
    xi = pair(CAC, CApD) + pair(CAD, CApC)
    return cl_from_xi_samples(xi, c, wgl, lmax)


# ------------------------------------------------------------------ tests
def _toy_spectra(lmax_in, rng):
    l = np.arange(lmax_in + 1)
    CAA = 1.0 / (1 + l) ** 2 * (1 + 0.3 * rng.uniform(size=lmax_in + 1))
    CBB = 1.0 / (1 + l) ** 1.5 * (1 + 0.3 * rng.uniform(size=lmax_in + 1))
    rho = 0.6 * (1 - l / (lmax_in + 2))
    CAB = rho * np.sqrt(CAA * CBB)
    return CAA, CBB, CAB


def test_typeB_exact(lmax=6, lmax_in=10, verbose=True):
    rng = np.random.default_rng(7)
    CAA, CBB, CAB = _toy_spectra(lmax_in, rng)
    out = []
    for name, w in (("plain product", w_plain), ("gradient product", w_grad)):
        h = Cl_product(CAA, CBB, CAB, CAB, lmax, lmax_in, w=w)
        r = Cl_product_realspace(CAA, CBB, CAB, CAB, lmax, w=w)
        err = np.max(np.abs(h[1:] - r[1:]) / np.abs(r[1:]))
        out.append(("type-B C_l, 3j formula vs real-space Wick+Legendre (%s, correlated maps)" % name, err))
    # uncorrelated version (second pairing absent)
    z = 0 * CAB
    for name, w in (("plain product", w_plain), ("gradient product", w_grad)):
        h = Cl_product(CAA, CBB, z, z, lmax, lmax_in, w=w)
        r = Cl_product_realspace(CAA, CBB, z, z, lmax, w=w)
        err = np.max(np.abs(h[1:] - r[1:]) / np.abs(r[1:]))
        out.append(("type-B C_l, 3j vs real-space (%s, uncorrelated)" % name, err))
    if verbose:
        for n, e in out:
            print("  %-80s max rel err = %.2e" % (n, e))
    return out


def correlated_realisation(sht, CAA, CBB, CAB, rng, lmax_in):
    """two real Gaussian maps with spectra CAA, CBB and cross CAB"""
    a = sht.random_alm(CAA, rng, lmax_in)
    b_ind = sht.random_alm(np.ones(lmax_in + 1), rng, lmax_in)
    b = {}
    for (l, m), v in a.items():
        rho = CAB[l] / np.sqrt(CAA[l] * CBB[l]) if CAA[l] * CBB[l] > 0 else 0.0
        b[(l, m)] = np.sqrt(CBB[l]) * (rho * v / np.sqrt(CAA[l]) + np.sqrt(max(1 - rho ** 2, 0)) * b_ind[(l, m)])
    return a, b


def test_typeB_montecarlo(nreal=400, lmax=4, lmax_in=6, verbose=True):
    """Monte Carlo over Gaussian realisations on the sphere: product maps, gradient products, monopole x map."""
    rng = np.random.default_rng(11)
    CAA, CBB, CAB = _toy_spectra(lmax_in, rng)
    sht = SHT(2 * lmax_in + 2)        # product of two band-limited maps is band-limited at 2 lmax_in
    accP = np.zeros(lmax + 1); accG = np.zeros(lmax + 1); accM = np.zeros(lmax + 1); accX = np.zeros(lmax + 1, dtype=complex)
    for r in range(nreal):
        a, b = correlated_realisation(sht, CAA, CBB, CAB, rng, lmax_in)
        A, B = sht.alm2map(a).real, sht.alm2map(b).real
        prod = sht.map2alm(A * B, lmax)
        accP += cl_from_alm(prod, lmax)
        gA, gB = sht.gradient(a).real, sht.gradient(b).real
        gprod = sht.map2alm(np.sum(gA * gB, axis=-1), lmax)
        accG += cl_from_alm(gprod, lmax)
        meanA = a[(0, 0)].real / np.sqrt(4 * np.pi)
        mono = {k: meanA * v for k, v in b.items() if k[0] <= lmax}
        accM += cl_from_alm(mono, lmax)
        accX += cross_cl_from_alm(prod, mono, lmax)
    accP /= nreal; accG /= nreal; accM /= nreal; accX /= nreal
    predP = Cl_product(CAA, CBB, CAB, CAB, lmax, lmax_in, w=w_plain)
    predG = Cl_product(CAA, CBB, CAB, CAB, lmax, lmax_in, w=w_grad)
    varA = CAA[0] / (4 * np.pi)
    predM = varA * CBB[:lmax + 1]
    predX = (CAA[0] / (4 * np.pi)) * CBB[:lmax + 1] + (CAB[0] / (4 * np.pi)) * CAB[:lmax + 1]   # (ii) with A=A, A'=B, X=A, Y=B
    out = []
    for name, mc, pr in (("MC plain product", accP, predP), ("MC gradient product", accG, predG),
                         ("MC <A>.B  (monopole x map)", accM, predM), ("MC <(AB)_lm (<A> b_lm)^*>", accX.real, predX)):
        rel = np.abs(mc[1:] - pr[1:]) / np.abs(pr[1:])
        out.append((name + " (l=1..%d), nreal=%d" % (lmax, nreal), rel.max()))
    if verbose:
        for n, e in out:
            print("  %-70s max rel dev = %.3f  (expect ~ few/sqrt(nreal))" % (n, e))
    return out


if __name__ == "__main__":
    print("type-B exact tests")
    test_typeB_exact()
    print("type-B Monte Carlo tests")
    test_typeB_montecarlo()
