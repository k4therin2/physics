"""
typeA.py -- linear line-of-sight maps (type A) : a_lm rules and C_l^(11).

Convention (task statement): background ray x(eta) = -n chi, chi = eta_o - eta, so every linear map reads
   Theta_A(n) = int d^3k/(2pi)^3 zeta_k int d eta sum_s g_s(k,eta) (i n.khat)^s e^{-i k.n chi}      (scalar part)
             + [ helicity pieces  g^{(lam)}(k,eta) Q_lam(n,khat) e^{-i k.n chi},  lam = +-1, +-2 ].
We write the angular monomials as (i n.khat)^s because i khat = grad/k, so a REAL map has REAL g_s.

Rayleigh (x = -n chi):   e^{-i k.n chi} = 4 pi sum_{lm} (-i)^l j_l(k chi) Y_lm(n) Y*_lm(khat).
   [ With the opposite sign x(eta) = +n chi one gets (+i)^l ; C_l is unchanged, cross-spectra between maps
     are unchanged as long as ONE convention is used for all maps (a_lm -> (-1)^l a_lm for every map). ]

RULE A-s (scalar monomials):  (i n.khat)^s e^{-i k.n chi}  -->  (-d/dx)^s j_l(x) ,  x = k chi,
   i.e.   a_lm = 4 pi (-i)^l int d^3k/(2pi)^3 zeta_k Y*_lm(khat) Delta_l(k),
          Delta_l(k) = int d eta sum_s g_s(k,eta) (-1)^s j_l^{(s)}(k chi(eta))
   (endpoint / delta(eta-eta_1) pieces: replace int d eta -> evaluation at eta_1).
   Equivalently (n.khat)^s --> (i d/dx)^s j_l, which for s=1 gives i j_l', s=2 gives -j_l''.

RULE A-lam (helicity, frame khat = z, e^{(+-1)} = (xhat +- i yhat)/sqrt2, e^{(+-2)}_ij = e^{(+-1)}_i e^{(+-1)}_j):
   Q_{+-1} = n.e^{(+-1)} = sin(th) e^{+- i phi}/sqrt2 ,   Q_{+-2} = e^{(+-2)}_ij n^i n^j = sin^2(th) e^{+-2 i phi}/2
   Q_{+-1} e^{-i x cos th} = sum_l  (-i)^l sqrt(4pi(2l+1)) [ -+ i sqrt(l(l+1)/2) j_l(x)/x ]        Y_{l,+-1}(n)
   Q_{+-2} e^{-i x cos th} = sum_l  (-i)^l sqrt(4pi(2l+1)) [ -(1/2) sqrt((l+2)!/(l-2)!) j_l(x)/x^2 ] Y_{l,+-2}(n)
   (the signs are fixed by the numerical tests below; see test_typeA()).
   General frame: Y_{l lam}(n)|_{k-frame} = sum_m D^l_{m lam}(phi_k, theta_k, 0) Y_lm(n), and
   D^l_{m lam}(phi_k,theta_k,0) = (-1)^lam sqrt(4pi/(2l+1)) [ _{-lam}Y_{lm}(khat) ]^*   (checked numerically)
   so that   a^{(lam)}_lm = 4 pi (-i)^l int d^3k/(2pi)^3 zeta^{(lam)}_k (-1)^lam [_{-lam}Y_lm(khat)]^* Delta^{(lam)}_l(k)
   with Delta^{(+-1)}_l = -+ i sqrt(l(l+1)/2) int d eta g^{(+-1)} j_l(x)/x ,   (sign: Y_{l1} = -sqrt(.) sin(th) P_l'(cos th) e^{i phi})
        Delta^{(+-2)}_l = -(1/2) sqrt((l+2)!/(l-2)!) int d eta g^{(+-2)} j_l(x)/x^2 .
   For C_l only |Delta^{(lam)}_l|^2 matters (sum_m |a_lm|^2 is frame independent).

C_l^(11) (any helicity, maps A and A'):   C^{AA'}_l = (2/pi) int k^2 dk P(k) Delta^A_l(k) Delta^{A'}_l(k)^* .
   Derivation: <a a'^*> = (4pi)^2 int d^3k/(2pi)^3 P |Y_lm(khat)|^2 Delta Delta'^*; sum_m |Y_lm|^2 = (2l+1)/(4pi);
   (4pi)^2 (2l+1)/(4pi) * 4pi/(2pi)^3 = (2l+1) * 2/pi.    [verified numerically in test_C11_normalisation()]
Sky monopole of a type-A map:  <A> := (1/4pi) int dOmega A = a_00/sqrt(4pi) = int d^3k/(2pi)^3 zeta_k Delta^A_0(k),
   <<A><A'>> = (1/(2 pi^2)) int k^2 dk P(k) Delta^A_0 Delta^{A'*}_0  = C^{AA'}_0/(4 pi).
"""
import numpy as np
from math import factorial, comb
from scipy.special import eval_legendre, spherical_jn
from harm_core import SHT, Y, jl_d, gaunt, w3j


# ------------------------------------------------------------------ analytic radial rules
def Delta_scalar_monomial(s, l, x):
    """radial function for (i n.khat)^s : (-1)^s j_l^{(s)}(x)"""
    return (-1) ** s * jl_d(l, x, s)


def Delta_helicity(lam, l, x):
    """radial function multiplying 4pi(-i)^l (-1)^lam [_{-lam}Y_lm(khat)]^* for Q_lam (see module docstring)."""
    x = np.asarray(x, dtype=float)
    if abs(lam) == 1:
        if l < 1:
            return 0.0 * x
        return -np.sign(lam) * 1j * np.sqrt(l * (l + 1) / 2.0) * spherical_jn(l, x) / x
    if abs(lam) == 2:
        if l < 2:
            return 0.0 * x
        return -0.5 * np.sqrt(factorial(l + 2) / factorial(l - 2)) * spherical_jn(l, x) / x ** 2
    raise ValueError


# ------------------------------------------------------------------ Wigner d / spin-weighted harmonics (independent defs)
def wigner_d(l, m, mp, beta):
    """d^l_{m m'}(beta), Wigner small-d, explicit sum (Wikipedia/Edmonds convention, d^l_{m m'} = <lm|e^{-i beta Jy}|lm'>)."""
    if abs(m) > l or abs(mp) > l:
        return 0.0
    pref = np.sqrt(factorial(l + m) * factorial(l - m) * factorial(l + mp) * factorial(l - mp))
    s = 0.0
    for k in range(0, 2 * l + 1):
        a, b, c, d = l + mp - k, k, m - mp + k, l - m - k
        if min(a, b, c, d) < 0:
            continue
        s += ((-1) ** (m - mp + k) * np.cos(beta / 2) ** (2 * l + mp - m - 2 * k) * np.sin(beta / 2) ** (m - mp + 2 * k)
              / (factorial(a) * factorial(b) * factorial(c) * factorial(d)))
    return pref * s


def wigner_D(l, m, mp, alpha, beta, gamma):
    """D^l_{m m'}(alpha,beta,gamma) = e^{-i m alpha} d^l_{m m'}(beta) e^{-i m' gamma}; rotation R = R_z(alpha)R_y(beta)R_z(gamma)
    acting on functions as (R f)(n) = f(R^{-1} n):  Y_{lm'}(R^{-1} n) = sum_m D^l_{m m'}(R) Y_lm(n)."""
    return np.exp(-1j * m * alpha) * wigner_d(l, m, mp, beta) * np.exp(-1j * mp * gamma)


def sYlm(s, l, m, theta, phi):
    """spin-weighted spherical harmonic, Goldberg et al. (1967) explicit formula (reduces to Y_lm for s=0)."""
    if abs(m) > l or abs(s) > l:
        return 0.0
    pref = ((-1) ** m * np.sqrt(factorial(l + m) * factorial(l - m) * (2 * l + 1)
                                / (4 * np.pi * factorial(l + s) * factorial(l - s))))
    tot = 0.0
    for r in range(0, l - s + 1):
        if r + s - m < 0 or r + s - m > l + s:
            continue
        tot += (comb(l - s, r) * comb(l + s, r + s - m) * (-1) ** (l - r - s)
                * np.tan(theta / 2) ** (-(2 * r + s - m)))
    return pref * np.sin(theta / 2) ** (2 * l) * np.exp(1j * m * phi) * tot


# ------------------------------------------------------------------ tests
def test_typeA(lmax=6, verbose=True):
    rng = np.random.default_rng(1)
    sht = SHT(lmax + 6, ntheta=40, nphi=48)   # generous quadrature (integrand not band-limited: e^{-i k.n chi})
    n = sht.n_vectors()
    report = []
    # (1) scalar monomial rule in a random frame
    for trial in range(2):
        khat = rng.normal(size=3); khat /= np.linalg.norm(khat)
        x = rng.uniform(2.0, 6.0)
        mu = n @ khat
        thk, phk = np.arccos(khat[2]), np.arctan2(khat[1], khat[0])
        for s in range(0, 4):
            f = (1j * mu) ** s * np.exp(-1j * x * mu)
            alm = sht.map2alm(f, lmax)
            err = 0.0
            for l in range(lmax + 1):
                for m in range(-l, l + 1):
                    pred = 4 * np.pi * (-1j) ** l * np.conj(Y(l, m, thk, phk)) * Delta_scalar_monomial(s, l, x)
                    err = max(err, abs(alm[(l, m)] - pred))
            report.append(("scalar rule s=%d x=%.2f" % (s, x), err))
    # (2) helicity rules in the k-frame
    th, ph = sht.TH, sht.PH
    for lam in (1, -1, 2, -2):
        x = 4.3
        if abs(lam) == 1:
            Q = np.sin(th) * np.exp(1j * np.sign(lam) * ph) / np.sqrt(2)
        else:
            Q = np.sin(th) ** 2 * np.exp(2j * np.sign(lam) * ph) / 2
        alm = sht.map2alm(Q * np.exp(-1j * x * np.cos(th)), lmax)
        err = 0.0
        for l in range(lmax + 1):
            for m in range(-l, l + 1):
                pred = (-1j) ** l * np.sqrt(4 * np.pi * (2 * l + 1)) * Delta_helicity(lam, l, x) if m == lam else 0.0
                err = max(err, abs(alm[(l, m)] - pred))
        report.append(("helicity rule lam=%+d (k-frame)" % lam, err))
    # (3) helicity in a random frame: polarisation vectors from the rotation R(phi_k,theta_k,0) z -> khat
    for lam in (1, -1, 2):
        thk, phk = rng.uniform(0.3, 2.8), rng.uniform(0, 2 * np.pi)
        khat = np.array([np.sin(thk) * np.cos(phk), np.sin(thk) * np.sin(phk), np.cos(thk)])
        e1 = np.array([np.cos(thk) * np.cos(phk), np.cos(thk) * np.sin(phk), -np.sin(thk)])   # R xhat
        e2 = np.array([-np.sin(phk), np.cos(phk), 0.0])                                        # R yhat
        epm = (e1 + 1j * np.sign(lam) * e2) / np.sqrt(2)
        x = 3.7
        if abs(lam) == 1:
            Q = n @ epm
        else:
            Q = (n @ epm) ** 2
        alm = sht.map2alm(Q * np.exp(-1j * x * (n @ khat)), lmax)
        errD, errS = 0.0, 0.0
        for l in range(lmax + 1):
            for m in range(-l, l + 1):
                predD = ((-1j) ** l * np.sqrt(4 * np.pi * (2 * l + 1)) * Delta_helicity(lam, l, x)
                         * wigner_D(l, m, lam, phk, thk, 0.0))
                predS = (4 * np.pi * (-1j) ** l * (-1) ** lam * np.conj(sYlm(-lam, l, m, thk, phk))
                         * Delta_helicity(lam, l, x))
                errD = max(errD, abs(alm[(l, m)] - predD))
                errS = max(errS, abs(alm[(l, m)] - predS))
        report.append(("helicity lam=%+d random frame, D-matrix form" % lam, errD))
        report.append(("helicity lam=%+d random frame, spin-weighted form" % lam, errS))
    # (4) sanity: sYlm with s=0 equals Y_lm
    err = max(abs(sYlm(0, l, m, 0.7, 1.9) - Y(l, m, 0.7, 1.9)) for l in range(4) for m in range(-l, l + 1))
    report.append(("sYlm(s=0) == Y_lm", err))
    if verbose:
        for name, e in report:
            print("  %-55s max abs err = %.2e" % (name, e))
    return report


def test_C11_normalisation(verbose=True):
    """Monochromatic check of C_l = (2/pi) int k^2 dk P |Delta_l|^2 via real-space correlation function:
    Theta = int d^3k/(2pi)^3 zeta_k f e^{-i k.n chi}, P = 2 pi^2 sigma^2 delta(k-k0)/k0^2 ->
    xi(n1.n2) = sigma^2 |f|^2 j_0(k0 chi |n1-n2|),  C_l = 2 pi int_{-1}^1 xi(c) P_l(c) dc  vs  4 pi sigma^2 |f|^2 j_l(k0 chi)^2."""
    x = 5.1
    c, w = np.polynomial.legendre.leggauss(200)
    xi = spherical_jn(0, x * np.sqrt(2 - 2 * c))
    errs = []
    for l in range(0, 8):
        Cl_num = 2 * np.pi * np.sum(w * xi * eval_legendre(l, c))
        Cl_an = 4 * np.pi * spherical_jn(l, x) ** 2
        errs.append(abs(Cl_num - Cl_an))
    if verbose:
        print("  C11 normalisation (monochromatic, real-space route) max abs err = %.2e" % max(errs))
    return max(errs)


if __name__ == "__main__":
    print("type-A tests")
    test_typeA()
    test_C11_normalisation()
