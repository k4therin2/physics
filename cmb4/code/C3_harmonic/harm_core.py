"""
harm_core.py -- shared angular-momentum tools for agent C3 (harmonic-space machinery).

Conventions (SPEC §5): f(x) = int d^3k/(2pi)^3 e^{ik.x} f_k ;  <zeta_k zeta_k'> = (2pi)^3 delta^3(k+k') P(k).
Spherical harmonics: Condon-Shortley phase (scipy.special.sph_harm_y == sympy Ynm, checked).
Gaunt integral  G^{l l1 l2}_{m m1 m2} := int dOmega Y*_{lm} Y_{l1 m1} Y_{l2 m2}
      = (-1)^m sqrt((2l+1)(2l1+1)(2l2+1)/(4pi)) (l l1 l2; 0 0 0) (l l1 l2; -m m1 m2).
Exact 3j symbols from sympy.physics.wigner (rational arithmetic), cached and converted to float.
"""
import functools
import numpy as np
from scipy.special import sph_harm_y, spherical_jn, eval_legendre
from sympy.physics.wigner import wigner_3j
from sympy import Rational, sqrt as ssqrt, pi as spi, Integer


@functools.lru_cache(maxsize=None)
def w3j(l1, l2, l3, m1, m2, m3):
    """Exact Wigner 3j symbol as a float (sympy exact evaluation, cached)."""
    if abs(m1) > l1 or abs(m2) > l2 or abs(m3) > l3 or m1 + m2 + m3 != 0:
        return 0.0
    if l3 < abs(l1 - l2) or l3 > l1 + l2:
        return 0.0
    return float(wigner_3j(l1, l2, l3, m1, m2, m3))


@functools.lru_cache(maxsize=None)
def w3j000(l1, l2, l3):
    return w3j(l1, l2, l3, 0, 0, 0)


def gaunt(l, m, l1, m1, l2, m2):
    """G^{l l1 l2}_{m m1 m2} = int dOmega Y*_{lm} Y_{l1m1} Y_{l2m2}  (float)."""
    if m != m1 + m2:
        return 0.0
    if (l + l1 + l2) % 2 or l < abs(l1 - l2) or l > l1 + l2:
        return 0.0
    return ((-1) ** m * np.sqrt((2 * l + 1) * (2 * l1 + 1) * (2 * l2 + 1) / (4 * np.pi))
            * w3j000(l, l1, l2) * w3j(l, l1, l2, -m, m1, m2))


def Y(l, m, theta, phi):
    """Y_lm(theta, phi), Condon-Shortley phase; theta = polar angle."""
    if abs(m) > l:
        return np.zeros(np.broadcast(theta, phi).shape, dtype=complex)
    return sph_harm_y(l, m, theta, phi)


def jl(l, x):
    return spherical_jn(l, x)


def jl_d(l, x, n):
    """n-th derivative of j_l w.r.t. its argument (n = 0,1,2,3)."""
    if n == 0:
        return spherical_jn(l, x)
    if n == 1:
        return spherical_jn(l, x, derivative=True)
    # use recursion j_l' = j_{l-1} - (l+1)/x j_l  applied symbolically via finite combos:
    # simplest robust route: numerical differentiation of the exact d/dx formula is avoided;
    # instead use the identity j_l'' = [ (l(l+1)/x^2 - 1) j_l - (2/x) j_l' ]  (spherical Bessel ODE)
    j0 = spherical_jn(l, x)
    j1 = spherical_jn(l, x, derivative=True)
    x = np.asarray(x, dtype=float)
    j2 = (l * (l + 1) / x ** 2 - 1.0) * j0 - 2.0 / x * j1
    if n == 2:
        return j2
    if n == 3:
        # differentiate the ODE once more
        return ((l * (l + 1) / x ** 2 - 1.0) * j1 - 2 * l * (l + 1) / x ** 3 * j0
                + 2.0 / x ** 2 * j1 - 2.0 / x * j2)
    raise ValueError("n<=3")


# ---------------------------------------------------------------- Legendre monomial expansion
@functools.lru_cache(maxsize=None)
def monomial_legendre_coeffs(a):
    """x^a = sum_{l'} alpha_{a l'} P_{l'}(x);  returns dict l' -> alpha (exact rationals as float).
    alpha_{a l'} = (2l'+1)/2 int_{-1}^1 x^a P_{l'}(x) dx  (nonzero only for l' <= a, l' = a mod 2)."""
    from sympy import legendre, symbols, integrate
    x = symbols('x')
    out = {}
    for lp in range(a % 2, a + 1, 2):
        val = Rational(2 * lp + 1, 2) * integrate(x ** a * legendre(lp, x), (x, -1, 1))
        out[lp] = float(val)
    return out


# ---------------------------------------------------------------- minimal exact SHT on a Gauss-Legendre grid
class SHT:
    """Spherical-harmonic transform on a Gauss-Legendre (theta) x uniform (phi) grid.
    Exact (to rounding) for functions band-limited at lmax if ntheta >= lmax+1 and nphi >= 2 lmax+1.
    Intended only for self-tests (small lmax)."""

    def __init__(self, lmax, ntheta=None, nphi=None):
        self.lmax = lmax
        self.nt = ntheta or (lmax + 2)
        self.np_ = nphi or (2 * lmax + 2)
        xg, wg = np.polynomial.legendre.leggauss(self.nt)
        self.theta = np.arccos(xg)[::-1]
        self.wtheta = wg[::-1]
        self.phi = 2 * np.pi * np.arange(self.np_) / self.np_
        self.wphi = 2 * np.pi / self.np_
        self.TH, self.PH = np.meshgrid(self.theta, self.phi, indexing='ij')
        # precompute Y_lm on grid: dict (l,m) -> array
        self.Ygrid = {}
        for l in range(lmax + 1):
            for m in range(-l, l + 1):
                self.Ygrid[(l, m)] = Y(l, m, self.TH, self.PH)
        self.W = self.wtheta[:, None] * self.wphi

    def n_vectors(self):
        """unit vectors n on the grid, shape (nt, nphi, 3)"""
        st, ct = np.sin(self.TH), np.cos(self.TH)
        return np.stack([st * np.cos(self.PH), st * np.sin(self.PH), ct], axis=-1)

    def integrate(self, f):
        return np.sum(self.W * f)

    def map2alm(self, f, lmax=None):
        lmax = self.lmax if lmax is None else lmax
        alm = {}
        for l in range(lmax + 1):
            for m in range(-l, l + 1):
                alm[(l, m)] = np.sum(self.W * np.conj(self.Ygrid[(l, m)]) * f)
        return alm

    def alm2map(self, alm):
        f = np.zeros_like(self.TH, dtype=complex)
        for (l, m), a in alm.items():
            if a != 0:
                f += a * self.Ygrid[(l, m)]
        return f

    def gradient(self, alm):
        """covariant gradient on the sphere of the map sum alm Y_lm, returned as 3-vectors in R^3
        (tangent to the sphere), shape (nt,nphi,3).  Uses grad Y_lm = dY/dtheta e_theta + (1/sin) dY/dphi e_phi
        with dY/dtheta from the ladder identity and dY/dphi = i m Y."""
        f_th = np.zeros_like(self.TH, dtype=complex)
        f_ph = np.zeros_like(self.TH, dtype=complex)
        st, ct = np.sin(self.TH), np.cos(self.TH)
        for (l, m), a in alm.items():
            if a == 0:
                continue
            Ylm = self.Ygrid[(l, m)]
            # dY_lm/dtheta = m cot(theta) Y_lm + sqrt((l-m)(l+m+1)) e^{-i phi} Y_{l,m+1}
            Ylm1 = self.Ygrid[(l, m + 1)] if m + 1 <= l else 0.0
            dth = m * ct / st * Ylm + np.sqrt((l - m) * (l + m + 1)) * np.exp(-1j * self.PH) * Ylm1
            f_th += a * dth
            f_ph += a * (1j * m) * Ylm / st
        e_th = np.stack([ct * np.cos(self.PH), ct * np.sin(self.PH), -st], axis=-1)
        e_ph = np.stack([-np.sin(self.PH), np.cos(self.PH), 0 * st], axis=-1)
        return f_th[..., None] * e_th + f_ph[..., None] * e_ph

    def random_alm(self, Cl, rng, lmax=None, real=True):
        """Gaussian realisation with spectrum Cl[l]; real map (a_{l,-m} = (-1)^m a_lm^*)."""
        lmax = self.lmax if lmax is None else lmax
        alm = {}
        for l in range(lmax + 1):
            s = np.sqrt(max(Cl[l], 0.0))
            alm[(l, 0)] = s * rng.normal()
            for m in range(1, l + 1):
                z = s * (rng.normal() + 1j * rng.normal()) / np.sqrt(2)
                alm[(l, m)] = z
                alm[(l, -m)] = (-1) ** m * np.conj(z)
        return alm


def cl_from_alm(alm, lmax):
    """(1/(2l+1)) sum_m |a_lm|^2"""
    return np.array([sum(abs(alm[(l, m)]) ** 2 for m in range(-l, l + 1)) / (2 * l + 1) for l in range(lmax + 1)])


def cross_cl_from_alm(alm, blm, lmax):
    return np.array([sum(alm[(l, m)] * np.conj(blm[(l, m)]) for m in range(-l, l + 1)) / (2 * l + 1)
                     for l in range(lmax + 1)])
