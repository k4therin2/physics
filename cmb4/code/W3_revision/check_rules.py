"""W3: numerical check of the harmonic rules rewritten in the observation-direction convention
(ray x = x_o + n chi, phase e^{+i k.n chi}, Rayleigh e^{+ik.n chi} = 4pi sum i^l j_l Y_lm(n) Y*_lm(khat)).
k-frame khat = z.  Rules checked: (A-s), (A-lambda) for lambda = +-1, +-2, (A-L).  Quadrature on the sphere."""
import numpy as np
from scipy.special import sph_harm_y, spherical_jn, eval_legendre
from sympy.physics.wigner import wigner_3j

def jl_deriv(l, x, s):
    # s-th derivative of j_l by finite recursion (s <= 2)
    if s == 0: return spherical_jn(l, x)
    if s == 1: return spherical_jn(l, x, derivative=True)
    if s == 2:  # j'' = -2/x j' - (1 - l(l+1)/x^2) j
        return -2/x*spherical_jn(l, x, True) - (1 - l*(l+1)/x**2)*spherical_jn(l, x)
nth, nph = 200, 400
th, wth = np.polynomial.legendre.leggauss(nth); th = np.arccos(th)
ph = np.linspace(0, 2*np.pi, nph, endpoint=False); wph = 2*np.pi/nph
TH, PH = np.meshgrid(th, ph, indexing='ij'); W = wth[:, None]*wph
n = np.array([np.sin(TH)*np.cos(PH), np.sin(TH)*np.sin(PH), np.cos(TH)])
Y = lambda l, m: sph_harm_y(l, m, TH, PH)
alm = lambda f, l, m: np.sum(W*np.conj(Y(l, m))*f)
x = 3.7; mu = n[2]; phase = np.exp(1j*x*mu)
err = 0
# (A-s): (i n.k)^s e^{ikn chi} -> a_lm = 4pi i^l j_l^{(s)}(x) Y*_lm(z)
for s in (0, 1, 2):
    f = (1j*mu)**s*phase
    for l in range(0, 6):
        pred = 4*np.pi*1j**l*jl_deriv(l, x, s)*np.conj(sph_harm_y(l, 0, 0.0, 0.0))
        err = max(err, abs(alm(f, l, 0) - pred))
print('(A-s) max err', err)
# (A-lambda): Q_lam e^{ikn chi} = sum_l i^l sqrt(4pi(2l+1)) Dhat_l Y_{l lam}(n),
# Dhat_{+-1} = +-i sqrt(l(l+1)/2) j_l/x,  Dhat_{+-2} = -1/2 sqrt((l+2)!/(l-2)!) j_l/x^2
e1 = {+1: np.array([1, 1j, 0])/np.sqrt(2), -1: np.array([1, -1j, 0])/np.sqrt(2)}
err = 0
from math import factorial
for lam in (1, -1, 2, -2):
    if abs(lam) == 1: Q = np.einsum('i,i...->...', e1[lam], n)
    else: Q = np.einsum('i,i...->...', e1[np.sign(lam)], n)**2
    for l in range(abs(lam), 7):
        if abs(lam) == 1: D = np.sign(lam)*1j*np.sqrt(l*(l+1)/2)*spherical_jn(l, x)/x
        else: D = -0.5*np.sqrt(factorial(l+2)/factorial(l-2))*spherical_jn(l, x)/x**2
        pred = 1j**l*np.sqrt(4*np.pi*(2*l+1))*D
        err = max(err, abs(alm(Q*phase, l, lam) - pred))
        for m in range(-l, l+1):           # other m must vanish
            if m != lam: err = max(err, abs(alm(Q*phase, l, m)))
print('(A-lambda) max err', err)
# (A-L): P_L(n.k) e^{ikn chi} = sum_l (2l+1) i^l Delta_l P_l(n.k),  Delta_l = sum_J (2J+1) (l L J;000)^2 i^{J-l} j_J(x)
err = 0
for L in (0, 1, 2, 3):
    f = eval_legendre(L, mu)*phase
    for l in range(0, 7):
        D = sum((2*J+1)*float(wigner_3j(l, L, J, 0, 0, 0))**2*1j**(J-l)*spherical_jn(J, x) for J in range(abs(l-L), l+L+1))
        pred = (2*l+1)*1j**l*D*np.sqrt(4*np.pi/(2*l+1))   # coefficient of Y_l0 = (2l+1) i^l D P_l -> sqrt(4pi/(2l+1))
        err = max(err, abs(alm(f, l, 0) - pred))
print('(A-L) max err', err)
# addition theorem used for xi(0): sum_l (2l+1) j_l(x)^2 = 1 and Parseval form sum (2l+1) D_l^X D_l^Y* = 1/2 int dmu x(mu) y(mu)*
print('sum (2l+1) j_l^2 - 1 =', sum((2*l+1)*spherical_jn(l, x)**2 for l in range(80)) - 1)
# M8 tail bound: |j_l(x)| <= x^l/(2l+1)!! <= (e x/(2l+1))^l  and |j_l'(x)| <= (l/x + 1) x^l/(2l+1)!!
import mpmath as mp
worst = 0
for X in (0.5, 3.0, 10.0, 40.0):
    for l in range(0, 160):
        jl = abs(mp.besselj(l + 0.5, X)*mp.sqrt(mp.pi/(2*X)))
        b1 = mp.mpf(X)**l/mp.fac2(2*l + 1); b2 = (mp.e*X/(2*l + 1))**l
        worst = max(worst, float(jl/b1), float(b1/b2))
print('max of |j_l|/bound1 and bound1/bound2 over grid (must be <= 1):', worst)
