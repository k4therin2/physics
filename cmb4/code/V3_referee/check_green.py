"""check_green.py -- V3: Green's functions (Wronskian sign), ghat -> 1, finiteness of S_Phi at eta->0 (my sources),
and the value K_SPhi(eta=0) vs C1's SPhi0_2nd."""
import sys, sympy as sp
from pw_tools import *
import pw_tools
C1 = '/home/user/physics/cmb4/code/C1_einstein'; sys.path.insert(0, C1)
import kernels_api as KA; K = KA.load()
from fourier import eta, k, q, p, mu
etap = sp.Symbol('etap', positive=True); kk = sp.Symbol('k', positive=True)
# Phi operator L = d^2 + (4/eta) d + k^2/3 ; particular solution Phi = int_0^eta G(eta,eta') S(eta') deta' requires
# G(eta,eta)=0, d_eta G|_{eta'=eta} = 1, L_eta G = 0.
G_Phi = K['second']['G_Phi'].subs(k, kk); G_h = K['second']['G_h'].subs(k, kk)
S = sp.Function('S')
Phi = sp.integrate(G_Phi * S(etap), (etap, 0, eta))
LPhi = sp.diff(Phi, eta, 2) + 4 / eta * sp.diff(Phi, eta) + kk**2 / 3 * Phi
print("L[int G_Phi S] - S  =", sp.simplify(LPhi - S(eta)), "  (Leibniz rule: 0 iff G(eta,eta)=0 and d_eta G|=1)")
h = sp.integrate(G_h * S(etap), (etap, 0, eta))
print("L_h[int G_h S] - S =", sp.simplify(sp.diff(h, eta, 2) + 2 / eta * sp.diff(h, eta) + kk**2 * h - S(eta)))
print("ghat = 3 u1 at eta->0:", sp.limit(K['second']['ghat'].subs(k, kk), eta, 0), "; ghat solves homogeneous ODE:",
      sp.simplify((lambda g_: sp.diff(g_, eta, 2) + 4 / eta * sp.diff(g_, eta) + kk**2 / 3 * g_)(K['second']['ghat'].subs(k, kk))) == 0)
# S_Phi finiteness at eta -> 0 from MY sources with symbolic (k,q,mu), Taylor-expanded transfer functions
E = E_mine()
kq, qq, m = sp.symbols('k q mu', positive=True)
def Tser(kmag):
    xv = kmag * t / sp.sqrt(3)
    return sp.series(-2 * (sp.sin(xv) - xv * sp.cos(xv)) / xv**3, t, 0, 10).removeO()
pw_tools.T_Phi = Tser
zq, zp = sp.symbols('zeta_q zeta_p')
kv, qv, pv = vec3(kq, qq, m); pm = sp.sqrt(kq**2 + qq**2 - 2 * kq * qq * m)
f = {}
for nm, v_ in pw_tools.mode(zq, qv, qq).items(): f[nm] = f.get(nm, 0) + v_
for nm, v_ in pw_tools.mode(zp, pv, pm).items(): f[nm] = f.get(nm, 0) + v_
src = lambda a, b: sp.expand(coeff_at_origin(substitute(E[('E', a, b, 2)], f), (zq, zp)))
S00 = src(0, 0); Sxx, Syy, Szz = src(1, 1), src(2, 2), src(3, 3)
STL = Szz - (Sxx + Syy) / 2; Str = Sxx + Syy + Szz
Q = STL / kq**2
SPhi = sp.expand(-(Str - S00) / 6 + sp.diff(Q, t) / t - kq**2 / 3 * Q)
ser = sp.series(SPhi, t, 0, 3).removeO()
print("S_Phi small-eta series (my sources): leading power present:", [n for n in range(-3, 3) if sp.simplify(ser.coeff(t, n)) != 0])
SPhi0 = sp.simplify(ser.coeff(t, 0))
C1_SPhi0 = K['third_ic']['SPhi0_2nd'].subs(p, pm).subs({k: kq, q: qq})
print("K_SPhi(eta=0) mine - C1:", sp.simplify(SPhi0 - C1_SPhi0))
print("K_SPhi(eta=0) =", sp.factor(SPhi0))
