"""thirteen.py (A3) -- assembly of the contracted third-order map and of C^(13).
Th3 = m3 - 3<m1> m2 - 3<m2> m1 + 6<m1>^2 m1   (l >= 1),  m1 = L1, m2 = L2 + L1^2, m3 = L3 + 3 L1 L2 + L1^3,
with L1 = list1, L2 = 2 list2, L3 = 6 list3 (C2's lists are the eps^n coefficients of Lambda).
coef_abc(k,q,-q;n): mask-7 component in the polarisation algebra; kappa(k;n) = (1/2) int d^3q/(2pi)^3 P(q) coef_abc.
Monochromatic (P = 2 pi^2 sigma^2 delta(k-k0)/k0^2, sigma = 1):  kappa = (1/2) <coef_abc>_{qhat},  C13_l = 4 pi Delta_l (Delta^kappa_l)*,
kappa(nu) = sum_l (2l+1) i^l Delta^kappa_l P_l(nu),  nu = n.khat."""
import numpy as np
from numpy.polynomial import legendre as Lg
from scipy.special import spherical_jn, eval_legendre
import tf
from grid import TGrid
from engine import Ctx, RayGrid, LT, amul, aadd, ascale, apart
from third import X3c

G = TGrid(); RAY = RayGrid(eta_o=G.E)

def Delta_l(l, k, r=RAY):
    chs = r.eta_o - r.eta_s
    v = (tf.tdt(k, r.eta_s) + tf.tphi(k, r.eta_s)) * spherical_jn(l, k * chs) - k * tf.tv(k, r.eta_s) * spherical_jn(l, k * chs, derivative=True) \
        + r.integ(2 * tf.tphi(k, r.t, 1) * spherical_jn(l, k * r.chi))
    if l == 0: v -= tf.tdt(k, r.eta_o) + tf.tphi(k, r.eta_o)
    if l == 1: v += k * tf.tv(k, r.eta_o) / 3
    return v

def sky_grid(Nt=24, Nphi=48):
    x, w = Lg.leggauss(Nt)
    ph = 2 * np.pi * (np.arange(Nphi) + 0.5) / Nphi
    th = np.arccos(x)
    nv = np.stack([np.outer(np.sin(th), np.cos(ph)).ravel(), np.outer(np.sin(th), np.sin(ph)).ravel(), np.repeat(x, Nphi)], 1)
    wv = np.repeat(w, Nphi) / (2 * Nphi)            # sum wv = 1 (sky average)
    return nv, wv

def M2_sky(k1v, k2v, form='L_split', sky=(24, 48), tadshift=None):
    """sky average of the m2 AB-coefficient for momentum pairs k1v[b], k2v[b] (arrays (B,3)) -> (B,)."""
    nv, wv = sky_grid(*sky)
    ctx = Ctx([k1v, k2v], nv, g=G, ray=RAY, tadshift=tadshift)
    L1, _ = ctx.termlist(LT[form][1], maxorder=1)
    L2, _ = ctx.termlist(LT[form][2], maxorder=2)
    m2 = aadd(ascale(L2, 2.0), amul(L1, L1))[3]
    return m2 @ wv

def coef_abc(k0, qm, mu, phis, nv, form='L_split', third=None, M2=None, D0=None, tadshift=None, split=True):
    """mask-7 coefficient of Th3 for k = k0 zhat, q = qm (sqrt(1-mu^2) cos phi, sqrt(1-mu^2) sin phi, mu), -q.
    third: dict of abc-coefficient profiles (X3c(mu)+X3c(-mu)) or None; M2 = (M2(k,q), M2(k,-q), M2(q,-q)); D0 = (Delta0(k0), Delta0(qm)).
    Returns dict channel -> (B, Nn) array."""
    phis = np.atleast_1d(phis)
    s = np.sqrt(1 - mu**2)
    qv = qm * np.stack([s * np.cos(phis), s * np.sin(phis), np.full_like(phis, mu)], 1)
    kv = np.array([[0.0, 0.0, k0]])
    ctx = Ctx([kv, qv, -qv], nv, g=G, ray=RAY, third=third, tadshift=tadshift)
    L1, _ = ctx.termlist(LT[form][1], maxorder=1)
    L2h, _ = ctx.termlist(LT[form][2], maxorder=2)
    L2 = ascale(L2h, 2.0)
    out = {}
    L3h, parts = ctx.termlist(LT[form][3], maxorder=3)
    for nm, el in parts.items():
        if 7 in el: out['L3: ' + nm] = 6 * el[7]
    out['3 L1 L2'] = 3 * amul(L1, L2).get(7, 0)
    out['L1^3'] = amul(amul(L1, L1), L1).get(7, 0)
    if M2 is not None:
        m1 = L1; m2 = aadd(L2, amul(L1, L1))
        mon1 = {1: D0[0], 2: D0[1], 4: D0[1]}
        mon2 = {3: M2[0], 5: M2[1], 6: M2[2]}
        out['-3<m1>m2'] = -3 * amul(mon1, m2).get(7, 0)
        out['-3<m2>m1'] = -3 * amul(mon2, m1).get(7, 0)
        out['6<m1>^2m1'] = 6 * amul(amul(mon1, mon1), m1).get(7, 0)
    out = {k_: np.broadcast_to(v, (len(phis), len(nv))) for k_, v in out.items()}
    return out, ctx

def nu_grid(Nn=48):
    nu, wnu = Lg.leggauss(Nn)
    nv = np.stack([np.sqrt(1 - nu**2), np.zeros_like(nu), nu], 1)
    return nu, wnu, nv

def project(kappa_nu, nu, wnu, lmax):
    """kappa(nu) = sum (2l+1) i^l D_l P_l(nu)  ->  D_l."""
    return np.array([np.sum(wnu * eval_legendre(l, nu) * kappa_nu) / (2 * 1j**l) for l in range(lmax + 1)])
