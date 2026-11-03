"""X1: compare Adam's first-order transfer function (Prop. Tl / eq. (54) of his report)
     T_l = -3 int_{x_*}^{x_o} s'(x) w_l(sqrt3 (x_o - x)) dx,  s = 2 sin x/x - cos x,  w_l = j_l/3 + j_l''
with our closed form Delta_l = A(x_*) j_l(k chi) - k T_V(x_*) j_l'(k chi) + I_l(k)   (C5 linear_cl.py, ISW via Si/Ci).
Conventions: Adam's t is conformal time (a = a1 t), x = k t/sqrt3 = our k eta/sqrt3; both use
P(k) = 2 pi^2 sigma^2 delta(k-k0)/k0^2 with <zeta^2> = sigma^2, and C11 = 4 pi sigma^2 T_l^2 (Adam (66)) = 4 pi sigma^2 Delta_l^2 (ours).
Adam's T_l is implemented here from his definitions (nothing imported from his repo)."""
import sys, os
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'C5_linear_Cl'))
import linear_cl as L

mp.mp.dps = 40
s3 = mp.sqrt(3)

def jl(l, z):
    return L._mp_jl(l, z)

def jl_pp(l, z):  # second derivative via the Bessel ODE: j'' = -2/z j' - (1 - l(l+1)/z^2) j
    z = mp.mpf(z)
    jp = mp.diff(lambda t: jl(l, t), z)
    return -2 / z * jp - (1 - l * (l + 1) / z**2) * jl(l, z)

def T_adam(l, k, ts, to):
    xs, xo = k * ts / s3, k * to / s3
    sp_ = lambda x: -2 * mp.sin(x) / x**2 + 2 * mp.cos(x) / x + mp.sin(x)   # s'(x)
    def w(y):
        if y == 0:
            return mp.mpf(2) / 15 if l == 2 else mp.mpf(0)
        return jl(l, y) / 3 + mp.diff(lambda t: jl(l, t), y, 2)
    n = int(xo - xs) * 2 + 4
    pts = [xs + (xo - xs) * i / n for i in range(n + 1)]
    return -3 * mp.quad(lambda x: sp_(x) * w(s3 * (xo - x)), pts)

def Delta_ours(l, k, es, eo):
    k, es, eo = mp.mpf(k), mp.mpf(es), mp.mpf(eo)
    xs = k * es / s3; c = k * (eo - es)
    j0 = lambda z: mp.sin(z) / z
    j1 = lambda z: mp.sin(z) / z**2 - mp.cos(z) / z
    A = -mp.cos(xs) + 2 * j0(xs) - 4 * j1(xs) / xs
    kTV = s3 * (mp.sin(xs) - 2 * j1(xs))
    jlp = mp.diff(lambda t: jl(l, t), c)
    em = A * jl(l, c) - kTV * jlp
    I = L.I_l_closed_mp(l, k, es, eo, dps=40)   # closed form (Si/Ci)
    Iq = L.I_l_quad_mp(l, k, es, eo, dps=40)    # quadrature, as a control of the closed form
    return em + I, em + Iq

if __name__ == '__main__':
    cfgs = [(1, 2, 8), (0.3, 2, 8), (3.0, 2, 8), (1, 1, 3), (2.5, 0.7, 5.2), (0.05, 2, 8), (10, 1, 4)]
    print('%-18s %3s %-28s %-28s %-10s %-10s' % ('(k0,eta*,eta_o)', 'l', 'Delta_ours (closed Si/Ci)', 'T_adam (eq. 54)', 'reldiff', 'cl-vs-quad'))
    worst = 0
    for (k, es, eo) in cfgs:
        for l in (2, 3, 4, 6):
            D, Dq = Delta_ours(l, k, es, eo)
            T = T_adam(l, mp.mpf(k), mp.mpf(es), mp.mpf(eo))
            rd = abs(D - T) / abs(T)
            worst = max(worst, rd)
            print('%-18s %3d %-28s %-28s %-10s %-10s' % ((k, es, eo), l, mp.nstr(D, 20), mp.nstr(T, 20), mp.nstr(rd, 3), mp.nstr(abs(D - Dq) / abs(T), 3)))
    print('worst relative difference Delta_ours vs T_adam:', mp.nstr(worst, 3))
