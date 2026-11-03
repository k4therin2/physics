"""X1 (c): late-time (x_* fixed) and subhorizon (k0 -> inf at fixed eta's) behaviour of our closed-form Delta_l (closed form only)."""
import mpmath as mp
import transfer_compare as T
L = T.L; s3 = mp.sqrt(3); mp.mp.dps = 40
def D(l, k, es, eo):
    k, es, eo = mp.mpf(k), mp.mpf(es), mp.mpf(eo); xs = k * es / s3; c = k * (eo - es)
    j0 = lambda z: mp.sin(z) / z; j1 = lambda z: mp.sin(z) / z**2 - mp.cos(z) / z
    A = -mp.cos(xs) + 2 * j0(xs) - 4 * j1(xs) / xs; kTV = s3 * (mp.sin(xs) - 2 * j1(xs))
    return A * T.jl(l, c) - kTV * mp.diff(lambda t: T.jl(l, t), c) + L.I_l_closed_mp(l, k, es, eo, dps=40)
print('(c1) late time x_* = 1 fixed: Delta*x_o')
for xo in (10, 100, 1000, 10000, 100000):
    print('   x_o=%-7d' % xo, '  '.join('l=%d: % .4f' % (l, D(l, 1, s3, s3 * xo) * xo) for l in (2, 3)))
print('(c2) k0 -> inf at eta_*=2, eta_o=8: (Delta - SW/Doppler lead)*x_o^2, Delta*k0*chi_*')
for k in (5, 20, 80, 320, 1280):
    for l in (2, 3):
        d = D(l, k, 2, 8); xs, xo, y = k * 2 / s3, k * 8 / s3, k * 6
        lead = -(mp.cos(xs) * T.jl(l, y) + s3 * mp.sin(xs) * mp.diff(lambda t: T.jl(l, t), y))
        print('   k0=%-5d l=%d  Delta=% .6e  lead=% .6e  (D-lead)*x_o^2=% .4f  Delta*k0*chi_*=% .4f' % (k, l, d, lead, (d - lead) * xo**2, d * y))
