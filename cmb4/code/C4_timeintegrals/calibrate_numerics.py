"""calibrate_numerics.py -- double-precision closed form vs ODE series vs mpmath reference, as a function of z=(k+q+p)eta/sqrt3."""
import sys, time, numpy as np, sympy as sp, mpmath as mp
sys.path.insert(0, '.')
import second_order_closed as SO
from second_order_closed import *
t0 = time.time()
kk, qq, mm = 1.0, 0.8, 0.3
pp = float(np.sqrt(kk**2 + qq**2 - 2 * kk * qq * mm))
Pl = SO._lam('P', P_expr); dPl = SO._lam('dP', sp.diff(P_expr, eta)); Il = SO._lam('Ih', h2_I_expr)
print("lambdified in %.0fs" % (time.time() - t0))
print("   z      closed(P)       series(P)      ref(P)          relerr closed   relerr series | h2: closed  series ref  rel  rel")
for z in [0.5, 1, 2, 3, 4, 5, 6, 8, 10, 14, 20, 30]:
    ee = z * np.sqrt(3) / (kk + qq + pp)
    ref = float(sp.re(sp.N(P_expr.subs({k: kk, q: qq, p: pp, eta: ee}), 40)))
    refh = float(sp.re(sp.N(h2_I_expr.subs({k: kk, q: qq, p: pp, eta: ee}), 40)))
    c = float(Pl(kk, qq, pp, ee)); s = float(SO._series_phi(kk, qq, pp, ee)[0])
    ch = float(Il(kk, qq, pp, ee)); sh = float(SO._series_h(kk, qq, pp, ee))
    print("%5.1f  %+.10e  %+.10e  %+.10e  %.1e  %.1e | %+.6e %+.6e %+.6e %.1e %.1e" % (z, c, s, ref, abs(c / ref - 1), abs(s / ref - 1), ch, sh, refh, abs(ch / refh - 1), abs(sh / refh - 1)))
# derivative check
ee = 3.0; refd = float(sp.re(sp.N(sp.diff(P_expr, eta).subs({k: kk, q: qq, p: pp, eta: ee}), 40)))
print("dP: closed %.12e series %.12e ref %.12e" % (dPl(kk, qq, pp, ee), SO._series_phi(kk, qq, pp, ee)[1], refd))
# full API vs sympy at a generic point
vals = {k: kk, q: qq, p: pp, eta: 2.5}
for nm, f, e in (('Phi2', Phi2, Phi2_expr), ('Psi2', Psi2, Psi2_expr), ('delta2', delta2, delta2_expr), ('deltaT2', deltaT2, deltaT2_expr), ('V2', V2, V2_expr), ('h2_I', h2_I, h2_I_expr)):
    a = f(kk, qq, mm, 2.5); b = complex(sp.N(e.subs(vals), 30))
    print("  %s: numeric %s  sympy %s" % (nm, a, b))
print("  omega2:", omega2(kk, qq, mm, 2.5), complex(sp.N(omega2(k, q, mu, eta).subs({k: kk, q: qq, mu: mm, eta: 2.5}), 20)))
print("  vT2:", vT2(kk, qq, mm, 2.5), complex(sp.N(vT2(k, q, mu, eta).subs({k: kk, q: qq, mu: mm, eta: 2.5}), 20)))
# vectorised timing
N = 20000
rng = np.random.default_rng(0)
K_ = rng.uniform(0.3, 2, N); Qv = rng.uniform(0.05, 3, N) * K_; M_ = rng.uniform(-1, 1, N); E_ = rng.uniform(0.1, 30, N) / K_
t1 = time.time(); v = Phi2(K_, Qv, M_, E_); print("Phi2 on %d points: %.2fs, finite: %s" % (N, time.time() - t1, np.all(np.isfinite(v))))
t1 = time.time(); v = h2_I(K_, Qv, M_, E_); print("h2_I on %d points: %.2fs, finite: %s" % (N, time.time() - t1, np.all(np.isfinite(v))))
