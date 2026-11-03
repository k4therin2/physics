"""check C2's conservation rules d_eta dT_n, d_eta v_n with A1's numeric jets (tests field naming/jets/fields)."""
import sys, numpy as np, sympy as sp
sys.path.insert(0, '.')
from pw_map import *
sys.path.insert(0, C2)
from fields import build
fd = build(verbose=False)
rules = fd['rules']
pw = PlaneWave()
etas = list(np.linspace(1.0, 7.5, 7))
tab = Tables(pw, etas)
X = np.linspace(0.1, 6.0, 9)[:, None] + 0 * np.array(etas)[None, :]
ctx = Ctx(tab, np.arange(len(etas))[None, :], X, {}, {})
for lhs, rhs in rules.items():
    nm = str(lhs.expr.func)
    if nm.startswith(ZERO_FIELDS): 
        r = planar_reduce(rhs)
        print(nm, 'zero field; reduced rhs ->', 'ZERO' if r == 0 else r); continue
    e = NpExpr(lhs - rhs)
    val = e(ctx)
    print('%-10s residual max %.2e   (|lhs| max %.2e)' % (nm, np.max(np.abs(val)), np.max(np.abs(NpExpr(lhs)(ctx)))))
Xs = np.linspace(0, 2*np.pi, 16, endpoint=False)[:, None] + 0*np.array(etas)[None, :]
ctx2 = Ctx(tab, np.arange(len(etas))[None, :], Xs, {}, {})
for lhs, rhs in rules.items():
    nm = str(lhs.expr.func)
    if nm not in ('dT2',): continue
    val = NpExpr(lhs - rhs)(ctx2)
    f = np.fft.rfft(val, axis=0) / 16
    print(nm, 'mode amplitudes (rows m=0..3) vs eta', etas); print(np.round(np.real(f[:4]), 6)); print(np.round(np.imag(f[:4]), 6))
    print(planar_reduce(rhs))
