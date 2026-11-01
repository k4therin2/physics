"""Convergence / error-budget study on the observable at eps = 0.03 (x_o = 0) and eps=0.03 at x_o = pi/2."""
import numpy as np, sys, time, json
from planar import PlanarRun
ic = int(sys.argv[1]) if len(sys.argv) > 1 else 4
eps = 0.03
base = dict(N=48, ic_order=ic, nmu=48, rtol=1e-13, atol=1e-16, t_i=0.01, cH=-0.5, verbose=False)
def run(**kw):
    p = dict(base); p.update(kw)
    t0 = time.time(); r = PlanarRun(eps, **p).run()
    return r, time.time() - t0
ref, dt = run()
print(f"reference: N=48 ic={ic} rtol=1e-13 cH=-0.5 t_i=0.01 nmu=48  ({dt:.1f}s)  t_o={ref.t_o:.14f} T_o={ref.T_o:.14e}")
aref = ref.alm(ref.Theta, 24)
print("a_l0 (l=0..8):", np.array2string(aref[:9], precision=12))
for kw in [dict(N=32), dict(N=64), dict(rtol=1e-12, atol=1e-15), dict(rtol=3e-14, atol=1e-17), dict(cH=0.0), dict(cH=-1.0), dict(t_i=0.005), dict(t_i=0.02)]:
    r, dt = run(**kw)
    print(f"{str(kw):40s} ({dt:5.1f}s) max|dTheta|={np.max(np.abs(r.Theta-ref.Theta)):.2e}  dt_o={r.t_o-ref.t_o:+.2e}  dT_o/T_o={r.T_o/ref.T_o-1:+.2e}  max|d a_l0|={np.max(np.abs(r.alm(r.Theta,24)-aref)):.2e}  exactS={r.exact_formula_residual:.1e}")
for nmu in [32, 64, 96]:
    r, dt = run(nmu=nmu)
    a = r.alm(r.Theta, 24)
    print(f"nmu={nmu:3d} ({dt:5.1f}s) max|d a_l0|(l<=24)={np.max(np.abs(a-aref)):.2e}  l<=20: {np.max(np.abs(a-aref)[:21]):.2e}  monopole diff={r.monopole/ref.monopole-1:+.2e}")
