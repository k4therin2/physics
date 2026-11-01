import numpy as np, time, sys
from planar import PlanarRun
ic = int(sys.argv[1]) if len(sys.argv) > 1 else 4
t0 = time.time()
r0 = PlanarRun(0.0, N=32, ic_order=ic, nmu=16).run()
print("FRW: t_o - t_obar =", r0.t_o - r0.t_obar, " max|T_obs/T_o - 1| =", np.max(np.abs(r0.ratio - 1)), " max|S| =", np.max(np.abs(r0.S)), f" ({time.time()-t0:.1f}s)")
t0 = time.time()
r = PlanarRun(0.01, N=32, ic_order=ic, nmu=16)
r.evolve_forward()
# constraints along the forward solution
for i in range(0, len(r.sol_forward.t), max(1, len(r.sol_forward.t)//8)):
    tt = r.sol_forward.t[i]; c = r.constraints(tt, r.sol_forward.y[:-1, i])
    print(f"  t={tt:8.4f}  |H|={c[0]:.2e} rel={c[1]:.2e}  |M|={c[2]:.2e} rel={c[3]:.2e}")
r.trace_rays()
print(f"({time.time()-t0:.1f}s)  Theta:", r.Theta[:4])
