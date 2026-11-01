"""Production amplitude scan.  Usage: python3 run_scan.py <t_obar> <x_o_over_pi> [tag]
Runs eps = +-{0.01,...,0.10}, saves output/scan_<tag>.npz with all maps and diagnostics."""
import numpy as np, sys, time, json
from planar import PlanarRun
t_obar = float(sys.argv[1]); x_o = float(sys.argv[2])*np.pi
tag = sys.argv[3] if len(sys.argv) > 3 else f"to{t_obar:g}_xo{sys.argv[2]}pi"
ic_order = int(sys.argv[4]) if len(sys.argv) > 4 else 6
t_i = float(sys.argv[5]) if len(sys.argv) > 5 else 0.01
eps_list = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10]
eps_all = []
for e in eps_list:
    eps_all += [e, -e]
P = dict(N=48, ic_order=ic_order, nmu=96, rtol=1e-13, atol=1e-16, t_i=t_i, cH=-0.5, t_star=2.0, t_obar=t_obar, x_o=x_o, verbose=True)
res = dict(eps=[], T_obs=[], T_o=[], t_o=[], t_init=[], t_e=[], x_e=[], S=[], Theta=[], monopole=[], alm=[], exactS_res=[], H=[], Hrel=[], M=[], Mrel=[], T_obs_null=[])
t0 = time.time()
for eps in eps_all:
    r = PlanarRun(eps, **P).run()
    c = r.constraints(r.t_o, r.y_o)
    res['eps'].append(eps); res['T_obs'].append(r.T_obs); res['T_o'].append(r.T_o); res['t_o'].append(r.t_o); res['t_init'].append(r.t_init)
    res['t_e'].append(r.t_e); res['x_e'].append(r.x_e); res['S'].append(r.S); res['Theta'].append(r.Theta); res['monopole'].append(r.monopole)
    res['alm'].append(r.alm(r.Theta, 30)); res['exactS_res'].append(r.exact_formula_residual); res['T_obs_null'].append(r.T_obs_null)
    res['H'].append(c[0]); res['Hrel'].append(c[1]); res['M'].append(c[2]); res['Mrel'].append(c[3])
    print(f"   done eps={eps:+.3f}  ({time.time()-t0:.0f}s)", flush=True)
out = {k: np.array(v) for k, v in res.items()}
out.update(mu=r.mu, wmu=r.wmu, params=json.dumps({k: (v if not isinstance(v, float) else v) for k, v in P.items()}), t_star=2.0, t_obar=t_obar, x_o=x_o, tag=tag)
np.savez(f"output/scan_{tag}.npz", **out)
print("saved", f"output/scan_{tag}.npz", f"total {time.time()-t0:.0f}s")
