"""Markdown summary of the extracted coefficients (a_l0 of Th1,Th2,Th3) for all tags -> output/summary.md"""
import numpy as np, json
lines = []
for tag in ['to8_xo0', 'to4_xo0', 'to8_xo0.5pi', 'to4_xo0.5pi']:
    c = np.load(f'output/coeff_{tag}.npz'); s = np.load(f'output/scan_{tag}.npz', allow_pickle=True)
    lines.append(f"\n### {tag}: t_* = {float(c['t_star'])}, t_obar = {float(c['t_obar'])}, x_o = {float(c['x_o']):.6f}")
    lines.append(f"runs: eps = +-{{0.01..0.10}}; max exact-formula residual {np.max(s['exactS_res']):.1e}; constraints (rel) H<{np.max(s['Hrel']):.1e}, M<{np.max(s['Mrel']):.1e}")
    lines.append("| eps | t_o | T_o | <T_obs> | a_10(Theta) |\n|---|---|---|---|---|")
    for i in np.argsort(s['eps']):
        lines.append(f"| {s['eps'][i]:+.2f} | {s['t_o'][i]:.12f} | {s['T_o'][i]:.12e} | {s['monopole'][i]:.12e} | {s['alm'][i][1]:+.3e} |")
    lines.append("\n| l | a_l0[Th1] | a_l0[Th2] | a_l0[Th3] |\n|---|---|---|---|")
    for l in range(0, 13):
        lines.append(f"| {l} | {c['alm_Th1'][l]:+.10e} | {c['alm_Th2'][l]:+.10e} | {c['alm_Th3'][l]:+.9e} |")
    lines.append(f"\nmap rms: Th1 {np.sqrt(0.5*np.sum(c['wmu']*c['Th1']**2)):.6e}, Th2 {np.sqrt(0.5*np.sum(c['wmu']*c['Th2']**2)):.6e}, Th3 {np.sqrt(0.5*np.sum(c['wmu']*c['Th3']**2)):.6e}; "
                 f"max|.|: Th1 {np.max(np.abs(c['Th1'])):.4e}, Th2 {np.max(np.abs(c['Th2'])):.4e}, Th3 {np.max(np.abs(c['Th3'])):.4e}")
    lines.append(f"extraction uncertainty (max over map, spread of fit subsets): Th1 {float(c['spread_Th1']):.1e}, Th2 {float(c['spread_Th2']):.1e}, Th3 {float(c['spread_Th3']):.1e}")
    lines.append(f"Th1 at mu=+1 (from +x), mu=0, mu=-1: {np.interp([1,0,-1], c['mu'], c['Th1'])}")
open('output/summary.md', 'w').write("\n".join(lines))
print("\n".join(lines))
