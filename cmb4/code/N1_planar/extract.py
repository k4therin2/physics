"""Taylor-coefficient extraction in eps from the amplitude scan.
Theta(eps, mu) = sum_n c_n(mu) eps^n,  SPEC: Theta = eps Th1 + eps^2/2 Th2 + eps^3/6 Th3  => Th1=c1, Th2=2c2, Th3=6c3.
Even/odd split: E(e)=(Theta(e)+Theta(-e))/2 = c2 e^2 + c4 e^4 + ..., O(e)=(Theta(e)-Theta(-e))/2 = c1 e + c3 e^3 + ...
Polynomial least-squares fits on several subsets of amplitudes; spread = truncation/noise uncertainty.
Usage: python3 extract.py <tag>"""
import numpy as np, sys, json
from scipy.special import eval_legendre
tag = sys.argv[1]
d = np.load(f"output/scan_{tag}.npz", allow_pickle=True)
eps, Theta, mu, w = d['eps'], d['Theta'], d['mu'], d['wmu']
pos = {e: Theta[i] for i, e in enumerate(eps) if e > 0}
neg = {-e: Theta[i] for i, e in enumerate(eps) if e < 0}
amps = sorted(pos)
def fit(amp_subset, parity, nterms):
    """parity 'odd': O(e)/e = c1 + c3 e^2 + ...; 'even': E(e)/e^2 = c2 + c4 e^2 + ...  Returns coefficient maps (nterms, nmu)."""
    A = np.array([[a**(2*j) for j in range(nterms)] for a in amp_subset])
    if parity == 'odd':
        Y = np.array([(pos[a] - neg[a])/2/a for a in amp_subset])
    else:
        Y = np.array([(pos[a] + neg[a])/2/a**2 for a in amp_subset])
    coef, *_ = np.linalg.lstsq(A, Y, rcond=None)
    return coef
subsets = {
    'A: e<=0.04, 4pts, 4terms': ([0.01, 0.02, 0.03, 0.04], 4),
    'B: e<=0.06, 6pts, 5terms': ([0.01, 0.02, 0.03, 0.04, 0.05, 0.06], 5),
    'C: e<=0.10, 8pts, 6terms': (amps, 6),
    'D: e<=0.08, 7pts, 6terms': ([0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08], 6),
    'E: e=.02-.08 even, 4pts,4terms': ([0.02, 0.04, 0.06, 0.08], 4),
    'F: e<=0.05, 5pts, 4terms': ([0.01, 0.02, 0.03, 0.04, 0.05], 4),
}
def alm(f, lmax=30):
    return np.array([np.sqrt(4*np.pi*(2*l + 1))*0.5*np.sum(w*f*eval_legendre(l, mu)) for l in range(lmax + 1)])
results = {}
for name, (sub, nt) in subsets.items():
    co = fit(sub, 'odd', nt); ce = fit(sub, 'even', nt)
    results[name] = dict(c1=co[0], c3=co[1], c5=co[2], c2=ce[0], c4=ce[1])
names = list(results)
good = [n for n in names if not n.startswith('E')]   # E (no small-eps points) is a deliberately weak fit, excluded from the uncertainty
ref = 'D: e<=0.08, 7pts, 6terms'
summary = {}
print(f"=== {tag} ===  (rms over sphere; spread = max over fit subsets A,B,C,F of max|diff| vs subset {ref})")
for cn, fac, lab in [('c1', 1, 'Th1 = c1'), ('c2', 2, 'Th2 = 2 c2'), ('c3', 6, 'Th3 = 6 c3'), ('c4', 24, 'Th4 = 24 c4 (indicative)')]:
    r = results[ref][cn]*fac
    rms = np.sqrt(0.5*np.sum(w*r**2))
    spread = max(np.max(np.abs(results[n][cn]*fac - r)) for n in good if n != ref)
    spreads = {n: float(np.max(np.abs(results[n][cn]*fac - r))) for n in names}
    a = alm(r)
    print(f"{lab:28s} rms={rms:.6e}  max|.|={np.max(np.abs(r)):.6e}  spread={spread:.2e}  a_l0(l=0..6)=", np.array2string(a[:7], precision=8))
    for n in names:
        print(f"      {n:34s} max|diff|={spreads[n]:.2e}")
    summary[lab] = dict(rms=float(rms), spread=float(spread), a_l0=a.tolist(), map=r.tolist())
np.savez(f"output/coeff_{tag}.npz", mu=mu, wmu=w, Th1=results[ref]['c1'], Th2=2*results[ref]['c2'], Th3=6*results[ref]['c3'], Th4=24*results[ref]['c4'],
         alm_Th1=alm(results[ref]['c1']), alm_Th2=alm(2*results[ref]['c2']), alm_Th3=alm(6*results[ref]['c3']), alm_Th4=alm(24*results[ref]['c4']),
         spread_Th1=max(np.max(np.abs(results[n]['c1'] - results[ref]['c1'])) for n in good),
         spread_Th2=2*max(np.max(np.abs(results[n]['c2'] - results[ref]['c2'])) for n in good),
         spread_Th3=6*max(np.max(np.abs(results[n]['c3'] - results[ref]['c3'])) for n in good),
         fit_subset=ref, t_star=float(d['t_star']), t_obar=float(d['t_obar']), x_o=float(d['x_o']))
json.dump(summary, open(f"output/coeff_{tag}.json", 'w'), indent=1)
# diagnostics
print("per-run: eps, t_o, T_o, monopole, dipole a_10, exact-formula residual, |H|rel, |M|rel")
for i in range(len(eps)):
    print(f"  {eps[i]:+.3f}  t_o={d['t_o'][i]:.13f}  T_o={d['T_o'][i]:.13e}  <T_obs>={d['monopole'][i]:.13e}  a10={d['alm'][i][1]:+.2e}  expS={d['exactS_res'][i]:.1e}  H={d['Hrel'][i]:.1e} M={d['Mrel'][i]:.1e}")
