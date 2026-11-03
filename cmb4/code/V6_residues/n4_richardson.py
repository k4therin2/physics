"""n4 (V6): residues of the monochromatic mu-integrand from A3's engine by polynomial (Richardson) extrapolation of
y_c,l(e) = e * Delta^kappa_mu,c,l(mu = 1-e), e = 0.02 ... 0.000625; also a fit allowing an e^{1/2} term (i.e. (1-mu)^{-1/2} in the
integrand) to bound it; comparison with the exact gauge responses n1 (the predicted residues)."""
import pickle, glob
import numpy as np
es = [0.02, 0.01, 0.005, 0.0025, 0.00125, 0.000625]
D = {e: pickle.load(open('n2_e%.6g_s+1.pkl' % e, 'rb')) for e in es}
n1 = pickle.load(open('n1_resp_mu0.31.pkl', 'rb'))['resp']
E = np.array(es); names = list(D[es[0]].keys())
def fit(y, deg, sqrt=False):
    cols = [E**j for j in range(deg + 1)]
    if sqrt: cols.insert(1, np.sqrt(E))
    A = np.stack(cols, 1)
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    return c
R = {}; Rlow = {}; S = {}
for nm in names:
    Y = np.array([e * D[e][nm].real for e in es])          # (6, 9)
    R[nm] = np.array([fit(Y[:, l], 4)[0] for l in range(9)])
    Rlow[nm] = np.array([fit(Y[:, l], 3)[0] for l in range(9)])
    S[nm] = np.array([fit(Y[:, l], 3, sqrt=True)[1] for l in range(9)])
big = [nm for nm in names if np.max(np.abs(R[nm])) > 1e-9]
print('channels with |residue| > 1e-9:')
for nm in big:
    print('  %-58s R  = %s' % (nm[:58], ' '.join('%+.10f' % x for x in R[nm])))
    print('  %-58s n1 = %s   |R-n1| max %.1e, deg3-deg4 %.1e' % ('', ' '.join('%+.10f' % x for x in n1[nm].real), np.max(np.abs(R[nm] - n1[nm].real)), np.max(np.abs(R[nm] - Rlow[nm]))))
others = max(np.max(np.abs(R[nm])) for nm in names if nm not in big)
print('max |residue| of all other channels: %.1e' % others)
tot = sum(R.values()); scale = max(np.max(np.abs(R[nm])) for nm in names)
print('SUM of residues (l=0..8):', ' '.join('%+.1e' % x for x in tot), '  max %.2e, relative to largest channel residue %.3f: %.1e' % (np.max(np.abs(tot)), scale, np.max(np.abs(tot)) / scale))
Stot = sum(S.values())
print('(1-mu)^(-1/2) coefficient of the TOTAL (fit R + S e^1/2 + poly):', ' '.join('%+.1e' % x for x in Stot))
print('   largest single-channel sqrt coefficient: %.1e' % max(np.max(np.abs(S[nm])) for nm in names))
# finite part F_l = lim [Delta_total - R_total/e]: total integrand
Ytot = np.array([sum(D[e][nm].real for nm in names) for e in es])
F = np.array([fit(Ytot[:, l], 4)[0] for l in range(9)])
Fs = np.array([fit(Ytot[:, l] - tot[l] / E, 3, sqrt=False)[0] for l in range(9)])
print('total integrand at e=0.000625:', np.round(Ytot[-1], 6))
print('finite limit F_l (l=0..8) = lim_{mu->1} Delta^kappa_mu,l:', np.round(Fs, 6))
pickle.dump({'R': R, 'S': S, 'sum': tot, 'F': Fs}, open('n4_results.pkl', 'wb'))
