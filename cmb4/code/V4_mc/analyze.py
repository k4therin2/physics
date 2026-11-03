"""analyze.py -- combine mc_<tag>.pkl files: ensemble means with standard errors (direction sets are the independent units;
within a set the R amplitude draws share the momentum directions), vs A2's C22 (cl22_eo8_c.pkl) and 4 pi Delta_l^2."""
import sys, os, pickle, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
A2 = pickle.load(open(os.path.join(HERE, '..', 'A2_second_order', 'cl22_eo8_c.pkl'), 'rb'))
sets = []; Ns = set()
for t in sys.argv[1:]:
    d = pickle.load(open(os.path.join(HERE, 'mc_%s.pkl' % t), 'rb')); sets += d['sets']; Ns.add(d['N'])
assert len(Ns) == 1; N = Ns.pop(); S = len(sets); R = len(sets[0]['C22'])
def est(key, fac=1.0):
    m = np.array([np.atleast_2d(s[key]).mean(0) for s in sets]) * fac             # per-set means (S, L)
    return m.mean(0), m.std(0, ddof=1) / np.sqrt(S)
C11, e11 = est('C11'); C22, e22 = est('C22'); C22o, e22o = est('C22off', N / (N - 1))
Wk, eWk = est('wickoff', N / (N - 1)) if 'wickoff' in sets[0] else (None, None)
mean2 = np.array([s['mean2'] for s in sets]).mean(0)
print('N=%d modes, %d direction sets x %d amplitude draws = %d realisations' % (N, S, R, S * R))
print('<zeta^2>/sigma^2 = %.4f +- %.4f ; max normalisation-term check %.1e' % (np.mean([s['zeta0sq'] for s in sets]),
      np.std([s['zeta0sq'] for s in sets]) / np.sqrt(S), max(s['normchk'] for s in sets)))
print(' l | C11 MC/s^2          4pi D^2 (A2)  dev/err | (1/4)C22 full  +- err      | (1/4)C22 off*N/(N-1) +- err | A2 (1/4)C22  | dev_off/err  ratio_off')
for l in range(1, 13):
    a2 = A2['C22'][l] / 4
    print('%2d | %.5e+-%.1e %.5e %+5.1f | %.5e +- %.2e | %.5e +- %.2e | %.5e | %+5.1f  %.4f' % (l, C11[l], e11[l], A2['C11'][l],
          (C11[l] - A2['C11'][l]) / e11[l], C22[l] / 4, e22[l] / 4, C22o[l] / 4, e22o[l] / 4, a2, (C22o[l] / 4 - a2) / (e22o[l] / 4), C22o[l] / 4 / a2))
print()
print(' l | Wick-in-amplitudes (directions MC only) (1/4)C22 +- err | ratio to A2 | dev/err')
for l in range(1, 13):
    a2 = A2['C22'][l] / 4
    print('%2d | %.5e +- %.2e | %.4f | %+5.1f' % (l, Wk[l] / 4, eWk[l] / 4, Wk[l] / 4 / a2, (Wk[l] / 4 - a2) / (eWk[l] / 4)))
ls = np.repeat(np.arange(13), 2 * np.arange(13) + 1)
print('mean map check: sum_m <Th2_lm>^2/(2l+1) / C22_l (should -> 0):', ['%.1e' % (np.sum(mean2[ls == l]**2) / (2 * l + 1) / C22[l]) for l in range(1, 11)])
