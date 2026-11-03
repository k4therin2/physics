"""X1: superhorizon ratios R_l(x) = [1/4 C22 + 1/3 Re C13]/C11 (sigma = 1) from our A2/A3 evaluators on eta_*=2, eta_o=8,
x = k0 eta_o in {0.4, 0.2, 0.1, 0.05, 0.025, 0.0125}; fits R = a + b x^2 + c x^4 and Richardson extrapolation.
Adam's prediction (adam_sh_pred.py): l=2: 1/4C22/C11 = 17/3 - 1187/360 x^2, 1/3ReC13/C11 = 4 - 257/60 x^2, sum 29/3 - 2729/360 x^2;
l=3: 1/4C22/C11 = 148/4 = 37, 1/3ReC13/C11 = 27/3 = 9, sum 46 (leading order only)."""
import pickle, glob, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); W4 = os.path.join(HERE, '..', 'W4_final')
def load(kind, x):
    for d in (HERE, W4):
        f = sorted(glob.glob(os.path.join(d, 'sh%s_%g_*.pkl' % (kind, x))))
        if f: return pickle.load(open(f[0], 'rb'))
xs = [0.4, 0.2, 0.1, 0.05, 0.025, 0.0125]
rows = {}
for x in xs:
    a, b = load('22', x), load('13', x)
    for l in (2, 3, 4):
        C11, Q22, Q13 = a['C11'][l], a['C22'][l] / 4, b['C13'][l].real / 3
        C11b = 4 * np.pi * b['Delta'][l]**2
        rows.setdefault(l, []).append((x, Q22 / C11, Q13 / C11, (Q22 + Q13) / C11, C11 / C11b - 1))
pred = {2: (17/3, -1187/360, 4, -257/60), 3: (37, None, 9, None)}
for l in (2, 3, 4):
    print('l=%d' % l)
    print('   x        1/4C22/C11      1/3ReC13/C11     sum/C11         (C11 A2 vs A3 Delta)-1')
    R = np.array(rows[l])
    for r in R: print('  %-7g  %-15.10f %-15.10f %-15.10f %.1e' % tuple(r))
    for col, nm in ((1, '22'), (2, '13'), (3, 'sum')):
        # Richardson on the three smallest x (ratio 2 steps, error x^2, x^4)
        f1, f2, f3 = R[-3, col], R[-2, col], R[-1, col]
        r1 = f2 + (f2 - f1) / 3; r2 = f3 + (f3 - f2) / 3; rich = r2 + (r2 - r1) / 15
        A = np.vstack([np.ones(len(xs)), R[:, 0]**2, R[:, 0]**4]).T
        coef = np.linalg.lstsq(A, R[:, col], rcond=None)[0]
        A4 = A[1:]; coef4 = np.linalg.lstsq(A4, R[1:, col], rcond=None)[0]
        print('   %-4s Richardson(3 smallest x) = %.6f   LSQ a+bx^2+cx^4 (6 pts): a=%.6f b=%.4f c=%.3f  (5 pts: a=%.6f b=%.4f)' %
              (nm, rich, coef[0], coef[1], coef[2], coef4[0], coef4[1]))
    if l in pred:
        p = pred[l]
        print('   Adam: 22: %.6f%s   13: %.6f%s   sum: %.6f%s' % (p[0], '' if p[1] is None else ' %+.4f x^2' % p[1], p[2],
              '' if p[3] is None else ' %+.4f x^2' % p[3], p[0] + p[2], '' if p[1] is None else ' %+.4f x^2' % (p[1] + p[3])))
