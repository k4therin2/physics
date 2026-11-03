"""analyze.py -- tables from cl22_<tag>.pkl: C11, C22, (1/4)C22, convergence between tags, group decomposition, K->0 integrand."""
import sys, os, pickle
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from typeB_check import CBB
H = os.path.dirname(os.path.abspath(__file__))
ld = lambda t: pickle.load(open(os.path.join(H, 'cl22_%s.pkl' % t), 'rb'))
main, ref = sys.argv[1], (sys.argv[2] if len(sys.argv) > 2 else None)
r = ld(main); rr = ld(ref) if ref else None
print('== %s: k0=%g eta_*=%g eta_o=%g nK=%d nth=%d   (units sigma^2 for C11, sigma^4 for C22)' % (main, r['k0'], r['eta_s'], r['eta_o'], r['nK'], r['nth']))
print(' l   C11              C22              C22/4            rel.diff vs %s   C22/4 / C11' % ref)
for l in range(2, 21):
    d = abs(r['C22'][l] / rr['C22'][l] - 1) if rr else np.nan
    print('%3d  %.9e  %.9e  %.9e  %.1e   %.3e' % (l, r['C11'][l], r['C22'][l], r['C22'][l] / 4, d, r['C22'][l] / 4 / r['C11'][l]))
gs = r['groups']; G = r['G']
print('\nGroup diagonal contributions G_gg(l) / C22(l)  (cross terms = remainder):')
print('  l ' + ' '.join('%9s' % g_[:9] for g_ in gs) + '  cross')
for l in (2, 3, 4, 6, 8, 10, 15, 20):
    dg = np.array([G[i, i, l] for i in range(len(gs))]) / r['C22'][l]
    print('%3d ' % l + ' '.join('%9.3f' % x for x in dg) + '  %7.3f' % (1 - dg.sum()))
for i, g_ in enumerate(gs): print('   col %d = %s' % (i, g_))
cb = CBB(r['k0'], r['eta_s'], r['eta_o'], 20)
i = gs.index('B: L1^2 - 2<m1>m1')
print('\nType-B 3j check, max rel diff l=2..20:', max(abs(cb[l] / G[i, i, l] - 1) for l in range(2, 21)), ' (l=2..12: %.1e)' % max(abs(cb[l] / G[i, i, l] - 1) for l in range(2, 13)))
Ks = sorted(r['integrand'])
print('\nK->0: sum_m |K_lm|^2 at the smallest K nodes (l = 0,1,2,3,4):')
for K in Ks[:4]:
    print('  K=%.4f  ' % K + ' '.join('%.3e' % x for x in r['integrand'][K][:5]))
