"""compare.py tag [N nmu] -- A1 map vs N1 coefficient maps (Th1, Th2, Th3) and a_l0; also T_o, <T_obs>."""
import sys, numpy as np, pickle, json, os
sys.path.insert(0, '.')
from pw_map import *
from n1ref import load, taylor, N1
tag = sys.argv[1] if len(sys.argv) > 1 else 'to8_xo0'
N = int(sys.argv[2]) if len(sys.argv) > 2 else 64
fields = sys.argv[3] if len(sys.argv) > 3 else None
xL = np.pi / 2 if '0.5pi' in tag else 0.0
etao = 8.0 if tag.startswith('to8') else 4.0
pw = PlaneWave(fields); d, r = load_c2(); mb = MapBuilder(pw, d, r, verbose=False)
res = mb.run(xL, etao, 2.0, nmu=96, N=N, verbose=False)
c = np.load(os.path.join(N1, 'coeff_%s.npz' % tag))
assert np.allclose(c['mu'], res['mu'])
out = {'tag': tag, 'N': N, 'rows': []}
print('== %s  (x_o=%.4f, etabar_o=%g, N=%d, fields=%s)' % (tag, xL, etao, N, fields or 'pw_fields.pkl'))
for k in range(3):
    A = alm(res['mu'], res['w'], res['Th'][k], 12)
    B = c['alm_Th%d' % (k + 1)][:13]
    dm = np.max(np.abs(res['Th'][k] - c['Th%d' % (k + 1)]))
    print('Th%d: max|map diff| = %.2e   (N1 extraction uncertainty %.1e)' % (k + 1, dm, float(c['spread_Th%d' % (k + 1)])))
    for l in range(0, 13):
        if abs(B[l]) > 1e-6 or abs(A[l]) > 1e-6:
            print('   l=%2d  A1 % .10e  N1 % .10e  diff % .2e' % (l, A[l], B[l], A[l] - B[l]))
        out['rows'].append((k + 1, l, float(A[l]), float(B[l])))
    out['Th%d_maxdiff' % (k + 1)] = float(dm)
n1 = load(tag)
cT = taylor(n1['eps'], n1['T_o'], 0.1, 4, 5); cM = taylor(n1['eps'], n1['monopole'], 0.1, 4, 5)
o = res['obs']; T = [float(x) for x in o['T_o']]
a1, a2, a3 = res['mono']
Mbar = [1, a1, a2 / 2]
Tobs = [T[0], T[1] + T[0] * a1, T[2] + T[1] * a1 + T[0] * a2 / 2]
print('T_o     coeffs A1', T, ' N1', [float(cT[j]) for j in range(3)])
print('<T_obs> coeffs A1', Tobs, ' N1', [float(cM[j]) for j in range(3)])
print('sky monopoles <m1>,<m2>,<m3> =', res['mono'])
print('Poisson-gauge clock: deta_o1=%.10f deta_o2=%.10f dxo1=%.3e dxo2=%.3e' % (o['deta_o1'], o['deta_o2'], o['dxo1'], o['dxo2']))
out.update(T_o=(T, [float(cT[j]) for j in range(3)]), Tobs=(Tobs, [float(cM[j]) for j in range(3)]), mono=[float(x) for x in res['mono']])
json.dump(out, open('compare_%s_N%d%s.json' % (tag, N, '' if fields is None else '_' + os.path.basename(fields)), 'w'), indent=1)
np.savez('map_%s_N%d.npz' % (tag, N), mu=res['mu'], w=res['w'], Th1=res['Th'][0], Th2=res['Th'][1], Th3=res['Th'][2],
         L1=res['L'][0], L2=res['L'][1], L3=res['L'][2])
