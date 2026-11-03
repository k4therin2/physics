import sys, time, numpy as np, pickle
sys.path.insert(0, '.')
from pw_map import *
from n1ref import load, taylor
tag = sys.argv[1] if len(sys.argv) > 1 else 'to8_xo0'
xL = {'to8_xo0': 0.0, 'to4_xo0': 0.0, 'to8_xo0.5pi': np.pi/2, 'to4_xo0.5pi': np.pi/2}[tag]
etao = 8.0 if tag.startswith('to8') else 4.0
N = int(sys.argv[2]) if len(sys.argv) > 2 else 48
pw = PlaneWave(); d, r = load_c2(); mb = MapBuilder(pw, d, r)
res = mb.run(xL, etao, 2.0, nmu=96, N=N)
n1 = load(tag)
assert np.allclose(n1['mu'], res['mu'])
cS = taylor(n1['eps'], -n1['S'], 0.08, 4, 4)
Lam = [cS[1], 2 * cS[2], 6 * cS[3]]
for k in range(3):
    diff = res['L'][k] - Lam[k]
    print('L%d: max|A1-N1| = %.3e   (max|L| %.3e)  a_l0 diff:' % (k + 1, np.max(np.abs(diff)), np.max(np.abs(Lam[k]))),
          np.round(alm(res['mu'], res['w'], diff, 6), 8))
print('obs', {k: (float(v) if not isinstance(v, tuple) else tuple(float(x) for x in v)) for k, v in res['obs'].items()})
pickle.dump(dict(res={k: v for k, v in res.items() if k not in ('obs',)}, Lam=Lam), open('cmp_%s.pkl' % tag, 'wb'))
