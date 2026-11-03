import sys, numpy as np
sys.path.insert(0, '.')
from pw_map import *
pw = PlaneWave(); d, r = load_c2()
mb = MapBuilder(pw, d, r)
res = mb.run(0.0, 8.0, 2.0, nmu=32, N=48)
d2 = dict(d); d2['L'] = d['L_direct']
mbd = MapBuilder(pw, d2, r)
resd = mbd.run(0.0, 8.0, 2.0, nmu=32, N=48, obs=res['obs'])
for k in range(3):
    print('L%d shear vs direct: max diff %.3e  (max %.3e)' % (k+1, np.max(np.abs(res['L'][k]-resd['L'][k])), np.max(np.abs(res['L'][k]))))
res96 = mbd.run(0.0, 8.0, 2.0, nmu=96, N=48, obs=res['obs'])
from n1ref import load, taylor
n1 = load('to8_xo0'); cS = taylor(n1['eps'], -n1['S'], 0.08, 4, 4)
Lam = [cS[1], 2 * cS[2], 6 * cS[3]]
for k in range(3):
    print('direct L%d - N1: max %.3e' % (k+1, np.max(np.abs(res96['L'][k]-Lam[k]))), np.round(alm(res96['mu'], res96['w'], res96['L'][k]-Lam[k], 6), 6))
# group breakdown at second order (both forms) at a few mu
for nm, rr in (('shear', mb.run(0.0, 8.0, 2.0, nmu=96, N=48, obs=res['obs'], verbose=False)), ('direct', res96)):
    print('==', nm)
    for (m, g), v in rr['groups'].items():
        if m == 2: print('  %-70s a20=% .6f a00=% .6f' % (g, alm(rr['mu'], rr['w'], v, 2)[2], alm(rr['mu'], rr['w'], v, 2)[0]))
print('N1 L2 a00,a20', alm(res96['mu'], res96['w'], Lam[1], 2))
