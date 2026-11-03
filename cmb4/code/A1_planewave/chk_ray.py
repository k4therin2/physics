import sys, numpy as np
sys.path.insert(0, '.')
from pw_map import *
pw = PlaneWave(); d, r = load_c2(); mb = MapBuilder(pw, d, r, verbose=False)
res = mb.run(0.0, 8.0, 2.0, nmu=8, N=48, verbose=False)
F = res['funs']; mu = res['mu']; etan = res['etan']
tab = Tables(pw, list(etan))
X = 0.0 + mu[:, None] * (8.0 - etan[None, :])
ctx = Ctx(tab, np.arange(49)[None, :], X, F, {})
P = ctx.jet('Phi1', 0, 0) + ctx.jet('Psi1', 0, 0)
n1 = mu[:, None]; n2 = np.sqrt(1 - mu**2)[:, None]
print('null cond first order max', np.max(np.abs(n1 * F['k1_1'] + n2 * F['k1_2'] + P)))
print('k1 sample', F['k1_1'][3, ::12], F['k1_2'][3, ::12])
