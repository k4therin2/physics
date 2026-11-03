import sys, time, numpy as np
sys.path.insert(0, '.')
from pw_map import *
t0=time.time()
pw = PlaneWave(); print('pw loaded %.1f'%(time.time()-t0))
d, r = load_c2()
mb = MapBuilder(pw, d, r)
obs = observer(pw, 0, 8, 0.01)
print({k: (v if not isinstance(v, tuple) else tuple(float(x) for x in v)) for k, v in obs.items()})
res = mb.run(0.0, 8.0, 2.0, nmu=48, N=48, obs=obs)
for k in range(3):
    print('Th%d a_l0:'%(k+1), alm(res['mu'], res['w'], res['Th'][k], 8))
