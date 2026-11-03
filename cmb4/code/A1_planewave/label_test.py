"""Effect of the label convention (x_L at eta->0 [default] vs SPEC-literal x_L at eta_init) for x_o = pi/2, etao = 8."""
import sys, numpy as np, json
sys.path.insert(0, '.')
from pw_map import *
pw = PlaneWave(); d, r = load_c2(); mb = MapBuilder(pw, d, r, verbose=False)
c = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'N1_planar', 'output', 'coeff_to8_xo0.5pi.npz'))
for lower0 in (True, False):
    res = mb.run(np.pi / 2, 8.0, 2.0, nmu=96, N=64, lower0=lower0, verbose=False)
    print('lower0=%s dxo1=%.12f' % (lower0, float(res['obs']['dxo1'])))
    for k in range(3):
        A = alm(res['mu'], res['w'], res['Th'][k], 3)
        print('  Th%d a_l0 (l=0..3) diff vs N1:' % (k + 1), A - c['alm_Th%d' % (k + 1)][:4])
