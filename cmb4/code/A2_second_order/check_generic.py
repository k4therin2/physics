"""check_generic.py -- CHECKS on the full K_2 for generic (non-collinear) momenta:
 (A) shear form == direct form pointwise in n (equivalence modulo 2nd-order conservation: tests Phi2,Psi2,dT2,V2,v2T,om2,h2 at K != 0);
 (B) symmetry K_2(k1,k2) = K_2(k2,k1);  (C) reality K_2(-k1,-k2;n) = K_2(k1,k2;n)^*;
 (D) constant mode (i) and gradient (ii): l>=2 power of K_2(k1 = s e, k2; n) as s -> 0 (full map vs channels)."""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kernel2 import Geometry, MapKernels
from sphere import Sky
g = Geometry(2.0, 8.0)
sky = Sky(40)
nv = sky.nv
k1 = np.array([0.7, -0.2, 0.4]); k2 = np.array([-0.3, 0.9, 0.5])
Kd, chd = MapKernels(g, k1, k2, nv, form='direct').K2()
Ks, chs = MapKernels(g, k1, k2, nv, form='shear').K2()
print('(A) max|K2_direct - K2_shear| / max|K2| =', np.max(np.abs(Kd - Ks)) / np.max(np.abs(Kd)))
Kt, _ = MapKernels(g, k2, k1, nv, form='shear').K2()
print('(B) symmetry  max|K2(k1,k2)-K2(k2,k1)| =', np.max(np.abs(Ks - Kt)))
Km, _ = MapKernels(g, -k1, -k2, nv, form='shear').K2()
print('(C) reality   max|K2(-k1,-k2)-conj K2(k1,k2)| =', np.max(np.abs(Km - np.conj(Ks))))
# (D)
L = 30
e = np.array([0.3, 0.5, -0.8]); e /= np.linalg.norm(e)
k2 = np.array([-0.3, 0.9, 0.5])
print('(D) squeezed k1 = s e, |k2| = %.3f ; l>=2 rms of the full K_2 and of selected channels; l=0,1 separately' % np.linalg.norm(k2))
for form in ('shear', 'direct'):
    for s in (1e-1, 3e-2, 1e-2, 3e-3, 1e-3):
        K, ch = MapKernels(g, s * e, k2, nv, form=form).K2()
        P = sky.power(K, L)
        big = max(ch, key=lambda c: np.sum(sky.power(ch[c], L)[2:]))
        Pc = sky.power(ch[big], L)
        print('  %-6s s=%.0e  full: l>=2 %.3e  l=1 %.3e  l=0 %.3e   | largest channel (%s): l>=2 %.3e' % (
            form, s, np.sqrt(np.sum(P[2:])), np.sqrt(P[1]), np.sqrt(P[0]), big[:40], np.sqrt(np.sum(Pc[2:]))))
