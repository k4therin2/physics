"""check_gradient.py -- CHECK (ii): O(k1) squeezed term vs V2's prediction (report_V2 §3, eta_i -> 0):
   K_2(k1,k2;n) = -D K_1(k2;n) - S_{k1} K_1(k2;n) + O(k1^2),  D = k2.d_k2,  S_p = (k2.d)(p.d) - (1/2)(p.k2) Lap_k2.
With D, S evaluated on the ENGINE's K_1 by finite differences (Hessian in 3D k-space).  Compared multipole-wise (l >= 1)
via the l-power of the residual  R = [K_2(s e,k2) + D K_1]/s + S_e K_1, isolated as the k1-odd part [K_2(s e) - K_2(-s e)]/(2s)."""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kernel2 import Geometry, MapKernels, Slot
from sphere import Sky
g = Geometry(2.0, 8.0); sky = Sky(40); nv = sky.nv; L = 30
e = np.array([0.3, 0.5, -0.8]); e /= np.linalg.norm(e)
for k2 in (np.array([-0.3, 0.9, 0.5]), np.array([1.2, -2.0, 0.7])):
    mk = MapKernels(g, k2, k2, nv)
    F = lambda kv: mk.K1(Slot(g, kv, nv))
    h = 2e-3 * np.linalg.norm(k2)
    F0 = F(k2); I3 = np.eye(3)
    grad = [(F(k2 + h * I3[i]) - F(k2 - h * I3[i])) / (2 * h) for i in range(3)]
    Hs = {}
    for i in range(3):
        Hs[(i, i)] = (F(k2 + h * I3[i]) - 2 * F0 + F(k2 - h * I3[i])) / h**2
        for j in range(i + 1, 3):
            Hs[(i, j)] = Hs[(j, i)] = (F(k2 + h * (I3[i] + I3[j])) - F(k2 + h * (I3[i] - I3[j])) - F(k2 - h * (I3[i] - I3[j])) + F(k2 - h * (I3[i] + I3[j]))) / (4 * h * h)
    DK = sum(k2[i] * grad[i] for i in range(3))
    # S_p F = k_i p_j H_ij + (p.grad)  [from (k.d)(p.d)F = k_i d_i(p_j d_j F) = k_i p_j H_ij] ; careful: (k.d)(p.d) acts as k_i p_j H_ij only
    SK = sum(k2[i] * e[j] * Hs[(i, j)] for i in range(3) for j in range(3)) - 0.5 * (e @ k2) * sum(Hs[(i, i)] for i in range(3))
    pred = -SK
    s = 2e-3                                      # odd part in k1 isolates O(k1) (O(k1^0) cancels, error O(s^2))
    O1x = (MapKernels(g, s * e, k2, nv).K2()[0] - MapKernels(g, -s * e, k2, nv).K2()[0]) / (2 * s)
    P_pred = sky.power(pred, L); P_res = sky.power(O1x - pred, L); P_num = sky.power(O1x, L)
    print('|k2|=%.3f  O(k1) coefficient (unit e):  l=1: num %.5e pred %.5e resid %.2e | l>=2 rms: num %.5e pred %.5e resid %.2e' % (
        np.linalg.norm(k2), np.sqrt(P_num[1]), np.sqrt(P_pred[1]), np.sqrt(P_res[1]), np.sqrt(P_num[2:].sum()), np.sqrt(P_pred[2:].sum()), np.sqrt(P_res[2:].sum())))
    print('     per l (2..8) resid/num:', ' '.join('%.1e' % np.sqrt(P_res[l] / P_num[l]) for l in range(2, 9)))
