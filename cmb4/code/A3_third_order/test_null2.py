"""test_null2.py (A3): second-order null condition along the ray for the engine's k2 functional:
k2_z = -Psi2 - Phi2 - 2 Dx1.grad(Psi1+Phi1) - om2_z + |k1|^2 + 4 Phi1 k1_z + h2_zz/2   (frame n = z, fields on the background ray)."""
import numpy as np
from engine import Ctx, aadd, amul, ascale
rng = np.random.default_rng(3)
nv = rng.normal(size=(5, 3)); nv /= np.linalg.norm(nv, axis=1)[:, None]
ctx = Ctx([[[0.7, -0.2, 0.4]], [[-0.3, 0.9, 0.5]], [[-0.7, 0.2, -0.4]]], nv)
F = lambda nm, o, comp=None, der=(0, 0, 0, 0): ctx.field(('f', nm, comp, o, der), 'ray')
rhs = aadd(ascale(F('Psi', 2), -1), ascale(F('Phi', 2), -1))
for j in range(3):
    der = [0, 0, 0, 0]; der[j + 1] = 1; der = tuple(der)
    rhs = aadd(rhs, ascale(amul(ctx.functional('Dx', 1, j), aadd(F('Psi', 1, None, der), F('Phi', 1, None, der))), -2))
rhs = aadd(rhs, ascale(F('om', 2, 2), -1))
for j in range(3):
    rhs = aadd(rhs, amul(ctx.functional('k', 1, j), ctx.functional('k', 1, j)))
rhs = aadd(rhs, ascale(amul(F('Phi', 1), ctx.functional('k', 1, 2)), 4))
rhs = aadd(rhs, ascale(aadd(F('h', 2, (0, 0)), F('h', 2, (1, 1))), -0.5))
k2 = ctx.functional('k', 2, 2)
for m in (3, 5, 6):
    print('mask', m, 'max|k2_z - null| =', np.max(np.abs(k2[m] - rhs[m])), ' scale', np.max(np.abs(rhs[m])))
k1 = ctx.functional('k', 1, 2)
print('first order: max|k1_z + Phi1 + Psi1| =', max(np.max(np.abs(k1[m] + F('Phi', 1)[m] + F('Psi', 1)[m])) for m in (1, 2, 4)))
