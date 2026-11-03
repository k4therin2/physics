import sys; sys.path.insert(0, '../C1_einstein')
import sympy as sp
from kernels_api import load
K = load(); K3 = K['third']
for nm in ('Q','S00','S0z','SPhi'):
    e = K3[nm]
    print(nm, e.free_symbols, {str(a) for a in e.atoms(sp.core.function.AppliedUndef)})
print(K3['delta3_expr']); print(K3['V3_expr'])
ic = K['third_ic']
print('C3', ic['C3']); print('Q3_0', ic['Q3_0']); print('S00_0', ic['S00_0']); print('KPhi2_small', ic['KPhi2_small'])
