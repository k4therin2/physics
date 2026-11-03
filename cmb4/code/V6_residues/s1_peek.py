import sys; sys.path.insert(0,'../C1_einstein')
from kernels_api import load
from fourier import eta,k,q,p,mu,TF
import sympy as sp
K=load()
for nm,v in K['second'].items():
    print(nm, type(v), (sp.count_ops(v) if isinstance(v,sp.Basic) else ''), (v.free_symbols if isinstance(v,sp.Basic) else ''))
for nm,v in K['third'].items():
    print('3:',nm, sp.count_ops(v) if isinstance(v,sp.Basic) else type(v))
print(K['third_ic'].keys(), K['third_ic']['tads'].keys())
print(K['second']['K_Q'])
print(K['second']['K_S0z'])
print(K['second']['K_Phi2_0'])
