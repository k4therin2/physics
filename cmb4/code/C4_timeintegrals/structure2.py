"""structure2.py -- special-function skeleton of the closed forms: coefficients of Cin(w eta), Si(w eta) in A, B (Phi2) and in
A_h, B_h (h2); test the sign-flip ("master coefficient") structure  c_{(sk k + sq q + sp p)} = c_master(sk k, sq q, sp p)."""
import sympy as sp, pickle, sys, time
sys.path.insert(0, '.'); sys.path.insert(0, '../C1_einstein')
from trig_integrate import *
from kernels_api import eta, k, q, p
t0 = time.time()
L = {'Cin': Cin}
d = pickle.load(open('phi2_closed.pkl', 'rb'))
print(list(d.keys()))
s3 = sp.sqrt(3)
def skeleton(e):
    e = sp.expand(e)
    out = {}
    for term in sp.Add.make_args(e):
        fs = list(term.atoms(Cin) | term.atoms(Si))
        if not fs: continue
        assert len(fs) == 1
        f = fs[0]; w = sp.expand(f.args[0] / eta)
        out.setdefault((f.func.__name__, w), 0)
        out[(f.func.__name__, w)] += term / f
    return {kk: sp.factor(sp.expand(v)) for kk, v in out.items()}
res = {}
for nm in ('A', 'B'):
    if nm not in d: continue
    e = sp.parse_expr(d[nm], local_dict=L)
    sk = skeleton(e); res[nm] = sk
    print(nm, "(%.0fs)" % (time.time() - t0))
    for (fn, w), c in sorted(sk.items(), key=str):
        print("   %s(%s eta): %s" % (fn, w, c))
pickle.dump({nm: {str(kk): sp.srepr(v) for kk, v in sk.items()} for nm, sk in res.items()}, open('phi2_skeleton.pkl', 'wb'))
