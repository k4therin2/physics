"""structure_phi2.py -- collect P = sum_w [a_w(eta) Cin(w eta) + b_w(eta) Si(w eta)] + R(eta); check reality; factor coefficients."""
import sympy as sp, pickle, sys, time
sys.path.insert(0, '.'); sys.path.insert(0, '../C1_einstein')
from trig_integrate import *
from kernels_api import load, eta, k, q, p
d = pickle.load(open('phi2_closed.pkl', 'rb'))
L = {'Cin': Cin}
A = sp.parse_expr(d['A'], local_dict=L); B = sp.parse_expr(d['B'], local_dict=L)
t0 = time.time()
def split(e):
    e = sp.expand(e)
    re_, im_ = e.as_real_imag()
    im_ = sp.expand(im_)
    return sp.expand(re_), im_
def collect_special(e):
    """-> (dict {Cin(arg): coeff}, dict {Si(arg): coeff}, remainder)"""
    e = sp.expand(e)
    cins = {f: 0 for f in e.atoms(Cin)}; sis = {f: 0 for f in e.atoms(Si)}; rem = 0
    for term in sp.Add.make_args(e):
        fs = [f for f in term.atoms(Cin) | term.atoms(Si)]
        if not fs:
            rem += term
        else:
            assert len(fs) == 1
            f = fs[0]
            (cins if isinstance(f, Cin) else sis)[f] += term / f
    return cins, sis, rem
out = {}
for nm, e in (('A', A), ('B', B)):
    re_, im_ = split(e)
    print(nm, "imag part (should be 0):", sp.simplify(im_) if im_ != 0 else 0)
    cins, sis, rem = collect_special(re_)
    print(nm, ": %d Cin, %d Si terms; remainder nops %d  (%.1fs)" % (len(cins), len(sis), sp.count_ops(rem), time.time() - t0))
    cf = {}
    for f, c in list(cins.items()) + list(sis.items()):
        c = sp.factor(sp.cancel(sp.expand(c)))
        cf[f] = c
        print("   ", f, ":", c)
    # the remainder: elementary; collect by trig functions of eta
    rem = sp.collect(sp.expand(rem), rem.atoms(sp.sin, sp.cos), sp.factor)
    out[nm] = (cf, rem)
pickle.dump({nm: ({sp.srepr(f): sp.srepr(c) for f, c in cf.items()}, sp.srepr(rem)) for nm, (cf, rem) in out.items()}, open('phi2_structure.pkl', 'wb'))
print("done %.1fs" % (time.time() - t0))
