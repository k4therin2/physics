"""debug_third.py -- V3: isolate the third-order discrepancy.
(1) validate my pair-field construction against MY second-order constraints (00, 0i, traceless-scalar) with abstract KPhi2;
(2) run C1's own contracted_kernel (exec of third_order.py lines 1-95, read-only) on single monomials and compare with
    my discrete three-mode realisation of the same monomials."""
import sys, os, random, sympy as sp
from pw_tools import *
from third_tools import *
C1 = '/home/user/physics/cmb4/code/C1_einstein'
E = E_mine()
z1, z2, z3 = sp.symbols('zeta1 zeta2 zeta3')
pt = rand_point(0); rng = random.Random(7)
# ---------------- (1) pair fields satisfy my second-order constraints (coefficient of zeta1 zeta3, pair (q,k))
f = three_mode_fields(pt, z1, z2, z3, amp=1)
def c2(a, b, zz):
    return coeff_at_origin(substitute(E[('E', a, b, 2)], f), zz)
for zz, lab in (((z1, z3), 'pair (q,k)'), ((z2, z3), 'pair (-q,k)')):
    r00 = ev(c2(0, 0, zz), rng, pt['eta']); r0 = [ev(c2(0, i, zz), rng, pt['eta']) for i in (1, 2, 3)]
    rTL = ev(c2(3, 3, zz) - (c2(1, 1, zz) + c2(2, 2, zz)) / 2, rng, pt['eta'])
    rxz = ev(c2(1, 3, zz), rng, pt['eta']); rxy = ev(c2(1, 2, zz), rng, pt['eta'])
    print(lab, "constraints 00:", sp.N(r00, 6), " 0i:", [sp.N(r_, 6) for r_ in r0], " TL:", sp.N(rTL, 6), " (xz, xy involve h2'' -> need not vanish:", sp.N(rxz, 6), sp.N(rxy, 6), ")")
# ---------------- (2) C1's contracted_kernel on monomials vs discrete realisation
src = open(os.path.join(C1, 'third_order.py')).read()
cut = src.index('\nN = 3\n')
ns = {}
cwd = os.getcwd(); os.chdir(C1)
try:
    exec(src[:cut], ns)
finally:
    os.chdir(cwd)
jet = ns['jet']
fX = three_mode_fields(pt, z1, z2, z3, amp=1)      # X2 itself (discrete amplitude 2K)
def realise(mono_jets):
    tot = 1
    for (fld, od, m) in mono_jets:
        e = fX.get('%s%d' % (fld, od), sp.Integer(0))
        for var, n in zip((t, x, y, z), m):
            if n: e = sp.diff(e, var, n)
        tot = tot * e
    return coeff_at_origin(sp.expand(tot), (z1, z2, z3))
tests = {
    'Phi1*Phi2': [('Phi', 1, (0, 0, 0, 0)), ('Phi', 2, (0, 0, 0, 0))],
    'd_z Phi1 * d_t delta2': [('Phi', 1, (0, 0, 0, 1)), ('delta', 2, (1, 0, 0, 0))],
    'v1_x * om2_x': [('v1', 1, (0, 0, 0, 0)), ('om1', 2, (0, 0, 0, 0))],
    'v1_z * v2_z': [('v3', 1, (0, 0, 0, 0)), ('v3', 2, (0, 0, 0, 0))],
    'Phi1 * h2_11': [('Phi', 1, (0, 0, 0, 0)), ('h11', 2, (0, 0, 0, 0))],
    'delta1 * d_z v2_x': [('delta', 1, (0, 0, 0, 0)), ('v1', 2, (0, 0, 0, 1))],
    'Phi1 * Psi2': [('Phi', 1, (0, 0, 0, 0)), ('Psi', 2, (0, 0, 0, 0))],
    'Phi1^3': [('Phi', 1, (0, 0, 0, 0))] * 3,
    'd_x v1_x * delta1 * d_t Psi1': [('v1', 1, (0, 1, 0, 0)), ('delta', 1, (0, 0, 0, 0)), ('Psi', 1, (1, 0, 0, 0))],
    'v1_z * v1_z * d_t Phi1': [('v3', 1, (0, 0, 0, 0)), ('v3', 1, (0, 0, 0, 0)), ('Phi', 1, (1, 0, 0, 0))],
    'Phi1 * v2_x': [('Phi', 1, (0, 0, 0, 0)), ('v1', 2, (0, 0, 0, 0))],
    'Phi1 * h2_13': [('Phi', 1, (0, 0, 0, 0)), ('h13', 2, (0, 0, 0, 0))],
}
for name, mj in tests.items():
    P = ns['const'](1)
    for (fld, od, m) in mj:
        P = P * jet(fld, od, m)
    kc1 = ns['on_shell'](ns['contracted_kernel'](P))
    mine = realise(mj)
    vm = ev(mine, rng, pt['eta'])
    vp = ev(c1_image_sum(kc1, pt, 1), rng, pt['eta']); vn = ev(c1_image_sum(kc1, pt, -1), rng, pt['eta'])
    print("%-32s mine=%s  C1(+s)=%s  C1(-s)=%s  -> match: %s" % (name, sp.N(vm, 10), sp.N(vp, 10), sp.N(vn, 10),
          '+' if abs(vm - vp) < 1e-20 else ('-' if abs(vm - vn) < 1e-20 else 'NONE')))
