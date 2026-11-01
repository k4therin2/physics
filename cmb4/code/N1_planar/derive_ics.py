"""Gradient expansion (regular solution) for the planar comoving-gauge system:
   Ah = e^{zeta}(1 + sum_n a_n(x) t^n), Bh = e^{zeta}(1 + sum_n b_n(x) t^n), NONLINEAR in zeta.
Solved order by order from the evolution equations; constraints verified order by order.
Output: ics_generated.pkl with a_n, b_n (n=2..NMAX) as expressions in zeta and its x-derivatives.
"""
import sympy as sp, pickle, sys, time
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 6
d = pickle.load(open('eqs_generated.pkl', 'rb'))
S = d['symbols']
t, Ah, Bh, alh, Ax, Bx, Axx, Bxx, pA, pB, pAx, pBx, ze, zx, zxx = [S[k] for k in
    ('t', 'Ah', 'Bh', 'alh', 'Ax', 'Bx', 'Axx', 'Bxx', 'pA', 'pB', 'pAx', 'pBx', 'zeta', 'zx', 'zxx')]
x = sp.Symbol('x', real=True)
zeta = sp.Function('zeta')(x)
an = [sp.Function(f'a{n}')(x) for n in range(NMAX + 1)]
bn = [sp.Function(f'b{n}')(x) for n in range(NMAX + 1)]
u = sum(an[n]*t**n for n in range(1, NMAX + 1))
v = sum(bn[n]*t**n for n in range(1, NMAX + 1))
Ah_e = sp.exp(zeta)*(1 + u)
Bh_e = sp.exp(zeta)*(1 + v)
alh_e = ((1 + u)*(1 + v)**2)**sp.Rational(1, 3)
pA_e = sp.diff(sp.log(Ah_e), t)/alh_e
pB_e = sp.diff(sp.log(Bh_e), t)/alh_e
sub = {Ah: Ah_e, Bh: Bh_e, alh: alh_e, Ax: sp.diff(Ah_e, x), Bx: sp.diff(Bh_e, x),
       Axx: sp.diff(Ah_e, x, 2), Bxx: sp.diff(Bh_e, x, 2), pA: pA_e, pB: pB_e,
       pAx: sp.diff(pA_e, x), pBx: sp.diff(pB_e, x), ze: zeta, zx: sp.diff(zeta, x), zxx: sp.diff(zeta, x, 2)}
EA = (sp.diff(pA_e, t) - d['pA_t'].subs(sub))*t**2     # order t^n equations
EB = (sp.diff(pB_e, t) - d['pB_t'].subs(sub))*t**2
Ham = d['Ham'].subs(sub)*t**4
Mom = d['Mom'].subs(sub)*t
def ser(e, n):
    return sp.series(e, t, 0, n + 1).removeO()
sol = {}
t0 = time.time()
EAs = sp.expand(ser(EA, NMAX))
EBs = sp.expand(ser(EB, NMAX))
print("series done", time.time() - t0, flush=True)
for n in range(1, NMAX + 1):
    cA = sp.simplify(EAs.coeff(t, n).subs(sol).doit())
    cB = sp.simplify(EBs.coeff(t, n).subs(sol).doit())
    s = sp.solve([cA, cB], [an[n], bn[n]], dict=True)
    if n == 1:
        print("order t^1 solve:", s)   # expect a1=b1 undetermined by evolution -> fixed =0 by constraints
        sol[an[1]] = 0; sol[bn[1]] = 0
        continue
    assert len(s) == 1, (n, s)
    sol[an[n]] = sp.simplify(s[0][an[n]])
    sol[bn[n]] = sp.simplify(s[0][bn[n]])
    print(f"n={n}: a{n} =", sol[an[n]]); print(f"      b{n} =", sol[bn[n]], flush=True)
# constraints check
Hs = sp.expand(ser(Ham, NMAX)).subs(sol).doit()
Ms = sp.expand(ser(Mom, NMAX)).subs(sol).doit()
for n in range(0, NMAX + 1):
    print(f"Ham coeff t^{n}:", sp.simplify(Hs.coeff(t, n)), "   Mom coeff t^{n}:", sp.simplify(Ms.coeff(t, n)))
# linear check: a2+2b2 should be zeta''/3 at linear order
lin = sp.simplify((sol[an[2]] + 2*sol[bn[2]]))
print("a2+2b2 =", lin)
pickle.dump({'sol': {str(k): v for k, v in sol.items()}, 'NMAX': NMAX}, open(f'ics_generated_{NMAX}.pkl', 'wb'))
print("saved", time.time() - t0)
