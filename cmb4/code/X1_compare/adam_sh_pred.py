"""Adam's superhorizon expansions (his Sec. 9.1-9.2, eqs. after (83)) specialised to our evaluation line
eta_* = 2, eta_o = 8, initial event t_I -> 0 (our evaluators use eta_i -> 0), as ratios to C11 in powers of x = k0 eta_o."""
import sympy as sp
k, x, i = sp.symbols('k x i')
o, e = sp.Integer(8), sp.Integer(2)
D = o**2 - e**2
C11 = 4*sp.pi/2025*k**4*D**2 - 4*sp.pi/637875*k**6*(o-e)**2*(o+e)*(59*o**3+59*o**2*e-106*o*e**2+114*e**3)
Q22 = 68*sp.pi/6075*k**4*D**2 - 8*sp.pi/637875*k**6*D*(105*i**4+719*o**4-1920*o**2*e**2+2560*o*e**3-1464*e**4)
Q13 = 16*sp.pi/2025*k**4*D**2 - 2*sp.pi/1913625*k**6*D*(2835*i**4+9973*o**4-25440*o**2*e**2+33920*o*e**3-21288*e**4)
for nm, Q in (('1/4 C22', Q22), ('1/3 ReC13', Q13), ('sum', Q22 + Q13)):
    r = sp.series((Q / C11).subs(i, 0).subs(k, x / 8), x, 0, 3).removeO()
    print('%-10s / C11 (i=0) = %s  = %s' % (nm, sp.nsimplify(r), sp.N(r, 8)))
print('C11 = ', sp.N(sp.expand(C11.subs(k, x/8)), 10))
