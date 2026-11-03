import sympy as sp, pickle
x, y, s = sp.symbols('x y s', positive=True)
Tc = lambda z: -2 * (sp.sin(z) - z * sp.cos(z)) / z**3
# A_Q = s^2 [2T^2 + (d(xT)/dx)^2];  eta = sqrt3 x/s;  eta^2 tau = -int_0^eta e'^2 A_Q de' = -(3 sqrt3/s^3) s^2 int_0^x y^2 [..] dy
f = sp.expand(sp.expand_trig(sp.simplify(y**2 * (2 * Tc(y)**2 + sp.diff(y * Tc(y), y)**2))))
print('integrand y^2[2T^2+(yT)\'^2] =', sp.simplify(f))
F = sp.integrate(f, y)
F = sp.simplify(F)
F0 = sp.limit(F, y, 0)
J = sp.simplify(F.subs(y, x) - F0)
print('J(x) = int_0^x =', J)
print('check dJ/dx - f(x):', sp.simplify(sp.diff(J, x) - f.subs(y, x)))
# tau = -(3 sqrt3 / s) J / eta^2 = -(3 sqrt3/s) J s^2/(3 x^2) = -(sqrt3 s/x^2)... careful: eta^2 = 3x^2/s^2
tau = sp.simplify(-(3 * sp.sqrt(3) / s**3) * s**2 * J / (3 * x**2 / s**2))
print('tau(x) =', tau)
print('tau series:', sp.series(tau, x, 0, 6))
pickle.dump({'tau_x': sp.srepr(tau), 'J': sp.srepr(J)}, open('s1_tau.pkl', 'wb'))
