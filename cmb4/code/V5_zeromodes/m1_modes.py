"""V5 step 4 (referee M1): homogeneous regular modes of the second-order tensor and vector equations (C1 linear operators, RD, a = eta).
Tensor:  h'' + (2/eta) h' + k^2 h = 0.   Vector (C1): ij-vector  omega' + (2/eta) omega = 0 (k != 0),
0i-transverse constraint  -(k^2/2) omega - (4/eta^2)(omega + v^T) = 0,  transverse Euler (radiation, conformal frame):
[(rhohat + Phat)(v^T + omega)]' = 0 with rhohat + Phat = const  =>  v^T + omega = C_v (conserved vorticity, relativistic Kelvin)."""
import sympy as sp
eta, k, Cv = sp.symbols('eta k C_v', positive=True)
h = sp.Function('h')
sol = sp.dsolve(h(eta).diff(eta, 2) + 2*h(eta).diff(eta)/eta + k**2*h(eta), h(eta))
print('tensor general solution:', sol)
j0 = sp.sin(k*eta)/(k*eta); y0 = sp.cos(k*eta)/(k*eta)
L = lambda f: sp.simplify(f.diff(eta, 2) + 2*f.diff(eta)/eta + k**2*f)
print('L[j0] =', L(j0), '; L[cos/(k eta)] =', L(y0))
print('j0(k eta) series:', sp.series(j0, eta, 0, 7), ' -> analytic (even power series), j0(0) = 1: REGULAR, not excluded by regularity')
print('cos(k eta)/(k eta) series:', sp.series(y0, eta, 0, 2), ' -> singular')
print('Frobenius indicial roots at eta=0: r(r-1) + 2r = 0 ->', sp.solve(sp.Symbol('r')*(sp.Symbol('r') - 1) + 2*sp.Symbol('r'), sp.Symbol('r')), '(r=0 regular, r=-1 singular)')
print('k = 0: h = A + B/eta; A = constant anisotropic shear = pure gauge (V5 run_all: shear generator => rotation of the sky).')
print('=> h_n(eta->0) = 0 must be DECLARED (no primordial tensors at any order); it also removes the k=0 constant-shear zero mode.')
om = sp.Function('omega')
solv = sp.dsolve(om(eta).diff(eta) + 2*om(eta)/eta, om(eta)); print('\nvector ij equation homogeneous solution:', solv, ' (singular unless C1 = 0)')
w = sp.Symbol('w'); vT = Cv - w
wsol = sp.solve(-(k**2/2)*w - 4/eta**2*(w + vT), w)[0]
print('constraint with Kelvin (v^T + omega = C_v): omega =', wsol, '; satisfies omega\' + 2 omega/eta = 0:', sp.simplify(sp.diff(wsol, eta) + 2*wsol/eta) == 0)
print('=> for k != 0 the only homogeneous vector mode is omega = -8 C_v/(k^2 eta^2), v^T = C_v - omega: singular, carried by primordial vorticity C_v;')
print('   regularity <=> C_v = 0 <=> omega = v^T = 0.  No regular homogeneous vector mode exists for k != 0 (C1/V3 statement confirmed).')
print('k = 0: the k^2 term drops: constraint gives only omega + v^T = 0 (C_v = 0), omega(eta) free = time-dependent translation c(eta):')
print('   omega = c\', v = -c\' (V5 gauge engine, run_all: translation generator) -> pure gauge; also linear-in-x vector modes = rotations w(eta) (pure gauge).')
