# Lead's independent check of the linear RD solution in Poisson gauge (fluid variables), SPEC conventions.
# Equations (standard, Poisson gauge, radiation, a=eta, H=1/eta, Psi=Phi):
#   Phi'' + 3H(1+c_s^2)Phi' + c_s^2 k^2 Phi = 0  with c_s^2 = 1/3  ->  Phi'' + (4/eta)Phi' + (k^2/3)Phi = 0
#   Poisson: k^2 Phi + 3H(Phi' + H Phi) = -(3/2) H^2 delta   (8piG a^2 rhobar = 3H^2)   [00 eq]
#   0i:  Phi' + H Phi = -(3/2) H^2 (1+w) * (theta_v/k^2)  ... we use  v_i = -i k_i V,  (Phi'+H Phi) = -(3/2)H^2 (4/3) V ... check signs by substitution below
import sympy as sp
eta,k=sp.symbols('eta k',positive=True)
x=k*eta/sp.sqrt(3)
TPhi=-2*(sp.sin(x)-x*sp.cos(x))/x**3
ode=sp.diff(TPhi,eta,2)+4/eta*sp.diff(TPhi,eta)+k**2/3*TPhi
print("ODE residual:",sp.simplify(ode))
print("TPhi(x->0):",sp.limit(TPhi,eta,0))
H=1/eta
# delta from 00 equation: delta = -(2/(3H^2)) [k^2 Phi + 3H(Phi' + H Phi)]
delta=sp.simplify(-(2/(3*H**2))*(k**2*TPhi+3*H*(sp.diff(TPhi,eta)+H*TPhi)))
dT=delta/4
print("delta_T(x->0)/Phi(0):",sp.limit(dT,eta,0)/sp.limit(TPhi,eta,0))
print("zeta = delta_T - Phi at eta->0:",sp.limit(dT-TPhi,eta,0))
print("delta_T series:",sp.series(dT,eta,0,4))
# velocity potential from 0i: (Phi' + H Phi) = -(3/2)H^2 (rho+P)/rho * (something). Standard: Phi'+H Psi = -4 pi G a^2 (rho+P) v  (v_i = d_i v)  => with 8piG a^2 rho = 3H^2, (rho+P)=4rho/3: Phi'+H Phi = -2 H^2 v
v=sp.simplify(-(sp.diff(TPhi,eta)+H*TPhi)/(2*H**2))
print("v (velocity potential, v_i = d_i v) series:",sp.series(v,eta,0,5))
# continuity check: delta' = -(4/3) nabla.v + 4 Phi'  -> in k-space delta' = (4/3) k^2 v + 4 Phi'
print("continuity residual:",sp.simplify(sp.diff(delta,eta)-(sp.Rational(4,3)*k**2*v+4*sp.diff(TPhi,eta))))
# Euler: v' = -(1/4) delta - Phi  (radiation, no Hubble friction)
print("Euler residual:",sp.simplify(sp.diff(v,eta)+delta/4+TPhi))
print("delta_T closed:",sp.simplify(dT))
print("v closed:",sp.simplify(v))
