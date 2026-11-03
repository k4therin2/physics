"""W3: observer clock and worldline with finite initial time eta_i (manuscript Sec. obsevent).

Definitions (conformal frame, a = eta, coefficient convention X = eps X1 + eps^2 X2/2):
  initial event:   1 + dT(eta_init, x_o(eta_init)) = eta_init/eta_i       (T_phys = T_i, Tbar ∝ 1/eta)
  worldline:       x_o(eta) = x_lab + int_{eta_init}^{eta} v(eta', x_o(eta')) deta'   (label = position at the initial event)
  clock:           Delta tau = int_{eta_init}^{eta_o} eta / uhat^0 deta  = (etabar_o^2 - eta_i^2)/2   (held fixed)
  1/uhat^0 = 1 + eps Psi1 + eps^2 [Psi2/2 - Psi1^2/2 - v1^2/2] + O(eps^3)
One spatial dimension suffices for the algebra (dot products are restored in the text).
Fields on the worldline are expanded about x_lab:  F(eta, x_o(eta)) = F(eta,x_lab) + eps X1(eta) dF + ...
Outputs: deta_o1, deta_o2 (coefficient convention eta_o = etabar_o + eps deta_o1 + eps^2 deta_o2/2), dxo1, dxo2,
with the O(eta_i^2) pieces displayed, and two checks: (a) eta_i -> 0 reproduces the manuscript's eta_i-free formulas,
(b) the exact dilation family gives eta_o = etabar_o/lambda for finite eta_i.
"""
import sympy as sp

eps, e, ei, eo, s = sp.symbols('epsilon eta eta_i etabar_o s', positive=True)
Psi1, Psi2, dT1, dT2, v1, v2, gPsi1, gdT1, gv1 = [sp.Function(n) for n in
    ('Psi1', 'Psi2', 'dT1', 'dT2', 'v1', 'v2', 'gradPsi1', 'graddT1', 'gradv1')]
I = sp.Integral

# ---- initial event eta_init = eta_i + eps e1 + eps^2 e2/2 (worldline at the label point there: X(eta_init) = 0 exactly)
e1, e2 = sp.symbols('e1 e2')
etainit = ei + eps*e1 + eps**2*e2/2
# dT along the worldline at eta_init: x_o(eta_init) = x_lab exactly, so no gradient term
dT_init = eps*dT1(etainit) + eps**2*dT2(etainit)/2
eq = sp.series(1 + dT_init - etainit/ei, eps, 0, 3).removeO()
sol = sp.solve([eq.coeff(eps, 1), eq.coeff(eps, 2)], [e1, e2], dict=True)[0]
print('e1 =', sol[e1]); print('e2 =', sp.simplify(sol[e2]))

# ---- worldline displacement X(eta) = eps X1 + eps^2 X2/2, X1(eta) = int_{eta_init}^{eta} v1 -> int_{eta_i}^{eta} v1 + O(eps)
# X(eta) = int_{eta_init}^{eta} v(eta', x_lab+X) deta';  expand eta_init:  int_{eta_init}^{eta} = int_{eta_i}^{eta} - eps e1 v1(eta_i) + ...
X1 = lambda E: I(v1(s), (s, ei, E))
X2half = lambda E: I(v2(s)/2 + X1(s)*gv1(s), (s, ei, E)) - sol[e1]*v1(ei)   # eps^2 coefficient
print('X2/2 (label at initial event) =', X2half(e))

# ---- clock
d1, d2 = sp.symbols('d1 d2')
etao = eo + eps*d1 + eps**2*d2/2
f1 = lambda E: Psi1(E)
f2 = lambda E: Psi2(E)/2 - Psi1(E)**2/2 - v1(E)**2/2 + X1(E)*gPsi1(E)
# Delta tau = int_{eta_init}^{eta_o} eta (1 + eps f1 + eps^2 f2) d eta ;  split  int_{eta_init}^{eta_o} = int_{eta_i}^{etabar_o} + [etabar_o, eta_o] - [eta_i, eta_init]
def endpoint(E0, dE):     # int_{E0}^{E0+dE} eta(1+eps f1+eps^2 f2) to O(eps^2), dE = O(eps)
    g0 = E0*(1 + eps*f1(E0)); gp = 1
    return g0*dE + gp*dE**2/2
dtau = (eo**2 - ei**2)/2 + eps*I(s*f1(s), (s, ei, eo)) + eps**2*I(s*f2(s), (s, ei, eo)) \
       + endpoint(eo, etao - eo) - endpoint(ei, etainit - ei)
cl = sp.expand(sp.series(sp.expand(dtau - (eo**2 - ei**2)/2), eps, 0, 3).removeO())
c1 = cl.coeff(eps, 1); s1 = sp.solve(c1, d1)[0]
c2 = cl.coeff(eps, 2).subs(d1, s1); s2 = sp.solve(c2, d2)[0]
s1 = s1.subs(sol); s2 = sp.expand(s2.subs(sol))
print('\ndeta_o1 =', sp.simplify(s1))
print('deta_o2 =', s2)

# separate the eta_i-dependent (boundary) pieces
B1 = sp.simplify(s1 + I(s*Psi1(s), (s, ei, eo))/eo)
print('\nO(eta_i^2) boundary piece of deta_o1:', B1)
# eta_i->0 form of the manuscript: deta_o2^(0) = -(2/eo){ deta_o1^2/2 + eo deta_o1 Psi1(eo) + int_0^eo eta[Psi2/2 - Psi1^2/2 - v1^2/2 + X1 dPsi1] }
D1 = sp.Symbol('D1')
ms = -2/eo*(D1**2/2 + eo*D1*Psi1(eo) + I(s*f2(s), (s, ei, eo)))
diff = sp.expand(s2 - ms.subs(D1, s1))
print('deta_o2 minus manuscript form (lower limits eta_i, D1 = full deta_o1):', sp.simplify(diff))

# ---- check (b): dilation family, constant fields: lambda = 1 + eps a1 + eps^2 a2/2,  Psi: 1+2Psi = lambda^4,
#      1+dT = 1/lambda, v=0.  Expected: eta_o = etabar_o/lambda through O(eps^2) for any eta_i.
a1, a2 = sp.symbols('a1 a2')
lam = 1 + eps*a1 + eps**2*a2/2
PsiS = sp.series((lam**4 - 1)/2, eps, 0, 3).removeO(); dTS = sp.series(1/lam - 1, eps, 0, 3).removeO()
rep = {Psi1: sp.Lambda(s, PsiS.coeff(eps, 1)), Psi2: sp.Lambda(s, 2*PsiS.coeff(eps, 2)),
       dT1: sp.Lambda(s, dTS.coeff(eps, 1)), dT2: sp.Lambda(s, 2*dTS.coeff(eps, 2)),
       v1: sp.Lambda(s, 0), v2: sp.Lambda(s, 0), gPsi1: sp.Lambda(s, 0), gv1: sp.Lambda(s, 0)}
def R(x):
    for f, l in rep.items(): x = x.replace(f, l)
    return x.doit()
o1 = sp.simplify(R(s1)); o2 = sp.simplify(R(s2))
exp_ = sp.series(eo/lam - eo, eps, 0, 3).removeO()
print('\ndilation: deta_o1 - expected =', sp.simplify(o1 - exp_.coeff(eps, 1)),
      '; deta_o2 - expected =', sp.simplify(o2 - 2*exp_.coeff(eps, 2)))
print('(eta_i kept symbolic: both differences must vanish identically)')

# ---- LaTeX for the manuscript
print('\nLATEX deta_o1:', sp.latex(sp.simplify(s1)))
print('LATEX deta_o2:', sp.latex(s2))
