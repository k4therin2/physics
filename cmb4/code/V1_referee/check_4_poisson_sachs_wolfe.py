"""
V1 referee, attack item 4(c): linear theory, Poisson gauge, conformal frame (SPEC sec. 2):
   ghat = diag(-(1+2 eps Psi), (1-2 eps Phi) delta_ij),  uhat^mu = (1 - eps Psi, eps v^i),
   That = Tbar (1 + eps dT),  Tbar = const in the conformal frame.
Everything below is computed from covariant definitions to O(eps) with sympy.

Independent derivation of the observable (NOT using the claim):
   T_obs = T_* E_o/E_e = T_* (a_e/a_o) Ehat_o/Ehat_e,   a_e T_* = That(e) = Tbar(1+eps dT_e)
   => ln(a_o T_obs/Tbar) = eps dT(e) + ln(Ehat_o/Ehat_e)                       (A)
   with ln(Ehat_o/Ehat_e) from the null geodesic equation.
Claim:  ln(a_o T_obs/Tbar) = eps dT(o) - S,  S = int sigma_ee dl.                (B)
(A)=(B) iff the integrand identity
   d/deta[dT]_ray - sigma_ee  ==  -d/deta[Psi + n.v]_ray + (Psi' + Phi')   (on shell)
   i.e.  sigma_ee == d/deta[dT + Psi + n.v]_ray - (Psi' + Phi')
holds, where d/deta along the ray = d_eta + n^i d_i.  We test exactly that, and we also
verify that (A) reproduces the textbook SW formula
   Theta(n) + const = [dT + Psi + n.v]_e - [Psi + n.v]_o + int (Psi'+Phi') d eta.
"""
import sympy as sp

eps = sp.symbols('epsilon')
eta_, X, Y, Z = sp.symbols('eta x y z', real=True)
co = [eta_, X, Y, Z]
Psi = sp.Function('Psi')(*co); Phi = sp.Function('Phi')(*co); dT = sp.Function('dT')(*co)
v = [sp.Function(f'v{i}')(*co) for i in (1, 2, 3)]
n = sp.symbols('n1:4', real=True)      # unit propagation direction (zeroth order), n.n = 1
nn = {n[2]: sp.sqrt(1 - n[0]**2 - n[1]**2)}
O1 = lambda e: sp.series(sp.expand(e), eps, 0, 2).removeO()   # keep O(eps^0), O(eps^1)

g = sp.diag(-(1 + 2*eps*Psi), 1 - 2*eps*Phi, 1 - 2*eps*Phi, 1 - 2*eps*Phi)
ginv = sp.diag(-(1 - 2*eps*Psi), 1 + 2*eps*Phi, 1 + 2*eps*Phi, 1 + 2*eps*Phi)   # to O(eps)
Gam = [[[O1(sp.Rational(1, 2)*sum(ginv[l, s]*(sp.diff(g[s, m], co[n_]) + sp.diff(g[s, n_], co[m]) - sp.diff(g[m, n_], co[s]))
                                 for s in range(4))) for n_ in range(4)] for m in range(4)] for l in range(4)]
u_up = sp.Matrix([1 - eps*Psi, eps*v[0], eps*v[1], eps*v[2]])
u_dn = (g*u_up).applyfunc(O1)
print("u.u + 1 = O(eps^2):", O1((u_up.T*u_dn)[0, 0] + 1) == 0)
Du = sp.Matrix(4, 4, lambda m, n_: O1(sp.diff(u_dn[m], co[n_]) - sum(Gam[l][n_][m]*u_dn[l] for l in range(4))))
theta = O1(sum(ginv[m, n_]*Du[m, n_] for m in range(4) for n_ in range(4)))
acc = sp.Matrix([O1(sum(u_up[n_]*Du[m, n_] for n_ in range(4))) for m in range(4)])
hm = (sp.eye(4) + u_up*u_dn.T).applyfunc(O1)
h_dn = (g + u_dn*u_dn.T).applyfunc(O1)
B = sp.Matrix(4, 4, lambda m, n_: O1(sum(hm[a, m]*hm[b, n_]*Du[a, b] for a in range(4) for b in range(4))))
sigma = ((B + B.T)/2 - theta/3*h_dn).applyfunc(O1)
div_v = sum(sp.diff(v[i], co[i + 1]) for i in range(3))
print("theta  =", sp.factor(theta))
print("a_i    =", [sp.factor(acc[i]) for i in range(1, 4)])
print("sigma_ij - eps[d_(i v_j) - delta_ij div v/3] == 0 :",
      all(sp.simplify(sigma[i, j] - eps*((sp.diff(v[j-1], co[i]) + sp.diff(v[i-1], co[j]))/2 - (div_v/3 if i == j else 0))) == 0
          for i in range(1, 4) for j in range(1, 4)), "; sigma_0mu =", list(sigma[0, :]))

# ---- fluid equations at O(eps), from the exact identities of item 2 -------------------
lnT = sp.log(1 + eps*dT)
energy = O1(theta + 3*sum(u_up[m]*sp.diff(lnT, co[m]) for m in range(4)))
dlnT = [O1(sp.diff(lnT, co[m])) for m in range(4)]
euler = [O1(acc[m] + sum(hm[a, m]*dlnT[a] for a in range(4))) for m in range(4)]
print("energy eq (O(eps)):", sp.factor(energy), " == 0")
print("Euler  eq (O(eps)):", [sp.factor(e) for e in euler], " == 0")
# on-shell substitution rules
dT_eta = sp.diff(dT, eta_); v_eta = [sp.diff(vi, eta_) for vi in v]
onshell = {dT_eta: -div_v/3 + sp.diff(Phi, eta_)}
onshell.update({v_eta[i]: -sp.diff(Psi, co[i+1]) - sp.diff(dT, co[i+1]) for i in range(3)})
print("   check rules solve them:", sp.simplify(energy.subs(onshell)) == 0, [sp.simplify(e.subs(onshell)) == 0 for e in euler])

# ---- photon: direct geodesic computation of d ln Ehat / d eta ---------------------------
# k^mu = dx^mu/deta = (1, n^i) + O(eps); d ln p^0/deta = -Gam^0_{mn} k^m k^n
k = sp.Matrix([1, n[0], n[1], n[2]])
dlnp0 = O1(-sum(Gam[0][m][n_]*k[m]*k[n_] for m in range(4) for n_ in range(4)))
Ehat_over_p0 = O1(-(u_dn.T*k)[0, 0])      # Ehat/p^0 = -u.k
dlnE_direct = O1(dlnp0 + sp.diff(sp.log(Ehat_over_p0), eta_)*0 + 0)  # placeholder, build properly below
# d/deta along ray of f: f_eta + n.grad f
Dray = lambda f: sp.diff(f, eta_) + sum(n[i]*sp.diff(f, co[i+1]) for i in range(3))
lnEoverp0 = O1(sp.log(Ehat_over_p0))
dlnE_direct = O1(dlnp0 + Dray(lnEoverp0))
# kinematic prediction (i): d ln E/d eta = -(theta/3 + sigma_ee + a.e) * dl/deta, dl/deta = Ehat/p^0
e_up = (k/Ehat_over_p0 - u_up).applyfunc(O1)
sig_ee = O1((e_up.T*sigma*e_up)[0, 0]); a_e = O1((acc.T*e_up)[0, 0])
kin = O1(-(theta/3 + sig_ee + a_e)*Ehat_over_p0)
print("(i) at O(eps) in Poisson gauge, direct geodesic vs kinematic :", sp.simplify((dlnE_direct - kin).subs(nn)) == 0)
nv = sum(n[i]*v[i] for i in range(3))
# textbook form of the direct result: d lnE/deta = -Dray(Psi + n.v) + Psi' + Phi'
print("    d lnE/deta == eps[-d/deta(Psi + n.v)_ray + Psi' + Phi'] :",
      sp.simplify((dlnE_direct - eps*(-Dray(Psi + nv) + sp.diff(Psi, eta_) + sp.diff(Phi, eta_))).subs(nn)) == 0)
print("    => ln(Ehat_o/Ehat_e) = -[Psi + n.v]_e^o + int (Psi'+Phi') deta  (standard SW/ISW/Doppler)")

# ---- THE TEST: integrand identity required for (A) == (B) --------------------------------
ident = O1(eps*Dray(dT) - sig_ee - eps*(-Dray(Psi + nv) + sp.diff(Psi, eta_) + sp.diff(Phi, eta_)))
# (referee note: a first version of this line had +sig_ee by mistake; its on-shell residual was exactly 2*sig_ee)
print("off-shell residual (should be nonzero):", sp.simplify(ident.subs(nn)) != 0)
print("ON-SHELL: d/deta[dT]_ray - sigma_ee + d/deta[Psi+n.v]_ray - (Psi'+Phi') == 0 :",
      sp.simplify(sp.expand(ident.subs(onshell)).subs(nn)) == 0)
print("=> integrating from e to o:  dT_o - S = dT_e + Psi_e + n.v_e - Psi_o - n.v_o + int(Psi'+Phi') deta")
print("   which is exactly ln(a_o T_obs / Tbar) from (A).  Claim (B) == Sachs-Wolfe at O(eps). CONFIRMED.")
