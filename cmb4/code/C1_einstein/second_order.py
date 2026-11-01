"""second_order.py -- L[X2] = S_X[fields1]: SVT decomposition at second order, Fourier kernels, Green's functions.
Equations are E^{(2)} := 2 * [eps^2] E  so that X2 enters linearly with unit weight and the sources are
2 x (quadratic in first-order fields).  Kernels K(q,p) are per zeta_q zeta_p with p = k - q,
X2_k = int d^3q/(2pi)^3 K_X(q, k-q; eta) zeta_q zeta_{k-q}."""
import sympy as sp, pickle, time
from tpoly import *
from fourier import *

t0 = time.time()
d = load_polys('E.pkl'); E = d['E']; cons = d['cons']
N = 2
def comp(a, b):
    P = E[a][b].order_part(N) * sp.factorial(N)
    return linear_part(P, N), quadratic_kernel(P)
lin, ker = {}, {}
for (a, b) in [(0, 0), (0, 1), (0, 2), (0, 3), (1, 1), (2, 2), (3, 3), (1, 2), (1, 3), (2, 3)]:
    lin[(a, b)], ker[(a, b)] = comp(a, b)
print("components done (%.1fs)" % (time.time() - t0))

Phi2, Psi2, de2 = amp('Phi', 2), amp('Psi', 2), amp('delta', 2)
v2x, v2y, v2z = [amp('v%d' % i, 2) for i in (1, 2, 3)]
omx, omy = amp('om1', 2), amp('om2', 2)
hp, hx = amp('h11', 2), amp('h12', 2)
Q2 = sp.Function('Q2')(eta)        # Q2 := Phi2 - Psi2
qperp2 = q**2 * (1 - mu**2)

# first-order transfer functions obey: T_Psi = T_Phi, T_delta, T_V from constraints, T_Phi'' from ODE
TPhi = lambda kk: TF['Phi'](kk, eta)
first_rules = {}
def on_shell(expr):
    """impose first-order relations on an expression containing T_X(kk, eta)."""
    expr = expr.doit()
    for kk in (q, p, k):
        TP = TPhi(kk)
        expr = expr.subs(TF['Psi'](kk, eta), TP)
        expr = expr.subs(TF['delta'](kk, eta), -sp.Rational(2, 3) * eta**2 * kk**2 * TP - 2 * eta * sp.Derivative(TP, eta) - 2 * TP)
        expr = expr.subs(TF['V'](kk, eta), -(eta**2 / 2) * (sp.Derivative(TP, eta) + TP / eta))
        expr = expr.doit()
        # eliminate second and higher eta-derivatives via ODE
        for nn in (4, 3, 2):
            d2 = -4 / eta * sp.Derivative(TP, eta) - kk**2 / 3 * TP
            dn = sp.diff(d2, eta, nn - 2)
            expr = expr.subs(sp.Derivative(TP, (eta, nn)), dn).doit()
        expr = expr.subs(sp.Derivative(TP, (eta, 2)), -4 / eta * sp.Derivative(TP, eta) - kk**2 / 3 * TPhi(kk)).doit()
    return sp.expand(expr)

def invariants(expr):
    """express via magnitudes k,q,p (mu = (k^2+q^2-p^2)/(2kq))."""
    expr = sp.expand(expr)
    expr = expr.subs(mu, (k**2 + q**2 - p**2) / (2 * k * q))
    return sp.simplify(expr)

# ---------------- traceless scalar:  (zz) - (xx+yy)/2
linTL = sp.expand(lin[(3, 3)] - (lin[(1, 1)] + lin[(2, 2)]) / 2)
kerTL = sp.expand(ker[(3, 3)] - (ker[(1, 1)] + ker[(2, 2)]) / 2)
print("TL linear part:", sp.factor(linTL))
assert sp.simplify(linTL + k**2 * (Phi2 - Psi2)) == 0
K_Q = on_shell(kerTL / k**2)            # Q2 = Phi2 - Psi2 = int K_Q zeta zeta
print("K_Q (raw, mu form):", K_Q)
# ---------------- 00 -> delta2
sol_de2 = sp.solve(lin[(0, 0)] + sp.Symbol('S00'), de2)[0]
print("delta2 =", sol_de2, "   with S00 = source kernel of 2*[eps^2]E_00")
K_S00 = on_shell(ker[(0, 0)])
# ---------------- 0i
print("lin 0z:", lin[(0, 3)]); print("lin 0x:", lin[(0, 1)]); print("lin 0y:", lin[(0, 2)])
sol_v2z = sp.solve(lin[(0, 3)] + sp.Symbol('S0z'), v2z)[0]
sol_v2x = sp.solve(lin[(0, 1)] + sp.Symbol('S0x'), v2x)[0]
K_S0z = on_shell(ker[(0, 3)]); K_S0x = on_shell(ker[(0, 1)]); K_S0y = on_shell(ker[(0, 2)])
assert K_S0y == 0
print("S0y = 0 (as required by reflection symmetry)")
# ---------------- vector ij: xz
print("lin xz:", lin[(1, 3)], " lin yz:", lin[(2, 3)])
K_Sxz = on_shell(ker[(1, 3)]); K_Syz = on_shell(ker[(2, 3)]); assert K_Syz == 0
# ---------------- tensor: xx - yy, xy
linT = sp.expand(lin[(1, 1)] - lin[(2, 2)]); kerT = on_shell(ker[(1, 1)] - ker[(2, 2)])
print("lin xx-yy:", linT); print("lin xy:", lin[(1, 2)])
K_Sxy = on_shell(ker[(1, 2)]); assert K_Sxy == 0
# ---------------- trace -> Phi2 ODE:  (trace - 00)/6 with Psi2 = Phi2 - Q2, delta2 eliminated
lin_tr = sp.expand(lin[(1, 1)] + lin[(2, 2)] + lin[(3, 3)])
ker_tr = on_shell(ker[(1, 1)] + ker[(2, 2)] + ker[(3, 3)])
odeL = sp.expand((lin_tr - lin[(0, 0)]) / 6).subs(Psi2, Phi2 - Q2).doit()
odeL = sp.expand(odeL)
print("Phi2 ODE linear part (with Q2):", odeL)
# check structure
target = sp.Derivative(Phi2, (eta, 2)) + 4 / eta * sp.Derivative(Phi2, eta) + k**2 / 3 * Phi2 - sp.Derivative(Q2, eta) / eta + k**2 / 3 * Q2
assert sp.simplify(odeL - target) == 0
print("VERIFIED: Phi2'' + 4Phi2'/eta + k^2 Phi2/3 = S_Phi,  S_Phi = -(S_tr - S_00)/6 + Q2'/eta - k^2 Q2/3")
K_SPhi = sp.expand(-(ker_tr - K_S00) / 6 + sp.diff(K_Q, eta) / eta - k**2 / 3 * K_Q)
K_SPhi = on_shell(K_SPhi)

# ---------------- covariant vector / tensor kernels
# vector: transverse source S_i(q,p) = P_ij(k) q_j * Kvec(q,p);   in frame x-component = qperp * Kvec
sq = sp.sqrt(1 - mu**2)
K_vec_0i = sp.expand(sp.cancel(K_S0x / (q * sq)))     # kernel of the transverse part of 2[eps^2]E_0i source
K_vec_ij = sp.expand(sp.cancel(K_Sxz / (q * sq)))     # kernel of the vector part of 2[eps^2]E_ij source (xz comp = ik_z (q_x stuff))
# tensor: (S_xx - S_yy) = qperp^2 * Kten ;  h2_ij = Lambda_ij,kl q_k q_l * (...)
K_ten = sp.expand(sp.cancel(kerT / qperp2))
for nm, ex in (('K_vec_0i', K_vec_0i), ('K_vec_ij', K_vec_ij), ('K_ten', K_ten)):
    assert not ex.has(sq), (nm, ex)
print("vector/tensor kernels free of sqrt(1-mu^2): OK")

# ---------------- omega2 ODE from xz;  v2_x from 0x
# lin xz = c1 * (om' + 2 om/eta) * (i k)/2 ... solve generally:
omp = sp.Derivative(omx, eta)
sol_omp = sp.solve(lin[(1, 3)] + sp.Symbol('Sxz'), omp)[0]
print("om2_x' =", sol_omp)

# ---------------- h2 ODE
sol_hpp = sp.solve(linT + sp.Symbol('Sxxyy'), sp.Derivative(hp, (eta, 2)))[0]
print("h_p'' =", sol_hpp)

res = dict(linTL=linTL, K_Q=K_Q, lin00=lin[(0, 0)], K_S00=K_S00, sol_de2=sol_de2,
           lin0z=lin[(0, 3)], lin0x=lin[(0, 1)], K_S0z=K_S0z, K_S0x=K_S0x, sol_v2z=sol_v2z, sol_v2x=sol_v2x,
           linxz=lin[(1, 3)], K_Sxz=K_Sxz, sol_omp=sol_omp, linT=linT, kerT=kerT, sol_hpp=sol_hpp,
           lin_tr=lin_tr, ker_tr=ker_tr, odeL=odeL, K_SPhi=K_SPhi, K_vec_0i=K_vec_0i, K_vec_ij=K_vec_ij, K_ten=K_ten)
pickle.dump({kk: sp.srepr(vv) for kk, vv in res.items()}, open('second_order_raw.pkl', 'wb'))
print("saved second_order_raw.pkl (%.1fs)" % (time.time() - t0))
