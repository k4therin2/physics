# Report A3 — third-order contracted map kernel K_3c and C_l^(13)
Agent A3 (relaunch 2026-11-03). Code: `cmb4/code/A3_third_order/`. Conventions: SPEC + C2 (ray xbar = n chi, e = -n).

## Milestone 1 — engine (DONE, validated at orders 1-2)
* `tf.py`: transfer functions (series branch for small x), C1 second-order kernels, C1 contracted third-order kernels
  K3c_{Q,S00,S0z,SPhi} lambdified with KPhi2, Kh, tadpoles as explicit arguments, C1 tadpoles, C3 amplitude.
* `grid.py`: global Chebyshev grid on (0, eta_o]; Green's integrals (G_Phi with small-K eta Taylor branch); second-order
  pair profiles (AB-coefficients) and the K=0 tadpole profiles.  CHECK: angle-averaged K->0 limit of the pair profiles
  converges to C1's tadpoles (K=0.005: Phi2, Psi2, dT2 agree to 1e-4 relative); tadpole delta2 constraint 1e-11.
* `engine.py`: polarisation engine.  zeta = a e^{ik0.x} + b e^{ik1.x} + c e^{ik2.x} with a^2=b^2=c^2=0; every quantity
  is an element of R[a,b,c]/(a^2,b^2,c^2); C2's term lists (shear and direct forms, n = z-hat split) are evaluated in a
  per-direction frame (e1,e2,n); ray functionals k1, Dx1, k2, Dx2 built by spectral cumulative integration;
  k2_init (Gram-Schmidt triad) evaluated in GLOBAL components with general n then rotated into the frame;
  observer clock/position at second order re-derived (SPEC 6, eta_i -> 0):
     deta_o1 = -(1/eta_o) int_0^eta_o eta Psi1,  deta_o2/2 = -(1/eta_o)[deta_o1^2/2 + deta_o1 eta_o Psi1(eta_o) + int_0^eta_o eta f2/2],
     f2/2 = Psi2/2 - Psi1^2/2 - v1^2/2 + X1.grad Psi1,  X1 = int_0^eta v1,  dxo2/2 = int_0^eta_o [v2/2 + X1.grad v1] + deta_o1 v1(eta_o).
  (checked against the exact dilation family by hand: d2/2 = eta_o(c^2/9 - a2/2) = exact series of eta_o/lambda.)
* CHECKS (`test_lowest.py`): K1 from L1 (shear and direct lists) = type-A expansion with Delta_l to 2e-15;
  L2 shear form = L2 direct form for a generic triple of momenta, all three pairs incl. the K=0 pair, to 1e-13.
* **BUG FOUND IN C2's `definitions` dict (affects any consumer, e.g. A1/A2):** it writes k1(eta) = k1_init + int_eta^eta_o G1
  and k2/2 = ... + int_eta^eta_o [G2/2 + Dx1.dG1]. The correct (and the one los.py itself uses via D k = G) is
  k(eta) = k(eta_o) - int_eta^eta_o G.  `test_sign.py`: with C2's sign, L2(shear) - L2(direct) = O(8) instead of 1e-13.
  (C2's L_terms lists themselves are fine.)

## Milestone 2 — third order: contracted fields + full L3 assembly (DONE, validated)
* `third.py`: channel (iii) fields from C1's K3c: Phi3c = int G_Phi(k) K3c_SPhi + C3 ghat, Psi3c = Phi3c - Q3c,
  delta3c, V3c = -(eta^2/4)[2(Phi3c' + Psi3c/eta) + S0z_c/(ik)] (C1 normalisation per P(q)); KPhi2/Kh supplied by the
  pair profiles (scalar second-order kernels checked symmetric to 1e-16), H''' from the h ODE.  Third-order delta_T:
  delta_T3 = delta3/4 - (9/16) delta1 delta2 + (21/64) delta1^3.
* Per configuration (k, q, -q) the abc coefficient of a third-order field is K3c(q) + K3c(-q) = X3c(mu) + X3c(-mu).
* `thirteen.py`: Th3 = m3 - 3<m1>m2 - 3<m2>m1 + 6<m1>^2 m1 with m3 = 6 list3 + 3 L1 (2 list2) + L1^3 in the algebra.
* **CHECK (decisive): C2's shear-form L3 and direct-form L3, evaluated with C1's contracted third-order fields and the engine's
  second-order functionals / clock / displacement / emission shift, agree to 2e-11 (relative 3e-14) for the phi-averaged
  coef_abc(k,q,-q;n) at mu = 0.31, -0.62 (`test_third.py frame`).**  This tests C1's third-order contracted solution
  (third-order conservation), C2's L3 lists and all of the engine's third-order objects at once.  Also `test_null2.py`:
  engine k2_z = second-order null condition to 1e-13.
* Triad subtlety found on the way: C2's direct-form *observer* terms at third order contain k2_init with the Gram-Schmidt
  triad of the n = z-hat frame (the 'split' substitution), whereas the canonical (SPEC) triad is Gram-Schmidt of the GLOBAL
  coordinate basis.  With the global triad in the ray functional and the split direct-form obs terms, direct - shear =
  11.32 (P_3 - P_1)(n.k) exactly (an observer-local mismatch, `dbg` in report history).  The shear-form L3 has no k2 at the
  observer, so PRODUCTION = shear form + global (canonical) triad.  The per-direction GS triad is not a rigid rotation and
  changes kappa at the 0.3% level; the canonical triad is rigid, which leaves Re C^(13) invariant (|a_lm|^2 sum invariant).
* Bookkeeping bug fixed: duplicate term names in C2 lists (engine channel dict now accumulates).
