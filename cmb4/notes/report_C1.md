# Report C1 — Einstein + fluid equations to O(eps^3), Poisson gauge, RD (a = eta)

Code: `cmb4/code/C1_einstein/`.  Conventions: SPEC.md §1, §4, §5 exactly.
All fields X = eps X1 + eps^2 X2/2 + eps^3 X3/6; the engine stores the coefficient of eps^n
(so an "order-n component" of an equation E is E^(n)/n!).

## Milestone 1 — ENGINE (status: COMPLETE, verified)
- `tpoly.py`: sparse truncated polynomials in "jet symbols" J(field, order, multi-index), coefficients
  rational x eta^n; eps-order carried by each jet, products truncated at eps^3 at multiplication time
  (no higher orders ever generated).  Series inverse / sqrt / log1p; derivatives d_eta, d_x, d_y, d_z;
  conversion to sympy (Derivative(Function)) for LaTeX; pickling with jet registry.
- `einstein.py`: g = eta^2 ghat with ghat_00 = -(1+2Psi), ghat_0i = omega_i, ghat_ij = (1-2Phi)delta_ij + h_ij
  (omega, h start at eps^2; all six h_ij kept independent here, TT imposed at the decomposition stage).
  Inverse metric as truncated Neumann series (checked g g^-1 = 1 to eps^3); Christoffels, Ricci, G_{mu nu}
  computed DIRECTLY from g (not via conformal formula); fluid 8 pi G T_{mu nu} = (3/eta^4)(1+delta)[(4/3) u u + (1/3) g],
  u^mu = eta^-1 uhat^0 (1, v^i), uhat^0 from the normalisation as a series (checked g u u = -1 to eps^3).
  Also nabla_mu T^{mu nu} as a series.
- CHECK: background (order-0) part of every E_{mu nu} and of the conservation equations vanishes identically.
- CHECK: first-order equations reproduce the standard linearised equations:
  E_00^(1) = 2 lap Phi1 - 6 (Phi1' + Psi1/eta)/eta - 3 delta1/eta^2,
  E_0i^(1) = 2 d_i(Phi1' + Psi1/eta) + 4 v1_i/eta^2,
  E_ij^(1) traceless: d_i d_j (Phi1 - Psi1),  trace: 2Phi1'' + ... (see E.tex).
- Outputs: `E.pkl` (engine format, incl. Gam, G, T, cons, uhat0), `E_sympy.pkl` (srepr strings of every
  component at orders 1,2,3), `E.tex` (LaTeX dumps).  Term counts: E_00: [6, 36, 127], E_0i: [3, 18, 83],
  E_ii: [9, 40, 147], E_ij (i!=j): [2, 19, 89] at orders [1,2,3].  Runtime < 1 s.

## Milestone 2 — FIRST ORDER (status: in progress)
Status: COMPLETE, verified (`first_order.py`, outputs `first_order.pkl/.tex`).  Frame k = k zhat.
- Traceless-scalar ij:  -k^2 (Phi1 - Psi1) = 0  => Psi1 = Phi1 (no anisotropic stress).
- 0i transverse: 4 v1_{x,y}/eta^2 = 0 => no transverse velocity; 0i longitudinal: v1_i = d_i V1,
  V1 = -(eta^2/2)(Phi1' + Psi1/eta)  (irrotational, curl v1 = 0 identically).
- 00: delta1 = -(2/3) eta^2 k^2 Phi1 - 2 eta Phi1' - 2 Psi1.
- trace ij (after eliminating delta1, Psi1):  Phi1'' + (4/eta) Phi1' + (k^2/3) Phi1 = 0   [VERIFIED].
- Conservation equations (nu = 0..3) vanish identically on the solution of the Einstein equations [VERIFIED].
- Regular solution, normalised with zeta_1 = lim (delta_T1 - Phi1), x = k eta/sqrt3:
    T_Phi(x)    = -2 (sin x - x cos x)/x^3                      -> -2/3
    T_delta(x)  = -4 cos x + 8 sin x/x + 8 cos x/x^2 - 8 sin x/x^3   (delta1 = zeta T_delta)
    T_deltaT(x) = T_delta/4 = -cos x + 2 sin x/x + 2 cos x/x^2 - 2 sin x/x^3  -> 1/3
    T_V(x)      = eta (x^2 sin x + 2 x cos x - 2 sin x)/x^3    -> eta/3   (v1_i = d_i V1, V1 = zeta T_V)
  all verified by substitution; singular second mode (cos x + x sin x)/x^3 ~ eta^-3 discarded.
  Agrees with lead's linear_check (T_Phi, delta_T -> zeta/3, v_pot -> eta zeta/3).

## Milestone 3 — SECOND ORDER (status: in progress)
