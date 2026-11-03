# report_V2 — squeezed-limit predictions from the adiabatic-mode transformation (independent route)

Status (2026-11-03): items (1), (2), (3), (4) DONE analytically; sympy checks PASS; numerical evaluator written
(code/V2_squeezed/squeezed_predictions.py). Not done: explicit O(qS) Poisson-gauge field kernels after gauge
restoration (recipe only, Sec. 2), and no comparison with other agents' K_2 (deliberately not read).
Code: `code/V2_squeezed/adiabatic_mode_rd.py` (script 1, inherited from the crashed run, re-run 2026-11-03: ALL CHECKS PASS),
`code/V2_squeezed/squeezed_predictions.py` (script 2, new).

Conventions: SPEC. Observer label point = origin; ray x(eta) = n(eta_o - eta); K_1(k,n) = F(k, mu) - F_0(k), mu = khat.n,
F_0 = sky monopole. zeta_L(x) = zeta_{k1} e^{i k1.x} = c + q.x + O(k1^2) with c = zeta_{k1}, q = i k1 zeta_{k1}.
Theta = eps Th1 + eps^2 Th2/2 + ..., Th_r = int zeta^r K_r  =>  the eps^2 long x short term of Theta is int int zeta_L zeta_S K_2.

## 1. Weinberg mode in RD Poisson gauge (script 1, PART A/B/D) — verified
Linearised G - 8piG T derived from scratch (a = eta, linear parametrisation). The finite-k RD solution
T_Phi = -2(sin x - x cos x)/x^3, x = k eta/sqrt3, T_dT = -(eta^2/6)(k^2 T_Phi + 3(T_Phi' + T_Phi/eta)/eta), T_v = -(eta^2/2)(T_Phi' + T_Phi/eta)
satisfies all components (residual 0). For zeta_L = c + q.x:
  Phi_L = Psi_L = -(2/3) zeta_L,  delta_T,L = zeta_L/3,  v_L^i = (eta/3) q^i   (exact solution, time independent),
and it is the pullback of FRW under
  eta~ = eta - (eta/3) zeta_L(x),   x~ = (1 + zeta_L) x - (x^2/2) q - (eta^2/6) q     (f = eta^2/6, f' = a^-2 int a^2 = eta/3),
fluid at rest and T = Tbar(eta~) in the tilde chart (all checked, O(c), O(q)). At eta -> 0 the spatial map is conformal,
d x~^2 = (1 + 2 zeta_L) dx^2 (checked). Constant mode to all orders in c: x~ = e^c x, eta~ = lambda(c) eta with
lambda^4 + lambda^2 e^{2c} = 2 keeps Poisson form, lambda = 1 - c/3 - 13c^2/54 + ...  (only needed for field-level checks).

## 2. Transformation of short-mode fields (script 1, PART C/E) — verified
Image of a Poisson short-mode solution under the map: cross terms beyond the argument shift F(eta~, x~):
  Phi: +(4/3) zeta_L Phi_S,   Psi: -(4/3) zeta_L Psi_S,   delta_T: +(1/3) zeta_L dT_S,
  v^i: -(4/3) zeta_L v_S^i + [q^i (x.v_S) - x^i (q.v_S)]  (local SCT rotation),
  ghat_0i: +(2 eta/3) q_i (Phi_S + Psi_S)  [O(qS) only; no g_0i at O(cS)].
Primordial field: on the uniform-density slice (eta -> 0) the invariant metric a(rho)^2 e^{2 zeta_nl(x)} dx^2 equals
a(rho)^2 e^{2 zeta~(x~)} dx~^2, and dx~^2 = e^{2 zeta_L} dx^2, hence EXACTLY  zeta~_S(x~) = zeta_S(x(x~))  (scalar relabelling).
Field-level predictions (constant mode, Poisson gauge, linear parametrisation; X2 cross kernel per unit c, k = k_S):
  F2_Phi(0,k) = (4/3)(1 - eta d_eta) T_Phi,   F2_Psi(0,k) = -(4/3)(1 + eta d_eta) T_Phi,
  F2_dT(0,k)  = (1/3) T_dT - (4/3) eta d_eta T_dT,   F2_vpot(0,k) = -(4/3)(1 + eta d_eta) T_v   (v^i = d_i vpot).
  These satisfy zeta_nl^(2) cross term = 0 at eta -> 0 (checked), and agree with CPV (49)-(52) for f = -1/3 after converting
  exponential -> linear parametrisation. O(qS): the longitudinal part of ghat_0i must be removed by a time shift
  alpha = -(4 eta/3) nabla^-2 (q.grad Phi_S); the transverse remainder is the squeezed vector mode
  omega_i = (4 eta/3) T_Phi(k eta) zeta_k [q_i - k_i (q.k)/k^2]  (recipe; not cross-checked against the 2nd-order vector eq.).

## 3. Observable: Theta[zeta_L + zeta_S] in terms of Theta[zeta_S]  (MAIN RESULT)
Key exact fact (SPEC Sec. 3): T_obs = T_o exp(-S), S = int sigma_ee dl. Every ingredient is covariant: shear tensor, null
geodesic, emission surface T_phys = T_*, tagged fluid worldline, proper-time clock from the T_phys = T_i event, monopole
normalisation; the tetrad enters only through its rotation part. Therefore, with the long mode removed by the map above,
   Theta[zeta_L + zeta_S](n) = Theta[zeta~_S ; observer at x~_o](R n),
   x~_o = -(eta_i^2/6) e^{-c} q  (label 0 at eta_init -> tilde position; fluid at rest in tilde chart; vanishes for eta_i -> 0),
   R = 1 + O(q^2, qS)  (antisymmetric Jacobian q_j x^i - x_j q^i vanishes on the worldline x_o = (eta^2-eta_i^2) q/6 || q; checked).
Note S[zeta_L] = 0 identically through O(q) (FRW has no shear), so <m1_L> = 0 and the normalisation produces no L x S term.
Plane-wave rule Theta^(1)[P(x) e^{ip.x}] = P(-i d_p) K_1(p) turns the relabelling into operators (sympy-verified):

### PREDICTION for lim K_2 (through O(k1)):
   K_2(k1, k2; n) = - D K_1(k2,n)  -  S_{k1} K_1(k2,n)  +  (eta_i^2/6)(k1.k2) K_1(k2,n)  +  O(k1^2)
   D = k2.d_{k2} = d/dln k2 (fixed khat2, n);   S_p = (k2.d_{k2})(p.d_{k2}) - (1/2)(p.k2) Lap_{k2}.
In variables K_1 = F(k, mu) (monopole-subtracted), a = k1.khat2, b = k1.n:
   S_{k1} F = (a k/2) F_kk - a F_k + (b - mu a) F_kmu + ((2 a mu - b)/k) F_mu - (a (1-mu^2)/(2k)) F_mumu .
Since RD has no scale other than eta, F = F(k eta_o, k eta_*, mu), so D = eta_o d_{eta_o} + eta_* d_{eta_*}:
   O(q^0):  Theta_LS = -c (eta_o d/d eta_o + eta_* d/d eta_*) Theta^(1)_S(n)      [operator on the map; exact for k1 = 0]
   O(q^1):  Theta_LS = Theta^(1)[ delta zeta ](n) - (eta_i^2/6) q.grad_{x_o} Theta^(1)_S(n),
            delta zeta(x) = -[(q.x) x - (x^2/2) q].grad zeta_S(x)   (SCT relabelling of the primordial field).
Real-space form of both:  Theta_LS(n) = Theta^(1)[ -xi_0.grad zeta_S ](n) + observer-label term, xi_0 = zeta_L x - x^2 q/2.

### Cancellations (what is absent)
O(q^0): the separate CPV-type pieces — amplitude rescaling (the "1" term), emission-time shift eta_e = eta_*(1 + c/3)
(from T_phys = T_*: 1 + delta_T,L = eta_e/eta_*), observer clock eta_o -> eta_o(1 + c/3), ray coordinate length
(1 - 4c/3)(eta_o - eta_e), field cross terms (Sec. 2), monopole normalisation — do NOT survive individually; they sum to the
pure relabelling zeta_S(x) -> zeta_S(e^{-c} x), i.e. -c d/dln k2. There is no leftover "1 + d/dln eta_* - 5 n.grad_n" (CPV 27).
O(q^1): observer velocity v_o = eta_o q/3 (aberration/Doppler), observer drift x_o = (eta_o^2 - eta_i^2) q/6, time-dependent
translation f q, emission-surface tilt, gradient lensing/time-delay, the O(qS) vector mode and the SCT rotation of v all cancel;
no dipolar modulation survives except (i) the SCT relabelling S_{k1} of the primordial Gaussian field and (ii) the trivial
eta_i^2 observer-label translation. Long mode alone: K_1(k,n) = O(k^2) (checked symbolically: monopole and dipole of K_1 vanish
at k^0, k^1; even before monopole subtraction, S_1[c] = 0).
IMPORTANT for SPEC 9(d): K_2(0,k2) is NOT zero; it is -dK_1/dln k2. "Pure gauge" means the map equals the map of the
RELABELLED short field; for fixed Gaussian zeta_S this is a physical dilation (Maldacena's relation is violated by an exactly
Gaussian zeta unless n_s = 1). Residual physical effects start at O(q^2) (MZ 5.13: the q^2 terms are curvature, not gauge).

### IR consequence for the one-loop C_l (derived, not numerically tested)
Soft region of C^(22) (two regions, k1 or k2 soft) and of C^(13) (q soft) with K_3(q,-q,k) = D^2 K_1 + O(q^2) (Sec. 4):
  [C^(22)/4 + Re C^(13)/3]_IR = sigma_L^2 x (1/(2l+1)) sum_m int_k P(k) [ |D K_lm|^2 + Re K*_lm D^2 K_lm ]
                              = sigma_L^2 x (1/(2l+1)) sum_m (1/2) int d^3k/(2pi)^3 |K_lm|^2 (3 + D)^2 P(k),
  sigma_L^2 = int_soft d^3q P(q)/(2pi)^3. For scale-invariant P ∝ k^-3, (3 + D) P = 0: the log IR divergence cancels between
  C^(22) and C^(13). O(k1) terms drop by parity; O(k1^2) is IR finite.

## 4. Third order: relations among K_3, K_2, K_1
Composite map (c then q): kernel = (1 + i S_q) e^{-cD}[e^{ik.x~_o} K], using e^{-cD} S_q = e^{c} S_q e^{-cD} (checked).
  (a) all orders in c:      K_{m+1}(0,...,0,k) = (-D)^m K_1(k);   in particular K_3(0,0,k) = D^2 K_1(k).
  (b) two long, c^2 and cq: K_3(ka, kb, k) = D^2 K_1 + (S_{ka} + S_{kb}) D K_1 - (eta_i^2/6) k.(ka + kb) (2 + D) K_1 + O(ka^2, kb^2, ka kb)
      => K_3(q, -q, k) = D^2 K_1(k) + O(q^2)   (the C^(13) soft limit).
  (c) one long x two short: K_3(k1, ka, kb) = -(D_a + D_b) K_2(ka,kb) - (S^(a)_{k1} + S^(b)_{k1}) K_2(ka,kb)
                                             + (eta_i^2/6) k1.(ka + kb) K_2(ka,kb) + O(k1^2).
  O(ka kb), O(q^2) pieces are physical (not fixed by the gauge argument).

## 5. Numerics (squeezed_predictions.py)
K_1 built from the SPEC 9(a) linear formula (emission Theta = dT + Psi - n.v at eta_*, + n.v_o, + 2 int Phi', minus monopole;
monopole of the local and ISW parts done analytically with j0, j1); derivatives in k, mu by sympy, eta-integrals by mpmath.
Functions: K1, DK1, SK1(a,b,...), K2_squeezed(k1vec,k2vec,n,eta_o,eta_*,eta_i) -> (O(k1^0) part, O(k1^1) part, K_1).
Self-checks: D K_1 vs eta-scaling finite difference; (k,mu) form of S_p vs direct 3D finite differences. Table below.
