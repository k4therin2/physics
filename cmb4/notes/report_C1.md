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
Status: COMPLETE, all three control checks PASS (`second_order.py`, `second_order_solve.py`, `bianchi.py`,
`gauge_check.py`, `position_sources.py`; outputs `second_order.pkl/.tex`, `position_sources.pkl/.tex`, `gauge_check.pkl`).
Notation: E^{(n)} := n! [eps^n] E_{mu nu} = Lin[X_n] + S^{(n)};  kernels per zeta_q zeta_p, p = k - q,
X2_k = int d^3q/(2pi)^3 K_X(q,k-q;eta) zeta_q zeta_{k-q};  frame k = k zhat, q = (q sqrt(1-mu^2),0,q mu);
invariant form uses magnitudes (k,q,p), mu = (k^2+q^2-p^2)/(2kq). T := T_Phi, T' := d_eta T_Phi.

### Linear operators (identical at every order) and sources
- Traceless scalar: Q_n := Phi_n - Psi_n = S_TL/k^2,  S_TL = S_zz - (S_xx+S_yy)/2  (position space:
  Q_n = -(3/2) nabla^-2 (d_i d_j nabla^-2 - delta_ij/3) S_ij).
- 00:  delta_n = (eta^2/3) S_00 - (2/3) eta^2 k^2 Phi_n - 2 eta Phi_n' - 2 Psi_n.
- 0i longitudinal: v_{n,i} = d_i V_n,  V_n = -(eta^2/4)[2(Phi_n' + Psi_n/eta) + S_0z/(ik)].
- 0i transverse (constraint relating nabla^2 omega to the transverse momentum density):
  (1/2) nabla^2 omega_i - (4/eta^2)(omega_i + v^T_i) = S^T_{0i}   [8 pi G a^2 (rho+P) = 4/eta^2].
- ij vector: omega_n' + (2/eta) omega_n = (2/(ik)) S_xz  -> regular solution omega = eta^-2 int eta'^2 S.
- ij TT:  h'' + (2/eta) h' + k^2 h = -2 (S_xx - S_yy)  for the (xx - yy) component; h2_ij = Lambda_ij,kl q_k q_l x (...).
- trace (delta_n eliminated by 00, Psi_n = Phi_n - Q_n):
  Phi_n'' + (4/eta) Phi_n' + (k^2/3) Phi_n = S_Phi := -(S_tr - S_00)/6 + Q_n'/eta - (k^2/3) Q_n.
Position-space S^{(2)}_{mu nu} (explicit derivatives, 6-20 terms each) and S^{(3)} (71-127 terms) in position_sources.tex/.pkl.

### Second-order kernels (abstract T form; closed forms obtained by inserting T_Phi)
- K_Q(q,p) = [ (13k^4-6k^2(p^2+q^2)+9(p^2-q^2)^2) T(q)T(p) + ... ]/(4k^4) (full form in second_order.tex); mu-form:
  K_Q = T(q)T(p)[-4+6 mu q/k+3q^2/k^2-9 mu^2 q^2/k^2] + eta[T(p)T'(q)+T(q)T'(p)][2 mu q/k + q^2/k^2 - 3 mu^2 q^2/k^2]
        + eta^2 T'(q)T'(p)[2 mu q/k + q^2/k^2 - 3 mu^2 q^2/k^2].
  Superhorizon: K_Q(0) = -(4/9)(4k^2 - 6 k q mu + 9 mu^2 q^2 - 3 q^2)/k^2  = (-4 + 9 F/k^2) Phi1Phi1-units, F = BMR's F.
- S_Phi kernel K_SPhi (second_order.tex), finite as eta->0 (no 1/eta^2 terms), so the particular integral vanishes at eta=0.
- VECTOR SECTOR (important finding): curl(uhat_i) = 0 is NOT implied at second order; what the equations imply
  (verified: the algebraic omega2 obtained from it satisfies the ij-vector ODE identically) is relativistic Kelvin:
  curl[(1+delta_T) uhat_i] = 0 (enthalpy-weighted velocity, h ∝ T ∝ rho^{1/4}).  Consequently omega2 is ALGEBRAIC:
    omega2_i(k) = int P_ij(k) q_j K_om(q,p) zeta zeta,
    K_om = (2i/k^2)[ eta (q^2-p^2) T(q)T(p) + eta^2 (q^2 T(q)T'(p) - p^2 T(p)T'(q)) ]   (antisymmetric in q<->p),
    v2^T_i = P_ij q_j K_vT,  K_vT = -(eta^2 k^2/8 + 1) K_om - (eta^2/4) S^T_{0i}-kernel  (explicit in second_order.pkl).
- TENSOR: h2_ij(k) = Lambda_ij,kl(k) q_k q_l int_0^eta G_h K_h,  K_h = 4[ 3 T(q)T(p) + eta(T'(q)T(p)+T(q)T'(p)) + eta^2 T'(q)T'(p) ]
        = 4[ 2 T(q)T(p) + (eta T'+T)_q (eta T'+T)_p ].
  CHECK (a): identical to Kohri-Terada (8) / Baumann et al (11)-(12) with w=1/3, Phi=Psi, h2^{ours} = h^{KT}: PASS (numerically exact).
- GREEN'S FUNCTIONS (verified L G = 0, G(eta,eta)=0, d_eta G|_{eta'=eta}=1):
  Phi: u1 = (sin x - x cos x)/x^3, u2 = (cos x + x sin x)/x^3, W = -3 sqrt3/(eta^4 k^3),
       G_Phi(eta,eta') = [u1(eta')u2(eta) - u1(eta)u2(eta')]/W(eta') ;  Phi2 = int_0^eta G_Phi S_Phi + C_Phi(k) ghat(k eta), ghat = 3u1 -> 1.
       (k->0: G_Phi -> (eta'/3)(1 - eta'^3/eta^3)).
  h:   G_h = (eta'/eta) sin(k(eta-eta'))/k.      omega: algebraic (or G_om = (eta'/eta)^2 for the ODE form).
  delta2, V2, Psi2: algebraic from the constraints (expressions saved as K_delta2_expr, K_v2z_expr).
- CHECK (b) Bianchi: engine verifies g^{mu a} nabla_a E_{mu nu} + g_{nu b} cons^b == 0 identically to O(eps^3)
  (bianchi.py).  Hence d/deta of the constraints (00, 0i) at any order = combination of the evolution equations,
  conservation, and lower-order equations; with first-order fields on shell the second-order constraints propagate.
- CHECK (c) constant zeta: exact family eta = A eta~, x = B x~: Psi=(A^4-1)/2, Phi=(1-A^2B^2)/2, delta=A^-4-1, v=0,
  zeta_nl = ln B; engine verifies it solves E=0 to O(eps^3) for arbitrary A,B.  With b1 = c, a1 = -c/3 (Psi1=Phi1 = -2c/3,
  delta1 = 4c/3), the k=0 second-order equations give ONLY delta2 = -2 Psi2 + 32c^2/9 (Phi2, Psi2 free: residual
  zero-mode gauge freedom of Poisson gauge).  zeta^(2)=0 fixes b2 = c^2.  The member with Psi2 = Phi2 has
  Phi2 = Psi2 = -8c^2/27, delta2 = 112c^2/27 (delta_T2 = 19c^2/27).  The angle-averaged k->0 limit of the finite-k
  regular solution instead picks <Phi2(0)> = -(8/9)c^2, <Psi2(0)> = +(8/9)c^2, delta2 = 16c^2/9: on the same gauge family.

## Milestone 4 — INITIAL CONDITIONS FROM GAUSSIAN ZETA (status: COMPLETE)
zeta_nl = (1/4) ln(1+delta) + (1/2) ln(1-2Phi) expanded by the engine:
  zeta^(2) = delta2/4 - delta1^2/4 - Phi2 - 2 Phi1^2,
  zeta^(3) = delta3/4 - (3/4) delta1 delta2 + delta1^3/2 - Phi3 - 6 Phi1 Phi2 - 8 Phi1^3.
lim eta^2 K_S00 = 32/3, K_Q(0) above, quadratic kernel of zeta^(2) -> -3 TT(0)... ; imposing zeta^(2)(eta->0) = 0:
  Phi2(k, eta->0) = C_Phi(k) with kernel
     K_Phi2(q,p; 0) = -(4/9) (2k^2 - 2 k q mu + 3 mu^2 q^2 - q^2)/k^2
                    = -(7k^4 - 2k^2(p^2+q^2) + 3(p^2-q^2)^2)/(9 k^4)
  position space (Phi1 = -2zeta/3): Phi2(0) = -2 Phi1^2 - nabla^-2(d_i Phi1 d_i Phi1) + 3 nabla^-4 d_i d_j (d_i Phi1 d_j Phi1),
  Psi2(0) - Phi2(0) = 4 Phi1^2 + 3 nabla^-2(dPhi dPhi) - 9 nabla^-4 d d(dPhi dPhi).
COMPARISON: EXACT agreement with Pitrou-Uzan-Bernardeau (3.12b,c) with R^(2) = -2 R^(1)2 (Gaussian Salopek-Bond zeta),
  after swapping their (Phi,Psi) names (PUB/BMR call the lapse Phi).  Bartolo-Matarrese-Riotto II (7.30): converting
  their exponential variables, BMR's Q2(0) = (-4 + 33 F/k^2) Phi1Phi1 vs ours (-4 + 9 F/k^2): the gradient terms (6F) and
  local term agree; the velocity-quadrupole term differs (BMR 27F, ours 3F) because BMR evaluate it with the late-time
  oscillating velocity (their 7.29) rather than the eta->0 velocity v = -(eta/2) grad Phi(0).  We keep ours (= PUB).

## Milestone 5 — THIRD ORDER (status: in progress)
Status: COMPLETE for the contracted ("13") form; position-space full sources COMPLETE; the un-contracted Fourier
kernel K3(k1,k2,k3) was NOT produced (not needed for C^(13)).  Scripts: `third_order.py`, `third_order_ic.py`,
`squeezed_check3.py`; outputs `third_order.pkl/.tex` (contracted sources), `third_order_ic.pkl/.tex` (tadpoles, C3),
`position_sources.pkl/.tex` (full S^{(3)}_{mu nu} in position space: 71-127 terms per component).

### Contracted construction (per P(q); NOT symmetrised under q -> -q, so only p = |k-q| appears; symmetrise by averaging with q->-q if desired)
  X3c(k,eta) = int d^3q/(2pi)^3 P(q) K3c_X(k,q,mu;eta),   <X3_k zeta_k'> = (2pi)^3 delta(k+k') P(k) X3c.
  cubic monomial  A1 B1 C1  ->  A(q)B(-q)C(k) + A(q)C(-q)B(k) + B(q)C(-q)A(k)
  bilinear A1 B2            ->  2 A(q) B2[(q1,q2)=(-q,k); total k-q]  +  A(k) <B2>(eta)
  Second-order fields entering: Q2, delta2, V2, v2^T, omega2 explicit (second_order.pkl), Phi2 = KPhi2(a,b,c,eta) and
  h2 = Kh(a,b,c,eta) kept as abstract symbols ((a,b,c) = (|q1|,|q2|,|q1+q2|)) with
    KPhi2(a,b,c,eta) = int_0^eta G_Phi(c;eta,eta') K_SPhi(a,b,c;eta') deta' + K_Phi2_0(a,b,c) ghat(c eta),
    Kh(a,b,c,eta)    = int_0^eta G_h(c;eta,eta') K_h(a,b,c;eta') deta'.
  TADPOLES (mean-field shifts <X2>): vector/tensor/velocity tadpoles vanish by isotropy; scalar tadpoles are the
  angle-averaged k->0 limits (this is a CONVENTION for the zero mode; Phi2,Psi2 at k=0 are gauge, see check (c)):
    tadQ2(q,eta) [closed form in third_order_ic.tex], tadQ2(0) = -16/9;  tadS00 (local), eta^2 tadS00 -> 32/3;
    tadPhi2(eta) = -8/9 + int_0^eta (eta'/3)(1-eta'^3/eta^3) tadSPhi deta'  (closed form; = -8/9 - 28 q^2 eta^2/405 + ...);
    taddelta2(eta) closed form, -> 16/9.   These agree with the k=0 gauge family (delta2 = -2 Psi2 + 32c^2/9).
### Equations for X3c (same linear operators): Q3c = K3c['Q'];  delta3c = (eta^2/3) K3c['S00'] - (2/3)eta^2 k^2 Phi3c - 2 eta Phi3c' - 2 Psi3c;
  V3c = -(eta^2/4)[2(Phi3c' + Psi3c/eta) + K3c['S0z']/(ik)];  Phi3c'' + 4Phi3c'/eta + k^2 Phi3c/3 = K3c['SPhi']  (830 terms),
  Phi3c = int_0^eta G_Phi K3c['SPhi'] + C3(k;q,mu) ghat(k eta).   Contracted vector/tensor sources (S0x, Sxz, Sxx-yy) are
  given too but integrate to zero over angles by isotropy (S0y, Syz, Sxy vanish identically: checked).
  K3c['SPhi'] is finite at eta -> 0 (checked with small-eta expansions), so the particular integral vanishes at eta=0.
### Third-order amplitude from zeta_nl^(3)(eta->0) = 0
  zeta^(3) = delta3/4 - (3/4) delta1 delta2 + delta1^3/2 - Phi3 - 6 Phi1 Phi2 - 8 Phi1^3  (engine), contracted (100 terms);
  with Q3c(0) = -16 q (44k^4 mu - 142 k^3 mu^2 q - 18 k^3 q + 108 k^2 mu^3 q^2 + 104 k^2 mu q^2 + 9 k mu^4 q^3 - 126 k mu^2 q^3 - 3 k q^3 + 24 mu q^4)/(9 k p^4),
  lim eta^2 K3c_S00 = -256 k q (2k^2 mu - k mu^2 q - 3kq + 2 mu q^2)/(3 p^4):
  C3(k;q,mu) = Phi3c(k,eta->0) = -16 (2k^5 + 12k^4 mu q - 42k^3 mu^2 q^2 - 10k^3 q^2 + 36k^2 mu^3 q^3 + 32k^2 mu q^3 + 3k mu^4 q^4 - 42k mu^2 q^4 + k q^4 + 8 mu q^5) / (9 k p^4),
  p^2 = k^2 - 2kq mu + q^2 (invariant form and angle average with logs in third_order_ic.tex; integrable 1/p^2-type behaviour at q->k).
  CHECK: squeezed limit C3(q->0) = -32/9 = 6 x [c^2 zeta_s](A^2B^2)(-2/3) from the exact dilation family (A^2B^2 = e^{4c/3}): PASS.

## Summary of status
COMPLETE+VERIFIED: engine; E_{mu nu} to eps^3 (background, Bianchi identity, exact gauge family); first order (ODE,
T_Phi, T_delta, T_V, zeta normalisation, conservation implied, irrotational); second order (all SVT equations, kernels,
Green's functions, KT tensor source, PUB ICs, vorticity-free structure [enthalpy-weighted], gauge check); third order
contracted sources and C3 (squeezed check).  PARTIAL: tadpole zero-mode convention (angle-averaged) is a choice to be
validated by the full observable's dilation invariance (SPEC 9d).  NOT DONE: un-contracted K3(k1,k2,k3); closed-form
evaluation of the Green's-function integrals (left as declared integrals, cf. SPEC §0).
How to use: `from kernels_api import *; K = load()`; K['second'][...] (K_Q, K_SPhi, K_om, K_vT, K_h, K_Phi2_0, G_Phi, G_h, ...),
K['third']['SPhi','Q','S00','S0z',...], K['third_ic']['C3','tads']; `closed_form()` inserts T_Phi.
