# Report A2 — second-order map kernel K_2 and C_l^(22)
Agent A2 (relaunch 2026-11-03). Code: `cmb4/code/A2_second_order/`. Conventions: SPEC.md + C2 (ray xbar = x_o + n chi,
chi = eta_o - eta, n = observation direction, fields on the ray carry e^{+i k.n chi}).

## Milestone 1 — engine repaired, plane-wave check (iv) vs N1: PASS
Reused from the crashed run: `transfer.py` (first-order transfer functions, C1 second-order kernels as numpy callables),
`kernel2.py` (numerical assembly of K_1, K_2 from C2's L_terms.pkl term lists, n-frame representation, Chebyshev time grid,
Gauss-Legendre Green's integrals).  Two defects found and fixed:
1. **Sign of the ray deflection (affects C2's `definitions`).**  The photon direction perturbation must be
   k1(eta) = k1_init - int_eta^{eta_o} G1 deta'   (since dk/deta = G and the ray is integrated backwards from the observer),
   i.e. k1_perp(eta) = i k_perp [T_V(eta_o) + 2 int_eta^{eta_o} T_Phi e^{i k.n chi'} deta'].  C2's `ray.py` line 93
   (`definitions['k1']`, and likewise `definitions['k2half']`, line 96-99) writes **+** int G1.  The L-term lists themselves
   are unaffected (they keep k1, Dx1 abstract), but any consumer of `definitions` (A3!) must flip the sign.  With the C2 sign
   the second-order maps disagree with N1 at O(1) and the shear and direct forms disagree with each other; with the
   corrected sign both agree with N1 (below).  The radial part n.k1 = -(Phi+Psi) is the null condition and is unaffected.
2. Psi2' was taken by finite differences of Q2 on the grid; replaced by the exact eta-derivative of C1's K_Q (needs T_Phi'').
Observer dipole in Delta_l: L1 contains +(n.v)_o, v_o = i k T_V(eta_o) zeta_k  =>  Delta_1 gets +k T_V(eta_o)/3 (j_1'(0) = 1/3);
confirmed by the N1 x_o = pi/2 dipole a_10[Th1] (2.0812147132, both).

CHECK (iv): zeta = cos(k(x + x_o)), k = 1 along x, eta_* = 2, eta_o = 4, 8, x_o = 0, pi/2 (N1's four configurations);
Th2(n) = sum_{s1,s2} c_s1 c_s2 K_2(s1 k, s2 k; n) incl. the K = 0 pairs (tensor zero mode h_ij = traceless part, scalar zero
modes only monopole/dipole).  a_l0 (l = 1..12) of Th1 and Th2 vs N1 (`check_planewave.py`, `check_planewave.log`):
max |diff| = 1.2e-8 (to8_xo0), 1.0e-9 (to4_xo0), 6.5e-9 (to8_xo0.5pi), 7.0e-10 (to4_xo0.5pi) — all within N1's quoted
uncertainties (3.2e-8, 9.7e-10, 9.2e-6, 5.4e-9).  E.g. to8_xo0: a_20[Th2] = -4.2224942187 (N1 -4.2224942193),
a_40 = +1.3785439772 (N1 +1.3785439781).  Shear and direct forms of L2 agree with each other to 1e-13.

## Milestone 2 — generic-momentum checks (iii), symmetry, form equivalence, constant mode (i)
`check_generic.py` / `.log`, `check_squeezed.py` / `.log`; eta_* = 2, eta_o = 8, sky grid 40x80 GL, l <= 30.
* (A) shear form vs direct form of L2 for generic non-collinear k1 = (0.7,-0.2,0.4), k2 = (-0.3,0.9,0.5): max rel. diff 6e-15
  pointwise in n.  The two forms differ by total derivatives and the 2nd-order conservation equations, and the shear form
  contains only v2 and h2 among the second-order fields while the direct form contains Phi2, Psi2, dT2, V2, v2T, om2, h2:
  this is a strong check of C1's K != 0 kernels (incl. the vector and tensor sectors, which the plane wave cannot probe).
* (B) symmetry K_2(k1,k2) = K_2(k2,k1): 9e-15.  (C) REALITY K_2(-k1,-k2;n) = K_2(k1,k2;n)^*: exact (0.0).
* (i) CONSTANT MODE — correction to the brief.  K_2(k1 -> 0, k2; n) does NOT vanish (neither channel-wise nor in total) and
  should not: zeta_nl = c + zeta_s is the coordinate rescaling x~ = e^{c} x of the configuration zeta~_s(x~) = zeta_s(e^{-c}x~)
  (early-time spatial metric e^{2 zeta} dx^2; observer at x = 0, its proper time and T_* are invariant), so
      lim_{k1->0} K_2(k1,k2;n) = - k2 . grad_{k2} K_1(k2;n)      (Maldacena-type consistency relation for the map),
  i.e. the constant mode is unobservable only up to relabelling the short wavenumbers (which leaves C_l invariant for a
  scale-free spectrum but not mode by mode).  Verified for the FULL map, all multipoles (incl. l = 0, 1):
  |K_2(s e, k2) - (-k dK_1/dk)| (rms over l>=2) = 3.18e-1, 3.18e-2, 3.18e-3 for s = 1e-2, 1e-3, 1e-4 (|k2| = 1.07; the
  dilation itself has rms 11.25), likewise for |k2| = 2.44 — the residual is exactly linear in s (the gradient term).
  Individual channels are O(10) at s -> 0 (e.g. [ray] 2nd:V2 alone: 11.0; [emi] time delay: 13.9) and cancel to the dilation.
* (ii) gradient O(k1) vs V2 (report_V2 §3, eta_i -> 0): prediction K_2 = -D K_1 - S_{k1} K_1 + O(k1^2), D = k.d_k,
  S_p = (k.d)(p.d) - (1/2)(p.k) Lap_k (SCT relabelling).  `check_gradient.py`/`.log`: k1-odd part [K_2(se)-K_2(-se)]/(2s)
  vs -S_e K_1 with S applied to the engine's K_1 by a finite-difference Hessian: l=1 and l>=2 agree to 1.5e-5 .. 6e-5
  (|k2| = 1.07) and 1.6e-4 .. 2.2e-4 (|k2| = 2.44) per multipole, at the level of the FD-Hessian error: PASS.
  So through O(k1) the squeezed limit of the full map is pure relabelling (dilation + SCT), as V2 predicts.

## Milestone 3 — explicit kernels (deliverables 1, 2)
### K_1, Delta_l, <m1>  (direct form of L1 = standard SW + Doppler + ISW incl. observer terms)
  K_1(k;n) = [T_dT + T_Phi - i k mu T_V]_{eta_*} e^{i k mu chi_*} - [T_dT + T_Phi - i k mu T_V]_{eta_o} + int_{eta_*}^{eta_o} 2 T_Phi'(k,eta) e^{i k mu chi} d eta,
  mu = n.khat, chi = eta_o - eta.  Rayleigh e^{i k.n chi} = 4pi sum i^l j_l(k chi) Y_lm(n) Y*_lm(khat), (i n.khat) -> +d/dx:
  K_1 = 4 pi sum_lm i^l Delta_l(k) Y_lm(n) Y*_lm(khat),
  Delta_l(k) = (T_dT + T_Phi)(k,eta_*) j_l(k chi_*) - k T_V(k,eta_*) j_l'(k chi_*) + int 2 T_Phi' j_l(k chi) d eta
               - delta_l0 (T_dT + T_Phi)(k,eta_o) + delta_l1 k T_V(k,eta_o)/3 .
  (verified: K_1 rebuilt from Delta_l, l<=40, equals the engine's L1 to 7e-15; observer dipole sign confirmed by N1 x_o=pi/2.)
  Sky monopole <m1> = int d^3k/(2pi)^3 zeta_k Delta_0(k);  <<m1>^2> = (1/2pi^2) int k^2 dk P Delta_0^2 (mono: sigma^2 Delta_0(k0)^2).
  Observer objects (x_o = 0, closed forms, x_o := k eta_o/sqrt3):
   clock      deta_o1(k) = -(1/eta_o) int_0^{eta_o} eta T_Phi d eta = (6/(k^2 eta_o)) (1 - sin x_o/x_o)    (-> eta_o/3: dilation c eta_o/3 OK)
              [from Delta tau = int a/uhat^0 = int eta (1 + Psi1) d eta fixed; the initial-surface term ~ eta_i^2 -> 0].  NOTE: this is the
              integral of Psi1 (lapse), not of delta_T1 as the brief says; the two differ for k eta_o ~ 1.
   displacement dx_o1 = i k X_o zeta_k,  X_o(k) = int_0^{eta_o} T_V d eta = (3/k^2)(2 sin x_o/x_o - cos x_o - 1)   (-> eta_o^2/6)
   emission   deta_e1 = eta_* T_dT(k,eta_*) e^{i k mu chi_*}.
### Normalisation chain
  Lambda = eps L1 + eps^2 L2/2, C2's list L[2] IS L2/2;  m2 = L2 + L1^2;  Theta2 = m2 - 2<m1> m1 (l >= 1);  hence
  **K_2(k1,k2;n) = 2 Sym[K_{L2/2}] + K_1(k1)K_1(k2) - Delta_0(k1) K_1(k2) - Delta_0(k2) K_1(k1)**,  Sym = (1/2)(k1 <-> k2).
### Channel kernels in named objects (`symbolic.py` -> `symbolic_kernels.pkl`, `symbolic_kernels.tex`)
Each channel is an explicit sympy expression in (k1,k2,mu,mu1=n.khat1,mu2,K,Kn=n.K, Pperp = k1k2(mu - mu1 mu2)), phases
E_a = e^{i k_a mu_a chi}, EK = E1E2, first-order transfer functions (TPhi, TdT, TV and eta-derivatives, slot a), ray functionals
Kf_a (aberration+deflection), Dp_a (transverse displacement), Dz_a (radial displacement = time delay + clock + displacement),
observer/emission objects (Xo, detao1, detae1) and C1's pair kernels (Phi2, Psi2, dT2, V2, K_vT, K_om, H = int G_h K_h).
Examples (ray terms are int_{eta_*}^{eta_o} d eta):
  ray time delay        -2 [Dz_1 TPhi_2 k2^2 mu2^2 E2 + (1<->2)]
  ray lensing           -2i Pperp [Dp_1 TPhi_2 k2 mu2 E2 + (1<->2)]      (= 2 grad_perp Dx1 . grad_perp (2Psi)-type, cf. C3 type D)
  ray aberr.+deflection  2 Pperp [Kf_2 TPhi_1 E1 + (1<->2)]
  ray local 1x1          4 (TPhi_1 TPhi_2' + TPhi_2 TPhi_1') E1E2
  ray 2nd-order          EK [ Phi2' - Psi2' + 2i Kn Psi2 - (1/2) h_nn'/H' * Hp + i Kn n.om2 ]   (h_nn, n.om2 via projectors above)
  endpoint 2nd-order     +-EK [ Psi2 - dT2 + i Kn V2 + n.(v2T + om2) ]  (obs: +, emi: -)
  endpoint 1x1, observer clock/displacement, emission shift, emission lensing/time delay/aberration: see .tex.
CHECK: every symbolic channel evaluated numerically equals the engine channel to <= 1.8e-15 (`check_symbolic.py`, `.log`).
C3 type of each channel:  B = L1^2, -2<m1>m1 (products of full LOS maps: finite 3j sums), and endpoint 1x1 products (shell maps
at eta_* / eta_o: finite 3j sums with Bessel shells);  C = ray local 1x1 (single-time product) and all second-order fields
(Phi2/Psi2, dT2, V2, v2T+om2 [helicity +-1], h2 [helicity +-2]);  D = lensing, time delay, aberration+deflection, clock,
observer displacement, emission-surface shift (nested two-time or endpoint x functional).
### K -> 0 (mu -> -1) analytically (`k0_limits.py`, `k0_limits.log`)
Leading behaviour as K = |k1+k2| -> 0 (p^2 = q^2 - 2Kqc + K^2, c = Khat.qhat):  K_Q, K_SPhi, K_Phi2_0 ~ K^-2 (coefficient
∝ (3c^2-1) q^2: Poisson-gauge d_i d_j nabla^-4 S_ij), hence Phi2, Psi2 ~ K^-2;  K_om ~ K^-1 and K_vT ~ K^-1 with EXACTLY
opposite coefficients, K_om + K_vT ~ K^+1;  K_S0z ~ K (so V2 ~ K^-2, v2 = grad V2 ~ K^-1);  K_S00, K_h ~ K^0.
Power counting: a term f(K) P_d(n) e^{i K.n chi} with f ~ K^-p and P_d a degree-d polynomial in n contributes to multipole l
through j_L(K chi), |L - l| <= d, i.e. ~ K^{l-d-p}.  All terms have d + p <= 2 (Psi2: d=0,p=2; d_z Psi2: d=1,p=1; n.om2,
n.v2T, d_z V2: d=1,p=1; h_nn, d_z om_z: p=0), so K_{2,lm} = O(1) for l >= 2 (finite, Khat-dependent limit) while l = 1, 0
pieces are O(K^-1), O(K^-2).  With dmu = K dK/k0^2 the mu-integrand of C^(22)_l is bounded at mu = -1 for every l >= 2:
**integrable, no cancellation needed** (the dipole, l = 1, would be log-divergent channel by channel).  Numerical confirmation
below (integrand near K -> 0).  In the monochromatic integral we integrate in K (analytic) instead of mu.

## Milestone 4 — C_l^(22): structure and monochromatic evaluation (deliverable 4)
### Assembly (factors explicit)
  C_l ⊃ eps^4 (1/4) C^(22)_l,   C^(22)_l = (1/(2l+1)) sum_m <Th2_lm Th2*_lm>_c ,   Th2 = int int zeta zeta K_2 (K_2 symmetric)
  => Wick: <zeta1 zeta2 zeta1'* zeta2'*> has 2 connected pairings, each giving |K_2|^2:
  C^(22)_l = (2/(2l+1)) sum_m int d^3k1 d^3k2/(2pi)^6 P(k1)P(k2) |K_{2,lm}(k1,k2)|^2
           = 1/(4 pi^4 (2l+1)) int k1^2dk1 k2^2dk2 P P int_{-1}^{1} dmu sum_m |K~_lm(k1,k2,mu)|^2      (isotropy, Khat = z frame),
  K~_lm = sum_c K~^c_lm over channels; hence C^(22)_l = sum_{c,c'} G_cc'(l) (channel Gram matrix, cross terms included).
  Per C3 type (with C2's sign convention i^L instead of C3's (-i)^L):
   B x B (L1^2 - 2<m1>m1):  C^{BB}_l = C^{m1^2}_l - 4<<m1>^2> C^(11)_l,  C^{m1^2}_l = 2 sum_{l1l2} (2l1+1)(2l2+1)/(4pi) (l l1 l2;000)^2 C^(11)_l1 C^(11)_l2
                     [the -4<<m1>^2>C^(11) = -8 (cross with m1^2) + 4 (|<m1>m1|^2)];
   C: K~^C_lm = int d eta K(k1,k2,K,eta) sum_{|L-l|<=d} i^L sqrt(4pi(2L+1)) j_L(K chi) <lm|Q|L0>  (finite 3j sum, one mu-integral);
   D: bipolar double Rayleigh sums (two times) or, equivalently, the same frame/mu reduction with the nested time integrals
      kept inside K~_lm;  B x C, C x D cross terms: C3 (C2).
  Implementation: the full K_2 (all channels, verified vs N1 and between forms) is projected on Y_lm by sky quadrature in the
  Khat = z frame for each mu node, which realises the C and D reductions exactly and the B channels as a cross-check
  against the finite 3j sum (`typeB_check.py`).
### Monochromatic P = 2 pi^2 sigma^2 delta(k-k0)/k0^2
  C^(22)_l = (sigma^4/(2l+1)) int_{-1}^{1} dmu sum_m |K~_lm(k0,k0,mu)|^2 ,   K = k0 sqrt(2(1+mu)),  frame k1,2 = k0(+-sin a, 0, cos a), cos a = K/2k0;
  B part: C^{BB}_l = (4 pi sigma^2)^2 2 sum (2l1+1)(2l2+1)/(4pi)(l l1 l2;000)^2 Delta_l1^2 Delta_l2^2 - 16 pi sigma^4 Delta_0^2 Delta_l^2  (all at k0).
### Numerical evaluation (labelled EVALUATION: k0 = 1; point 1: eta_* = 2, eta_o = 8; point 2: eta_* = 2, eta_o = 4)
`cl22.py` (runs `cl22_<tag>.pkl/.log`), tables `analyze.py` -> `analyze_eo8.log`, `analyze_eo4.log`.  Units sigma^2 (C11), sigma^4 (C22).
```
 l | eta_o=8: C11/s^2      C22/s^4          (1/4)C22/s^4    | eta_o=4: C11/s^2      (1/4)C22/s^4
 2 | 9.546006061e-01  5.703017520e+01  1.425754380e+01 | 8.613350843e-03  3.993324543e+00
 3 | 6.250779591e-01  1.151253897e+01  2.878134742e+00 | 8.052093896e-03  6.124256484e-01
 4 | 1.750712875e-01  7.463675108e+00  1.865918777e+00 | 1.116078495e-03  5.287309571e-02
 5 | 1.933453283e-02  4.520572234e+00  1.130143058e+00 | 7.175773518e-05  7.780704277e-03
 6 | 1.770349639e-04  2.088249491e+00  5.220623728e-01 | 2.780479034e-06  1.672196841e-03
 7 | 3.174580556e-04  1.496197895e+00  3.740494738e-01 | 7.318262641e-08  2.619102919e-04
 8 | 2.304634771e-04  1.047895782e+00  2.619739455e-01 | 1.403749744e-09  2.798933064e-05
 9 | 6.045085961e-05  4.967381754e-01  1.241845438e-01 | 2.057008620e-11  2.143254526e-06
10 | 9.679528515e-06  1.603390998e-01  4.008477495e-02 | 2.382758584e-13  1.234364265e-07
11 | 1.105691361e-06  3.789202401e-02  9.473006002e-03 | 2.239446102e-15  5.557386369e-09
12 | 9.699082749e-08  8.010458273e-03  2.002614568e-03 | 1.743328105e-17  2.015547066e-10
13 | 6.823169755e-09  2.205827422e-03  5.514568556e-04 | 1.143109189e-19  6.028691151e-12
14 | 3.962943779e-10  7.797776236e-04  1.949444059e-04 | 6.402183756e-22  1.515307705e-13
15 | 1.940877947e-11  2.540234250e-04  6.350585624e-05 | 3.099033167e-24  3.249712209e-15  *
16 | 8.146618828e-13  6.723015045e-05  1.680753761e-05 | 1.309717407e-26  6.022021672e-17  *
17 | 2.968887381e-14  1.441356324e-05  3.603390811e-06 | 4.875242191e-29  9.746048238e-19  *
18 | 9.494896918e-16  2.557045655e-06  6.392614138e-07 | 1.610736317e-31  1.390051124e-20  *
19 | 2.688913155e-17  3.834413557e-07  9.586033894e-08 | 4.755727176e-34  1.758712758e-22  *
20 | 6.795223637e-19  4.946249047e-08  1.236562262e-08 | 1.262423388e-36  2.002188565e-24  *
```
Convergence: eta_o=8: (nK,nth) = (24,40), (40,48), (56,60) agree to <= 1.5e-11 relative for l <= 20 (<= 1e-13 for l <= 14).
eta_o=4: (24,40) vs (40,48) agree to <= 4e-11 for l <= 12, 1e-9 for l = 13,14; entries marked * (l >= 15) are below the
double-precision cancellation floor (~1e-16 of the largest channel, values < 1e-14) and are NOT reliable (rel. spread up to 2e-3).
Type-B cross-check: the B-channel diagonal G_BB from the mu-integral equals the finite 3j sum to 4e-12 (l<=12), 5e-9 (l<=20) at
eta_o = 8, and 3e-10 (l<=12) at eta_o = 4.
K -> 0 (mu -> -1): the total integrand sum_m |K_lm|^2 tends to finite values (eta_o=8, K = 0.0018: l=2: 1.335, l=4: 63.1, odd l
-> 0 as K^2 by the k1 <-> -k2 parity at K = 0); even l=0 is finite (213.5) in the total, i.e. the K^-2 channel pieces cancel in the
monopole as well — consistent with the analytic power counting.
Structure of the result: strong cancellations between channel groups (diagonal group contributions / C22 at l=2,4,8 for eta_o=8:
lensing 0.99, 3.6, 3.1; aberration+deflection 0.51, 1.85, 3.2; time delay 0.20, 0.22, 0.32; Phi2/Psi2 0.16, 0.30, 0.12;
B (L1^2-2<m1>m1) 0.025, 0.10, 0.016; cross terms -0.96, -5.1, -5.9): the transverse displacement (lensing) and direction
(aberration+deflection) remappings largely cancel, as expected for a remapping of a map that is itself generated on the ray.
Vector (om2+v2T) and observer clock/displacement diagonals are <1e-3 (the observer terms are l<=1 except through products).
Remark: in this monochromatic toy the first-order spectrum is band-limited (k0 chi_* = 6) whereas (22) reaches K <= 2k0, so
(1/4)C22/C11 grows rapidly with l (15 at l=2, ~1e3 at l=6-8 for eta_o=8).

## STATUS
DONE+VERIFIED: K_1/Delta_l/<m1> with observer terms; K_2 = all channels (engine + named-object sympy kernels, agreement 2e-15);
checks (iv) plane wave vs N1 (4 configs, within N1 errors), form equivalence, symmetry, reality (iii), constant mode (i:
consistency relation -k dK_1/dk, NOT zero — brief corrected), gradient (ii) vs V2 (PASS at FD accuracy); C^(22) assembly with
factors 1/4 and Wick 2, type-B 3j check, K->0 integrability analytic + numeric; monochromatic l=2..20 at two points.
FINDINGS for others: C2 ray.definitions k1/k2 sign (also found by A1/lead); clock shift = integral of Psi1 (lapse), closed form above.
NOT DONE: closed-form (non-numerical) expressions of the C-channel mu-integrand (the l-projection is done by sky quadrature,
equivalent to C3's finite 3j sums); per-channel LaTeX for the B/1x1 groups beyond symbolic_kernels.tex; A1 per-group cross-check.
