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
* (ii) gradient O(k1): V2's report does not exist; the O(k1) residual is recorded (linear in s, coefficient 31.8 (l>=2 rms)
  at |k2| = 1.07), no independent prediction compared — NOT DONE.
