# Report V3 — hostile referee of C1 (second/third order), relaunch 2026-11-03

Code: `cmb4/code/V3_referee/` (logs of every run in `logs/`). C1's files were not modified (read and imported only).
Arbiter: sympy exact, or 30-digit numerics at random rational (k, q, mu, eta).

## 0. Independent engine (basis of attacks A, B, D)
- `indep_einstein.py`: a separate plain-sympy computation of E_{mu nu} = G - 8piG T for g = eta^2 ghat (SPEC §4),
  coded independently of C1's `tpoly` engine. `compare_engine.py`: **term-by-term identical to C1's E.pkl for all 10
  components at orders 1, 2, 3 and for u_mu**.
- `check_dilation.py`: the exact dilation family (Psi=(A^4-1)/2, Phi=(1-A^2B^2)/2, delta=A^-4-1) solves my E to O(eps^3)
  for arbitrary A(eps), B(eps); also checked non-perturbatively (G - 8piG T = 0 exactly), zeta_nl = ln B.
- Caveat: the brief asked for a *different method* (ADM / exponential parametrisation). What I have is the same method
  coded independently. To compensate, each sector was checked against independent literature (PUB App. A scalar
  sources, Kohri–Terada tensor source, Lu et al. vector constraint), see below.

## Verdicts

### Claim (1) linear — CONFIRMED
`check_first.py`: Psi1 = Phi1; v_x = v_y = 0; V1 = -(eta^2/2)(Phi1' + Phi1/eta); delta1 = -(2/3)eta^2k^2Phi1 - 2eta Phi1' - 2Phi1;
Phi'' + 4Phi'/eta + k^2Phi/3 = 0, solved by T_Phi = -2(sin x - x cos x)/x^3, T_Phi(0) = -2/3; lim(delta1/4 - Phi1) = 1;
all ten E^(1) vanish on the closed-form solution.

### Claim (2) operator structure + tensor kernel — CONFIRMED
- `check_ops.py` (TT/transversality imposed, k || z), at **both orders 2 and 3** from my E: Q_n = S_TL/k^2;
  delta_n = eta^2 S00/3 - (2/3)eta^2k^2Phi_n - 2eta Phi_n' - 2Psi_n; V_n = -(eta^2/4)[2(Phi_n'+Psi_n/eta) + S0z/(ik)];
  0x: k^2 om/2 + 4(om+v^T)/eta^2 + S0x = 0; xz: -(ik/2)(om'+2om/eta) + Sxz = 0; h'' + 2h'/eta + k^2h = -2(Sxx-Syy);
  trace => Phi'' + 4Phi'/eta + k^2Phi/3 = -(S_tr - S00)/6 + Q'/eta - k^2Q/3. All exact. (The "False" lines in
  `check_second.py` part 1 came only from not imposing TT/transversality; superseded by check_ops.)
- `check_second.py` part 2: K_S00, K_S0z, K_Q, K_SPhi, K_h from discrete plane waves in my E vs C1's closed forms,
  4 random points: agreement to 1e-164; S0y = Syz = Sxy = 0.
- `check_KT.py`: C1's K_h vs Kohri–Terada (20) f_RD transcribed from the paper; mapping h^KT = h2 (KT g_ij ⊃ h/2, ours
  ⊃ eps^2 h2/2), w = 1/3: **K_h = 4 f_RD(q/k, p/k, k eta)**, relative diff 5e-29 at 6 random points.
- `check_green.py`: G_Phi, G_h satisfy L[int G S] = S (Wronskian sign right); ghat -> 1; S_Phi^(2) finite at eta -> 0
  (only eta^0, eta^2), K_SPhi(0) mine = C1.
- PUB (A.1)–(A.8) from raw text (`check_pub2.py`): S1 = -K_S00 (00), S3 = 2K_Q (traceless), S4 = -K_S0z/(ik) (0i, after
  identifying PUB's velocity as v^I = a u^I = uhat^0 v^i_ours, i.e. V2_PUB = V2 - 2 d^-1(Psi1 dV1)): exact (15 digits) at 3
  points. (The crashed run's `check_pub.py` had a sign error, (DD)^-1 d_I d_J -> -k_I k_J/k^4; fixed in check_pub2.)
  PUB's trace source S2 (A.6) does NOT match: `fit_pub_S2.py` finds exactly
  S2_PUB + S_tr/6 = -2Phi^2/eta^2 - 4Phi Phi'/eta + Phi'^2/2 + (4/3)Phi lap Phi + (4/3)(dPhi)^2 + (eta/3)dPhi'.dPhi + (eta^2/6)(dPhi')^2.
  Printed (A.6) is visibly garbled (8H(Phi-Psi)Psi' appears twice); our S_tr is fixed by two engines that agree term by
  term and satisfy the exact dilation family. I attribute this to PUB (A.6) as printed, not to C1 (open literature note).

### Claim (3) omega2 algebraic — CONFIRMED
- `check_vector.py`: my 0x constraint + relativistic Kelvin curl[(1+delta_T)uhat_i] = 0 give omega, v^T identical to C1's
  K_om, K_vT (|diff| ~ 1e-164, 4 points); this omega satisfies my ij-vector evolution equation identically.
- Kelvin-free argument: K_om = O(eta) as eta -> 0 while the homogeneous solution of om' + 2om/eta = src is eta^-2, so the
  algebraic omega is the unique regular solution of the evolution equation.
- Printed K_om = (2i/k^2)[eta(q^2-p^2)T_qT_p + eta^2(q^2T_qT'_p - p^2T_pT'_q)] equals the pickled kernel.
- Lu–Ananda–Clarkson–Maartens 0812.1349 eq. (20) (RD, w = 1/3), with ds^2 ⊃ -a^2 S_i dx^i deta (ghat_0i = -S/2), S of X2
  type: omega2 = -S has kernel **identical** to C1's K_om (symbolic ratio 1).

### Claim (4) zeta^(2), Phi2(k,0) — CONFIRMED
- zeta^(2) = delta2/4 - delta1^2/4 - Phi2 - 2Phi1^2 (`check_ic2.py`, `check_zeta3.py`; equal to SPEC's ln(1+delta_T) form).
- Independent eta -> 0 limit of my 00 constraint and traceless source (Taylor-expanded transfer functions):
  lim eta^2 S00 = 32/3, Q2(0) = -(4/9)(4k^2 - 6kq mu + 9mu^2q^2 - 3q^2)/k^2; zeta^(2)(0) = 0 gives
  Phi2(k,0) = -(4/9)(2k^2 - 2kq mu + 3mu^2q^2 - q^2)/k^2 = -(7k^4 - 2k^2(p^2+q^2) + 3(p^2-q^2)^2)/(9k^4); difference from C1: 0.
- PUB (3.12b,c) from raw text (paper.txt l.457–462), R^(1) = -zeta, R^(2) = -2R^(1)2 (Gaussian R_M, PUB (3.19)), PUB Psi =
  our Phi: (3.12b) - mine = 0; (3.12c) vs -Q2(0): 0.
- Zero mode: angle-averaged k -> 0 limits Phi2 -> -8/9, Q2 -> -16/9 reproduced; they are the dilation-family member
  A = e^{-c/3}, B = e^c (Psi2 = 8/9, Phi2 = -8/9, delta2 = 16/9), satisfy the k=0 00 constraint and zeta^(2) = 0.

### Claim (5) third order — CONFIRMED (with a convention caveat on the zero mode, and a harmless erratum)
zeta^(3) expansion exact; contraction rules confirmed by brute-force Wick; contracted sources K3c match my E^(3) exactly;
C3 re-derived fully independently (my E^(3), my small-eta series, my tadpoles, my zeta expansion) matches C1's C3 at 3 random
points to 30 digits; squeezed limit -32/9 re-derived from the separate-universe/dilation argument. Caveat: C3 depends on the
(gauge) zero-mode convention for the second-order tadpoles; C1's angle-averaged choice is consistent, but downstream assembly
must use the same tadpoles together with the observer clock/normalisation (SPEC 9d). Details §5.

## 5. Third order
- `check_C3_indep.py` (new, the main independent test): discrete three-mode realisation in MY E^(3), every first/second-order
  field replaced by its exact small-eta Taylor series (Phi2 = C(0) + eta^2[S_Phi(0)/10 - c^2 C(0)/30], h2 = (8/9)eta^2 Lambda qq,
  omega/V/delta from my constraints), tadpoles computed HERE as angle-averaged k -> 0 limits (tadQ2 = -16/9 + 88q^2eta^2/405,
  tadPhi2 = -8/9 - 28q^2eta^2/405, taddelta2 = 16/9 - 472q^2eta^2/405), zeta_nl expanded HERE; zeta^(3)(eta->0) = 0 solved for
  the image-summed amplitude.  Result: C3sym(mine) = C3(k,q,mu) + C3(k,q,-mu) of C1 **exactly (|diff| = 0 at 30 digits) at
  3 random points**; unchanged with C1's (wrong) 3/10 coefficient (log C3i_a310).
- Zero-mode convention test (logs C3i_shiftPhi/Psi): shifting the tadpole <Phi2> by +1 shifts C3sym by exactly +8 (4 per
  image); shifting <Psi2> by +1 (with delta2 -> delta2 - 2, as the k=0 00 constraint requires) changes nothing.  The zero-mode
  equations leave (Phi2, Psi2) free (Poisson-gauge residual dilations); zeta^(2)_{k=0} = 0 forces Phi2 = -Psi2/2 - 4/9 (C1's
  -8/9, 8/9 satisfies it).  Along this allowed one-parameter family C3 changes by -2 dPsi2 per image, so **C3 is a
  gauge-convention-dependent intermediate**, not unique.  C1's choice (continuity in k, = dilation member A = e^{-c/3},
  B = e^c) is internally consistent; the final C_l is only correct if A3/assembly uses the same tadpoles together with the
  matching observer clock and normalisation (SPEC 9d).  Not an error, but must be carried forward explicitly.
- Squeezed limit, independent argument: a constant long mode c in the dilation family A = e^{-c/3}, B = e^c gives, at eta -> 0,
  Phi_total = Phi_bg(c) + A^2B^2 Phi_s = Phi_bg + e^{4c/3}(-2/3)zeta_s with zeta_total = c + zeta_s(x/B) (linear => Gaussian);
  the c^2 zeta_s coefficient is (8/9)(-2/3), times 3! = -32/9 per P(q). C1's C3 at q = 0 is -16*2k^5/(9k^5) = -32/9,
  mu-independent. Since C1's rational C3 agrees with my independent value at random points, the limit is confirmed.
- `check_zeta3.py`: zeta^(3) = delta3/4 - (3/4)delta1 delta2 + delta1^3/2 - Phi3 - 6Phi1Phi2 - 8Phi1^3: exact.
- `check_wick.py`: brute-force lattice Wick contraction confirms C1's cubic (3-term) and bilinear
  (2A(q)B2(-q,k) + A(k)<B2>) contraction rules.
- `check_third.py`: C1's contracted sources K3c (S00, S0z, Q, Sxx-Syy, S0x, Sxz) vs my E^(3) on an explicit three-mode
  realisation (zeta1 e^{iq.x}, zeta2 e^{-iq.x}, zeta3 e^{ik.x}) with explicit second-order pair fields: the identity
  6[z1z2z3]E^(3) = K3c(mu) + K3c(-mu; p -> |k+q|) holds exactly at 3 points (worst diff 0). (`debug_third.py` from the crashed
  run was chasing the image s -> -s of transverse components; resolved, no discrepancy.)
- ERRATUM (minor, harmless) in C1 `third_order_ic.py` l.72: KPhi2_small uses Phi2 = C + (3/10)S_Phi(0)eta^2 - ...; the regular
  particular solution of Phi'' + 4Phi'/eta = S0 is S0 eta^2/10. (C1's own tadpole -8/9 - 28q^2eta^2/405 does use 1/10.)
  `check_C3.py`, `check_SPhi3_finite.py`: C3, Q3c(0), lim eta^2 S00c, the zeta^(3) contracted limit and the absence of
  eta^-4..eta^-1 terms in K3c[SPhi] are all independent of this coefficient; only the O(1) value K3c[SPhi](eta=0) (feeding the
  O(eta^2) term of Phi3c) changes.
