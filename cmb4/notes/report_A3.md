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

## Milestone 3 — collinear points mu -> +-1 of the monochromatic loop (DONE, numerical; `diag_mu.py`, logs diag_mu_[pm].log)
Monochromatic |q| = |k| = k0: at mu -> 1 the pair (k,-q) has |k-q| = k0 sqrt(2(1-mu)) -> 0 (and mu -> -1 likewise for (k,q)).
Every second-order Poisson-inverted object of that pair ~ 1/|k-q|^2 (e.g. C1's C3 -> 16/(9(1-mu)), Phi3c(eta_o) ~ -0.117/(1-mu),
V3c, delta3c similarly), so INDIVIDUAL channels of the mu-integrand diverge like 1/(1-mu) (log-divergent mu-integral):
at 1-mu = 1e-4 the channels 'time delay 2nd (radial Dx2)' (-2804), 'observer clock 2nd' (+2410), 'local fields:3' (channel iii;
-610 at l=3), 'emission shift 2nd' (+295), 'local 2+1' (+272) ... are O(10^3).  THE TOTAL IS FINITE: per-mu Delta^kappa_l
(l=2: -6.639, -6.294, -6.259 at 1-mu = 1e-2, 1e-3, 1e-4) converges linearly in (1-mu) (i.e. analytically in p^2), identically
at mu -> -1 (integrand even in mu).  Interpretation: the 1/p^2 parts of a long (p -> 0) second-order mode form, at leading
order, a K=0 homogeneous solution, i.e. a pure-gauge mode of the residual K=0 Poisson-gauge freedom (time reparametrisation),
which the full observable does not see (consistency relation for a long second-order mode); gradient pieces O(1/p) are
integrable anyway.  Hence the monochromatic C^(13) is finite and the mu-integral needs no regularisation; plain Gauss-Legendre.
sqrt3 resonances (|k-+q| = 2k0/sqrt3, mu = +-1/3): at finite eta all Green's integrals are entire in the frequencies
(cf. C4: Cin, Si entire), so the integrand is analytic there; nothing to do at finite eta_o.

## Milestone 4 — IR (two long modes, q -> 0)  (DONE numerically; `ir.py`, ir.log)
Brief's criterion "O(q^2 P(q))" is superseded (lead/A2): a constant zeta = c is the same universe in coordinates x -> e^c x,
so Th[c + zeta_s] = Th[zeta_s](e^c x) and, with Th3 = int zeta^3 K3, K3(0,0,k;n) = (k.d_k)^2 K1(k;n); the qhat-averaged
q-integrand kappa~(q) := 3 K3sym(k,q,-q) averaged = (1/2)<coef_abc>_qhat therefore tends to 3 (k d_k)^2 K1, i.e.
Delta^kappa~_l(k; q -> 0) = 3 (k d/dk)^2 Delta_l(k)   (O(q) gradient/SCT terms are odd in qhat and average out).
Numbers (k = 1, eta_* = 2, eta_o = 8; prediction 3(k d/dk)^2 Delta_l, l=1..4 = 25.748, -16.360, -10.235, 0.335):
  q = 0.10: 25.523, -15.785, -10.342, -0.218 ;   q = 0.05: 25.685, -16.210, -10.255, 0.197
  q = 0.025: 25.732, -16.322, -10.240, 0.301
  deviations shrink by 4.0 per halving (l=2: 0.575 -> 0.151 -> 0.038; l=4: 0.553 -> 0.138 -> 0.034):
  kappa~(q) = 3(k d_k)^2 K1 + O(q^2).  PASS.  This is exactly V2's prediction K_3(q,-q,k) = D^2 K_1(k) + O(q^2) (D = k.d_k,
  per unit c^2; kappa~ = 3 K3sym(k,q,-q) qhat-averaged), checked here for the FULL assembled kernel incl. normalisation.
Consequence: IR piece of 1/3 Re C13 = sigma^2_IR (2/pi) int k^2 dk P Delta_l (k d_k)^2 Delta_l, which with A2's (22) piece
sigma^2_IR (2/pi) int k^2 P (k d_k Delta_l)^2 forms (2/pi) int dlnk (k^3 P) (1/2)(d_lnk)^2 Delta_l^2 = (1/2)(2/pi) int dlnk Delta_l^2 d^2_lnk(k^3P)
(integration by parts), sigma^2_IR = int_IR d^3q/(2pi)^3 P(q).  For P ~ k^{n_s-4}: sigma^2_IR ~ int dq q^{n_s-2} diverges at q -> 0
for n_s <= 1, but its coefficient is ~ (n_s-1)^2 x (finite): exactly zero for n_s = 1; for n_s != 1 the total is IR-finite only
for n_s > 1 (n_s < 1: log/power IR sensitivity proportional to (n_s-1)^2, the usual dilation-induced IR sensitivity of a
non-scale-invariant spectrum).  O(q^2) remainder: int dq q^{n_s} -> converges for n_s > -1.  Monochromatic spectrum: no IR issue.

## Milestone 5 — tadpoles / zero-mode gauge (DONE; `gauge.py`, gauge_lin.log, gauge_cub.log)
The K=0 second-order modes (tadpoles <Phi2>, <Psi2>, <delta2>, C1's angle-averaged convention) enter kappa through: L3
'local 2+1' (tadpole x first order), radial k2 and Dx2 (time delay 2nd), the observer clock deta_o2 (lapse tadpole on the
worldline), the emission shift deta_e2 (dT2 tadpole), 3 L1 L2 and the normalisation -3<m2>m1 (<m2> tadpole), and channel
(iii) (C1's sources A1(k)<B2> and C3).  The tadpoles are NOT pure gauge as a whole (their eta-dependence ~ q^2 eta^2 is the
backreaction of the short modes on the homogeneous expansion), so their total coefficient is not expected to vanish; what must
vanish is the dependence on the residual K=0 Poisson-gauge freedom eta -> eta + eps^2 tau(eta)/2:
  d<Psi2> = tau' + tau/eta, d<Phi2> = -tau/eta, d<delta2> = -4 tau/eta (zeta_nl unchanged), and the exact transformation of the
  third-order fields dPhi3 = 3 tau T' + 6 (tau/eta) T, dPsi3 = 3 tau T' + 6(tau' + tau/eta) T, dV3 = 3 tau T_V' + 3 tau' T_V.
RESULTS: (1) C1's K3c sources with shifted tadpoles, solved with the same Green's function (+ the shift of the regular-mode
amplitude dC3 = dPhi3(0)), reproduce dPhi3c, dPsi3c, dV3c to 1e-13: C1's tadpole bookkeeping is exactly gauge covariant.
(2) Full assembled kappa (all tadpole routes incl. normalisation) for tau = 0.1 eta: unchanged to 6e-14 (scale 4e2) while
individual channels change by O(1-4) (time delay 2nd 3.6, clock 2nd 2.3, channel iii 1.35, local 2+1 1.18, radial k2 1.18,
emission 2nd 0.15).  tau ~ eta^3 (time-dependent gauge): see gauge_cub.log.
Note: in the shear form L2 contains NO tadpole field (the shear of a homogeneous perturbation vanishes), so 3 L1 L2 and
-3<m2>m1 involve the (q,-q) pair only through products of first-order objects (variance terms, channel i); the tadpole fields
act only at third order (routes listed above).  So check (a) holds in the precise sense: the gauge-dependent tadpole content
cancels identically among those routes; the remaining (gauge-invariant) tadpole dependence is physical backreaction and is kept.

## Milestone 6 — monochromatic C^(13) (in progress; `mono.py`, `combine.py`)
Definitions (explicit factors): Th3 = int zeta^3 K3 (K3 symmetric); polarisation coefficient coef_abc(k,q,-q;n) = 6 K3(k,q,-q;n);
K_3c(k;n) := int d^3q/(2pi)^3 P(q) K3(k,q,-q;n) = kappa/3, kappa = <dTh3/dzeta_k> = (1/2) int d^3q/(2pi)^3 P coef_abc;
kappa(nu) = sum_l (2l+1) i^l Delta^kappa_l P_l(nu) => Delta^(3c)_l = Delta^kappa_l / 3 and
   C13_l = 3 (2/pi) int k^2 dk P Delta_l Delta^(3c)*_l = (2/pi) int k^2 dk P Delta_l Delta^kappa*_l   (normalisation terms included in Th3).
Monochromatic: kappa(k0) = (1/2) int_{-1}^{1} dmu (1/2)<coef_abc>_phi  (|q| = k0), C13_l = 4 pi sigma^4 Delta_l Delta^kappa*_l
(= 12 pi sigma^2 Delta_l Delta^(3c)*_l of the brief).  <m2> enters via sky averages M2(k,q), M2(k,-q) (24x48 sky quadrature per mu)
and M2_tad = <m2>_(q,-q) = -7.46928 sigma^2-units; <m1> via Delta_0.  Im C13 ~ 1e-14 (reality OK).
Resolution lesson: the phi-integral needs Nphi >= 16 (phases e^{i q_perp.n_perp chi}, q chi_* = 6): at mu = 0.5,
Delta^kappa_2 = -0.970 (Nphi=8), -2.82625 (16), -2.82632 (32).  First runs A (16 mu x 12 phi) / B (10 mu x 8 phi) are
under-resolved in phi (30% spread at l=2) and are superseded by runs C (16 x 16, Nn=32) and D (12 x 20).

## Referee R1 points B1/B2 (addressed)
### B1 — mu -> +-1 behaviour of the full contracted kernel (monochromatic shell q = k)  [`residues.py`, residues.log]
Agreed with the referee's input: C3|_{q=k} = 16/(9(1-mu)) + O(1) (our tf.C3: 1.7759/(1-mu) at 1-mu=1e-4), and every object of the
long pair p = |k-q| -> 0 carries 1/p^2 = 1/(2k^2(1-mu)).  Channel residues R_c,l = lim_{mu->1} (1-mu) Delta^kappa_mu,c,l(mu)
(Delta^kappa_mu = Legendre coefficient of the per-mu integrand (1/2)<coef_abc>_phi; k = q = 1, eta_* = 2, eta_o = 8; from 1-mu = 1e-6):
| channel | l=0 | l=1 | l=2 | l=3 | l=4 | l=5 | l=6 |
|---|---|---|---|---|---|---|---|
| (iii) local fields:3 [v3 = grad V3c, incl. C3 ghat] | -0.34416 | +0.00803 | +0.01165 | -0.05832 | -0.07581 | -0.05898 | -0.03418 |
| (ii) time delay 2nd (radial Dx2 . grad I1) | +0.35486 | -0.05097 | -0.28054 | +0.05950 | +0.08132 | +0.05308 | +0.02098 |
| (ii) observer clock 2nd (deta_o2 I1|_o) | 0 | 0 | +0.23933 | 0 | 0 | 0 | 0 |
| (ii) emission-surface shift 2nd (deta_e2 I1|_e) | -0.01075 | +0.04295 | +0.02954 | -0.00119 | -0.00551 | +0.00590 | +0.01319 |
| (ii) local 2+1 | +0.00997 | -0.04053 | +0.02473 | +0.01970 | +0.00444 | -0.00261 | -0.00378 |
| (ii) radial k2 | -0.01001 | +0.04054 | -0.02473 | -0.01970 | -0.00444 | +0.00261 | +0.00379 |
| all other channels (lensing, aberration, transverse k2/Dx2, products, normalisation) | <= 4e-5 | | | | | | |
| **SUM** | -5e-5 | +3e-5 | -6e-6 | -2e-5 | -3e-6 | -1e-5 | +4e-7 |
The sum vanishes to the accuracy of the extraction (|R| ~ 1e-6 x finite part).  The six singular channels are EXACTLY the channels
that respond to the K=0 time-reparametrisation gauge mode in the gauge test (Milestone 5): the 1/p^2 part of the long second-order
mode is, at leading order, a homogeneous (K=0) solution, i.e. residual Poisson-gauge freedom, to which the observable is blind.
Hence (analytic leading behaviour): Delta^kappa_mu,l(mu) = F_l + G_l (1-mu) + ..., NO 1/(1-mu) and NO (1-mu)^{-1/2} term
(the O(1/p) gradient/SCT pieces are odd under the azimuthal average around k and drop); numerically F_l (l=0..6) =
-52.7137, +28.7453, -6.2562, -18.9101, -3.2841, -13.2461, +0.3736, and the approach is linear in (1-mu) (differences 1e-4 -> 1e-5
-> 1e-6 shrink x10).  mu -> -1 identical (integrand even in mu: the configuration (k,q,-q) is symmetric under q -> -q).
=> the monochromatic C^(13) is FINITE; no narrow-shell regularisation is needed (a narrow shell gives the same limit).
### B2 — tadpoles
The k->0 limit of K_Phi2 indeed does not exist pointwise (quadrupole q^2(3mu^2-1)/k^2); C1's angle average removes the quadrupole
(= the constant anisotropic shear zero mode, which averages to zero over qhat; in our assembly the tadpole multiplies first-order
objects of the free leg k only, so only the qhat-averaged tadpole can ever enter kappa: the quadrupole drops identically, not by
a convention).  The remaining isotropic tadpole has 3 functions (Phi, Psi, delta)(eta) - 2 K=0 constraints = 1 free function =
the time reparametrisation tau(eta) (+ the dilation constant fixed by zeta^(2)=0).  We showed the full kappa is invariant under
arbitrary tau (tau = 0.1 eta and tau = 0.1 eta^3: |d kappa| = 6e-14, 3e-13), C1's K3c covariant to 1e-13.  The particular
(sourced) part of the tadpole is not a gauge mode (source = qhat-averaged K=0 stress of the short modes) and is NOT expected to
have zero coefficient; it is the physical homogeneous backreaction and is kept.  (V5's constant-shear proof covers the quadrupole.)
Machine-generated channel listing: `code/A3_third_order/A3_channels.tex` (every L1, L2/2, L3/6 group with its expression, the
functionals, observer quantities and the polarisation contraction rule for channels i/ii/iii).

## Milestone 6 results — monochromatic numbers (DONE; runs C = 16 mu x 16 phi x 32 nu, D = 12 x 20 x 32; `combine.py C D`)
EVALUATION at k0 = 1, eta_* = 2, eta_o = 8, sigma = 1 (C11 in sigma^2, C22, C13 in sigma^4). C11 = 4 pi Delta_l^2 agrees with A2 to 3e-14.
C22 from A2 (cl22_eo8_b.pkl). Run C vs D: |rel diff| of 1/3 Re C13 <= 8e-5 (l <= 10), 3e-4 (l=12), 2e-3 (l=14); l >= 16 unreliable
(C13 there < 1e-9 and irrelevant next to C22).  Im C13 ~ 1e-14.
| l | C11 | 1/4 C22 | 1/3 Re C13 | 1/4 C22 + 1/3 Re C13 | ratio to C11 |
|---|---|---|---|---|---|
| 2 | 9.546006e-01 | 1.425754e+01 | -5.39801e+00 | +8.85953e+00 | +9.28 |
| 3 | 6.250780e-01 | 2.878135e+00 | -2.35152e+01 | -2.06370e+01 | -33.0 |
| 4 | 1.750713e-01 | 1.865919e+00 | -9.06394e+00 | -7.19802e+00 | -41.1 |
| 5 | 1.933453e-02 | 1.130143e+00 | -1.37218e+00 | -2.42041e-01 | -12.5 |
| 6 | 1.770350e-04 | 5.220624e-01 | -1.58643e-02 | +5.06198e-01 | 2.9e3 |
| 7 | 3.174581e-04 | 3.740495e-01 | -3.37230e-02 | +3.40326e-01 | 1.1e3 |
| 8 | 2.304635e-04 | 2.619739e-01 | -2.79364e-02 | +2.34038e-01 | 1.0e3 |
| 9 | 6.045086e-05 | 1.241845e-01 | -8.22586e-03 | +1.15959e-01 | 1.9e3 |
| 10 | 9.679529e-06 | 4.008477e-02 | -1.46307e-03 | +3.86217e-02 | 4.0e3 |
| 11 | 1.105691e-06 | 9.473006e-03 | -1.8408e-04 | +9.28892e-03 | 8.4e3 |
| 12 | 9.699083e-08 | 2.002615e-03 | -1.7666e-05 | +1.98495e-03 | 2.0e4 |
| 13 | 6.823170e-09 | 5.514569e-04 | -1.351e-06 | +5.50106e-04 | 8.1e4 |
| 14 | 3.962944e-10 | 1.949444e-04 | -8.49e-08 | +1.94860e-04 | 4.9e5 |
| 15-20 | | (A2 table) | < 5e-9 | = 1/4 C22 to 1e-4 | |
Delta^kappa_l (l=0..10): -96.795, +53.949, -4.6756, -25.171, -18.333, -8.3515, -1.0090, +1.6018, +1.5573, +0.89536, +0.39797.
Largest channel contributions to Delta^kappa_2 (run C): channel iii (v3) -22.47, k1-aberration/deflection x 2nd-order fields +12.54,
radial k2 -11.40, observer clock 2nd +11.14, time delay x 2nd +7.00, lensing x 2nd -6.76, transverse k2 +5.52, ... normalisation
-3<m2>m1 +3.96, 3L1L2 -4.00, 6<m1>^2 m1 small; strong cancellations (sum -4.68).  Interpretation caveat (as A2): the monochromatic
first-order spectrum is band-limited (k0 chi_* = 6) while the loop reaches |k+-q| <= 2k0, so for l >= 6 C11 is tiny and the
eps^4 terms dominate by orders of magnitude; at l = 2-5 Re C13/3 partially cancels 1/4 C22 (l=2: 38%; l=3,4 C13 dominates).
