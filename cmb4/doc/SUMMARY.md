# CMB temperature power spectrum through fourth order in a radiation-only universe — summary of results

Full manuscript: `cmb4/doc/main.tex` (compiled `main.pdf`). Code: `cmb4/code/`. Verification ledger: appendix of the manuscript.

## 1. The result that organises everything (exact, nonperturbative)
For a perfect radiation fluid (P = rho/3, rho ∝ T^4) in any spacetime, photons emitted isotropically at T_* from the surface
T_phys = T_* and observed by a comoving observer satisfy

    T_obs(n) = T_o * exp[ - S(n) ],     S(n) = ∫_e^o sigma_{mu nu} e^mu e^nu dl ,

i.e. the local fluid temperature times the exponential of minus the fluid SHEAR integrated along the ray (e = propagation
direction in the fluid frame, dl = proper length in that frame). Proof: the kinematic redshift identity
d ln E/dl = -(theta/3 + sigma_ee + a.e) combined with the two radiation-fluid identities theta = -3 u.d ln T and
a_mu = -h_mu^nu d_nu ln T, which turn theta/3 + a.e into -d ln T/dl along the ray. Checked in Kasner (1e-31), symbolically to third order, and nonperturbatively by an exact planar-symmetric GR solver (residual 1.4e-14 at
Theta_rms = 0.08).

Consequences: T_o drops out of the normalised map, 1 + Theta = e^{-S}/<e^{-S}>_sky; Sachs-Wolfe, ISW, Doppler, lensing,
time delay and the emission-surface location are all inside S (the emission surface enters only as the endpoint of the
integral); the observer's frame enters only through a rotation of the sky, which cannot affect C_l. For first-order scalar
modes the Fermi-Walker frame does not rotate relative to the canonical triad at O(eps) (vorticity, curl omega and Thomas
precession are O(eps^2)), so the rotation is rigid, R = 1 + eps^2 J_2/2, and C^(22) and Re C^(13) are even SEPARATELY the same
in both frames (only a generic O(eps) rotation would mix them). All orders are computed in the global canonical triad; a
per-direction triad is not a rigid rotation and must not be used at third order.

## 2. Bookkeeping ("the sweet spot")
* Normalisation is pure algebra: Lambda = -S = eps L1 + eps^2 L2/2 + eps^3 L3/6, m1 = L1, m2 = L2 + L1^2,
  m3 = L3 + 3 L1 L2 + L1^3; for l >= 1, Theta2 = m2 - 2<m1> m1, Theta3 = m3 - 3<m1> m2 - 3<m2> m1 + 6<m1>^2 m1
  (<.> = sky monopole, a random variable).
* Gaussian primordial variable: zeta_nl = lim_{eta->0}[ln(1+delta_T) + (1/2) ln(1-2Phi)] in Poisson gauge (= Maldacena's
  uniform-density zeta); "no higher-order primordial modes" = zeta_nl^(2)(0) = zeta_nl^(3)(0) = 0 and h_n(0) = 0; "regular"
  = analytic at eta = 0. Initial data at finite eta_i are the restriction of this solution (constraints automatic).
* Named objects: T_Phi(x) = -2(sin x - x cos x)/x^3 (x = k eta/sqrt3), T_delta, T_V; second-order kernels K_Q, K_SPhi,
  K_om (vector mode is ALGEBRAIC via relativistic Kelvin), K_h (= Kohri-Terada); Green's functions G_Phi, G_h; ALL second-order
  fields in closed form (elementary + Si + Cin, one special-function coefficient D(k,q,p)); IC kernel
  K_Phi2(0) (= Pitrou-Uzan-Bernardeau); third-order contracted kernels K3c and amplitude C_3; map functionals L1 (3 terms),
  L2 (34 terms), L3 (359 terms, 32 named groups); harmonic rules (types A-D).

## 3. First order (closed form)
Delta_l(k) = A(x_*) j_l(k chi_*) - k T_V(x_*) j_l'(k chi_*) + I_l(k) + [observer terms, l <= 1],
A = T_dT + T_Phi = -cos x + 2 j_0(x) - 4 j_1(x)/x, k T_V = sqrt3 (sin x - 2 j_1(x)), I_l = ISW integral in closed form for
every l (Si/Ci). C_l^(11) = (2/pi) ∫ k^2 dk P |Delta_l|^2. RD Sachs-Wolfe amplitude -zeta/3; scale-invariant plateau
l(l+1)C_l/(2 pi A_s) rises from 0.377 (l=2) to 0.914 (l -> inf) because the ISW is not negligible in RD.

## 4. Second and third order
* Second-order map: 15 channels (products incl. normalisation; second-order fields Phi2/Psi2, delta_T2, V2, omega2+v2T
  (helicity 1), h2 (helicity 2); nested lensing, time delay, aberration, emission shift, clock shift
  delta eta_o1 = -(1/eta_o) ∫ eta Psi1, observer displacement). C_l^(22) = (2/(2l+1)) sum_m ∫∫ P P |K_{2,lm}|^2.
* Squeezed limits (two independent routes agree): K_2(k1->0,k2) = -k.d_k K_1 (a constant zeta is a DILATION, not zero),
  O(k1) = special conformal relabelling; at third order K_3(q,-q,k) -> (k.d_k)^2 K_1. Infrared part of
  C22/4 + Re C13/3 for a power law = (1/2)(n_s-1)^2 sigma_L^2 C11: vanishes for n_s = 1, IR-sensitive for red tilt.
* Third order: kappa(k;n) = <dTheta3/d zeta_k> assembled from the self-contracted cubic products, first x second-order
  fields (tadpole and cross-contracted), and third-order fields; C13_l = (2/pi) ∫ k^2 dk P Delta_l Delta^kappa*_l.
* Zero modes (PROVEN): the only residual k = 0 freedom of Poisson gauge is the time reparametrisation eta -> eta + F(eta)
  (dilation, shear, rotation, translation are fixed by zeta^(2)(0) = 0, h(0) = 0 and the label). The divergent 1/k^2 (3mu^2-1)
  term of K_Phi2(0) is a quadrupolar/dipolar time-reparametrisation mode (not a constant shear); the q-angle average at fixed k
  removes all negative powers of k, so the tadpole limit exists and equals the angle-averaged values. The map is invariant to
  second order under every zero mode (rigid rotation + total derivative); at third order invariance is inferred from the same
  identity and confirmed numerically (assembled kappa unchanged to 3e-13 under tau(eta)).

## 5. Monochromatic spectrum, P = 2 pi^2 sigma^2 delta(k-k0)/k0^2 (variance <zeta^2> = sigma^2 held fixed)
C11 = 4 pi sigma^2 Delta_l(k0)^2 (closed form). C22 collapses to one integral over K = k0 sqrt(2(1+mu)); C13 to one
mu-integral. Collinear poles (PROVEN to cancel): at mu -> 1 the long second-order mode of two nearly opposite short modes has a
1/p^2 part that is exactly the K=0 time reparametrisation with (eta^2 tau)' = -eta^2 A_Q, A_Q = s^2[2T^2 + (T + eta T')^2]
(closed-form tau); the C_3 pole is the same shift. The six singular channel residues (clock, emission, time delay, v3, local 2+1,
radial k2) are explicit and cancel identically for every l by a total-derivative identity; numerically the residues agree to
9e-11 and sum to <= 8e-10. So C13 is finite, with no (1-mu)^(-1/2) term. Final result (boxed in the manuscript):
  C_l = eps^2 C11 + eps^4 [C22/4 + Re C13/3],  C11 = 4 pi sigma^2 Delta_l^2,
  C22/4 = sigma^4/(4(2l+1)) ∫ dmu sum_m |K~_{2,lm}(k0,k0,mu)|^2,  Re C13/3 = (4 pi/3) sigma^4 Delta_l Re Delta^kappa_l,
a function of (k0 eta_*, k0 eta_o) only at fixed sigma^2; evaluated at one point (no k0 scan). Numbers at k0 = 1, eta_* = 2,
eta_o = 8 (units sigma^2 / sigma^4):

    l :  C11       1/4 C22    1/3 Re C13    1/4 C22 + 1/3 Re C13
    2 :  0.9546    14.2575    -5.3980       +8.8595
    3 :  0.6251     2.8781   -23.5152      -20.6370
    4 :  0.1751     1.8659    -9.0639       -7.1980
    6 :  1.77e-4    0.5221    -1.59e-2      +0.5062
    10:  9.68e-6    4.01e-2   -1.46e-3      +3.86e-2

Positivity through this order requires eps^2 sigma^2 < 0.030 (l = 3), 0.024 (l = 4). Superhorizon check (k0 eta_o = 0.2, 0.1,
0.05): the O(eps^4) coefficient vanishes as (k0 eta_o)^(2l), like C11 (no uncancelled gauge piece); its ratio to C11 tends to a
finite constant (9.6 sigma^2 at l = 2, 46 sigma^2 at l = 3).

## 6. Verification (highlights)
Exact formula: algebraic, Kasner, planar GR solver. Einstein equations: Bianchi identity to eps^3; independent referee;
Kohri-Terada tensor source, Lu et al. vector constraint, PUB initial conditions reproduced. END-TO-END: the complete
analytic chain reproduces the exact planar solver's Theta1, Theta2, Theta3 maps for a single plane wave in four
configurations within the solver's noise (a_20[Theta3] to 7e-8), in a different gauge. Harmonic machinery reproduces the
full-sky lensed C_l (Hu 2000) exactly. C22: independent Monte Carlo (180,000 realisations of a 20-mode shell field, K_2 as a black
box) agrees with the quadrature within 1.1 sigma at 0.3-0.7% for l = 1..11 (l = 12: 2.0 sigma); it tests the assembly, Wick
factors and normalisation, not K_2 itself. The normalisation formulas (E2), (E3) reproduce the corresponding channels of kappa
to 1e-6 / 2e-5.

## 7. Unresolved / limitations (honest list, see manuscript Sec. "Unresolved steps")
UV behaviour of the loop for power-law spectra (probed numerically to q = 4k only; perfect fluid has no damping: expect a log
divergence at n_s = 1 like <zeta^2> itself); third-order Green's integrals and Delta^kappa_l left as explicit one-dimensional
integrals (nothing remains as an integral at second order); third-order zero-mode invariance inferred, not run separately;
the harmonic formula T^(ii) untested numerically; IR formula for tilted spectra derived, not tested; C13 only at one point; no general second-order gauge-transformation test
(only k = 0 family, residual tau(eta), and the different-gauge planar comparison); single-plane-wave test probes only
axisymmetric parts.
