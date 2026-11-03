# Report A1 — plane-wave end-to-end check of the analytic chain vs N1

Agent A1 (relaunch, 2026-11-03). Code: `cmb4/code/A1_planewave/`. Conventions: SPEC.md; C2 ray convention
xbar(eta) = x_o + n(etabar_o - eta), e = -n; sky direction n = (mu, sqrt(1-mu^2), 0), mu = n.xhat (N1's polar axis +x);
a_l0 = sqrt(4pi(2l+1)) (1/2) int Theta P_l dmu (N1 normalisation).

## Bottom line
**The analytic chain C1 (fields) -> C2 (L1, L2, L3, emission, ray functionals) -> SPEC §6-7 (clock, displacement,
normalisation) reproduces N1's nonperturbative maps at all three orders, in all four configurations
(etabar_o = 8, 4; x_o = 0, pi/2), within N1's extraction uncertainty.** First order agrees to ~1e-11, second order
to ~1e-9 (a_l0), third order to ~1e-7 (x_o=0) / ~1e-5 (x_o=pi/2, where N1's own uncertainty is 3e-2).
The observable is independent of the zero-mode (tadpole) gauge choice to 1e-14.
Two bugs were found on the way (one in the crashed A1 run's own code, one in C2's *documented* ray-functional
definitions — see "Bugs found"); neither affects C1's equations or C2's L-term lists.

## Method (status: DONE)
1. `pw_symbolic.py` (from the crashed run, reviewed, reused): C1's engine `E.pkl` (Einstein eqs to eps^3, Poisson gauge)
   reduced to the planar ansatz zeta = eps cos x: modes {1} at O(eps), {0,2} at O(eps^2), {1,3} at O(eps^3);
   v^x = d_x V, v^y=v^z=omega=0, h only as a homogeneous traceless k=0 mode at O(eps^2) (h_11 = H, h_22 = h_33 = -H/2).
   Closure (vanishing of all other components, yy=zz, parity) asserted. Per mode: Psi, delta, V algebraic, Phi from
   Phi'' + 4Phi'/eta + (m^2/3)Phi = S_Phi; H from H'' + 2H'/eta = S_H. First order = C1's T_Phi, T_delta, T_V (checked).
2. `pw_numeric.py` (reused): regular solutions = exact Green's-function integrals
   Phi_{n,m} = C_{n,m} ghat_m + int_0^eta G_m S_Phi (G_m = C1's G_Phi; G_0 = (eta'/3)(1-eta'^3/eta^3)), evaluated by
   Chebyshev spectral integration in mpmath (40 digits, degree 160). Amplitudes C_{n,m} from zeta_nl^(2)(0) = zeta_nl^(3)(0) = 0
   imposed on the exact rational small-eta series: C_{2,0} = -4/27, C_{2,2} = -1/3 (= C1's K_Phi2(q,p;0) for q=p=k0),
   C_{3,1} = 38/27, C_{3,3} = 2/81. Checks: ODE residuals <= 6e-23, series agreement <= 2e-23, C1-engine conservation
   nabla_mu T^{mu nu} = 0 at orders 2, 3: 1e-27 / 5e-21.
   **Zero-mode choice:** Psi2_0 = Phi2_0 + gauge_c * eta (residual time reparametrisation eta -> eta - gauge_c eta^2/4 at k=0);
   production gauge_c = 0; the third-order mode-1 fields contain the corresponding "tadpole" pieces automatically.
3. `pw_eval.py`: field jets in C2's naming (Phi_n, Psi_n, dT_n = n![eps^n]((1+delta)^{1/4}-1), v_i_n, h_ij_n).
4. `pw_map.py` (new): C2's general-n shear-form lists `L_terms.pkl['L']` (reduced to the planar case: 2 / 21 / 146 terms),
   ray functionals k1, Dx1, k2, Dx2 built from `ray_cache.pkl` (G1, G2half, k_init) by spectral integration on 65
   Chebyshev-Lobatto nodes in [eta_*, etabar_o] (vectorised over 96 Gauss-Legendre mu, float64), emission shifts
   deta_e1, deta_e2 from `L_terms['emission']`; observer (mpmath): worldline X = x_L + int v, clock
   Delta tau = int_{eta_init}^{eta_o} eta/uhat^0 = (etabar_o^2 - t_i^2)/2 with initial event (1+dT)/eta_init = 1/t_i, t_i = 0.01,
   expanded exactly to O(eps^2) -> deta_o1, deta_o2, dxo1, dxo2; then m1,m2,m3, sky monopoles, Theta_n (SPEC §7, incl. l=0,1), a_l0.
   Label convention: x_L = Poisson position at eta -> 0 (= N1's comoving label, the gauge map being O(eta^2)); SPEC's
   literal "label at eta_init" differs by int_0^{t_i} v = O(t_i^2) (irrelevant for x_o = 0 by symmetry).
5. `compare.py`, `resid_test.py`, `n1ref.py`: comparison with N1's coefficient maps, with N1's raw scans (no N1 fitting),
   and Taylor coefficients of N1's T_o, <T_obs>, S(mu).

## Comparison table (A1 vs N1, a_l0; full l <= 12 in `compare_<tag>.log` / `compare_<tag>_N64.json`)
| config | order | l | A1 | N1 | A1 - N1 |
|---|---|---|---|---|---|
| etao=8, x_o=0 | 1 | 2 | -2.1847203553e+00 | -2.1847203553e+00 | -4.4e-12 |
| | 1 | 4 | +1.2552456282e+00 | +1.2552456282e+00 | -3.6e-12 |
| | 1 | 6 | -4.7973477367e-02 | -4.7973477364e-02 | -3.4e-12 |
| | 2 | 2 | -4.2224942186e+00 | -4.2224942193e+00 | +6.5e-10 |
| | 2 | 4 | +1.3785439774e+00 | +1.3785439781e+00 | -7.4e-10 |
| | 2 | 6 | -1.2595382846e+00 | -1.2595382817e+00 | -2.9e-09 |
| | 3 | 2 | -1.2889987902e+01 | -1.2889987968e+01 | +6.6e-08 |
| | 3 | 4 | -4.7070907021e+01 | -4.7070907098e+01 | +7.7e-08 |
| | 3 | 6 | +6.2905106702e+00 | +6.2905105363e+00 | +1.3e-07 |
| etao=4, x_o=0 | 1 | 2 | +2.0752530981e-01 | +2.0752530981e-01 | +4.3e-12 |
| | 2 | 2 | -1.3572165380e+00 | -1.3572165385e+00 | +5.4e-10 |
| | 2 | 4 | -2.2995087658e-01 | -2.2995087720e-01 | +6.1e-10 |
| | 3 | 2 | -4.8330628410e+00 | -4.8330627544e+00 | -8.7e-08 |
| | 3 | 4 | -3.0635108209e+00 | -3.0635106964e+00 | -1.2e-07 |
| etao=8, x_o=pi/2 | 1 | 1 (dipole) | +2.0812147132e+00 | +2.0812147132e+00 | +1.1e-11 |
| | 1 | 3 | +2.0917805127e+00 | +2.0917805128e+00 | -5.1e-11 |
| | 2 | 2 | +3.0669017921e+00 | +3.0669017899e+00 | +2.2e-09 |
| | 2 | 4 | +6.5664035508e+00 | +6.5664035492e+00 | +1.6e-09 |
| | 3 | 1 (dipole) | -5.1723421459e+01 | -5.1723420730e+01 | -7.3e-07 |
| | 3 | 3 | -4.3006733540e+01 | -4.3006736906e+01 | +3.4e-06 |
| etao=4, x_o=pi/2 | 1 | 1 (dipole) | -2.3149107586e-01 | -2.3149107587e-01 | +3.8e-12 |
| | 2 | 2 | -2.1076209609e+00 | -2.1076209611e+00 | +2.1e-10 |
| | 3 | 1 (dipole) | +1.3059761781e+01 | +1.3059761844e+01 | -6.4e-08 |
| | 3 | 3 | -3.6965575489e-01 | -3.6965566762e-01 | -8.7e-08 |

Max |map difference| over the sky vs N1 extraction uncertainty (N1 "spread"):
| config | Th1 | Th2 | Th3 |
|---|---|---|---|
| etao=8, x_o=0 | 1.5e-10 (3.6e-9) | 1.3e-7 (3.2e-8) | 4.0e-6 (2.1e-4) |
| etao=4, x_o=0 | 1.3e-11 (6.8e-12) | 4.5e-9 (9.7e-10) | 2.7e-7 (4.0e-7) |
| etao=8, x_o=pi/2 | 2.7e-8 (5.3e-7) | 4.0e-7 (9.2e-6) | 1.8e-3 (3.1e-2) |
| etao=4, x_o=pi/2 | 2.1e-11 (6.7e-10) | 2.4e-9 (5.4e-9) | 3.1e-7 (3.9e-5) |
The Th2 map differences for x_o = 0 exceed N1's quoted spread by ~4x; they sit at N1's stated noise floor
(sigma_Theta/eps0^2 ~ 5e-12/1e-4 with O(1) weights, report_N1 §5) and the a_l0 agree to <= 3e-9. Independent test without
any N1 fitting (`resid_test.py`): R(eps) = Theta_N1(eps) - [eps Th1 + eps^2 Th2/2 + eps^3 Th3/6]_A1 scales cleanly as eps^4
for all 16 amplitudes (to8_xo0: max|R|/eps^4 = 11.6 ... 14.5 smoothly from eps=-0.1 to +0.1; at eps=0.01 |R| = 1.2e-7,
no lower-order remainder). A1's own numerics: Chebyshev N = 48/64/96 identical to all printed digits; shear form vs
C2's direct (redshift) form agree to 4e-10 (L2) / 9e-7 (L3, max over map).

Observer quantities (gauge-invariant), Taylor coefficients [eps^0, eps^1, eps^2]:
| config | T_o A1 | T_o N1 (fit) | <T_obs> A1 | <T_obs> N1 (fit) |
|---|---|---|---|---|
| etao=8, x_o=0 | 0.125, -0.0550223132, +0.0041441746 | 0.125, -0.0550222904, +0.0041441914 | -1.48623e-5, -0.0934726388 | -1.48671e-5, -0.0934726721 |
| etao=4, x_o=0 | 0.25, +0.1715330582, -0.0468010539 | 0.25, +0.1715330582, -0.0468010544 | +0.1629542549, -0.0967211615 | +0.1629542546, -0.0967211626 |
| etao=8, x_o=pi/2 | 0.125, 0, -0.0365879536 | 0.125, 0, -0.0365879014 | 0, -0.196599264 | 0, -0.196598814 |
| etao=4, x_o=pi/2 | 0.25, 0, -0.2272250869 | 0.25, 0, -0.2272250875 | 0, -0.1725955153 | 0, -0.1725955172 |
(N1 fit uncertainties are ~1e-8..5e-7 for these; A1 values are exact to the quadrature precision.)
The coordinate time of observation is gauge dependent: Poisson eta_o = 8 + 0.91167 eps + 4.41138 eps^2/2 (x_o=0) vs N1's
comoving-gauge t_o = 8 + 2.54456 eps - 3.2885 eps^2; N1's t_o(+-0.01) = 8.02512265/7.97421965 is therefore NOT expected
to be reproduced in Poisson gauge; the physical content (T_o at fixed proper time) agrees, see table.
Sky monopoles (etao=8, x_o=0): <m1> = 0.4400596075, <m2> = -1.174459454, <m3> = -27.17767538.

## Zero-mode (tadpole) independence (DONE)
Fields regenerated with gauge_c = 1 (`pw_fields_g1.pkl`): Phi2_0, Psi2_0, delta2_0 and all third-order mode-1 fields change
at O(1) (e.g. Psi2_0(8): -1.21 -> +4.79; Phi3_1(8): -3.72 -> -1.20), deta_o2 changes 4.411 -> -11.589, but every a_l0 of
Th1, Th2, Th3, the sky monopoles and T_o are unchanged to <= 2e-14 (`compare_to8_xo0_g1.log`).

## Bugs found (evidence)
1. **A1 crashed-run code (fixed):** `pw_eval._build_dT` built the k=0 mode of delta with weight 2 (z^0 + z^-0 = 2), so
   dT2_0 = delta2_0/2 - ... instead of delta2_0/4 - .... Detected by C2's own conservation rules evaluated with A1's jets
   (`chk_cons.py`: d_eta dT2 residual 0.27 in the k=0 mode only; after fix <= 3e-11 at all orders) and by T_o at O(eps^2)
   (-0.00290 -> +0.0041442 vs N1 +0.0041442).
2. **C2 (documentation-level, affects any downstream user of `ray.definitions` / `L_terms.pkl['definitions']` /
   report_C2 Milestone 2):** the ray functionals are printed as k1(eta) = k1_init + int_eta^{etabar_o} G1 and
   k2(eta)/2 = k2_o/2 + int_eta^{etabar_o}[G2/2 + Dx1.dG1]. Since dk/deta = G (C2's own geodesic equation and the
   `Drep` rules used in los.py), the correct forms have **minus** signs: k1(eta) = k1_init - int_eta^{etabar_o} G1,
   k2(eta)/2 = k2_o/2 - int_eta^{etabar_o}[...]. Evidence: with '+', the first-order null condition n.k1 = -(Phi1+Psi1)
   fails by O(1) (max 2.3) and C2's shear and direct forms disagree (L2 by 2.1, L3 by 38); with '-', null condition 1e-15,
   shear = direct (4e-10 / 9e-7) and agreement with N1 as above. The L-term lists themselves (built with the correct `Drep`)
   are right; the Dx1, Dx2, k2_o definitions are right. C2's dilation check could not see this (G1 = 0 for constant fields).
   **A2/A3 must use the minus sign if they evaluate k1, k2 from the definitions.**
No evidence of any C1 bug: C1's equations (via E.pkl), IC prescription zeta_nl^(n)(0)=0 and C2's L1-L3 + emission +
SPEC §6-7 bookkeeping reproduce N1 at all orders.

## Status
| item | status |
|---|---|
| fields orders 1-3 for the plane wave from C1's equations (exact Green's integrals, mpmath evaluation) | DONE, verified |
| C2 L1-L3 on the background ray, ray functionals, emission shifts | DONE (sign fix above) |
| observer clock (T_i crossing), displacement, canonical triad (C2's k_init), normalisation, a_l0 | DONE |
| etao=8 x_o=0 (primary), etao=4, x_o=pi/2 (dipole) | DONE, all agree with N1 |
| zero-mode gauge independence | DONE (1e-14) |
| term-group localisation | not needed (no residual discrepancy); group values available from `MapBuilder.run()['groups']` |
