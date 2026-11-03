# Report C5 — first-order transfer function Delta_l(k), C_l^(11), sky monopole

Code: `cmb4/code/C5_linear_Cl/linear_cl.py` (run `python3 linear_cl.py`; log `test_log.txt`).
Conventions: SPEC.md; C2's ray convention x(eta) = x_o + n (eta_o - eta), n = observation direction, observer x_o = 0,
chi = eta_o - eta, x = k eta/sqrt3, s_* = k eta_*, S = k eta_o.

## Milestone 0 — reconstruction of the crashed run (status: code + full test suite PASS, written up here)
The previous run (2026-11-01) produced `linear_cl.py`, `run.log`, `test_log.txt` (identical; ALL PASS, 64 s) but no
report.  Re-read and re-run; content:

### Linear map and transfer functions (from C1 §2, re-verified)
Theta1(n) = [dT1 + Psi1 - n.v1]_e - [same]_o + int_{eta_*}^{eta_o}(Phi1'+Psi1') d eta  (SPEC §9a).
Phi1 = Psi1 = zeta T_Phi, T_Phi = -2 j_1(x)/x;  dT1 = zeta T_dT, T_dT = -cos x + 2 j_0(x) + T_Phi;
v1 = grad V1, k T_V = sqrt3 (sin x - 2 j_1(x));  2 d_eta T_Phi = (4k/sqrt3) j_2(x)/x.
Checked against C1's pickles and the first-order ODE, 00/0i constraints, continuity: residual 0.
Emission source A := T_dT + T_Phi -> -1/3 (RD Sachs-Wolfe), kT_V/x -> 1/sqrt3.

### Rayleigh rules (C2 convention)
e^{i k chi n.khat} = 4pi sum_lm i^l j_l(k chi) Y_lm(n) Y*_lm(khat);  n.khat e^{...} -> 4pi i^l (-i j_l') Y*_lm, i.e.
(i n.khat)^s -> (d/dz)^s j_l(z).  Observer terms carry no phase (chi = 0): the observer monopole is pure l = 0
(j_l(0) = delta_l0) and the observer Doppler n.v_o = i k mu T_V(x_o) zeta is pure l = 1 (j_l'(0) = delta_l1/3).
Verified numerically (Gauss-Legendre sphere quadrature, l <= 6): max error 2.9e-14.

### Result
a_lm = 4pi i^l int d^3k/(2pi)^3 zeta_k Y*_lm(khat) Delta_l(k),
  Delta_l(k) = G_l(eta_*) - G_l(eta_o) + I_l(k),
  G_l(eta) = A(x) j_l(k chi) - kT_V(x) j_l'(k chi),
  G_l(eta_o) = A(x_o) delta_l0 - kT_V(x_o) delta_l1 / 3,
  I_l(k) = int_{eta_*}^{eta_o} 2 T_Phi' j_l(k(eta_o-eta)) d eta = (4/sqrt3) int_{s_*}^{S} [j_2(s/sqrt3)/(s/sqrt3)] j_l(S-s) ds.
C_l^(11) = (2/pi) int k^2 dk P(k) Delta_l(k)^2 (from <zeta zeta>, Y addition theorem, i^l(-i)^l = 1, Delta_l real).
Monochromatic P = 2pi^2 sigma^2 delta(k-k0)/k0^2: C_l^(11) = 4pi sigma^2 Delta_l(k0)^2.
Prefactor, observer terms and ISW verified by a REAL-SPACE check: xi(n1.n2) = sigma^2 <K K*>_khat computed directly
from the kernel vs sum_l (2l+1)/(4pi) C_l P_l: agreement 6e-15 at 5 angles.

### ISW closed form (general l algorithm, explicit l = 0..6)
j_l(u) = C_l(u) sin u + D_l(u) cos u (Laurent polynomials); j_2(x)/x likewise; product-to-sum gives phases
cos/sin(a_± s - S), a_± = 1 ± 1/sqrt3 (never resonant); partial fractions in 1/(s^p (S-s)^q); antiderivatives of
cos(a t)/t^m, sin(a t)/t^m by recursion down to Ci, Si.  The observer endpoint u = S - s -> 0 is taken as a finite
part; all poles and log u cancel (checked symbolically/numerically).  Result: I_l = rational x trig in (s_*, S) +
Ci/Si at the six arguments a_± s_*, a_± S, a_± (S - s_*).  Term counts l = 0..6: 106, 163, 228, 307, 400, 495, 598.
Closed form vs mpmath quadrature (5 parameter sets): max rel. err 2e-18.  eta_* -> 0 closed forms also checked.

### Limits (as in the log)
- small k: M(k) = Delta_0 = -k^4 (eta_o-eta_*)^3 (eta_o+3 eta_*)/810 + O(k^6),
  Delta_1 = 2 k^3 (eta_o-eta_*)^2 (eta_o+2eta_*)/405 + O(k^5): constant and gradient zeta modes are invisible
  (dilation/SPEC 9d consistent at first order).
- SW plateau: with Delta_l -> -(1/3) j_l(k chi) and P = 2pi^2 A_s/k^3: l(l+1)C_l/2pi = A_s/9 (RD analogue of A_s/25);
  int dy/y j_l^2 = 1/(2l(l+1)) verified.
- k -> infinity: Delta_l -> -cos x_* j_l(k chi) - sqrt3 sin x_* j_l'(k chi) + O(1/(k eta_*)) (verified k = 200, 2000).
- Example (sigma^2=1, eta_*=0.5, eta_o=2, k0=4): Delta_0..4 = 0.44006, -0.33896, 0.27562, 0.22303, 0.11803.

## Milestone 1 — general-l structure of the ISW closed form (status: COMPLETE, verified l = 0..6)
Script `isw_structure.py` (log `isw_structure.log`), explicit forms `isw_closed_forms.tex` (l = 0..6).
THEOREM (all l).  With s_* = k eta_*, S = k eta_o, w = S - s_* = k(eta_o - eta_*), a_± = 1 ± 1/sqrt3,
Cin(z) = int_0^z (1 - cos t) dt/t = gamma + ln z - Ci(z) (entire), J(s) = j_2(s/sqrt3)/(s/sqrt3), Y(s) = y_2(s/sqrt3)/(s/sqrt3):
    I_l = E_l(s_*,S) + A_l(S) [Ci(a_+S) - Ci(a_+s_*) - Ci(a_-S) + Ci(a_-s_*)] + B_l(S) [same with Si]
          + C_l(S) [Cin(a_+w) - Cin(a_-w)] + D_l(S) [Si(a_+w) + Si(a_-w)],
    A_l = (2/sqrt3) Res_{s=0} Y(s) y_l(S-s) = sqrt3 [3 y_l^(3)(S) + y_l^(1)(S)],
    B_l = -(2/sqrt3) Res_{s=0} Y(s) j_l(S-s) = -sqrt3 [3 j_l^(3)(S) + j_l^(1)(S)],
    C_l = -(2/sqrt3) Res_{u=0} y_l(u) Y(S-u),     D_l = -(2/sqrt3) Res_{u=0} y_l(u) J(S-u),
    E_l elementary: Laurent polynomials in (s_*, S, w) times sin/cos of a_± s_*, a_± S, S - s_*, s_*/sqrt3, S/sqrt3.
(Y = -27/s^4 - 3/(2s^2) + regular gives the derivative form of A, B; y_l has a pole of order l+1, so C_l, D_l involve
derivatives of y_2(x)/x, j_2(x)/x up to order l at x = S/sqrt3.)
Derivation: integrand = (1/sqrt3)(g + gbar)(s/sqrt3)(h + hbar)(S - s), g = h_2^(1)(x)/x, h = h_l^(1).  A phase term
e^{ibs}R(s) contributes to Ci/Si with coefficient Res[e^{ibs}R] at each pole (higher poles feed in through the recursion
int e^{ibs}/s^m -> (ib)^{m-1}/(m-1)! int e^{ibs}/s); poles: s = 0 (order 4, from y_2/x) and u = S - s = 0 (order l+1, from y_l).
Since j_l and j_2(x)/x are entire, the four Hankel products pair up: A^- = -A^+, B^- = -B^+, C^- = -C^+, D^- = D^+, leaving the
real residue forms above.  Re Res_{u=0}[j_l(u) g] = 0 is exactly the vanishing of the ln w coefficient (regularity at the
observer: no finite-part ambiguity).
CHECKS (40-digit numerics at 3 points, l = 0..6): every Ci/Si coefficient of the brute-force closed form equals the residue
formula (diff < 1e-169); ln w coefficient = 0 symbolically; E_l contains no Ci, Si, log, gamma.  CAUTION: sympy.residue gave
WRONG answers for these expressions; residues are computed by explicit Laurent x Taylor products.
Example l = 0:  C_0 = (2/sqrt3) Y(S), D_0 = (2/sqrt3) J(S),
    A_0 = -2 sqrt3 (S^3 sin S + 4S^2 cos S - 9S sin S - 9 cos S)/S^4,  B_0 = 2 sqrt3 (S^3 cos S - 4S^2 sin S - 9S cos S + 9 sin S)/S^4;
  E_0 has 56 terms; E_l (l = 1..6): 99, 156, 221, 294, 375, 464 terms (all in isw_closed_forms.tex).

## Milestone 2 — monochromatic C_l^(11), limits, sky monopole, LaTeX (status: COMPLETE)
New code in `linear_cl.py` (section 6b; tests [7], [8]), `plateau_isw.py` (log `plateau_isw.log`); LaTeX section
`linear_cl.tex` (+ `isw_coeffs.tex`, `isw_closed_forms.tex`); it compiles standalone against doc/macros.tex (only the external
ref sec:observable is unresolved outside main.tex).  NOT inserted into doc/ (other agents' files) — the lead should \input it.

### Monochromatic (P = 2 pi^2 sigma^2 delta(k-k0)/k0^2, <zeta^2> = sigma^2)
C_l^(11) = 4 pi sigma^2 [ A(x0*) j_l(k0 chi_*) - k0T_V(x0*) j_l'(k0 chi_*) - delta_l0 A(x0o) + delta_l1 k0T_V(x0o)/3 + I_l(k0) ]^2,
A(x) = -cos x + 2 j_0(x) - 4 j_1(x)/x, kT_V = sqrt3 (sin x - 2 j_1(x)), I_l from Milestone 1 (C11_mono_sym(l) gives it in sympy).

### k0 eta_* -> 0 (Sachs-Wolfe)
RD SW amplitude derived: superhorizon Phi' = 0 and the 00 constraint give delta_1 = -2 Psi_1, delta_T1 = -Psi_1/2;
zeta = delta_T1 - Phi_1 => Phi_1 = -2zeta/3, (delta_T + Psi)_e = Psi/2 = -zeta/3 (matter era: -zeta/5); emission Doppler ~ x -> 0.
At fixed S = k eta_o: Delta_l -> -(1/3) j_l(S) + I_l^0(S) (l >= 2), checked to 1e-10 at k eta_* = 3e-5.
Scale-invariant P = 2pi^2 A_s/k^3 => exact eta_o-independent plateau l(l+1)C_l/(2pi A_s) = 2l(l+1) int dS/S Delta_l(S)^2.
SW alone: 1/9 (RD analogue of 1/25).  BUT in RD Phi decays after horizon entry and every observed mode (S >~ l) has entered,
so ISW (int 2 T_Phi' = +4/3) is not negligible: SW+ISW = 0.3772, 0.4474, 0.5003, 0.5757, 0.6649, 0.7629, 0.8290, 0.8770 for
l = 2, 3, 4, 6, 10, 20, 40, 100.  Large-l limit (flat sky): int_0^1 dmu |Ktilde(mu)|^2 = 0.913769 (ISW alone 0.865111), with the
CLOSED FORM Ktilde(mu) = 1 - 2b^2 - b(1-b^2) ln|(1+b)/(1-b)| + i pi b (1-b^2) theta(1-b), b = sqrt3 mu (Fourier transform of the
line-of-sight source -delta(s)/3 + 2 dT_Phi/ds, via j_1(x)/x = (1/4) int_{-1}^1 (1-t^2) e^{ixt} dt; checked numerically to 3e-4,
limited by truncation).  Convergence to the asymptote is slow (~1/l).

### k0 -> infinity
Delta_l -> -cos x0* j_l(k0 chi_*) - sqrt3 sin x0* j_l'(k0 chi_*) + delta_l0 cos x0o + delta_l1 sin(x0o)/sqrt3, corrections
O(1/(k eta_*)) (checked k = 300, 3000 for l = 0, 1, 2, 7; for l >= 2 the correction is in fact O(1/k^2)).  For fixed l >= 2:
C_l -> 4pi sigma^2 [cos x0* sin(k0chi_* - l pi/2) + sqrt3 sin x0* cos(k0chi_* - l pi/2)]^2/(k0 chi_*)^2 -> 0;
observer terms keep C_0 -> 4pi sigma^2 cos^2 x0o, C_1 -> (4pi/3) sigma^2 sin^2 x0o.

### Sky monopole
<m1> = a_00/sqrt(4pi) = int d^3k/(2pi)^3 zeta_k M(k),
  M(k) = Delta_0(k) = A(x_*) j_0(k chi_*) + kT_V(x_*) j_1(k chi_*) - A(x_o) + I_0(k)     (observer monopole included; M_explicit()),
  M(k) = -k^4 chi_*^3 (eta_o + 3 eta_*)/810 + O(k^6)  (IR safe).
<<m1>^2> = (1/2pi^2) int k^2 P M^2 dk = C_0^(11)/(4pi); monochromatic sigma^2 M(k0)^2 (example eta_*=0.5, eta_o=2, k0=4: M = 0.440060,
<<m1>^2> = 0.193652).
<<m1> a*_lm> = 4pi (-i)^l int k^2dk/(2pi)^3 P M Delta_l int dOmega_k Y_lm(khat) = 0 for l >= 1 (int Y_lm = sqrt(4pi) delta_l0),
for any isotropic P.  Numerically: kernel sky-average = M(k0) for every khat (5e-15); |<<m1>a*_lm>| < 1e-14 (l = 1..3);
<<m1>^2> = C_0/4pi (4e-16).  Caveat: first-order (Gaussian) statement only; at second order <m1> multiplies m1_lm in Th2_lm.

## Status summary
COMPLETE + VERIFIED: Delta_l(k) incl. observer terms (Rayleigh rules, observer Doppler l=1 only, observer monopole l=0 only);
C_l^(11) prefactor (analytic + real-space check); ISW closed form for all l (residue theorem for the Ci/Si/Cin coefficients,
elementary remainder by algorithm), explicit l = 0..6; monochromatic C_l^(11); limits k0 eta_* -> 0 (RD SW amplitude -1/3, plateau
1/9 for SW alone, SW+ISW table and closed-form large-l asymptote 0.913769) and k0 -> infinity; small-k cancellations; sky monopole
M(k), <<m1>^2>, <<m1>a_lm> = 0; LaTeX section.  Full suite `python3 linear_cl.py`: ALL PASS (141 s, run.log = test_log.txt).
NOT DONE / caveats: E_l is given explicitly only for l <= 6 (algorithm general); the large-l plateau asymptote is a numerical
integral of a closed-form |Ktilde|^2 (no closed-form number); the LaTeX section is not \input into doc/main.tex (lead's file).
