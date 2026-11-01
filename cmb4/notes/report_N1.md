# N1 report: nonperturbative planar GR + radiation numerics (independent arbiter)

Code: `cmb4/code/N1_planar/` (README.md there documents formats). Status: **running report, updated per milestone.**

## 1. Formulation (milestone 1 — DONE)
**Gauge.** Planar symmetry (fields of (t,x); y,z Killing; y-z rotation symmetry => B_yy = B_zz, no tensor
polarisation can be sourced by a single plane wave zeta(x)). Comoving-orthogonal (fluid-rest) gauge:
```
ds^2 = -alpha^2 dt^2 + A^2 dx^2 + B^2 (dy^2+dz^2),   u^mu = (1/alpha, 0,0,0)   (irrotational => allowed)
```
Derived with sympy (`derive_equations.py`), not copied: conservation nabla_mu T^{mu nu}=0 gives
* Euler  (nu=x): d_x ln alpha = -(1/4) d_x ln rho = -d_x ln T  =>  alpha T = F(t); the residual relabelling t->f(t) sets F=1.
* Energy (nu=t): d_t ln rho = -(4/3) d_t ln(AB^2)  =>  T (AB^2)^{1/3} = C(x); the regular solution (A,B -> t e^{zeta}, T -> 1/t as t->0) gives C = e^{zeta}.

Hence EXACTLY (no approximation):  `alpha = e^{-zeta(x)} (A B^2)^{1/3}`,  `8piG rho = 3/alpha^4`,  T_phys = 1/alpha.
The whole system reduces to the two metric functions A,B with the lapse fixed algebraically (a Bona-Masso
f=1/3 slicing in disguise: gauge speed = c_s = 1/sqrt3; the "sound wave" is carried by the metric).
Units/conventions = SPEC §1-2: k=1 (period 2pi), a(t)=t, conformal temperature That = tT, That_bar = 1,
8piG=1, 8piG a^4 rhobar = 3.  Hatted variables A = t Ah, B = t Bh, alpha = t alh, alh = e^{-zeta}(Ah Bh^2)^{1/3},
momenta pA = (1/alh) d_t ln Ah, pB = (1/alh) d_t ln Bh (= minus mixed extrinsic curvature of ghat; zero in background).

**Equations** (from G_xx, G_yy; G_tt, G_tx are the constraints; see `derive_equations.py` output / `eqs_generated.pkl`):
```
d_t Ah = alh Ah pA ,  d_t Bh = alh Bh pB
d_t pA = -alh pA^2 - alh pA pB + alh pB^2/2 - 5pA/(3t) + 2pB/(3t) + (alh^2-1)/(2 alh^3 t^2)
         + alh[ (zx^2 - zxx)/Ah^2 - 4 Bx zx/(3Ah^2 Bh) + 5Bxx/(3Ah^2 Bh) - 13Bx^2/(18Ah^2Bh^2) + Ax zx/(3Ah^3) - 11AxBx/(9Ah^3Bh) + Axx/(3Ah^3) - 5Ax^2/(9Ah^4) ]
d_t pB = -3alh pB^2/2 + pA/(3t) - 4pB/(3t) + (alh^2-1)/(2 alh^3 t^2) + alh[ -Bx zx/(Ah^2 Bh) + 7Bx^2/(6Ah^2Bh^2) + AxBx/(3Ah^3Bh) ]
Ham/alpha^2 = [3(alh^2-1) + 2alh^3 pB t + alh^4 pB^2 t^2]/(alh^4 t^4) + 2(pA+pB)/(alh t^3) + 2pApB/t^2 - 2Bxx/(Ah^2Bh t^2) - Bx^2/(Ah^2Bh^2t^2) + 2AxBx/(Ah^3Bh t^2)
Mom = -2 alh pBx - 2zx/t + 2Bx alh (pA-pB)/Bh + 4Bx/(3Bh t) + 2Ax/(3Ah t)
```
Shear (verified symbolically): sigma_xx = (2/3)A^2 s, sigma_yy = sigma_zz = -(1/3)B^2 s, s = (d_t lnA - d_t lnB)/alpha = (pA-pB)/t,
acceleration a_x = d_x ln alpha, theta = (d_t lnA + 2 d_t lnB)/alpha.

**Hyperbolicity** (`hyperbolicity.py`): linearised about FRW, principal symbol in scaled variables (a,b,pA/kappa,pB/kappa)
has eigenvalues +-i/sqrt3 (sound) and a Jordan block at 0 => plain ADM here is only weakly hyperbolic (standard).
Adding cH t^2 (Ham/alpha^2) to the pB equation (parity-preserving) with cH<0, cH != -1/6 gives 4 distinct eigenvalues
+-i/sqrt3, +-i sqrt(2|cH|): strongly hyperbolic. Production uses cH = -1/2 (constraint mode at light speed).
Lower-order terms (1/t, 1/t^2) give kappa-bounded eigenvalue real parts (no exponential ill-posedness) for both.
Numerically cH=0, -1/2, -1 agree to 2.4e-12 in Theta (constraint violations stay ~1e-13 either way).

**Numerics.** Fourier pseudo-spectral in x (N=48 modes; N=32..64 agree to 2e-12), derivatives by FFT, products pointwise;
DOP853 (scipy) rtol=1e-13, atol=1e-16; rtol 1e-12 / 3e-14 change Theta by 6e-12 / 3e-12 (double-precision floor ~3e-12).

## 2. Initial data (milestone 2 — DONE)
Gradient expansion of the regular solution, nonlinear in zeta (`derive_ics.py`, sympy, zeta(x) generic):
Ah = e^{zeta}(1 + sum a_n t^n), Bh = e^{zeta}(1 + sum b_n t^n). Indicial analysis: det ∝ (n+1)^2 n (n-1): n=0 is the
constant (zeta) mode, n=1 is a would-be mode killed by the constraints (a1=b1=0), n=-1 the excluded decaying mode;
odd orders vanish. Results (e.g.)
```
a2 = e^{-2zeta}(4 zeta'' - zeta'^2)/18 ,  b2 = e^{-2zeta}(zeta'' + 2 zeta'^2)/18 ,   a2 + 2 b2 = e^{-2zeta}(2zeta'' + zeta'^2)/6
a4 = e^{-4zeta}(33 z1^4 - 12 z1^2 z2 - 38 z1 z3 + 4 z2^2 + 10 z4)/1080,  b4 = e^{-4zeta}(-8 z1^4 - 2 z1^2 z2 + 7 z1 z3 - z2^2 + z4)/1080
a6, b6: see ics6.log / ics_generated_6.pkl (z_n = d^n zeta/dx^n)
```
Hamiltonian and momentum constraints verified to vanish order by order through the order kept.
Linear check of the expansion: delta_c = alh^{-4}-1 -> -(4/3)(a2+2b2)t^2 -> (4/9) k^2 t^2 zeta = -(2/3)k^2 eta^2 Psi_1 with Psi_1 -> -(2/3)zeta. OK.
Start t_i = 0.01 (k eta_i = 1e-2). t_i = 0.005 / 0.02 with O(t^4) data change Theta by 4e-12 / 3e-12 and t_o by -3e-12 / +4e-11.
Production uses O(t^6) data at t_i = 0.01 (truncation << 1e-12).

**Linear-limit check vs SPEC §5** (`test_linear.py`): Poisson-gauge potentials obtained from my variables through the
linear gauge transformation (derived in `gauge_transform.py`: xi^x = (hA-hB)/k, xi^0 = -(hA'-hB')/k^2,
Phi_P = -hB - (hA'-hB')/(k^2 t), Psi_P = Psi_c + [t(hA''-hB'') + (hA'-hB')]/(k^2 t)), and the gauge-invariant
comoving density contrast delta_c = alh^{-4}-1, compared with T_Phi(x) = -2(sin x - x cos x)/x^3, x = t/sqrt3, and
delta_c = 4 zeta (sin x - x cos x)/x over 0.01 <= t <= 12:
 max|Phi_P/eps - T_Phi| = 1.8e-6 (eps=1e-3), 1.8e-8 (eps=1e-4); Psi: 2.9e-5 / 2.9e-7; delta_c: 4.7e-5 / 4.7e-7.
Residuals scale exactly ∝ eps, i.e. they are the genuine second-order corrections: linear limit verified
(including Phi_P = Psi_P and the sign/normalisation T_Phi(0) = -2/3).

## 3. Observer, emission, ray tracing (milestone 3-4 — DONE, validated)
* Observer: fluid element x = x_o (fixed label). Initial event: T_phys(t_init, x_o) = T_i = 1/t_i (t_init from the
  series solution by root finding; t_init - t_i ~ 1e-10 eps). Clock: Delta tau = (t_obar^2 - t_i^2)/2 (physical proper
  time, = background value); tau integrated as an extra ODE (d tau/dt = alpha(t,x_o)); t_o found by event + Newton
  polish with re-integration (residual < 1e-15).
* Tetrad: u, e_x = d_x/A, e_y = d_y/B, e_z = d_z/B (already orthonormal; FW irrelevant: map is axisymmetric about x).
  Sky direction n = (mu, sqrt(1-mu^2), 0): photon momentum p = E_o(u - n^a e_a), E_o = 1, i.e. p_x = -mu A, p_y = -sqrt(1-mu^2) B,
  p_t = -alpha at the observer. mu = n_x > 0 <=> photon arrives FROM the +x side (travels in -x direction).
* Rays: Hamilton's equations with p_y conserved, parametrised by t, integrated BACKWARDS together with the PDE
  (the system is time-reversible, so the fields are available exactly at the ray's time; x interpolation is spectral).
  Emission event along each ray: alpha(t,x_r) = t_* (<=> T_phys = T_* = 1/t_*); bracketed on integrator steps,
  brentq on dense output, then Newton with exact re-integration (|alpha_e - t_*| < 5e-16). Null constraint check:
  E from p_t vs from (p_x,p_y) agree to 4e-15.  T_obs = T_* E_o/E_e, E = -p_t/alpha.
* Shear integral S = int_e^o sigma_{mu nu} e^mu e^nu dl, dl = alpha dt, sigma_ee = s (e_x^2 - 1/3), integrated along each ray.
* Sky: Gauss-Legendre in mu (nmu=96 in production; 64 gives a_l0 to 5e-13 for l<=24, 32 is NOT enough), map is a
  function of mu only (axisymmetry), a_l0 = sqrt(4pi(2l+1)) (1/2) int Theta P_l dmu (real Y_l0 normalisation, polar axis = +x).
  For x_o = 0 the map is even in mu (reflection symmetry) => odd l vanish (dipole ~ 2e-13 numerically).
* Validation: eps=0 (FRW): T_obs/T_o - 1 = 2e-16 for all directions, S = 0, t_o = t_obar exactly.

## 4. BONUS CHECK: exact formula T_obs(n) = T_o exp(-S(n)) (SPEC §3)  — VERIFIED NONPERTURBATIVELY
max_n |T_obs/T_o - exp(-S)| = 7.5e-15 (eps=0.01), 1.4e-14 (eps=0.1, Theta_rms = 0.077), 8e-15 (eps=0.03).
This ties together the geodesic energies, the local fluid temperatures at emission and observation, and the shear
integral: the lead's exact formula holds at the double-precision level at fully nonlinear amplitude.

## 5. Production runs / extraction  — IN PROGRESS
(see below when done)

## 5. Production runs and eps-Taylor extraction (milestones 5-7 — DONE)
Parameters (all runs): k=1, t_i=0.01, t_* = 2 (x_* = k eta_*/sqrt3 = 1.155), t_obar = 8 (eta_o/eta_* = 4) and 4 (ratio 2),
observer x_o = 0 (primary; sits at the maximum of zeta) and x_o = pi/2 (generic: zeta(x_o)=0, zeta' maximal, non-zero
dipole), N=48, nmu=96, O(t^6) ICs, DOP853 rtol 1e-13, cH=-1/2.  Amplitudes eps = +-{0.01,0.02,0.03,0.04,0.05,0.06,0.08,0.10}
(16 runs per configuration, ~8 s each).  Every run: exact-formula residual max|T_obs/T_o - e^{-S}| <= 7.5e-14,
constraint violations (relative to largest term) < 1e-9 at t_o (absolute ~1e-13), |alpha_e - t_*| < 3e-15, null condition 1e-14.
Symmetry checks (x_o=0, eps=0.05): Theta(-eps; x_o=0) vs Theta(+eps; x_o=pi): 5.4e-12; parity Theta(mu)=Theta(-mu): 5.5e-12.
For x_o = pi/2: Theta(-eps, mu) = Theta(+eps, -mu) holds exactly (t_o, T_o, monopole identical for +-eps) => odd-eps
coefficients are odd in mu (odd l only), even-eps coefficients even in mu.

**Extraction.** E(e) = [Theta(e)+Theta(-e)]/2 = c2 e^2 + c4 e^4 + ..., O(e) = [Theta(e)-Theta(-e)]/2 = c1 e + c3 e^3 + ...;
least-squares polynomial fits in e^2 on several amplitude subsets (A: e<=0.04/4 terms, B: e<=0.06/5, C: e<=0.10/6,
D: e<=0.08/6 [reference], F: e<=0.05/4). SPEC normalisation Th1 = c1, Th2 = 2c2, Th3 = 6c3. Uncertainty := max over the
map of the spread between subsets (conservative; the three best subsets B,C,D typically agree 5-20x better).
Full tables (a_l0 for l<=12 per configuration, t_o(eps), T_o(eps), <T_obs>(eps)): `code/N1_planar/output/summary.md`;
maps/a_l0 (l<=30): `output/coeff_<tag>.npz/.json`; raw runs: `output/scan_<tag>.npz` (formats in README.md).

| configuration | Th1 rms / max | Th2 rms / max | Th3 rms / max | uncertainty Th1 / Th2 / Th3 (abs, max over map) |
|---|---|---|---|---|
| to8_xo0  (ratio 4, x_o=0)    | 0.7111 / 1.076 | 1.309 / 2.708 | 14.21 / 37.23 | 3.6e-9 / 3.2e-8 / 2.1e-4 (best subsets B,C,D: 1.2e-5) |
| to4_xo0  (ratio 2, x_o=0)    | 0.06503 / 0.09907 | 0.3955 / 1.258 | 1.644 / 6.816 | 6.8e-12 / 9.7e-10 / 4.0e-7 |
| to8_xo0.5pi (ratio 4, x_o=pi/2) | 0.8428 / 2.106 | 2.145 / 4.871 | 24.13 / 38.20 | 5.3e-7 / 9.2e-6 / 3.1e-2 (B,C,D: ~8e-3) |
| to4_xo0.5pi (ratio 2, x_o=pi/2) | 0.09388 / 0.2649 | 0.6417 / 1.658 | 4.078 / 6.319 | 6.7e-10 / 5.4e-9 / 3.9e-5 |

Primary configuration to8_xo0, a_l0 (real Y_l0, polar axis +x; odd l vanish by parity, numerically <2e-11 for Th1):
```
 l    a_l0[Th1]            a_l0[Th2]            a_l0[Th3]
 2   -2.1847203553e+00    -4.2224942193e+00    -1.288998797e+01
 4   +1.2552456282e+00    +1.3785439781e+00    -4.707090710e+01
 6   -4.7973477364e-02    -1.2595382817e+00    +6.290510536e+00
 8   -6.2592963739e-02    -6.4106435777e-02    +6.924507646e+00
10   +1.4257282317e-02    -1.1039960997e-01    -5.605800149e+00
12   -1.5571675075e-03    +3.8560168349e-01    +3.831186883e+00
Th1(mu=+1) = Th1(mu=-1) = -0.4221462,  Th1(mu=0) = +1.07625651
```
Observer clock/temperature (to8_xo0): t_o(+-0.01) = 8.0251226534 / 7.9742196483 (O(eps) shift because zeta(x_o)=eps != 0),
T_o(+-0.01) = 0.12445094278 / 0.12554985919, <T_obs> = 0.12499073752 / 0.12499056880.  For x_o = pi/2 the shifts are O(eps^2).
Th4 (indicative, 24 c4) is also stored; its uncertainty is ~1e-4 (ratio 2) to ~4e-3 (ratio 4, x_o=0).

**Error budget (observable T_obs(n)/T_o, single run).** Double-precision floor: changing N (32/48/64), rtol (1e-12..3e-14),
cH (0,-1/2,-1), t_i (0.005/0.01/0.02), IC order (4/6) changes Theta by 2e-12 .. 8e-12 (max over sky) and t_o by <5e-11
(relative <1e-11); exact e^{-S} identity 1e-14. => T_obs/T_o accurate to ~1e-11 relative (target 1e-9 met by 2 orders).
Angular: nmu=64 already gives a_l0 (l<=24) to 5e-13; production nmu=96.
Extraction: noise floor sigma_Theta ~ 5e-12 propagates as ~sigma/e0^n (e0 = 0.01) with O(1) weights: Th1 ~1e-9, Th2 ~1e-7,
Th3 ~1e-5 for the ratio-4 runs (consistent with the subset spreads); truncation (c_{2m+1} e^{2m}) is controlled by the
subset comparison (the ratio-4 maps have fast-growing Taylor coefficients, c4 ~ 3.5, so e<=0.08 is near the useful limit).

## 6. What is verified (summary)
1. Equations derived symbolically; conservation identities and shear verified; hyperbolicity analysed (strongly hyperbolic variant used).
2. FRW limit exact (2e-16). Constraints ~1e-13 absolute throughout.
3. Linear limit = SPEC's T_Phi (Phi_P = Psi_P = zeta T_Phi, delta_c) with residuals ∝ eps (2nd-order), to 1e-8 at eps=1e-4.
4. Nonlinear regular initial data (gradient expansion) consistent with constraints to O(t^6); IC-order independence 8e-12.
5. Exact observable identity T_obs = T_o e^{-S} (SPEC §3) verified to 1e-14 at eps up to 0.1 (Theta_rms 0.08) — nonperturbative.
6. Symmetries (zeta sign flip = half-period shift; parity) at 5e-12.
7. Th1, Th2, Th3 maps/a_l0 for 4 configurations with quantified uncertainties (Th2 to ~1e-8..1e-5, Th3 to ~1e-6..1e-3 relative).

Caveats: single plane wave => only the m=0 (axisymmetric) part of the kernels is probed (for x_o=0 only even l;
x_o=pi/2 adds odd l and the dipole); eta_o/eta_* is 2 and 4 (not CMB-like 50); monopole normalisation and the
observer-clock shift are included exactly, so the analytic comparison must include SPEC §6-7 bookkeeping.
