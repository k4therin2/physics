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
