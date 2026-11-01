# N1_planar — exact planar GR + radiation-fluid solver, ray tracing, and eps-Taylor extraction

Independent nonperturbative numerical arbiter for the analytic O(eps^2), O(eps^3) CMB-map calculation (SPEC.md).
Physics summary and all derivations: `../../notes/report_N1.md`.

## Files
| file | purpose |
|---|---|
| `derive_equations.py` | sympy: Einstein + fluid equations in planar comoving-orthogonal gauge -> `eqs_generated.pkl` |
| `derive_ics.py N` | sympy: gradient expansion (regular solution, nonlinear in zeta) to O(t^N) -> `ics_generated_N.pkl` |
| `hyperbolicity.py` | principal-symbol analysis (plain ADM weakly hyperbolic; cH-modified strongly hyperbolic) |
| `gauge_transform.py` | linear gauge transformation of this gauge -> Poisson gauge (for the linear check) |
| `gen_numeric.py` | lambdifies the sympy results (with exact rewrites avoiding background cancellations) |
| `planar.py` | the solver: `PlanarRun(eps, ...).run()` (spectral + DOP853 evolution, observer clock, ray tracing) |
| `test_smoke.py`, `test_linear.py`, `test_convergence.py`, `test_symmetry.py` | validations (logs `*.log`) |
| `run_scan.py t_obar x_o/pi tag [ic_order t_i]` | production amplitude scan -> `output/scan_<tag>.npz` |
| `extract.py tag` | even/odd polynomial extraction of Taylor coefficients -> `output/coeff_<tag>.npz/.json` |

Regenerate the pickles (gitignored) with `python3 derive_equations.py; python3 derive_ics.py 6` (~1 min + ~20 min).

## Conventions (all SPEC-compatible)
* Units: k = 1 (zeta period 2 pi), a(t) = t with t the coordinate time of the comoving-orthogonal gauge
  (= conformal time in the background), conformal temperature That = t T_phys, That_bar = 1, 8 pi G = 1,
  so 8 pi G a^4 rhobar = 3.  T_phys = 1/alpha exactly (alpha = lapse).
* Primordial curvature: zeta(x) = eps cos(x), defined by g_ij -> a^2 e^{2 zeta} delta_ij as t->0 in comoving
  coordinates (no eps^2 corrections, no decaying/tensor modes). Fourier: zeta_k = eps (2pi)^3 [delta^3(k - k0 x^) + delta^3(k + k0 x^)]/2 with k0 = 1.
* Emission surface: T_phys = T_* = 1/t_star, t_star = 2 (k eta_* = 2, x_* = k eta_*/sqrt3 = 1.155).
* Observer: fluid element with label x_o (0, or pi/2 for the generic runs), initial event T_phys = T_i = 1/t_i, t_i = 0.01,
  observation after proper time Delta tau = (t_obar^2 - t_i^2)/2 ; t_obar = 8 (eta_o/eta_* = 4) and t_obar = 4 (ratio 2).
  The actual coordinate time t_o(eps) and the local temperature T_o(eps) are stored.
* Sky direction n: unit vector in the observer's fluid-frame triad (e_x = d_x/A, e_y = d_y/B, e_z = d_z/B);
  the photon ARRIVES from direction n, i.e. its spatial momentum is along -n. The map depends only on
  mu = n.e_x = cos(theta) (axisymmetry about x). mu > 0 : photon came from the +x side.
* Theta(n) = T_obs(n)/<T_obs>_sky - 1 with <.>_sky = (1/2) int dmu (exact Gauss-Legendre). Dipole kept.
* a_l0 = sqrt(4 pi (2l+1)) (1/2) int_{-1}^{1} f(mu) P_l(mu) dmu  (= int dOmega Y_l0 f, real Y_l0, polar axis +x; m != 0 vanish).
* Taylor convention: Theta = eps Th1 + eps^2/2 Th2 + eps^3/6 Th3  (Th1 = c1, Th2 = 2 c2, Th3 = 6 c3 with Theta = sum c_n eps^n).
* S(n) = int_e^o sigma_{mu nu} e^mu e^nu dl (fluid shear in fluid frame along the ray, dl = E dlambda = alpha dt);
  exact claim T_obs/T_o = exp(-S) (SPEC §3) is checked in every run (`exactS_res`).

## Output formats
`output/scan_<tag>.npz` (one entry per run, arrays stacked along axis 0 over `eps`):
`eps`, `T_obs` (nmu), `T_obs_null` (same, E_e from null condition), `T_o`, `t_o`, `t_init`, `t_e` (nmu, emission times),
`x_e` (nmu, emission label), `S` (nmu), `Theta` (nmu), `monopole` (<T_obs>), `alm` (31: a_l0 of Theta, l=0..30),
`exactS_res` (max|T_obs/T_o - e^{-S}|), `H`,`Hrel`,`M`,`Mrel` (constraint violations at t_o, absolute and relative to the
largest term), plus `mu`, `wmu` (Gauss-Legendre nodes/weights), `t_star`, `t_obar`, `x_o`, `params` (json).

`output/coeff_<tag>.npz`: `mu`, `wmu`, maps `Th1`, `Th2`, `Th3`, `Th4` (SPEC normalisation), `alm_Th1..4` (l=0..30),
`spread_Th1..3` (max over fit subsets = extraction uncertainty), `fit_subset`, `t_star`, `t_obar`, `x_o`.
`output/coeff_<tag>.json`: same summary in text (rms, spread, a_l0 lists, maps).
Tags: `to8_xo0` (t_obar=8, x_o=0, O(t^6) ICs, PRIMARY), `to4_xo0`, `to8_xo0.5pi`, `to4_xo0.5pi`, `to8_xo0_ic4` (O(t^4) ICs, for the IC-error estimate).
