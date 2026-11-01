# Report C4 — closed-form time integrals of the second-order regular solutions (RD, Poisson gauge, a = eta)

Code: `cmb4/code/C4_timeintegrals/` (C1's files untouched; C1 kernels loaded through `kernels_api`).
Conventions: SPEC.md §1, §4, §5 and report_C1.md (X = eps X1 + eps^2 X2/2; X2_k = int d^3q/(2pi)^3 K_X(q,k-q;eta) zeta_q zeta_{k-q};
k = k zhat, q = (q sqrt(1-mu^2), 0, q mu), p = |k-q|;  x_k := k eta/sqrt3).

## Milestone 1 — ENGINE (`trig_integrate.py`)  [COMPLETE, verified]
Integrands are sums of  c(k,q,p) eta'^n {sin,cos}(w eta'),  n in [-7, 0],  w in {(+-k +-q +-p)/sqrt3} (Phi2) or {k +- (q +- p)/sqrt3} (h2).
Term-by-term:  trig -> exp,  I_n(w) = int_eps^eta t^n e^{iwt} dt by the recursion I_n = [t^{n+1}e^{iwt}]/(n+1) - (iw/(n+1)) I_{n+1},
I_{-1} = -Cin(w eta) + i Si(w eta) + ln(eta/eps);  the eps-singular pieces (eps^m, ln eps) are collected and REQUIRED to cancel
(checked symbolically: they do, for Phi2 and h2 — this is the regularity of the total integrand at eta'=0).
Special functions used:  Si(z)  and  Cin(z) := EulerGamma + ln|z| - Ci(|z|) = int_0^z (1-cos t)/t dt.
KEY POINT: Cin is entire and even, so the closed forms are manifestly analytic in all frequency combinations: NO absolute values,
NO principal values, NO branch choices.  Rewriting Cin -> Ci reproduces the Kohri-Terada form with its Ci(|w| x) and
log|...| terms (the logs are the ln|w| parts of Cin), and shows where the printed KT (22) needs |.| in the log.
Implemented as a sympy Function `Cin` (derivative, evalf, series, rewrite(Ci)) + numpy evaluator (series for |z|<0.5, sici above).

## Milestone 2 — TENSOR CONTROL CASE  h2 vs Kohri-Terada (22)  [COMPLETE: PASS]
I_h(k,q,p;eta) := int_0^eta G_h K_h deta'  (h2_ij = Lambda_ij,kl q_k q_l I_h;  k-frame: h2_xx = -h2_yy = q^2(1-mu^2) I_h/2, h2_xy = 0;
plus amplitude e^+_ij h2_ij = q^2(1-mu^2) I_h/sqrt2).  Closed form (`build_h2.py`, `h2_closed.pkl`):
   I_h = [sin(k eta) A_h(eta) - cos(k eta) B_h(eta)]/eta,  A_h, B_h = sum over 8 frequencies k +- (q +- p)/sqrt3 of
   [polynomial(1/eta) x trig] + c_w [Cin(w eta), Si(w eta)];  powers eta'^-5..eta'^-1 only (no 1/w poles anywhere).
   Relation:  I_h = (4/k^2) I_RD^{KT}(u = p/k, v = q/k, x = k eta).
CHECK: numerically identical to KT (22) [with |.| in the log] at 12 random points covering u+v < sqrt3 and > sqrt3: rel. diff <= 1e-15
  (with exact rational inputs: 25 digits).  Quadrature arbiter (Gauss-Legendre, 60 digits, 40 panels): agreement to 25 digits.
SUBTLETY (resonance u+v = sqrt3, i.e. q+p = sqrt3 k): in Ci-form the pair  -Ci(|1-(u+v)/sqrt3| x) + log|3-(u+v)^2|  is 0*inf at the
  resonance; the limit is  -EulerGamma - log(x) - log|3-(u-v)^2| + ... ; in Cin-form nothing happens (Cin(0)=0).  KT (22) as printed
  omits the |.| in the log (needed for u+v > sqrt3); their (25) has it.  Exact resonant form built independently by substituting
  p = sqrt3 k - q in the integrand (`build_collinear.py`): it contains an explicit log(eta) term (from the zero-frequency 1/eta' term)
  and agrees with quadrature to 40 digits and with the generic Cin-form evaluated 1e-12 away from the resonance.
  The other resonance |q-p| = sqrt3 k is OUTSIDE the physical domain (|q-p| <= k < sqrt3 k): no Theta(|u-v|-sqrt3) can appear.
NOTE: literature_extracts.md §3.1 quotes "I_RD ~ x^2/2 for small x"; with KT's f_RD(0)=4/3 the correct small-x limit is 2x^2/9
  (ours: I_h -> 8 eta^2/9, i.e. (4/k^2)(2x^2/9)); the printed (22) itself is fine.
Small-eta series (closed form, cross-checked with the ODE Frobenius series):  I_h = (8/9) eta^2 - (2/135)(3k^2+p^2+q^2) eta^4 + O(eta^6).
Late time (`asymptotics.py`):  I_h -> (1/eta)[...]  = (4/k^2) x KT (25) EXACTLY (numerically identical incl. the Theta(u+v-sqrt3) pi cos x term).

## Milestone 3 — SCALAR  Phi2  [closed form COMPLETE, verified; see below for the remaining items]
Phi2 = P + K_Phi2_0(k,q,p) ghat(k eta),  P := int_0^eta G_Phi S_Phi = u2(eta) A(eta) - u1(eta) B(eta),
   A = int_0^eta (u1/W) S_Phi,  B = int_0^eta (u2/W) S_Phi,  u1/W = -eta' sin x' + (k/sqrt3) eta'^2 cos x',  u2/W = -eta' cos x' - (k/sqrt3) eta'^2 sin x'.
Integrand: 200 terms -> 64 (n,w) monomials per integral, n in [-7,0], w in {(+-k +-q +-p)/sqrt3} (8 frequencies); all eps-singular
pieces cancel symbolically.  A, B = sum_w [a_w Cin(w eta) + b_w Si(w eta)] + elementary (trig x Laurent polynomial in eta).
Only the n = 0 monomials give 1/w factors: (e^{iw eta}-1)/(iw)  -> removable 0/0 at the collinear configurations
   w = (k-q-p)/sqrt3 = 0 (mu=1, q<k),  (k+p-q)/sqrt3 = 0 (mu=1, q>k),  (k+q-p)/sqrt3 = 0 (mu=-1)   [these are the ONLY degenerate
   combinations for Phi2: the sound speed enters both G_Phi and the source, so no sqrt3-resonances; signs of the four combinations
   are fixed by the triangle inequalities, so Si(w eta) and Cin(w eta) need no |.| anywhere inside the domain].
   Exact collinear limiting forms built by substituting p = k-q, q-k, k+q in the integrand (`collinear_closed.pkl`): zero-frequency
   terms occur only at n = -4, -2, 0 => NO log(eta) in the collinear Phi2 (unlike the h2 resonance).
CHECK (quadrature arbiter, 60 digits): 8 random points + sqrt3-point: |rel. diff| <= 3.5e-40; 4 collinear points vs limiting forms:
   <= 2.5e-40; generic Cin-form evaluated 1e-12 off the collinear point agrees to 1e-11 (continuity).
