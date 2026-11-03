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

## Milestone 4 — ASYMPTOTICS AND SMALL-eta SERIES  [COMPLETE, verified]  (relaunch 2026-11-03)
### eta -> infinity (`asymptotics.py` -> `asymptotics.pkl`; `asymptotics_compact.py`; check `probe_asym2.py`)
Cin(z) -> EulerGamma + ln|z| - sin z/z + cos z/z^2 + 2 sin z/z^3 - 6 cos z/z^4,  Si(z) -> (pi/2)sgn z - cos z/z - sin z/z^2 + ...,
with the frequency signs fixed by the triangle inequalities.  The "sign error" suspected by the crashed run is NOT present in the
saved expansion: re-checked with a 80-digit mpmath evaluation of the full closed forms at k eta = 1e2, 1e3, 1e4, 1e5 for three
configurations, (full - asym) x^4 (Phi2, kept eta^-2, eta^-3) and (full - asym) x^3 (h2, kept eta^-1, eta^-2) stay bounded
(no growth over three decades) => the expansions are correct through the kept orders.  [The earlier log only looked bad because
for q/k = 0.1 the next coefficient is ~1e4 (1/(p^3 q^3) prefactors) — the expansion is non-uniform as q -> 0 and near q+p = sqrt3 k.]
LEADING TERMS (x = k eta/sqrt3):
  Phi2 -> (1/eta^2) { C_L [ Lambda cos x + pi sin x ] + C_0 cos x } + O(eta^-3),
     Lambda = ln[ (k+q+p)(q+p-k) / ((k+q-p)(k+p-q)) ],
     C_L = -3 (k^2+p^2+q^2) R / (8 k^4 p^3 q^3),
     R = 3k^6 - 3k^4(p^2+q^2) - 3k^2(p^4+q^4) - 2k^2 p^2 q^2 + 3(p^6+q^6) - 3 p^2 q^2 (p^2+q^2),
     C_0 = 3[ -3k^6 - 3k^4(p^2+q^2) + 3k^2(p^4+q^4) + 2k^2 p^2 q^2 + 3(p^6+q^6) - 3p^2q^2(p^2+q^2) ] / (2 k^4 p^2 q^2)
     (C_0 includes the homogeneous piece K_Phi2_0 ghat -> -9 K_Phi2_0 cos x/(k eta)^2).  No Theta function: for Phi2 all four
     frequency signs are fixed in the physical domain; Lambda is real and finite except at the collinear edges (log divergence of
     the *asymptotic* coefficient only; the exact closed form is finite there, see Milestone 3).
  I_h -> (4/k^2) x KT(25)  (sin(k eta) log term + Theta(q+p-sqrt3 k) pi cos(k eta) term), identical to Kohri-Terada (25).
### eta -> 0 (`series_small.py` -> `series_small_eta.pkl`, log `series_small.log`)
Frobenius recursion a_n n(n+3) = s_{n-2} - (k^2/3) a_{n-2} on the exact source series + C1 constraints; all fields regular:
  Phi2   = -(7k^4 - 2k^2(p^2+q^2) + 3(p^2-q^2)^2)/(9k^4)
           + eta^2 [-21k^4(p^2+q^2) + 10k^2(p^4+q^4) - 28k^2p^2q^2 + 3(p^6+q^6) - 3p^2q^2(p^2+q^2)]/(270 k^4)
           - eta^4 [-150k^6(p^2+q^2) - 65k^4(p^4+q^4) - 398k^4p^2q^2 + 80k^2(p^6+q^6) - 144k^2p^2q^2(p^2+q^2)
                    + 15(p^8+q^8) + 36p^2q^2(p^4+q^4) - 102p^4q^4]/(113400 k^4) + O(eta^6)
  I_h    = (8/9) eta^2 - (2/135)(3k^2+p^2+q^2) eta^4 + (45k^4+15k^2(p^2+q^2)+5(p^4+q^4)+22p^2q^2) eta^6/42525 + O(eta^8)
  Psi2, Q2, delta2, deltaT2 to eta^4 and V2 to eta^5: in series_small.log (eta^0 terms reproduce PUB/C1 ICs; zeta^(2)(0) = 0 checked).
CHECK: closed forms (mpmath 120 digits) minus the series through eta^6: (closed - series)/eta^8 -> constant for eta = 1e-1..1e-3
  (two configurations) for both P and I_h.

## Milestone 5 — ALL SECOND-ORDER FIELDS: END-TO-END EINSTEIN CHECK  [COMPLETE: PASS]
`verify_einstein.py` (log `verify_einstein.log`): first-order plane-wave pair zeta = A e^{i q1.x} + B e^{i q2.x} plus the C4
second-order fields (Phi2, Psi2 = Phi2 - Q2, delta2, v2 = grad V2 + v2^T, omega2, h2_ij from I_h; coefficient of A B e^{ik.x},
summed over both orderings) inserted into C1's FULL second-order Einstein equations E^(2)_{mu nu} (10 components, from
E_sympy.pkl, i.e. independent of the SVT decomposition used to derive the kernels) and the 4 second-order conservation equations;
eta-derivatives of the closed forms by mpmath.diff at 50 digits.  Three points (k eta = 2.5, 9.6, 20.1; q/k from 0.19 to 1.6):
  max |residual| / max|individual term| = 5.5e-50 over all 14 equations x 3 points.
(This simultaneously checks the closed-form Phi2, its eta-derivatives, all constraint-derived fields, the vector fields and h2.)

## Milestone 6 — COMPACT STRUCTURE + LaTeX  [COMPLETE]
`structure5.py` (-> `phi2_compact.pkl`, log `structure5.log`), `structure4.py` (h2 part -> `h2_compact.pkl`; its Phi2 stage was too
slow and was superseded by structure5), `make_tex.py` -> `C4_closed_forms.tex/.pdf` (16 pp; structure in main text, full forms in
the appendix, all formulas generated from the pickles).
KEY STRUCTURAL RESULT: P = u2(x_k) A(eta) - u1(x_k) B(eta) and the whole special-function content of Phi2 is governed by ONE coefficient
   D(k,q,p) = (k^2+p^2+q^2) R / (8 k^2 p^3 q^3)   (R as in the asymptotics):
   A = D [Si(w1 eta) - Si(w2 eta) - Si(w3 eta) - Si(w4 eta)] + sum_j [C^A_j cos(w_j eta) + S^A_j sin(w_j eta)],
   B = -D [Cin(w1 eta) - Cin(w2 eta) - Cin(w3 eta) + Cin(w4 eta)] + sum_j [C^B_j cos + S^B_j sin] + c^B_0,
   w1..4 = (k+q+p), (k+q-p), (k-q+p), (q+p-k) over sqrt3 (all >= 0 in the triangle), C_j, S_j Laurent polynomials eta^-6..eta^0
   (eta^0 terms carry the removable 1/w_j collinear poles).  No ln(eta) survives.  Rebuild from the compact form = stored P (30 digits).
   (Consistency: the eta->infinity coefficient C_L = -(3/k^2) D ... matches Milestone 4.)
