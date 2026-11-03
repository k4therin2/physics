# Referee report R1 (cold read): "The CMB temperature angular power spectrum through fourth order in a radiation-only universe"

Scope: Secs. 1-7, 9 and Apps. A (kernels), A2 (L-term lists), B (ledger) of `cmb4/doc/main.pdf`, freshly compiled
(69 pp., no undefined references). I skipped Secs. 8, 10 and 11 (assembly, monochromatic result, unresolved), which other
authors are still writing. Only the manuscript and the problem statement were read. I did my own re-derivations by hand and in
sympy, in `cmb4/code/R1_referee/` (`check1.py`, `check2.py`, `check_tensor.py`, `check_dirmap.py`, `check_hel.py`,
`check_smallk.py`, `check_flatsky.py`).

Overall: the core of the paper is strong. The exact reformulation T_obs = T_o e^{-S} (Thm 2.3) is the paper's best
idea, and the setup is careful. I found one false claim that matters for the monochromatic (13) term, plus several gaps of
logic and completeness. As written, the document is a project log, not a paper.

---------------------------------------------------------------------------------------------------------------------
## A. What I re-derived and confirmed (no action needed)

| item | result |
|---|---|
| Lemma 2.1 (kinematic redshift), Lemma 2.2 (radiation identities), Thm 2.3 proof, conformal-invariance remark, Cor. 2.4 | correct by hand, including signs and the vorticity drop-out |
| Linear kinematics (2.13), linear fluid eqs (2.14), on-shell identity (2.15), L_1 (2.16), including the observer monopole and the sign of the observer dipole +n.v_o | correct |
| S as an eta-integral (2.12) | correct |
| zeta^(1,2,3) expansions (3.4)-(3.6) | sympy: exact |
| K_Q(0) and K_Phi2(0): equality of both printed forms, K_Phi2 = K_Q/3 - 8/27, the "reader may verify" line, position-space forms (3.9)-(3.10) | sympy: exact |
| C_3: (k,q,mu) form = invariant form; lim_{q->0} C_3 = -32/9; the separate-universe coefficient 6 x 1/2 (4c/3)^2 x (-2/3) | exact |
| Dilation family: a_1 = -c/3, the a_2 = 0 member (-8/9, +8/9, 16/9), the Psi_2 = Phi_2 member (-8/27, 112/27), delta_2 = -2 Psi_2 + 32c^2/9 | exact |
| First-order transfer functions (4.3): ODE, 00 constraint, V, continuity, Euler, small-x limits; the j_n forms in Sec. 5b; 2 dT/d eta = (4k/sqrt3) j_2/x | exact |
| W, G_Phi (both forms), its k->0 limit, G_h: L G = 0, G(eta,eta) = 0, dG = 1 | exact |
| **Tensor source normalisation**: second-order E_xx - E_yy computed from scratch for two plane waves with the linear solution, then -2c/q_x^2 = K_h = 4[2TT + (eta T' + T)(eta T' + T)] with the doc's L_h h = -2 Lambda S and convolution convention | exact (`check_tensor.py`) |
| m_n from e^Lambda; Theta^(2), Theta^(3) for l >= 1 (5.2), including the l = 0 constants | exact |
| Canonical Gram-Schmidt triad (6.4) and second-order direction map (6.6) for n = z (e_3^0, e_3^z, e_3^x, e_1^z, k_2^perp, k_2^z) | exact (`check_dirmap.py`) |
| k_1,init (6.5); clock delta eta_o1, delta eta_o2 (6.11)-(6.12), including 1/u^0 to second order; emission shifts (6.13); eta_o = eta_o/lambda and delta eta_e2 = (2c^2/9 - a_2) eta_* for the constant mode | correct by hand |
| Rayleigh rule (A-s), C^{AA'} = (2/pi) int k^2 P Delta Delta*, sky-monopole rules, C11 prefactor in Prop. 5b.1 | correct |
| Helicity rules (A-lambda), lambda = 1, 2, l = 2..4 | quadrature agrees to 1e-15 |
| Loop prefactors 2 and 3 (5.3); (C1) reduction 1/(4 pi^4 (2l+1)); (B1); (E2), (E3) prefactors (independently re-derived by Wick counting, agree) | correct |
| SW amplitude A(0) = -1/3; ISW integral of 2T' = 4/3; SW-only plateau 1/9; Delta_0 = -k^4 chi_*^3 (eta_o + 3 eta_*)/810, Delta_1 = 2k^3 chi_*^2 (eta_o + 2 eta_*)/405 | exact (`check_smallk.py`) |
| Flat-sky K~(mu) and int_0^1 abs(K~)^2 = 0.913769 | confirmed |
| Lensing control: Theta^(2) = 2 grad phi . grad Theta, Theta^(3) = 3 grad phi grad phi : grad grad Theta, (13) = -3 l(l+1) R C_l | correct |

---------------------------------------------------------------------------------------------------------------------
## B. Findings ranked by severity

### BLOCKING

**B1. The C_3 "integrability" claim is false at q = k (Sec. 3.3, after Eq. 3.12), and this breaks the monochromatic (13) term.**
At q = k, `C_3(k;k,mu) = -4(3mu^3 + 39mu^2 - 45mu + 7) / (9(mu - 1))`. Its expansion is 16/(9(1-mu)) - 56/3 + O(1-mu). This is a
simple pole in mu at mu = 1 (p = abs(k - q) -> 0), so int dmu C_3 diverges logarithmically. The numerator vanishes only to first
order in p^2 against the p^{-4} denominator, so C_3 ~ 1/p^2 there.

- For a smooth P(q), the 3-d integral over d^3q is still convergent, because d^3p/p^2 is integrable.
- For the monochromatic spectrum that the problem requires, q = k = k_0 exactly. Then Phi_3c(k_0, eta) contains
  sigma^2 g(k_0 eta) <C_3(k_0;k_0,mu)>_mu, which is infinite unless the particular integral or other terms cancel the pole.
- Delicate point (2) of Sec. 5.8 mentions mu -> +/-1 but does not resolve it.

Required:
- (a) Correct the statement.
- (b) Determine the 1/(1-mu) coefficient of the full contracted map Theta^(3)_c at q = k (all K_3c, the tadpole and the
  normalisation terms) and show analytically that it cancels in C^(13). Otherwise, define the monochromatic limit as the
  limit of a narrow spectrum and show that limit exists.

Geometric hint: at q = k and mu -> 1, p-hat is perpendicular to k-hat. The (3mu_p^2 - 1) angular structure of the 1/p^2
Poisson terms therefore does not average away.

**B2. The k -> 0 limit of the second-order scalar kernels does not exist (Sec. 3.4, Eq. 3.13 ff.), so the tadpole convention exchanges a limit with an angular average without justification.**
At p = -q + k, `K_Phi2(0) = -[8 - 8 q mu/k + 4 q^2 (3mu^2 - 1)/k^2]/9`. This diverges like q^2/k^2. It is not merely
"direction dependent (proportional to 3mu^2 - 1)", as the text says. The tadpoles tad_X are defined as "angle-averaged k -> 0
limits", which is the angle average of a divergent quantity: the divergent terms average to zero only if the average is taken
before the limit.

- The claim that the convention drops out of C^(13) is supported only numerically, for one plane wave (A1-3). That test also
  has only collinear k-vectors and no vector or tensor modes.
- The problem asks for an analytic justification of every exchange of limits.

Required: an analytic proof that Theta^(3) (clock + emission + tadpoles + normalisation) is independent of the zero-mode
assignment (Phi_2,0, Psi_2,0) along the residual family. The dilation computation C2-9 almost does this already: promote it to a
displayed derivation with a general k = 0 mode.

### MAJOR

**M1. "No independent higher-order tensor modes" is assumed, not derived (Def. 3.2, Sec. 3.2).**
The text says regularity "discards all homogeneous vector and tensor modes". This is false for tensors: the homogeneous solution
h = j_0(k eta) = sin(k eta)/(k eta) of L_h is analytic at eta = 0. The uniform-density argument in Sec. 3.1 then uses
"h_ij -> 0 as eta -> 0" to identify zeta_nl with Maldacena's zeta, which is circular. Gaussianity of zeta_nl fixes the scalars
only (gamma_ij in e^{2zeta}(e^gamma) is a separate variable).

Fix: declare as part of the initial-data definition that h_n,ij(eta -> 0) = 0 for n = 2, 3 (no primordial tensors, gamma = 0 at
every order), alongside zeta^(2)(0) = zeta^(3)(0) = 0. Vector modes are indeed singular (proportional to eta^-2) and are removed by
regularity.

**M2. Initial time eta_i is not kept, contrary to the brief (Secs. 2.4, 6.3).**
The clock, Eq. (6.11)-(6.12), and the observer label drop O(eta_i^2) terms. The label is moved to eta -> 0 rather than to the
physical initial event. The text says these terms are "kept in the numerical test" but gives no formula. The problem requires
finite eta_i, distinct initial, emission and observation times, and all boundary contributions.

- Give the O(eta_i^2) terms analytically, or prove that they vanish as eta_i -> 0 uniformly and state eta_i -> 0 as a declared limit.
- The claim "irrelevant for C_l" in Sec. 6.3 and A1-5 is unproved. An Eulerian-selected and a Lagrangian-selected fluid element
  sample the field with different weights at O(eps^2), so homogeneity alone does not make them equivalent.
- The observer label is part of the definition of the observable. State it in Sec. 2.4, together with the frame of isotropic
  emission (fluid rest frame), which is also an assumption.

**M3. Status statements contradict each other, and "verified" is over-claimed.**
- The abstract still says that Secs. 6, 7 and 9 are placeholders, but they are written.
- The abstract says "every result quoted here has been verified", yet (E2), (E3) and T^(ii) are explicitly unverified (C3-10).
- Sec. 4 marks the kernel time integrals for Phi_2 and h_2 as unevaluated, "being derived". Sec. 9.2 and C4-1/2 say closed
  forms exist and are verified, but they are not printed anywhere.
- Sec. 5.4 says the closed form of C^(11) is in Sec. 8; it is in Sec. 5b.

Fix: one status table, as in the roadmap of Sec. C below.

**M4. Two direction conventions are used side by side.**
- Sec. 5b uses n-hat, e^{+ik.n chi} and i^l.
- Secs. 5.3-5.7 use e-hat = -n-hat, e^{-ik.e chi} and (-i)^l. Sec. 5.3 also tells the reader to "read n-hat as e-hat" inside
  Eq. (5.3), which defines Theta^(r)(n-hat).

I checked that the two conventions give the same transfer function Delta_l, because a^n = (-1)^l a^e is absorbed by i^l versus
(-i)^l. A reader cannot see this, however, and C^(13) mixes Delta_l from Sec. 5b with Delta^(3c) from Sec. 5.7.

Fix: adopt n-hat with i^l throughout (it is the physical convention of the problem), rewrite the rules (A-s), (A-L) and (E1), and
state explicitly that Delta_l is convention independent.

**M5. The document reads as a project log.**
It refers to "agent C1/C2/A1/N1/V3", "SPEC §5/§7/§9a/§9d" (undefined), "report C2", ".pkl" loaders, and a "documentation slip"
remark (Rem. 6.1) about a sign printed in someone else's code. A referee cannot evaluate SPEC references.

Fix: remove all agent names and SPEC references from the main text. Give equations in their final correct form; Rem. 6.1 should
simply be deleted. Move provenance to the ledger and refer to it by check ID only.

**M6. The third-order map is never displayed.**
Sec. 7 gives only group names and term counts (Table 7.1). No explicit formula for L_3 or for the contracted map Theta^(3)_c
appears in the main text. Even App. A2 lists only groups of at most ten terms. C^(13) cannot be verified independently from the
paper.

Minimum fix: display the groups that survive the self-contraction, as named-object formulas: the tadpole group, the variance x
map group, and the cross-contracted group. Also give the q-contracted map in the (A-L) Legendre form (h_L) that Sec. 5.7 needs.
Similarly, the second-order kernel K_2(k1, k2; n) is never written in Fourier space.

**M7. Gauge invariance at fixed physical labels is argued but never tested perturbatively.**
- The exact formula is manifestly covariant.
- The dilation check covers only the k = 0 residual gauge.
- The N1 comparison is actually a gauge-invariance test (comoving-orthogonal gauge against the Poisson-gauge chain), but it is not
  presented as one.

Fix: (i) say this explicitly in Sec. 9; (ii) add a perturbative check under a general second-order time/space gauge
transformation of the assembled L_2, or prove that the labels used (T_phys surfaces, proper time, fluid worldline, Fermi-Walker
transport) make each L_n gauge invariant term by term.

**M8. Analytic treatment of infinite sums and delta-limits.**
- Sec. 5.8 truncates the l1, l2 sums with the numerical heuristic l <~ k0 chi + 10 (k0 chi)^{1/3} + 10.
- xi^{XY}(0) is written as an infinite l' sum.

Both are avoidable analytically. xi^{XY}(0) is a same-point real-space correlator, int d^3k/(2pi)^3 P x(k) y(k) averaged over
angles, so for a monochromatic spectrum it is a closed form; this is the addition theorem, sum (2l+1) j_l^2 = 1. Do this, or give
an analytic tail bound. The problem forbids numerics replacing the derivation.

**M9. Theorem 2.3 needs one more explicit hypothesis, and Fermi-Walker transport is not constructed.**
- State "photons emitted isotropically in the fluid rest frame" in the problem setup, since the problem leaves the frame implicit.
- The FW tetrad is never constructed. The argument is that C_l is invariant under the realisation-dependent rotation R = e^{eps A},
  but (22) and (13) are separately frame dependent. Give the one-line proof that sum_m abs(a_lm)^2 is invariant order by order,
  and that the eps^4 coefficient is (1/4)abs(a2)^2 + (1/3) Re a1 a3*. Then say explicitly that the problem's "specified initial
  orientation + FW" is replaced by the canonical triad without loss.

**M10. The meaning of < >_c is never stated.**
C^(22) includes <<m1>^2> C^(11), which is a product of two 2-point functions. This makes sense only if < >_c means subtracting
<a><a*> (zero for l >= 1), not "connected in zeta". State this explicitly, since the problem uses < >_c.

### MINOR

- m1. Sec. 3.4: a_1 is undefined (it means A = e^{eps a_1}). Sec. 6.3 and C2-9 use lambda = 1 + eps a_1 + eps^2 a_2/2 (polynomial)
  and mu = e^{eps c} for A and B, so a_2 means different things in different places. Unify the notation.
- m2. Notation clashes:
  - g-hat(x) = 3u_1 (homogeneous solution) clashes with the conformal metric g-hat. Rename it, e.g. U(x).
  - S is the shear integral, the source S^(n), and k eta_o in Sec. 5b.
  - A is the dilation amplitude, the FW generator, type-A maps, the SW source A(x), the coefficients A_l, and A_s.
  - M is e^Lambda, M(k) (monopole transfer function), and M_2.
  - Lambda is -S, the TT projector, and a summation index in (A-L).
  - lambda is the affine parameter, the helicity, and the time rescaling.
  - omega is the vorticity tensor and the metric vector.
  - sigma is the shear and sigma^2 the variance.
  - mu is k.q, mu = e^{eps c}, and the m index.
  - K is the kernels, abs(k1 + k2), and the line-of-sight source K(s).
  - h is the tensor mode, the projector h_munu, and the Hankel function h^(1).
  - Rename where cheap, e.g. sigma_zeta^2, Pi_{ij,kl}, M_0(k) -> Delta_0.
- m3. Sec. 5.6: "F the coupling of Hu". Hu's F_{l l1 l2} is not symmetric in l1 and l2, so the correlated cross term needs
  F_{l l1 l2} F_{l l2 l1}; it reduces to F^2 only with the symmetric w_grad of (B1). Write the formula with w_grad.
- m4. Sec. 5.8 (13)(iii) gives two prefactors (12 pi sigma^2 versus 4 pi sigma^2) and "K_3 x Q_3-coefficients". Define Delta^(3c)
  once, unambiguously.
- m5. Sec. 4.1 0i equation: say that the factor 4/eta^2 = 8 pi G a^2 (rho + P) (it is stated only in 4.2).
- m6. Sec. 2.4: the eta_init, eta_i and t_i notations are used interchangeably (also "T-bar_i" versus T_i).
- m7. Eq. (6.17), L2 direct form: the observer uses B_1' (partial) and the emitter uses D B_1. This is correct, because the observer
  displacement is separate and the emitter moves along the ray, but say why in one clause.
- m8. Table 4.1 lists K^0_Phi2 among the "second-order kernels" but it is defined in Sec. 3. Fine, but add an equation reference.
- m9. Sec. 9.2 mentions an erratum about "3/10 vs 1/10 in a helper". This is irrelevant to the paper; drop it.
- m10. main.tex comments still say [PLACEHOLDER] for written sections. The author line "(collaborative derivation; see ledger)"
  is unusual.
- m11. The "Relation to existing work" paragraph should move to an introduction, ahead of the technical setup.

---------------------------------------------------------------------------------------------------------------------
## C. Completeness against the problem statement

| requirement | where addressed | verdict |
|---|---|---|
| Model (flat, radiation, Lambda = 0, rho = alpha T^4) | Sec. 1 | yes |
| Specify the Gaussian nonlinear variable | Def. 3.1-3.2 (zeta_nl = Salopek-Bond/Maldacena zeta, Poisson-gauge limit) | yes |
| No independent higher-order modes | Def. 3.2, Sec. 3.2 | scalars yes; **tensor modes assumed, not derived (M1)**; k = 0 zero mode is a convention (B2) |
| Constraints at every order; do not set higher orders to zero | Sec. 3.2 (regular solution + Bianchi propagation), 3.3 (X_2(0), C_3) | yes; v_2(0) not displayed |
| Emission on T_phys = T_* | Thm 2.3 (T_e = T_* exactly), Eq. (2.11), (6.13) | yes |
| Observer on a fluid worldline, fixed proper time from a physical initial event | Sec. 2.4, 6.3 | partly: O(eta_i^2) dropped and label moved to eta -> 0 (M2) |
| FW tetrad, initial orientation | Sec. 2 (iii) | argued by invariance; not constructed (M9) |
| Start from E_o/E_e; normalise by sky mean; dipole after normalising | Thm 2.3, Cor. 2.4, Sec. 5.1 | yes, cleanly |
| Theta^(1,2,3) with all metric/fluid/propagation/emission/observer/normalisation terms, vectors/tensors, products | Secs. 2.5, 5b, 6, 7, App. A2 | Theta^(1) and Theta^(2) explicit; **Theta^(3) only as term lists in code (M6)**; omega_2, h_2 included |
| a_lm, C_l^(rs), interference terms, normalisations fixed | Secs. 1, 5 | machinery yes; assembled C^(22), C^(13) in Sec. 8 (not reviewed) |
| Analytic; evaluated versus unevaluated marked | \unevaluated markers | inconsistent status (M3); heuristic l-truncation (M8) |
| Finite general times; initial/emission/observation distinguished; boundary terms | Sec. 1 (four times), endpoints Sec. 6.4 | yes except eta_i (M2) |
| Verify constraints, fluid normalisation, null/geodesic, gauge invariance | C1-1/5, C2-5/6, A1-4, dilation, N1 | constraints, normalisation, null yes; gauge invariance only k = 0 + implicit N1 (M7) |
| Monochromatic only after the full observable; what is held fixed | Sec. 5.8 (sigma^2 fixed) | stated; **(13) monochromatic limit divergent at mu -> 1 unless cancellation shown (B1)** |
| Justify exchanges of limits | partial (eta_i -> 0 versus Green's integrals, Sec. 3.2) | missing for k -> 0 versus angle average (B2), delta-function limit (B1), l-sums (M8) |
| Identify unresolved steps | Sec. 11 (not reviewed), Sec. 9.6 "what checks do not cover" | Sec. 9.6 is good; merge it into the status table |

---------------------------------------------------------------------------------------------------------------------
## D. Clarity for the target reader (postdoc in GR/linear perturbation theory)

Where a postdoc gets lost:

1. **The headline result is buried.** Theorem 2.3 is the conceptual heart, but the abstract and Sec. 1 do not say early enough what
   it buys: no separate SW, Doppler, ISW, lensing or intrinsic terms; T_o drops out; only sigma along the ray matters. Put a
   one-paragraph "why this works" in the introduction.
2. **Section order.** Harmonic machinery (Sec. 5) precedes the maps it acts on (Secs. 6-7), and the linear section (5b) sits
   between them. Proposed order: 1 Intro/roadmap; 2 Setup; 3 Exact observable; 4 Initial data; 5 Field solution; 6 Linear map and
   C^(11) as a worked example of the whole pipeline; 7 Second-order map; 8 Third-order map; 9 Harmonic machinery and loop formulas;
   10 Assembly; 11 Monochromatic; 12 Checks (summary table only; details to the appendix ledger); 13 Unresolved.
3. **Sec. 5b is the ideal worked example.** Rewrite it to run explicitly through the same pipeline used later: shear -> L_1 ->
   m_1 -> Theta^(1) -> Delta_l -> C^(11). Then Secs. 7-8 read as the same steps at higher order.
4. **Over-expanded:** the term counts in tables (34/229/359/2484 terms), agent provenance, numeric agreement digits in the
   main text, and the long plane-wave tables. Move these to the ledger.
   **Too abstract:** Sec. 7 (no formulas); Sec. 5.7 items (ii) and (iii) (objects K^C, h_L, M_2 never instantiated); the
   contracted rules (4.13) without a worked example. Add one fully worked contraction, e.g. the tadpole of delta_T2 inside
   delta eta_e2 producing <delta_T2> x Theta^(1)-like map, to show the reader what a (13) contribution looks like.
5. **Definitions to move earlier:**
   - the coefficient convention X_n/n! together with E^(n) = n![eps^n]E (good, but repeat it in the map sections where 1/2 L_2,
     1/6 L_3 appear);
   - the meaning of < >_c (M10);
   - the observer label and emission frame (M2, M9);
   - the n versus e convention (M4), all in a single conventions box in Sec. 1.
6. **Named objects are good.** The "sweet spot" is reached in Secs. 2, 3, 4 and 6 (Table 4.1 of named objects, Table 6.1 of named
   L_2 pieces). Extend the same treatment to Theta^(3): name the survivors (variance x map, tadpole, cross, genuine third order)
   with one formula each in terms of named objects.

### Proposed 2-page roadmap/summary (to place right after the abstract)

**Page 1: what is computed and how.**
- (R1) The target formula, eq. (1.x), with its explicit conventions box: Fourier, P(k), X_n/n!, n versus e, Y_lm phase, < >_c, and
  the monochromatic definition with sigma^2 fixed.
- (R2) Three structural results, one displayed line each:
  - (a) T_obs = T_o e^{-S} and Theta = e^{-S}/<e^{-S}> - 1 (exact);
  - (b) zeta_nl = eps zeta_1 Gaussian, h_n(0) = 0, regular data: X_n(eta_i) are fixed non-local functionals;
  - (c) Theta^(n)_lm in terms of m_n and sky monopoles (5.2).
- (R3) A pipeline diagram: zeta_1 -> {T_X, K_X, K_3c, C_3} -> ray functionals {k_n, Delta x_n, delta eta_o, delta eta_e} ->
  L_n -> m_n -> Theta^(n) -> {Delta_l, K_2,lm, h_L} -> C^(11), C^(22), C^(13) -> monochromatic.
  Annotate each arrow with its section number and the type of object (closed form / declared integral / machine list).

**Page 2: status and requirements.**
- (R4) A requirement -> section table: Sec. C above, condensed.
- (R5) A status table with columns object | evaluated in closed form | declared integral (complete integrand, limits, nesting) |
  verified by | open. Rows: T_X, K_X, the Phi_2 and h_2 time integrals, K_3c, C_3, Phi_3c, L_1, L_2, L_3, Delta_l, C^(11),
  K_2,lm, (E2), (E3), T^(ii), C^(22), C^(13), mono.
- (R6) The three known delicate points, stated as open problems with where they are treated:
  - zero-mode/tadpole independence;
  - the mu -> 1 pole in the (13) loop at q = k;
  - K -> 0 in the (22) loop.
- (R7) A one-paragraph summary of the checks: exact solutions, dilation, the non-perturbative plane-wave comparison to 1e-7 at
  third order (which also counts as a gauge-invariance test), and the lensing control.

---------------------------------------------------------------------------------------------------------------------
## E. Typos / small notation

- Sec. 2.2: "dl = E d lambda" is called the proper-length element. It is the proper length measured in the u-frame; say so.
- Sec. 4.2: "Kelvin theorem" uses u-hat_j with a lower index; specify u-hat_j = g-hat_{j mu} u-hat^mu = omega_j + v_j + ...
- Sec. 5.3: "in this subsection and the next n-hat ... read as e-hat": remove after M4.
- Sec. 5.8: "(effectively finite: ...)": replace (M8).
- Sec. 6, Table 6.1: footnote dagger placement is confusing ("(d) and (e) form one group").
- Sec. 9.3 bullet "V3 finds that the image-summed C_3 shifts by 8": "image-summed" is undefined.
- App. A, v_2,z: "S0z" is typeset as a plain identifier; use S_{0z}.
- Ledger C2-9: "mu = e^{eps c}" should be B.
