# Report W1 — settled sections of the final document (cmb4/doc/)

Status: DONE for the assigned sections. No TeX engine in the container (`which pdflatex` empty); syntactic sanity
checked with a custom script (brace/environment/$ balance, \left/\right, \input existence, all \ref/\eqref labels and
\cite keys resolved): PASS. Compile with `pdflatex main` (x2) when available; packages: amsmath, amssymb, amsthm,
geometry, xcolor, hyperref, breqn (dmath for machine-generated kernels).

## Files (all in cmb4/doc/)
- `main.tex` — article skeleton; \input of every section; placeholders clearly marked.
- `macros.tex` — notation macros, `\todo{}` (red "[To be supplied: ...]"), `\unevaluated{}` (blue "[unevaluated integral: ...]").
- `sec_setup.tex` (Sec. 1) — model, a = eta units, conformal variables, Fourier/statistics, expansion convention X = eps X1 + eps^2 X2/2 + ...,
  E^{(n)} = n! [eps^n] E, Poisson gauge (linear parametrisation, dictionary to BMR exponentials), direction convention
  (n = observation direction, e = -n propagation, background ray x = x_o + n chi, D = d_eta - n.grad), the four times, a_lm/C_l definitions, literature paragraph.
- `sec_observable.tex` (Sec. 2) — conformal frame; Lemma 1 (kinematic redshift), Lemma 2 (radiation identities, both P=rho/3 and rho∝T^4 needed),
  Theorem T_obs = T_o e^{-S}, remark on conformal invariance (S = Shat), Corollary 1+Theta = e^{-S}/<e^{-S}>; consequences (i)-(iv) incl. V1's C6
  (C22, C13 separately frame dependent → canonical triad at all orders); validity remarks (V1 C1-C3); observer/clock/emission/triad/direction-map
  definitions (SPEC §6) and S in eta-parametrisation; linear order: sigma_ee = D(delta_T+Psi-n.v) - (Phi'+Psi') on shell → L1 = SW+Doppler+ISW with observer terms.
- `sec_primordial.tex` (Sec. 3) — zeta_nl definition, uniform-density-slice argument (= Salopek–Bond/Maldacena zeta), Gaussianity = zeta^(2)(0)=zeta^(3)(0)=0,
  regularity, finite-eta_i data as restriction of the regular solution; zeta^(2), zeta^(3) expansions; derivation of K_Phi2^0 from zeta^(2)(0)=0
  (reader-checkable one-liner), position-space form, PUB agreement, BMR difference; C3(k;q,mu) boxed; exact dilation family (A,B), k=0 residual
  freedom, angle-averaged member = {A = e^{-c/3}, B = e^{c}} (I verified by hand that C1's numbers Phi2=-8/9, Psi2=+8/9, delta2=16/9 lie on it),
  tadpole convention and why it must cancel; squeezed C3 check interpreted as long-mode modulation.
- `sec_einstein.tex` (Sec. 4) — table of named objects; first order (equations, T_Phi, T_delta, T_V); operator structure (universal L_Phi, L_h, L_omega,
  SVT-projected sources, the ten equations); relativistic Kelvin (algebraic omega2); kernels K_Q, K_om, K_h (mu-form) + pointers; Green's functions
  (G_Phi rewritten as (sqrt3/k)(x'/x^3)[(1+xx')sin(x-x') - (x-x')cos(x-x')], verified against C1's form and L G = 0, G(eta,eta)=0, dG=1);
  declared integrals K_Phi2, K_h; third order contracted construction, rules, tadpoles, equations for X3c, C3; summary of what is settled.
- `sec_harmonic.tex` (Sec. 5) — SPEC §7 bookkeeping (Theta1..3 boxed), loop formulas with C3's confirmed prefactors, direction-convention statement
  (rules written with e-hat = -n-hat = C3's "n"; a_lm differ by (-1)^l, C_l unchanged), type A (A-s, A-lambda, C11, monopoles), type B (B1 both pairings,
  monopole×map (i)-(iii), 1/4 C22 formula), type D, type C (K_2lm, frame form, (C1),(C2)), (13): (E0), (A-L), (E1), structures (i)-(iii), (E2),(E3);
  lensing control statement; monochromatic definitions (P = 2pi^2 sigma^2 delta(k-k0)/k0^2, sigma^2 fixed) and collapsed rules; delicate points (1)-(4).
- `sec_map2.tex, sec_map3.tex, sec_assembly.tex, sec_checks.tex, sec_mono.tex, sec_unresolved.tex` — PLACEHOLDERS (one \todo each saying what is expected).
- `app_kernels.tex` — machine generated from C1 pickles via kernels_api (script in scratchpad; symbols T_q, \dot T_q): K_Q, K_S00, K_S0z, K_SPhi, K_om, K_vT, K_h,
  K_Phi2^0, K_Q(0), zeta2 quadratic kernel, delta2/v2z algebraic forms, u1, u2, ghat, G_Phi, G_h, G_om, C3 (mu and invariant), Q3c(0), lim eta^2 K3c_S00, five tadpoles.
  NOTE: C1's report wrote the invariant K_Q without its overall minus sign; the appendix uses the machine form (verified against the mu-form).
- `app_ledger.tex` — verification ledger: V1-1..6, C2-1..6 (milestones 1-2, marked), C1-1..14, C3-1..10, N1-1..2, pending list.
- `refs.tex` — thebibliography, 25 entries (all cited).

## Conventions I fixed in the document (later agents please respect)
- n-hat = observation direction, ray x = x_o + n chi (SPEC §6 / C2). Harmonic rules in Sec. 5 are stated for maps written via e-hat = -n-hat
  (identical to C3's code, where this vector is called n). Translation: a_lm -> (-1)^l a_lm at every order, all C_l unchanged.
- T_q := T_Phi(q eta/sqrt3), \dot T_q := d_eta T_q in all kernels; p = |k - q|; mu = k.q-hat.
- "Declared integrals": K_Phi2(k,q,p;eta) = int_0^eta G_Phi K_SPhi + K_Phi2^0 ghat, K_h = int_0^eta G_h K_h, and the third-order int G_Phi K3c_SPhi.

## What remains for later sections
Secs. 6-11 (second/third-order map from C2, assembly incl. C^(11) closed form from C5, checks incl. dilation/N1/V3, closed-form time integrals from C4,
monochromatic result, unresolved list). Sec. 5.6 has a \todo for the explicit h_L, K^C, M_2 of our map.
