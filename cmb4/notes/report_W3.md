# Report W3 — revision of the manuscript per referee report R1

Agent W3 (2026-11-03). Owns: main.tex, sec_setup, sec_roadmap (new), sec_observable, sec_primordial, sec_einstein,
sec_harmonic, sec_linear, sec_map2, sec_map3, sec_checks, app_kernels, app_lterms, app_ledger, app_closedforms (new), refs, macros.
Not touched: sec_assembly, sec_mono, sec_unresolved, app_k2channels.

## Milestone 1 — code (code/W3_revision/)  [DONE]
* clock_eta_i.py (+ .log): observer clock with finite eta_i, label at the initial event. Result
  deta_o1 = [eta_i^2 dT1(eta_i) - int_{eta_i}^{etabar_o} eta Psi1]/etabar_o;
  deta_o2 = -(2/etabar_o){deta_o1^2/2 + etabar_o deta_o1 Psi1|_o + int_{eta_i}^{etabar_o} eta[Psi2/2 - Psi1^2/2 - v1^2/2 + dx1.grad Psi1]}
            + (eta_i^2/etabar_o)[dT2 + dT1^2 + 2 Psi1 dT1 + 2 eta_i dT1 dT1']_{(eta_i, x_lab)};
  label at the initial event: dxo2/2 gains -eta_i dT1(eta_i) v1(eta_i) = O(eta_i^2).  Identical to A1's pw_map.observer
  (which keeps eta_i); dilation family eta_o = etabar_o/lambda holds for finite eta_i (symbolic).
* check_rules.py (+ .log): harmonic rules (A-s), (A-lambda), (A-L) rewritten in the n-hat convention verified by sphere
  quadrature (<= 6e-14); new (A-lambda) helicity factors: Dhat_{+-1} = +- i sqrt(l(l+1)/2) j_l/x (sign flipped relative to the
  e-hat convention), Dhat_{+-2} unchanged; bound |j_l(x)| <= x^l/(2l+1)!! <= (e x/(2l+1))^l checked on a grid.
* gen_lterms_full.py: extends W2a's generator; writes doc/app_lterms.tex with ALL terms (L1 3, L2/2 34, direct L2/2 76, L3/6 359).

## Milestone 2 — revision of the text  [DONE; pdflatex x2: 0 errors, 0 undefined refs]
* main.tex: R1's order (Intro/Setup, Roadmap, Observable, Primordial, Einstein, Linear [worked example], Map2, Map3, Harmonic,
  Assembly, Mono, Checks, Unresolved; apps kernels, closedforms, lterms, k2channels, ledger+reproducibility); new abstract (status
  consistent; analytic/numerical/pending listed); author line; tikz, multirow.
* sec_setup: new Introduction ("why this works", relation to existing work, organisation); setup with assumptions (O1) emission isotropic
  in fluid rest frame, (O2) observer label at initial event + FW frame, (O3) normalisation; <.>_c defined (M10); single direction
  convention + Rayleigh eq (M4); four times; notation table (renames).
* sec_roadmap (new, ~2 pp): conventions box, 3 structural results, TikZ pipeline, requirement->location table, status table,
  delicate points; \todo for (13)/combination numbers.
* sec_observable: affine param varsigma, vorticity varpi, projector perp; worldline/clock with label at the initial event (eq:worldline);
  new subsection sec:FW (FW construction, rotation, unitarity proof, eps^4 coefficient, l=2 example X-6) (M9); normalisation subsection
  moved here from sec_harmonic (label kept).
* sec_primordial (M1): Maldacena's gamma_ij is independent; Def. gauss now includes h_n(0)=0 and no primordial vectors (eq:notensor);
  regularity list (j0 tensor mode regular, NOT removed); eq:greens_split (int_0 = int_{eta_i} + homogeneous data, exact at finite
  eta_i); v2(0); C3 pole at q=k corrected (eq:C3pole) + \todo hook (A3); k->0 limit of K_Phi2 does not exist (eq:KPhi20_smallk) +
  \todo hook "zero-mode invariance" (V5); dilation notation lambda (polynomial, a_1,a_2), beta = e^{eps c}.
* sec_einstein: provenance removed; U(x); Pi projector; 0i factor; Table row for K^0; closed forms stated (App. closedforms).
* sec_harmonic: all rules rewritten in the n-hat convention (A-s, A-lambda [Dhat_{+-1} sign flips], B1, type D, type C, A-L, E1, E2/E3);
  closed form of xi(0) (eq:xi0, addition theorem) and analytic l-tail bound (eq:jbound, eq:tailbound) (M8); w_grad symmetric lensing
  formula (m3); Theta3c and Delta^(3c) defined once (m4); M_2 -> \bar K_2.
* sec_linear: pipeline steps (1)-(5) as worked example; S->s_o, A->A_SW, K(s)->kappa, M(k)->Delta_0; SPEC refs removed.
* sec_map2: provenance removed; sign remark rewritten (label rem:signslip kept because sec_assembly cites it); finite-eta_i clock and
  displacement (eq:dxo, eq:deta_o1, eq:deta_o2) with O(eta_i^2) terms displayed, eta_i->0 statement and what the numerics used (M2);
  m7 clause; table dagger; pointer to Fourier K_2 (sec:asm-K2).
* sec_map3: references full appendix; tadpole group as named formula (eq:tadgroup); worked contraction (emission-shift tadpole).
* sec_checks: rewritten without agent names; dilation = k=0 gauge test; plane wave vs comoving-gauge solver = finite-k gauge test (M7);
  \todo hooks for zero-mode invariance and mu->1; erratum dropped.
* app_lterms: generated, ALL terms (L1 3, L2/2 34, direct 76, L3/6 359) (M6).  app_closedforms (new, from C4; C4 complete): Phi2
  special-function skeleton (coefficient D), I_h complete (skeleton + 9 trig monomials, re-verified 2e-29 by code/W3_revision/closedforms_tex.py),
  limits, verification.  app_ledger: IDs renamed by check category (X, R, EF, CF, FO, SQ, H, NP, PW), no agent names; new R-8, H-11, H-12,
  CF-4; new Reproducibility appendix (app:repro) listing code directories.

## For the other writer (sec_assembly / sec_mono / sec_unresolved, not edited by W3)
* sec_assembly l.12-13 says Secs. typeA-typeC use e-hat and (-i)^l: now stale (rules are in the n-hat convention, i^l); delete the sentence.
* sec_assembly uses lambda for the dilation time rescaling (consistent with the new notation) and cites rem:signslip (kept).
* sec_assembly "observer label at the origin ... label = position as eta->0": this is the eta_i->0 limit of the label (O2); OK.
* sec_mono/assembly: agent names in \todo notes; "(13)(iii)" prefactor is now unambiguous: Delta^(3c) includes the factor 3 (eq:Delta3c).
* M_2 renamed \bar K_2 in (E2)/(E3).
