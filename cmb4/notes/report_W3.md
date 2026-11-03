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
