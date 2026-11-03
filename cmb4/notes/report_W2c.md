# Report W2c — completion of the manuscript: (13) assembly, monochromatic (13), unresolved steps

Owns: sec_assembly, sec_mono, sec_unresolved, app_k2channels, app_k3channels (new), roadmap (13) todo, B1/B2 hooks in
sec_primordial/sec_checks, ledger entries for A3/V4/V5.

## Status at start (2026-11-03)
* report_V5.md: ABSENT (zero-mode invariance proof still running) -> hook to be left if still absent at the end.
* report_V4.md: only the Method section (MC of C22 running; mc_n20.* being written) -> re-check before finishing.

## Milestone 1 — app_k3channels.tex generator  [DONE]
* doc/gen_k3.py: typographic conversion of code/A3_third_order/A3_channels.tex (sympy function notation -> v^x_1, h_{2,xz},
  k^x_1, Delta x^z_2, delta eta_{o1}, ...), no algebra; header rewritten in prose (contraction rule, Theta3 assembly).

## Milestone 2 — sec_assembly, sec_mono  [DONE, not yet compiled]
* sec_assembly: stale e-hat/(-i)^l sentence deleted; agent/code names stripped; new 9.6 (kappa definition, factors, eq:asm-clock3
  second-order observer objects, classes i/ii/iii with compact kernels eq:asm-classi..iii, kappa_norm eq:asm-kappanorm, triad Lemma
  lem:asm-triad), 9.7 (C13 split, validation list, tadpole convention = B2 resolution in the gauge sense, pending analytic proof),
  9.8 (Theorem thm:asm-squeezed3 + table, eq:asm-IRtot, eq:asm-IRns: soft = (1/2)(n_s-1)^2 sigma_L^2 C11, UV table, numerical only).
* sec_mono: (13) collapse eq:mono-C13, residue table l=0..6, F_l, structural explanation, numerical-only caveat, evaluation settings,
  Delta^kappa_l; Table tab:mono-comb (l=2..12), remarks; hook sec:mono-mc for V4.

## Milestone 3 — unresolved section, roadmap/hooks, ledger, compile  [DONE]
* sec_unresolved: 10 numbered items (residues, UV, IR tilted, closed forms vs declared integrals, E2/E3/T(ii), gauge, end-to-end coverage,
  tadpole convention [HOOK-V5], limit exchanges, model limitations).
* roadmap: (13)/combination todo replaced by l=2,4,10 numbers; status rows (K3c, kappa/C13, zero-mode) and two delicate-point bullets updated.
* sec_primordial / sec_checks: B1/B2 todos replaced by resolution text; analytic zero-mode proof left as \pending with %% HOOK-V5.
* app_ledger: new category L; entries R-9..R-12, EF-22, SQ-2, SQ-3, H-13, L-1..L-4 (L-4 = V4 hook); repro list split A2/A3/V4; pending
  paragraph updated.  app_k2channels intro: code/pkl names removed.  main.tex: \input{app_k3channels}; abstract status sentence updated
  (minimal edit, W3-owned file -- flagged).
* pdflatex x2: 0 errors, 0 undefined references; no overfull boxes in the new material.
