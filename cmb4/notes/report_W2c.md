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
