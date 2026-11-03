# Report W2b — writer: assembly (second-order half) and monochromatic sections
Files (cmb4/doc/): sec_assembly.tex, sec_mono.tex, app_k2channels.tex (new; \input in main.tex after app_kernels).
Generator for app_k2channels.tex: scratchpad gen_k2.py/gen_app.py (loads A2's symbolic_kernels.pkl, renames symbols, and
verifies symbolically that each compact channel formula quoted in sec_assembly equals A2's machine kernel: residual 0 for
all 15 channels).

## Milestone 1 — sec_assembly.tex (DONE, compiles)
- 8.1 K_1, Delta_l, <m1>, clock (lapse integral, closed form re-derived by hand), X_o (re-derived), emission shift.
- 8.2 K_2 bookkeeping (boxed), Fourier rules, J_a/U_a (direct-form integrand/boundary function), functionals F, D_perp, D_r
  (A2's sign), second-order pair kernels with projectors.
- 8.3 channels B/C/D with compact verified kernels (new compact forms: emission shift = -deta_e (J + DU), clock = deta_o (J + dU/deta),
  displacement = dx_o.grad U; endpoint 1x1 brackets are NOT opposite at obs/emi: the P_perp T_V T_V term has the same sign).
  Remark: obs channels have degree <=1 in n (clock, displacement, obs second-order fields => l<=1), obs 1x1 degree 2 (l<=2).
- 8.4 squeezed theorem (V2 derivation + A2 numerics), two paragraphs, IR formula for the (22) soft region
  (1/4)C22|soft = sigma_L^2 (2/pi) int k^2 dk P (k d_k Delta_l)^2 (derived here from V2's limit; consistent with V2's combined formula).
- 8.5 C22 general formula, radial functions per type (incl. Legendre radial functions F_l, D_perp l, D_r l of the functionals), checks.
- 8.6-8.8 \todo placeholders for A3 (K_3c, C13, IR structure).
