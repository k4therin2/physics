# X1: comparison with Adam's independent attempt (saved by the lead from the agent's hand-back)

I couldn't write `/home/user/physics/cmb4/notes/report_X1_compare.md`: the harness blocks subagents from writing report files. Everything else is committed (4f627db) in `/home/user/physics/cmb4/code/X1_compare/`, which holds the scripts, logs and pkl files behind every number below. The parent can write the note from this summary. The two calculations agree everywhere I could test except one term: the O(k0^6) part of C13.

## Numeric cross-checks
**(a) First order.** I wrote Adam's T_l (his eq. 54, his p.25) directly from his definitions and compared it with our closed-form Δ_l. They are the same function with the same sign: agreement to ~1e-38 at 28 points (l = 2–6, k0η from 0.05 to 40). The worst case, 3e-21, is where our Si/Ci closed form itself loses digits.
- Symbolically, our small-k series of C_2^(11) through k0^4 and k0^6 minus his p.40 formula gives exactly 0.
- His leading T_l coefficients (eq. 80, l = 2–5) are reproduced exactly.
- Our C11 matches all 18 C11 entries in his App. A tables to his 3 printed digits.
- Conventions need no mapping: his t is conformal time, x = kη/√3 as ours, and his P(k) is identical to ours with ⟨ζ²⟩ = σ² held fixed. C11 = 4πσ²T_l² is the same as our 4πσ²Δ_l².

**(b) Superhorizon fourth order.** I reran our A2/A3 evaluators at k0η_o = 0.4, 0.2, 0.1, 0.05, 0.025, 0.0125 and Richardson-extrapolated:

| l | 1/4 C22/C11 | 1/3 Re C13/C11 | total/C11 | Adam |
|---|---|---|---|---|
| 2 | 5.666667 | 4.000000 | **9.666667** | 17/3, 4, **29/3** |
| 3 | 37.000000 | 9.000000 | **46.000000** | 37, 9, **46** |
| 4 | 252.257 | 16.000000 | 268.257 | no prediction |

- Both totals and both 22/13 splits agree to 7 digits, from two codes in different gauges. This is a strong independent confirmation of both C22 machineries and of the leading C13.
- The 13 part equals l²σ², exactly Adam's half of the observer dilation (D²/2 with D = 2l).
- **Next order (k0^6):** C11 agrees exactly, and 1/4 C22 agrees to 1e-4 (b22 = −3.2972 vs his −1187/360), which includes his endpoint and clock pieces.
- **Disagreement:** for 1/3 Re C13 we get b13 = −4.3796 (k0η_o)² against his −257/60 = −4.2833. We are lower by 0.0962 (k0η_o)²σ²·C11, about 2.2% of the k0^6 coefficient.
- **Not numerical noise.** The leading term is exact to 1e-7, the offset scales as x² from 0.4 down to 0.0125, and doubling the resolution leaves it unchanged.
- **Where it lives.** Repeating at η_* = 0.5, 1, 4, 6 gives the same offset (−0.0961 to −0.0963), while b22 still agrees each time. The offset does not depend on emission time, so it is an observer-local factor in Θ^(3): most likely the second-order clock and its K = 0 tadpole/homogeneous-lapse content, not the ray or emission.
- Our κ is invariant under the residual τ(η) to 1e-13, so the tadpole gauge choice is not the cause.
- Adam's two printed k0^6 forms agree with each other. I could not resolve which side is wrong. Next step: compute the clock contribution of the (q,−q) contraction analytically at O(k0^6) in both frameworks.

**(c) Late time vs subhorizon.** These are different limits. On our closed form:
- **Late time (x_* = 1 fixed, x_o up to 1e5):** Δ_l·x_o stays bounded, consistent with his O(ln x_o/x_o) bound.
- **Subhorizon (k0 → ∞ at fixed times):** (Δ_l − SW/Doppler form)·x_o² stays O(1–10). This confirms his Prop. 30, which is literally our `Delta_l_large_k` formula.

**(d) Inside the horizon.** I ran our evaluators on his trajectory x_* = x_o/3. We use η_i → 0; he uses x_I = x_o/6.
- **1/4 C22:** agrees to ≤ 0.13% at x_o = 0.15–4 (e.g. 8.31298 vs 8.324 at x_o = 4; 2.61082 vs 2.611 at x_o = 2).
- **Total C_2^(4):** differs by 1–5%:

| x_o | 0.3 | 0.6 | 1 | 2 | 3 | 4 |
|---|---|---|---|---|---|---|
| ours | 2.633e-3 | 0.01221 | −0.1473 | 1.991 | −1.600 | −32.34 |
| Adam | 2.64e-3 | 0.0127 | −0.141 | 2.010 | −1.57 | −32.2 |

- At x_o ≤ 1 the differences match the superhorizon 13 offset quantitatively (predicted −0.026 and −0.10 in the ratio at x_o = 0.3 and 0.6; observed −0.027 and −0.10).
- His qualitative statements are reproduced: 22 and 13 are comparable and partly cancel, and the sign pattern in x_o is the same.
- Our default A3 resolution (16, 8, 24) was off by 1.8% at x_o = 4; the values above use converged (32, 16, 40) runs.
- Our published point (k0η_* = 2, k0η_o = 8) is not on his trajectory, so it can't be compared directly.

## Setup, organisation, depth
- **Gaussian variable:** the same object. His ζ = (1/6) ln det γ in the temperature frame (Def. 6, p.7); ours ζ_nl (§5). The dilation is fixed by ⟨g⟩ = 0 for him and ζ^(2)(0) = 0 for us. Both keep the forced K = 0 mode.
- **Gauges differ:** his temperature frame (comoving, uniform-temperature slicing) vs our Poisson gauge.
- **Adam has the exact shear formula.** His Lemma 2 (p.5), T_obs = T_fluid(o) e^{−I}, is our Thm 4.3 in coordinates. Ours is covariant and holds with vorticity; his assumes an irrotational fluid. He has no Ellis/Dunsby attribution.
- His "triangular scheme", adjoints and nine-moment reduction reorganise the integrals. By his own account (§11) he never evaluates C13 at finite k0; his finite-k0 numbers are labelled "not evidence".
- **Adam has, we lack:**
  - exact-in-ζ superhorizon theorems (curvature quadrupole and octupole) with a Gaussian resummation;
  - spectrum bounds [8, 29/3]m0 and [18, 46]m0;
  - a closed-form O(k0^6) C_2;
  - the log-map split (Thm 18) and the observer-dilation decomposition (Prop 19);
  - T_l rebuilt from its superhorizon series;
  - finite t_I throughout.
- **We have, he lacks:**
  - an evaluated C13 and full spectrum at a finite-k0 point;
  - the planar non-perturbative GR end-to-end check at orders 1–3;
  - an analytic collinear-residue cancellation;
  - the zero-mode classification;
  - a Monte Carlo of C22;
  - results up to l = 14.

## What each could borrow
- **We from Adam:**
  - his k0^6 formula as an analytic test of κ, to settle the discrepancy;
  - the observer-dilation reading of 29/3 = 8 + 5/3, which should replace our "9.6σ², expected for non-scale-invariant" sentence in §12.7;
  - the log-map split;
  - finite η_i in the evaluators.
- **Adam from us:**
  - a covariant proof of Lemma 2 and the Dunsby attribution;
  - a non-perturbative end-to-end check;
  - the residue analysis;
  - the zero-mode classification;
  - an evaluated, converged C13.

Files are in `/home/user/physics/cmb4/code/X1_compare/`:
- `transfer_compare.py`, `series_compare.py`, `adam_table_c11.py`, `latetime2.py`: first order and limits
- `superhorizon_x1.py`, `sh_analyse.py`, `sh_ratio.py`, `adam_sh_pred.py`, `adam_b_of_e.py`, `adam_consistency.py`: superhorizon fourth order
- `adam_config.py` with `ad22.log` and `ad13_*.log`: his App. A trajectory