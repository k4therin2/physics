# Referee report R2 (final cold read)
(Saved by the lead from R2's hand-back. Checks: cmb4/code/R2_referee/check_r2.py, check_kq.py.)

Scope: fresh compile (102 pp, no undefined refs); every section and appendix; problem statement; report_R1.md. HOOK-V4/V5 not counted.
VERDICT: ready after the listed edits; no blocking items.
R1 items: all addressed (B1, B2 settled, analytic proofs listed as open; M1, M2, M4-M6, M8-M10 done; M3 done apart from stale
leftovers; M7 argued; general gauge test declared not done).

Re-derived and correct: normalisation series m2, m3, Theta2, Theta3; K_Q(0), K_Phi2(0) forms; zeta^(2)(0)=0; small-k form (5.13);
C_3 pole at q=k: 16/(9(1-mu)) - 56/3; C_3(q->0) = -32/9; finite-eta_i clock delta eta_o1, delta eta_o2 (eta_init and 1/u^0 at 2nd
order); delta x_o2; delta eta_e1,2; closed forms delta eta_o(k), X_o(k); dilation clock check; L2 endpoint bookkeeping; worked tadpole
contraction; kappa = 3 int P K_3 = 1/2 int P coef_abc; class (i)/(ii) rules; kappa_norm (4 terms); C13 = (2/pi) int k^2 P Delta Delta^kappa*;
(E2)/(E3) prefactors and rule (A-L) by independent Wick count; K2(0,k) = -D K1, K3(0,0,k) = D^2 K1; soft (22)+(13) = 1/2 (n_s-1)^2
sigma_L^2 C11; residual tau-gauge shifts; triad lemma; K->0 power counting; 1/3 Re C13 = 4 pi Delta_l Delta^kappa_l/3 reproduces tables.

MAJOR
M1. Frame-dependence claim overstated (Sec 4.3(iii), 4.5, roadmap row, 12.6, ledger X-6): for first-order scalar modes the FW frame
    does not rotate relative to the canonical triad at O(eps) (vorticity, curl omega, Thomas precession all O(eps^2)); R = 1 + eps^2 J2/2
    rigid => C22 and Re C13 SEPARATELY invariant for this observer. Fix: one-line J1 = 0 argument; merge lemma into 4.5; restate.
M2. Inconsistent status statements: Sec 6 "declared convention" vs 5.4 "not a convention"; 10.8 delicate point (2) "pending", (3) stale;
    Sec 14 item 9(e) "justified" vs item 8 "proof missing"; Sec 11 intro "when available, K_3c"; zero-mode discussion in six places.
    Fix: one authoritative statement in 11.6 with pointers; delete 10.8 "Delicate points".
M3. Residue table: "sum" row equals (1-mu) F_l at 1-mu = 1e-6 (residues extracted including the finite part); quoted "<= 1.5e-4" is an
    artefact; l=0 rows add to -9e-5 not -5e-5. Fix: subtract finite part; update 5.3, 12.5, 13.6, 14.1, ledger L-1. Structural
    explanation: at fixed p-hat the O(p^-2) piece is spatially constant => K=0 homogeneous solution = tau-mode (+ dilation).
M4. 47 overfull boxes (up to 612 pt): App A tadpoles, invariant C_3, E2/E3, monochromatic rules, (6.1), (5.11), Sec 7. Fix with
    (1/denominator) x [sum] or multline.
M5. "Theta^(3) fully explicit" has an unstated triad caveat: App C L3 lists and k2,init display are in the n=z split; per-direction
    triad changes kappa by 0.3%. Fix: state k2,init must come from the global triad with general n^a; print/describe general-n form.
M6. Monochromatic O(eps^4) result thin: (13) at one point; no k0-dependence or asymptotic limit; no displayed exact expression.
    Fix: boxed result in named objects; "depends only on k0 eta_*, k0 eta_o; sigma^2 fixed; no k0 scan"; ideally check the O(eps^4)
    coefficient vanishes as k0 eta_o -> 0.
MINOR: abstract UV "established numerically" overstated (probed to q=4k); add IR result; C22 two points, C13 one. (E2)/(E3) unused:
    demote or test against kappa_norm channels (+3.96 at l=2). Sec 11.7 definite sign only for scale-invariant. Sec 12.6(b) bound at
    l=4 is 0.024 not 0.03. Inconsistent numbers: ISW 5e-18 vs 2e-18; large-k 300/3000 vs 200/2000; solver 7.5e-15/1.4e-14 vs 7e-15.
    Only one \unevaluated marker. Phi2 elementary part and K_3c unprinted. Sec 10.6 "assembled in Sec 11" not displayed. Notation
    clashes: h2, U, E_a, mu, calligraphic D, c, P, bold omega, R. Log leftovers: EF-3, EF-12, EF-6 "A,B", FO-4 "M(k0)", X-1, "Working draft".
Completeness: met except partly: evaluated-vs-unevaluated (unprinted closed forms, result integral not displayed), monochromatic
    k0-dependence/limits (M6), limit exchanges (mu->1 and soft region numerical, UV unjustified; all declared).
Clarity edits (3-4 h): (1) roadmap "Main results" block; (2) merge 4.5 with triad lemma; (3) single zero-mode discussion; move 5.3
    collinear "Resolution" to 12.5; (4) slim Sec 11 (numbers to ledger; split 2nd/3rd order); (5) move Sec 12 "Numerical reliability"
    to ledger; (6) delete 10.8; (7) notation clashes and log leftovers.
