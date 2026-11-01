# Lead notes: CMB temperature power spectrum through O(eps^4), radiation-only GR universe
(Persistent continuity file. Update after every milestone. Read this first after any compaction.)

## Task (from user upload /root/.claude/uploads/.../46e6b817-message.txt)
Compute C_l = eps^2 C^(11) + eps^4 [ 1/4 C^(22) + 1/3 Re C^(13) ] for T_obs = T_* E_o/E_e, radiation-only flat GR,
Gaussian nonlinear primordial variable, observer = tagged fluid worldline at fixed proper time, tetrad FW-transported,
normalize by monopole then drop dipole. Keep analytic; declare unevaluated integrals; monochromatic k0 only at end.
Deliverable must be understandable by a physics postdoc. Also: research Anthropic's recent (Oct 2026) physics/math writeups.

## Working directories
- cmb4/notes/NOTES.md  (this file)
- cmb4/lit/            literature extracts
- cmb4/code/           sympy derivations & checks
- cmb4/doc/            final write-up (LaTeX + md)

## Status log
- [t0] Repo empty. sympy installing. Launching research agents (Anthropic writeups; literature).
- [t1] KEY RESULT (lead, verified at linear order by hand): T_obs(n) = T_o exp[-S(n)], S = int_e^o sigma_ee dl
  (shear of the fluid along the ray). Exact for radiation fluid. See SPEC.md §3. Plan: Poisson gauge evaluation
  of S to 3rd order; zeta ICs via zeta_nl = ln(1+delta_T)+1/2 ln(1-2Phi) at eta->0 (SPEC §5).
- [t1] Agents launched: research(Anthropic), literature. Next: C1 (Einstein eqs to 3rd order), C2 (observable S to
  3rd order + checks), C3 (harmonic machinery).
- [t2] Anthropic-writeups research done -> cmb4/lit/anthropic_writeups.md. Adopted lessons: (1) orchestrator does
  no heavy math, briefs = memos with fixed conventions + control cases; (2) separate generator/verifier agents,
  blind hostile referee for each claim; (3) CAS/numerics as arbiter; (4) state on disk, ledger of failed approaches;
  (5) treat "done" as a claim to refute.
## Ledger of failed approaches / decisions
- Decided: Poisson gauge for dynamics (literature-checkable, no secular growth); exact e^{-S} formula for observable.
- [t3] Literature extracted -> cmb4/lit/literature_extracts.md (33 papers in cmb4/lit/raw). Key: NO one-loop C_l exists
  in literature; (13) term never computed; only Mirbabayi-Zaldarriaga (MZ 1409.4777) keeps observer terms. MZ: constant+gradient
  zeta unobservable once x_o, v_o, proper time, sky monopole included (their 4.21-4.25) -> our dilation check (SPEC 9d) + gradient check.
  CPV 1109.1822: time shift eps = -eta zeta_L/3 in RD. Kohri-Terada RD kernel I_RD(u,v,x) with Si/Ci (their eq 22). PUB 1003.0481
  RD superhorizon ICs: Psi2 = 2/3 R2 + Psi^2 - Delta^-1[dPsi dPsi] + 3 Delta^-2 dd[dPsi dPsi]; Phi2 = Psi2 + 4Psi^2 + 3Delta^-1[..] - 9 Delta^-2[..].
- [t3] Lead's independent linear check (cmb4/code/lead/linear_check.py): T_Phi=-2(sin x - x cos x)/x^3, delta_T->zeta/3, v_pot->eta zeta/3 (v_i=d_i v),
  continuity+Euler verified. Use as reference against C1.
- [t3] Running: C1 (Einstein 3rd order), C2 (observable S), C3 (harmonics), V1 (referee of exact formula).
- [t4] Launched V2 (squeezed-limit predictions via MZ/CPV transformations; independent route) and N1 (exact planar 1+1
  numerical GR+radiation + ray tracing; nonperturbative arbiter; extracts Th2, Th3 maps for single plane wave).
  Agents in flight: C1, C2, C3, V1, V2, N1. Waiting. Next after results: assemble kernels (C1 x C2 -> K1,K2,K3c),
  feed to C3 formulas, run checks (dilation/gradient from V2; planar numerics from N1), then write-up + simplifier + cold-read referee.
- [t5] V1 referee: exact formula CONFIRMED (Kasner+radiation test to 1e-31; linear SW reproduced). Caveats: needs both P=rho/3 and
  rho∝T^4; e = last crossing of T=T_*; C^(22), Re C^(13) separately frame-dependent, only 1/4C22+1/3ReC13 invariant -> compute all
  orders in ONE frame (canonical triad). report_V1.md.
- [t6] C1 DONE (report_C1.md; code C1_einstein, loader kernels_api.load()). Engine verified (Bianchi to eps^3; exact dilation family
  Psi=(A^4-1)/2, Phi=(1-A^2B^2)/2, delta=A^-4-1, zeta=ln B solves eqs). 1st order = lead's. 2nd order: tensor source = Kohri-Terada
  (PASS); omega2 ALGEBRAIC via relativistic Kelvin curl[(1+delta_T)u_i]=0: K_om=(2i/k^2)[eta(q^2-p^2)TT+eta^2(q^2 T_q T'_p - p^2 T_p T'_q)];
  Phi2 superhorizon kernel = Pitrou-Uzan-Bernardeau 3.12 (PASS; BMR II 7.30 differs by their velocity quadrupole evaluation).
  3rd order: contracted kernels K3c (Q3,S00,S0z,SPhi) + 3rd-order IC amplitude C3 (squeezed limit = dilation family PASS).
  CAVEAT: scalar k=0 tadpoles (<Phi2>,<Psi2>) are residual zero-mode gauge; C1 used angle-averaged k->0 limit; must drop out of the
  observable (dilation check in assembly). Green's integrals left as integrals: G_Phi=[u1(eta')u2(eta)-u1(eta)u2(eta')]/W(eta'),
  u1=(sin x - x cos x)/x^3, u2=(cos x + x sin x)/x^3, W=-3sqrt3/(eta^4 k^3); G_h=(eta'/eta) sin(k(eta-eta'))/k.
- [t7] Launched V3 (referee of C1 2nd/3rd order, independent ADM re-derivation), C4 (closed-form Green's integrals, Kohri-Terada
  control), C5 (first-order Delta_l, C_l^(11) closed form, monochromatic, sky monopole <m1>).
  In flight: C2, C3, V2, N1, V3, C4, C5. Blocked on C2: assembly A1 (map kernels K1,K2,K3c), dilation check of tadpoles.
- [t8] C3 DONE (report_C3.md; code C3_harmonic). Type A rule (i n.k)^s -> (-d/dx)^s j_l; helicity rules; type B (product) C_l with
  both Wick pairings + w_grad for grad.grad; lensing control = Hu (62)-(63) EXACT; SPEC §8 prefactors 2 and 3 confirmed;
  type C (convolution kernel) C^(22) = 1/(4pi^4(2l+1)) int k1^2 k2^2 P P dmu sum_m |K~_lm|^2 in K-hat = z frame; (13) q-reduction.
  Monochromatic: C11 = 4 pi sigma^2 |Delta_l(k0)|^2; type C -> single mu integral; delicate mu->-1 (K->0) for 1/K^2 kernels.
  NOTE: C3 used ray x = -n chi (source at -n chi), SPEC §6 implies +... sign: a_lm -> (-1)^l; harmless if consistent. ASSEMBLY MUST
  USE ONE CONVENTION: ray x(eta) = x_o - n (eta_o - eta) (SPEC) i.e. source at -n chi. OK consistent with C3.
- [t9] Launched W1 (draft settled sections of doc: setup, exact observable, zeta ICs, Einstein solution, harmonics). Still blocked on C2 for assembly.
- [t10] N1 DONE (report_N1.md; output in code/N1_planar/output, summary.md). Exact formula T_obs/T_o = exp(-S) verified
  NONPERTURBATIVELY to 7e-15 (eps=0.1, Theta_rms=0.08). Reference maps for single plane wave zeta = eps cos(kx), k=1, t_*=2,
  t_o=8 (and 4), x_o=0 (and pi/2): a_20[Th1]=-2.18472036, a_40=+1.25524563; a_20[Th2]=-4.22249422, a_40=+1.37854398;
  a_20[Th3]=-12.890, a_40=-47.071 (unc. Th2 3e-8, Th3 2e-4). Clock shift t_o(±0.01)=8.02512265/7.97421965 for x_o=0.
  Gauge used: comoving-orthogonal with alpha = e^{-zeta}(AB^2)^{1/3} exactly. USE THESE AS ARBITER for assembled K1,K2,K3.
