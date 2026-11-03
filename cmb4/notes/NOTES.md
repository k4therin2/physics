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
- [t11] C2 DONE (report_C2.md; L_terms.pkl term dicts; shear form + direct form). Exact formula verified symbolically to 3rd order
  (both routes). L1 = SW formula. Dilation: full map zero through 3rd order incl. clock + emission (a2 tadpole free -> drops out).
  Boubekeur et al structure reproduced. CONVENTION FIX: background ray xbar(eta) = x_o + n (eta_o - eta), e = -n, n = observation
  direction (C2). C3/C5 used source at -n chi => their a_lm differ by (-1)^l: harmless for C_l but ASSEMBLY must use C2's convention
  and re-derive Rayleigh rule signs: e^{+ik.n chi} = 4pi sum i^l j_l Y_lm(n) Y*_lm(k).
- [t11] NEXT: A1 = plane-wave end-to-end check vs N1 (Th2, Th3 maps); A2 = Fourier kernels K2, K3c + C_l assembly + monochromatic.
- [t12] Launched A1 (plane-wave vs N1), A2 (K2 + C22 + mono), A3 (K3c + C13 + mono). In flight also: V2, V3, C4, C5, W1. Blocked: W2 (results sections), cold-read referee, simplifier.
- [t13] W1 DONE: doc skeleton cmb4/doc/main.tex with sec_setup, sec_observable, sec_primordial, sec_einstein, sec_harmonic, app_kernels, app_ledger, refs; placeholders sec_map2/map3/assembly/checks/mono/unresolved. No pdflatex in container (try pip install? or tectonic). NOTE W1: C1 report omitted overall sign of invariant K_Q; appendix uses machine form.
- [2026-11-03] All 7 in-flight agents (A1,A2,A3,V2,V3,C4,C5) died on usage-credit limit on 11-01. User asked to use Opus 5.5 for
  subagents. Briefs saved persistently in cmb4/notes/briefs/. All seven RELAUNCHED on Opus 5.5 (resume from partial files).
  Blocked after them: W2 (results sections sec_map2/map3/assembly/checks/mono/unresolved), cold-read referee, simplifier.
  pdflatex now installed; doc compiles (29 pp).
- [11-03] V3 DONE: all 5 C1 claims CONFIRMED (K_h = 4 f_RD Kohri-Terada to 5e-29; omega2 = Lu et al eq 20; Phi2(0) = PUB 3.12; C3
  re-derived independently to 30 digits). Erratum: third_order_ic.py line 72 coefficient 3/10 should be 1/10 (affects only O(1) value of
  K3c[SPhi](0), not C3). CAVEAT: tadpole zero-mode family Phi2 = -Psi2/2 - 4/9; C3 shifts by -2 dPsi2 per image; assembly must use the
  SAME tadpoles with matching clock/normalisation (dilation check). PUB (A.6) trace source garbled in paper (literature note).
- [11-03] C5 DONE (report_C5.md; code C5_linear_Cl; linear_cl.tex to \input into doc). Delta_l = A(x_*) j_l(w) - k T_V(x_*) j_l'(w)
  - d_l0 A(x_o) + (1/3) d_l1 k T_V(x_o) + I_l(k); A = T_dT + T_Phi = -cos x + 2 j_0(x) - 4 j_1(x)/x; k T_V = sqrt3 (sin x - 2 j_1(x));
  ISW I_l closed form for every l (Ci/Si with a± = 1 ± 1/sqrt3; coefficients as residues; verified 1e-169 / 5e-18). RD SW amplitude
  -zeta/3. Monochromatic C11 = 4 pi sigma^2 Delta_l(k0)^2. Sky monopole M = Delta_0; <<m1>^2> = C_0/4pi; <<m1> a_lm>=0 for l>=1.
  Scale-invariant plateau l(l+1)C_l/(2pi A_s): SW alone 1/9; with ISW 0.377 (l=2) -> 0.914 (l->inf).
- [11-03] *** A1 DONE: FULL ANALYTIC CHAIN MATCHES EXACT SOLVER N1 AT ORDERS 1,2,3 in 4 configs (eta_o=8,4; x_o=0,pi/2) within N1
  noise (a_20[Th3] diff 7e-8; residual Theta_N1 - series scales as eps^4). Tadpole independence confirmed (2e-14). Bugs found:
  (1) C2 DOC sign: k1(eta) = k1_init - int_eta^{eta_o} G1 (minus, not plus); L-term lists correct. A2/A3 informed. (2) crashed-run bug
  in A1's own code fixed. (3) SPEC §6 wording: observer label = position as eta->0: x_o(eta) = x_label + int_0^eta v. N1's eta_o shift
  is in N1's gauge, not Poisson. Report: report_A1.md; code A1_planewave (pw_symbolic, pw_numeric, pw_map).
- [11-03] A2 DONE (report_A2.md; symbolic_kernels.pkl/.tex: all K_2 channels, B/C/D types). Plane-wave vs N1 PASS (Th2 a_l0 to 1e-8);
  shear vs direct form at generic non-collinear momenta 6e-15 (tests vector/tensor kernels at K!=0). CORRECTION TO SPEC 9d: constant
  mode limit K_2(k1->0,k2) = -k.d_k K_1 (DILATION of first-order kernel), not zero; O(k1) = V2's SCT -S_{k1}K_1. IR of the loop:
  dilation-response structure (1/2)(k d_k)^2 of the C11 integrand -> vanishes for scale-invariant P after k-integration; IR-sensitive
  for tilted P. Clock shift delta eta_o1 = -(1/eta_o) int_0^{eta_o} eta Psi1 (lapse, not delta_T). Observer displacement X_o =
  (3/k^2)(2 sin x_o/x_o - cos x_o - 1). mu->-1 (K->0): Phi2,Psi2,V2 ~ K^-2, omega2+v2T ~ K; K_{2,lm} finite for l>=2 (dipole log).
  Monochromatic numbers (k0=1, eta_*=2, eta_o=8; sigma units): l=2: C11=0.95460, 1/4 C22=14.2575; l=4: 0.17507, 1.86592; l=10:
  9.68e-6, 4.0085e-2. A3 informed of the corrected IR criterion.
- [11-03] V2 DONE (independent route; report_V2.md): K_2(k1->0,k2;n) = -D K_1 - S_{k1} K_1 + O(k1^2), D = k.d_k, S_p = (k.d_k)(p.d_k)
  - (1/2)(p.k) nabla_k^2; agrees with A2's numerics. Third order: K_3(q,-q,k) = D^2 K_1 + O(q^2); K_{m+1}(0..0,k) = (-D)^m K_1.
  IR: [C22/4 + Re C13/3]_IR = sigma_L^2 (1/2) int |K_lm|^2 (3+D)^2 P -> total derivative, vanishes for P ∝ k^-3. "Gaussian zeta breaks
  Maldacena's consistency relation unless n_s = 1" (comment to include). A3 informed.
- [11-03] W2a DONE (sec_map2, sec_map3, sec_checks, app_lterms, ledger). Lead fixed sec_observable: label at eta->0; ledger ref V1-6.
- [11-03] W2b DONE: sec_assembly (9.1-9.5, 15 K_2 channels verified vs A2 machine kernels), sec_mono 2nd-order half, app_k2channels. Placeholders 9.6-9.8 + mono (13)/combination for A3. Generator scripts copied to doc/.
- [11-03] R1 cold-read DONE (report_R1.md): Secs 1-5, linear, maps re-derived OK. BLOCKING: B1 C_3 has simple pole at mu=1 (q=k, p->0):
  monochromatic (13) mu-integral log-divergent unless 1/p^2 poles cancel in full K_3c (forwarded to A3). B2 tadpole convention swaps
  limit/angle-average: k->0 limit of K_Phi2(0) ~ q^2(3mu^2-1)/k^2 (constant anisotropic shear residual gauge) -> V5 launched to prove
  zero-mode invariance analytically. MAJOR: M1 declare h_n(0)=0 (regular j_0 mode not excluded by regularity); M2 keep eta_i terms in
  clock/label; state fluid-rest-frame emission; M3 stale status statements; M4 unify direction convention (n vs e); M5 remove project-log
  language (agent names, SPEC refs, pkl); M6 display Theta^(3) fully (appendix) and K_2 Fourier form; M7 present N1 comparison as gauge
  test; M8 analytic l-cutoff; M9 FW tetrad proof; M10 define <>_c. Proposed reorder + 2-page roadmap in report.
- [11-03] Launched V5 (zero-mode gauge invariance proof, M1 tensor-mode statement), W3 (revision per R1 M-items, reorder, roadmap, full L3 appendix, eta_i restored). In flight: A3, C4, V4, V5, W3. Blocked: W2c (third-order sections, unresolved), R2 final referee, final summary.
- [11-03] A3 DONE (report_A3.md; code A3_third_order; A3_channels.tex). kappa := <dTh3/dzeta_k> = (1/2) int P coef_abc, K_3c = kappa/3;
  C13_l = (2/pi) int k^2 P Delta_l Delta^kappa*_l; mono: 4 pi sigma^4 Delta_l Delta^kappa*_l. Shear vs direct L3 with C1's contracted
  fields: 2e-11. Tadpoles: kappa invariant under residual K=0 gauge eta->eta+tau(eta) to 1e-13 (B2 resolved in gauge sense; only
  qhat-averaged tadpole can enter). B1: six channels have 1/(1-mu) poles (residues v3 +0.0117, Dx2 time delay -0.2805, clock +0.2393,
  emission +0.0295, local 2+1 +0.0247, radial k2 -0.0247; sum ~1e-5 = 0 within accuracy) -> monochromatic C13 FINITE; integrand finite
  at mu->±1. IR: qhat-avg integrand -> 3(k d_k)^2 K1 = V2's D^2 K1; combined IR = (2/pi) int dlnk (1/2) Delta_l^2 d^2_{lnk}(k^3 P) -> 0 for n_s=1.
  UV: numerically O(q^0) oscillatory up to q=4k (not proven). Mono numbers (k0=1,eta_*=2,eta_o=8, sigma^4): 1/3ReC13: l=2 -5.39801,
  3 -23.5152, 4 -9.06394, 6 -1.586e-2, 10 -1.463e-3; 1/4C22+1/3ReC13: l=2 +8.85953, 3 -20.6370, 4 -7.19802, 6 +0.506198, 10 +3.8622e-2.
  TRIAD: A3 uses shear form + canonical GLOBAL (rigid) triad; Theta2 is triad-independent (first-order direction map has no rotation
  for scalar perturbations; O(eps^2) direction corrections enter S only at O(eps^3)); a rigid O(eps^2) rotation drops out of Re C13
  identically (a*(-i L a) imaginary) -> combination consistent. Not done: eta_o=4 for C13; closed-form Delta^kappa_l; UV proof.
- [11-03] W3 DONE (revision per R1: reorder, roadmap, intro, eta_i restored, app_closedforms from C4, full L3 appendix, FW lemma; 85 pp). Launched W2c (third-order sections, unresolved, roadmap fill) and V6 (analytic collinear residue cancellation). In flight: W2c, V4, V5, V6, (C4 finishing). Blocked: R2 final cold read, final user summary.
- [11-03] W2c DONE: manuscript complete (102 pp) except HOOK-V4 (MC C22) and HOOK-V5 (zero-mode proof). Lead confirms W2c's out-of-ownership edits (abstract, roadmap status rows/bullets, ledger closing). Next: R2 final cold read after V4/V5/V6 land; then final summary.
- [11-03] V6 DONE: collinear 1/(1-mu) residue cancellation PROVEN analytically (report_V6.md). Long second-order mode pole is the
  K=0 time reparametrisation eta -> eta + tau, (eta^2 tau)' = -eta^2 A_Q, A_Q = s^2[2T^2 + (T + eta T')^2]; pole parts Phi2 = -tau/eta,
  Psi2 = tau' + tau/eta, delta2 = -4 tau/eta; C_3 pole = same shift. Residues (tt = tau/2k^2): clock 3tt I1|_o, emission -3tt I1|_e,
  time delay + v3 = -3 int d/deta[tt I1] (total derivative) cancels clock+emission; local 2+1 = -radial k2. Sum zero for every n.
  Numerically: residues to 9e-11, sum <= 8e-10 (2.3e-9 relative). No (1-mu)^{-1/2} term; finite limits F_l match A3 to 1e-4.
  -> B1 RESOLVED ANALYTICALLY. Manuscript must be updated: "numerical only" -> proven (W4 final pass).
- [11-03] V5 DONE (report saved by lead): residual k=0 freedom = time reparametrisation eta -> eta + F(eta) ONLY (dilation/shear/
  rotation/translation fixed by zeta^(2)=0, h(0)=0, label). Lead's constant-shear hypothesis REFUTED: the 1/k^2(3mu^2-1) term is a
  quadratic/dipolar time-reparametrisation mode (satisfies Psi = -(G'+G/eta), delta = 4Phi, V = eta Phi); angle average at fixed k
  removes all negative powers of k => tadpole limit exists = C1's values. Map invariant to 2nd order (rigid rotation + total
  derivative). B2 RESOLVED. M1: declare h_n(0)=0; regular <=> omega = v^T = 0 for k!=0.
- [11-03] R2 DONE (report saved by lead): READY AFTER EDITS, no blocking. M1 FW frame: no O(eps) rotation for scalar modes => C22 and
  Re C13 separately invariant for this observer (restate); M2 stale/inconsistent status statements; M3 residue table extraction
  artefact (now superseded by V6 analytic proof); M4 overfull boxes; M5 triad caveat for printed L3 (n=z split); M6 monochromatic
  result thin (box formula; k0 eta_o -> 0 check). Next: W4 final editor pass (R2 + V5 + V6 + V4 hook).
- [11-03] V4 DONE: MC (180k realisations) agrees with A2's C22 within 1.1 sigma at 0.3-0.7% for l=1..11 (l=12 +2 sigma); rules out wrong Wick factor (50%) and missing -2<m1>m1 (6 sigma at l=3). report_V4.md. Added *.pkl ignore for V4.
- [11-03] C4 DONE (report_C4.md; second_order_closed.py; C4_closed_forms.tex 16pp): ALL second-order fields in closed form (Si, Cin, elementary; D = (k^2+p^2+q^2)R/(8k^2p^3q^3); w_i = (k+q+p, k+q-p, k-q+p, q+p-k)/sqrt3); full Einstein-eq residual 5.5e-50; quadrature agreement 1e-31..1e-44; small-eta series; eta->inf asymptotics correct (no sign error). W3 already created app_closedforms from C4.
