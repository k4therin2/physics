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
