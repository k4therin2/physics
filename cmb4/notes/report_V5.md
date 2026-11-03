# V5 report: zero-mode (k = 0) freedoms of second-order Poisson gauge (referee items B2 and M1)
(Saved by the lead from V5's hand-back; code and logs committed in cmb4/code/V5_zeromodes/, commit d13da15.)

## Step 1: residual freedom (gauge.py, zeta_check.py)
Each candidate applied as an exact coordinate change to FRW + generic first-order fields; all satisfy d.omega=0, h TT.
- Time reparametrisation eta -> eta + F(eta), F = O(eta): Psi = F' + F/eta, Phi = -F/eta, delta = -4F/eta. Leaves zeta_nl, h(0),
  observer label unchanged: GENUINE residual freedom.
- Dilation, constant shear, constant rotation, constant translation: each changes zeta^(2) (e.g. by x.grad zeta_1), or h(0), or the
  label; fixed by the initial-data conditions: NOT free.
- Time-dependent translations c(eta), rotations w(eta) with c(0)=w(0)=0: pure gauge. eta -> eta + c0: Psi ∝ c0/eta singular, excluded.
LEAD'S HYPOTHESIS REFUTED: the 1/k^2 (3mu^2-1) term of K_Phi2(0) is NOT a constant shear (kphi2_limit.py):
- K_Q^(-2) = -(3mu^2-1) q^2 K_h(q,q)/4 (same function as the tensor source); K_SPhi^(-2) = (K_Q^(-2))'/eta; at eta=0, Q_d = 3 Phi_d.
- The 1/k^2 and 1/k parts of (Phi2, Psi2, delta2, V2) satisfy exactly the time-reparametrisation relations Psi = -(G' + G/eta),
  delta = 4 Phi, V = eta Phi: a quadratic/dipolar time-reparametrisation mode; its induced homogeneous shear vanishes at eta=0.
- The q-angle average at fixed k removes every negative power of k: the tadpole limit exists and equals C1's angle-averaged value
  (-16/9 for K_Q, -8/9 for K_Phi2(0), 32/3 for eta^2 K_S00); the 1/k part of omega2 averages to zero. B2 RESOLVED.

## Step 2: map Theta to second order unchanged (run_all.py, zm_invariance.py)
C2's line-of-sight code rebuilt with first-order omega, h allowed (reproduces C2's L1, L2 exactly). Each generator superposed on
arbitrary first-order fields with observer event, emission surface, triad and ray functionals transformed. Result: Delta L1 = 0;
Delta L2 = rigid sky rotation + exact total derivative int D(F.integrand) whose boundary terms cancel. Holds for: time rescaling and
general cubic F(eta) (no rotation); dilation and c(eta) (no change); constant shear (rotation omega = (-S23, S13, -S12)/2 = Gram-Schmidt
prediction); w(eta) (rotation w(eta_o)); quadratic mode (rotation from g(eta_o)C); dipolar mode.
Controls: h1 = -S without cross terms fails (as it should); two deliberately wrong ray functionals fail; C2's dilation check still zero
through third order. Bugs found: sympy solve silently dropped equations (fixed; residuals checked explicitly); C2's documented sign
k1 = k1_init + int G1 is wrong (minus used).

## Step 3: physical tadpole combination
zeta^(2)(k=0) = 0 with the 00 constraint delta2 = -2 Psi2 + 32c^2/9 at eta->0 fixes the dilation direction b2. Remaining: the
time-reparametrisation family F(eta), pure gauge; constant member (dPsi2, dPhi2, ddelta2) = (2a, -a, -4a). C1's two tadpole members
differ by a = 16/27 along this family: the choice is pure gauge. Third order compensates (reparametrised solution solves the field
equations, keeps zeta^(3)=0 and h(0)=0; regular solution unique). Theta3 invariance follows from the same multilinear identity (no
separate third-order run).

## Step 4: referee M1 (m1_modes.py)
Tensors: j_0(k eta) analytic and regular (indicial roots 0, -1): regularity does not exclude it; h_n(eta->0)=0 must be declared as
initial data (also removes the constant-shear zero mode at k=0). Vectors k!=0: only homogeneous mode omega = -8 C_v/(k^2 eta^2),
singular, carried by primordial vorticity: regular <=> omega = v^T = 0. Vectors k=0: only omega = -v = c'(eta), pure gauge.
