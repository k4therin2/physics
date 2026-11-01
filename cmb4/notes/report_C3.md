# Report C3 — harmonic-space machinery for C_l^(11), C_l^(22), C_l^(13)
(agent C3; code in `cmb4/code/C3_harmonic/`; run `python3 run_all_tests.py` there)

Status log
- [M1] Conventions fixed; type-A rules (scalar monomials, helicity ±1, ±2) derived and verified to 1e-14; C^(11) normalisation verified.
- [M2] Type-B product formula (both Wick pairings, monopole×map terms) verified exactly (1e-12, real-space route) and by Monte Carlo;
  lensing control case: Hu (62)-(63) reproduced from SPEC (22)/(13) bookkeeping (exact algebra + MC); SPEC §8 prefactors confirmed;
  type-D (nested/unequal-time) reduction and its Fourier↔angular-gradient identity verified (3e-12).

## 0. Conventions used throughout (SPEC §1, §5, §7, §8)
- Fourier: f(x) = ∫ d³k/(2π)³ e^{ik·x} f_k, ⟨ζ_k ζ_k'⟩ = (2π)³ δ³(k+k') P(k). Hence ⟨ζ²(x)⟩ = ∫ d³k/(2π)³ P = (1/2π²)∫k²dk P(k).
- Observer at x = 0; background ray x(η) = −n χ, χ = η_o − η, n = observation direction (as in the task text; this is the sign
  under which the other agents' kernels carry e^{−ik·nχ}). NOTE: SPEC §6 (p̂ = E(û − n^a e_a), photon moves along −n) puts the
  source at x = +nχ. The two choices differ only by n → −n, i.e. a_lm → (−1)^l a_lm for EVERY map; all C_l and all
  cross-spectra below are unchanged provided one convention is used consistently. I use x = −nχ everywhere.
- Y_lm with Condon–Shortley phase (scipy.sph_harm_y = sympy Ynm, checked). Gaunt G^{l l1 l2}_{m m1 m2} := ∫dΩ Y*_lm Y_{l1m1} Y_{l2m2}
  = (−1)^m √((2l+1)(2l1+1)(2l2+1)/4π) (l l1 l2;0 0 0)(l l1 l2; −m m1 m2).
- a_lm = ∫dΩ Y*_lm Θ, C_l^{XY} = (1/(2l+1)) Σ_m ⟨a^X_lm a^{Y*}_lm⟩. Θ = εTh1 + ε²Th2/2 + ε³Th3/6,
  C_l = ε² C^(11) + ε⁴[¼C^(22) + ⅓ Re C^(13)].
- Rayleigh (x = −nχ): e^{−ik·nχ} = 4π Σ_lm (−i)^l j_l(kχ) Y_lm(n) Y*_lm(k̂).

## 1. Type A (linear line-of-sight maps) — `typeA.py`
Definition. Θ_A(n) = ∫ d³k/(2π)³ ζ_k ∫_{η_*}^{η_o} dη Σ_s g_s(k,η) (i n·k̂)^s e^{−ik·nχ(η)}  (+ helicity pieces, below).
Monomials are written as (i n·k̂)^s because i k̂ = ∇/k: a real map has real g_s (the Doppler term n·v has v_k ∝ i k̂ ζ_k).
δ(η−η_1) pieces (emission/observer endpoint terms) are included by replacing ∫dη → evaluation at η_1.

RULE A-s (scalar):  (i n·k̂)^s e^{−ik·nχ} → (−d/dx)^s j_l(x), x = kχ.  Equivalently (n·k̂) → i j_l', (n·k̂)² → −j_l'', (n·k̂)³ → −i j_l'''.
   a_lm = 4π (−i)^l ∫ d³k/(2π)³ ζ_k Y*_lm(k̂) Δ_l(k),    Δ_l(k) = ∫dη Σ_s g_s(k,η) (−1)^s j_l^{(s)}(kχ(η)).
   [With x = +nχ the prefactor is 4π i^l and (n·k̂) → −i j_l'.]

RULE A-λ (helicity). Frame k̂ = ẑ, e^{(±1)} = (x̂ ± iŷ)/√2, e^{(±2)}_ij = e^{(±1)}_i e^{(±1)}_j, so
   Q_{±1} = n·e^{(±1)} = sinθ e^{±iφ}/√2,   Q_{±2} = e^{(±2)}_ij n^i n^j = sin²θ e^{±2iφ}/2:
   Q_λ e^{−ik·nχ} = Σ_l (−i)^l √(4π(2l+1)) Δ̂^{(λ)}_l(x) Y_{lλ}(n)   with
   Δ̂^{(±1)}_l(x) = ∓ i √(l(l+1)/2) j_l(x)/x,      Δ̂^{(±2)}_l(x) = −½ √((l+2)!/(l−2)!) j_l(x)/x².
   (Derivation: ∫dμ (1−μ²)P_l'(μ) e^{−ixμ} = 2(−i)^{l−1} l(l+1) j_l/x via (1−μ²)P_l' = l(l+1)(P_{l−1}−P_{l+1})/(2l+1) and
   j_{l−1}+j_{l+1} = (2l+1)j_l/x, together with Y_{l1} = −√((2l+1)/(4πl(l+1))) sinθ P_l' e^{iφ}; the ±2 case likewise.
   Signs fixed and checked numerically to 1e-14.)
   General frame: Y_{lλ}(n)|_{k-frame} = Σ_m D^l_{mλ}(φ_k,θ_k,0) Y_lm(n),  D^l_{mλ}(φ_k,θ_k,0) = (−1)^λ √(4π/(2l+1)) [_{−λ}Y_lm(k̂)]*,
   hence  a^{(λ)}_lm = 4π(−i)^l ∫ d³k/(2π)³ ζ^{(λ)}_k (−1)^λ [_{−λ}Y_lm(k̂)]* Δ^{(λ)}_l(k),  Δ^{(λ)}_l = ∫dη g^{(λ)} Δ̂^{(λ)}_l(kχ).
   (D-matrix convention D^l_{mm'}(α,β,γ) = e^{−imα} d^l_{mm'}(β) e^{−im'γ}, Y_{lm'}(R⁻¹n) = Σ_m D^l_{mm'}(R) Y_lm(n); spin-weighted
   harmonics from the Goldberg explicit formula; both verified numerically.)
   For spectra only |Δ^{(λ)}_l|² matters; scalar/vector/tensor pieces of the same k do not cross-correlate (different m in the k-frame).

C^(11):  C^{AA'}_l = (2/π) ∫ k²dk P(k) Δ^A_l(k) Δ^{A'}_l(k)*   (any helicity; Δ real for real scalar maps).
   Derivation: (4π)² Σ_m|Y_lm(k̂)|² /(2π)³ ×4π = (2l+1)(2/π). Verified independently by the real-space route
   (monochromatic: ξ(n1·n2) = σ²|f|² j_0(k0χ|n1−n2|) ⇒ C_l = 4πσ²|f|² j_l(k0χ)², 1e-14).
Sky monopole of a type-A map (a random variable):  ⟨A⟩_sky = a_00/√(4π) = ∫ d³k/(2π)³ ζ_k Δ^A_0(k);
   ⟨⟨A⟩⟨A'⟩⟩ = (1/2π²) ∫ k²dk P Δ^A_0 Δ^{A'}_0 = C^{AA'}_0/(4π);   ⟨⟨A⟩ a^{A'}_lm⟩ = 0 for l ≥ 1 (isotropy).


## 2. Type B (products of two maps at the same n) — `typeB.py`
a^B_lm = Σ a_{l1m1} a'_{l2m2} G^{l l1 l2}_{m m1 m2}; "weighted" products a^B_lm = Σ w(l,l1,l2) G a a' cover ∇A·∇A' with
w_grad = ½[l1(l1+1)+l2(l2+1)−l(l+1)] (Hu 47/50). Wick with C^{XY}_l, ⟨a^X_{lm}a^Y_{l'm'}⟩ = (−1)^m C^{XY}_l δ_{ll'}δ_{m,−m'};
Σ_{m m1 m2}|G|² = (2l+1)(2l1+1)(2l2+1)/(4π)(l l1 l2;000)². For B = AA', B' = CD and l ≥ 1:

  (B1)  C^{BB'}_l = Σ_{l1 l2} (2l1+1)(2l2+1)/(4π) (l l1 l2;0 0 0)² w(l,l1,l2) w'(l,l1,l2) [ C^{AC}_{l1} C^{A'D}_{l2} + C^{AD}_{l1} C^{A'C}_{l2} ]

(parity l+l1+l2 even; the third pairing ⟨AA'⟩⟨CD⟩ feeds only l = 0.) Special case B = B' = AA':
  C^{BB}_l = Σ (2l1+1)(2l2+1)/(4π)(000)² [C^{AA}_{l1}C^{A'A'}_{l2} + C^{AA'}_{l1}C^{A'A}_{l2}]  — the two Wick pairings; NOT a factor 2
  on a single term unless A = A'. For A = A': C^{(A²)}_l = 2Σ(...)C^{AA}_{l1}C^{AA}_{l2}.
Monopole×map terms (SPEC §7), writing ⟨X⟩ := x_00/√(4π) and treating ⟨X⟩·1 as an l = 0 map (C^{⟨X⟩1,Y}_l = 4π⟨⟨X⟩⟨Y⟩⟩δ_{l0}),
(B1) collapses to (l ≥ 1):
  (i)   (1/(2l+1))Σ_m⟨|⟨X⟩ y_lm|²⟩ = ⟨⟨X⟩²⟩ C^{YY}_l
  (ii)  (1/(2l+1))Σ_m⟨(AA')_lm (⟨X⟩ y_lm)*⟩ = ⟨⟨A⟩⟨X⟩⟩ C^{A'Y}_l + ⟨⟨A'⟩⟨X⟩⟩ C^{AY}_l
  (iii) (1/(2l+1))Σ_m⟨(⟨X⟩y_lm)(⟨Z⟩w_lm)*⟩ = ⟨⟨X⟩⟨Z⟩⟩ C^{YW}_l
with ⟨⟨X⟩⟨Y⟩⟩ = (1/2π²)∫k²dk P Δ^X_0 Δ^Y_0 for type-A maps. Hence for Th2 = m2 − 2⟨m1⟩m1 with m2 = type-B:
  ¼C^(22)_l = ¼[C^{m2m2}_l − 4 Re(ii)(A,A'→ m2's factors; X=Y=m1) + 4⟨⟨m1⟩²⟩C^(11)_l].
Checks: (B1) vs exact real-space Wick + Legendre transform (ξ_B = ξ_{AC}ξ_{A'D} + ξ_{AD}ξ_{A'C}; for the gradient weight
ξ = ξ1''ξ2''(1−c²)² − c(1−c²)(ξ1''ξ2'+ξ1'ξ2'') + (1+c²)ξ1'ξ2' per pairing) agree to 1e-12 for correlated maps, both weights.
Monte Carlo (400 Gaussian realisations, exact band-limited SHT on a Gauss–Legendre grid): plain product, gradient product,
⟨X⟩·map and the cross term (ii) agree with (B1)/(i)/(ii) within 3–7% (≈ statistical error), excluding any factor-2 error.

## 3. CONTROL CASE: lensed C_l to O(C^φ) (Hu 2000; Lewis–Challinor 2006) — `lensing.py`
Θ̃ = Θ + ∇φ·∇Θ + ½∇^iφ∇^jφ∇_i∇_jΘ. With Θ = εTh1 + ε²Th2/2 + ε³Th3/6: Th2 = 2∇φ·∇Θ, Th3 = 3∇φ∇φ:∇∇Θ, so
¼C^(22) = ⟨|δ1_lm|²⟩ and ⅓Re C^(13) = 2Re⟨Θ_lm δ2*_lm⟩ — exactly the (first)² + 2(zeroth×second) split of LC (4.11). Bookkeeping OK.
(22): (B1) with w_grad ⇒ ¼C^(22)_l = Σ_{l1l2} F²_{ll1l2}/(2l+1) [C^{φφ}_{l1}C^{ΘΘ}_{l2} + C^{φΘ}_{l1}C^{Θφ}_{l2}],
      F = ½[l1(l1+1)+l2(l2+1)−l(l+1)]√((2l+1)(2l1+1)(2l2+1)/4π)(l l1 l2;000) = Hu (52); first term = Hu (62). The second term is the
      extension for ⟨Θφ⟩ ≠ 0 (relevant for us: SW/ISW correlate with the deflection potential).
(13): ⟨Θ_lm δ2*_lm⟩: pairing ⟨∇^iφ∇^jφ⟩_{same point} = g^{ij}R, R = ½⟨|∇φ|²⟩ = ½Σ l(l+1)(2l+1)/(4π) C^{φφ}_l, ∇²Y_lm = −l(l+1)Y_lm
      ⇒ ⅓Re C^(13)_l = −l(l+1) R C^{ΘΘ}_l = Hu (63). The two other pairings contain ⟨∇_i∇_jφ ∇_kΘ⟩ (no parity-even isotropic
      rank-3 tangent tensor ⇒ 0) and ∇_i⟨∇φ·∇Θ⟩ = 0, so they vanish even when ⟨Θφ⟩ ≠ 0.
Numerical verification: (a) (B1,w_grad) = ΣF²/(2l+1)(...) to 2e-16 (algebraic identity); (b) Monte Carlo of the Taylor-expanded
remapping on the sphere (δ2 built from first gradients: δ2 = ½[∇φ·∇(∇φ·∇Θ) − ½∇|∇φ|²·∇Θ]), with correlated φ,Θ:
⟨|δ1|²⟩ vs (22) formula: 2.3% ; 2Re⟨Θδ2*⟩ vs −l(l+1)RC_l: 4.6% (nreal = 400, consistent with MC noise ≈ 1/√nreal·O(1)).
SPEC §8 re-derived: C^(22) factor 2 (two connected pairings, K_2 symmetric) — correct; C^(13) factor 3 (three pairings k↔q_i,
remaining q_j = −q_k) — correct; SPEC's K*_{1,lm}K_{3,lm} is the complex conjugate of ⟨a^(1)a^(3)*⟩ — harmless since only Re C^(13)
enters, but the un-conjugated form is C^(13) = (3/(2l+1))Σ_m ∫d³k/(2π)³ P K_{1,lm}(k)[∫d³q/(2π)³ P(q) K_{3,lm}(k,q,−q)]*.
The "tadpole" pairing q_j = −q_k inside Th3 is the entire (13) term and must be kept (it is the K=0 piece of the 3-field kernel).

## 4. Type D (nested / unequal-time products) — `lensing.py` docstring
Θ_D(n) = ∫dη1 ∫_{η1}^{η_o}dη2 W(η1,η2) (δ_ij − n_in_j) ∂_iX(η1,−nχ1) ∂_jY(η2,−nχ2).
Angular-gradient identity: (δ−nn)_{ij}∂_iX∂_jY|_{x=−nχ} = (1/(χ1χ2)) ∇_nX̃·∇_nỸ, X̃(n) := X(η1,−nχ1) (a "shell map", whose a_lm is
4π(−i)^l ∫ζ_k Y*_lm(k̂) x(k,η1) j_l(kχ1)). Therefore a^D_lm = ∫dη1dη2 W/(χ1χ2) Σ w_grad G x̃_{l1m1}(η1) ỹ_{l2m2}(η2) and, by (B1),
  C^{DD}_l = Σ_{l1l2}(2l1+1)(2l2+1)/(4π)(000)² w_grad² ∫dη1dη2dη1'dη2' WW'/(χ1χ2χ1'χ2') [C^{xx}_{l1}(η1,η1')C^{yy}_{l2}(η2,η2') + C^{xy}_{l1}(η1,η2')C^{yx}_{l2}(η2,η1')],
  C^{xy}_l(η,η') := (2/π)∫k²dk P x(k,η) y(k,η')* j_l(kχ)j_l(kχ')  (unequal-time shell spectra). For CMB lensing X = Θ at η_* only
  (δ-function in η1) and the η2-integral is the lensing potential: reduces to §3. Radial contraction n·∂X n·∂Y: n·∂X(−nχ) = −∂_χX̃ ⇒
  replace j_l(kχ) → −k j_l'(kχ) in the shell a_lm, plain Gaunt (w = 1). Mixed radial/transverse pieces: same with one factor each.
Check: the Fourier form of the transverse contraction, Q = −k1k2[μ − (n·k̂1)(n·k̂2)] with e^{−ik1·nχ1−ik2·nχ2}, projected numerically,
equals the Gaunt/shell-map route to 3e-12 (random k1,k2,χ1≠χ2). So type D = two-time type C, and the l-dependent weight is exact.
