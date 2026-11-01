# Report C3 — harmonic-space machinery for C_l^(11), C_l^(22), C_l^(13)
(agent C3; code in `cmb4/code/C3_harmonic/`; run `python3 run_all_tests.py` there)

Status log
- [M1] Conventions fixed; type-A rules (scalar monomials, helicity ±1, ±2) derived and verified to 1e-14; C^(11) normalisation verified.
- [M2] Type-B product formula (both Wick pairings, monopole×map terms) verified exactly (1e-12, real-space route) and by Monte Carlo;
  lensing control case: Hu (62)-(63) reproduced from SPEC (22)/(13) bookkeeping (exact algebra + MC); SPEC §8 prefactors confirmed;
  type-D (nested/unequal-time) reduction and its Fourier↔angular-gradient identity verified (3e-12).
- [M3] Type-C kernels (finite 3j sums + μ-integral), C^(22) and C×B formulas verified to 1e-15 (incl. type-C ↔ type-B consistency);
  (13) reductions: generalised rule (A-L) and q-angular reduction (E1) verified to 1e-15; monochromatic collapse written out;
  full test suite `run_all_tests.py` (log in `test_log.txt`). DONE.

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

## 5. Type C (second-order convolution kernels) — `typeC.py`
Θ_C(n) = ∫d³k1d³k2/(2π)^6 ζ_{k1}ζ_{k2} ∫dη K(k1,k2,K,η) Q(n;k̂1,k̂2,K̂) e^{−iK·nχ}, K = k1+k2,
Q = Σ q_{abc}(k1,k2,μ)(n·k̂1)^a(n·k̂2)^b(n·K̂)^c (degree d; transverse/TT projectors contracted with n are of this form, e.g.
n^iP_ij(K̂)k̂1^j = (n·k̂1) − (n·K̂)(K̂·k̂1)). Define K_{2,lm}(k1,k2) := ∫dΩ Y*_lm ∫dη K Q e^{−iK·nχ} so a_lm = ∫∫ζζK_{2,lm}.
General frame:  K_{2,lm}(k1,k2) = ∫dη K Σ_{L=|l−d|}^{l+d} Σ_M 4π(−i)^L j_L(Kχ) Y*_{LM}(K̂) ⟨lm|Q|LM⟩,
  ⟨lm|Q|LM⟩ = ∫dΩ Y*_lm Q Y_LM: a finite 3j sum, generated by applying (a+b+c) times the coupling rule
  (n·v̂)Y_{lm} = (4π/3)Σ_μ Y*_{1μ}(v̂) Σ_{l'=l±1,m'} G^{l' l 1}_{m' m μ} Y_{l'm'}.  The L-sum is finite (|L−l| ≤ d) for polynomial Q.
Frame K̂ = ẑ, k1,k2 in the x–z plane (cosα1 = (k1+k2μ)/K, sinα1 = k2√(1−μ²)/K; cosα2 = (k2+k1μ)/K, sinα2 = −k1√(1−μ²)/K):
  K̃_lm(k1,k2,μ) = ∫dη K Σ_L (−i)^L √(4π(2L+1)) j_L(Kχ) ⟨lm|Q|L0⟩,  |m| ≤ d.
Isotropy: Σ_m|K_{2,lm}|² depends only on (k1,k2,μ) ⇒ ∫dΩ1dΩ2 → 8π²∫dμ, and with SPEC §8 (K_2 symmetrised):
  (C1)  C^(22)_l = (2/(2l+1))Σ_m∫∫d³k1d³k2/(2π)^6 PP|K_{2,lm}|² = 1/(4π⁴(2l+1)) ∫k1²dk1 k2²dk2 P(k1)P(k2) ∫_{−1}^{1}dμ Σ_m |K̃_lm(k1,k2,μ)|².
  (C2)  (1/(2l+1))Σ_m⟨C_lm B*_lm⟩ = 1/(4π⁴(2l+1)) ∫k1²dk1k2²dk2 PP ∫dμ Σ_m K̃^C_lm K̃^{B*}_lm,
        K̃^B_lm = ½Σ G^{l l1 l2}_{m m1 m2}[Ã_{l1m1}(k1)Ã'_{l2m2}(k2) + (A↔A')], Ã_{l1m1}(k) = 4π(−i)^{l1}Y*_{l1m1}(k̂)Δ^A_{l1}(k) (frame angles).
  Non-symmetric kernels must be symmetrised first (|K|² → K_lm(k1,k2)[K_lm(k1,k2)+K_lm(k2,k1)]*/2 is NOT enough for the m-sum; use ½(K+K^T)).
Two-time kernels (type D as type C) use the double Rayleigh expansion (bipolar) instead: K_{2,lm} = Σ_{l1l2}(4π)²(−i)^{l1+l2}∫∫dη1dη2 K j_{l1}(k1χ1)j_{l2}(k2χ2)
  Σ_{m1m2}Y*_{l1m1}(k̂1)Y*_{l2m2}(k̂2)⟨lm|Q|l1m1;l2m2⟩ (infinite but Bessel-truncated sums); the frame/μ-integral reduction (C1) is unchanged.
Checks: K_{2,lm} finite sum vs direct 2D n-quadrature (random k1,k2, two time nodes, Q of degree 3): 1e-14; Σ_m|K_lm|² random frame vs
K̂=ẑ frame: 1e-15; (C1) with K = ½[f(k1)g(k2)+g(k1)f(k2)], Q = 1 (i.e. the product map A·A') reproduces (B1) with monochromatic spectra
to 1e-15 [equivalently ∫dμ j_l(x√(2(1+μ)))² = 2Σ(2l1+1)(2l2+1)(l l1 l2;000)² j_{l1}(x)²j_{l2}(x)²]; (C2) vs (B1): 4e-15.

## 6. The (13) term — `thirteen.py`
C^(13)_l = (1/(2l+1))Σ_m⟨Th1_lm Th3*_lm⟩, Th3 = m3 − 3⟨m1⟩m2 − 3⟨m2⟩m1 + 6⟨m1⟩²m1 ⇒
  (E0) C^(13)_l = T3_l − 3T_{⟨1⟩2,l} − 3T_{⟨2⟩1,l} + 6⟨⟨m1⟩²⟩C^(11)_l,  ⟨⟨m1⟩²⟩ = (1/2π²)∫k²dk P Δ_0².
RULE (A-L) (generalised type-A rule, verified 3e-14): h(k,η)P_L(n·k̂)e^{−ik·nχ} ⇒ Δ_l = ∫dη h Σ_{Λ=|l−L|}^{l+L}(2Λ+1)(l L Λ;000)²(−i)^{Λ−l}j_Λ(kχ)
  (L=0: j_l; L=1: i j_l'; …). Every q-contracted kernel is a function of n·k̂ only, hence a Legendre series of this type.
q-ANGULAR REDUCTION (E1) (verified 1e-15, with and without Δχ):
  ∫dΩ_q/(4π) g(q̂·k̂)(n·q̂)^b e^{−iq·nΔχ} = Σ_{L'}(−i)^{L'}(2L'+1)j_{L'}(qΔχ) Σ_{L''}β_{L'bL''} g_{L''} P_{L''}(n·k̂),
  g_{L''} = ½∫dμ P_{L''}(μ)g(μ), x^bP_{L'} = Σβ P_{L''}; Δχ = 0 ⇒ only L' = 0: a finite polynomial of degree b with Legendre-moment coefficients.
  An extra (n·k̂)^a is absorbed by re-expanding ν^aP_{L''}(ν). Then (A-L).
(iii) Third-order kernel K3(k1,k2,k3,η)Q3 e^{−i(k1+k2+k3)·nχ} (symmetric): K3c(k,η;n) = ∫d³q/(2π)³P(q)K3(k,q,−q,η)Q3(n;k̂,q̂,−q̂)e^{−ik·nχ}
  = [∫q²dq/(2π²)P(q)Σ_L h_L(k,q,η)P_L(n·k̂)]e^{−ik·nχ}, h_L from (E1) with Δχ = 0 (depends on |k| only ⇒ isotropy manifest), and
  T3^{(iii)}_l = 3(2/π)∫k²dk P(k) Δ^{m1}_l(k) Δ^{(3c)}_l(k)*,  Δ^{(3c)}_l(k) = ∫q²dq/(2π²)P(q)∫dη Σ_L h_L Σ_Λ(2Λ+1)(l L Λ;000)²(−i)^{Λ−l}j_Λ(kχ).
  (factor 3 = three pairings for symmetric K3; otherwise 3K3(k,q,−q) → K3(k,q,−q)+K3(q,k,−q)+K3(q,−q,k).)
(i) m3 ⊃ ABC (three type-A maps): T3^{(i)}_l = C^{m1A}_l ξ^{BC}(0) + C^{m1B}_l ξ^{AC}(0) + C^{m1C}_l ξ^{AB}(0),
  ξ^{XY}(0) = ⟨X(n)Y(n)⟩ = Σ_{l'}(2l'+1)/(4π)C^{XY}_{l'} = (1/2π²)∫k²dk P ∫dη dη' x y* j_0(k|χ−χ'|) (for j_l-type maps);
  with gradient contractions ξ^{∇B∇C}(0) = Σ l'(l'+1)(2l'+1)/(4π)C^{BC}_{l'} and ∇B∇C:∇∇A → −l(l+1)R_{BC}C^{m1A}_l (as in §3).
(ii) m3 ⊃ A·C (type-A × type-C): T3^{(ii)}_l = C^{m1A}_l⟨⟨C⟩⟩ + 2(2/π)∫k²dk P Δ^{m1}_l(k) Δ^{AC}_l(k)*,
  ⟨⟨C⟩⟩ = ∫d³k/(2π)³P(k)K^C(k,−k;n) (K=0 tadpole mean), Δ^{AC}_l = type-A transfer function of ∫d³q/(2π)³P(q)K^A(q;n)K^C(k,−q;n)
  = ∫q²dq/(2π²)P(q)∫dη_Adη_C g_A K(k,q,|k−q|) × [(E1) with Δχ = χ_A − χ_C applied to Q_A(n,q̂)Q_C(n;k̂,−q̂)] e^{−ik·nχ_C} → (A-L).
Normalisation pieces (M2(q1,q2) := K_{2,00}(q1,q2)/√(4π) the sky-average kernel of m2):
  (E2) T_{⟨2⟩1,l} = ⟨⟨m2⟩⟩C^(11)_l + (1/π³)∫k²dk k'²dk' PP Δ_l(k)Δ_l(k')*∫dμ P_l(μ)M2(kk̂,−k'k̂')*,  ⟨⟨m2⟩⟩ = (1/2π²)∫k²dk P M2(k,−k);
  (E3) T_{⟨1⟩2,l} = 2(1/2π²)∫q²dq P(q)Δ_0(q)(2/π)∫k²dk P(k)Δ_l(k)Δ̄^{m2}_l(k;q)*, Δ̄_l(k;q) = ½∫dμ K_{2,l0}(kẑ,q)/((−i)^l√(4π(2l+1)))
  (q̂-averaged kernel; for type-B m2 = AA' it is ½[Δ^A_l(k)Δ^{A'}_0(q)+Δ^{A'}_l(k)Δ^A_0(q)], recovering §2 (ii) — consistency check of (E3)).
  The same objects give ¼C^(22) = ¼[C^{m2m2} − 4Re T^{(22)}_{2⟨1⟩1} + 4⟨⟨m1⟩²⟩C^(11)], T^{(22)}_{2⟨1⟩1} = (E3)* with the (22)-m2.

## 7. Monochromatic spectrum — `mono.py`
P(k) = 2π²σ²δ(k−k0)/k0² ⇔ ⟨ζ²⟩ = σ² (checked). Rules: ∫k²dk P F → 2π²σ²F(k0); ∫d³k/(2π)³ P F → σ²⟨F(k0k̂)⟩_{k̂}.
  C^(11): C^{AA'}_l = 4πσ² Δ^A_l(k0)Δ^{A'}_l(k0)*;  ⟨⟨A⟩⟨A'⟩⟩ = σ²Δ^A_0Δ^{A'}_0.
  Type B: C^{BB'}_l = (4πσ²)²Σ_{l1l2}(2l1+1)(2l2+1)/(4π)(000)² w w'[Δ^A_{l1}Δ^C_{l1}Δ^{A'}_{l2}Δ^D_{l2} + Δ^A_{l1}Δ^D_{l1}Δ^{A'}_{l2}Δ^C_{l2}] (all at k0;
     effectively finite: truncate at l ≈ k0χ_max + 10(k0χ_max)^{1/3} + 10).
  Type C: C^(22)_l = (σ⁴/(2l+1))∫_{−1}^{1}dμ Σ_{|m|≤d}|K̃_lm(k0,k0,μ)|², K = k0√(2(1+μ)), frame cosα1 = cosα2 = √((1+μ)/2).
  C×B: (σ⁴/(2l+1))∫dμ Σ_m K̃^C_lm K̃^{B*}_lm.
  (13)(iii): T3 = 12πσ²Δ_l(k0)Δ^{(3c)}_l(k0)*, Δ^{(3c)}_l(k0) = σ²Σ_L h_L(k0,k0)Σ_Λ(2Λ+1)(l L Λ;000)²(−i)^{Λ−l}j_Λ(k0χ), h_L = Legendre moments
     ½∫dμ P_L(μ)[K3(k0,k0,μ)×Q3-coefficients], |k±q| = k0√(2(1±μ)).  (i): ξ^{XY}(0) = σ²Σ_{l'}(2l'+1)Δ^X_{l'}Δ^{Y*}_{l'}.
     (ii): ⟨⟨C⟩⟩ = σ²K^C(k0k̂,−k0k̂;n); Δ^{AC}_l(k0) = σ²×[(E1) with q = k0, g(μ) = K(k0,k0,k0√(2(1−μ)))g_A …].
  Normalisation: ⟨⟨m1⟩²⟩ = σ²Δ_0(k0)²; ⟨⟨m2⟩⟩ = σ²M2(k0k̂,−k0k̂); (E2) cross = 4πσ⁴|Δ_l|²∫dμ P_l(μ)M2(k0k̂,−k0k̂')*; (E3) = 8πσ⁴Δ_0Δ_lΔ̄_l(k0;k0)*.
DELICATE POINTS: (1) μ → −1 in type-C (22): K → 0, dμ = KdK/k0², j_L(Kχ) → δ_{L0}; a kernel ~K^{−p} converges only for p < 1, so any
1/K² (Poisson inversion) must be accompanied by vanishing numerators (TT/transverse projections of k1⊗k2 at k2 = −k1 do vanish) — check the
K → 0 limit of every kernel analytically; Gauss–Legendre nodes never hit μ = −1 and can hide a divergence (test nmu-convergence).
(2) μ → ±1 in the (13) loops: |k∓q| → 0 inside K3(k,q,−q), K(k,q,|k−q|). (3) K = 0 tadpoles ⟨⟨m2⟩⟩, ⟨⟨C⟩⟩ require kernels exactly at K = 0;
they enter only via Th3 = m3 − 3⟨m1⟩m2 − 3⟨m2⟩m1 + 6⟨m1⟩²m1 where the sky-mean normalisation cancels the direction-independent mean of e^Λ —
individual tadpoles may be IR-singular while the combination is finite; keep symbolic until combined. (SPEC §9(d) constrains k1 → 0 of K_2,
a different limit.) (4) ⟨ζ_{k1}ζ_{k2}⟩ inside Th2 (K = 0) feeds only l = 0 — consistent with monopole removal.

## 8. What is verified (all in `run_all_tests.py`; deterministic unless marked MC)
- Rule A-s (s = 0..3), rule A-λ (λ = ±1, ±2; k-frame and random frame via D-matrix and spin-weighted forms): 1e-14.
- C^(11) normalisation (2/π) via the real-space route: 1e-14.
- (B1) both pairings, plain and gradient weights, correlated maps, vs exact real-space Wick + Legendre transform: 1e-12; MC (400 realisations,
  exact band-limited SHT): plain, gradient, ⟨X⟩·map, (ii) cross term within 3–7% (statistical).
- Lensing control: ¼C^(22) = ΣF²/(2l+1)C^φC^Θ exactly (2e-16); MC of the Taylor-expanded remap: (22) 2.3%, (13) = −l(l+1)RC_l 4.6% (correlated φ,Θ);
  SPEC §8 factors 2 and 3 re-derived and confirmed (SPEC's K*_1K_3 is the conjugate of ⟨a^(1)a^(3)*⟩; irrelevant for Re).
- Type-D Fourier ↔ angular-gradient identity (two radii): 3e-12.
- Type-C: K_{2,lm} finite sums vs quadrature 1e-14; isotropy reduction 1e-15; (C1) ≡ (B1) for product kernels 1e-15; (C2) 4e-15.
- (A-L) 3e-14; (E1) 1e-15 (Δχ = 0 and ≠ 0, with extra (n·k̂)^a). Monochromatic normalisation 6e-7 (Gaussian δ approximation).
NOT verified numerically (formulas only, derived by Wick counting consistent with the verified ones): (E2), (E3), T3^{(ii)} as written; the
collapse of (E3) to §2(ii) for type-B m2 is an analytic consistency check that passes.
