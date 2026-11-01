# Report C3 — harmonic-space machinery for C_l^(11), C_l^(22), C_l^(13)
(agent C3; code in `cmb4/code/C3_harmonic/`; run `python3 run_all_tests.py` there)

Status log
- [M1] Conventions fixed; type-A rules (scalar monomials, helicity ±1, ±2) derived and verified to 1e-14; C^(11) normalisation verified.

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

(sections 2–6 follow in later milestones)
