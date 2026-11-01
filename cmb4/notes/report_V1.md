# Referee report V1 — exact shear formula for T_obs(n)

**Verdict: CONFIRMED WITH CAVEATS.** No error found in the chain (i)–(iv). Every step was
re-derived with explicit indices and checked with sympy (scripts in `code/V1_referee/`), and the
formula was tested against three exact solutions and against linear Poisson-gauge theory. The
caveats (sec. 7) are conditions of validity / sharpenings of the surrounding statements in
SPEC sec. 3, not errors.

Conventions: signature (-+++), ∇_b u_a = σ_ab + ω_ab + (θ/3) h_ab − u_b a_a, a_a = u^b∇_b u_a,
h = g + uu, θ = ∇·u, σ_ab = h_(a^c h_b)^d ∇_d u_c − (θ/3)h_ab (as in the claim). E = −u·p.
p = E(u+e): p·p = 0 ⇒ e·e = 1 given u·e = 0. dl = E dλ.

Scripts (all PASS):
- `check_1_pointwise_identity.py` — items 1, 3 (generic metric, no symmetry, exact)
- `check_2_fluid_identities.py` — item 2 (generic metric, general ρ, P; thermodynamic algebra)
- `check_3_bianchi_milne.py` — items 4a, 4b (FRW, Bianchi I/Kasner test fluid, Milne)
- `check_4_poisson_sachs_wolfe.py` — item 4c (Poisson gauge, Sachs–Wolfe/ISW/Doppler)
- `check_5_tetrad_rotation.py` — item 6 (realisation-dependent rotation, C_l pieces)

---

## 1. Step (i): d ln E/dl = −(θ/3 + σ_ee + a·e)

Affine null geodesic p^ν∇_ν p^μ = 0. Then
dE/dλ = p^ν∇_ν(−u_μ p^μ) = −p^μ p^ν ∇_ν u_μ − u_μ (p^ν∇_ν p^μ) = −p^μ p^ν ∇_ν u_μ.
Insert the decomposition with p = E(u+e):
- σ: p^μp^ν σ_μν = E² σ_μν e^μe^ν (σ ⊥ u) = E² σ_ee.
- ω: p^μp^ν ω_μν = 0 identically (antisymmetric; symmetric contraction). **Vorticity drops with no
  assumption** — verified explicitly (`omega_{mn} p^m p^n == 0: True`).
- θ: (θ/3) h_μν p^μp^ν = (θ/3)E² h_ee = (θ/3)E².
- a: −p^μp^ν u_ν a_μ = −E²(u·(u+e))(a·(u+e)) = −E²(−1)(a·e) = E² a·e.
Hence dE/dλ = −E²(θ/3 + σ_ee + a·e), and with d/dl = (1/E)d/dλ:
**d ln E/dl = −(θ/3 + σ_ee + a·e).** Sign of each term as claimed.

Check 1 builds a metric g = G0 + G1_α x^α with random rational coefficients (Lorentzian, no
symmetry), a generic unit timelike u(x), generic T(x), generic null p at x=0, and verifies the
geodesic-derived dE/dλ against −p p ∇u, the decomposition, tracelessness/orthogonality, and the
final identity, all exactly: all `True`.

## 2. Step (ii): fluid identities, general EOS

With T^{μν} = (ρ+P)u^μu^ν + P g^{μν}, at a generic point (check 2, exact):
- u_ν ∇_μT^{μν} = −[u·∂ρ + (ρ+P)θ]  ⇒ energy: u·∂ρ + (ρ+P)θ = 0.
- h^ν_λ ∇_μT^{μλ} = (ρ+P)a^ν + h^{νμ}∂_μP ⇒ Euler: a_μ = −h_μ^ν ∂_νP/(ρ+P).
No Einstein equations used; the claim holds for a test fluid in a fixed spacetime as well.

Barotropic ρ(T), P(T):
- a = −h ∂ ln T ⇔ T P′(T) = ρ+P.
- θ = −3 u·∂ ln T ⇔ T ρ′(T) = 3(ρ+P).
Solving the pair: ρ′ = 3P′ ⇒ ρ = 3P + c, and T P′ = 4P + c ⇒ P = C T⁴ − c/4. With c = 0 this is
exactly {P = ρ/3, ρ = αT⁴}. Neither assumption alone suffices: P = ρ/3 alone needs ρ ∝ T⁴ for the
Euler form; ρ = αT⁴ with P = wρ gives a = −[4w/(1+w)] h∂lnT and θ = −[4/(1+w)] u·∂lnT, which
reduce to the claim iff w = 1/3. **No hidden chemical-potential / particle-number assumption**:
T is *defined* by ρ = αT⁴, and the Gibbs–Duhem-like relation ρ+P = T dP/dT is a consequence,
not an input. (A constant c ≠ 0 is a Λ-like shift and also works; not needed.)

## 3. Step (iii)

e^μ h_μ^ν = e^μ(δ_μ^ν + u_μu^ν) = e^ν + (e·u)u^ν = e^ν (checked). So
θ/3 + a·e = −u·∂lnT − e^μ h_μ^ν ∂_ν lnT = −(u+e)·∂lnT = −(1/E) p·∂lnT = −(1/E) d lnT/dλ = −d lnT/dl.
Then d lnE/dl = d lnT/dl − σ_ee, i.e. d/dλ[lnE − lnT] = −E σ_ee, and since p is future directed, λ
increases from e to o: **ln(E_o/E_e) = ln(T_o/T_e) − ∫_e^o σ_ee dl.**
Check 1 verifies the *off-shell* identity
[d lnE/dl − d lnT/dl + σ_ee] = −[(θ + 3u·∂lnT)/3 + (a + h·∂lnT)·e]
exactly at a generic point, so the claimed relation holds iff the two fluid identities hold
(the residual is generically nonzero off shell, −0.93 in the random instance).

## 4. Step (iv)

T_e = T_* on the emission surface; Liouville with f_e = f_Planck(E_e/T_*) gives f_o =
f_Planck(E_o /(T_* E_o/E_e)), a Planck spectrum in each direction (E_o/E_e is independent of the
photon's energy by linearity of the geodesic equation in p). Hence
**T_obs(n) = T_* E_o/E_e = T(o) exp[−∫_e^o σ_ee dl].**

## 5. Exact-solution tests (check 3)

(a) FRW radiation (Bianchi I with equal exponents): σ = 0 symbolically and T_* E_o/E_e / T_o = 1
for all p_i, symbolically.

(b) Bianchi I, ds² = −dt² + Σ t^{2q_i} dx_i², radiation test fluid u = ∂_t (a = 0, Euler holds since
T = T(t); θ + 3 d lnT/dt = 0 verified). p_i = const, E = (Σ p_i²/a_i²)^{1/2}, e^i = p_i/(a_i²E),
σ_ij = a_i²(H_i − θ/3)δ_ij, dl = dt. Symbolically d/dt[lnE − lnT] + σ_ee = 0 (with e^i(t) varying
along the ray). Integrated numerically for Kasner (2/3,2/3,−1/3) and non-Kasner (1/2,1/3,1/5),
three directions each, t_e = 1 → t_o = 5: |T_obs^exact − T_o e^{−S}| ≲ 2×10⁻³¹ (quadrature
precision), with genuine direction dependence (e.g. 0.342 vs 1.386 vs 1.233 for the Kasner case).

(c) Milne flow u = x/τ in Minkowski: θ = 3/τ, σ = 0, a = 0 symbolically; with p constant and
x_o = x_e + λp, E_o/E_e · τ_o/τ_e = 1 exactly ⇒ T_obs = T_* τ_e/τ_o = T(o), isotropic, as predicted.

## 6. Linear theory, Poisson gauge (check 4)

Conformal frame (SPEC sec. 2): ĝ = diag(−(1+2εΨ), (1−2εΦ)δ_ij), û = (1−εΨ, εv^i), T̂ = T̄(1+ε δ_T),
T̄ = const. Computed from the covariant definitions to O(ε):
θ̂ = ε(∂·v − 3Φ′), â_i = ε(v_i′ + ∂_iΨ), σ̂_ij = ε[∂_(i v_j) − δ_ij ∂·v/3], σ̂_0μ = 0 (Φ′ drops out).
Fluid identities ⇒ δ_T′ = −∂·v/3 + Φ′ and v_i′ = −∂_iΨ − ∂_iδ_T.

Independent observable: a_e T_* = T̂(e) = T̄(1+εδ_T(e)) (emission-surface shift included), so
ln(a_o T_obs/T̄) = ε δ_T(e) + ln(Ê_o/Ê_e).  Direct geodesic computation (Γ̂^0_00 = εΨ′,
Γ̂^0_0i = ε∂_iΨ, Γ̂^0_ij = −εΦ′δ_ij; Ê = p̂⁰(1 + εΨ − ε n·v)):
d lnÊ/dη = ε[−D(Ψ + n·v) + Ψ′ + Φ′], D := ∂_η + n^i∂_i (along the ray),
⇒ ln(a_o T_obs/T̄) = ε[δ_T + Ψ + n·v]_e − ε[Ψ + n·v]_o + ε∫(Ψ′+Φ′)dη — the textbook
Sachs–Wolfe + Doppler + ISW formula (n = propagation direction).

Claim: ln(a_o T_obs/T̄) = ε δ_T(o) − S, S = ∫σ̂_ee dη at O(ε). Equality ⇔ the integrand identity
**σ̂_ee = ε[D(δ_T + Ψ + n·v) − (Ψ′ + Φ′)] on shell.** Hand proof: n^in^j∂_iv_j = D(n·v) − n·v′ =
D(n·v) + n·∇Ψ + n·∇δ_T; n·∇δ_T = Dδ_T + ∂·v/3 − Φ′; n·∇Ψ = DΨ − Ψ′; subtract ∂·v/3. Sympy: off
shell the residual is nonzero, on shell it is identically 0. (Referee note: a first version of my
test had the sign of σ_ee wrong and reported a residual of exactly 2σ_ee, which is what the correct
identity predicts for the wrong sign; corrected and documented in the script.) Also verified (i)
again in this gauge: geodesic vs kinematic d lnÊ/dη agree.

## 7. Items 5–6 and caveats

**C1 (smoothness).** The proof uses the strong (differential) form of ∇·T = 0 and C¹ fields u, T
along the ray. At a shock (which a P = ρ/3 fluid can form nonperturbatively) Rankine–Hugoniot
replaces (ii), σ_ee becomes a distribution multiplying discontinuous E, and the formula is not
defined across the shock. Irrelevant for the perturbative programme.

**C2 (which emission event).** e must be the actual emission event: the *last* crossing of
T = T_* along the past ray (fluid opaque for T > T_*), and the ray must intersect T = T_*. If the
past ray crosses T = T_* several times, S must start at the last crossing; the identity of sec. 3
holds between any two points so nothing else changes. Nonperturbatively the T = T_* hypersurface
need not be spacelike everywhere; still fine. Perturbatively a single crossing.

**C3 (model vs physics).** Isotropic Planck emission in the fluid frame is self-consistent at
T = T_* (perfect fluid ⇔ isotropic in the rest frame; the photons *are* the fluid). After T_*
the "fluid" (u, T) in the formula is the model's continuation (SPEC sec. 0 defines it so); real
free-streaming radiation is not a perfect fluid and would source Φ ≠ Ψ. The claim is a theorem
in the stated model; nothing is "lost" in the derivation. The test photons' distribution differs
from the fluid's after T_* — that difference is exactly e^{−S}.

**C4 (direction label).** e is the propagation direction; the observer sees direction −e(o).
n → −n maps a_lm → (−1)^l a_lm, so C_l is unaffected; just fix the convention.

**C5 (T(o) and normalisation).** T(o) is a scalar at o, independent of n; 1+Θ = e^{−S}/⟨e^{−S}⟩
with ⟨·⟩ the solid-angle average in the observer's (= fluid) rest frame: exact, T(o) cancels.

**C6 (tetrad orientation, sharpened).** The FW tetrad is a realisation-dependent rotation
R(ε) = exp(εA) of the canonical triad. Check 5 (l = 2, exact spin-2 generators): Σ_m|a_lm|² is
invariant through O(ε⁴) and a^(1) is invariant, but a^(2) → a^(2) + 2A a^(1), so **C^(22)_l and
Re C^(13)_l are NOT separately invariant** (example: Σ|a₂|² = 321/25 vs 329107/1225 − 8798√6/1225),
while ¼C^(22) + ⅓Re C^(13) is (107/60 both). Hence Θ₁, Θ₂, Θ₃ must be computed in one and the
same frame; the SPEC's statement is correct for the total C_l only. The dipole power Σ_m|a_1m|²
is likewise invariant, so "remove dipole after normalising" is frame independent.

**C7.** The identities also hold for ρ = 3P + c, P = CT⁴ − c/4 (Λ-like constant); the claim's
hypotheses are sufficient and, within {P = wρ, ρ ∝ T⁴}, necessary (w = 1/3).
