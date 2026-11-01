# Report C2 — observable T_obs = T_o exp[-S] through third order (Poisson gauge)
Agent C2. Code: `cmb4/code/C2_observable/`. Conventions: SPEC.md. sympy 1.14.
Status legend: [DONE] [VERIFIED] [PARTIAL] [NOT DONE].

## Conventions fixed by C2 (please read)
* Series storage: every quantity is a truncated series `Ser` with coefficients c_n = X_n/n! (c1 = X1, c2 = X2/2, c3 = X3/6).
* Fields: sympy Functions of (eta,x,y,z): `Psi1,Psi2,Psi3, Phi1.., dT1.. (= delta_T), v_i_n (v^i), om_i_n (omega_i, n>=2),
  h_ij_n (n>=2, h_33 := -h_11-h_22 imposed; d_i h_ij = 0 = d_i om_i NOT imposed, left to the user)`.
* Direction: n = sky direction the observer looks at; propagation direction e = -n; background ray
  **x̄(eta) = x̄_o + n (eta_o^bar - eta)**, dx̄/deta = -n. (The brief wrote x_o - n(eta_o-eta); that has dx/deta=+n and is
  inconsistent with SPEC §6 (phat_o ∝ uhat_o - n^a e_a) and §9a (Doppler -n.v); I follow the SPEC.)
* Total derivative along the background ray: D := d_eta - n^i d_i.

## Milestone 1 — exact formula verified independently  [DONE, VERIFIED]  (`check_identity.py`)
(A) Kinematic identity, exact in a generic (even non-null) vector p, no equations of motion used:
    p^a p^b nabla_a u_b = sigma_pp + theta (p.p + E^2)/3 + E a.p ,  E = -u.p.  ZERO at orders 1,2,3.
    (Hence dE/dlambda = -p p nabla u = -(...), i.e. d ln Ehat/d lhat = -(theta/3 + sigma_ee + a.e) on the null cone.)
(B) Radiation identities, using ONLY nabla_mu That^{mu nu} = 0 (solved order by order for d_eta dT_n, d_eta v_n^i, n=1,2,3;
    no Einstein equations): theta + 3 u.d ln That = 0 and a_mu + h_mu^nu d_nu ln That = 0.  ZERO at orders 1,2,3.
(C) The two routes compared along the ray. Route (i): ln(Ehat_o/Ehat_e) = -int Gam^0_{kk} deta + [ln(-u.k)]_e^o
    (integrating dp^0/dlambda = -Gam^0 p p for p^0, E = p^0(-u.k)). Route (ii): ln(That_o/That_e) - S with
    S = int sigma_kk/(-u.k) deta. Writing the endpoint terms of both as integrals of total derivatives along the SAME
    perturbed ray, (i) - (ii) = int deta of
      R := -Gam^0_kk + d/deta[ln(-u.k) - ln That] + sigma_kk/(-u.k) + theta g_kk/(3(-u.k)),
    with k^i = -n^i + eps k1^i(eta) + eps^2 k2^i/2 + eps^3 k3^i/6 arbitrary functions, dk^i/deta substituted from the
    non-affine geodesic equation dk^i/deta = -Gam^i_kk + Gam^0_kk k^i, and d/deta = d_eta + k^i d_i on the fields.
    R == 0 at orders 1, 2 AND 3 modulo the conservation rules (no null condition, no |n|=1 needed; the theta g_kk
    term is the only trace of non-nullness and vanishes on the cone). Term counts before reduction: 14 / 70 / 340.
    => the integrated Ehat_o/Ehat_e agree term by term at every order once both are expanded about the same
    background ray (any common expansion of a vanishing integrand vanishes). Honest caveat: the comparison is done at the
    integrand level on the perturbed ray; the expansion about the background ray (milestone 3) is applied afterwards to the
    sigma-form only, and is checked separately (dilation check, milestone 4).
(D) Linear: with I1 := sigma1_kk (k=(1,-n)) = n^i n^j d_i v_j - d.v/3, verified  -I1 = -D[dT1+Psi1-n.v1] + Phi1' + Psi1'
    mod first-order conservation and |n|=1, hence
      L1 = -S1 = [dT + Psi - n.v]_e - [dT + Psi - n.v]_o + int_{eta_*}^{eta_o} (Phi'+Psi') deta   (SPEC §9a), observer
    terms included (the -[dT+Psi]_o piece is a monopole, +n.v_o the observer dipole).

## Milestone 2 — the ray through second order  [DONE, VERIFIED]  (`ray.py`, cache `ray_cache.pkl`)
Geodesic eq. in the eta parametrisation (exact): dk^i/deta = G^i = -Gam^i_kk + Gam^0_kk k^i, k=(1,k^i). Background ghat is
Minkowski, so the background ray is exactly straight, k^i = -n^i.
* Canonical triad: Gram-Schmidt of (d_1,d_2,d_3) wrt ghat and orthogonal to uhat (order of GS: 1,2,3). Verified
  g(e_a,e_b)=delta_ab, g(e_a,u)=0 through order 3. Explicit (second order, c2 = X2/2 convention):
  e_1^0 = eps v_1 + eps^2[ v_1_2/2 + om_1/2 - (Phi1+2Psi1) v_1 ],  e_1^1 = 1 + eps Phi1 + eps^2[ 3Phi1^2/2 + Phi2/2 - h_11/4 + v_1^2/2 ],
  e_1^2 = eps^2 v_1 v_2, e_1^3 = eps^2 v_1 v_3;  e_2: first order (v_2, 0, Phi1, 0), second order (x2):
  (-2Phi1 v_2 - 4Psi1 v_2 + om_2 + v_2_2, -h_12, 3Phi1^2 + Phi2 - h_22/2 + v_2^2, 2 v_2 v_3); e_3 analogous (GS asymmetry only in h_ab, v_a v_b).
* Direction map k^i(eta_o) = (u^i - n^a e_a^i)/(u^0 - n^a e_a^0), verified null through order 3 (with |n|=1):
    k1_init^i = v^i - n^i (n.v) - n^i (Phi1 + Psi1)                      (aberration v_perp; radial |k| = 1+Phi+Psi)
    k2_init^i/2 (displayed for n = z-hat, perp = 1,2, radial = 3):
       perp:   v_perp_2/2 + v_perp (n.v) + h_{perp z}/2
       radial: -Psi2/2 - Phi2/2 - 3Phi1^2/2 - Phi1 Psi1 + Psi1^2/2 + h_zz/4 - om_z/2 + v_perp^2/2
  (general n: `ray.k_init[i].c[2]`, 23-26 terms each). All fields at the observer EVENT (eta_o, x_o).
* First-order geodesic RHS (|n|=1):  G1^i = -d_i(Phi+Psi) + 2 n^i n.d(Phi+Psi) - n^i (Phi'+Psi')
  (transverse: -d_perp(Phi+Psi) = lensing deflection; radial n.G1 = -D(Phi+Psi), consistent with the null condition n.k1 = -(Phi+Psi)).
* Second order RHS G2/2 = [Gam2 at background point] + [Gam1 x k1 cross terms] (82-85 terms, `ray.G2half`, contains k1_i(eta)),
  plus the positional shift Dx1^i d_i G1 added in the functional definition.
* Ray functionals (all in terms of fields on the BACKGROUND ray; `ray.definitions` holds them as sympy Integral equations):
    k1(eta)    = k1_init|_(etabar_o,xbar_o) + int_eta^{etabar_o} G1 deta'
    Dx1(eta)   = dx_o1 + n deta_o1 - int_eta^{etabar_o} k1 deta'                  (Dx := x(eta) - xbar(eta))
    k2(eta)/2  = k2_init/2 + deta_o1 d_eta k1_init + dx_o1.d k1_init - deta_o1 G1(etabar_o) + int_eta^{etabar_o} [G2/2 + Dx1.d G1] deta'
    Dx2(eta)/2 = dx_o2/2 + n deta_o2/2 - int_eta^{etabar_o} k2/2 deta' - deta_o1 k1(etabar_o)
  (deta_o, dx_o = perturbation of the observer event, kept symbolic deta_o1, deta_o2, dxo1_i, dxo2_i.)
  Decomposition: n.Dx1 = radial (time-delay-like) part, Dx1 - n(n.Dx1) = transverse (lensing) part.
