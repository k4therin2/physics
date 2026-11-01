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

## Milestone 3 — S through third order, named line-of-sight terms  [DONE]  (`los.py` -> `L_terms.pkl`, log `los_run.log`)
### Method
S = int_{eta_e}^{eta_o} F deta, F(eta) = I(fields at (eta, x(eta)), k(eta)), I = sigma_kk/(-u.k) expanded to eps^3 with
k^i = -n^i + eps k1^i(eta) + eps^2 k2^i(eta)/2 (abstract functions). Then
 (1) position: F = I + Dx.dI + (1/2) Dx Dx ddI with Dx = eps Dx1(eta) + eps^2 Dx2(eta)/2 (fields now on the background ray);
 (2) endpoints: int_{eta_e}^{eta_o} = int_{eta_*}^{etabar_o} + sum_j [deta_o^j/j! D^{j-1}F]_obar - [deta_e^j/j! D^{j-1}F]_*
     with D = total derivative along the ray in the background representation (d_eta - n.d on fields, dk1/deta -> G1,
     dDx1/deta -> k1, ...), deta_o = eps deta_o1 + eps^2 deta_o2/2 (symbolic), deta_e solved from 1+dT(eta_e, x(eta_e)) = eta_e/eta_*:
       deta_e1 = eta_* dT1|_*
       deta_e2 = eta_* [ dT2 + 2 deta_e1 d_eta dT1 + 2 (Dx1(eta_*) - n deta_e1).d dT1 ]|_*      (x(eta_e) = x(eta_*) - n deta_e + ...)
 (3) at 'obs' the functionals take their observer values (k1 -> k1_init, Dx1 -> dx_o1 + n deta_o1, ...; milestone 2).
### Representation of the deliverable (pickle `L_terms.pkl`)
dict with: `L[m]` (m=1,2,3) = list of terms for the coefficient of eps^m in Lambda = -S (i.e. L1, L2/2, L3/6); each term is a dict
{name, where in {'ray','obs','emi'}, tags, expr, nterms}. 'ray' = integrate expr over eta from eta_* to etabar_o with fields at
(eta, xbar(eta)); 'obs' = fields at (etabar_o, xbar_o); 'emi' = fields at (eta_*, xbar_*). `L_split[m]`: same with n = z-hat
(k1_3 -> -(Phi1+Psi1)), tags lensing/time delay/aberration. `L_direct`, `L_direct_split`: the equivalent redshift ("direct") form
Lambda = -int Gam^0_kk + [ln(-u.k) - ln(1+dT)]_e^o (equal to -S modulo the conservation equations; see check T2 below).
`emission` = {deta_e1, deta_e2}; `definitions` = Integral forms of k1, Dx1, k2/2, Dx2/2; `m` = (m1,m2,m3); `Theta` = Theta_n in terms of
L_n and sky monopoles <m_n>. |n| = 1 is not imposed in the general-n expressions (reduce with ray.unit_rules if desired).
Sizes: L1: 15 terms; L2/2: 229 terms (general n) / 34 terms (n=z, 7 groups); L3/6: 2484 terms (general n) / 359 terms (n=z, 32 groups).
### L1  (shear form):  L1 = -int_{eta_*}^{etabar_o} I1 deta,  I1 = n^i n^j d_i v_j - d.v/3   (= standard formula after IBP, milestone 1D)
### L2  (shear form, n = z-hat, perp = x,y; all signs as they enter L2/2):
  L2/2 = int deta { [second-order shear]  -(1/4) h'_zz - (1/2) I1[v2]
                    + [quadratic local]    v_perp.(v_perp'/3 - d_z v_perp - d_perp v_z) + v_z (d_perp.v_perp/3 - 2 v_z'/3 - 2 d_z v_z/3)
                    + [direction]          k1_perp.(d_z v_perp + d_perp v_z)                     (aberration at o + deflection, propagated)
                    + [lensing]            - Dx1_perp . d_perp I1
                    + [time delay]         - Dx1_z d_z I1 }
         - deta_o1 I1|_o                       [observer clock]
         + deta_e1 I1|_e                       [emission-surface shift]
  (no explicit Psi1, Phi1 in the local quadratic terms: the (1-2Phi) of sigma_ij, the 1/(-u.k) = 1 - Psi - n.v factor and the radial
  null-condition piece n.k1 = -(Phi+Psi) cancel exactly, as they must for sigma_ee with e a unit vector; v' may be eliminated with the
  Euler equation v' = -d(Psi + dT) at first order if desired - not done, no equations of motion are used in L.)
  Observer displacement dx_o1 enters L2 only through Dx1 (lensing/time delay) and k2; it is a monopole/dipole-type term.
### L2  (direct/redshift form, n = z-hat; literature-friendly; coefficient of eps^2):
  ray:  (1/2)(Phi2' - Psi2') + d_z Psi2 + (1/4)h'_zz... written as -(1/2)Gam^0_kk[2nd-order fields] = (1/2)(Phi2'-Psi2') + d_z Psi2 + (1/4)(h'_11+h'_22) + (1/2) d_z om_z
        + 2 Phi1 Phi1' + 2 Phi1 d_z Psi1 + 2 Psi1 Psi1' - 2 Psi1 d_z Psi1               [quadratic local, incl. radial k1]
        - 2 k1_perp . d_perp Psi1                                                       [transverse direction]
        + Dx1_perp . d_perp (Phi1' - Psi1' + 2 d_z Psi1)  + Dx1_z d_z (Phi1' - Psi1' + 2 d_z Psi1)   [lensing, time delay of the ISW/SW integrand]
  obs:  (1/2)(Psi2 - dT2 + om_z + v2_z) - Phi1 v_z - Psi1^2 - Psi1 v_z + dT1^2/2 - v_perp^2/2
        + deta_o1 [ -(Gam^0_kk)_1 + d_eta(Psi1 - dT1 + v_z) ] + dx_o1 . d(Psi1 - dT1 + v_z)
  emi:  -(same local terms) + deta_e1 (Gam^0_kk)_1 - deta_e1 D(Psi1 - dT1 + v_z) - Dx1(eta_*).d(Psi1 - dT1 + v_z) + k1_perp(eta_*).v_perp
  where (Gam^0_kk)_1 = Psi1' - 2 d_z Psi1 - Phi1', D = d_eta - d_z.  The boundary function is ln(-u.k) - ln(1+dT) =
  Psi + n.v - dT + (second order: Psi^2... see pickle), i.e. the standard SW+Doppler combination; the emission-surface shift
  deta_e1 = eta_* dT1 reproduces the "temperature at the background time" bookkeeping of the standard derivations.
### L3 (shear form, n = z-hat): 32 named groups (359 terms); list in `los_run.log` / `L_terms.pkl['L_split'][3]`:
  ray: third-order shear (fields:3, 5 terms); second x first local (57); cubic local (55); k1_perp x (1+1) (22), k1 x 2nd-order (6), k1^2 (8),
       transverse k2 (4), radial k2 (3); lensing Dx1 x (1+1) (50), Dx1 x 2nd (10), Dx1_perp^2 (9), Dx2_perp (6); time delay Dx1_z x (1+1) (26),
       Dx1_z x 2nd (5), Dx1_z^2 (3), Dx2_z (3); mixed Dx1_perp Dx1_z (6), k1 Dx1_perp (8), k1 Dx1_z (4).
  obs: deta_o2 (3), deta_o1 x (1+1) (6), deta_o1 x 2nd (5), deta_o1 dx_o1 (9), deta_o1^2 (3+6).
  emi: deta_e2 (3), deta_e1 x (1+1) (10), deta_e1 x 2nd (5), deta_e1 x lensing (6), deta_e1 x time delay (3), deta_e1 x k1 (4), deta_e1^2 (6).
### m_n and monopole subtraction (SPEC §7) — verified algebraically from exp/ratio series:
  m1 = L1, m2 = L2 + L1^2, m3 = L3 + 3 L1 L2 + L1^3;
  Theta1 = m1 - <m1>;  Theta2 = m2 - <m2> - 2<m1>(m1 - <m1>);
  Theta3 = m3 - <m3> - 3<m1>(m2 - <m2>) - 3<m2>(m1 - <m1>) + 6<m1>^2 (m1 - <m1>)   => for l>=1 exactly SPEC §7.

## Check 4(b) — dilation  [DONE, VERIFIED]  (`check_dilation.py`, `dilation_run.log`)
Exact FRW in rescaled coordinates eta = lam etat, x = mu xt (a = etat):  ghat_00 = -lam^4, ghat_ij = lam^2 mu^2 delta_ij,
1+dT = 1/lam, v = 0, uhat^0 = lam^-2, zeta_nl = ln mu EXACTLY.  With lam = 1 + eps a1 + eps^2 a2/2, mu = exp(eps b1 + eps^2 b2/2):
 Psi1 = 2a1, Phi1 = -a1 - b1, dT1 = -a1  => Psi1 = Phi1 = -2c/3, dT1 = c/3 fixes a1 = -c/3, b1 = c  (dT1 = c/3 as in SPEC §5).
 zeta_nl^(2) = b2 => b2 = 0; a2 is FREE (residual constant time rescaling). Second-order constants:
 Psi2 = 2a2 + 2c^2/3,  Phi2 = -a2 - 7c^2/9,  dT2 = -a2 + 2c^2/9  (third order also generated: Psi3 = -6a2 c - 4c^3/9, Phi3 = -5a2 c - 2c^3/3, dT3 = -2a2 c + 2c^3/9).
Results (all symbolic in c, a2):
 * sigma_mu nu == 0 identically through order 3 for these constant fields; uhat^0 = lam^-2 reproduced through order 3.
 * L1, L2, L3 from the assembled term lists (shear AND direct forms, with deta_e1, deta_e2 substituted by their solutions
   deta_e1 = c eta_*/3, deta_e2 = (-a2 + 2c^2/9) eta_*, which equal the exact eta_*(1/lam - 1) series) are ZERO at orders 1,2,3 => Lambda = 0.
 * Observer clock: Delta tau = lam^2 (eta_o^2 - eta_init^2)/2 with eta_init = eta_i/lam, fixed to its background value gives
   eta_o = etabar_o/lam through second order (any a2), hence T_o/Tbar(etabar_o) = (1+dT)(etabar_o/eta_o) = 1 + O(eps^3):
   the observed monopole is unperturbed. Emission: eta_e = eta_*/lam, a_e/a_o = eta_*/etabar_o exactly.
 => T_obs(n) = T_o e^{-S} = Tbar(etabar_o): Theta = 0 through second order (in fact to all orders), INCLUDING clock and emission surface.
 Honest remark: in the exact e^{-S} formulation the dilation check is automatically trivial for the anisotropy (S depends only on
 derivatives of the fields) — the clock and emission-surface machinery is only exercised by the monopole T_o (verified above) and by the
 direct form (route i), whose assembled L1..L3 also vanish. The nontrivial content is thus (a) S has no undifferentiated-field terms
 (true by construction: sigma is built from nabla u), (b) T_o = Tbar(etabar_o) at fixed proper time.
