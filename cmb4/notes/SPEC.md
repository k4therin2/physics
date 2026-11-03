# MASTER SPEC — CMB temperature power spectrum to O(eps^4), radiation-only GR
(Lead: this file is the single source of truth for conventions. Every agent reads it first.)

## 0. The problem (verbatim essentials)
Flat FRW, GR, Lambda=0, radiation-only perfect fluid P=rho/3, rho = alpha T_phys^4.
First order: regular adiabatic scalar perturbations. Nonlinear primordial variable zeta Gaussian,
no independent higher-order primordial modes. Photons emitted isotropically (Planck, T_*) on the
physical hypersurface T_phys = T_*, propagate collisionlessly. Observer = tagged fluid worldline,
measures at fixed elapsed proper time from a physical initial event; tetrad FW-transported.
T_obs(n) = T_* E_o/E_e, E = -u.p. Theta(n) = T_obs(n)/[(1/4pi) int dOmega T_obs] - 1. Dipole removed
only after normalising. Theta = eps Th1 + eps^2/2 Th2 + eps^3/6 Th3.
a_lm^(r) = int dOmega Y*_lm Th^(r); C_l^(rs) = 1/(2l+1) sum_m <a^(r) a^(s)*>_c.
C_l = eps^2 C^(11) + eps^4 [ 1/4 C^(22) + 1/3 Re C^(13) ] + O(eps^6).
Everything analytic; declare unevaluated integrals; monochromatic spectrum at k0 only at the end.

## 1. Units, background
Signature (-+++). Conformal time eta. Choose a(eta) = eta (a_i/eta_i = 1). H = a'/a = 1/eta.
Friedmann: 3H^2 = 8 pi G a^2 rhobar  =>  8 pi G a^4 rhobar = 3, i.e. 8 pi G rhohat_bar = 3 where
rhohat = a^4 rho is the CONFORMAL energy density (constant in background). Tbar(eta) = Tbar_i eta_i/eta.
Conformal temperature That = a T_phys ; That_bar = const.

## 2. Conformal frame (exact, nonperturbative)
g_{mu nu} = a(eta)^2 ghat_{mu nu}.  Null geodesics of g = null geodesics of ghat (phat^mu = a^2 p^mu).
Fluid: uhat^mu = a u^mu (unit in ghat). Radiation stress tensor is traceless, hence conformally
covariant: That^{mu nu} = a^6 T^{mu nu} = (rhohat+Phat) uhat uhat + Phat ghat^{mu nu},
rhohat = a^4 rho = alpha That^4, Phat = rhohat/3, and nablahat_mu That^{mu nu} = 0 EXACTLY.
Energies: E = -u.p = Ehat/a with Ehat = -ghat_{mu nu} uhat^mu phat^nu.
Consequences (exact): acceleration ahat_mu = -hhat_mu^nu d_nu ln That ; expansion thetahat = -3 uhat.d ln That.

## 3. THE EXACT OBSERVABLE (lead's key result; every agent must respect it)
For any null geodesic with affine tangent phat and the fluid field uhat:
   d ln Ehat / d lhat = -( thetahat/3 + sigmahat_{mu nu} e^mu e^nu + ahat.e ),
 where phat = Ehat (uhat + e), e = unit propagation direction in fluid frame, d lhat = Ehat d lambdahat.
Using the radiation-fluid identities,  thetahat/3 + ahat.e = -(1/Ehat) d ln That/d lambdahat, so
   ln(Ehat_o/Ehat_e) = ln(That_o/That_e) - S ,   S(n) := int_e^o sigmahat_{mu nu} e^mu e^nu d lhat .
With That_e = a_e T_*, That_o = a_o T_o:
   T_obs(n) = T_o * exp[ - S(n) ]                                              (EXACT)
T_o (local fluid temperature at observer) is direction independent, so
   1 + Theta(n) = exp[-S(n)] / < exp[-S] >_sky .
S is first order in eps (background shear vanishes). All Sachs-Wolfe/ISW/Doppler/emission-surface/
lensing/time-delay physics is contained in S. Conformal rescaling drops out: S = int sigma_ee dl in the
physical metric too. NOTE: the observer's tetrad orientation (FW transport) is a rigid rotation of
the sky => drops out of C_l exactly (sum_m |a_lm|^2 invariant realisation by realisation). Use the
canonical triad (Gram-Schmidt of coordinate basis wrt ghat, orthogonal to uhat).

Parametrising the ray by eta, with k^mu = dx^mu/d eta = (1, dx^i/deta):
   S(n) = int_{eta_e}^{eta_o} d eta  [ sigmahat_{mu nu} k^mu k^nu ] / [ -uhat_mu k^mu ]   along x^i(eta).
(sigmahat is orthogonal to uhat so sigma_{mu nu} p p = Ehat^2 sigma_ee ; Ehat/phat^0 = -uhat.k.)

## 4. Gauge and field parametrisation (Poisson gauge in the conformal metric)
ghat_00 = -(1 + 2 Psi),  ghat_0i = omega_i (d_i omega_i = 0),  ghat_ij = (1 - 2 Phi) delta_ij + h_ij (h_ii = 0 = d_i h_ij).
Expansions: X = eps X1 + eps^2 X2/2 + eps^3 X3/6 for X in {Psi, Phi, omega_i, h_ij, delta, v^i};
omega and h start at order eps^2 (no primordial vectors/tensors). NO exponential parametrisation.
Fluid: rho = rhobar(eta)(1+delta);  1 + delta_T := (1+delta)^{1/4}.  uhat^mu = uhat^0 (1, v^i),
uhat^0 fixed by ghat uhat uhat = -1. Do NOT impose irrotationality by hand; verify it (curl of uhat_i = 0).
Christoffels, Einstein tensor: G_{mu nu}[g] with g = eta^2 ghat (a = eta). Equations:
G_{mu nu} = 8 pi G T_{mu nu}, with 8 pi G a^4 rhobar = 3. Also nabla_mu T^{mu nu} = 0 (redundant; use as check).

## 5. Primordial variable and initial data
zeta_nl(x) := lim_{eta->0} [ ln(1 + delta_T(eta,x)) + (1/2) ln(1 - 2 Phi(eta,x)) ]   (Poisson-gauge fields,
 nonlinear delta_T and Phi). Justification: the uniform-density slice through Poisson time eta_P = eta + T has
 1+T/eta = 1+delta_T, its induced metric at eta->0 is a^2 (1+T/eta)^2 (1-2Phi) delta_ij + O(eta^2), tensors
 vanish there for regular solutions, so zeta_nl is Maldacena's zeta (uniform-density curvature perturbation).
REQUIREMENT: zeta_nl = eps zeta_1 exactly, zeta_1 Gaussian:  zeta_nl^(2) = zeta_nl^(3) = 0 at eta->0.
"Regular": solution analytic at eta=0 (no decaying modes) at every order; no homogeneous tensor/vector modes.
Initial data at finite eta_i := values of this regular solution at eta_i (satisfies constraints by
construction; keep eta_i explicit in lower limits of Green's integrals where asked, note eta_i->0 reproduces).
Expected linear solution (to be DERIVED, not assumed): with x = k eta / sqrt3,
   Phi_1 = Psi_1 = zeta_1 * T_Phi(x), T_Phi(x) = -2 (sin x - x cos x)/x^3  (-> -2/3 as x->0),
   delta_T1 -> -Phi_1/2 superhorizon, v_1 -> -(eta/2) grad Phi-type.
Fourier: f(x) = int d^3k/(2pi)^3 e^{ik.x} f_k ; <zeta_k zeta_k'> = (2pi)^3 delta^3(k+k') P(k).

## 6. Observer and emission
Observer worldline: tagged fluid element, x_o(eta) = x_o^label + int_0^eta v^i(eta', x_o(eta')) d eta'  (label = position as eta->0; A1 2026-11-03).
Clock: proper time from the crossing of the initial-data surface T_phys = T_i (Poisson-time eta_init(x_o),
1+delta_T(eta_init,x_o) = eta_init/eta_i): Delta tau = int_{eta_init}^{eta_o} a(eta)/uhat^0 d eta  fixed
=> eta_o = eta_o^bar + delta eta_o[zeta] (direction independent).
Emission event: ray meets T_phys = T_*  <=>  1 + delta_T(eta_e, x(eta_e)) = eta_e/eta_*  (since Tbar ∝ 1/eta).
Direction map: phat_o = Ehat_o ( uhat_o - n^a e_a ), canonical triad e_a; k^i(eta_o) = phat^i/phat^0.

## 7. Normalisation bookkeeping (exact algebra)
Lambda := -S = eps L1 + eps^2 L2/2 + eps^3 L3/6.  M := e^Lambda = 1 + eps m1 + eps^2 m2/2 + eps^3 m3/6 with
 m1 = L1, m2 = L2 + L1^2, m3 = L3 + 3 L1 L2 + L1^3.   Mbar := (1/4pi) int dOmega M.
For l >= 1:  a_lm[Theta] = a_lm[M]/Mbar  =>
 Th1_lm = m1_lm ;  Th2_lm = m2_lm - 2 <m1> m1_lm ;  Th3_lm = m3_lm - 3<m1> m2_lm - 3<m2> m1_lm + 6 <m1>^2 m1_lm,
 where <m> := sky monopole (1/4pi) int dOmega m  (a random variable correlated with the map).

## 8. Loop structure
Th^(r)(n) = int Pi d^3k_i/(2pi)^3 zeta_{k_1}..zeta_{k_r} K_r(k_1..k_r; n), K_r symmetric.
C^(22)_l = (2/(2l+1)) sum_m int d^3k1 d^3k2/(2pi)^6 P(k1)P(k2) |K_{2,lm}(k1,k2)|^2  (l>=1),
C^(13)_l = (3/(2l+1)) sum_m int d^3k/(2pi)^3 P(k) K*_{1,lm}(k) [ int d^3q/(2pi)^3 P(q) K_{3,lm}(k,q,-q) ].
(Agents: re-derive signs/conjugations; the Wick factors 2 and 3 come from the pairings.)

## 9. Mandatory checks
(a) S_1 reproduces the standard linear formula Th1 = [delta_T + Psi - n.v]_e - [same]_o + int (Phi'+Psi') d eta.
(b) Exact formula vs direct second-order integration of the geodesic energy (independent derivation).
(c) Constraint consistency at each order (Bianchi), vorticity-free, uhat normalisation.
(d) Constant zeta_1 = c is pure gauge: the FULL map (with observer clock, normalisation) must be unchanged
    => K_2(k1->0,k2) and K_3 squeezed limits are fixed; K_2(k1=0,.)=0 structure. (dilation check)
(e) Literature: Kohri-Terada RD induced-tensor kernel; Bartolo et al second-order Poisson-gauge RD ICs.
