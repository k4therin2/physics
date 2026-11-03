# Report V4 — Monte-Carlo referee of the monochromatic C_l^(22) assembly
Agent V4 (2026-11-03). Code: `cmb4/code/V4_mc/` (`mc.py`, `analyze.py`). Uses only A2's kernel evaluator
(`kernel2.MapKernels(geom,k1,k2,nv).K2()/K1()`, `Delta_l`) as a black box; independent of A2's `cl22.py` assembly
(no Khat=z frame, no mu/K quadrature, no Wick formula, no complex-Y m-sum identity).

## Method
Field zeta(x) = sum_{j=1}^N (sigma/sqrt N)[a_j cos(k_j.x) + b_j sin(k_j.x)], |k_j| = k0 = 1, isotropic random directions,
a_j,b_j ~ N(0,1).  zeta = sum_j [c_j e^{ik_j.x} + c.c.], c_j = sigma(a_j - i b_j)/(2 sqrt N); in the SPEC convention
zeta_k = (2pi)^3 sum_q c_q delta^3(k-q), q in {+-k_j}, so
  Th1(n) = sum_q c_q K_1(q;n),   Th2(n) = sum_{q,q'} c_q c_q' K_2(q,q';n) = sum_{i,j} 2Re[c_i c_j K_2(k_i,k_j) + c_i c_j^* K_2(k_i,-k_j)]
(all ordered pairs, incl. i=j and the K=0 pairs (k,-k), which the engine treats with its zero-mode prescription).
A2's K_2 already contains -Delta_0(k1)K_1(k2) - Delta_0(k2)K_1(k1); checked per realisation that its pair sum equals
-2<m1>Th1 with <m1> = sum_q c_q Delta_0(k0) (rel. diff <= 2e-15), i.e. the normalisation term is included exactly once.
a_lm: real orthonormal harmonics, Gauss-Legendre (24) x uniform phi (48) quadrature, l <= 12. Ensemble: S direction sets, each
with R amplitude draws (pair-kernel a_lm computed once per set). C_l = (1/(2l+1)) sum_m a_lm^2 averaged; standard errors from
the scatter of the S per-set means (sets are the independent units). Estimators:
 * full: physical Th2 of the finite-mode field; E = C22 (1 + O(1/N)) because i=j pairs (k2 = +-k1) are not isotropic pairs.
 * off: Th2 without i=j terms; exactly E = (N-1)/N C22 (distinct modes = independent isotropic pairs; the q,-q pairing absent);
   reported x N/(N-1): unbiased.
 * Wick-in-amplitudes: exact amplitude average of 'off' for the given directions (MC over directions only), diagnostic.
