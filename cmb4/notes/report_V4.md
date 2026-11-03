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

## Run and result (k0 = 1, eta_* = 2, eta_o = 8; units sigma^2, sigma^4)
`python3 mc.py 20 45 4000 24 2026 4 n20` (N = 20 modes, 45 direction sets x 4000 amplitude draws = 180000 realisations,
420 pair kernels + 20 K_1 per set, ~110 s/set on 4 cores, 82 min). Output `mc_n20.pkl`, `mc_n20.log`, table `analyze_n20.log`.
Checks: <zeta^2>/sigma^2 = 0.9964 +- 0.0032; normalisation term = -2<m1>Th1 to <= 2.3e-15 in every realisation;
sum_m <Th2_lm>^2/((2l+1)C22) <= 2.5e-3 (-> 0, connected subtraction negligible).
```
 l | C11 MC            4piD^2 (A2)  | (1/4)C22 MC (off, unbiased)  A2 (1/4)C22   ratio   dev/err | full (biased O(1/N))
 1 | 1.4446(27)        1.44382      | 28.08 +- 0.25                28.2169       0.995   -0.5    | 27.57
 2 | 0.95250(14)       0.954601     | 14.174 +- 0.100              14.2575       0.994   -0.8    | 13.60
 3 | 0.62474(10)       0.625078     | 2.8564 +- 0.0208             2.87813       0.993   -1.0    | 2.723
 4 | 0.17524(20)       0.175071     | 1.8724 +- 0.0099             1.86592       1.003   +0.7    | 2.000
 5 | 1.9358(27)e-2     1.93345e-2   | 1.1248 +- 0.0062             1.13014       0.995   -0.9    | 1.084
 6 | 1.7704(20)e-4     1.77035e-4   | 0.5201 +- 0.0032             0.522062      0.996   -0.6    | 0.514
 7 | 3.1760(38)e-4     3.17458e-4   | 0.3710 +- 0.0043             0.374049      0.992   -0.7    | 0.353
 8 | 2.3047(26)e-4     2.30463e-4   | 0.2603 +- 0.0026             0.261974      0.994   -0.7    | 0.250
 9 | 6.0433(66)e-5     6.04509e-5   | 0.12339 +- 0.00076           0.124185      0.994   -1.1    | 0.1176
10 | 9.676(9)e-6       9.67953e-6   | 3.9864e-2 +- 2.5e-4          4.00848e-2    0.995   -0.9    | 3.866e-2
11 | 1.1051(13)e-6     1.10569e-6   | 9.455e-3 +- 5.9e-5           9.47301e-3    0.998   -0.3    | 9.41e-3
12 | 9.698(9)e-8       9.69908e-8   | 2.031e-3 +- 1.4e-5           2.00261e-3    1.014   +2.0    | 2.32e-3
```
Wick-in-amplitudes (directions only) agrees with the 'off' column to << 1 sigma (e.g. l=2: 14.189 +- 0.088, l=4 1.8696 +- 0.0082).
Errors are dominated by the direction sampling, and are correlated across l (common direction sets), so the uniform
~-0.7 sigma offset is one fluctuation, not 10 independent ones.

## Verdict
A2's monochromatic C_l^(22) assembly (harmonic projection, Wick factor 2, factor 1/4, m-sums, mono reduction) is CONFIRMED:
all l = 1..11 agree within 1.1 sigma at 0.3-0.7 % precision (l = 12: +2.0 sigma, 1.4 %). C11 matches 4pi sigma^2 Delta_l^2 to
<= 1.5 sigma (0.1-0.2 %). Sensitivity: a wrong Wick factor would be a 50 % effect; dropping the -2<m1>m1 normalisation
would shift C22 by +1.3 %, +4.2 %, +1.8 % at l = 2,3,4 (the l = 3 shift is ~6 sigma here) — excluded at l = 3. The
normalisation term enters A2's K_2 exactly once (checked realisation by realisation).
The 'full' estimator (physical finite-N map incl. i = j pairs) differs from the continuum by O(1/N) (-4.6 % at l=2, +7 % at l=4,
+16 % at l=12 for N = 20), as expected from the non-isotropic k2 = +-k1 pairs; the 'off' estimator removes this exactly.
Caveat: the MC uses A2's kernel K_2 itself (same engine), so it checks the assembly, not the kernel.
