"""Run every self-test of the C3 harmonic machinery and print a summary (deterministic tests ~1e-12; MC tests ~ few %)."""
import time
import typeA, typeB, lensing, typeC, thirteen, mono

t0 = time.time()
print("== type A (Rayleigh / helicity rules, C11 normalisation)")
typeA.test_typeA(); typeA.test_C11_normalisation()
print("== type B (products, Wick pairings, monopole terms)")
typeB.test_typeB_exact(); typeB.test_typeB_montecarlo()
print("== lensing control case (Hu 2000) and type-D identity")
lensing.test_lensing_control(); lensing.test_typeD_identity()
print("== type C (convolution kernels): matrix elements, isotropy reduction, C22 vs type-B, C x B")
typeC.test_typeC_matrix_elements(); typeC.test_typeC_vs_typeB()
print("== (13) reductions")
thirteen.test_rule_AL(); thirteen.test_qavg_reduction()
print("== monochromatic normalisation")
mono.test_mono_normalisation()
print("all tests done in %.0f s" % (time.time() - t0))
