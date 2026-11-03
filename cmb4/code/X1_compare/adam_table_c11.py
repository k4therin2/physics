"""X1: our closed-form C11_2 = 4 pi Delta_2^2 (sigma=1) at Adam's App. A configurations vs his tabulated C2^(11)."""
import mpmath as mp
from transfer_compare import Delta_ours
s3 = mp.sqrt(3)
fixed_ratio = [(0.5, 0.00243), (1, 0.0258), (1.5, 0.0596), (2, 0.0423), (3, 0.116), (4, 0.958), (6, 0.0458), (8, 0.00059), (12, 0.00036), (16, 0.0533)]
fixed_em = [(2, 0.0340), (4, 0.946), (8, 0.0014), (12, 0.0302), (16, 0.00091)]
fr13 = [(0.15, 2.209e-5), (0.3, 3.414e-4), (0.6, 0.00475)]
print('fixed ratio x_* = x_o/3 (k0 = 1, eta = sqrt3 x):')
for xo, a in fixed_ratio + fr13:
    D, _ = Delta_ours(2, 1, s3 * xo / 3, s3 * xo)
    print('  x_o=%-5g ours C11=%.6g  Adam %.6g  rel %.1e' % (xo, 4 * mp.pi * D**2, a, (4 * mp.pi * D**2 - a) / a))
print('fixed emission x_* = 1:')
for xo, a in fixed_em:
    D, _ = Delta_ours(2, 1, s3, s3 * xo)
    print('  x_o=%-5g ours C11=%.6g  Adam %.6g  rel %.1e' % (xo, 4 * mp.pi * D**2, a, (4 * mp.pi * D**2 - a) / a))
