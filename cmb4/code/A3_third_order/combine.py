"""combine.py (A3): eps^4 coefficient 1/4 C22 + 1/3 Re C13 (monochromatic, k0=1, eta_*=2, eta_o=8; units sigma^4),
C22 from A2 (code/A2_second_order/cl22_eo8_b.pkl), C13 from mono_<tag>.pkl; convergence vs a second run."""
import pickle, sys, numpy as np
a = pickle.load(open('mono_%s.pkl' % sys.argv[1], 'rb')); b = pickle.load(open('mono_%s.pkl' % sys.argv[2], 'rb'))
A2 = pickle.load(open('../A2_second_order/cl22_eo8_b.pkl', 'rb'))
C22 = np.asarray(A2['C22']); C11a2 = np.asarray(A2['C11'])
print('C11 check (A3 4pi Delta^2 vs A2):', np.max(np.abs(4 * np.pi * a['Delta'][2:21]**2 - C11a2[2:21]) / C11a2[2:21]))
print(' l  C11/s^2          1/4 C22/s^4      1/3 Re C13/s^4   (run diff rel)  sum/s^4           sum/C11')
rows = []
for l in range(2, 21):
    c13 = a['C13'][l].real / 3; c13b = b['C13'][l].real / 3
    s = C22[l] / 4 + c13
    rows.append((l, C11a2[l], C22[l] / 4, c13, abs(c13 - c13b) / abs(c13), s))
    print('%2d  %.9e  %.9e  %+.9e  (%.1e)  %+.9e  %+.3e' % (l, C11a2[l], C22[l] / 4, c13, abs(c13 - c13b) / abs(c13), s, s / C11a2[l]))
pickle.dump(rows, open('combined_eo8.pkl', 'wb'))
