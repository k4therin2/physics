"""n3 (V6): explicit closed expressions for the six channel residues as 1-D line-of-sight integrals of the first-order shear integrand
of the free leg (k = k zhat, x_ray = n chi, chi = eta_o - eta, nu = n.khat):
    I1(eta, nu) = -k^2 (nu^2 - 1/3) T_V(eta) e^{i k nu chi},   tt := tau/(2k^2)  (tau from taufun),
compared with the channel responses of n1 (= the residues).  Delta_l = (1/(2 i^l)) int dnu P_l(nu) R(nu)."""
import sys, pickle
sys.path.insert(0, '../A3_third_order')
import numpy as np
from numpy.polynomial import legendre as Lg
from scipy.special import eval_legendre
import tf
from taufun import tau as taufn
k, es, eo = 1.0, 2.0, 8.0
xe, we = Lg.leggauss(400); t = es + (xe + 1) * (eo - es) / 2; wt = we * (eo - es) / 2
nu, wnu = Lg.leggauss(64)
tt = lambda e, d=0: taufn(k, e, d) / (2 * k**2)
def I1(e, d_eta_TV=0):
    e = np.atleast_1d(e)
    return -k**2 * (nu[:, None]**2 - 1 / 3) * tf.tv(k, e, d_eta_TV)[None, :] * np.exp(1j * k * nu[:, None] * (eo - e)[None, :])
Bo = (tt(np.array([eo])) * I1(eo))[:, 0]
Bs = (tt(np.array([es])) * I1(es))[:, 0]
Ja = (I1(t, 1) * tt(t)) @ wt                     # tau~ * dT_V/deta part
Jb = (I1(t) * tt(t, 1)) @ wt                     # tau~' I1
Jc = (1j * k * nu[:, None] * I1(t) * tt(t)) @ wt # tau~ (n.grad) I1
Jd = (I1(t) * (tt(t) / t)) @ wt                  # tau~/eta I1
def proj(R, L=8): return np.array([np.sum(wnu * eval_legendre(l, nu) * R) / (2 * 1j**l) for l in range(L + 1)]).real
n1 = pickle.load(open('n1_resp_mu0.31.pkl', 'rb'))['resp']
nm = {'v3': 'L3: [ray] local | fields:3', 'TD': 'L3: [ray] time delay 2nd (radial Dx2) | fields:1', 'clock': 'L3: [obs] observer clock 2nd | fields:1',
      'emi': 'L3: [emi] emission-surface shift 2nd | fields:1', 'loc21': 'L3: [ray] local | fields:2+1', 'k2r': 'L3: [ray] radial k2 | fields:1'}
pred = {'v3': -3 * (Ja + Jb), 'TD': 3 * Jc, 'clock': 3 * Bo, 'emi': -3 * Bs}
for c, R in pred.items():
    a, b = proj(R), n1[nm[c]].real
    print('%-6s explicit %s\n       n1       %s   max diff %.1e' % (c, np.round(a, 8), np.round(b, 8), np.max(np.abs(a - b))))
# local 2+1 / radial k2: least-squares in the basis {Jb, Jd, Ja, Jc}
basis = {'Jb': proj(Jb), 'Jd': proj(Jd), 'Ja': proj(Ja), 'Jc': proj(Jc)}
A = np.stack(list(basis.values()), 1); y = n1[nm['loc21']].real
coef, *_ = np.linalg.lstsq(A, y, rcond=None)
print('local 2+1 fit coefficients', dict(zip(basis, np.round(coef, 10))), ' resid %.1e' % np.max(np.abs(A @ coef - y)))
print('radial k2 + local 2+1 :', np.max(np.abs(n1[nm['loc21']] + n1[nm['k2r']])))
tot = sum(pred.values())
print('explicit sum v3+TD+clock+emi (should vanish identically: total derivative):', np.max(np.abs(proj(tot))), np.max(np.abs(tot)))
