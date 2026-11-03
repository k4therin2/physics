"""n1 (V6): channel-by-channel response of A3's assembled kappa to the residual K=0 gauge mode with tau_A3(eta) = tau(eta)/(2 k^2)
(tau = pole gauge function, taufun.py; 1/p^2 = 1/(2k^2(1-mu)) at q = k), injected exactly as in A3 gauge.py (zero-pair tadpole shift
+ third-order shift).  By s1-s3 this response IS the 1/(1-mu) residue of the monochromatic integrand, channel by channel.
Output: Legendre coefficients l = 0..8 of each channel's response and of the total."""
import sys, pickle, time
sys.path.insert(0, '../A3_third_order')
import numpy as np
import tf
from thirteen import *
from taufun import tau as taufn
K0 = 1.0
TAU = (lambda t: taufn(K0, t) / (2 * K0**2), lambda t: taufn(K0, t, 1) / (2 * K0**2), lambda t: taufn(K0, t, 2) / (2 * K0**2))
def tadshift(t, q1m):
    ta, dta, ddta = (f(t)[None, :] for f in TAU)
    dPsi = dta + ta / t; dPhi = -ta / t; dde = -4 * ta / t
    ddPsi = ddta + dta / t - ta / t**2; ddPhi = -dta / t + ta / t**2
    dddPhi = -ddta / t + 2 * dta / t**2 - 2 * ta / t**3
    ddde = -4 * dta / t + 4 * ta / t**2
    return {'Phi': 2 * dPhi, 'dPhi': 2 * ddPhi, 'ddPhi': 2 * dddPhi, 'Psi': 2 * dPsi, 'dPsi': 2 * ddPsi, 'delta': 2 * dde,
            'Q': 2 * (dPhi - dPsi), 'dT': 2 * dde / 4, 'ddT': 2 * ddde / 4}
def third_shift(k):
    t = G.t; ta, dta, ddta = (f(t) for f in TAU)
    T, dT_, ddT = (tf.tphi(k, t, d) for d in (0, 1, 2))
    V, dV, ddV = (tf.tv(k, t, d) for d in (0, 1, 2))
    de, dde = tf.tdelta(k, t), tf.tdelta(k, t, 1)
    dPhi3 = 3 * ta * dT_ + 6 * ta / t * T
    ddPhi3 = 3 * dta * dT_ + 3 * ta * ddT + 6 * (dta / t - ta / t**2) * T + 6 * ta / t * dT_
    dPsi3 = 3 * ta * dT_ + 6 * (dta + ta / t) * T
    ddPsi3 = 3 * dta * dT_ + 3 * ta * ddT + 6 * (ddta + dta / t - ta / t**2) * T + 6 * (dta + ta / t) * dT_
    dV3 = 3 * ta * dV + 3 * dta * V
    ddV3 = 3 * dta * dV + 3 * ta * ddV + 3 * ddta * V + 3 * dta * dV
    dde3 = 3 * ta * dde - 12 * ta / t * de
    s = {'Phi': dPhi3, 'dPhi': ddPhi3, 'Psi': dPsi3, 'dPsi': ddPsi3, 'V': dV3, 'dV': ddV3, 'dT3lin': dde3 / 4}
    return {k_: 2 * v for k_, v in s.items()}
if __name__ == '__main__':
    LMAX = 8
    mu = float(sys.argv[1]) if len(sys.argv) > 1 else 0.31
    nu, wnu, nv = nu_grid(32)
    phis = 2 * np.pi * (np.arange(16) + 0.5) / 16
    D0 = (Delta_l(0, K0), Delta_l(0, K0))
    res = {}
    for shifted in (False, True):
        t0 = time.time()
        ts = tadshift if shifted else None
        th = third_profiles(K0, K0, mu)
        if shifted:
            sh = third_shift(K0); th = {k_: th[k_] + sh[k_] for k_ in th}
        s = np.sqrt(1 - mu**2); kv = np.array([[0, 0, K0]]); qv = K0 * np.array([[s, 0, mu]])
        m2 = M2_sky(np.vstack([kv, kv]), np.vstack([qv, -qv]), tadshift=ts)
        qq = K0 * np.array([[0.3, 0.4, np.sqrt(0.75)]])
        mt = M2_sky(qq, -qq, tadshift=ts)[0]
        ch, _ = coef_abc(K0, K0, mu, phis, nv, third=th, M2=(m2[0], m2[1], mt), D0=D0, tadshift=ts)
        res[shifted] = {k_: project(0.5 * v.mean(0), nu, wnu, LMAX) for k_, v in ch.items()}
        print('shifted=%s M2tad=%s (%.0fs)' % (shifted, mt, time.time() - t0), flush=True)
    resp = {k_: res[True][k_] - res[False][k_] for k_ in res[False]}
    tot = sum(resp.values()); scale = max(np.max(np.abs(v)) for v in resp.values())
    print('channel responses (Delta_l, l=0..%d; real parts; imag max %.1e):' % (LMAX, max(np.max(np.abs(v.imag)) for v in resp.values())))
    for k_, v in sorted(resp.items(), key=lambda kv: -np.max(np.abs(kv[1]))):
        if np.max(np.abs(v)) > 1e-12: print('  %-62s %s' % (k_[:62], ' '.join('%+.8f' % x for x in v.real)))
    print('  %-62s %s' % ('TOTAL', ' '.join('%+.1e' % x for x in tot.real)))
    print('max |total| = %.2e  (largest channel %.3f)' % (np.max(np.abs(tot)), scale))
    pickle.dump({'mu': mu, 'resp': resp}, open('n1_resp_mu%.2f.pkl' % mu, 'wb'))
