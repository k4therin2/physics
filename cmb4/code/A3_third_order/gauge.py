"""gauge.py (A3): zero-mode (tadpole) gauge test.  Residual K=0 Poisson-gauge freedom eta -> eta + eps^2 tau(eta)/2
(tau = second-order homogeneous, per P(q): tau~).  Exact transformation (a = eta, rho_bar ~ eta^-4):
  d tadPsi2 = tau' + tau/eta,  d tadPhi2 = -tau/eta,  d taddelta2 = -4 tau/eta     (zeta_nl unchanged for regular tau)
  third order (C1 normalisation, per P(q), slot-k transfer functions T_X):
  dPhi3c = 3 tau T_Phi' + 6 (tau/eta) T_Phi,  dPsi3c = 3 tau T_Phi' + 6 (tau' + tau/eta) T_Phi,
  dV3c = 3 tau T_V' + 3 tau' T_V,  ddelta3c = 3 tau T_delta' - 12 (tau/eta) T_delta.
Test: kappa(nu) at fixed mu with (tadpoles + d) in every K=0 pair object (L3 local/functionals/clock/emission/<m2>) and
third-order fields + 2 d X3c must equal the unshifted kappa.   Second part: C1's K3c sources with shifted tadpoles,
solved by the same Green's function, must reproduce dPhi3c up to the homogeneous mode fixed by dPhi3c(eta->0)."""
import numpy as np, sys, time
import tf
from thirteen import *
from third import X3c

def taus(kind, a=1.0):
    if kind == 'lin': return (lambda t: a * t, lambda t: a + 0 * t, lambda t: 0 * t)
    if kind == 'cub': return (lambda t: a * t**3, lambda t: 3 * a * t**2, lambda t: 6 * a * t)
TAU = taus(sys.argv[1] if len(sys.argv) > 1 else 'lin', 0.1)
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
    k0 = q0 = 1.0; mu = 0.31
    nu, wnu, nv = nu_grid(8)
    phis = 2 * np.pi * (np.arange(8) + 0.5) / 8
    D0 = (Delta_l(0, k0), Delta_l(0, q0))
    res = {}
    for shifted in (False, True):
        ts = tadshift if shifted else None
        th = third_profiles(k0, q0, mu)
        if shifted:
            sh = third_shift(k0); th = {k_: th[k_] + sh[k_] for k_ in th}
        s = np.sqrt(1 - mu**2); kv = np.array([[0, 0, k0]]); qv = q0 * np.array([[s, 0, mu]])
        m2 = M2_sky(np.vstack([kv, kv]), np.vstack([qv, -qv]), tadshift=ts)
        qq = q0 * np.array([[0.3, 0.4, np.sqrt(0.75)]])
        mt = M2_sky(qq, -qq, tadshift=ts)[0]
        ch, _ = coef_abc(k0, q0, mu, phis, nv, third=th, M2=(m2[0], m2[1], mt), D0=D0, tadshift=ts)
        res[shifted] = (0.5 * sum(ch.values()).mean(0), {k_: 0.5 * v.mean(0) for k_, v in ch.items()}, mt, m2)
        print('shifted=%s  M2tad=%s  kappa_mu(nu) =' % (shifted, mt), np.round(res[shifted][0], 6), flush=True)
    d = res[True][0] - res[False][0]
    print('max |kappa(shifted) - kappa| =', np.max(np.abs(d)), '  scale', np.max(np.abs(res[False][0])))
    print('largest channel changes (cancel in the sum):')
    chg = {k_: np.max(np.abs(res[True][1][k_] - res[False][1][k_])) for k_ in res[False][1]}
    for k_, v in sorted(chg.items(), key=lambda kv: -kv[1])[:8]: print('   %-70s %.4g' % (k_, v))
    # part 2: C1 K3c sources with shifted tadpoles
    t = G.t
    ta, dta, ddta = (f(t) for f in TAU)
    base = {nm: (lambda d, nm=nm: tf.TAD[nm](q0, t, d)) for nm in ('tadPhi2', 'tadQ2', 'taddelta2')}
    dPh = -ta / t; dPs = dta + ta / t
    sh = {'tadPhi2': [dPh, -dta / t + ta / t**2, -ddta / t + 2 * dta / t**2 - 2 * ta / t**3],
          'tadQ2': [dPh - dPs, (-dta / t + ta / t**2) - (ddta + dta / t - ta / t**2)], 'taddelta2': [-4 * ta / t]}
    shifted = {nm: (lambda d, nm=nm: base[nm](d) + sh[nm][d]) for nm in base}
    mus = np.array([0.31])
    dPhi3_expected = (3 * ta * tf.tphi(k0, t, 1) + 6 * ta / t * tf.tphi(k0, t))
    C3shift = 6 * TAU[0](1e-8) / 1e-8 * tf.tphi(k0, 1e-8)    # dPhi3c(eta -> 0)
    X0 = X3c(G, k0, q0, mus); X1 = X3c(G, k0, q0, mus, tads=shifted, C3shift=np.array([C3shift]))
    for nm, exp in (('Phi', dPhi3_expected), ('Psi', 3 * ta * tf.tphi(k0, t, 1) + 6 * (dta + ta / t) * tf.tphi(k0, t)),
                    ('V', 3 * ta * tf.tv(k0, t, 1) + 3 * dta * tf.tv(k0, t))):
        dd = (X1[nm] - X0[nm])[0]
        print('C1 K3c with shifted tadpoles: d%s3c - expected: max %.3e (scale %.3e)' % (nm, np.max(np.abs(dd - exp)[t > 0.3]), np.max(np.abs(exp))))
