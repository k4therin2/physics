"""check_symbolic.py -- numerical evaluation of the named-object kernels of symbolic.py vs the engine (kernel2) channel by channel."""
import os, sys, pickle
import numpy as np
import sympy as sp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from symbolic import channel_kernels, k1, k2, mu, mu1, mu2, K, Kn, etas, etao, Pperp, E, EK, second, S, ETA
from kernel2 import Geometry, MapKernels, Pair
from transfer import tphi, tdt, tv
g = Geometry(2.0, 8.0)
rng = np.random.default_rng(3)
nv = rng.normal(size=(6, 3)); nv /= np.linalg.norm(nv, axis=1)[:, None]
kv1 = np.array([0.7, -0.2, 0.4]); kv2 = np.array([-0.3, 0.9, 0.5])
mk = MapKernels(g, kv1, kv2, nv, form='direct')
Ktot, ch_eng = mk.K2()
pr = mk.pair
q1, q2 = np.linalg.norm(kv1), np.linalg.norm(kv2); Kvec = kv1 + kv2
geo = {k1: q1, k2: q2, mu: kv1 @ kv2 / q1 / q2, K: np.linalg.norm(Kvec), etas: g.eta_s, etao: g.eta_o}
m1 = nv @ kv1 / q1; m2 = nv @ kv2 / q2
def values(where):
    T = {'ray': g.t, 'obs': np.array([g.eta_o]), 'emi': np.array([g.eta_s])}[where]
    chi = g.eta_o - T
    v = dict(geo)
    v[mu1] = m1[:, None]; v[mu2] = m2[:, None]; v[Kn] = (q1 * m1 + q2 * m2)[:, None]
    v[Pperp] = (q1 * q2 * (geo[mu] - m1 * m2))[:, None]
    for s, (q, mm) in ((1, (q1, m1)), (2, (q2, m2))):
        Es = np.exp(1j * q * mm[:, None] * chi[None, :]); v[E[s]] = Es
        Eray = np.exp(1j * q * mm[:, None] * g.chi[None, :])
        Kf = tv(q, g.eta_o) + 2 * g.cum(tphi(q, g.t)[None, :] * Eray)
        Xo = np.sum(g.ww * tv(q, g.tw)); deo = -np.sum(g.ww * g.tw * tphi(q, g.tw)) / g.eta_o
        Dp = Xo - g.cum(Kf); Dz = 1j * q * mm[:, None] * Xo + deo + 2 * g.cum(tphi(q, g.t)[None, :] * Eray)
        pick = {'ray': slice(None), 'obs': slice(-1, None), 'emi': slice(0, 1)}[where]
        v[S('Kf', s)] = Kf[:, pick]; v[S('Dp', s)] = Dp[:, pick]; v[S('Dz', s)] = Dz[:, pick]
        v[S('Xo', s)] = Xo; v[S('detao1', s)] = deo
        v[S('detae1', s)] = (g.eta_s * tdt(q, g.eta_s) * np.exp(1j * q * mm * (g.eta_o - g.eta_s)))[:, None]
        for nm, f in (('TPhi', tphi), ('TdT', tdt), ('TV', tv)):
            v[S(nm, s)] = f(q, T, 0)[None, :]; v[S(nm + 'p', s)] = f(q, T, 1)[None, :]
    v[EK] = v[E[1]] * v[E[2]]
    pick = {'ray': slice(None), 'obs': slice(-1, None), 'emi': slice(0, 1)}[where]
    for nm, arr in (('Phi2', pr.Phi2), ('Phi2p', pr.dPhi2), ('Psi2', pr.Psi2), ('Psi2p', pr.dPsi2), ('dT2', pr.dT2), ('V2', pr.V2),
                    ('KvT12', pr.KvT), ('Kom12', pr.Kom), ('H', pr.H), ('Hp', pr.dH)):
        v[second[nm]] = arr[None, pick]
    from transfer import KER
    v[second['KvT21']] = KER['K_vT'](pr.K, q2, q1, g.t)[None, pick]
    v[second['Kom21']] = KER['K_om'](pr.K, q2, q1, g.t)[None, pick]
    return v
chs = channel_kernels()
worst = 0
for name, (where, expr) in chs.items():
    v = values(where)
    syms = sorted(expr.free_symbols, key=str)
    f = sp.lambdify(syms, expr, 'numpy')
    val = f(*[v[s_] for s_ in syms]) * np.ones((len(nv), 1))
    val = g.integ(val) if where == 'ray' else val[:, 0]
    eng = ch_eng.get(name)
    if eng is None:   # engine splits 2nd-order v into V2/v2T sub-channels
        eng = sum(vv for kk, vv in ch_eng.items() if kk.startswith('[%s] 2nd:' % where))
    d = np.max(np.abs(val - eng)) / max(1e-12, np.max(np.abs(eng)))
    worst = max(worst, d)
    print('%-70s rel.diff %.1e' % (name, d))
print('WORST', worst)
