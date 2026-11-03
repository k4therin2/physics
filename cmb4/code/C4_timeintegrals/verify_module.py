"""verify_module.py -- accuracy of the fast double-precision evaluators of second_order_closed.py against a 50-digit
reference: the closed forms evaluated in mpmath with exact (mpf) inputs.  The reference itself is checked against the
direct quadrature arbiter int_0^eta G S deta' at a subset of points.  Includes near-collinear configurations."""
import sys, time, numpy as np, sympy as sp, mpmath as mp, pickle
sys.path.insert(0, '.'); sys.path.insert(0, '../C1_einstein')
import second_order_closed as SO
from trig_integrate import Cin
from kernels_api import load, closed_form, eta, k, q, p
t0 = time.time()
mp.mp.dps = 50
def cin_mp(z):
    z = abs(z); return mp.euler + mp.log(z) - mp.ci(z) if z != 0 else mp.mpf(0)
NS = {'Cin': cin_mp, 'Si': mp.si}
lam = lambda e: sp.lambdify((k, q, p, eta), e, modules=[NS, 'mpmath'])
S = load()['second']
ref = {nm: lam(getattr(SO, nm + '_expr')) for nm in ('Phi2', 'dPhi2', 'Psi2', 'delta2', 'deltaT2', 'V2', 'h2_I')}
ref['omega2'] = lam(SO.K_om); ref['vT2'] = lam(SO.K_vT)
etap = sp.Symbol('etap', positive=True)
intP = sp.lambdify((k, q, p, eta, etap), sp.expand(S['G_Phi'] * closed_form(S['K_SPhi']).subs(eta, etap)), 'mpmath')
intPd = sp.lambdify((k, q, p, eta, etap), sp.expand(sp.diff(S['G_Phi'], eta) * closed_form(S['K_SPhi']).subs(eta, etap)), 'mpmath')
intH = sp.lambdify((k, q, p, eta, etap), sp.expand(S['G_h'] * closed_form(S['K_h']).subs(eta, etap)), 'mpmath')
fK0 = lam(S['K_Phi2_0']); fgh = lam(S['ghat']); fdgh = lam(sp.diff(S['ghat'], eta))
print("lambdified (%.0fs)" % (time.time() - t0), flush=True)
rng = np.random.default_rng(11)
pts = []
for i in range(40):
    kk = rng.uniform(0.2, 3); qq = kk * np.exp(rng.uniform(np.log(0.03), np.log(20))); mm = rng.uniform(-1, 1)
    z = np.exp(rng.uniform(np.log(0.05), np.log(300)))
    pp = np.sqrt(kk**2 + qq**2 - 2 * kk * qq * mm); ee = z * np.sqrt(3) / (kk + qq + pp)
    pts.append((kk, qq, mm, ee))
for mm in (1 - 1e-6, 1 - 1e-10, -1 + 1e-9, 1.0, -1.0):
    for qq, ee in ((0.4, 6.0), (2.5, 9.0), (0.7, 1.2)):
        pts.append((1.0, qq, mm, ee))
pts.append((1.0, np.sqrt(3) / 2, 1 / np.sqrt(3), 9.0))   # q+p = sqrt3 k
worst = {}
nq = 0
for (kk, qq, mm, ee) in pts:
    K, Q, M, T = [mp.mpf(float(v)) for v in (kk, qq, mm, ee)]
    Pm = mp.sqrt(K**2 + Q**2 - 2 * K * Q * M)
    pn = float(Pm)
    collinear = abs(abs(mm) - 1) < 1e-3
    for nm in ('Phi2', 'dPhi2', 'Psi2', 'delta2', 'deltaT2', 'V2', 'h2_I', 'omega2', 'vT2'):
        if collinear and abs(abs(mm) - 1) < 1e-12 and nm in ('Phi2', 'dPhi2', 'Psi2', 'delta2', 'deltaT2', 'V2'):
            # exactly collinear: reference from the quadrature arbiter (closed form has removable 0/0)
            continue
        r = ref[nm](K, Q, Pm, T)
        if nm in ('omega2', 'vT2'): r = Q * mp.sqrt(1 - M**2) * r
        r = complex(r)
        v = complex(np.asarray(getattr(SO, nm)(kk, qq, mm, ee)).ravel()[0])
        err = abs(v - r) / max(abs(r), 1e-300)
        key = nm + (' (collinear)' if collinear else '')
        worst[key] = max(worst.get(key, 0), err)
        if err > 1e-11: print("   large err %s at k=%g q=%g mu=%.12g eta=%g: %g (num %s ref %s)" % (nm, kk, qq, mm, ee, err, v, r))
    if nq < 10 and not collinear:      # quadrature arbiter for the reference
        nq += 1
        a = mp.quad(lambda t: intP(K, Q, Pm, T, t), mp.linspace(0, T, 30), method='gauss-legendre')
        ad = mp.quad(lambda t: intPd(K, Q, Pm, T, t), mp.linspace(0, T, 30), method='gauss-legendre')
        b = mp.quad(lambda t: intH(K, Q, Pm, T, t), mp.linspace(0, T, 30), method='gauss-legendre')
        c0 = fK0(K, Q, Pm, T)
        e1 = abs(a + c0 * fgh(K, Q, Pm, T) - ref['Phi2'](K, Q, Pm, T)) / abs(ref['Phi2'](K, Q, Pm, T))
        e2 = abs(ad + c0 * fdgh(K, Q, Pm, T) - ref['dPhi2'](K, Q, Pm, T)) / abs(ref['dPhi2'](K, Q, Pm, T))
        e3 = abs(b - ref['h2_I'](K, Q, Pm, T)) / abs(ref['h2_I'](K, Q, Pm, T))
        print("  arbiter k=%.3g q/k=%.3g mu=%.3g k eta=%.4g: Phi2 %s dPhi2 %s I_h %s" % (kk, qq / kk, mm, kk * ee, mp.nstr(e1, 3), mp.nstr(e2, 3), mp.nstr(e3, 3)), flush=True)
        worst['ARBITER vs 50-digit closed form'] = max(worst.get('ARBITER vs 50-digit closed form', 0), float(max(e1, e2, e3)))
# exactly collinear: module vs quadrature
for (qq, mm, ee) in ((0.4, 1.0, 6.0), (2.5, 1.0, 9.0), (0.7, -1.0, 1.2), (0.4, -1.0, 6.0)):
    K, Q, M, T = mp.mpf(1), mp.mpf(qq), mp.mpf(mm), mp.mpf(ee); Pm = abs(K - Q * M) if abs(M) == 1 else None
    a = mp.quad(lambda t: intP(K, Q, Pm, T, t), mp.linspace(0, T, 30), method='gauss-legendre') + fK0(K, Q, Pm, T) * fgh(K, Q, Pm, T)
    v = float(SO.Phi2(1.0, qq, mm, ee)); err = abs(v - a) / abs(a)
    worst['Phi2 exactly collinear vs quadrature'] = max(worst.get('Phi2 exactly collinear vs quadrature', 0), float(err))
    print("  exactly collinear q=%g mu=%g: module %r quad %s rel %s" % (qq, mm, v, mp.nstr(a, 17), mp.nstr(err, 3)))
print("WORST relative errors of the double-precision evaluators:")
for kk_, v in worst.items(): print("   %-40s %.2e" % (kk_, v))
print("done %.0fs" % (time.time() - t0))
