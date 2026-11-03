"""verify_einstein.py -- END-TO-END CHECK of all second-order closed forms: insert the first-order plane-wave solution
zeta(x) = A e^{i q1.x} + B e^{i q2.x} and the C4 second-order fields (coefficient of A B e^{i k.x}, k = q1 + q2,
symmetrised over the two orderings) into C1's full second-order Einstein equations E^{(2)}_{mu nu} (all 10 components)
and the second-order energy-momentum conservation equations (4 components), evaluated numerically (mpmath, 50 digits;
eta-derivatives of the closed forms by mpmath.diff).  Residuals are compared with the size of the individual terms."""
import sympy as sp, mpmath as mp, pickle, sys, time
sys.path.insert(0, '.'); sys.path.insert(0, '../C1_einstein')
from trig_integrate import *
from kernels_api import load, closed_form, eta, k, q, p
t0 = time.time()
DEBUG = len(sys.argv) > 1
S = load()['second']; L = {'Cin': Cin}
mp.mp.dps = 50
def cin_mp(z):
    z = abs(z); return mp.euler + mp.log(z) - mp.ci(z) if z != 0 else mp.mpf(0)
NS = {'Cin': cin_mp, 'Si': mp.si}
P = sp.parse_expr(pickle.load(open('phi2_closed.pkl', 'rb'))['P'], local_dict=L)
Ih = sp.parse_expr(pickle.load(open('h2_closed.pkl', 'rb'))['I_h'], local_dict=L)
lam = lambda e: sp.lambdify((k, q, p, eta), e, modules=[NS, 'mpmath'])
fP = lam(P); fI = lam(Ih)
kern = {nm: lam(closed_form(S[nm])) for nm in ('K_Q', 'K_S00', 'K_S0z', 'K_om', 'K_vT')}
fK0 = lam(S['K_Phi2_0']); fgh = lam(S['ghat'])
print("lambdified (%.0fs)" % (time.time() - t0))
E = pickle.load(open('../C1_einstein/E_sympy.pkl', 'rb'))
X, Y, Z = sp.symbols('x y z')
A, B = sp.symbols('A B')
s3 = sp.sqrt(3)
def TPhi(kk): xx = kk * eta / s3; return -2 * (sp.sin(xx) - xx * sp.cos(xx)) / xx**3

def run(kk, qq, pp, ee):
    mu = (kk**2 + qq**2 - pp**2) / (2 * kk * qq)
    sn = sp.sqrt(1 - mu**2)
    q1 = (qq * sn, 0, qq * mu); kv = (0, 0, kk); q2 = tuple(kv[i] - q1[i] for i in range(3))
    E1 = sp.exp(sp.I * (q1[0] * X + q1[2] * Z)); E2 = sp.exp(sp.I * (q2[0] * X + q2[2] * Z)); E12 = sp.exp(sp.I * kk * Z)
    # first order (exact)
    def first(mag):
        T = TPhi(mag); dT = sp.diff(T, eta)
        return T, -2 * eta**2 * mag**2 * T / 3 - 2 * eta * dT - 2 * T, -eta * (eta * dT + T) / 2
    Tq, Dq, Vq = first(qq); Tp, Dp, Vp = first(pp)
    sub = {}
    args = (eta, X, Y, Z)
    sub['Phi1'] = A * Tq * E1 + B * Tp * E2; sub['Psi1'] = sub['Phi1']
    sub['delta1'] = A * Dq * E1 + B * Dp * E2
    for i, nm in enumerate(('v11', 'v21', 'v31')):
        sub[nm] = sp.I * (A * q1[i] * Vq * E1 + B * q2[i] * Vp * E2)
    # second order: undefined amplitude functions of eta
    names2 = ['Phi2', 'Psi2', 'delta2', 'v12', 'v22', 'v32', 'om12', 'om22', 'om32', 'h112', 'h122', 'h132', 'h222', 'h232', 'h332']
    F = {nm: sp.Function('F_' + nm)(eta) for nm in names2}
    for nm in names2: sub[nm] = A * B * F[nm] * E12
    # numeric amplitudes (symmetrised over orderings (q,p) and (p,q))
    kq, kp_, pp_ = [sp.N(v, 60) for v in (kk, qq, pp)]
    K, Qm, Pm = mp.mpf(str(kq)), mp.mpf(str(kp_)), mp.mpf(str(pp_))
    def sym(f): return lambda t: f(K, Qm, Pm, t) + f(K, Pm, Qm, t)
    Phi2 = lambda t: mp.re(sym(fP)(t)) + 2 * fK0(K, Qm, Pm, t) * fgh(K, Qm, Pm, t)
    Q2 = sym(kern['K_Q']); S00 = sym(kern['K_S00']); S0z = sym(kern['K_S0z'])
    dPhi2 = lambda t: mp.diff(Phi2, t)
    Psi2 = lambda t: Phi2(t) - Q2(t)
    delta2 = lambda t: t**2 / 3 * S00(t) - mp.mpf(2) / 3 * t**2 * K**2 * Phi2(t) - 2 * t * dPhi2(t) - 2 * Psi2(t)
    vz = lambda t: t * (-S0z(t) * t - 2j * t * K * dPhi2(t) - 2j * K * Psi2(t)) / 4          # C1 K_v2z_expr = i k V2
    # vector/tensor: omega_i = P_ij q1_j [K_om(q,p) - K_om(p,q)]... use explicit: P_ij q2_j = -P_ij q1_j
    qx = mp.mpf(str(sp.N(q1[0], 60)))
    omx = lambda t: qx * (kern['K_om'](K, Qm, Pm, t) - kern['K_om'](K, Pm, Qm, t))
    vTx = lambda t: qx * (kern['K_vT'](K, Qm, Pm, t) - kern['K_vT'](K, Pm, Qm, t))
    hxx = lambda t: qx**2 / 2 * mp.re(fI(K, Qm, Pm, t) + fI(K, Pm, Qm, t))
    num = {'Phi2': Phi2, 'Psi2': Psi2, 'delta2': delta2, 'v12': vTx, 'v22': lambda t: 0, 'v32': vz,
           'om12': omx, 'om22': lambda t: 0, 'om32': lambda t: 0,
           'h112': hxx, 'h122': lambda t: 0, 'h132': lambda t: 0, 'h222': lambda t: -hxx(t), 'h232': lambda t: 0, 'h332': lambda t: 0}
    te = mp.mpf(str(sp.N(ee, 60)))
    out = []
    for key in [('E', a, b, 2) for a in range(4) for b in range(a, 4)] + [('cons', a, 2) for a in range(4)]:
        e = sp.parse_expr(E[key]).xreplace({sp.Symbol('eta'): eta})
        fs = {f for f in e.atoms(sp.Function) if isinstance(f, sp.core.function.AppliedUndef)}
        e = e.xreplace({f: sp.Function('G_' + str(f.func))(*args) for f in fs})
        rep = {sp.Function('G_' + str(f.func))(*args): sub[str(f.func)] for f in fs}
        e = e.subs(rep).doit()
        e = sp.expand(e)
        c = sp.expand(e.coeff(A, 1).coeff(B, 1) / E12)
        c = sp.expand(sp.powsimp(c))
        # numeric evaluation: F's and their eta-derivatives -> symbols, then mpmath lambdify per term
        atoms = sorted([d for d in c.atoms(sp.Derivative)], key=lambda z: -z.derivative_count)
        atoms += [f for f in c.atoms(sp.Function) if str(f.func).startswith('F_')]
        syms = tuple(sp.symbols('a0:%d' % len(atoms))) if atoms else ()
        vals = []
        for a in atoms:
            if isinstance(a, sp.Derivative):
                vals.append(mp.diff(num[str(a.expr.func)[2:]], te, a.derivative_count))
            else:
                vals.append(num[str(a.func)[2:]](te))
        cs = c
        for a, s_ in zip(atoms, syms): cs = cs.xreplace({a: s_})
        assert not [f for f in cs.atoms(sp.Function) if isinstance(f, sp.core.function.AppliedUndef)]
        tot = mp.mpc(0); mx = mp.mpf(0); lin = mp.mpc(0)
        for t in sp.Add.make_args(cs):
            vv = sp.lambdify((eta,) + tuple(syms), t, 'mpmath')(te, *vals)
            tot += vv; mx = max(mx, abs(vv))
            if t.free_symbols & set(syms): lin += vv
        if DEBUG: print('      lin %s  src %s' % (mp.nstr(lin, 12), mp.nstr(tot - lin, 12)))
        out.append((key, tot, mx))
        print("   %-16s residual %s   max term %s   ratio %s" % (str(key), mp.nstr(abs(tot), 5), mp.nstr(mx, 5), mp.nstr(abs(tot) / mx if mx else 0, 3)), flush=True)
    return out

pts = [(sp.Integer(1), sp.Rational(4, 5), sp.Rational(3, 4), sp.Rational(5, 2)),
       (sp.Rational(13, 10), sp.Rational(21, 10), sp.Rational(17, 10), sp.Rational(37, 5)),
       (sp.Rational(7, 4), sp.Rational(1, 3), sp.Rational(8, 5), sp.Rational(23, 2))]
worst = 0
for kk, qq, pp, ee in pts:
    print("point k=%s q=%s p=%s eta=%s (k eta=%s):" % (kk, qq, pp, ee, sp.N(kk * ee, 4)), flush=True)
    for key, tot, mx in run(kk, qq, pp, ee):
        worst = max(worst, abs(tot) / mx if mx else 0)
print("WORST relative residual (|sum| / max|term|): %s   (%.0fs)" % (mp.nstr(worst, 3), time.time() - t0))
