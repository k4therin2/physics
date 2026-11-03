"""pw_symbolic.py -- Stage A of A1: plane-wave (planar) reduction of C1's Einstein+fluid equations.

Ansatz (zeta = eps cos x, k = 1, all fields functions of (eta, x) only, y-z rotation symmetry, x -> -x parity):
  scalars  X_n(eta,x) = sum_{m in M_n} X_n_m(eta) cos(m x),   M_1={1}, M_2={0,2}, M_3={1,3}
  velocity v^x_n = d_x V_n,  V_n = sum_m V_n_m(eta) cos(m x);  v^y = v^z = 0;  omega_i = 0;
  tensor:  only the k=0 homogeneous traceless mode at second order: h_11 = H2_0(eta), h_22 = h_33 = -H2_0/2
           (mode m != 0 tensors vanish by transversality + y-z symmetry; verified below that the ansatz closes).
  Psi_1 = Phi_1 (first order), Psi_2_0 = Phi_2_0 + gf(eta) (k=0 residual time-reparametrisation gauge; gf = gauge_c*eta).
The engine's E_{mu nu} (coefficient of eps^n, including 1/n!) is converted with a custom jet callback into a Laurent
polynomial in z = e^{ix}; the coefficient of z^m gives the mode-m equation.  Per mode we solve, exactly as C1 does:
  traceless (E_xx - E_yy)  -> Psi (m != 0)   or the ODE for H (m = 0)
  00                       -> delta           0x -> V            trace -> ODE  Phi'' + 4 Phi'/eta + (m^2/3) Phi = S_Phi
Output: pw_eqs.pkl with, per (n, m): 'Psi','delta','V' (algebraic in terms of Phi_n_m, Phi_n_m' and lower orders),
'SPhi' (ODE source), 'SH' (m=0), plus zeta^(n), uhat^0, cons^nu as z-polynomials, and the first-order closed forms.
"""
import sympy as sp, sys, os, pickle, time
HERE = os.path.dirname(os.path.abspath(__file__))
C1 = os.path.join(HERE, '..', 'C1_einstein')
sys.path.insert(0, C1)
from tpoly import load_polys, REG, Poly, field_series, const

t0 = time.time()
eta = sp.Symbol('eta', positive=True)
z = sp.Symbol('z')
gauge_c = sp.Symbol('gauge_c')
I = sp.I
EVEN = {'Phi', 'Psi', 'delta', 'h11', 'h22', 'h33'}
ODD = {'v1'}
ZERO = {'v2', 'v3', 'om1', 'om2', 'om3', 'h12', 'h13', 'h23'}
MODES = {1: [1], 2: [0, 2], 3: [1, 3]}

def amp(field, n, m):
    return sp.Function('%s%d_%d' % (field, n, m))(eta)

def trig_z(m, kind, nx):
    """d_x^nx cos(mx) or sin(mx) as a Laurent polynomial in z = e^{ix}."""
    if kind == 'cos':
        return ((I * m) ** nx * z ** m + (-I * m) ** nx * z ** (-m)) / 2
    return ((I * m) ** nx * z ** m - (-I * m) ** nx * z ** (-m)) / (2 * I)

def gf_expr():
    return gauge_c * eta

def jet_value(field, n, mi):
    """Fourier image of d^mi X_n for the ansatz (mi = (nt, nx, ny, nz))."""
    nt, nx, ny, nz = mi
    if ny or nz or field in ZERO:
        return sp.Integer(0)
    if field in EVEN:
        if field == 'Psi':
            if n == 1:
                base = {1: amp('Phi', 1, 1)}
            elif n == 2:
                base = {0: amp('Phi', 2, 0) + gf_expr(), 2: amp('Psi', 2, 2)}
            else:
                base = {m: amp('Psi', n, m) for m in MODES[n]}
        elif field in ('h11', 'h22', 'h33'):
            if n != 2:
                return sp.Integer(0)
            H = amp('H', 2, 0)
            base = {0: H if field == 'h11' else -H / 2}
        else:
            base = {m: amp(field, n, m) for m in MODES[n]}
        tot = 0
        for m, A in base.items():
            tot += sp.diff(A, eta, nt) * trig_z(m, 'cos', nx)
        return tot
    # v1 = d_x V
    tot = 0
    for m in MODES[n]:
        if m == 0:
            continue
        tot += sp.diff(amp('V', n, m), eta, nt) * trig_z(m, 'cos', nx + 1)
    return tot

def callback(field, od, mi):
    return jet_value(field, od, mi)

def to_z(poly, n):
    """coefficient of eps^n of a Poly as expanded Laurent polynomial in z (times n! so that it is the n-th order equation E^{(n)})."""
    ex = poly.order_part(n).to_sympy(eta=eta, funcs=callback)
    return sp.expand(ex * sp.factorial(n))

def modes_of(ex):
    """dict m -> coefficient of z^m (m >= 0 kept, with the real-amplitude normalisation: A_m = 2 c_m (m>0), A_0 = c_0;
    for odd components the coefficient c_m of z^m equals B_m/(2i))."""
    ex = sp.expand(ex)
    out = {}
    if ex == 0:
        return out
    P = sp.Poly(ex * z ** 8, z)
    for (e,), c in P.terms():
        out[e - 8] = c
    return out

print('loading E.pkl ...')
d = load_polys(os.path.join(C1, 'E.pkl'))
E = d['E']; cons = d['cons']; uhat0 = d['uhat0']
comp = {'00': (0, 0), '0x': (0, 1), '0y': (0, 2), '0z': (0, 3), 'xx': (1, 1), 'yy': (2, 2), 'zz': (3, 3), 'xy': (1, 2), 'xz': (1, 3), 'yz': (2, 3)}

def unknowns(n):
    names = []
    for m in MODES[n]:
        names += [amp('Phi', n, m), amp('delta', n, m)]
        if m != 0:
            names += [amp('Psi', n, m), amp('V', n, m)]
        else:
            names += [amp('H', n, m)]
    return names

def collect_mode_eqs(n):
    """Return dict name -> {m: expr} of E^{(n)} components (z-coefficients)."""
    out = {}
    for nm, (a, b) in comp.items():
        ex = to_z(E[a][b], n)
        out[nm] = modes_of(ex)
        print('  E_%s order %d: modes %s  (%.1fs)' % (nm, n, sorted(out[nm].keys()), time.time() - t0))
    return out

def lin_solve(ex, unk):
    """solve linear expression ex == 0 for the AppliedUndef unk (treated as a symbol)."""
    s = sp.Dummy()
    ex2 = ex.xreplace({unk: s})
    sol = sp.solve(ex2, s)
    assert len(sol) == 1, (unk, sol)
    return sp.expand(sol[0])

def ode_form(ex, P, order_coeffs):
    """ex linear in P, P', P''; return (alpha, source) with ex = alpha*(P'' + c1 P' + c0 P) - alpha*source."""
    s0, s1, s2 = sp.Dummy(), sp.Dummy(), sp.Dummy()
    ex2 = ex.xreplace({sp.Derivative(P, (eta, 2)): s2}).xreplace({sp.Derivative(P, eta): s1}).xreplace({P: s0})
    ex2 = sp.expand(ex2)
    a2 = ex2.coeff(s2); a1 = ex2.coeff(s1); a0 = ex2.coeff(s0)
    rest = sp.expand(ex2 - a2 * s2 - a1 * s1 - a0 * s0)
    assert not rest.has(s0) and not rest.has(s1) and not rest.has(s2)
    c1, c0 = order_coeffs
    assert sp.simplify(a1 / a2 - c1) == 0, (sp.simplify(a1 / a2), c1)
    assert sp.simplify(a0 / a2 - c0) == 0, (sp.simplify(a0 / a2), c0)
    return a2, sp.expand(-rest / a2)

results = {'MODES': MODES}
for n in (1, 2, 3):
    print('=== order', n)
    Em = collect_mode_eqs(n)
    # closure checks: components that must vanish identically
    for nm in ('0y', '0z', 'xy', 'xz', 'yz'):
        assert all(sp.expand(v) == 0 for v in Em[nm].values()), ('nonzero', nm, n, Em[nm])
    assert all(sp.expand(Em['yy'].get(m, 0) - Em['zz'].get(m, 0)) == 0 for m in set(Em['yy']) | set(Em['zz'])), 'yy != zz'
    # parity: even components have no odd modes etc.
    for nm in ('00', 'xx', 'yy', 'zz', '0x'):
        for m, v in Em[nm].items():
            if (m % 2) != (n % 2):
                assert sp.expand(v) == 0, ('parity violation', nm, n, m)
    # odd component 0x: coefficient of z^m and z^-m must be opposite; 00: equal
    for m in MODES[n]:
        if m > 0:
            assert sp.expand(Em['0x'].get(m, 0) + Em['0x'].get(-m, 0)) == 0
            assert sp.expand(Em['00'].get(m, 0) - Em['00'].get(-m, 0)) == 0
    print('  closure/parity checks passed (%.1fs)' % (time.time() - t0))
    for m in MODES[n]:
        P = amp('Phi', n, m)
        e00 = Em['00'].get(m, 0); e0x = Em['0x'].get(m, 0)
        eTL = sp.expand(Em['xx'].get(m, 0) - Em['yy'].get(m, 0)); etr = sp.expand(Em['xx'].get(m, 0) + Em['yy'].get(m, 0) + Em['zz'].get(m, 0))
        res = {}
        if m != 0:
            S = amp('Psi', n, m); D = amp('delta', n, m); V = amp('V', n, m)
            if n == 1:
                Ssol = P     # Psi1 = Phi1 built into the ansatz; check traceless eq vanishes identically
                assert sp.expand(eTL) == 0, eTL
            else:
                Ssol = lin_solve(eTL, S)
                assert not Ssol.has(sp.Derivative(P, eta)), 'Psi solution should be algebraic in Phi'
            Dsol = lin_solve(e00.xreplace({S: Ssol}) if n > 1 else e00, D)
            Vsol = lin_solve(e0x.xreplace({S: Ssol}) if n > 1 else e0x, V)
            tr = etr
            if n > 1:
                tr = tr.xreplace({sp.Derivative(S, (eta, 2)): sp.diff(Ssol, eta, 2), sp.Derivative(S, eta): sp.diff(Ssol, eta), S: Ssol})
            tr = tr.xreplace({sp.Derivative(D, eta): sp.diff(Dsol, eta), D: Dsol})
            tr = tr.xreplace({sp.Derivative(V, eta): sp.diff(Vsol, eta), V: Vsol})
            tr = sp.expand(tr.doit())
            assert not tr.has(S) and not tr.has(D) and not tr.has(V)
            alpha, src = ode_form(tr, P, (4 / eta, sp.Rational(m * m, 3)))
            res = {'Psi': Ssol, 'delta': Dsol, 'V': Vsol, 'SPhi': src, 'alpha': alpha}
            print('  mode %d: Phi ODE operator verified (alpha = %s); SPhi %d terms; Psi %d, delta %d, V %d terms' % (
                m, alpha, len(sp.Add.make_args(src)), len(sp.Add.make_args(Ssol)), len(sp.Add.make_args(Dsol)), len(sp.Add.make_args(Vsol))))
        else:
            D = amp('delta', n, 0); H = amp('H', n, 0)
            assert sp.expand(e0x) == 0, ('0x zero mode nonzero', e0x)
            # H ODE from the traceless eq.
            assert not eTL.has(P) and not eTL.has(D), 'traceless zero mode should not contain Phi/delta'
            alphaH, srcH = ode_form(eTL, H, (2 / eta, 0))
            Dsol = lin_solve(e00, D)
            assert not e00.has(H), '00 should not contain H'
            tr = etr.xreplace({sp.Derivative(D, eta): sp.diff(Dsol, eta), D: Dsol})
            tr = sp.expand(tr.doit())
            assert not tr.has(H), 'trace should not contain H'
            alpha, src = ode_form(tr, P, (4 / eta, 0))
            res = {'delta': Dsol, 'SPhi': src, 'alpha': alpha, 'SH': srcH, 'alphaH': alphaH, 'Psi': P + gf_expr()}
            print('  mode 0: Phi ODE verified (alpha = %s), H ODE verified (alphaH = %s); SPhi %d terms, SH %d terms' % (
                alpha, alphaH, len(sp.Add.make_args(src)), len(sp.Add.make_args(srcH))))
        results[(n, m)] = res
    # conservation equations as z-polynomials (for numerical verification later)
    results[('cons', n)] = [modes_of(to_z(cons[nu], n)) for nu in range(4)]
    results[('uhat0', n)] = modes_of(to_z(uhat0, n) / sp.factorial(n))      # coefficient of eps^n (no n!)
    print('  cons/uhat0 extracted (%.1fs)' % (time.time() - t0))

# zeta_nl = (1/4) log(1+delta) + (1/2) log(1-2Phi): order parts as z-polynomials (times n!)
Phi = field_series('Phi'); delta = field_series('delta')
zeta_nl = (const(1) + delta).log1p_series() / 4 + (const(1) - Phi * 2).log1p_series() / 2
for n in (1, 2, 3):
    results[('zeta', n)] = modes_of(to_z(zeta_nl, n))
    print('zeta^(%d) modes:' % n, {m: len(sp.Add.make_args(v)) for m, v in results[('zeta', n)].items()})

# first-order closed forms (derived here from the mode-1 equations: Phi1_1 = T_Phi solves the ODE with zero source)
X = eta / sp.sqrt(3)
TPhi = -2 * (sp.sin(X) - X * sp.cos(X)) / X ** 3
assert sp.simplify(results[(1, 1)]['SPhi']) == 0
assert sp.simplify(sp.diff(TPhi, eta, 2) + 4 / eta * sp.diff(TPhi, eta) + TPhi / 3) == 0
P1 = amp('Phi', 1, 1)
Tdel = sp.simplify(results[(1, 1)]['delta'].xreplace({sp.Derivative(P1, eta): sp.diff(TPhi, eta)}).xreplace({P1: TPhi}).doit())
TV = sp.simplify(results[(1, 1)]['V'].xreplace({sp.Derivative(P1, eta): sp.diff(TPhi, eta)}).xreplace({P1: TPhi}).doit())
print('T_Phi =', TPhi); print('T_delta =', Tdel); print('T_V =', TV)
# compare with C1's first_order.pkl
fo = {kk: (sp.parse_expr(vv) if isinstance(vv, str) else vv) for kk, vv in pickle.load(open(os.path.join(C1, 'first_order.pkl'), 'rb')).items()}
eta_c1 = sp.Symbol('eta', positive=True)
for mine, theirs in ((Tdel, fo['T_delta_eta']), (TV, fo['T_V_eta'])):
    diff = sp.simplify(mine - theirs.subs({sp.Symbol('k', positive=True): 1, eta_c1: eta}))
    print('   difference to C1:', diff)
    assert diff == 0
results['first_closed'] = {'Phi': TPhi, 'delta': Tdel, 'V': TV}
results['zeta1_check'] = sp.limit(Tdel / 4 - TPhi, eta, 0)
print('lim (delta_T1 - Phi1) =', results['zeta1_check'])

with open(os.path.join(HERE, 'pw_eqs.pkl'), 'wb') as fh:
    def enc(v):
        if isinstance(v, sp.Basic): return sp.srepr(v)
        if isinstance(v, dict): return {kk: enc(vv) for kk, vv in v.items()}
        if isinstance(v, list): return [enc(vv) for vv in v]
        return v
    pickle.dump({kk: enc(vv) for kk, vv in results.items()}, fh)
print('saved pw_eqs.pkl (%.1fs)' % (time.time() - t0))
