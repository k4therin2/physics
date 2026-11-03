"""pw_numeric.py -- Stage B of A1: solve the plane-wave mode equations of pw_symbolic.py.
(1) Exact small-eta (Frobenius/Taylor) series of every mode amplitude, rational coefficients, used to impose the
    regular initial conditions zeta_nl^(2)(eta->0) = zeta_nl^(3)(eta->0) = 0 exactly (amplitudes C_{n,m} of the
    regular homogeneous solution ghat_m(eta) = 3 u1(m eta/sqrt3)).
(2) High-precision numerical evaluation of the exact Green's-function solutions
      Phi_{n,m}(eta) = C_{n,m} ghat_m(eta) + int_0^eta G_m(eta,eta') S_Phi(eta') deta',
      G_m = [u1(eta')u2(eta) - u1(eta)u2(eta')]/W(eta')  (m>0),  G_0 = (eta'/3)(1 - eta'^3/eta^3),
      H_{2,0}(eta) = int_0^eta (eta' - eta'^2/eta) S_H(eta') deta',
    by Chebyshev spectral integration in mpmath (dps = 40, degree NCH), on [0, L].
(3) Checks: ODE residuals, agreement with the exact series near eta = 0, the engine's conservation equations
    nabla_mu T^{mu nu} = 0 at orders 2 and 3 (not used in the solution), and the eta->0 limits.
Output: pw_fields.pkl (Chebyshev coefficients of all mode functions, constants, series), usable via load_fields().
"""
import sympy as sp, pickle, os, sys, time, random
from mpmath import mp, mpf, sin, cos, sqrt, matrix, exp, mpc
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from cheb import Cheb, Fun, nodes

eta = sp.Symbol('eta', positive=True)
gauge_c = sp.Symbol('gauge_c')
MODES = {1: [1], 2: [0, 2], 3: [1, 3]}
L_ETA = 8          # fields needed on [0, etabar_o] with etabar_o <= 8
NCH = int(os.environ.get('A1_NCH', 160))
DPS = int(os.environ.get('A1_DPS', 40))
GAUGE_C = sp.Rational(os.environ.get('A1_GAUGE_C', '0'))
KSER = int(os.environ.get('A1_KSER', 16))
OUT = os.environ.get('A1_FIELDS_OUT', os.path.join(HERE, 'pw_fields.pkl'))

def load_eqs():
    R = pickle.load(open(os.path.join(HERE, 'pw_eqs.pkl'), 'rb'))
    def dec(v):
        if isinstance(v, str): return sp.parse_expr(v)
        if isinstance(v, dict): return {kk: dec(vv) for kk, vv in v.items()}
        if isinstance(v, list): return [dec(vv) for vv in v]
        return v
    return {kk: dec(vv) for kk, vv in R.items()}

def amp(field, n, m):
    return sp.Function('%s%d_%d' % (field, n, m))(eta)

# ------------------------------------------------------------------ evaluators
def jets_in(expr):
    """set of (funcname, derivative order) appearing in expr."""
    out = set()
    for d in expr.atoms(sp.Derivative):
        out.add((str(d.expr.func), sum(c for v, c in d.variable_count)))
    for f in expr.atoms(sp.core.function.AppliedUndef):
        out.add((str(f.func), 0))
    return out

def make_evaluator(expr, registry, extra_syms=()):
    """expr -> callable(eta_value, *extra) using Funs in registry[name] (Fun objects)."""
    expr = expr.xreplace({gauge_c: GAUGE_C})
    js = sorted(jets_in(expr))
    syms = {j: sp.Symbol('J_%s_%d' % j) for j in js}
    rep = {}
    for d in expr.atoms(sp.Derivative):
        rep[d] = syms[(str(d.expr.func), sum(c for v, c in d.variable_count))]
    e2 = expr.xreplace(rep)
    rep2 = {f: syms[(str(f.func), 0)] for f in e2.atoms(sp.core.function.AppliedUndef)}
    e2 = e2.xreplace(rep2)
    args = [eta] + [syms[j] for j in js] + list(extra_syms)
    f = sp.lambdify(args, e2, modules='mpmath')
    def ev(x, *extra):
        vals = [registry[name].val(k, x) for name, k in js]
        return f(x, *vals, *extra)
    return ev

# ------------------------------------------------------------------ exact series
def ser_of(expr, K):
    """Taylor series of a closed-form expression in eta to O(eta^K) as a polynomial (exact)."""
    return sp.expand(sp.series(expr, eta, 0, K + 1).removeO())

def laurent_trunc(expr, K):
    """expand a rational expression in eta into a Laurent polynomial and keep exponents <= K."""
    expr = sp.expand(expr)
    out = 0
    for t in sp.Add.make_args(expr):
        c, e = t.as_coeff_exponent(eta)
        if e <= K:
            out += t
    return out

def subs_series(expr, sers):
    """substitute mode functions (and derivatives) by their polynomial series."""
    rep = {}
    for d in expr.atoms(sp.Derivative):
        nm = str(d.expr.func); k = sum(c for v, c in d.variable_count)
        rep[d] = sp.diff(sers[nm], eta, k)
    e2 = expr.xreplace(rep)
    rep2 = {f: sers[str(f.func)] for f in e2.atoms(sp.core.function.AppliedUndef) if str(f.func) in sers}
    return e2.xreplace(rep2)

def particular_series(src, m, K, kind='Phi'):
    """regular particular solution of P'' + 4P'/eta + (m^2/3)P = src (kind Phi) or H'' + 2H'/eta = src (kind H),
    src a Laurent polynomial with exponents >= 0; returns polynomial up to eta^K."""
    src = sp.expand(src)
    coeffs = {}
    for t in sp.Add.make_args(src):
        if t == 0: continue
        c, e = t.as_coeff_exponent(eta)
        assert e >= -1, ('source too singular at eta=0', t)   # e=-1 arises for the gauge mode gf=c*eta; solution ~ eta, regular
        coeffs[int(e)] = coeffs.get(int(e), 0) + c
    p = {}
    for j in range(-1, K - 1):
        s = coeffs.get(j, 0)
        if kind == 'Phi':
            p[j] = (s - sp.Rational(m * m, 3) * p.get(j - 2, 0)) / ((j + 2) * (j + 5))
        else:
            p[j] = s / ((j + 2) * (j + 3))
    return sp.expand(sum(p[j] * eta ** (j + 2) for j in p))

def ghat_series(m, K):
    if m == 0:
        return sp.Integer(1)
    X = m * eta / sp.sqrt(3)
    return ser_of(3 * (sp.sin(X) - X * sp.cos(X)) / X ** 3, K)

# ------------------------------------------------------------------ main
def build(verbose=True):
    mp.dps = DPS
    t0 = time.time()
    R = load_eqs()
    first = R['first_closed']
    series = {}      # name -> polynomial in eta (exact)
    consts = {}
    K1 = KSER
    series['Phi1_1'] = ser_of(first['Phi'], K1)
    series['delta1_1'] = ser_of(first['delta'], K1)
    series['V1_1'] = ser_of(first['V'], K1)
    if verbose: print('first-order series:', series['Phi1_1'], '|', series['delta1_1'], '|', series['V1_1'])

    # ---- numerical first-order Funs on [0, L] from the closed forms (evaluated at high precision)
    xs = nodes(0, L_ETA, NCH)
    reg = {}
    with mp.workdps(DPS + 40):
        fl = {nm: sp.lambdify(eta, first[key], modules='mpmath') for nm, key in (('Phi1_1', 'Phi'), ('delta1_1', 'delta'), ('V1_1', 'V'))}
        vals = {nm: [f(x) for x in xs] for nm, f in fl.items()}
    for nm in vals:
        reg[nm] = Fun(Cheb.fit(0, L_ETA, [mpf(v) for v in vals[nm]]))
    # accuracy of the Chebyshev representation incl. derivatives vs closed forms
    with mp.workdps(DPS + 40):
        dfl = {nm: [sp.lambdify(eta, sp.diff(first[key], eta, k), modules='mpmath') for k in range(4)] for nm, key in (('Phi1_1', 'Phi'), ('delta1_1', 'delta'), ('V1_1', 'V'))}
        err = 0
        for _ in range(10):
            x = mpf(random.random()) * L_ETA
            for nm in dfl:
                for k in range(4):
                    err = max(err, abs(reg[nm].val(k, x) - dfl[nm][k](x)))
    if verbose: print('first-order Chebyshev representation: max error (f..f\'\'\') = %.1e  tail %.1e (%.1fs)' % (err, reg['Phi1_1'].ch.tail(), time.time() - t0))
    reg['Psi1_1'] = reg['Phi1_1']
    series['Psi1_1'] = series['Phi1_1']

    # Green's function ingredients
    def u12(m):
        X = m * eta / sp.sqrt(3)
        u1 = (sp.sin(X) - X * sp.cos(X)) / X ** 3; u2 = (sp.cos(X) + X * sp.sin(X)) / X ** 3
        W = sp.simplify(u1 * sp.diff(u2, eta) - sp.diff(u1, eta) * u2)
        assert sp.simplify(W + 3 * sp.sqrt(3) / (eta ** 4 * m ** 3)) == 0, W
        return u1, u2, W

    checks = {}
    for n in (2, 3):
        Kn = KSER - 2 * (n - 1)
        for m in MODES[n]:
            res = R[(n, m)]
            P = amp('Phi', n, m)
            # ---------- exact series
            src_ser = laurent_trunc(subs_series(res['SPhi'].xreplace({gauge_c: GAUGE_C}), series), Kn - 2)
            Pp = particular_series(src_ser, m, Kn, 'Phi')
            C = sp.Symbol('C_%d_%d' % (n, m))
            Pser = sp.expand(C * ghat_series(m, Kn) + Pp)
            sers_m = dict(series); sers_m['Phi%d_%d' % (n, m)] = Pser
            if m == 0:
                sH = laurent_trunc(subs_series(res['SH'], series), Kn - 2)
                sers_m['H%d_%d' % (n, m)] = particular_series(sH, 0, Kn, 'H')
            for key, nm in (('Psi', 'Psi'), ('delta', 'delta'), ('V', 'V')):
                if key in res and (m != 0 or key != 'V'):
                    sers_m['%s%d_%d' % (nm, n, m)] = laurent_trunc(subs_series(res[key].xreplace({gauge_c: GAUGE_C}), sers_m), Kn)
            zeta_m = laurent_trunc(subs_series(R[('zeta', n)][m], sers_m), Kn)
            z0 = zeta_m.subs(eta, 0)     # constant term (negative powers would signal an error)
            assert all(t.as_coeff_exponent(eta)[1] >= 0 for t in sp.Add.make_args(zeta_m)), zeta_m
            Csol = sp.solve(z0, C)[0]
            consts[(n, m)] = Csol
            for kk in list(sers_m):
                sers_m[kk] = sp.expand(sers_m[kk].subs(C, Csol))
            series.update({kk: vv for kk, vv in sers_m.items() if kk not in series})
            zeta_check = sp.expand(zeta_m.subs(C, Csol))
            if verbose:
                print('order %d mode %d: C = %s ;  zeta^(%d)_%d(eta) series after imposing = %s' % (n, m, Csol, n, m, sp.series(zeta_check, eta, 0, 5)))
                print('   Phi%d_%d series: %s' % (n, m, sp.series(series['Phi%d_%d' % (n, m)], eta, 0, 5)))
            # ---------- numerical Green's solution
            evS = make_evaluator(res['SPhi'], reg)
            Svals = [evS(x) for x in xs]
            if m > 0:
                u1, u2, W = u12(m)
                with mp.workdps(DPS + 30):
                    f_u1 = sp.lambdify(eta, u1, 'mpmath'); f_u2 = sp.lambdify(eta, u2, 'mpmath'); f_W = sp.lambdify(eta, W, 'mpmath')
                    u1v = [mpf(f_u1(x)) for x in xs]; u2v = [mpf(f_u2(x)) for x in xs]; Wv = [mpf(f_W(x)) for x in xs]
                A = Cheb.fit(0, L_ETA, [u1v[i] * Svals[i] / Wv[i] for i in range(len(xs))]).antideriv(0)
                B = Cheb.fit(0, L_ETA, [u2v[i] * Svals[i] / Wv[i] for i in range(len(xs))]).antideriv(0)
                Cm = mpf(sp.N(Csol, DPS + 5))
                Pvals = [u2v[i] * A(xs[i]) - u1v[i] * B(xs[i]) + Cm * 3 * u1v[i] for i in range(len(xs))]
            else:
                A0 = Cheb.fit(0, L_ETA, [xs[i] * Svals[i] for i in range(len(xs))]).antideriv(0)
                B0 = Cheb.fit(0, L_ETA, [xs[i] ** 4 * Svals[i] for i in range(len(xs))]).antideriv(0)
                Cm = mpf(sp.N(Csol, DPS + 5))
                Pvals = [(A0(x) - B0(x) / x ** 3) / 3 + Cm for x in xs]
            reg['Phi%d_%d' % (n, m)] = Fun(Cheb.fit(0, L_ETA, Pvals))
            if m == 0:
                evH = make_evaluator(res['SH'], reg)
                Hv = [evH(x) for x in xs]
                AH = Cheb.fit(0, L_ETA, [xs[i] * Hv[i] for i in range(len(xs))]).antideriv(0)
                BH = Cheb.fit(0, L_ETA, [xs[i] ** 2 * Hv[i] for i in range(len(xs))]).antideriv(0)
                reg['H%d_%d' % (n, m)] = Fun(Cheb.fit(0, L_ETA, [AH(x) - BH(x) / x for x in xs]))
            for key, nm in (('Psi', 'Psi'), ('delta', 'delta'), ('V', 'V')):
                if key in res and (m != 0 or key != 'V'):
                    ev = make_evaluator(res[key], reg)
                    reg['%s%d_%d' % (nm, n, m)] = Fun(Cheb.fit(0, L_ETA, [ev(x) for x in xs]))
            # ---------- checks: ODE residual, series agreement, Chebyshev tails
            Pf = reg['Phi%d_%d' % (n, m)]
            resid = 0
            for _ in range(8):
                x = mpf('0.3') + mpf(random.random()) * (L_ETA - mpf('0.3'))
                resid = max(resid, abs(Pf.val(2, x) + 4 * Pf.val(1, x) / x + mpf(m * m) / 3 * Pf(x) - evS(x)))
            if m == 0:
                Hf = reg['H%d_%d' % (n, m)]
                for _ in range(8):
                    x = mpf('0.3') + mpf(random.random()) * (L_ETA - mpf('0.3'))
                    resid = max(resid, abs(Hf.val(2, x) + 2 * Hf.val(1, x) / x - evH(x)))
            sdiff = 0
            for kk in [k_ for k_ in sers_m if k_.endswith('%d_%d' % (n, m))]:
                f_ser = sp.lambdify(eta, series[kk], 'mpmath')
                for xv in (mpf('0.02'), mpf('0.05')):
                    sdiff = max(sdiff, abs(reg[kk](xv) - f_ser(xv)))
            tails = max(reg[kk].ch.tail() for kk in sers_m if kk.endswith('%d_%d' % (n, m)))
            checks[(n, m)] = dict(ode_resid=resid, series_diff=sdiff, cheb_tail=tails)
            if verbose: print('   numerics: ODE residual %.1e, |num - series| at eta<=0.05: %.1e, Chebyshev tail %.1e (%.1fs)' % (resid, sdiff, tails, time.time() - t0))

    # ---------- conservation equations (engine, orders 2,3) at random (eta, x): must vanish
    zsym = sp.Symbol('z')
    cons_res = {}
    for n in (2, 3):
        worst = 0
        for nu in range(4):
            for m, ex in R[('cons', n)][nu].items():
                ev = make_evaluator(ex, reg)
                for _ in range(5):
                    x = mpf('0.5') + mpf(random.random()) * 7
                    worst = max(worst, abs(ev(x)))
        cons_res[n] = worst
        if verbose: print('conservation equations nabla_mu T^{mu nu} at order %d (all modes, random eta): max |.| = %.1e' % (n, worst))

    # ---------- save
    out = {'cheb': {kk: (str(f.ch.a), str(f.ch.b), [str(c) for c in f.ch.c]) for kk, f in reg.items()},
           'consts': {kk: sp.srepr(v) for kk, v in consts.items()}, 'series': {kk: sp.srepr(v) for kk, v in series.items()},
           'checks': checks, 'cons_res': cons_res, 'DPS': DPS, 'NCH': NCH, 'GAUGE_C': sp.srepr(GAUGE_C), 'L': L_ETA,
           'uhat0': {n: {m: sp.srepr(v) for m, v in R[('uhat0', n)].items()} for n in (1, 2, 3)},
           'first_closed': {kk: sp.srepr(v) for kk, v in first.items()}}
    with open(OUT, 'wb') as fh:
        pickle.dump(out, fh)
    if verbose: print('saved', OUT, '(%.1fs)' % (time.time() - t0))
    return reg, consts, series, checks

def load_fields(fn=None):
    fn = fn or OUT
    d = pickle.load(open(fn, 'rb'))
    reg = {kk: Fun(Cheb(mpf(a), mpf(b), [mpf(c) for c in cs])) for kk, (a, b, cs) in d['cheb'].items()}
    consts = {kk: sp.parse_expr(v) for kk, v in d['consts'].items()}
    series = {kk: sp.parse_expr(v) for kk, v in d['series'].items()}
    uh = {n: {m: sp.parse_expr(v) for m, v in dd.items()} for n, dd in d['uhat0'].items()}
    return dict(reg=reg, consts=consts, series=series, uhat0=uh, meta=d)

if __name__ == '__main__':
    build()
