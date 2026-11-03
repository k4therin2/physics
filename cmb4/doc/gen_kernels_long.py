"""gen_kernels_long.py (W4): typesets the long machine-generated expressions of App. kernels as
(prefactor) x [sum of terms], broken over lines (multline*).  Pure typesetting of C1's stored sympy expressions; prints LaTeX blocks.
Tadpoles are written in x = q eta/sqrt3."""
import pickle, re, sympy as sp
C1 = '/home/user/physics/cmb4/code/C1_einstein/'
so = pickle.load(open(C1 + 'second_order.pkl', 'rb')); ic = pickle.load(open(C1 + 'third_order_ic.pkl', 'rb'))
S = lambda v: sp.sympify(v) if isinstance(v, str) else v
def vis(t):  # crude visible width of a latex term
    t = re.sub(r'\\(left|right|cdot)', '', t); t = re.sub(r'\\[a-zA-Z]+', 'x', t); t = re.sub(r'[{}^_ ]', '', t)
    return len(t)
def tidy(t):
    for f in ('sin', 'cos'):
        t = t.replace('\\%s^{2}{\\left(x \\right)}' % f, '\\%s^2x' % f).replace('\\%s{\\left(x \\right)}' % f, '\\%s x' % f)
    return t
def block(lhs, expr, names=None, width=50):
    expr = sp.together(sp.expand(expr))
    num, den = sp.fraction(expr)
    num = sp.expand(num)
    c, num2 = sp.factor_terms(num).as_coeff_Mul() if False else (1, num)
    terms = sp.Add.make_args(num)
    g = sp.gcd_list([sp.Rational(sp.Poly(t, *t.free_symbols).LC()) if False else sp.factor_terms(t).as_coeff_Mul()[0] for t in terms])
    pre = sp.nsimplify(g) / den
    terms = sorted(terms, key=lambda t: (-sp.Poly(t, x, sp.sin(x), sp.cos(x)).degree(x) if t.has(x) else 0, str(t)))
    lt = [tidy(sp.latex(sp.expand(t / g), symbol_names=names or {})) for t in terms]
    lines, cur, w = [], '', 0
    for i, t in enumerate(lt):
        s = t if (i == 0 or t.startswith('-')) else '+' + t
        if w + vis(s) > width and cur:
            lines.append(cur); cur, w = '', 0
        cur += s; w += vis(s)
    lines.append(cur)
    body = ('\\\\\n').join(lines)
    return '\\begin{multline*}\n%s=%s\\Big[%s\\Big]\n\\end{multline*}' % (lhs, sp.latex(pre, symbol_names=names or {}), body)
eta, q, k, p, x = sp.symbols('eta q k p x', positive=True)
out = {}
tads = ic['tads']
for key, lab in [('tadQ2', r'\mathrm{tad}_{Q_2}'), ('tadS00', r'\mathrm{tad}_{S_{00}}'), ('tadSPhi', r'\mathrm{tad}_{S_\Phi}'),
                 ('tadPhi2', r'\mathrm{tad}_{\Phi_2}'), ('taddelta2', r'\mathrm{tad}_{\delta_2}')]:
    e = S(tads[key])
    e = e.subs({sp.Symbol(str(s), positive=True): s for s in []})
    syms = {s.name: s for s in e.free_symbols}
    e = e.subs(syms['eta'], sp.sqrt(3) * x / syms['q'])
    e = sp.expand_trig(sp.simplify(e))
    e = e.rewrite(sp.cos) if False else e
    e0 = S(tads[key]); s0 = {v.name: v for v in e0.free_symbols}
    for qq, ee in ((0.7, 1.3), (2.1, 0.4)):
        a = complex(e0.subs({s0['q']: qq, s0['eta']: ee})); b = complex(e.subs({syms['q']: qq, x: qq * ee / 3**0.5}))
        assert abs(a - b) < 1e-9 * (1 + abs(a)), (key, a, b)
    out[key] = block(lab, sp.expand(e))
e = S(ic['C3_inv']); out['C3inv'] = block(r'C_3', e)
e = S(so['K_SPhi']); out['KSPhi'] = None
open('/tmp/claude-0/-home-user-physics/fb9ccdc3-37b3-5823-8617-3fdb5adcb1f8/scratchpad/kl.pkl', 'wb').write(pickle.dumps(out))
for k_, v in out.items():
    if v: print('%%%%', k_); print(v)
# ---- K_SPhi grouped by the four transfer-function structures
e = S(so['K_SPhi'])
Tq_, Tp_, dTq, dTp = sp.symbols('T_q T_p Tdq Tdp')
syms = {s.name: s for s in e.free_symbols}
ee, qq, pp = syms['eta'], syms['q'], syms['p']
fq = sp.Function('T_Phi')(qq, ee); fp = sp.Function('T_Phi')(pp, ee)
e2 = e.subs({sp.Derivative(fq, ee): dTq, sp.Derivative(fp, ee): dTp}).subs({fq: Tq_, fp: Tp_})
e2 = sp.expand(e2)
assert not e2.has(sp.Function('T_Phi'))
nm = {Tq_: 'T_q', Tp_: 'T_p', dTq: r'\dot T_q', dTp: r'\dot T_p'}
rows = []
for mono, lab in [(Tq_ * Tp_, r'T_qT_p'), (Tp_ * dTq, r'T_p\dot T_q'), (Tq_ * dTp, r'T_q\dot T_p'), (dTq * dTp, r'\dot T_q\dot T_p')]:
    c = sp.Poly(e2, Tq_, Tp_, dTq, dTp).coeff_monomial(mono)
    c = sp.factor(sp.together(c))
    rows.append((lab, sp.latex(c)))
chk = sum(sp.Poly(e2, Tq_, Tp_, dTq, dTp).coeff_monomial(m) * m for m in [Tq_*Tp_, Tp_*dTq, Tq_*dTp, dTq*dTp])
assert sp.simplify(chk - e2) == 0
blk = '\\begin{align*}\nK_{S_\\Phi}={}&' + '\\\\\n&+'.join(r'%s\times\Big[%s\Big]' % (lab, c) for lab, c in rows) + '\n\\end{align*}'
print('%%%% KSPhi'); print(blk)
out['KSPhi'] = blk
open('/tmp/claude-0/-home-user-physics/fb9ccdc3-37b3-5823-8617-3fdb5adcb1f8/scratchpad/kl.pkl', 'wb').write(pickle.dumps(out))
