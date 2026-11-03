"""gen_series_small.py (W4): typesets C4's small-eta series (code/C4_timeintegrals/series_small_eta.pkl) for App. closedforms,
each coefficient as (1/denominator) x [polynomial], polynomials written in s = p^2+q^2 and r = p^2 q^2 where shorter."""
import pickle, sympy as sp
d = pickle.load(open('/home/user/physics/cmb4/code/C4_timeintegrals/series_small_eta.pkl', 'rb'))
S = lambda v: sp.sympify(v) if isinstance(v, str) else v
rows = []
for key, lab, nmax in [('Phi2', r'\Phi_2', 4), ('I_h', r'\mathcal K_h', 4), ('Psi2', r'\Psi_2', 4), ('delta2', r'\delta_2', 4), ('V2', r'V_2', 5)]:
    dd = {n: S(v) for n, v in d[key].items()}
    k, p, q = [sp.Symbol(n_, positive=True) for n_ in 'kpq']
    sigma, rho = sp.symbols('s r')
    terms = []
    for n in range(0, nmax + 1):
        c = sp.together(dd.get(n, 0))
        if c == 0: continue
        num, den = sp.fraction(sp.factor(c))
        num = sp.expand(num)
        # rewrite in s, r
        P = sp.Poly(num, p, q)
        try:
            from sympy.polys.polyfuncs import symmetrize
            sym, rem, defs = symmetrize(num.subs({p: sp.sqrt(sigma), q: sp.sqrt(rho)}), formal=True) if False else (None, None, None)
        except Exception:
            pass
        u, v = sp.symbols('u v')
        nn = sp.expand(num.subs({p**2: u, q**2: v}).subs({p: sp.sqrt(u), q: sp.sqrt(v)}))
        from sympy.polys.polyfuncs import symmetrize
        sm, rem, dfs = symmetrize(nn, u, v, formal=True)
        if rem == 0:
            sm = sp.expand(sm.subs({dfs[0][0]: sigma, dfs[1][0]: rho}))
            num_l = sp.latex(sp.collect(sm, k), symbol_names={sigma: 's', rho: 'r'})
        else:
            num_l = sp.latex(num)
        den_l = sp.latex(den)
        terms.append((n, num_l, den_l))
    assert len(terms) > 0
    body = []
    for n, nl, dl in terms:
        body.append(r'\frac{\eta^{%d}}{%s}\Big[%s\Big]' % (n, dl, nl) if n else r'\frac{1}{%s}\Big[%s\Big]' % (dl, nl))
    rows.append(lab + '&=' + ('\\\\\n&\\quad+').join(body) + r'+O(\eta^{%d})' % (nmax + 2))
print('\\begin{align*}\n' + '\\\\\n'.join(rows) + '\n\\end{align*}')
