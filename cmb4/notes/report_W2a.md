# Report W2a — writer: second-order map, third-order map, checks, L-term appendix, ledger

Agent W2a (2026-11-03). Files written (doc/ only): `sec_map2.tex`, `sec_map3.tex`, `sec_checks.tex`, `app_lterms.tex` (new),
`app_ledger.tex` (extended), one `\input{app_lterms}` line in `main.tex` (after app_kernels). Not touched: sec_assembly,
sec_mono, sec_unresolved (A2/A3/W2b). Document compiles (pdflatex x2, no errors; remaining undefined refs belong to other
agents' in-progress files: `tab:mono-groups`, `sec:asm-*`).

## Milestones
1. [DONE] sec_map2 "The second-order map": expansion scheme; geodesic eq. in eta; G1; canonical triad to 2nd order (from
   ray_cache.pkl, general index form); direction map k1_init (general n), k2_init (n=z); ray functionals k1, Dx1, k2/2, Dx2/2
   with the CORRECT minus sign (Remark rem:signslip, documenting A1's finding); observer worldline (label at eta->0) with
   dxo1, dxo2; clock deta_o1, deta_o2 (expanded from SPEC §6 by me; identical to A1's pw_map.observer, and checked against the
   dilation family: eta_o = etabar_o/lambda through 2nd order); emission deta_e1, deta_e2; endpoint formula; L2 as 7 named
   pieces (table with meanings and counts 34 / 229), explicit n=z integrand; direct (SW) form explicit n=z (76 terms, from
   L_direct_split, transcribed and checked against the pickle term by term); equivalence (C2 T2 symbolic, A1 numeric, Boubekeur);
   m2, Theta2; rough kernel-type classification (B/C/D) for the assembly.
2. [DONE] sec_map3: ingredients, L3 group table (32 groups / 359 terms n=z; 24 / 2484 general), third-order shear explicit,
   which fields appear (shear form: no omega3, Phi3, Psi3, dT at all except via deta_e; direct form has them — checked in pickle),
   m3, Theta3, qualitative (13) survival classes (variance x map; tadpole; cross-contracted; third-order fields), pointing to
   sec:assembly.
3. [DONE] app_lterms: machine-generated from L_terms.pkl (generator script below): all 34 n=z terms of L2/2 (shear), the 76-term
   direct form, group tables for L2 general n (229) and L3 (n=z and general n), explicit L3 groups with <= 10 terms.
4. [DONE] sec_checks: (1) exact formula (V1, N1, C2), (2) Einstein (C1, V3, KT/PUB/Lu, C4 closed forms, A1 fields),
   (3) dilation (exact family, C2 symbolic, C3 squeezed, A1 tadpole independence 2e-14, V3 zero-mode dependence of C3),
   (4) A1 vs N1 full a_l0 table (24 rows) + map-difference vs N1 uncertainty table + eps^4 residual test + T_o,
   (5) C5 first-order controls, (6) C3 harmonic controls, gaps (axisymmetric/collinear only, no vectors, ratios 2 and 4,
   loops/monochromatic only internally consistent, V2 comparison pending, PUB A.6).
5. [DONE] app_ledger: new C2-7..C2-11, V3-1..7, C4-1..3, C5-1..4, V2-1, N1-2 (sharpened), N1-3, A1-1..7; pending list updated.

## Notes for the lead (not fixed: files not mine)
* sec_observable.tex §observer: worldline written as int_{eta_init}^eta v; A1 showed the label must be at eta->0 (SPEC already
  updated). sec_map2 states the eta->0 convention and remarks on the O(eta_i^2) difference.
* sec_observable.tex consequence (iii) cites "V1-5" for the tetrad-rotation check; the ledger entry is V1-6.
* sec_einstein "Summary" says tadpole-convention irrelevance "to be shown in Sec. checks": sec_checks now cites A1's numerical
  zero-mode independence (plane wave) — a full analytic statement for C_l still belongs to the assembly.
* The clock formula deta_o2 (eq:deta_o2) is my expansion of SPEC §6 (dropping O(eta_i^2)); it matches A1's implementation
  (pw_map.observer, which keeps the eta_i terms) and the dilation family; no separate agent derived it symbolically.

## Generator for app_lterms.tex
Run from code/C2_observable: `python3 gen_lterms.py L_terms.pkl body.tex`, then app_lterms.tex = hand-written header
(section, intro, `\begingroup\sloppy`) + body.tex + `\endgroup`.
```python
"""Generate doc/app_lterms_gen.tex from C2's L_terms.pkl (agent W2a).  Run: python3 gen_lterms.py <pkl> <out.tex>"""
import sys, re, pickle
import sympy as sp
from sympy.printing.latex import LatexPrinter

AX = {'1': 'x', '2': 'y', '3': 'z'}
VAR = {'eta': r'\eta', 'x': 'x', 'y': 'y', 'z': 'z'}

def fname(f):
    nm = f.func.__name__
    args = f.args
    arg_star = (len(args) == 1 and str(args[0]) == 'eta_star')
    m = re.fullmatch(r'(Psi|Phi)(\d)', nm)
    if m: return '\\%s_{%s}' % (m.group(1), m.group(2))
    m = re.fullmatch(r'dT(\d)', nm)
    if m: return r'\delta_{T,%s}' % m.group(1)
    m = re.fullmatch(r'v_(\d)_(\d)', nm)
    if m: return 'v^{%s}_{%s}' % (AX[m.group(1)], m.group(2))
    m = re.fullmatch(r'om_(\d)_(\d)', nm)
    if m: return r'\omega_{%s,%s}' % (m.group(2), AX[m.group(1)])
    m = re.fullmatch(r'h_(\d)(\d)_(\d)', nm)
    if m: return 'h_{%s,%s%s}' % (m.group(3), AX[m.group(1)], AX[m.group(2)])
    m = re.fullmatch(r'(k|Dx)(\d)_(\d)', nm)
    if m:
        base = 'k' if m.group(1) == 'k' else r'\Delta x'
        s = '%s^{%s}_{%s}' % (base, AX[m.group(3)], m.group(2))
        return s + (r'(\eta_*)' if arg_star else '')
    raise ValueError(nm)

class P(LatexPrinter):
    def _print_Function(self, e, exp=None):
        s = fname(e)
        return '(%s)^{%s}' % (s, exp) if exp else s
    def _print_Derivative(self, e):
        f = e.expr
        base = fname(f)
        vs = [(str(v), c) for v, c in e.variable_count]
        if vs == [('eta', 1)]:
            return '{%s}\'' % base
        ops = ''
        for v, c in vs:
            ops += r'\partial_{%s}' % VAR[v] + ('^{%d}' % c if c > 1 else '')
        return ops + base
    def _print_Symbol(self, e):
        s = str(e)
        m = re.fullmatch(r'deta_([oe])(\d)', s)
        if m: return r'\delta\eta_{%s%s}' % (m.group(1), m.group(2))
        m = re.fullmatch(r'dxo(\d)_(\d)', s)
        if m: return r'\delta x^{%s}_{o%s}' % (AX[m.group(2)], m.group(1))
        if s == 'eta_star': return r'\eta_*'
        m = re.fullmatch(r'n(\d)', s)
        if m: return 'n^{%s}' % AX[m.group(1)]
        return LatexPrinter._print_Symbol(self, e)

pr = P(dict(mul_symbol=None))
def term_tex(t):
    c, rest = t.as_coeff_Mul()
    sgn = '-' if c < 0 else '+'
    c = abs(c)
    factors = sp.Mul.make_args(rest)
    # order: symbols, ray functionals, fields
    def key(f):
        b = f.base if isinstance(f, sp.Pow) else f
        if isinstance(b, sp.Symbol): return (0, str(b))
        if isinstance(b, sp.Function) and re.match(r'(k|Dx)\d_', b.func.__name__): return (1, str(b))
        if isinstance(b, sp.Derivative): return (3, str(b))
        return (2, str(b))
    factors = sorted(factors, key=key)
    body = ' '.join(pr.doprint(f) for f in factors)
    if c == 1: cs = ''
    elif c.q == 1: cs = str(c.p)
    else: cs = r'\tfrac{%d}{%d}' % (c.p, c.q)
    return sgn, cs + body

def expr_tex(e, per_line=3):
    terms = sorted(sp.Add.make_args(sp.expand(e)), key=lambda t: sp.default_sort_key(t))
    out = []
    for i, t in enumerate(terms):
        s, b = term_tex(t)
        if i == 0 and s == '+': s = ''
        out.append(s + b)
    lines = [' '.join(out[i:i + per_line]) for i in range(0, len(out), per_line)]
    return lines

def esc(s):
    return s.replace('_', r'\_').replace('^', r'\^{}')

def main(pkl, out):
    d = pickle.load(open(pkl, 'rb'))
    W = []
    w = W.append
    w('% app_lterms_gen.tex -- MACHINE GENERATED by notes-side script gen_lterms.py (agent W2a) from')
    w('% code/C2_observable/L_terms.pkl.  Do not edit by hand.')
    def group_table(lst, caption, label):
        w(r'\begin{table}[htbp]\centering\small')
        w(r'\caption{%s}\label{%s}' % (caption, label))
        w(r'\begin{tabular}{@{}rlp{0.62\textwidth}lr@{}}')
        w(r'\hline \# & where & content & fields & terms\\ \hline')
        tot = 0
        for i, t in enumerate(lst, 1):
            name = t['name']
            m = re.match(r'\[(\w+)\]\s*(.*?)\s*\|\s*fields:(\S+)', name)
            where, content, fo = m.group(1), m.group(2), m.group(3)
            tag = t['tags'].split('|')[0].strip()
            tot += t['nterms']
            w(r'%d & %s & %s%s & %s & %d\\' % (i, where, content.replace('^2', r'$^2$').replace('&', r'\&'),
              (r' [\code{%s}]' % esc(tag)) if tag not in ('', '-') else '', fo, t['nterms']))
        w(r'\hline & & total & & %d\\ \hline' % tot)
        w(r'\end{tabular}')
        w(r'\end{table}')
    # ---------------- L2, n = z, explicit
    w(r'\subsection{$\tfrac12L_2$, shear form, $\nhat=\hat{\mathbf z}$: all 34 terms}')
    w(r'\label{app:L2z}')
    w(r'Each group is listed as \code{[where] name | fields:orders} exactly as stored in \code{L\_terms.pkl[\textquotesingle L\_split\textquotesingle][2]}; '
      r'``ray\textquotedblright\ groups are integrated $\int_{\eta_*}^{\bar\eta_o}\dd\eta$ with all fields at $(\eta,\bar\vx(\eta))$, ``obs\textquotedblright\ groups are evaluated at '
      r'$(\bar\eta_o,\bar\vx_o)$, ``emi\textquotedblright\ groups at $(\eta_*,\bar\vx(\eta_*))$.  Field orders refer to the explicit fields; the ray functionals '
      r'$k_1^i,\Delta x_1^i$ and the shifts $\delta\eta_{o1},\delta\eta_{e1}$ are first order.  Primes are $\partial_\eta$; all $v$ are $v^i_1$ unless indexed otherwise.')
    tot = 0
    for t in d['L_split'][2]:
        tot += t['nterms']
        w(r'\paragraph{\code{%s}} (%d terms)' % (esc(t['name']), t['nterms']))
        w(r'\begin{align*}')
        lines = expr_tex(t['expr'])
        w((r'\\' + '\n&\\quad ').join(['&' + lines[0]] + lines[1:]))
        w(r'\end{align*}')
    w(r'\noindent Total: %d terms ($\nhat=\hat{\mathbf z}$).  The general-$\nhat$ list \code{L[2]} (no $|\nhat|=1$ imposed) has %d terms in %d groups:'
      % (tot, sum(t['nterms'] for t in d['L'][2]), len(d['L'][2])))
    group_table(d['L'][2], r'Groups of $\tfrac12L_2$ for general $\nhat$ (\code{L\_terms.pkl[\textquotesingle L\textquotesingle][2]}).', 'tab:L2gen')
    # ---------------- direct form L2 n=z groups
    w(r'\subsection{$\tfrac12L_2$, direct (redshift) form, $\nhat=\hat{\mathbf z}$}')
    w(r'\label{app:L2direct}')
    tot = 0
    for t in d['L_direct_split'][2]:
        tot += t['nterms']
        w(r'\paragraph{\code{%s}} (%d terms)' % (esc(t['name']), t['nterms']))
        w(r'\begin{align*}')
        lines = expr_tex(t['expr'])
        w((r'\\' + '\n&\\quad ').join(['&' + lines[0]] + lines[1:]))
        w(r'\end{align*}')
    w(r'\noindent Total: %d terms ($\nhat=\hat{\mathbf z}$); general $\nhat$: %d terms in %d groups.' %
      (tot, sum(t['nterms'] for t in d['L_direct'][2]), len(d['L_direct'][2])))
    # ---------------- L3 groups
    w(r'\subsection{$\tfrac16L_3$: group lists}')
    w(r'\label{app:L3groups}')
    group_table(d['L_split'][3], r'Groups of $\tfrac16L_3$, shear form, $\nhat=\hat{\mathbf z}$ (\code{L\_split[3]}); bracketed tags as stored.', 'tab:L3z')
    group_table(d['L'][3], r'Groups of $\tfrac16L_3$, shear form, general $\nhat$ (\code{L[3]}).', 'tab:L3gen')
    w(r'\paragraph{Explicit short groups.}  The groups of Table~\ref{tab:L3z} with at most ten terms ($\nhat=\hat{\mathbf z}$):')
    for t in d['L_split'][3]:
        if t['nterms'] > 10: continue
        w(r'\paragraph{\code{%s}} (%d terms)' % (esc(t['name']), t['nterms']))
        w(r'\begin{align*}')
        lines = expr_tex(t['expr'], per_line=2)
        w((r'\\' + '\n&\\quad ').join(['&' + lines[0]] + lines[1:]))
        w(r'\end{align*}')
    w(r'\noindent The direct form of $\tfrac16L_3$ has %d terms ($\nhat=\hat{\mathbf z}$, %d groups) and %d terms for general $\nhat$.' %
      (sum(t['nterms'] for t in d['L_direct_split'][3]), len(d['L_direct_split'][3]), sum(t['nterms'] for t in d['L_direct'][3])))
    open(out, 'w').write('\n'.join(W) + '\n')

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
```
