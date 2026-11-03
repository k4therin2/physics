"""W3: generate doc/app_lterms.tex (complete machine-generated term lists of L1, L2/2, L3/6) from
code/C2_observable/L_terms.pkl.  Extends the generator reproduced in notes/report_W2a.md:
all 359 terms of L3/6 (n = z), all 34 of L2/2 (n = z), the 76-term direct form of L2/2, with readable group titles.
Run:  python3 gen_lterms_full.py ../C2_observable/L_terms.pkl ../../doc/app_lterms.tex"""
import sys, re, pickle
import sympy as sp
from sympy.printing.latex import LatexPrinter

AX = {'1': 'x', '2': 'y', '3': 'z'}
VAR = {'eta': r'\eta', 'x': 'x', 'y': 'y', 'z': 'z'}

def fname(f):
    nm = f.func.__name__; args = f.args
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
        base = fname(e.expr)
        vs = [(str(v), c) for v, c in e.variable_count]
        if vs == [('eta', 1)]: return "{%s}'" % base
        if vs == [('eta', 2)]: return "{%s}''" % base
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
    sgn = '-' if c < 0 else '+'; c = abs(c)
    factors = sp.Mul.make_args(rest)
    def key(f):
        b = f.base if isinstance(f, sp.Pow) else f
        if isinstance(b, sp.Symbol): return (0, str(b))
        if isinstance(b, sp.Function) and re.match(r'(k|Dx)\d_', b.func.__name__): return (1, str(b))
        if isinstance(b, sp.Derivative): return (3, str(b))
        return (2, str(b))
    body = ' '.join(pr.doprint(f) for f in sorted(factors, key=key))
    if c == 1: cs = ''
    elif c.q == 1: cs = str(c.p)
    else: cs = r'\tfrac{%d}{%d}' % (c.p, c.q)
    return sgn, cs + body

def expr_lines(e, per_line):
    terms = sorted(sp.Add.make_args(sp.expand(e)), key=sp.default_sort_key)
    out = []
    for i, t in enumerate(terms):
        s, b = term_tex(t)
        if i == 0 and s == '+': s = ''
        out.append(s + b)
    return [' '.join(out[i:i + per_line]) for i in range(0, len(out), per_line)]

PIECE = [  # raw description fragment -> readable
    ('transverse direction k1 (aberration+deflection)', r'direction $\vk_{1\perp}$'),
    ('lensing 2nd (transverse Dx2)', r'lensing $\Delta\vx_{2\perp}$'),
    ('time delay 2nd (radial Dx2)', r'time delay $\Delta x^z_2$'),
    ('lensing (transverse Dx1)', r'lensing $\Delta\vx_{1\perp}$'),
    ('time delay (radial Dx1)', r'time delay $\Delta x^z_1$'),
    ('transverse k2', r'direction $\vk_{2\perp}$'),
    ('radial k2', r'direction $k^z_2$'),
    ('observer clock 2nd', r'clock $\delta\eta_{o2}$'),
    ('observer clock', r'clock $\delta\eta_{o1}$'),
    ('observer displacement', r'observer displacement $\delta\vx_{o1}$'),
    ('emission-surface shift 2nd', r'emission shift $\delta\eta_{e2}$'),
    ('emission-surface shift', r'emission shift $\delta\eta_{e1}$'),
    ('local', 'local'),
]
WHERE = {'ray': r'ray, $\int_{\eta_*}^{\bar\eta_o}\!\dd\eta$', 'obs': r'observer, at $(\bar\eta_o,\bar\vx_o)$',
         'emi': r'emission, at $(\eta_*,\bar\vx_*)$'}
def readable(name):
    m = re.match(r'\[(\w+)\]\s*(.*?)\s*\|\s*fields:(\S+)', name)
    where, desc, fo = m.group(1), m.group(2), m.group(3)
    parts = []
    for p in desc.split(';'):
        p = p.strip(); sq = p.endswith('^2'); p = p[:-2] if sq else p
        for a, b in PIECE:
            if p == a: p = b; break
        else: raise ValueError(p)
        parts.append(p + (r' (squared)' if sq else ''))
    return where, r' $\times$ '.join(parts), fo.replace('+', r'{+}')

def block(w, lst, per_line, counter):
    tot = 0
    for t in lst:
        where, title, fo = readable(t['name'])
        counter[0] += 1
        tot += t['nterms']
        w(r'\paragraph{(%d) %s; %s; explicit fields of order %s} (%d terms)' % (counter[0], WHERE[where], title, fo, t['nterms']))
        lines = expr_lines(t['expr'], per_line)
        w(r'\begin{align*}')
        w((r'\\' + '\n&\\quad ').join(['&' + lines[0]] + lines[1:]))
        w(r'\end{align*}')
    return tot

def group_table(w, lst, caption, label):
    w(r'\begin{table}[htbp]\centering\small')
    w(r'\caption{%s}\label{%s}' % (caption, label))
    w(r'\begin{tabular}{@{}rlp{0.6\textwidth}lr@{}}')
    w(r'\hline \# & where & content & fields & terms\\ \hline')
    tot = 0
    for i, t in enumerate(lst, 1):
        where, title, fo = readable(t['name'])
        tot += t['nterms']
        w(r'%d & %s & %s & %s & %d\\' % (i, where, title, fo, t['nterms']))
    w(r'\hline & & total & & %d\\ \hline' % tot)
    w(r'\end{tabular}')
    w(r'\end{table}')

HEADER = r"""\section{Complete line-of-sight term lists of \texorpdfstring{$L_1$, $L_2$}{L1, L2} and \texorpdfstring{$L_3$}{L3}}
\label{app:lterms}
{\small\begingroup\sloppy\allowdisplaybreaks
This appendix displays the full perturbative expansion of $\Lambda=-S$ of \eqref{eq:S_eta} through third order for the
direction $\nhat=\hat{\mathbf z}$, as produced by the symbolic expansion of Secs.~\ref{sec:map2}--\ref{sec:map3}
(generator and code: App.~\ref{app:repro}).  It is \emph{complete}: together with the field kernels of
Sec.~\ref{sec:einstein} and App.~\ref{app:kernels}, the ray functionals \eqref{eq:k1def}--\eqref{eq:Dx2def}, the observer
shifts \eqref{eq:deta_o1}--\eqref{eq:dxo} and the emission shifts \eqref{eq:deta_e}, it defines $\Theta^{(3)}$ via
\eqref{eq:theta3}.  The listed expressions are the coefficients of $\eps^m$ in $\Lambda$, i.e.\ $L_1$, $\tfrac12L_2$,
$\tfrac16L_3$.  Notation: $v^i_n$, $\Psi_n$, $\Phi_n$, $\delta_{T,n}$, $\omega_{n,i}$, $h_{n,ij}$ are the order-$n$
Poisson-gauge fields of \eqref{eq:poisson}--\eqref{eq:fluidvars}, with $h_{n,zz}\equiv-h_{n,xx}-h_{n,yy}$ substituted;
$k^i_1,k^i_2,\Delta x^i_1,\Delta x^i_2$ are the ray functionals (functions of $\eta$ on the ray; ``$(\eta_*)$'' marks their
value at the emission end); primes are $\partial_\eta$.  ``Ray'' groups are integrated over $\eta\in[\eta_*,\bar\eta_o]$ with
all fields at $(\eta,\bar\vx(\eta))$; ``observer'' groups are evaluated at $(\bar\eta_o,\bar\vx_o)$; ``emission'' groups at
$(\eta_*,\bar\vx_*)$.  For $\nhat=\hat{\mathbf z}$ the radial direction functional has been eliminated with the null condition
$k^z_1=-(\Phi_1+\Psi_1)$ and $|\nhat|=1$ is used.  A general direction follows by rotation, since $\Lambda$ is a scalar
functional of the fields and $\nhat$ (Sec.~\ref{sec:map3}); the lists for symbolic $\nhat$ (without $|\nhat|=1$ imposed) are
longer ($L_1$: 15, $\tfrac12L_2$: 229, $\tfrac16L_3$: 2484 terms) and are available in machine-readable form (App.~\ref{app:repro}).
No equation of motion has been used in the lists.
"""

def main(pkl, out):
    d = pickle.load(open(pkl, 'rb'))
    W = []; w = W.append
    w('% app_lterms.tex -- MACHINE GENERATED by code/W3_revision/gen_lterms_full.py from code/C2_observable/L_terms.pkl.')
    w('% Do not edit by hand; rerun the generator.')
    w(HEADER)
    # L1
    w(r'\subsection{$L_1$, $\nhat=\hat{\mathbf z}$}')
    c = [0]
    n1 = block(w, d['L_split'][1], 3, c)
    w(r'\noindent i.e.\ $L_1=-\int_{\eta_*}^{\bar\eta_o}I_1\,\dd\eta$ with $I_1$ of \eqref{eq:I1}; on shell this equals \eqref{eq:L1}.')
    # L2 shear
    w(r'\subsection{$\tfrac12L_2$, shear form, $\nhat=\hat{\mathbf z}$: all terms}')
    w(r'\label{app:L2z}')
    c = [0]
    n2 = block(w, d['L_split'][2], 3, c)
    w(r'\noindent Total: %d terms; the named pieces (a)--(g) of Table~\ref{tab:L2pieces} are groups (1), (5), (4)/(2), (3) of this list.' % n2
      if False else r'\noindent Total: %d terms.' % n2)
    group_table(w, d['L'][2], r'Groups of $\tfrac12L_2$ for symbolic $\nhat$ (term counts).', 'tab:L2gen')
    # L2 direct
    w(r'\subsection{$\tfrac12L_2$, direct (redshift) form, $\nhat=\hat{\mathbf z}$}')
    w(r'\label{app:L2direct}')
    c = [0]
    nd = block(w, d['L_direct_split'][2], 3, c)
    w(r'\noindent Total: %d terms; for symbolic $\nhat$: %d terms.' % (nd, sum(t['nterms'] for t in d['L_direct'][2])))
    # L3
    w(r'\subsection{$\tfrac16L_3$, shear form, $\nhat=\hat{\mathbf z}$: all %d terms}' % sum(t['nterms'] for t in d['L_split'][3]))
    w(r'\label{app:L3groups}\label{app:L3full}')
    group_table(w, d['L_split'][3], r'The 32 groups of $\tfrac16L_3$ (shear form, $\nhat=\hat{\mathbf z}$), numbered as in the list below.', 'tab:L3z')
    group_table(w, d['L'][3], r'Groups of $\tfrac16L_3$ for symbolic $\nhat$ (term counts).', 'tab:L3gen')
    c = [0]
    order = sorted(d['L_split'][3], key=lambda t: ({'ray': 0, 'obs': 1, 'emi': 2}[t['where']]))
    n3 = block(w, d['L_split'][3], 2, c)
    w(r'\noindent Total: %d terms.  The direct form of $\tfrac16L_3$ has %d terms ($\nhat=\hat{\mathbf z}$, %d groups) and %d for symbolic $\nhat$; it agrees'
      r' with the shear form modulo the conservation equations (check R-4 of App.~\ref{app:ledger}).' %
      (n3, sum(t['nterms'] for t in d['L_direct_split'][3]), len(d['L_direct_split'][3]), sum(t['nterms'] for t in d['L_direct'][3])))
    w(r'\endgroup}')
    open(out, 'w').write('\n'.join(W) + '\n')
    print('L1 %d, L2/2 %d, direct %d, L3/6 %d terms written' % (n1, n2, nd, n3))

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
