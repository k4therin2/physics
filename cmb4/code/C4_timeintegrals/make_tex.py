"""make_tex.py -- generates C4_closed_forms.tex (structure in the main text, full forms in the appendix) from the pickles."""
import sympy as sp, pickle, sys
sys.path.insert(0, '.'); sys.path.insert(0, '../C1_einstein')
from trig_integrate import *
from kernels_api import load, eta, k, q, p
L = {'Cin': Cin}; s3 = sp.sqrt(3)
S = load()['second']
Tq, Tp, dTq, dTp = sp.symbols(r'T_q T_p \dot{T}_q \dot{T}_p')
def abstract(e):
    e = e.doit()
    e = e.subs({sp.Derivative(sp.Function('T_Phi')(q, eta), eta): dTq, sp.Derivative(sp.Function('T_Phi')(p, eta), eta): dTp})
    e = e.subs({sp.Function('T_Phi')(q, eta): Tq, sp.Function('T_Phi')(p, eta): Tp})
    return sp.expand(e)
def lx(e): return sp.latex(e).replace(r'\operatorname{Cin}', r'\mathrm{Cin}').replace(r'\operatorname{Si}', r'\mathrm{Si}')
def kernel_tex(e):
    e = abstract(e)
    P = sp.Poly(e, Tq, Tp, dTq, dTp)
    parts = []
    for mono, c in P.terms():
        m = sp.Mul(*[s**n for s, n in zip((Tq, Tp, dTq, dTp), mono)])
        parts.append(r"\left[%s\right] %s" % (lx(sp.factor(c)), lx(m)))
    return r" \\ &\quad + ".join(parts)
ser = {nm: {n: sp.parse_expr(c) for n, c in d.items()} for nm, d in pickle.load(open('series_small_eta.pkl', 'rb')).items()}
def series_tex(nm, nmax):
    out = []
    for n in range(nmax + 1):
        c = ser[nm].get(n, 0)
        if c != 0: out.append(r"\left[%s\right]\eta^{%d}" % (lx(sp.factor(c)), n))
    return r" \\ &\quad + ".join(out) + r" + O(\eta^{%d})" % (nmax + 2)
comp = {}
for nm in ('phi2', 'h2'):
    try:
        d = pickle.load(open(nm + '_compact.pkl', 'rb'))
        comp[nm] = ({kk: sp.parse_expr(v, local_dict=L) for kk, v in d['sf'].items()},
                    {sp.parse_expr(a, local_dict=L): sp.parse_expr(b, local_dict=L) for a, b in d['trig'].items()})
    except FileNotFoundError:
        pass
xk, xq, xp = sp.symbols('x_k x_q x_p')
def trig_tex(trig, xmap):
    lines = []
    for m, c in sorted(trig.items(), key=lambda t: str(t[0])):
        mm = m.subs(xmap)
        c = sp.expand(c)
        pw = sorted({sp.Poly(sp.expand(t * eta**30), eta).degree() - 30 for t in sp.Add.make_args(c)})
        terms = []
        for n in sorted(pw):
            cn = sp.expand(c * eta**(-n)).subs(eta, 0) if False else sp.Add(*[t for t in sp.Add.make_args(c) if sp.Poly(sp.expand(t * eta**30), eta).degree() - 30 == n]) / eta**n
            terms.append(r"\left[%s\right]\eta^{%d}" % (lx(sp.factor(sp.expand(cn))), n))
        lines.append(r"%s\Big\{ %s \Big\}" % (lx(mm), r" \\ &\qquad + ".join(terms)))
    return r" \\ &+ ".join(lines)

doc = r"""\documentclass[10pt]{article}
\usepackage[margin=0.8in]{geometry}
\usepackage{amsmath,amssymb}
\allowdisplaybreaks
\title{C4: closed-form time integrals of the second-order radiation-era solutions\\ (Poisson gauge, $a\propto\eta$)}
\date{2026-11-03}
\begin{document}\maketitle
\section{Setup and result in one paragraph}
Conventions of SPEC.md and report\_C1: $X=\epsilon X_1+\epsilon^2X_2/2$, $X_{2,\mathbf k}=\int\frac{d^3q}{(2\pi)^3}K_X(q,|\mathbf k-\mathbf q|;\eta)\,\zeta_{\mathbf q}\zeta_{\mathbf k-\mathbf q}$,
$p=|\mathbf k-\mathbf q|$, $x_a=a\eta/\sqrt3$, $T_\Phi(x)=-2(\sin x-x\cos x)/x^3$.
The only non-algebraic second-order time integrals are
\begin{align}
\Phi_2&=P(\eta)+K_{\Phi_2,0}\,\hat g(k\eta),\qquad P=\int_0^\eta G_\Phi(\eta,\eta')S_\Phi(\eta')d\eta',\qquad \hat g=3u_1=-\tfrac32T_\Phi(x_k),\\
h_{2,ij}&=\Lambda_{ij,kl}q_kq_l\,I_h,\qquad I_h=\int_0^\eta \frac{\eta'}{\eta}\frac{\sin k(\eta-\eta')}{k}K_h(\eta')d\eta' = \frac{4}{k^2}I_{\rm RD}^{\rm KT}\!\left(\tfrac pk,\tfrac qk,k\eta\right).
\end{align}
Both are obtained in CLOSED FORM: elementary functions plus $\mathrm{Si}(w\eta)$ and the entire function
$\mathrm{Cin}(z)=\int_0^z(1-\cos t)\,dt/t=\gamma_E+\ln|z|-\mathrm{Ci}(|z|)$ at the frequencies
\[\Phi_2:\ w_1=\tfrac{k+q+p}{\sqrt3},\ w_2=\tfrac{k+q-p}{\sqrt3},\ w_3=\tfrac{k-q+p}{\sqrt3},\ w_4=\tfrac{q+p-k}{\sqrt3}\ (\text{all}\ \ge0);\qquad
h_2:\ k\pm\tfrac{q\pm p}{\sqrt3}.\]
All other second-order fields are algebraic in $\Phi_2,\Phi_2'$ and first-order transfer functions (Sec.~\ref{sec:alg}); hence
\textbf{nothing remains as an unevaluated integral at second order}. Because Cin is entire and even, the forms need no
$|\cdot|$, principal values or branch choices; rewriting $\mathrm{Cin}\to\mathrm{Ci}$ reproduces Kohri--Terada (KT) eq.~(22) with
its $\log|\ldots|$ terms.

\section{Method}
The integrands are finite sums $\sum c_{n,w}(k,q,p)\,\eta'^n e^{iw\eta'}$, $-7\le n\le0$. Each term is integrated by
$I_n=[t^{n+1}e^{iwt}]/(n+1)-\frac{iw}{n+1}I_{n+1}$, $I_{-1}=-\mathrm{Cin}(w\eta)+i\,\mathrm{Si}(w\eta)+\ln(\eta/\varepsilon)$, on $[\varepsilon,\eta]$.
All $\varepsilon^{-m}$ and $\ln\varepsilon$ pieces cancel identically (symbolic check) --- this is the regularity of the total
integrand at $\eta'=0$. The source of $\Phi_2$ is
\begin{align}
S_\Phi&=K_{S\Phi}= %s
\end{align}
with $T_a\equiv T_\Phi(x_a)$, $\dot T_a\equiv\partial_\eta T_\Phi(x_a)$, and $G_\Phi=[u_1(\eta')u_2(\eta)-u_1(\eta)u_2(\eta')]/W(\eta')$,
$u_1=(\sin x-x\cos x)/x^3$, $u_2=(\cos x+x\sin x)/x^3$, $W=-3\sqrt3/(\eta^4k^3)$.

\section{Structure of $\Phi_2$}
\begin{align}
P&=u_2(x_k)\,\mathcal A(\eta)-u_1(x_k)\,\mathcal B(\eta),\qquad
\mathcal A=\int_0^\eta\frac{u_1}{W}S_\Phi,\quad \mathcal B=\int_0^\eta\frac{u_2}{W}S_\Phi .
\end{align}
The special-function skeleton is governed by ONE coefficient,
\begin{align}
\mathcal D(k,q,p)&=\frac{(k^2+p^2+q^2)\,R}{8k^2p^3q^3},\quad
R=3k^6-3k^4(p^2+q^2)-3k^2(p^4+q^4)-2k^2p^2q^2+3(p^6+q^6)-3p^2q^2(p^2+q^2),\\
\mathcal A&\supset \mathcal D\,[\mathrm{Si}(w_1\eta)-\mathrm{Si}(w_2\eta)-\mathrm{Si}(w_3\eta)-\mathrm{Si}(w_4\eta)],\qquad
\mathcal B\supset %s
\end{align}
the rest of $\mathcal A,\mathcal B$ being trigonometric polynomials times Laurent polynomials in $\eta$ (full form: Appendix~\ref{app:phi2}).
Only the $n=0$ monomials produce $1/w$ factors; they become removable $0/0$ at the collinear edges
$w_2=0$ ($\mu=-1$), $w_3=0$ ($\mu=1,q>k$), $w_4=0$ ($\mu=1,q<k$). Exact collinear limiting forms are obtained by substituting
$p=k+q,\ q-k,\ k-q$ in the integrand before integrating (\texttt{collinear\_closed.pkl}); zero-frequency terms occur only at
$n=-4,-2,0$, so no $\ln\eta$ appears in collinear $\Phi_2$. No $\sqrt3$-resonance exists for $\Phi_2$ (sound speed enters both $G_\Phi$ and the source).

\section{Tensor control: $h_2$}
$I_h=[\sin(k\eta)A_h(\eta)-\cos(k\eta)B_h(\eta)]/\eta$ with $A_h,B_h$ = Cin/Si at $k\pm(q\pm p)/\sqrt3$ plus elementary terms
(Appendix~\ref{app:h2}). It equals $(4/k^2)\times$KT(22) [with $|\cdot|$ inside its log] to $10^{-15}$ (double) / 25 digits
(exact inputs). Resonance $q+p=\sqrt3k$: harmless in Cin form; the exact resonant form (separately integrated) contains
an explicit $\ln\eta$. The other resonance $|q-p|=\sqrt3k$ lies outside the physical domain.

\section{All other fields (algebraic)}\label{sec:alg}
\begin{align}
\Psi_2&=\Phi_2-Q_2,\qquad Q_2=K_Q= %s\\
\delta_2&=\tfrac{\eta^2}{3}K_{S00}-\tfrac23\eta^2k^2\Phi_2-2\eta\Phi_2'-2\Psi_2,\qquad
\delta_{T2}=\tfrac14\delta_2-\tfrac3{16}T_\delta(q)T_\delta(p),\\
V_2&=-\tfrac{\eta^2}{4}\Big[2\big(\Phi_2'+\Psi_2/\eta\big)+\frac{K_{S0z}}{ik}\Big],\quad v_{2,i}=\partial_iV_2+v^T_{2,i},\\
\omega_{2,i}&=P_{ij}(\hat k)q_j\,K_\omega,\qquad K_\omega= %s,\\
v^T_{2,i}&=P_{ij}(\hat k)q_j\,K_{vT}\quad(\text{explicit in second\_order.pkl; } K_{vT}=-(\eta^2k^2/8+1)K_\omega-\tfrac{\eta^2}{4}K_{S0\perp}).
\end{align}
$K_{S00},K_{S0z}$ are C1's (closed in $T_\Phi$); $\Phi_2'$ is the closed-form derivative ($G_\Phi(\eta,\eta)=0$).

\section{Limits}
\subsection*{$\eta\to0$ (Frobenius + constraints; checked against the closed forms at 120 digits)}
\begin{align}
\Phi_2&=%s\\
I_h&=%s\\
\Psi_2&=%s\\
\delta_2&=%s\\
V_2&=%s
\end{align}
\subsection*{$\eta\to\infty$}
With $\Lambda=\ln\frac{(k+q+p)(q+p-k)}{(k+q-p)(k+p-q)}$:
\begin{align}
\Phi_2&\to\frac{1}{\eta^2}\Big\{C_L\big[\Lambda\cos x_k+\pi\sin x_k\big]+C_0\cos x_k\Big\}+O(\eta^{-3}),\qquad C_L=-\frac{3(k^2+p^2+q^2)R}{8k^4p^3q^3},\\
C_0&=\frac{3\left[-3k^6-3k^4(p^2+q^2)+3k^2(p^4+q^4)+2k^2p^2q^2+3(p^6+q^6)-3p^2q^2(p^2+q^2)\right]}{2k^4p^2q^2},\\
I_h&\to\frac{4}{k^2}\,{\rm KT}(25)\quad(\sin k\eta\ \text{log term}+\Theta(q+p-\sqrt3k)\,\pi\cos k\eta\ \text{term}).
\end{align}
Next orders ($\eta^{-3}$ for $\Phi_2$, $\eta^{-2}$ for $I_h$) in \texttt{asymptotics.pkl}; verified up to $k\eta=10^5$ (remainder
scales as the first dropped power). The expansion is non-uniform as $q\to0$ or $p\to0$ and near the collinear edges / $q+p=\sqrt3k$.

\section{Verification}
\begin{itemize}
\item Quadrature arbiter (Gauss--Legendre, 50--60 digits) vs closed forms: $P$ at 11 points incl.\ collinear and $\sqrt3$ points:
$\le3.5\times10^{-40}$; $I_h$: 25 digits; $\Phi_2,\Phi_2',I_h$ at 10 random points in \texttt{verify\_module.py}.
\item END-TO-END: all C4 fields inserted into C1's full second-order Einstein tensor (10 components) and the 4 conservation
equations for a plane-wave pair: residual/max term $\le5.5\times10^{-50}$ at 3 points ($k\eta=2.5,\,9.6,\,20.1$).
\item ODE residuals by substitution; small-$\eta$ series = Frobenius series; KT (22),(25) controls.
\item Fast double-precision evaluators (\texttt{second\_order\_closed.py}) vs 50-digit reference: see report.
\end{itemize}
\appendix
\section{Full form of $P$ (particular integral of $\Phi_2$)}\label{app:phi2}
%s
\section{Full form of $I_h$}\label{app:h2}
%s
\end{document}
"""
def full_tex(nm, xmap, skel_label):
    if nm not in comp: return "(not generated)"
    sf, trig = comp[nm]
    s = (r"Frequencies: $w_1=(k+q+p)/\sqrt3,\ w_2=(k+q-p)/\sqrt3,\ w_3=(k-q+p)/\sqrt3,\ w_4=(q+p-k)/\sqrt3$. " if nm == 'phi2' else
         r"Frequencies: $w_1=k+(q+p)/\sqrt3,\ w_2=k+(q-p)/\sqrt3,\ w_3=k-(q-p)/\sqrt3,\ w_4=k-(q+p)/\sqrt3$ (only $w_4$ can vanish). ")
    s += r"Special-function part (function, frequency):\begin{align*}" + \
        r" \\ ".join(r"&[%s,\ w_{%d}]:\ %s" % (f, j + 1, lx(c)) for (f, j), c in [((kk.split(',')[0].strip("('\" "), int(kk.split(',')[1].strip(")' "))), v) for kk, v in sf.items()]) + r"\end{align*}"
    s += r"Elementary part ($x_a=a\eta/\sqrt3$%s):\begin{align*}&%s\end{align*}" % (skel_label, trig_tex(trig, xmap))
    return s
x = lambda a: a * eta / s3
xmap_phi = {sp.sin(x(k)): sp.sin(xk), sp.cos(x(k)): sp.cos(xk), sp.sin(x(q)): sp.sin(xq), sp.cos(x(q)): sp.cos(xq), sp.sin(x(p)): sp.sin(xp), sp.cos(x(p)): sp.cos(xp)}
yk = sp.Symbol(r'k\eta')
xmap_h = dict(xmap_phi); xmap_h.update({sp.sin(k * eta): sp.sin(yk), sp.cos(k * eta): sp.cos(yk)})
Bsk = r"-\tfrac{\mathcal D}{1}\,[\ldots]\ \text{(see Appendix)}"
if 'phi2' in comp:
    pass
tex = doc % (kernel_tex(S['K_SPhi']), Bsk, kernel_tex(S['K_Q']), kernel_tex(S['K_om']),
             series_tex('Phi2', 4), series_tex('I_h', 4), series_tex('Psi2', 4), series_tex('delta2', 4), series_tex('V2', 5),
             full_tex('phi2', xmap_phi, ''), full_tex('h2', xmap_h, r', and $\sin k\eta,\cos k\eta$'))
open('C4_closed_forms.tex', 'w').write(tex)
print("written C4_closed_forms.tex, %d chars" % len(tex))
