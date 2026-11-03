"""make_tex.py (A3): machine-generated LaTeX appendix listing every channel of the contracted third-order kernel.
Th3 = 6 L3/6 + 3 L1 (2 L2/2) + L1^3 - 3<m1>m2 - 3<m2>m1 + 6<m1>^2 m1, with C2's lists (shear form, n = z-hat split).
Each L3 group is printed as its integrand/endpoint expression; the contraction rule (polarisation) turns it into K_3c."""
import pickle, sympy as sp
from engine import LT
L = LT['L_split']
out = [r"\section*{Appendix (A3): channels of the contracted third-order map kernel}",
r"""Conventions: $\Lambda=-S=\epsilon L_1+\epsilon^2L_2/2+\epsilon^3L_3/6$; the lists below are $L_1$, $L_2/2$, $L_3/6$ (C2, shear form, frame
$\hat n=\hat z$, fields on the background ray $\bar x=\hat n\chi$ ('ray', integrated $\int_{\eta_*}^{\bar\eta_o}d\eta$), at the observer ('obs') or at
emission ('emi')). Ray functionals: $k_1=k_{1,\rm init}-\int_\eta^{\eta_o}G_1$, $\Delta x_1=\delta x_{o1}+\hat n\delta\eta_{o1}-\int_\eta^{\eta_o}k_1$,
$k_2/2=k_{2o}/2-\int_\eta^{\eta_o}[G_2/2+\Delta x_1\cdot\partial G_1]$, $\Delta x_2/2=\delta x_{o2}/2+\hat n\delta\eta_{o2}/2-\int_\eta^{\eta_o}k_2/2-\delta\eta_{o1}k_1(\eta_o)$;
observer: $\delta\eta_{o1}=-\eta_o^{-1}\int_0^{\eta_o}\eta\Psi_1$, $\delta\eta_{o2}/2=-\eta_o^{-1}[\delta\eta_{o1}^2/2+\delta\eta_{o1}\eta_o\Psi_1(\eta_o)+\int_0^{\eta_o}\eta(\Psi_2/2-\Psi_1^2/2-v_1^2/2+X_1\cdot\partial\Psi_1)]$,
$\delta x_{o2}/2=\int_0^{\eta_o}[v_2/2+X_1\cdot\partial v_1]+\delta\eta_{o1}v_1(\eta_o)$. Contraction: $\kappa(k;\hat n)=\tfrac12\int\frac{d^3q}{(2\pi)^3}P(q)\,{\rm coef}_{abc}$ where
coef$_{abc}$ is the $abc$ component of each expression for $\zeta=a e^{ik\cdot x}+b e^{iq\cdot x}+c e^{-iq\cdot x}$ ($a^2=b^2=c^2=0$): a cubic first-order monomial
$ABC\to\sum_{\rm perm}A[k_{\pi1}]B[k_{\pi2}]C[k_{\pi3}]$ (channel i); $A_1B_2\to A[k]B_2^{AB}(q,-q)+A[q]B_2^{AB}(k,-q)+A[-q]B_2^{AB}(k,q)$,
$B_2^{AB}$ = second-order kernel $K_B(q_1,q_2)+K_B(q_2,q_1)$, $(q,-q)$ = tadpole (channel ii); $X_3\to K^{3c}_X(q)+K^{3c}_X(-q)$ (C1, channel iii).
Th$_3$ = $6\times$[list 3] + $3L_1L_2$ + $L_1^3$ $-3\langle m_1\rangle m_2-3\langle m_2\rangle m_1+6\langle m_1\rangle^2m_1$."""]
for m in (1, 2, 3):
    out.append(r"\subsection*{List %d (%d groups, %d terms)}" % (m, len(L[m]), sum(t['nterms'] for t in L[m])))
    for t in L[m]:
        out.append(r"\paragraph{%s} \texttt{where=%s}" % (t['name'].replace('_', r'\_').replace('^', r'\^{}').replace('|', r'$|$'), t['where']))
        out.append(r"\begin{dmath*}" + sp.latex(t['expr']) + r"\end{dmath*}")
open('A3_channels.tex', 'w').write('\n'.join(out))
print('wrote A3_channels.tex', sum(len(s) for s in out), 'chars')
