"""Generate app_k3channels.tex from code/A3_third_order/A3_channels.tex.
Purely typographic substitution (no algebra): sympy function notation -> the document's symbols."""
import re
SRC='/home/user/physics/cmb4/code/A3_third_order/A3_channels.tex'
OUT='/home/user/physics/cmb4/doc/app_k3channels.tex'
lines=open(SRC).read().splitlines()
start=next(i for i,l in enumerate(lines) if l.startswith('\\subsection*{List 1'))
body=lines[start:]
C={'1':'x','2':'y','3':'z'}
def math(s):
    s=s.replace('{\\left(\\eta,x,y,z \\right)}','')
    s=s.replace('{\\left(\\eta \\right)}','')
    s=s.replace('{\\left(\\eta_{star} \\right)}','(\\eta_*)')
    s=s.replace('\\eta_{star}','\\eta_*')
    # derivatives
    s=s.replace('\\frac{\\partial}{\\partial \\eta}','\\partial_\\eta ')
    for a in 'xyz': s=s.replace('\\frac{\\partial}{\\partial %s}'%a,'\\partial_%s '%a)
    def d2(m):
        toks=re.findall(r'\\partial (\\eta|[xyz])(?:\^\{(\d)\})?',m.group(1))
        return ''.join('\\partial_{%s}'%v+('^{%s}'%p if p else '') for v,p in toks)+' '
    s=re.sub(r'\\frac\{\\partial\^\{\d\}\}\{((?:\\partial (?:\\eta|[xyz])(?:\^\{\d\})?)+)\}',d2,s)
    # fields
    s=re.sub(r'\bv_\{([123]) ([123])\}',lambda m:'{v^{%s}_{%s}}'%(C[m.group(1)],m.group(2)),s)
    s=re.sub(r'\bh_\{([123])([123]) ([123])\}',lambda m:'h_{%s,%s%s}'%(m.group(3),C[m.group(1)],C[m.group(2)]),s)
    s=re.sub(r'\\operatorname\{om\}_\{([123]) ([123])\}',lambda m:'\\omega_{%s,%s}'%(m.group(2),C[m.group(1)]),s)
    s=re.sub(r'\\operatorname\{Dx\}_\{([12]) ([123])\}',lambda m:'{\\Delta x^{%s}_{%s}}'%(C[m.group(2)],m.group(1)),s)
    s=re.sub(r'\bk_\{([12]) ([123])\}',lambda m:'{k^{%s}_{%s}}'%(C[m.group(2)],m.group(1)),s)
    s=re.sub(r'\bdxo_\{([12]) ([123])\}',lambda m:'{\\delta x^{%s}_{o%s}}'%(C[m.group(2)],m.group(1)),s)
    s=re.sub(r'\bdeta_\{([oe])([12])\}',lambda m:'\\delta\\eta_{%s%s}'%(m.group(1),m.group(2)),s)
    s=re.sub(r'\\(Phi|Psi)_\{([123])\}',lambda m:'\\%s_{%s}'%(m.group(1),m.group(2)),s)
    return s
names={'1':'$L_1$','2':'$\\tfrac12L_2$','3':'$\\tfrac16L_3$'}
out=[]
nterms={}
for l in body:
    m=re.match(r'\\subsection\*\{List (\d) \((\d+) groups, (\d+) terms\)\}',l)
    if m:
        out.append('\\subsection*{Term list of %s (%s groups, %s terms)}'%(names[m.group(1)],m.group(2),m.group(3))); continue
    m=re.match(r'\\paragraph\{\[(ray|obs|emi)\] (.*) \$\|\$ fields:([0-9+]*)\} \\texttt\{where=\w+\}',l)
    if m:
        where={'ray':'ray, $\\int_{\\eta_*}^{\\eta_o}\\dd\\eta$','obs':'observer, $\\eta=\\eta_o$, $\\chi=0$','emi':'emission, $\\eta=\\eta_*$'}[m.group(1)]
        lab=m.group(2).replace('\\^{}2','${}^2$').replace('k1','$\\vk_1$').replace('k2','$\\vk_2$').replace('Dx1','$\\Delta\\vx_1$').replace('Dx2','$\\Delta\\vx_2$')
        fo=m.group(3).replace('+','{+}')
        out.append('\\paragraph{%s; fields of order %s [%s].}'%(lab,fo,where)); continue
    if l.startswith('\\begin{dmath*}'):
        out.append(math(l)); continue
    raise SystemExit('unparsed line: '+l[:80])
hdr=r'''\section{Machine-generated channels of the contracted third-order kernel}
\label{app:k3channels}
{\sloppy
This appendix lists the term groups from which the contracted third-order kernel $\kappa$ of Sec.~\ref{sec:asm-K3c} is assembled, exactly as
produced by the symbolic expansion of $S$ (Sec.~\ref{sec:map3}, App.~\ref{app:L3full}) in the shear form, written in the frame $\nhat=\hat{\mathbf z}$
($x,y$ transverse, $z$ along the observation direction; spatial components are superscripts or second subscripts, perturbative order the
first subscript: $v^x_1$, $h_{2,xz}$, $\omega_{2,z}$, $k^x_1$, $\Delta x^z_2$).  Only the typography has been changed (no algebra).  The three
lists are $L_1$, $\tfrac12L_2$, $\tfrac16L_3$ of $\Lambda=-S=\eps L_1+\tfrac{\eps^2}2L_2+\tfrac{\eps^3}6L_3+\dots$; ``ray'' groups are integrated
along the background ray $\bar\vx=\nhat\chi$, ``observer'' groups are evaluated at $(\eta_o,\mathbf 0)$ and ``emission'' groups at $(\eta_*,\nhat\chi_*)$.
Ray functionals and observer quantities are those of Sec.~\ref{sec:ray} and Sec.~\ref{sec:obsevent} in the limit $\eta_i\to0$:
$k_1=k_{1,\rm init}-\int_\eta^{\eta_o}G_1$, $\Delta\vx_1=\delta\vx_{o1}+\nhat\,\delta\eta_{o1}-\int_\eta^{\eta_o}\vk_1$,
$\tfrac12k_2=\tfrac12k_{2,\rm init}-\int_\eta^{\eta_o}[\tfrac12G_2+\Delta\vx_1\cdot\nabla G_1]$ (plus the transport terms of \eqref{eq:k2def}),
$\tfrac12\Delta\vx_2=\tfrac12\delta\vx_{o2}+\tfrac12\nhat\,\delta\eta_{o2}-\delta\eta_{o1}\vk_1(\eta_o)-\int_\eta^{\eta_o}\tfrac12\vk_2$; the observer clock and
displacement are \eqref{eq:asm-clock3}; $\delta\eta_{e1},\delta\eta_{e2}$ are the emission-time shifts \eqref{eq:deta_e}.

\emph{Contraction.}  For $\zeta=a\,e^{i\vk\cdot\vx}+b\,e^{i\vq\cdot\vx}+c\,e^{-i\vq\cdot\vx}$ with $a^2=b^2=c^2=0$, ${\rm coef}_{abc}$ denotes the
coefficient of $abc$ in an expression; then $\kappa(\vk;\nhat)=\tfrac12\int\!\frac{\dd^3q}{(2\pi)^3}P(q)\,{\rm coef}_{abc}$ \eqref{eq:kappa-def}, applied
to $\Theta^{(3)}=6\,[\tfrac16L_3]+3L_1L_2+L_1^3-3\sky{m_1}m_2-3\sky{m_2}m_1+6\sky{m_1}^2m_1$.  In terms of kernels: a product of three first-order
objects $ABC\to\sum_\pi A[\vk_{\pi1}]B[\vk_{\pi2}]C[\vk_{\pi3}]$ over the permutations of $(\vk,\vq,-\vq)$ (class (i));
a first$\times$second-order product $A_1B_2\to A[\vk]\,B^{\rm s}_2(\vq,-\vq)+A[\vq]\,B^{\rm s}_2(\vk,-\vq)+A[-\vq]\,B^{\rm s}_2(\vk,\vq)$, with
$B^{\rm s}_2(\vq_1,\vq_2)=K_B(\vq_1,\vq_2)+K_B(\vq_2,\vq_1)$ and $(\vq,-\vq)$ the tadpole (class (ii)); a third-order field
$X_3\to K_{3c,X}(\vq)+K_{3c,X}(-\vq)$, the contracted solution \eqref{eq:third_eqs} (class (iii)).  Every group below has been evaluated
numerically by the engine of Sec.~\ref{sec:asm-K3c}; the shear-form lists reproduce the direct-form lists at third order (Sec.~\ref{sec:asm-C13}).\par}

'''
open(OUT,'w').write(hdr+'\n'.join(out)+'\n')
print('written',OUT,len(out),'lines')
