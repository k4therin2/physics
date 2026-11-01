"""position_sources.py -- position-space sources: E^{(n)} = n! [eps^n] E = Lin[X_n] + S^{(n)}; dump S^{(n)} (n=2,3)."""
import sympy as sp, pickle
from tpoly import *
d = load_polys('E.pkl'); E = d['E']
names = {(0,0):'00',(0,1):'0x',(0,2):'0y',(0,3):'0z',(1,1):'xx',(1,2):'xy',(1,3):'xz',(2,2):'yy',(2,3):'yz',(3,3):'zz'}
out = {}
for n in (2, 3):
    for (a, b), nm in names.items():
        P = E[a][b].order_part(n) * sp.factorial(n)
        lin = Poly({kk: v for kk, v in P.d.items() if len(kk[2]) == 1})
        src = P - lin
        out[(n, nm, 'lin')] = lin; out[(n, nm, 'src')] = src
save_polys(out, 'position_sources.pkl')
with open('position_sources.tex', 'w') as fh:
    fh.write("% E^{(n)}_{mu nu} = Lin[X_n] + S^{(n)}_{mu nu};  Lin is the first-order operator acting on the order-n fields.\n")
    fh.write("% Convention: X = eps X1 + eps^2 X2/2 + eps^3 X3/6;  S^{(2)} is quadratic in X1, S^{(3)} = bilinear(X1,X2) + cubic(X1).\n")
    for n in (2, 3):
        for (a, b), nm in names.items():
            for part in ('lin', 'src'):
                ex = out[(n, nm, part)].to_sympy()
                fh.write("\n\\subsection*{$n=%d$, component %s, %s (%d terms)}\n\\begin{dmath}\n%s\n\\end{dmath}\n" % (n, nm, part, len(sp.Add.make_args(ex)), sp.latex(ex)))
    fh.write("""
\\section*{Projections (Fourier, $\\hat k = \\hat z$)}
Traceless scalar: $-k^2 Q_n + [S_{zz} - (S_{xx}+S_{yy})/2] = 0$, $Q_n = \\Phi_n - \\Psi_n$; position space
$Q_n = -\\tfrac32 \\nabla^{-2}(\\partial_i\\partial_j\\nabla^{-2} - \\delta_{ij}/3) S^{(n)}_{ij}$.
00: $\\delta_n = \\tfrac{\\eta^2}{3} S^{(n)}_{00} - \\tfrac23\\eta^2 k^2\\Phi_n - 2\\eta\\Phi_n' - 2\\Psi_n$.
0i longitudinal: $v_{n,z} = -\\tfrac{\\eta^2}{4}[2ik(\\Phi_n' + \\Psi_n/\\eta) + S^{(n)}_{0z}]$.
0i transverse: $\\tfrac12\\nabla^2\\omega_{n,i} - 4(\\omega_{n,i} + v^T_{n,i})/\\eta^2 = S^{(n)T}_{0i}$ (in Fourier: $-k^2\\omega/2 - 4(\\omega+v^T)/\\eta^2 = S^T_{0i}$).
ij vector: $-\\tfrac{ik}{2}(\\omega_{n,x}' + 2\\omega_{n,x}/\\eta) + S^{(n)}_{xz} = 0$.
TT: $h_{n,p}'' + 2h_{n,p}'/\\eta + k^2 h_{n,p} = -(S^{(n)}_{xx}-S^{(n)}_{yy})/2 \\cdot 2$, i.e. $(h_{xx}-h_{yy})'' + \\dots = -2(S_{xx}-S_{yy})$ (the TT operator is $\\tfrac12(h''+2\\mathcal H h' - \\nabla^2 h)_{ij}$).
Trace: $\\Phi_n'' + 4\\Phi_n'/\\eta + k^2\\Phi_n/3 = -(S^{(n)}_{xx}+S^{(n)}_{yy}+S^{(n)}_{zz} - S^{(n)}_{00})/6 + Q_n'/\\eta - k^2 Q_n/3$.
""")
print({kk: v.nterms() for kk, v in out.items() if kk[2] == 'src'})
