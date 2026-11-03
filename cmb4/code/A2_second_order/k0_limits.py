"""k0_limits.py -- ANALYTIC K -> 0 behaviour (K = |k1+k2|, i.e. mu -> -1 for |k1| = |k2|) of C1's second-order kernels.
p = |K - q| with p^2 = q^2 - 2 K q c + K^2 (c = Khat.qhat); T_Phi kept abstract and Taylor-expanded about q (symbols
T0 = T_Phi(q), T1 = dT/dk(q), T2 = d^2T/dk^2(q), and eta-derivatives Te0 = dT/deta(q), Te1 = d^2T/dk deta (q), ...).
Output: leading power of K and coefficient for each kernel (k0_limits.log)."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'C1_einstein'))
import sympy as sp
from kernels_api import load
from fourier import eta, k, q, p, TF
Ker = load()['second']
c, Ks = sp.symbols('c K', real=True)
T = TF['Phi']
NT = 4
Tk = sp.symbols('T0:%d' % NT); Tek = sp.symbols('Te0:%d' % NT); Teek = sp.symbols('Tee0:%d' % NT)
def expand_T(expr):
    """replace T(p), dT/deta(p), d2T/deta2(p) by Taylor polynomials around q in (p - q); T(q).. by T0.."""
    dp = sp.Symbol('dp')
    reps = {}
    reps[sp.Derivative(T(p, eta), (eta, 2))] = sum(Teek[j] * dp**j / sp.factorial(j) for j in range(NT))
    reps[sp.Derivative(T(p, eta), eta)] = sum(Tek[j] * dp**j / sp.factorial(j) for j in range(NT))
    reps[T(p, eta)] = sum(Tk[j] * dp**j / sp.factorial(j) for j in range(NT))
    reps2 = {sp.Derivative(T(q, eta), (eta, 2)): Teek[0], sp.Derivative(T(q, eta), eta): Tek[0], T(q, eta): Tk[0]}
    e = expr.doit()
    e = e.subs(reps).subs(reps2)
    return e, dp
out = []
pK = sp.sqrt(q**2 - 2 * Ks * q * c + Ks**2)
for nm in ['K_Q', 'K_om', 'K_vT', 'K_om+K_vT', 'K_S0z', 'K_SPhi', 'K_S00', 'K_Phi2_0', 'K_h']:
    if nm == 'K_om+K_vT':
        expr = Ker['K_om'] + Ker['K_vT']
    else:
        expr = Ker[nm]
    e, dp = expand_T(expr)
    e = e.subs(dp, pK - q).subs(p, pK).subs(k, Ks)
    s = sp.series(e, Ks, 0, 2).removeO()
    s = sp.expand(s)
    lead = None
    for pw in range(-3, 2):
        co = sp.simplify(s.coeff(Ks, pw))
        if co != 0:
            lead = (pw, sp.factor(co)); break
    line = '%-10s  leading K^%s :  %s' % (nm, lead[0] if lead else None, lead[1] if lead else 0)
    print(line); out.append(line)
open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'k0_limits.log'), 'w').write('\n'.join(out) + '\n')
