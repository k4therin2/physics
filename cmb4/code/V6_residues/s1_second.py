"""s1_second.py (V6): 1/P^2 poles of C1's second-order pair kernels for |q1| = |q2| = s, total momentum P -> 0.
Kernel symbols (C1 invariant form): k = |q1+q2| (total), q = |q1|, p = |q2|.  We set k -> P, q = p = s and Laurent-expand in P.
Abstract T(eta) := T_Phi(s, eta), first-order ODE T'' = -4T'/eta - s^2 T/3 used to reduce."""
import sys; sys.path.insert(0, '../C1_einstein')
import sympy as sp, pickle
from kernels_api import load
from fourier import eta, k, q, p, mu, TF
K = load(); K2 = K['second']
P, s = sp.symbols('P s', positive=True)
T = sp.Function('T')
def red(e):
    """replace T_Phi(s,eta) by T(eta) and reduce T''."""
    e = e.subs({k: P, q: s, p: s}).doit()
    e = e.subs(TF['Phi'](s, eta), T(eta)).doit()
    e = e.subs(sp.Derivative(T(eta), (eta, 2)), -4 * T(eta).diff(eta) / eta - s**2 * T(eta) / 3)
    return sp.expand(e)
def laurent(e, nmax=0):
    e = red(e)
    e = sp.expand(e)
    out = {}
    for n in range(-6, nmax + 1):
        c = sp.simplify(e.coeff(P, n)) if n != 0 else sp.simplify(e.subs(P, 0) if not e.has(1/P) else sp.expand(e).as_independent(P)[0])
        if c != 0: out[n] = c
    return out
res = {}
for nm in ('K_Q', 'K_S00', 'K_S0z', 'K_SPhi', 'K_om', 'K_vT', 'K_h', 'K_Phi2_0'):
    e = red(K2[nm])
    poly = sp.Poly(sp.expand(e * P**6), P)
    terms = {m[0] - 6: sp.factor(c) for m, c in zip(poly.monoms(), poly.coeffs())}
    neg = {n: c for n, c in terms.items() if n < 0}
    res[nm] = terms
    print(nm, ' negative powers of P:', {n: c for n, c in neg.items()})
Tp = T(eta).diff(eta)
A_Q = sp.factor(res['K_Q'].get(-2, 0))
print('A_Q =', A_Q, '   check = s^2[2T^2 + (T + eta T\')^2]:', sp.simplify(A_Q - s**2 * (2 * T(eta)**2 + (T(eta) + eta * Tp)**2)))
# pole of S_Phi vs A_Q'/eta
dAQ = red(sp.diff(A_Q, eta))
print('pole(K_SPhi) - A_Q\'/eta =', sp.simplify(res['K_SPhi'].get(-2, 0) - dAQ / eta))
# K_Phi2_0 pole vs A_Q(0)/3 (closed form T(0) = -2/3, T'(0) = 0)
print('K_Phi2_0 pole =', res['K_Phi2_0'].get(-2, 0), '  A_Q(0)/3 =', sp.simplify(A_Q.subs(Tp, 0).subs(T(eta), sp.Rational(-2, 3)) / 3))
pickle.dump({'A_Q': sp.srepr(A_Q), 'poles': {nm: {n: sp.srepr(c) for n, c in d.items() if n < 0} for nm, d in res.items()}}, open('s1_poles.pkl', 'wb'))
# closed-form tau: (eta^2 tau)' = -eta^2 A_Q, tau regular
x = sp.Symbol('x', positive=True)
Tc = -2 * (sp.sin(x) - x * sp.cos(x)) / x**3          # x = s eta/sqrt3
eta_x = sp.sqrt(3) * x / s
AQc = s**2 * (2 * Tc**2 + (sp.diff(x * Tc, x))**2)       # (T + eta T') = d(eta T)/d eta = d(x T)/dx
integrand = sp.simplify(eta_x**2 * AQc * sp.sqrt(3) / s)   # d eta = sqrt3/s dx
I = sp.integrate(sp.expand(sp.expand_trig(integrand)), (x, 0, x))
I = sp.simplify(I)
tau = sp.simplify(-I / eta_x**2)
print('tau(x) closed form (x = s eta/sqrt3):', tau)
print('series tau:', sp.series(tau, x, 0, 6))
pickle.dump({'tau_x': sp.srepr(tau)}, open('s1_tau.pkl', 'wb'))
