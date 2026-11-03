"""s3_gauge_check.py (V6): the 1/p^2 poles of C1's contracted third-order sources equal (exactly, symbolically) the sources generated
by the residual K=0 Poisson-gauge transformation eta -> eta + eps^2 tau/2 acting on the first-order mode k, with
tau = the second-order pair pole function: KPhi2(q,k,p) = -tau/(eta p^2) + O(1), tau' = -A_Q - 2 tau/eta, A_Q = k^2 [3T^2 + 2 eta T T' + eta^2 T'^2].
Expected shifts (A3 gauge.py formula, x c):  dPhi3 = c[3 tau T' + 6 tau T/eta], dPsi3 = c[3 tau T' + 6 (tau' + tau/eta) T],
dV3 = c[3 tau T_V' + 3 tau' T_V], ddelta3 = c[3 tau T_delta' - 12 tau T_delta/eta], T_V = -(eta^2/2)(T' + T/eta),
T_delta = -(2/3) eta^2 k^2 T - 2 eta T' - 2 T.   Expected c = 2 (bilinear A1 B2 -> 2 A(q) B2_pair + A(k) <B2>)."""
import sys, pickle
sys.path.insert(0, '../C1_einstein')
import sympy as sp
from fourier import eta, k
P = pickle.load(open('s2_poles.pkl', 'rb'))
poles = {nm: sp.sympify(sp.parse_expr(d[-2])) if False else eval('sp.' + 'S(0)') for nm, d in P.items()}
poles = {nm: sp.parse_expr(d[-2]) if -2 in d else 0 for nm, d in P.items()}
poles = {nm: (eval(d[-2], sp.__dict__.copy()) if -2 in d else 0) for nm, d in P.items()}
tau = sp.Function('tau')(eta); T = sp.Function('T')(eta)
Tk, dTk, ddTk = sp.symbols('Tk dTk ddTk'); a0, a1, a2 = sp.symbols('a0 a1 a2')
A = -tau / eta
rep = {a0: A, a1: sp.diff(A, eta), a2: sp.diff(A, eta, 2), Tk: T, dTk: T.diff(eta), ddTk: T.diff(eta, 2)}
AQ = k**2 * (3 * T**2 + 2 * eta * T * T.diff(eta) + eta**2 * T.diff(eta)**2)
def reduce(e):
    e = e.doit()
    for _ in range(4):
        e = e.subs(sp.Derivative(tau, (eta, 3)), sp.diff(-AQ - 2 * tau / eta, eta, 2))
        e = e.subs(sp.Derivative(tau, (eta, 2)), sp.diff(-AQ - 2 * tau / eta, eta))
        e = e.subs(sp.Derivative(tau, eta), -AQ - 2 * tau / eta)
        e = e.doit()
        for n in (4, 3, 2):
            e = e.subs(sp.Derivative(T, (eta, n)), sp.diff(-4 * T.diff(eta) / eta - k**2 * T / 3, eta, n - 2))
        e = e.doit()
    return sp.simplify(sp.expand(e))
c = sp.Symbol('c')
Tp = T.diff(eta); TV = -(eta**2 / 2) * (Tp + T / eta); Td = -sp.Rational(2, 3) * eta**2 * k**2 * T - 2 * eta * Tp - 2 * T
dPhi3 = c * (3 * tau * Tp + 6 * tau * T / eta)
dPsi3 = c * (3 * tau * Tp + 6 * (tau.diff(eta) + tau / eta) * T)
dV3 = c * (3 * tau * sp.diff(TV, eta) + 3 * tau.diff(eta) * TV)
dde3 = c * (3 * tau * sp.diff(Td, eta) - 12 * tau * Td / eta)
expect = {
    'Q': dPhi3 - dPsi3,
    'SPhi': sp.diff(dPhi3, eta, 2) + 4 * sp.diff(dPhi3, eta) / eta + k**2 * dPhi3 / 3,
    'S00': (3 / eta**2) * (dde3 + sp.Rational(2, 3) * eta**2 * k**2 * dPhi3 + 2 * eta * sp.diff(dPhi3, eta) + 2 * dPsi3),
    'S0z': sp.I * k * (-4 * dV3 / eta**2 - 2 * (sp.diff(dPhi3, eta) + dPsi3 / eta)),
}
for nm in ('Q', 'SPhi', 'S00', 'S0z'):
    got = reduce(poles[nm].xreplace(rep))
    ex = reduce(expect[nm])
    sol = sp.solve(sp.expand(got - ex), c, dict=True)
    diff2 = reduce(got - ex.subs(c, 2))
    print('%-5s pole(K3c) - expected(c=2) = %s ;  solve for c: %s' % (nm, diff2, sol))
# C3 pole: C3(k; q=k, mu) ~ 16/(9(1-mu)) = (32 k^2/9) / p^2  vs  c * lim_{eta->0} [3 tau T' + 6 tau T/eta]
x = sp.Symbol('x', positive=True)
tau_x = eval(pickle.load(open('s1_tau.pkl', 'rb'))['tau_x'], sp.__dict__.copy())
s = sp.Symbol('s', positive=True)
tau_eta = tau_x.subs(s, k).subs(x, k * eta / sp.sqrt(3))
Tc = -2 * (sp.sin(x) - x * sp.cos(x)) / x**3
Tc_eta = Tc.subs(x, k * eta / sp.sqrt(3))
lim = sp.limit((3 * tau_eta * sp.diff(Tc_eta, eta) + 6 * tau_eta * Tc_eta / eta), eta, 0)
print('lim_{eta->0} [3 tau T\' + 6 tau T/eta] =', lim, ' ; x2 =', 2 * lim, ' ; C3 pole coefficient (32/9)k^2 =', sp.Rational(32, 9) * k**2)
# closed-form tau satisfies its ODE with A_Q closed form
AQc = k**2 * (3 * Tc_eta**2 + 2 * eta * Tc_eta * sp.diff(Tc_eta, eta) + eta**2 * sp.diff(Tc_eta, eta)**2)
print('closed-form tau ODE residual:', sp.simplify(sp.diff(tau_eta, eta) + 2 * tau_eta / eta + AQc))
