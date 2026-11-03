"""V5 step 1b / referee B2: the k -> 0 (squeezed, p = |k - q| -> q) structure of C1's second-order kernels.
(i)  leading 1/k^2 and 1/k coefficients of K_Q, K_SPhi, K_Phi2(0), delta2, v2, omega2 at fixed q, GENERIC transfer function T;
(ii) the 1/k^2 parts obey EXACTLY the relations of a time-reparametrisation gauge mode (eta_old = eta + G(eta) P(x)):
       Phi_d = G/eta,  Psi_d = -(G' + G/eta),  delta_d = 4 Phi_d (i.e. -2 Psi_d at eta->0),  V_d = eta Phi_d,  Q_d = eta Phi_d' + 3 Phi_d;
(iii) the q-angle average at FIXED k annihilates every 1/k^2 and 1/k coefficient, so lim_{k->0} int d^3q P(q) K(q, k-q) exists
     and equals the angle-averaged k->0 kernel (C1's tadpole convention).  C1's files are only read."""
import sympy as sp, sys
sys.path.insert(0, '../C1_einstein')
from kernels_api import load, closed_form, T_Phi_closed
from fourier import eta, k, q, p, mu, TF, I
K = load(); K2 = K['second']
T = TF['Phi']
NT = 4
def squeeze(expr, order=0):
    """p -> sqrt(k^2 - 2 k q mu + q^2), T(p) Taylor-expanded about q; Laurent series in k up to k^order."""
    e = expr.doit()
    dp = sp.Symbol('dp')
    reps = {}
    for a in e.atoms(sp.Derivative):
        if a.expr == T(p, eta):
            reps[a] = sum(dp**j/sp.factorial(j)*sp.Derivative(T(q, eta), *( [(q, j)] if j else []), *a.variables) if j else sp.Derivative(T(q, eta), *a.variables) for j in range(NT))
    e = e.xreplace(reps)
    e = e.xreplace({T(p, eta): sum(dp**j/sp.factorial(j)*(sp.Derivative(T(q, eta), (q, j)) if j else T(q, eta)) for j in range(NT))})
    e = e.subs(p, sp.sqrt(k**2 - 2*k*q*mu + q**2)).subs(dp, sp.sqrt(k**2 - 2*k*q*mu + q**2) - q)
    ats = [a for a in e.atoms(sp.Derivative) if a.expr == T(q, eta)] + [T(q, eta)]
    rep = {a: sp.Symbol('tsym%d' % i) for i, a in enumerate(ats)}
    e = sp.expand(e.xreplace(rep))
    r = sp.expand(sp.series(e, k, 0, order + 1).removeO())
    return sp.expand(r.xreplace({v: u for u, v in rep.items()}))
RULE = -4*sp.Derivative(T(q, eta), eta)/eta - q**2*T(q, eta)/3      # first-order Phi ODE: T'' = -4T'/eta - q^2 T/3
def ode_reduce(e):
    e = sp.expand(e.doit())
    for _ in range(6):
        reps = {}
        for a in e.atoms(sp.Derivative):
            if a.expr != T(q, eta): continue
            vc = dict(a.variable_count)
            m, j = vc.get(eta, 0), vc.get(q, 0)
            if m >= 2:
                r = RULE
                if m > 2: r = sp.diff(r, eta, m - 2)
                if j: r = sp.diff(r, q, j)
                reps[a] = r
        if not reps: break
        e = sp.expand(e.xreplace(reps).doit())
    return e
def coeffs(ser):
    return {j: sp.simplify(ser.coeff(k, j)) for j in (-2, -1, 0)}
def angavg(e): return sp.simplify(sp.integrate(sp.expand(e), (mu, -1, 1))/2)

print('=== (i) Laurent coefficients at fixed q (k -> 0) ===')
KQ = coeffs(squeeze(K2['K_Q']))
Kh_qq = sp.simplify(K2['K_h'].doit().subs(p, q))
print('K_Q  k^-2 coefficient:', sp.factor(KQ[-2]))
print('   check K_Q^(-2) = -(3 mu^2 - 1) q^2 K_h(q,q)/4 :', sp.simplify(KQ[-2] + (3*mu**2 - 1)*q**2*Kh_qq/4) == 0)
print('K_Q  k^-1 coefficient:', sp.factor(KQ[-1]))
KS = coeffs(squeeze(K2['K_SPhi']))
print('K_SPhi k^-2 coefficient equals d_eta(K_Q^(-2))/eta (mod first-order ODE):', ode_reduce(KS[-2] - sp.diff(KQ[-2], eta)/eta) == 0)
print('K_SPhi k^-1 coefficient equals d_eta(K_Q^(-1))/eta (mod first-order ODE):', ode_reduce(KS[-1] - sp.diff(KQ[-1], eta)/eta) == 0)
KP0 = coeffs(squeeze(K2['K_Phi2_0']))
print('K_Phi2(0) Laurent:', {j: sp.factor(v) for j, v in KP0.items()})
# eta -> 0 values of Q^(-2): T(q,0) = -2/3, eta T' -> 0 (closed form)
Qd0 = sp.limit(closed_form(KQ[-2]), eta, 0)
print('Q_d(eta->0) =', sp.factor(Qd0), ';  3 Phi_d(0) =', sp.factor(3*KP0[-2]), ' -> equal:', sp.simplify(Qd0 - 3*KP0[-2]) == 0)
Qd1 = sp.limit(closed_form(KQ[-1]), eta, 0)
print('Q^(-1)(eta->0) =', sp.factor(Qd1), ';  3 Phi^(-1)(0) =', sp.factor(3*KP0[-1]), ' -> equal:', sp.simplify(Qd1 - 3*KP0[-1]) == 0)
print('S00 kernel: k^-2, k^-1 coefficients:', [sp.simplify(c) for c in list(coeffs(squeeze(K2['K_S00'])).values())[:2]])
S0z = coeffs(squeeze(K2['K_S0z'], order=1))
print('S0z kernel: k^-2, k^-1 coefficients:', S0z[-2], S0z[-1])

print('\n=== (ii) gauge-mode relations of the divergent parts ===')
print('Let Phi_d solve Phi_d\'\' + 4 Phi_d\'/eta = K_SPhi^(-n) = (K_Q^(-n))\'/eta with Phi_d(0) = K_Phi2(0)^(-n) (regular, unique).')
print('Then r := eta Phi_d\' + 3 Phi_d - Q_d has r\' = eta(Phi_d\'\' + 4Phi_d\'/eta - Q_d\'/eta) = 0 and r(0) = 3Phi_d(0) - Q_d(0) = 0 (checked above):')
print('   => Q_d = eta Phi_d\' + 3 Phi_d for all eta, i.e. with G := eta Phi_d:  Phi_d = G/eta, Psi_d = Phi_d - Q_d = -(G\' + G/eta).')
Ph, Qs, S0 = sp.Function('Phi_d')(eta), sp.Function('Q_d')(eta), sp.Symbol('S0z_d')
Psid = Ph - (eta*sp.diff(Ph, eta) + 3*Ph)
dd = K2['K_delta2_expr'].subs({sp.Function('Phi2')(eta): Ph, sp.Function('Psi2')(eta): Psid, sp.Symbol('S00'): 0}).doit()
dd = sp.expand(dd).coeff(k, 0) if False else sp.expand(dd)
print('delta2 divergent part (S00 has no divergent part; -2 eta^2 k^2 Phi/3 is O(k^0)):', sp.simplify(dd + 2*eta**2*k**2*Ph/3), ' (= 4 Phi_d = -4 G/eta... time-reparametrisation value)')
vz = K2['K_v2z_expr'].subs({sp.Function('Phi2')(eta): Ph, sp.Function('Psi2')(eta): Psid, sp.Symbol('S0z'): 0}).doit()
print('v2_z divergent part (S0z regular):', sp.simplify(vz), '  = i k (eta Phi_d)  => V_d = eta Phi_d = G')
G = sp.Function('G')(eta)
eng = dict(Phi=-G/eta, Psi=sp.diff(G, eta) + G/eta, dT=-G/eta)   # engine (run_all.log, tidal/time-reparam generator with F = G)
print('engine (gauge.py, eta_old = eta + F(eta) P(x), F = -G): Psi + (eta Phi)\' + Phi =', sp.simplify(eng['Psi'] + sp.diff(eta*eng['Phi'], eta) + eng['Phi']).subs(G, -G))

print('\n=== (iii) q-angle average at fixed k (the tadpole limit) ===')
for name in ('K_Q', 'K_SPhi', 'K_Phi2_0', 'K_S00'):
    c = coeffs(squeeze(K2[name]))
    print('%-9s <k^-2 coeff>_mu = %s, <k^-1 coeff>_mu = %s, <k^0 coeff>_mu = %s' % (name, angavg(c[-2]), angavg(c[-1]),
          sp.factor(sp.limit(closed_form(angavg(c[0])*(eta**2 if name == 'K_S00' else 1)), eta, 0)) if name != 'K_Phi2_0' else sp.factor(angavg(c[0]))))
print('   (K_S00 entry is eta^2 x <k^0>, C1: eta^2 tadS00 -> 32/3; K_SPhi entry is its eta->0 limit)')
print('Exact statement for K_Phi2(0): int dOmega_q/(4pi) K_Phi2(0)(q, k - q) at FIXED k =', angavg(K2['K_Phi2_0'].subs(p, sp.sqrt(k**2 - 2*k*q*mu + q**2))), '(k-independent)')
# vector: omega2_i(k) = P_ij(khat) q_j K_om ; leading 1/k term ~ mu ; q-angle average of P_ij q_j mu f(q) = P_ij khat_j q f/3 = 0
Kom = coeffs(squeeze(K2['K_om']))
print('K_om Laurent: k^-2:', sp.simplify(Kom[-2]), '; k^-1:', sp.factor(Kom[-1]))
print('   omega2 ~ P_ij(khat) q_j K_om: the 1/k term is odd in mu -> int dOmega_q q_j mu = (4pi/3) q khat_j, and P_ij khat_j = 0: ensemble tadpole of omega2 vanishes at every k.')
