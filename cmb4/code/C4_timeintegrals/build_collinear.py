"""build_collinear.py -- exact limiting forms at the degenerate frequency combinations.
Phi2:  (±k±q±p)/sqrt3 = 0  <=>  collinear:  (i) q+p=k (mu=1, q<k), (ii) q=k+p (mu=1, q>k), (iii) p=k+q (mu=-1).
h2:    k ± (q±p)/sqrt3 = 0  <=>  q+p = sqrt3 k  (the only one inside the physical domain |k-q|<=p<=k+q).
The substitution is made in the INTEGRAND before integrating, so the w=0 terms are integrated exactly."""
import sympy as sp, pickle, sys, time
sys.path.insert(0, '.'); sys.path.insert(0, '../C1_einstein')
from trig_integrate import *
from kernels_api import load, closed_form, eta, k, q, p
K = load(); S = K['second']
etap = sp.Symbol('etap', positive=True)
u1, u2 = S['u1'], S['u2']; W = -3 * sp.sqrt(3) / (eta**4 * k**3)
SPhi = closed_form(S['K_SPhi']).subs(eta, etap)
Kh = closed_form(S['K_h']).subs(eta, etap)
out = {}
cases = {'Phi2_q+p=k': (p, k - q), 'Phi2_q=k+p': (p, q - k), 'Phi2_p=k+q': (p, k + q)}
for nm, (sym, val) in cases.items():
    t0 = time.time()
    src = sp.expand(SPhi.subs(sym, val))
    A, tA = integrate_0_to_T(sp.expand((u1 / W).subs(eta, etap) * src), etap, eta)
    B, tB = integrate_0_to_T(sp.expand((u2 / W).subs(eta, etap) * src), etap, eta)
    Pc = sp.expand((u2 * A - u1 * B).rewrite(sp.cos))
    Pc = sp.expand(sp.re(Pc)) if Pc.has(sp.I) else Pc
    out[nm] = Pc
    print(nm, "done (%.0fs): zero-frequency terms:" % (time.time() - t0), [(n, w) for (n, w) in tA if w == 0])
t0 = time.time()
src = sp.expand(Kh.subs(p, sp.sqrt(3) * k - q))
A, tA = integrate_0_to_T(sp.expand(etap * sp.cos(k * etap) * src / k), etap, eta)
B, tB = integrate_0_to_T(sp.expand(etap * sp.sin(k * etap) * src / k), etap, eta)
Ih = sp.expand((sp.sin(k * eta) / eta * A - sp.cos(k * eta) / eta * B).rewrite(sp.cos))
Ih = sp.expand(sp.re(Ih)) if Ih.has(sp.I) else Ih
out['h2_q+p=sqrt3k'] = Ih
print("h2 resonance done (%.0fs): zero-frequency terms:" % (time.time() - t0), [(n, w) for (n, w) in tA if w == 0], " has log:", Ih.has(sp.log))
pickle.dump({kk: sp.srepr(vv) for kk, vv in out.items()}, open('collinear_closed.pkl', 'wb'))
print("saved collinear_closed.pkl")
