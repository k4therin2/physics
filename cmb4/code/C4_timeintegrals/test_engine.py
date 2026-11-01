import sympy as sp, mpmath as mp, sys
sys.path.insert(0, '.')
from trig_integrate import *
t, T = sp.symbols('t T', positive=True)
a, b = sp.symbols('a b', positive=True)
# 1. (cos a t - cos b t)/t
r, _ = integrate_0_to_T((sp.cos(a*t) - sp.cos(b*t))/t, t, T)
print("test1:", sp.simplify(r))
# 2. regular combination with t^-3: (sin(a t) - a t cos(a t))/t^3  (this is T_Phi-like)  -> int_0^T
f = (sp.sin(a*t) - a*t*sp.cos(a*t))/t**3
r, _ = integrate_0_to_T(f, t, T)
r = sp.expand(r.rewrite(sp.cos))
print("test2:", sp.simplify(r))
num = mp.quad(lambda s: (mp.sin(2*s) - 2*s*mp.cos(2*s))/s**3, [0, 1.3])
val = sp.N(r.subs({a: 2, T: sp.Rational(13, 10)}), 20)
print("   numeric check:", num, val)
# 3. product with three frequencies and t^-5: (sin at - a t cos at)(sin bt - bt cos bt) cos(ct) /t^5 ; c symbol
c = sp.Symbol('c', positive=True)
f = (sp.sin(a*t) - a*t*sp.cos(a*t))*(sp.sin(b*t) - b*t*sp.cos(b*t))*sp.cos(c*t)/t**5
r, terms = integrate_0_to_T(f, t, T)
r = sp.expand(r.rewrite(sp.cos))
print("test3 nterms:", len(terms), "frequencies:", sorted(set(w for (n, w) in terms), key=str))
vals = {a: 1.3, b: 0.7, c: 2.1, T: 2.4}
num = mp.quad(lambda s: (mp.sin(1.3*s) - 1.3*s*mp.cos(1.3*s))*(mp.sin(0.7*s) - 0.7*s*mp.cos(0.7*s))*mp.cos(2.1*s)/s**5, [0, 2.4])
val = sp.N(r.subs(vals), 20)
print("   numeric check:", num, val, "imag:", sp.N(sp.im(r.subs(vals))))
