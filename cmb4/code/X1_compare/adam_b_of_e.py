"""Adam's O(x^2) coefficients b22(e), b13(e) of [1/4C22]/C11 and [1/3ReC13]/C11 (x = k0 eta_o, eta_o = 8, i = 0) and the conversion
of a ratio-coefficient difference db into a k0^6 coefficient difference of 1/3 Re C13 (sigma=1)."""
import sympy as sp
k, x, e = sp.symbols('k x e', positive=True)
o, i = sp.Integer(8), 0
D = o**2 - e**2
C11 = 4*sp.pi/2025*k**4*D**2 - 4*sp.pi/637875*k**6*(o-e)**2*(o+e)*(59*o**3+59*o**2*e-106*o*e**2+114*e**3)
Q22 = 68*sp.pi/6075*k**4*D**2 - 8*sp.pi/637875*k**6*D*(105*i**4+719*o**4-1920*o**2*e**2+2560*o*e**3-1464*e**4)
Q13 = 16*sp.pi/2025*k**4*D**2 - 2*sp.pi/1913625*k**6*D*(2835*i**4+9973*o**4-25440*o**2*e**2+33920*o*e**3-21288*e**4)
def b(Q):
    return sp.simplify(sp.series((Q / C11).subs(k, x / 8), x, 0, 3).removeO().coeff(x, 2))
B22, B13 = b(Q22), b(Q13)
def b22(ev): return float(B22.subs(e, ev))
def b13(ev): return float(B13.subs(e, ev))
def dQ13_k6(db, ev):   # Q13 = C11*(4 + b x^2) with x = 8k: delta(k^6 coeff) = db * 64 * C11_k4coeff
    return float(db * 64 * 4 * sp.pi / 2025 * D.subs(e, ev)**2)
if __name__ == '__main__':
    for ev in (0.5, 1, 2, 4, 6):
        print('e=%-4g Adam b22=%.6f b13=%.6f' % (ev, b22(ev), b13(ev)))
